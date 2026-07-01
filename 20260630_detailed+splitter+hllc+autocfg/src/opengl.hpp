/*
  Copyright © Cambridge Numerical Solutions Ltd 2013
*/
#pragma once
// Comment me to disable OpenGL output
//#define GLOUTPUT
#ifdef GLOUTPUT

#include <stdio.h>
#include <stdlib.h>
#include <GL/glew.h>
#include <GL/freeglut.h>
#include <cutil.h>
#include <cuda_runtime_api.h>
#include <cuda_gl_interop.h>
#include <thrust/device_vector.h>
#include <thrust/scan.h>
#include <vector>
#include <cutil.h>
#include "render.hpp"


// OpenGL pixel buffer object and texture //
GLuint gl_PBO, gl_Tex;

real* plot,* plot_device;
unsigned int* cmap_rgba, *plot_rgba; // rgba arrays for plotting
unsigned int* plot_rgba_data;
// fun stuff
int cursor = GLUT_CURSOR_CROSSHAIR;cubehelix
int fullscreen = 0;

int WriteBitmap(const char* name, unsigned char* buff, unsigned width, unsigned height, bool flip = false);

__global__ void numSegmentsKernel(real* data, int* numSegments, int ni, int nj, int width, int height, real value) {
	const int nni = ni - 1, nnj = nj - 1;
	const int offset = blockDim.x * blockIdx.x + threadIdx.x;
	const int i = offset % nni, j = offset / nni;

	if (i < nni && j < nnj) {
		if (i < width - 1 && j < height - 1) {
			int signature = (data[j * ni + i] < value) * 8 + (data[j * ni + i + 1] < value) * 4 + (data[(j + 1) * ni + i + 1] < value) * 2 + (data[(j + 1) * ni + i] < value);
			if (signature == 0xA || signature == 0x5) {
				numSegments[offset] = 2;
			} else if (signature != 0xF && signature != 0x0) {
				numSegments[offset] = 1;
			} else {
				numSegments[offset] = 0;
			}
		} else {
			numSegments[offset] = 0;
		}
	}
}
__global__ void generateSegmentsKernel(real* data, int*  segmentIDs, int ni, int nj, int width, int height, real value, float* segments) {
	const int nni = ni - 1, nnj = nj - 1;
	const int offset = blockDim.x * blockIdx.x + threadIdx.x;
	const int i = offset % nni, j = offset / nni;

	if (i < nni && j < nnj) {
		if (i < width - 1 && j < height - 1) {
			int signature = (data[j * ni + i] < value) * 8 + (data[j * ni + i + 1] < value) * 4 + (data[(j + 1) * ni + i + 1] < value) * 2 + (data[(j + 1) * ni + i] < value);
			const int id = segmentIDs[offset];
			switch (signature) {
				case 0x1:
				case 0xE:
					segments[4 * id + 0] = i;
					segments[4 * id + 1] = j + 0.5;
					segments[4 * id + 2] = i + 0.5;
					segments[4 * id + 3] = j + 1;
					break;
				case 0x2:
				case 0xD:
					segments[4 * id + 0] = i + 0.5;
					segments[4 * id + 1] = j + 1;
					segments[4 * id + 2] = i + 1;
					segments[4 * id + 3] = j + 0.5;
					break;
				case 0x3:cubehelix
				case 0xC:
					segments[4 * id + 0] = i;
					segments[4 * id + 1] = j + 0.5;
					segments[4 * id + 2] = i + 1;
					segments[4 * id + 3] = j + 0.5;
					break;
				case 0x4:
				case 0xB:
					segments[4 * id + 0] = i + 0.5;
					segments[4 * id + 1] = j;
					segments[4 * id + 2] = i + 1;
					segments[4 * id + 3] = j + 0.5;
					break;
				case 0x6:
				case 0x9:
					segments[4 * id + 0] = i + 0.5;
					segments[4 * id + 1] = j;
					segments[4 * id + 2] = i + 0.5;
					segments[4 * id + 3] = j + 1;
					break;
				case 0x7:
				case 0x8:
					segments[4 * id + 0] = i;
					segments[4 * id + 1] = j + 0.5;
					segments[4 * id + 2] = i + 0.5;
					segments[4 * id + 3] = j;
					break;
				case 0x5:
				case 0xA:
					segments[4 * id + 0] = i;
					segments[4 * id + 1] = j + 0.5;
					segments[4 * id + 2] = i + 0.5;
					segments[4 * id + 3] = j;
					segments[4 * (id + 1) + 0] = i + 0.5;
					segments[4 * (id + 1) + 1] = j + 1;
					segments[4 * (id + 1) + 2] = i + 1;
					segments[4 * (id + 1) + 3] = j + 0.5;
					break;
				case 0x0:
				case 0xF:
					break;
			}
		}
	}
}

// MARCHING SQUARES
std::vector<float> getMarchingSquares(real* data, int ni, int nj, int width, int height, real value) {
	int* numSegments,* segmentIDs;
	int nni = ni - 1, nnj = nj - 1;
	cudaMalloc(&numSegments, sizeof(int) * nni * nnj);
	cudaMalloc(&segmentIDs, sizeof(int) * nni * nnj);

	dim3 blockDim(256, 1, 1);
	dim3 gridDim((nni * nnj + 255) / 256, 1, 1);

	numSegmentsKernel<<<gridDim, blockDim>>>(data, numSegments, ni, nj, width, height, value);
	cudaThreadSynchronize();
	checkForError();

	thrust::exclusive_scan(thrust::device_ptr<int>(numSegments), thrust::device_ptr<int>(numSegments + nni * nnj), thrust::device_ptr<int>(segmentIDs));
	cudaThreadSynchronize();

	cudaFree(numSegments);

	int maxSegmentID;
	cudaMemcpy(&maxSegmentID, segmentIDs + nni * nnj - 1, sizeof(int), cudaMemcpyDeviceToHost);

	float* segments;
	cudaMalloc(&segments, sizeof(float) * 4 * (maxSegmentID));
	generateSegmentsKernel<<<gridDim, blockDim>>>(data, segmentIDs, ni, nj, width, height, value, segments);
	cudaThreadSynchronize();

	std::vector<float> segments_h(4 * (maxSegmentID));
	cudaMemcpy(&segments_h[0], segments, sizeof(float) * 4 * (maxSegmentID), cudaMemcpyDeviceToHost);
	cudaFree(segments);
	cudaFree(segmentIDs);
	return segments_h;
}

class OpenGLOutputter {
public:
	static OpenGLOutputter* instance;
	int plotVariable;

public:
	int argc;
	char** argv;
	int width, height;
	const int ni, nj;
	int viewportWidth, viewportHeight, textHeight;
	bool paused, showBoundaries;
	ColourMode colourMode;
	bool autoscale;
	limits lastLimits;

	OpenGLOutputter(const OpenGLOutputter& copy);

	// sizes for the progress bar
	const int textHeightPx;
	const typename Mesh<GPU>::type& grid;

	dim3 blockDim, gridDim;

public:
	typename Mesh<GPU>::type* gridToRender;

	OpenGLOutputter(int argc_, char** argv_, typename Mesh<GPU>::type& grid) : argc(argc_), argv(argv_), ni(grid.activeNx()), nj(grid.activeNy()), gridToRender(NULL), grid(grid),
	textHeightPx(30),
	plotVariable(SCHLIEREN), colourMode(GREYSCALE), paused(!true), showBoundaries(false), autoscale(true) {
		instance = this;
		width = ni;

		blockDim = dim3(16, 16, 1);
		gridDim  = dim3((ni + 15) / 16, (nj + 15) / 16, 1);
	}

	static OpenGLOutputter& getInstance() {
		return *instance;
	}

	void operator()() {
		// Are we fullscreen?
		fullscreen = false;// params->fullscreen;
		// Window size
		height = nj;
		float window_xsize = width;// params->N[0] * params->dx[0];
		float window_ysize = height;// params->N[1] * params->dx[1];
		int window_xpixels, window_ypixels;
		int window_maxpixels = 1200;
		if (window_xsize >= window_ysize) {
			window_xpixels = window_maxpixels;
			window_ypixels = (int)((window_ysize / window_xsize) * window_maxpixels);
		} else {
			window_ypixels = window_maxpixels;
			window_xpixels = (int)((window_xsize / window_ysize) * window_maxpixels);
		}
		viewportWidth = window_xpixels;
		textHeight = textHeightPx * height / window_ypixels;
		viewportHeight = window_ypixels + textHeightPx;
		window_ypixels = viewportHeight;
		height += textHeight;

		plot = new real[ni * nj];
		cudaMalloc(&plot_device, sizeof(real) * ni * nj);
		plot_rgba = new unsigned int[ni * nj];

		//
		// Iinitialise OpenGL display - use glut
		//
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);

		/*glutInitWindowSize(ni, nj);      // Window of ni x nj pixels*/
		glutInitWindowSize(window_xpixels, window_ypixels);
		glutInitWindowPosition(50, 50);  // position
		glutCreateWindow("MUSCL solver");       // title

		// Check for OpenGL extension support
		std::cout << "OpenGL: Loading extensions: " << glewGetErrorString(glewInit()) << std::endl;
		if (!glewIsSupported(
					"GL_VERSION_2_0 "
					"GL_ARB_pixel_buffer_object "
					"GL_EXT_framebuffer_object "
					)) {
			std::cerr << "ERROR: Support for necessary OpenGL extensions missing." << std::endl;
			return;
		}

		// Set up view
		glClearColor(0.0, 0.0, 0.0, 0.0);
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glOrtho(0, width, height, 0, -1, 1);
		glMatrixMode(GL_MODELVIEW);
		glDisable(GL_DEPTH_TEST);
		glShadeModel(GL_FLAT);

		// Create texture which we use to display the result and bind to gl_Tex
		glEnable(GL_TEXTURE_2D);
		glGenTextures(1, &gl_Tex);                     // Generate 2D texture
		glBindTexture(GL_TEXTURE_2D, gl_Tex);          // bind to gl_Tex
		// texture properties:
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, ni, nj, 0,
				GL_RGBA, GL_UNSIGNED_BYTE, NULL);


		// Create pixel buffer object and bind to gl_PBO. We store the data we want to
		// plot in memory on the graphics card - in a "pixel buffer". We can then
		// copy this to the texture defined above and send it to the screen
		glGenBuffers(1, &gl_PBO);
		glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, gl_PBO);
		glBufferData(GL_PIXEL_UNPACK_BUFFER_ARB, sizeof(float) * ni * nj, NULL, GL_STREAM_COPY);
		CUDA_SAFE_CALL(cudaGLRegisterBufferObject(gl_PBO));

		// Set the call-back functions and start the glut loop
		glutDisplayFunc(dispatchDisplay);
		glutReshapeFunc(dispatchResize);
		glutIdleFunc(dispatchIdle);

		glutKeyboardFunc(dispatchKeyboard);
		glutIgnoreKeyRepeat(1);
		glutSetCursor(cursor);

		glutMainLoop();
		return;
	}

	void draw(typename Mesh<GPU>::type& grid) {
		// For plotting, map the plot_rgba_data array to the
		// gl_PBO pixel buffer
		glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, gl_PBO);
		CUDA_SAFE_CALL(cudaGLMapBufferObject((void**)&plot_rgba_data, gl_PBO));

		if (autoscale) {
			lastLimits = renderData(grid, plot_device, plot_rgba_data, ni, nj, plotVariable, colourMode);
		} else {
			renderData(grid, plot_device, plot_rgba_data, ni, nj, plotVariable, colourMode, lastLimits.min, lastLimits.max);
		}

		cudaGLUnmapBufferObject(gl_PBO);
		glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, 0);
		glutPostRedisplay();
	}

	void dispatchDraw(typename Mesh<GPU>::type& grid) {
		gridToRender = &grid;
	}

	static void dispatchDisplay() {
		//std::cout << "dispatching display" << std::endl;
		OpenGLOutputter::getInstance().display();
	}
	static void dispatchIdle() {
		OpenGLOutputter::getInstance().idle();
	}
	static void dispatchResize(int w, int h) {
		OpenGLOutputter::getInstance().resize(w, h);
	}
	static void dispatchKeyboard(unsigned char key, int x, int y) {
		OpenGLOutputter::getInstance().keyboard(key, x, y);
	}
	void drawText(float x, float y, std::string text) {
		glPushMatrix();
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glColor4f(1.0, 1.0, 1.0, 1.0);
		glRasterPos2f(x, y * height / viewportHeight);
		//glScalef(0.075, 0.075, 1.0);
		//glColor3f(1.0, 1.0, 1.0);
		glutBitmapString(GLUT_BITMAP_HELVETICA_12, (const unsigned char*) text.c_str());
		glPopMatrix();
	}
	void drawBoundaryRegion(int variable, real threshold) {
		renderKernel<<<gridDim, blockDim>>>(grid, plot_device, ni, nj, variable);
		cudaThreadSynchronize();
		checkForError();

		std::vector<float> boundary = getMarchingSquares(plot_device, ni, nj, ni, nj, threshold);

		glBegin(GL_LINES);
		for (int i = 0; i < boundary.size(); i += 2) {
			glVertex2f(boundary[i], boundary[i + 1] + height - nj);
		}
		glEnd();
	}

	std::string getVariableName(int variable) {
		switch (variable) {
			case RHO: return "Density";
			case PRESSUREFOO: return "Pressure";
			case XVELOCITYPLOT: return "X Velocity";
			case YVELOCITYPLOT: return "Y Velocity";
			case VELOCITYPLOT: return "|Velocity|";
			case SOUNDSPEED: return "Sound speed";
			case SCHLIEREN: return "Schlieren";
			case VORTICITY: return "Vorticity";
			case POSVORTICITY: return "+ve vorticity";
			case NEGVORTICITY: return "-ve vorticity";
			case DENSITY + DISPLAYTYPES: return "Density";
			case XMOMENTUM + DISPLAYTYPES: return "X Momentum";
			case YMOMENTUM + DISPLAYTYPES: return "Y Momentum";
			case ENERGY + DISPLAYTYPES: return "Energy";
#ifdef GHOST
			case PHI + DISPLAYTYPES: return "Phi";
			case FRACTION + DISPLAYTYPES: return "Fraction";
#endif
			default: return "";
		}
	}

	void display() {
		//std::cout << "display" << std::endl;
		glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, gl_PBO);
		// Copy the pixel buffer to the texture, ready to display
		glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, ni, nj, GL_RGBA, GL_UNSIGNED_BYTE, 0);

		// Render one quad to the screen and colour it using our texture
		// i.e. plot our plotvar data to the screen
		glEnable(GL_TEXTURE_2D);
		glClear(GL_COLOR_BUFFER_BIT);
		glColor3f(1.0, 1.0, 1.0);
		glBegin(GL_QUADS);
		glTexCoord2f(0.0, 0.0);
		glVertex2f(0.0, height - nj);
		glTexCoord2f(1.0, 0.0);
		glVertex2f(width, height - nj);
		glTexCoord2f(1.0, 1.0);
		glVertex2f(width, height);
		glTexCoord2f(0.0, 1.0);
		glVertex2f(0.0, height);
		glEnd();
		glDisable(GL_TEXTURE_2D);
		glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, 0);

#ifdef GHOST
		if (showBoundaries) {
			glLineWidth(1.0);
			glColor3f(0.0, 0.0, 0.0);
Gx=0ρD1∂Y1∂x...ρDNS−1∂YNS−1∂x−τxx−τxy−uτxx−vτxy+qx000000Gy=0ρD1∂Y1∂y...ρDNS−1∂YNS−1∂y−τxy−τyy−uτxy−vτyy+qy000000S=Np
			drawBoundaryRegion(PHI + DISPLAYTYPES, 0.0);
		}
#endif

		std::stringstream plotting;
		plotting << "Plotting " << getVariableName(plotVariable) << ": " << std::setprecision(3) << lastLimits.min << " to " << lastLimits.max;
		drawText(10, 15, plotting.str());

		std::stringstream time;
		time << "Time: " << std::fixed << std::setprecision(8) << grid.time() << " s" << std::endl;
		drawText(ni / 2, 15, time.str());

		glFlush();
		glutSwapBuffers();
		static int frame = 0;
		//if (frame % 10 == 0) dumpBuffer(frame);
		//std::cout << "ufoo" << std::endl;
		frame++;
	}
	void keyboard(unsigned char key, int x, int y) {
//		std::cout << key << std::endl;
		if (key == '-') {
			if (--plotVariable < 0) plotVariable = 0;
		} else if (key == '+') {
			if (++plotVariable >= Cell::length() + DISPLAYTYPES) plotVariable = Cell::length() + DISPLAYTYPES - 1;
		} else if (key == '1') {
			colourMode = GREYSCALE;
		} else if (key == '2') {
			colourMode = HUE;
		} else if (key == '3') {
			colourMode = CUBEHELIX;
		} else if (key == ' ') {
			paused = !paused;
		} else if (key == 'b') {
			showBoundaries = !showBoundaries;
		} else if (key == 'a') {
			autoscale = !autoscale;
		}
	}

	void idle() {
//		boost::mutex::scoped_lock(MEIN);
		if (gridToRender != NULL) {
			draw(*gridToRender);
			if (!paused) {
				gridToRender = NULL;
			}
		}
	}
	void resize(int w, int h) {
		h -= textHeightPx;
		height = nj;
		if (w * height < h * width) {
			h = (float) height * w / width;
		} else {
			w = (float) width * h / height;
		}
		viewportWidth = w; viewportHeight = h + textHeightPx;
		textHeight = textHeightPx * (nj) / viewportHeight;
		height = nj + textHeight;
		glViewport(0, 0, viewportWidth, viewportHeight);
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glOrtho(0, width, height, 0, -1, 1);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
	}
	void dumpBuffer(int n) {
		glPixelStorei(GL_PACK_ALIGNMENT, 1);
		glReadBuffer(GL_BACK_LEFT);
		unsigned char* image = new unsigned char[3 * viewportWidth * viewportHeight];
		glReadPixels(0, 0, viewportWidth, viewportHeight, GL_RGB, GL_UNSIGNED_BYTE, image);
		std::stringstream filename;
		filename << "frames/frame" << std::setfill('0') << std::setw(8) << n << ".png";
		WriteBitmap(filename.str().c_str(), image, viewportWidth, viewportHeight, true);
		delete[] image;
	}
};

OpenGLOutputter* OpenGLOutputter::instance = NULL;


#endif

