#pragma once
#include "ghost.hpp"
#include "png.h"
#include <iostream>
#include <vector>

enum DisplayTypes {
  RHO,
  PRESSUREFOO,
#ifdef REACTIVE
  PMAXPLOT,
  TEMPERATUREPLOT,
  YCH4PLOT,    // CH4 mass fraction
  YO2PLOT,     // O2 mass fraction
  YH2OPLOT,    // H2O mass fraction
  YCO2PLOT,    // CO2 mass fraction
  YOHPLOT,     // OH mass fraction
  YCOPLOT,     // CO mass fraction
  YN2PLOT,     // N2 mass fraction
#endif
  XVELOCITYPLOT,
  YVELOCITYPLOT,
  VELOCITYPLOT,
  SCHLIEREN,
  SOUNDSPEED,
  VORTICITY,
  POSVORTICITY,
  NEGVORTICITY,
  DISPLAYTYPES
};
enum ColourMode { GREYSCALE, HUE, CUBEHELIX, COLOURMODES };

__device__ unsigned int cubehelix(real v) {
  const real gamma = 1;
  const real rgb = 0.5;
  const real rots = -1.5;
  const real hue = 1.2;

  real dg = pow(v, gamma);
  real a = hue * dg * (1 - dg) / 2;
  real phi = 2 * M_PI * (rgb / 3 + rots * v);

  real cosphi = a * cos(phi);
  real sinphi = a * sin(phi);

  real r = dg - 0.14861 * cosphi + 1.78277 * sinphi;
  real g = dg - 0.29227 * cosphi - 0.90649 * sinphi;
  real b = dg + 1.97249 * cosphi;
  if (r < 0)
    r = 0;
  if (r > 1)
    r = 1;
  if (g < 0)
    g = 0;
  if (g > 1)
    g = 1;
  if (b < 0)
    b = 0;
  if (b > 1)
    b = 1;

  return (int)(r * 255) | ((int)(g * 255) << 8) | ((int)(b * 255) << 16);
}

__device__ unsigned int hue(real v) {
  const real S = 1, V = 1;
  const real C = V * S;
  const real H = 240 * (1 - v);
  const real Hp = H / 60.0;
  const real X = C * (1 - fabs(fmod(Hp, (real)2) - 1));
  real rp, gp, bp;
  if (0 <= Hp && Hp < 1) {
    rp = C;
    gp = X;
    bp = 0;
  } else if (1 <= Hp && Hp < 2) {
    rp = X;
    gp = C;
    bp = 0;
  } else if (2 <= Hp && Hp < 3) {
    rp = 0;
    gp = C;
    bp = X;
  } else if (3 <= Hp && Hp < 4) {
    rp = 0;
    gp = X;
    bp = C;
  } else if (4 <= Hp && Hp < 5) {
    rp = X;
    gp = 0;
    bp = C;
  } else if (5 <= Hp && Hp < 6) {
    rp = C;
    gp = 0;
    bp = X;
  } else {
    rp = 0;
    gp = 0;
    bp = 0;
  }
  const real m = V - C;
  const real r = rp + m, g = gp + m, b = bp + m;

  return (int)(r * 255) | ((int)(g * 255) << 8) | ((int)(b * 255) << 16);
}

__device__ unsigned int greyscale(real v) {
  unsigned int n = 255 * v;
  return n | (n << 8) | (n << 16);
}

__global__ void renderKernel(const typename Mesh<GPU>::type grid, real *data,
                             const int ni, const int nj, int plotVariable,
                             const real DV = 0) {
  const int overlap = 0;
  const int i = (blockDim.x - overlap) * blockIdx.x + threadIdx.x,
            j = (blockDim.y - overlap) * blockIdx.y + threadIdx.y;
  const int gi = i, gj = j;

  if (i < ni && j < nj) {
    real p[NUMBER_VARIABLES];
    conservativeToPrimitive(grid(gi, gj), p);

#ifdef GHOST
    if (plotVariable != PHI + DISPLAYTYPES &&
        plotVariable != FRACTION + DISPLAYTYPES && grid(gi, gj)[PHI] > 0) {
      data[ni * j + i] = NAN;
    } else
#endif
#ifdef REACTIVE
    /*for (int n_x = 1; n_x < (number-1); n_x++){
                   if (grid.x(gi) >= (300 + 0*space) && grid.x(gi) <= (300 +
    0*space + length_x))
                   {
                           data[ni * j + i] = NAN;
                           }
    }*/

    /*
     for (int n_x = 1; n_x < (number-1); n_x++){

            if (grid.x(gi) >= (300 + n_x*space) && grid.x(gi) <= (300 +
     n_x*space + length_x)){ data[ni * j + i] = NAN;
            }
            } else if (grid.x(gi) >= (300 + (number-1)*space) && grid.x(gi) <=
     (300 + (number-1)*space + length_x))
                    {
                            data[ni * j + i] = NAN;
                            } else */
    /*
        if ((grid.x(gi) >= X31 && grid.x(gi) <= X32 && grid.y(gj) >= Y31 &&
     grid.y(gj) <= Y32) || (grid.x(gi) >= X41 && grid.x(gi) <= X42 && grid.y(gj)
     >= Y41 && grid.y(gj) <= Y42) || (grid.x(gi) >= X51 && grid.x(gi) <= X52 &&
     grid.y(gj) >= Y51 && grid.y(gj) <= Y52)|| (grid.x(gi) >= X61 && grid.x(gi)
     <= X62 && grid.y(gj) >= Y61 && grid.y(gj) <= Y62)|| (grid.x(gi) >= X71 &&
     grid.x(gi) <= X72 && grid.y(gj) >= Y71 && grid.y(gj) <= Y72)|| (grid.x(gi)
     >= X81 && grid.x(gi) <= X82 && grid.y(gj) >= Y81 && grid.y(gj) <= Y82)||
     (grid.x(gi) >= X91 && grid.x(gi) <= X92 && grid.y(gj) >= Y91 && grid.y(gj)
     <= Y92)|| (grid.x(gi) >= X101 && grid.x(gi) <= X102 && grid.y(gj) >= Y101
     && grid.y(gj) <= Y102)|| (grid.x(gi) >= X103 && grid.x(gi) <= X104 &&
     grid.y(gj) >= Y103 && grid.y(gj) <= Y104) || (grid.x(gi) >= X105 &&
     grid.x(gi) <= X106 && grid.y(gj) >= Y105 && grid.y(gj) <= Y106)||
     (grid.x(gi) >= X107 && grid.x(gi) <= X108 && grid.y(gj) >= Y107 &&
     grid.y(gj) <= Y108) || (grid.x(gi) >= X109 && grid.x(gi) <= X110 &&
     grid.y(gj) >= Y109 && grid.y(gj) <= Y110) || (grid.x(gi) >= X111 &&
     grid.x(gi) <= X112 && grid.y(gj) >= Y111 && grid.y(gj) <= Y112) ||
     (grid.x(gi) >= X113 && grid.x(gi) <= X114 && grid.y(gj) >= Y113	&&
     grid.y(gj) <= Y114) || (grid.x(gi) >= X115 && grid.x(gi) <= X116 &&
     grid.y(gj) >= Y115 && grid.y(gj) <= Y116) || (grid.x(gi) >= X117 &&
     grid.x(gi) <= X118 && grid.y(gj) >= Y117 && grid.y(gj) <= Y118) ||
     (grid.x(gi) >= X119 && grid.x(gi) <= X120 && grid.y(gj) >= Y119 &&
     grid.y(gj) <= Y120) || (grid.x(gi) >= X121 && grid.x(gi) <= X122 &&
     grid.y(gj) >= Y121 && grid.y(gj) <= Y122) ) {
     //if ((grid.x(gi) >= X11 && grid.x(gi) <= X12 && grid.y(gj) >= Y11 &&
     grid.y(gj) <= Y12) || (grid.x(gi) >= X21 && grid.x(gi) <= X22 && grid.y(gj)
     >= Y21 && grid.y(gj) <= Y22) || (grid.x(gi) >= X31 && grid.x(gi) <= X32 &&
     grid.y(gj) >= Y31 && grid.y(gj) <= Y32) || (grid.x(gi) >= X41 && grid.x(gi)
     <= X42 && grid.y(gj) >= Y41 && grid.y(gj) <= Y42) || (grid.x(gi) >= X51 &&
     grid.x(gi) <= X52 && grid.y(gj) >= Y51 && grid.y(gj) <= Y52)|| (grid.x(gi)
     >= X61 && grid.x(gi) <= X62 && grid.y(gj) >= Y61 && grid.y(gj) <= Y62)) {
           //  if ((grid.x(gi) <= XCORNER && grid.y(gj) >= YCORNER)) {
          //   if (grid.x(gi) <= XCORNER && grid.y(gj) >= YCORNER) {
              data[ni * j + i] = NAN;
        //      } else if ( grid.x(gi) >= 1.2*XCORNER && grid.y(gj) >= YCORNER){
        //    data[ni * j + i] = NAN;
         } else
    */

#endif

        if (plotVariable >= DISPLAYTYPES) {
      data[ni * j + i] = grid(gi, gj)[plotVariable - DISPLAYTYPES];
    } else if (plotVariable == RHO) {
      data[ni * j + i] = grid(gi, gj)[DENSITY];
    } else if (plotVariable == PRESSUREFOO) {
      data[ni * j + i] = p[PRESSURE];
    } else if (plotVariable == PMAXPLOT) {
      data[ni * j + i] = grid(gi, gj)[PMAX];
    } else if (plotVariable == TEMPERATUREPLOT) {
      Cell c = grid(gi, gj);
      real Y[NUM_SPECIES];
      getMassFractions(c, Y);
      real rho = c[DENSITY];
      real e_kin = 0.5 * (c[XVELOCITY]*c[XVELOCITY] + c[YVELOCITY]*c[YVELOCITY]);
      real e_int = (c[ENERGY]/(rho+1e-30)) - e_kin;
      real R_mix = mixtureRgas(Y);
      data[ni * j + i] = (GAMMA_DEFAULT - 1.0) * e_int / (R_mix + 1e-30);
    } else if (plotVariable == YCH4PLOT) {
      data[ni * j + i] = p[SPECIES_INDEX(12)];  // CH4 index 12
    } else if (plotVariable == YO2PLOT) {
      data[ni * j + i] = p[SPECIES_INDEX(3)];   // O2 index 3
    } else if (plotVariable == YH2OPLOT) {
      data[ni * j + i] = p[SPECIES_INDEX(6)];   // H2O index 6
    } else if (plotVariable == YCO2PLOT) {
      data[ni * j + i] = p[SPECIES_INDEX(9)];   // CO2 index 9
    } else if (plotVariable == YOHPLOT) {
      data[ni * j + i] = p[SPECIES_INDEX(4)];   // OH index 4
    } else if (plotVariable == YCOPLOT) {
      data[ni * j + i] = p[SPECIES_INDEX(8)];   // CO index 8
    } else if (plotVariable == YN2PLOT) {
      data[ni * j + i] = p[SPECIES_INDEX(NUM_SPECIES - 1)];  // N2 (last species)
    } else if (plotVariable == XVELOCITYPLOT) {
      data[ni * j + i] = grid(gi, gj)[XMOMENTUM] / grid(gi, gj)[DENSITY];
    } else if (plotVariable == YVELOCITYPLOT) {
      data[ni * j + i] = grid(gi, gj)[YMOMENTUM] / grid(gi, gj)[DENSITY];
    } else if (plotVariable == VELOCITYPLOT) {
      data[ni * j + i] = sqrt(pow(grid(gi, gj)[XMOMENTUM]/grid(gi, gj)[DENSITY], 2) +
                              pow(grid(gi, gj)[YMOMENTUM]/grid(gi, gj)[DENSITY], 2));
    } else if (plotVariable == SCHLIEREN) {
      const real drhodx =
          (grid(gi + 1, gj)[DENSITY] - grid(gi - 1, gj)[DENSITY]) /
          (2.0 * grid.dx());
      const real drhody =
          (grid(gi, gj + 1)[DENSITY] - grid(gi, gj - 1)[DENSITY]) /
          (2.0 * grid.dy());

      // data[ni * j + i] = exp(-0.0003 * sqrt(drhodx * drhodx + drhody *
      // drhody));  // 无量纲原始值
      data[ni * j + i] = exp(-0.005 * sqrt(drhodx * drhodx + drhody * drhody));
    } else if (plotVariable == SOUNDSPEED) {

      // temperature
      data[ni * j + i] = p[PRESSURE] / p[DENSITY];

      // data[ni * j + i] = soundSpeed(grid(gi, gj));

    } else if (plotVariable == VORTICITY || plotVariable == POSVORTICITY ||
               plotVariable == NEGVORTICITY) {
      if (i > 0 && j > 0 && i < ni - 1 && j < nj - 1) {
        const real uPlus =
            grid(gi, gj + 1)[XMOMENTUM] / grid(i, j + 1)[DENSITY];
        const real uMinus =
            grid(gi, gj - 1)[XMOMENTUM] / grid(i, j - 1)[DENSITY];

        const real vPlus =
            grid(gi + 1, gj)[YMOMENTUM] / grid(gi + 1, gj)[DENSITY];
        const real vMinus =
            grid(gi - 1, gj)[YMOMENTUM] / grid(gi - 1, gj)[DENSITY];

        // second order approximation to du/dz and dv/dx
        data[ni * j + i] = (uPlus - uMinus) / (2 * grid.dy()) -
                           (vPlus - vMinus) / (2 * grid.dx());
        if (plotVariable == POSVORTICITY) {
          data[ni * j + i] = fmax(0, data[ni * j + i]);
        } else if (plotVariable == NEGVORTICITY) {
          data[ni * j + i] = fmin(0, data[ni * j + i]);
        }
      } else {
        data[ni * j + i] = 0;
      }
    }
  }

  if (Check3 > 0){
    if (i < ni && j < nj) {
    real p[NUMBER_VARIABLES];
    conservativeToPrimitive(grid(gi, gj), p);
			
	if (pow(grid.x(gi) - (l_bump * cos(theta_C * PI / 180) + A_x),2) + pow(grid.y(gj) - (AB_y - sin(theta_C * PI / 180) * l_bump),2) >= pow(r_bump,2)){
		//if (pow(grid.x(gi) - 40,2) + pow(grid.y(gj) - 40,2) < pow(r_bump,2)){
			data[ni * j + i] = NAN;
			}
		
	}
  }
}

struct limits {
  real min, max;
};

template <int blockSize>
__device__ void getLimitsKernelLastWarp(volatile limits *data, int tid) {
#pragma unroll
  for (int k = 32; k >= 1; k >>= 1) {
    if (blockSize >= k * 2) {
      data[tid].max = max(data[tid].max, data[tid + k].max);
      data[tid].min = min(data[tid].min, data[tid + k].min);
    }
  }
}

template <int blockSize>
__global__ void getLimitsKernel(const real *plot, const int ni, const int nj,
                                const int width, const int height,
                                Grid2D<limits, GPU> output) {
  // tid is offset to run from 0...blockSize in 1D
  const int tid = blockDim.x * threadIdx.y + threadIdx.x;
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  const int j = blockIdx.y * blockDim.y + threadIdx.y;

  __shared__ limits data[blockSize];

  data[tid].min = 1e30;
  data[tid].max = -1e30;

  if (i < width && j < height && i < ni && j < nj) {
    data[tid].min = plot[j * ni + i];
    data[tid].max = plot[j * ni + i];
  }
  __syncthreads();

#pragma unroll
  for (int k = blockSize; k > 32; k >>= 1) {
    if (blockSize >= k * 2) {
      if (tid < k) {
        data[tid].max = max(data[tid].max, data[tid + k].max);
        data[tid].min = min(data[tid].min, data[tid + k].min);
        __syncthreads();
      }
    }
  }
  if (tid < 32)
    getLimitsKernelLastWarp<blockSize>(data, tid);

  output(blockIdx.x, blockIdx.y, 0).min = data[0].min;
  output(blockIdx.x, blockIdx.y, 0).max = data[0].max;
}

template <typename T>
__global__ void normaliseKernel(const real *const data, const int ni,
                                const int nj, const limits l, T *const plot,
                                ColourMode colour) {
  const int i = blockDim.x * blockIdx.x + threadIdx.x;
  const int j = blockDim.y * blockIdx.y + threadIdx.y;

  unsigned int value;

  if (i < ni && j < nj) {
    if (isnan(data[j * ni + i])) {
      value = 0x808080;
    } else {
      real frac = (data[j * ni + i] - l.min) / (l.max - l.min);
      if (frac < 0)
        frac = 0;
      if (frac > 0.999)
        frac = 0.999;
      if (colour == GREYSCALE) {
        value = greyscale(frac);
      } else if (colour == HUE) {
        value = hue(frac);
      } else if (colour == CUBEHELIX) {
        value = cubehelix(frac);
      }
    }
    if (sizeof(T) == 1) {
      int offset = 3 * (j * ni + i);
      plot[offset + 0] = (value >> 0) & 0xFF;
      plot[offset + 1] = (value >> 8) & 0xFF;
      plot[offset + 2] = (value >> 16) & 0xFF;
    } else {
      plot[j * ni + i] = value;
    }
  }
}

template <typename T>
limits renderData(typename Mesh<GPU>::type &grid, real *plot_data, T *plot_rgba,
                  int ni, int nj, int plotVariable, ColourMode colourMode,
                  real lMin = 0, real lMax = 0) {
  dim3 blockDim = dim3(16, 16, 1);
  dim3 gridDim = dim3((ni + 15) / 16, (nj + 15) / 16, 1);

  renderKernel<<<gridDim, blockDim>>>(grid, plot_data, ni, nj, plotVariable);
  cudaThreadSynchronize();
  checkForError();

  limits l;
  if (lMin == 0 && lMax == 0) {
    Grid2D<limits, GPU> outputGPU(gridDim.x, gridDim.y);
    getLimitsKernel<256>
        <<<gridDim, blockDim>>>(plot_data, ni, nj, ni, nj, outputGPU);
    cudaThreadSynchronize();
    checkForError();

    Grid2D<limits, CPU> output(outputGPU);
    outputGPU.free();

    l.min = 1e100;
    l.max = -1e100;
    for (int i = 0; i < output.Nx(); i++) {
      for (int j = 0; j < output.Ny(); j++) {
        l.min = min(l.min, output(i, j, 0).min);
        l.max = max(l.max, output(i, j, 0).max);
      }
    }
    output.free();
  } else {
    l.min = lMin;
    l.max = lMax;
  }

  normaliseKernel<<<gridDim, blockDim>>>(plot_data, ni, nj, l, plot_rgba,
                                         colourMode);
  cudaThreadSynchronize();
  checkForError();

  return l;
}

int WriteBitmap(const char *name, unsigned char *buff, unsigned width,
                unsigned height, bool flip) {
  FILE *fp = fopen(name, "wb");
  if (!fp) {
    return -1;
  }

  png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
                                                (png_voidp)NULL, NULL, NULL);
  if (!png_ptr) {
    return -1;
  }

  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr) {
    png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
    return -1;
  }

  png_init_io(png_ptr, fp);

  png_set_IHDR(png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_RGB,
               PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
               PNG_FILTER_TYPE_DEFAULT);

  png_write_info(png_ptr, info_ptr);

  png_byte *image = (png_byte *)buff;

  unsigned k;
  png_bytep *row_pointers = new png_bytep[height];
  for (k = 0; k < height; k++) {
    row_pointers[k] = image + (flip ? (height - k - 1) : k) * width * 3;
  }

  png_write_image(png_ptr, row_pointers);
  png_write_end(png_ptr, info_ptr);
  png_destroy_write_struct(&png_ptr, &info_ptr);
  delete[] row_pointers;
  fclose(fp);
  return 0;
}

limits saveFrame(typename Mesh<GPU>::type &grid, int plotVariable,
                 ColourMode colourMode, const char *filename, real lMin = 0.0,
                 real lMax = 0.0, real *temp = 0) {
  int ni = grid.activeNx(), nj = grid.activeNy();
  real *plot_data;
  unsigned char *plot_rgba;

  if (temp == 0) {
    cudaMalloc(&plot_data, ni * nj * sizeof(real));
    cudaMalloc(&plot_rgba, ni * nj * sizeof(unsigned char) * 3);
  } else {
    plot_data = temp;
    plot_rgba = (unsigned char *)(temp + ni * nj);
  }

  limits l = renderData(grid, plot_data, plot_rgba, ni, nj, plotVariable,
                        colourMode, lMin, lMax);

  unsigned char *plot_rgba_host = new unsigned char[ni * nj * 3];
  cudaMemcpy(plot_rgba_host, plot_rgba, ni * nj * 3 * sizeof(unsigned char),
             cudaMemcpyDeviceToHost);
  WriteBitmap(filename, plot_rgba_host, ni, nj, false);

  if (temp == 0) {
    cudaFree(plot_data);
    cudaFree(plot_rgba);
  }
  delete[] plot_rgba_host;

  return l;
}
