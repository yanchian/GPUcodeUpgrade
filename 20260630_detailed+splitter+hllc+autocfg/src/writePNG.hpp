/* Utility routines for saving plots of the current simulation to PNG bitmap */
/* files */

#include "png.h"

void getColor(double value, double minimum, double maximum, float* colour) {
	float normalised = (value - minimum) / (maximum - minimum);
	if (normalised > 1.) {
		normalised = 1.;
	}
	if (normalised < 0.) {
		normalised = 0.;
	}
	for (int i = 0; i < 3; i++) {
		colour[i] = normalised;
	}
}

int WriteBitmap(const char* name, unsigned char* buff, unsigned width, unsigned height) {
	FILE* fp = fopen(name, "wb");
	if (!fp) {
		return -1;
	}

	png_structp png_ptr = png_create_write_struct
	                      (PNG_LIBPNG_VER_STRING, (png_voidp)NULL, NULL, NULL);
	if (!png_ptr) {
		return -1;
	}

	png_infop info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr) {
		png_destroy_write_struct(&png_ptr,
		                         (png_infopp)NULL);
		return -1;
	}

	png_init_io(png_ptr, fp);

	png_set_IHDR(png_ptr, info_ptr, width, height, 8,
	             PNG_COLOR_TYPE_RGB,
	             PNG_INTERLACE_NONE,
	             PNG_COMPRESSION_TYPE_DEFAULT,
	             PNG_FILTER_TYPE_DEFAULT);

	png_write_info(png_ptr, info_ptr);

	png_byte* image = (png_byte*)buff;

	unsigned k;
	png_bytep* row_pointers = new png_bytep[height];
	for (k = 0; k < height; k++) {
		row_pointers[k] = image + (k) * width * 3;
	}

	png_write_image(png_ptr, row_pointers);
	png_write_end(png_ptr, info_ptr);
	png_destroy_write_struct(&png_ptr, &info_ptr);
	delete[] row_pointers;
	fclose(fp);
	return 0;
}


