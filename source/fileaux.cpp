#include "fileaux.hh"
#include <sys/stat.h>

#include <string>

#define SYSTEM_VERBOSE_QUIET true

int system_1(char *str) {

	std::string sys_call_str(str);

#ifdef SYSTEM_VERBOSE_QUIET

#ifdef _WIN32
	sys_call_str.append(" > nul");
#endif
#ifdef __unix__
	sys_call_str.append(" > /dev/null");
#endif

#endif

	return system(sys_call_str.c_str());

}

int system_1(std::string sys_call_str) {

#ifdef SYSTEM_VERBOSE_QUIET

#ifdef _WIN32
	sys_call_str.append(" > nul");
#endif
#ifdef __unix__
	sys_call_str.append(" > /dev/null");
#endif

#endif

	return system(sys_call_str.c_str());

}

long aux_GetFileSize(char* filename)
{
	struct stat stat_buf;
	int rc = stat(filename, &stat_buf);
	return rc == 0 ? stat_buf.st_size : -1;
}

long aux_GetFileSize(const char* filename)
{
	struct stat stat_buf;
	int rc = stat(filename, &stat_buf);
	return rc == 0 ? stat_buf.st_size : -1;
}

bool auxReadRAW(const char *filename, const int height, const int width, const int ncomp, unsigned short *img) {

	FILE *filept = fopen(filename, "rb");

	if (filept == nullptr) {
		printf("%s does not exist\n", filename);
		return false;
	}

	//img = new unsigned short[width*height * ncomp]();

	unsigned short *Image16bit = new unsigned short[width*height * ncomp]();

	/*--< Read 16bit ppm image from filept >--*/
	int nread = static_cast<int>(fread(Image16bit, sizeof(unsigned short), width*height * ncomp, filept));

	if (nread != width*height * ncomp)
	{
		fprintf(stderr, "READ ERROR aux_read16ppm() %s\n", filename);
		delete[](img);
		delete[](Image16bit);
		return false;
	}

	fclose(filept);

	int i = 0;

	for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; y++) {

			int red, green, blue;

			red = Image16bit[(x + y*width) * ncomp];
			if (ncomp == 3) {
				green = Image16bit[(x + y*width) * ncomp + 1];
				blue = Image16bit[(x + y*width) * ncomp + 2];
			}

			// Exhange upper 8bit and lower 8bit for Intel x86
			red = ((red & 0x00ff) << 8) | ((red & 0xff00) >> 8);
			if (ncomp == 3) {
				green = ((green & 0x00ff) << 8) | ((green & 0xff00) >> 8);
				blue = ((blue & 0x00ff) << 8) | ((blue & 0xff00) >> 8);
			}

			img[i] = red;
			if (ncomp == 3) {
				img[i + height*width] = green;
				img[i + 2 * height*width] = blue;
			}

			i++;

		}
	}

	delete[](Image16bit);

	return true;

}
