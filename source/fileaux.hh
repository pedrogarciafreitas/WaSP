#ifndef FILEAUX_HH
#define FILEAUX_HH

#include <string>

#ifdef __unix__
#define _popen popen
#define _pclose pclose
#endif

int system_1(char *str);
int system_1(std::string sys_call_str);
long aux_GetFileSize(char* filename);
long aux_GetFileSize(const char* filename);

bool auxReadRAW(const char *filename, const int height, const int width, const int ncomp, unsigned short *img);

#endif