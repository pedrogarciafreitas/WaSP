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
