#ifndef CONFIG_LOADER_HH
#define CONFIG_LOADER_HH

#include "view.hh"

struct global_parameters {

	int32_t n_views_total;
	bool YUV_TRANSFORM;
	bool YUV_RATIO_SEARCH;
	bool STD_SEARCH;

};

global_parameters* load_global_parameters(const char *config_file);

view* load_config_and_init_LF(const char *config_file);

#endif // !CONFIG_LOADER_HH
