#include <cstdio>
#include <cmath>

#include "configLoader.hh"

global_parameters* load_global_parameters(const char *config_file) {

	global_parameters *global_params = new global_parameters;

	FILE *filept;

	filept = fopen(config_file, "rb");

	int n_views_total, yuv_transform_s, yuv_ratio_search_s, std_search_s;

	fread(&n_views_total, sizeof(int32_t), 1, filept); /*reading*/
	fread(&yuv_transform_s, sizeof(int32_t), 1, filept); /*reading*/
	fread(&yuv_ratio_search_s, sizeof(int32_t), 1, filept); /*reading*/
	fread(&std_search_s, sizeof(int32_t), 1, filept); /*reading*/

	global_params->n_views_total = n_views_total;
	global_params->YUV_TRANSFORM = yuv_transform_s > 0 ? true : false;
	global_params->YUV_RATIO_SEARCH = yuv_ratio_search_s > 0 ? true : false;
	global_params->STD_SEARCH = std_search_s > 0 ? true : false;

	fclose(filept);

	return global_params;
}

view* load_config_and_init_LF(const char *config_file) {

	global_parameters *global_params = load_global_parameters(config_file);

	view *LF = new view[global_params->n_views_total]();

	//unsigned short MINIMUM_DEPTH = 0;

	FILE *filept;

	filept = fopen(config_file, "rb");

	fseek(filept, 4 * sizeof(int32_t), SEEK_SET);

	for (int ii = 0; ii < global_params->n_views_total; ii++) {

		view *SAI = LF + ii;

		initView(SAI);

		SAI->i_order = ii;

		fread(&(SAI->r), sizeof(int), 1, filept); /*reading*/
		fread(&(SAI->c), sizeof(int), 1, filept); /*reading*/
												  /* find number of rows and cols */

		//maxR = (SAI->r + 1 > maxR) ? SAI->r + 1 : maxR;
		//maxC = (SAI->c + 1 > maxC) ? SAI->c + 1 : maxC;

		int xx = 0, yy = 0;

		fread(&xx, sizeof(int), 1, filept); /*reading*/
		fread(&yy, sizeof(int), 1, filept); /*reading*/

		if (abs(xx) > 0) {
			SAI->has_x_displacement = true;
			SAI->x = static_cast<float>(xx) / 100000;
		}

		if (abs(yy) > 0) {
			SAI->has_y_displacement = true;
			SAI->y = -static_cast<float>(yy) / 100000;
		}

		int rate_color, rate_depth;

		fread(&rate_color, sizeof(int), 1, filept); /*reading*/
		fread(&rate_depth, sizeof(int), 1, filept); /*reading*/

		SAI->residual_rate_color = static_cast<float>(rate_color) / 100000.0f;
		SAI->residual_rate_depth = static_cast<float>(rate_depth) / 100000.0f;

		fread(&SAI->Ms, sizeof(int), 1, filept); /*reading*/
		fread(&SAI->NNt, sizeof(int), 1, filept); /*reading*/

		int sigma = 0;

		fread(&sigma, sizeof(int), 1, filept); /*reading*/
		SAI->sigma = static_cast<float>(sigma) / 100000.0f;

		int tmpminv;

		fread(&tmpminv, sizeof(int), 1, filept); /*reading, if we have negative inverse depth,
												 for example in lenslet, we need to subtract min_inv_d
												 from the inverse depth maps*/
		//if (ii == 0) {
		//	if (tmpminv > 0) {
		//		MINIMUM_DEPTH = tmpminv;
		//	}
		//}

		//if (MINIMUM_DEPTH > 0) {
		//	//SAI->has_min_inv_depth = true;
		//	SAI->min_inv_d = (int)MINIMUM_DEPTH;
		//}

		if (ii > 0) {
			SAI->min_inv_d = LF->min_inv_d;
		}
		else {
			SAI->min_inv_d = tmpminv > 0 ? tmpminv : 0;
		}

		fread(&(SAI->n_references), sizeof(int), 1, filept); /*reading*/

		if (SAI->n_references > 0) {

			SAI->has_color_references = true;

			SAI->references = new int[SAI->n_references]();

			fread(SAI->references, sizeof(int), SAI->n_references, filept); /*reading*/

		}

		if (SAI->n_references > 0) {
			SAI->color_reference_views = new view*[SAI->n_references]();
			for (int ik = 0; ik < SAI->n_references; ik++) {
				SAI->color_reference_views[ik] = LF + SAI->references[ik];
			}
		}

		fread(&(SAI->n_depth_references), sizeof(int), 1, filept); /*reading*/

		if (SAI->n_depth_references > 0) {

			SAI->has_depth_references = true;

			SAI->depth_references = new int[SAI->n_depth_references]();

			fread(SAI->depth_references, sizeof(int), SAI->n_depth_references, filept); /*reading*/

		}

		if (SAI->n_depth_references > 0) {
			SAI->depth_reference_views = new view*[SAI->n_depth_references]();
			for (int ik = 0; ik < SAI->n_depth_references; ik++) {
				SAI->depth_reference_views[ik] = LF + SAI->depth_references[ik];
			}
		}

		fread(&SAI->has_segmentation, sizeof(int), 1, filept); // new,13.06.18 /*reading*/

		SAI->yuv_ratio_search = global_params->YUV_RATIO_SEARCH;
		SAI->yuv_transform = global_params->YUV_TRANSFORM;

	}

	fclose(filept);

	return LF;

}