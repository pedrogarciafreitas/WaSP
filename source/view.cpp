#include "view.hh"

#define nullptr 0

void initView(view* view)
{

	view->color = nullptr;
	view->depth = nullptr;
	view->segmentation = nullptr;

	view->r = 0;
	view->c = 0;

	view->min_inv_d = 0;
	view->n_references = 0;

	view->references = nullptr;
	view->depth_references = nullptr;

	view->merge_weights = nullptr;
	view->sparse_weights = nullptr;
	view->sparse_mask = nullptr;
	view->merge_weights_float = nullptr;
	view->number_of_pixels_per_region = nullptr;

	view->bmask = nullptr;
	view->seg_vp = nullptr;

	view->NNt = 0;
	view->Ms = 0;

	view->has_segmentation = 0;
	view->maxL = 0;

	view->DM_ROW = nullptr;
	view->DM_COL = nullptr;

	view->i_order = 0;

	view->stdd = 0.0;

	view->use_median = false;

	view->yuv_transform = true;

	view->has_color_residual = false;
	view->has_depth_residual = false;
	view->use_global_sparse = false;

	view->has_color_references = false;
	view->has_depth_references = false;
	//view->has_min_inv_depth = false;

	view->has_x_displacement = false;
	view->has_y_displacement = false;

	view->has_chrominance = false;

	view->depth_file_exist = false;

	view->final_psnr = 0.0;
	view->warp_psnr = 0.0;
	view->merge_psnr = 0.0;
	view->sparse_psnr = 0.0;

}