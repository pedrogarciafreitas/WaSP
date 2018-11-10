#include "view.hh"
#include "ppm.hh"
#include <cstdio>

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
	
	view->ncomp = 0;

}

bool loadColor(view* SAI) {

	unloadColor(SAI);

	return aux_read16PGMPPM(SAI->path_out_ppm, SAI->nc, SAI->nr, SAI->ncomp, SAI->color);
}

bool loadInverseDepth(view* SAI) {

	unloadInverseDepth(SAI);

	return aux_read16PGMPPM(SAI->path_out_pgm, SAI->nc, SAI->nr, SAI->ncomp, SAI->depth);
}

bool loadLabels(view* SAI) {

	unloadLabels(SAI);

	SAI->label_im = new int32_t[SAI->nr*SAI->nc]();

	FILE *tmpfile_im_labels = fopen(SAI->path_label_im, "rb");
	size_t nread = fread(SAI->label_im, sizeof(int32_t), SAI->nr*SAI->nc, tmpfile_im_labels);
	fclose(tmpfile_im_labels);

	if (nread < SAI->nr*SAI->nc * 4) {
		return false;
	}
	else {
		return true;
	}

}

void unloadLabels(view *SAI) {
	if (SAI->label_im != nullptr) {
		delete[](SAI->label_im);
		SAI->label_im = nullptr;
	}
}

void unloadColor(view* SAI) {
	if (SAI->color != nullptr) {
		delete[](SAI->color);
		SAI->color = nullptr;
	}
}

void unloadInverseDepth(view* SAI) {
	if (SAI->depth != nullptr) {
		delete[](SAI->depth);
		SAI->depth = nullptr;
	}
}

int32_t *loadWarpedLabelIm(view *SAI, view *ref_view) {

	int32_t *warpedLabelIm = new int32_t[SAI->nr*SAI->nc]();

	FILE *tmpfile_im_warped_labels;
	
	char warped_label_path[1024];
	sprintf(warped_label_path, "%s%c%03d_%03d%s%03d_%03d%s", SAI->output_dir, '/', SAI->c, SAI->r, "_im_labels_warped_to_",ref_view->c,ref_view->r,".int32");

	tmpfile_im_warped_labels = fopen(warped_label_path, "rb");
	fread(warpedLabelIm, sizeof(int), SAI->nr*SAI->nc, tmpfile_im_warped_labels);
	fclose(tmpfile_im_warped_labels);

}

bool writeWarpedLabelIm(view *SAI, view *ref_view, const int32_t *warpedLabelIm) {

	FILE *tmpfile_im_warped_labels;

	char warped_label_path[1024];
	sprintf(warped_label_path, "%s%c%03d_%03d%s%03d_%03d%s", SAI->output_dir, '/', ref_view->c, ref_view->r, "_im_labels_warped_to_", SAI->c, SAI->r, ".int32");

	tmpfile_im_warped_labels = fopen(warped_label_path, "wb");
	size_t nwritten = fwrite(warpedLabelIm, sizeof(int), SAI->nr*SAI->nc, tmpfile_im_warped_labels);
	fclose(tmpfile_im_warped_labels);

	if (nwritten < SAI->nr*SAI->nc * 4) {
		return false;
	}
	else {
		return true;
	}


}