#include "view.hh"
#include "ppm.hh"
#include <cstdio>

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

	view->sigma = 0.0;

	view->use_median = false;

	view->yuv_transform = false;

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

	view->yuv_ratio_search = false;

	view->original_color_view = nullptr;

	view->warped_color_views = nullptr;
	view->warped_depth_views = nullptr;
	view->occlusion_masks = nullptr;

	view->color_reference_views = nullptr; // parents
	view->depth_reference_views = nullptr;

}

void setViewFilePaths(view* SAI, const char *output_dir, const char *input_dir) {

	sprintf(SAI->output_dir, "%s", output_dir);
	sprintf(SAI->input_dir, "%s", input_dir);

	sprintf(SAI->path_input_ppm, "%s%c%03d_%03d%s", SAI->input_dir, '/', SAI->c, SAI->r, ".ppm");
	sprintf(SAI->path_input_pgm, "%s%c%03d_%03d%s", SAI->input_dir, '/', SAI->c, SAI->r, ".pgm");
	sprintf(SAI->path_input_seg, "%s%c%03d_%03d%s", SAI->input_dir, '/', SAI->c, SAI->r, "_segmentation.pgm");

	sprintf(SAI->path_out_ppm, "%s%c%03d_%03d%s", SAI->output_dir, '/', SAI->c, SAI->r, ".ppm");
	sprintf(SAI->path_out_pgm, "%s%c%03d_%03d%s", SAI->output_dir, '/', SAI->c, SAI->r, ".pgm");
	sprintf(SAI->path_label_im, "%s%c%03d_%03d%s", SAI->output_dir, '/', SAI->c, SAI->r, "_im_labels.int32");
}

void initializeWarpingArrays(view* SAI) {

	/* holds partial warped views for ii */
	SAI->warped_color_views = new unsigned short*[SAI->n_references]();
	SAI->warped_depth_views = new unsigned short*[SAI->n_references]();
	SAI->occlusion_masks = new float*[SAI->n_references]();

	for (int ij = 0; ij < SAI->n_references; ij++) {
		SAI->warped_color_views[ij] = new unsigned short[SAI->nr*SAI->nc * 3]();
		SAI->warped_depth_views[ij] = new unsigned short[SAI->nr*SAI->nc]();
		SAI->occlusion_masks[ij] = new float[SAI->nr*SAI->nc]();
	}

}

void deinitializeWarpingArrays(view* SAI) {
	for (int ij = 0; ij < SAI->n_references; ij++)
	{
		delete[](SAI->warped_color_views[ij]);
		delete[](SAI->warped_depth_views[ij]);
		delete[](SAI->occlusion_masks[ij]);
	}

	delete[](SAI->warped_color_views);
	delete[](SAI->warped_depth_views);
	delete[](SAI->occlusion_masks);
}

bool loadOriginalColor(view* SAI) {
	return aux_read16PGMPPM(SAI->path_input_ppm, SAI->nc, SAI->nr, SAI->ncomp, SAI->original_color_view);
}

void unloadOriginalColor(view *SAI) {
	if (SAI->original_color_view != nullptr) {
		delete[](SAI->original_color_view);
		SAI->original_color_view = nullptr;
	}
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

int getNB(view *SAI) {
	/* returns the number of elements in the binary occlusion mask matrix */
	return (1 << SAI->n_references)*SAI->n_references;
}

void setBMask(view *SAI)
{

	/* sets the binary mask used to derive view availability in each of the MMM classes,
	size of the binary mask is [MMM x n_references] */

	int MMM = 1 << SAI->n_references;

	bool *bmask = new bool[MMM * SAI->n_references]();

	SAI->bmask = bmask;

	for (int ij = 0; ij < MMM; ij++) {

		int uu = ij;

		for (int ik = SAI->n_references - 1; ik >= 0; ik--) {

			if ( uu / (1 << ik) > 0)
			{
				uu = uu - (1 << ik);
				bmask[ij + ik * MMM] = 1;
			}

		}
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

	return warpedLabelIm;

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

void cleanView(view *SAI) {

	if (SAI->color != nullptr)
		delete[](SAI->color);
	if (SAI->depth != nullptr)
		delete[](SAI->depth);
	if (SAI->references != nullptr)
		delete[](SAI->references);
	if (SAI->depth_references != nullptr)
		delete[](SAI->depth_references);
	if (SAI->merge_weights != nullptr)
		delete[](SAI->merge_weights);
	if (SAI->sparse_weights != nullptr)
		delete[](SAI->sparse_weights);
	if (SAI->bmask != nullptr)
		delete[](SAI->bmask);
	if (SAI->seg_vp != nullptr)
		delete[](SAI->seg_vp);
	if (SAI->sparse_mask != nullptr)
		delete[](SAI->sparse_mask);
	if (SAI->segmentation != nullptr)
		delete[](SAI->segmentation);
	if (SAI->label_im != nullptr) {
		delete[](SAI->label_im);
	}

}