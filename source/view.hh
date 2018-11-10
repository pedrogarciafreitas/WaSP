#ifndef VIEW_HH
#define VIEW_HH

#include <stdint.h>

#define MEDFILT_DEPTH false

struct view{

	unsigned short *color;
	unsigned short *depth;

	unsigned short *segmentation;

	int r, c; // SAI subscript

	int nr, nc; // image height, width

	float y, x; // camera displacement

	int min_inv_d; // needed only if inverse depth has negative values, [0,max]-mind = [-mind,max-mind]

	int n_references, n_depth_references;

	int *references, *depth_references; /* depth references not necessarily the same as color references
													  we can have, for example, depth warping only from the externally obtained depth but we still
													  warp color from neighbors that don't have depth provided. We don't want to propagate depth errors from
													  badly warped depth views, thus we restrict depth warping to some high quality subset (usually meaning the 
													  externally obtained high quality depth maps)*/

	signed short *merge_weights;
	int32_t *sparse_weights;

	unsigned char *sparse_mask;

	float *merge_weights_float;

	int *number_of_pixels_per_region;

	bool *bmask; /* view mask for merging weights */
	unsigned short *seg_vp; /* class segmentation, used for view merging weights */
	int NB;

	float residual_rate_color;
	float residual_rate_depth;

	float stdd;

	int NNt, Ms; //for global sparse, NNt defines the neighborhood size [ -NNt:NNt,-NNt:NNt ], Ms is the filter order

	int has_segmentation;
	int maxL; // number of regions in segmentation

	int ****region_displacements; /* region displacement vectors [iar][iac][iR][xy], e.g., [13][13][25][2], for 13x13 angular views with 25 regions for segmentation */

	char path_input_pgm[1024], path_input_ppm[1024], path_input_seg[1024];
	char path_out_pgm[1024], path_out_ppm[1024];
	char path_label_im[1024];

	float *DM_ROW, *DM_COL; /* for lenslet with region displacement vectors */

	int i_order; /* view position in encoding configuration */

	bool use_median; //use median merging or not

	bool yuv_transform;

	bool has_color_residual, has_depth_residual, use_global_sparse;
	bool has_color_references, has_depth_references;
	//bool has_min_inv_depth;

	bool has_x_displacement, has_y_displacement;

	bool has_chrominance;

	bool depth_file_exist;

	char ppm_residual_path[1024];
	char jp2_residual_path_jp2[1024];

	char pgm_residual_Y_path[1024];
	char jp2_residual_Y_path_jp2[1024];
	char pgm_residual_Cb_path[1024];
	char jp2_residual_Cb_path_jp2[1024];
	char pgm_residual_Cr_path[1024];
	char jp2_residual_Cr_path_jp2[1024];

	char *ycbcr_pgm_names[3];
	char *ycbcr_jp2_names[3];

	char pgm_residual_depth_path[1024];
	char jp2_residual_depth_path_jp2[1024];

	double final_psnr;
	double warp_psnr;
	double merge_psnr;
	double sparse_psnr;

	int ncomp;

	/* for regions */
	int32_t *label_im;
	int nregions;
	int *reg_histogram;

	char output_dir[1024];
	char input_dir[1024];


};

void initView(view* view);

bool loadColor(view* SAI);
bool loadInverseDepth(view* SAI);
void unloadColor(view* SAI);
void unloadInverseDepth(view* SAI);

bool loadLabels(view* SAI);
void unloadLabels(view *SAI);

int32_t *loadWarpedLabelIm(view *SAI, view *ref_view);
bool writeWarpedLabelIm(view *SAI, view *ref_view, const int32_t *warpedLabelIm);

#endif