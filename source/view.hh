#ifndef VIEW_HH
#define VIEW_HH

#include <stdint.h>
#include <vector>

struct MV_REGION {
	unsigned int iR;
	signed char dy;
	signed char dx;
};

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

	char path_input_pgm[1024], path_input_ppm[1024], path_input_seg[1024];
	char path_out_pgm[1024], path_out_ppm[1024];
	char path_label_im[1024];

	float *DM_ROW, *DM_COL; /* for lenslet with region displacement vectors */

	int i_order; /* view position in encoding configuration */

	bool use_median; //use median merging or not

	bool yuv_transform;

	bool has_color_residual, has_depth_residual, use_global_sparse;
	bool has_color_references, has_depth_references;

	bool has_x_displacement, has_y_displacement; /* camera displacement information flag */

	/* for regions */
	int *label_im;
	int nregions;
	int *reg_histogram;

	bool use_region_sparse;
	bool use_motion_vectors;

	std::vector< std::vector< unsigned char > > region_Regr;
	std::vector< std::vector< int32_t > > region_Theta;

	std::vector< int > mv_regions;

	std::vector< std::pair< unsigned short, std::vector< MV_REGION > > > mv_views; /* first one is view index, second contains vector of regions for that view */

};

void initView(view* view);


#endif