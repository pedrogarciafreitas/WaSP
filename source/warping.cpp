#include <cstdlib>
#include <cmath>
#include <cstdio>

#include "ppm.hh"
#include "warping.hh"
#include "bitdepth.hh"

#define SAVE_PARTIAL_WARPED_VIEWS false

void getDisparity(const float y0, const float y1, const float x0, const float x1, float *inverse_depth_view_0, const int nr, const int nc,
	float *DM_COL, float *DM_ROW) {

	float ddx = y0 - y1;
	float ddy = x0 - x1;

//#pragma omp parallel
	for (int ij = 0; ij < nr*nc; ij++)
	{
		float disp = inverse_depth_view_0[ij];
		DM_COL[ij] = disp*ddx;
		DM_ROW[ij] = disp*ddy;
	}

}

void warp_all_references_to_current_view(view *SAI) {

	for (int ij = 0; ij < SAI->n_references; ij++)
	{

		view *ref_view = SAI->color_reference_views[ij];

		loadColor(ref_view);
		loadInverseDepth(ref_view);

		/* FORWARD warp color AND depth */
		warpView0_to_View1(ref_view, SAI, 3, ij);

		unloadColor(ref_view);
		unloadInverseDepth(ref_view);

		if (SAVE_PARTIAL_WARPED_VIEWS) {

			char tmp_str[1024];

			sprintf(tmp_str, "%s/%03d_%03d%s%03d_%03d%s", SAI->output_dir, (ref_view)->c, (ref_view)->r, "_warped_to_", SAI->c, SAI->r, ".ppm");
			aux_write16PGMPPM(tmp_str, SAI->nc, SAI->nr, 3, SAI->warped_color_views[ij]);

			sprintf(tmp_str, "%s/%03d_%03d%s%03d_%03d%s", SAI->output_dir, (ref_view)->c, (ref_view)->r, "_warped_to_", SAI->c, SAI->r, ".pgm");
			aux_write16PGMPPM(tmp_str, SAI->nc, SAI->nr, 1, SAI->warped_depth_views[ij]);

		}

	}
}

void warpView0_to_View1(view *view0, view *view1, const int ncomp, const int ref_idx)
{

	/*this function forward warps from view0 to view1 for both color and depth*/

	float ddy = view0->y - view1->y;
	float ddx = view0->x - view1->x;

	unsigned short *AA1 = view0->color;
	unsigned short *DD1 = view0->depth;

	float *DispTarg = view1->occlusion_masks[ref_idx];
	unsigned short *warpedColor = view1->warped_color_views[ref_idx];
	unsigned short *warpedDepth = view1->warped_depth_views[ref_idx];

	for (int ij = 0; ij < view0->nr*view0->nc; ij++) {
		DispTarg[ij] = INIT_DISPARITY_VALUE;
	}

	for (int ij = 0; ij < view0->nr*view0->nc; ij++)
	{
		float disp = ((float)DD1[ij] - (float)view0->min_inv_d) / (float)(1 << D_DEPTH);// pow(2, D_DEPTH);
		float DM_COL = disp*ddx;
		float DM_ROW = disp*ddy;

		int iy = ij % view0->nr; //row
		int ix = (ij - iy) / view0->nr; //col

		int ixnew = ix + (int)floor(DM_COL + 0.5);
		int iynew = iy + (int)floor(DM_ROW + 0.5);

		if (ixnew >= 0 && ixnew < view0->nc && iynew >= 0 && iynew < view0->nr) {
			int indnew = iynew + ixnew*view0->nr;
			if (DispTarg[indnew] < disp) { /* Bug fix from VM1.0. Since we use also negative disparity,
										   but the largest positive value still represents nearest pixel,
										   we should compare against disp not disp0.*/
				DispTarg[indnew] = disp;
				warpedDepth[indnew] = DD1[ij];

				for (int icomp = 0; icomp < ncomp; icomp++) {
					warpedColor[indnew + view0->nr*view0->nc*icomp] = AA1[ij+view0->nr*view0->nc*icomp];
				}
				
			}
		}

	}

}