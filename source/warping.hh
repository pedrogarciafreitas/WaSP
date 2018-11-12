#ifndef WARPING_HH
#define WARPING_HH

#include "view.hh"

#define INIT_DISPARITY_VALUE -10000.0 /* Bug fix from VM1.0. Introduced because -1.0 is no longer a good initial value since we can have negative disparity as well.*/

void warp_all_references_to_current_view(view *SAI);

void warpView0_to_View1(view *view0, view *view1, const int ncomp, const int ref_idx);

void getDisparity(const float y0, const float y1, const float x1, const float x2, float *inverse_depth_view_0, const int nr, const int nc,
	float *DM_COL, float *DM_ROW);

template <class T>
void warp_0_to_1( const float *DM_ROW_0, const float *DM_COL_0, const T *im0, const int nr, const int nc, const int ncomp,
	T *warped_im0, float *warped_inverse_depth0, const float *inverse_depth0) {

	for (int ij = 0; ij < nr*nc; ij++) {
		warped_inverse_depth0[ij] = INIT_DISPARITY_VALUE;
	}

	for (int ij = 0; ij < nr*nc; ij++)
	{
		int iy = ij % nr; //row
		int ix = (ij - iy) / nr; //col

		int ixnew = ix + (int)floor(DM_COL_0[ij] + 0.5);
		int iynew = iy + (int)floor(DM_ROW_0[ij] + 0.5);

		if (ixnew >= 0 && ixnew < nc && iynew >= 0 && iynew < nr) {

			float disp = inverse_depth0[ij];

			int indnew = iynew + ixnew*nr;

			if (warped_inverse_depth0[indnew] < disp) {
				warped_inverse_depth0[indnew] = disp;
				for (int icomp = 0; icomp < ncomp; icomp++) {
					warped_im0[indnew + nr*nc*icomp] = im0[ij + nr*nc*icomp];
				}
			}
		}
	}
}

#endif