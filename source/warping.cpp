#include "warping.hh"
#include "bitdepth.hh"

#include <cstdlib>
#include <cmath>

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

void warpView0_to_View1(view *view0, view *view1, unsigned short *&warpedColor, unsigned short *&warpedDepth, float *&DispTarg, int ncomp)
{

	/*this function forward warps from view0 to view1 for both color and depth*/

	float ddy = view0->y - view1->y;
	float ddx = view0->x - view1->x;

	unsigned short *AA1 = view0->color;
	unsigned short *DD1 = view0->depth;

	warpedColor = new unsigned short[view0->nr*view0->nc * 3]();
	warpedDepth = new unsigned short[view0->nr*view0->nc]();
	DispTarg = new float[view0->nr*view0->nc]();

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