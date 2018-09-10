#include "warping.hh"
#include "bitdepth.hh"
#include "motioncompensation.hh"

#include <cstdlib>
#include <cmath>


void warpSubscripts_from_View0_to_View1_int(view *view0, view *view1, int iy, int ix, int &iynew, int &ixnew) {

	float ddy = view0->y - view1->y;
	float ddx = view0->x - view1->x;

	unsigned short *DD1 = view0->depth;

	int ij = view0->nr*ix + iy;

	if (ij >= 0 && ij < view0->nr*view0->nc) {

		float disp = ((float)DD1[ij] - (float)view0->min_inv_d) / (float)(1 << D_DEPTH);
		float DM_COL = disp*ddx;
		float DM_ROW = -disp*ddy;

		ixnew = ix + (int)floor(DM_COL + 0.5);
		iynew = iy + (int)floor(DM_ROW + 0.5);

	}
	else {
		ixnew = -1;
		iynew = -1;
	}

}

void warpView0_to_View1(view *view0, view *view1, unsigned short *&warpedColor, unsigned short *&warpedDepth, float *&DispTarg)
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

	float *DM_COL_MV = new float[view0->nr*view0->nc]();
	float *DM_ROW_MV = new float[view0->nr*view0->nc]();

	if (MOTION_VECTORS) {

		getMotionVectorsView0_to_View1(view0, view1);

		int i_v = findMVIndex(view0, view1);

		if (i_v > -1) {

			std::vector<MV_REGION> mv_regions_final = view0->mv_views.at(i_v).second;

			for (int ij = 0; ij < view0->nr*view0->nc; ij++) {
				int ik = *(view0->label_im + ij);
				for (int iR = 0; iR < mv_regions_final.size(); iR++) {
					if (mv_regions_final.at(iR).iR == ik) {
						*(DM_COL_MV + ij) += mv_regions_final.at(iR).dx;
						*(DM_ROW_MV + ij) += mv_regions_final.at(iR).dy;
						break;
					}
				}
			}

		}
	}

	for (int ij = 0; ij < view0->nr*view0->nc; ij++)
	{

		float disp = ((float)DD1[ij] - (float)view0->min_inv_d) / (float)(1 << D_DEPTH);// pow(2, D_DEPTH);
		float DM_COL = disp*ddx + *(DM_COL_MV+ij);
		float DM_ROW = -disp*ddy + *(DM_ROW_MV + ij);

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
				warpedColor[indnew] = AA1[ij];
				warpedColor[indnew + view0->nr*view0->nc] = AA1[ij + view0->nr*view0->nc];
				warpedColor[indnew + 2 * view0->nr*view0->nc] = AA1[ij + 2 * view0->nr*view0->nc];
				warpedDepth[indnew] = DD1[ij];
			}
		}

	}

}