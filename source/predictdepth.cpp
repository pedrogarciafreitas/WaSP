#include "predictdepth.hh"
#include "warping.hh"
#include "medianfilter.hh"
#include "ppm.hh"
#include "inpainting.hh"
#include "warping.hh"

#include <ctime>
#include <vector>

void predictDepth(view* SAI, view *LF)
{
	/* forward warp depth */
	if (SAI->has_depth_references) {
		/* currently we forward warp the depth from the N (for HDCA N = 5, lenslet maybe 1?) references */
		initializeWarpingArraysInverseDepth(SAI);

		for (int ij = 0; ij < SAI->n_depth_references; ij++)
		{
			view *ref_view = LF + SAI->depth_references[ij];

			int tmp_w, tmp_r, tmp_ncomp;

			aux_read16PGMPPM(ref_view->path_out_pgm, tmp_w, tmp_r, tmp_ncomp, ref_view->depth);
			//aux_read16PGMPPM(ref_view->path_out_ppm, tmp_w, tmp_r, tmp_ncomp, ref_view->color);

			warpView0_to_View1(ref_view, SAI, 1, ij,0);

			delete[](ref_view->depth);
			delete[](ref_view->color);

			ref_view->depth = nullptr;
			ref_view->color = nullptr;
		}

		/* merge depth using median*/

		int startt = clock();

#pragma omp parallel for
		for (int ij = 0; ij < SAI->nr*SAI->nc; ij++) {
			std::vector<unsigned short> depth_values;
			for (int uu = 0; uu < SAI->n_depth_references; uu++) {
				//for (int uu = 0; uu < SAI->n_references; uu++) {
				unsigned short *pp = SAI->warped_depth_views[uu];
				float *pf = SAI->occlusion_masks[uu];
				if (*(pf + ij) > INIT_DISPARITY_VALUE) {
					depth_values.push_back(*(pp + ij));
				}
			}
			if (depth_values.size() > 0)
				SAI->depth[ij] = getMedian(depth_values);
		}

		//std::cout << "time elapsed in depth merging\t" << (int)clock() - startt << "\n";

		/* hole filling for depth */
		holefilling(SAI->depth, 1, SAI->nr, SAI->nc, 0);

		deinitializeWarpingArraysInverseDepth(SAI);

	}

}
