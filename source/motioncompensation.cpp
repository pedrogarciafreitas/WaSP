#include <vector>
#include <algorithm>
#include <numeric>

#include "motioncompensation.hh"
#include "ppm.hh"
#include "warping.hh"
#include "medianfilter.hh"
#include "ycbcr.hh"
#include "bitdepth.hh"

void readLabelIm(view *view0){
	if (view0->label_im == NULL) {
		FILE *tmpfile_im_labels = fopen(view0->path_label_im, "rb");
		view0->label_im = new int[view0->nr*view0->nc]();
		fread(view0->label_im, sizeof(int), view0->nr*view0->nc, tmpfile_im_labels);
		fclose(tmpfile_im_labels);
	}
}

int findMVIndex(view *view0, view* view1) {
	int i_v = -1;
	for (int ii = 0; ii < view0->mv_views.size(); ii++) {
		if (view0->mv_views.at(ii).first == view1->i_order) {
			i_v = ii;
			break;
		}
	}
	return i_v;
}

void sortRegionsBySize(view *view0) {

	std::vector<std::pair<int, int>> reg_sz;

	for (int iR = 0; iR < view0->nregions; iR++) {

		std::pair<int, int> tmp_reg;
		tmp_reg.first = view0->reg_histogram[iR];
		tmp_reg.second = iR;

		reg_sz.push_back(tmp_reg);

	}

	std::sort(reg_sz.begin(), reg_sz.end());

	for (int iR = static_cast<int>(reg_sz.size()) - 1; iR >= 0; iR--) {

		if (view0->reg_histogram[iR] >= MV_REG_MIN_SIZE && view0->mv_regions.size() <= MV_MAX_REGS) {
			view0->mv_regions.push_back(iR);
		}

		if (view0->mv_regions.size() == MV_MAX_REGS) {
			break;
		}

	}

}

//
//void getMotionVectorsView0_to_View1(view *view0, view *view1) {
//
//	if (findMVIndex(view0, view1) >-1) {
//		printf("Motion vectors already available ...\n");
//		return;
//	}
//	else {
//		printf("Motion vector esimation between views %i (%i,%i) and %i (%i,%i) ...\n", view0->i_order,
//			view0->r, view0->c, view1->i_order, view1->r, view1->c);
//	}
//
//	int nregions = view0->nregions; /* number of regions in view0->label_im */
//
//	int ncomp1;
//
//	unsigned short *im0;
//	aux_read16PGMPPM(view0->path_input_ppm, view0->nc, view0->nr, ncomp1, im0);
//
//	unsigned short *im1;
//	aux_read16PGMPPM(view1->path_input_ppm, view1->nc, view1->nr, ncomp1, im1);
//
//	readLabelIm(view0);
//
//	int *segp = view0->label_im;
//
//	int **search_radiuses = new int*[nregions]();
//	for (int ik = 0; ik < nregions; ik++) {
//		search_radiuses[ik] = new int[2]();
//	}
//
//	std::vector<int> *sr_y = new std::vector<int>[nregions]();
//	std::vector<int> *sr_x = new std::vector<int>[nregions]();
//
//	for (int ijk = 0; ijk < view0->nr*view0->nc; ijk++) {
//
//		int ik = *(segp + ijk);
//
//		int iy = ijk % view0->nr; //row
//		int ix = (ijk - iy) / view0->nr; //col
//
//		int iy1 = iy;
//		int ix1 = ix;
//
//		warpSubscripts_from_View0_to_View1_int(view0, view1, iy, ix, iy1, ix1);
//
//		int dy = abs(iy - iy1);
//		int dx = abs(ix - ix1);
//
//		sr_y[ik].push_back(dy);
//		sr_x[ik].push_back(dx);
//
//	}
//
//	for (int ik = 0; ik < nregions; ik++) {
//		if ( sr_y[ik].size() > 0 ) {
//			search_radiuses[ik][0] = getMedian(sr_y[ik]);
//			search_radiuses[ik][1] = getMedian(sr_x[ik]);
//		}
//	}
//
//	std::vector<int> *mv_y = new std::vector<int>[nregions]();
//	std::vector<int> *mv_x = new std::vector<int>[nregions]();
//
//	for (int ijk = 0; ijk < view0->nr*view0->nc; ijk++) {
//
//		int ik = *(segp + ijk);
//
//		int search_radius_y = search_radiuses[ik][0] / 2;
//		int search_radius_x = search_radiuses[ik][1] / 2;
//
//		search_radius_y = search_radius_y < MIN_RADIUS ? MIN_RADIUS : search_radius_y;
//		search_radius_x = search_radius_x < MIN_RADIUS ? MIN_RADIUS : search_radius_x;
//
//		float best_match_score = FLT_MAX;
//		
//		int best_r = 0;
//		int best_c = 0;
//
//		int iy = ijk % view0->nr; //row
//		int ix = (ijk - iy) / view0->nr; //col
//
//		int iy1;
//		int ix1;
//
//		warpSubscripts_from_View0_to_View1_int(view0, view1, iy, ix, iy1, ix1);
//
//		for (int isr = -search_radius_y; isr <= search_radius_y; isr++) {
//			for (int isc = -search_radius_x; isc <= search_radius_x; isc++) {
//
//				bool extrm = abs(isr) == abs(search_radius_y) && abs(isc) == abs(search_radius_x);
//
//				int iy11 = iy1 + isr;
//				int ix11 = ix1 + isc;
//
//				if (iy11 >= 0 && iy11 < view0->nr && ix11 >= 0 && ix11 < view0->nc) {
//
//					int ijk1 = ix11*view0->nr + iy11;
//
//					float match_score = 0;
//
//					for (int ic = 0; ic < 3; ic++) {
//						int offc = ic*view0->nr*view0->nc;
//						match_score += abs(((float)*(im0 + ijk + offc) - (float)*(im1 + ijk1 + offc)));
//					}
//
//					if (match_score < best_match_score) {
//						best_match_score = match_score;
//						best_r = isr;
//						best_c = isc;
//					}
//
//				}
//			}
//		}
//
//		mv_y[ik].push_back(best_r);
//		mv_x[ik].push_back(best_c);
//	}
//
//	std::vector< MV_REGION > mv_regions_0;
//
//	/* find best matching displacement */
//
//	for (int iki = 0; iki < view0->mv_regions.size(); iki++) {
//
//		int ik = view0->mv_regions.at(iki);
//
//		int dy = getMode(mv_y[ik]);
//		int dx = getMode(mv_x[ik]);
//
//		if (abs(dy) > 0 || abs(dx) > 0) {
//
//			MV_REGION mv_region;
//			mv_region.iR = static_cast<unsigned int>(ik);
//			mv_region.dx = static_cast<signed char>(dx);
//			mv_region.dy = static_cast<signed char>(dy);
//
//			mv_regions_0.push_back(mv_region);
//		}
//
//	}
//
//	if (mv_regions_0.size() > 0) {
//
//		std::pair< unsigned short, std::vector<MV_REGION> > tmp_mv_reg;
//		tmp_mv_reg.first = static_cast<unsigned short>(view1->i_order);
//		tmp_mv_reg.second = mv_regions_0;
//
//		view0->mv_views.push_back(tmp_mv_reg);
//
//		for (int iR = 0; iR < mv_regions_0.size(); iR++) {
//			printf("iR=%i\tdy=%i\tdx=%i\n", mv_regions_0.at(iR).iR, mv_regions_0.at(iR).dy, mv_regions_0.at(iR).dx);
//		}
//
//	}
//
//	printf("Regions found: %i\n", static_cast<int>(mv_regions_0.size()));
//
//	delete[](im1);
//	delete[](im0);
//
//	delete[](mv_y);
//	delete[](mv_x);
//
//	delete[](sr_y);
//	delete[](sr_x);
//
//	for (int ik = 0; ik < nregions; ik++) {
//		delete[](search_radiuses[ik]);
//	}
//	delete[](search_radiuses);
//
//	if (view0->label_im != NULL) {
//		delete[](view0->label_im);
//		view0->label_im = NULL;
//	}
//}

void getMotionVectorsView0_to_View1(view *view0, view *view1) {

	if (findMVIndex(view0, view1) >-1 ) {
		printf("Motion vectors already available ...\n");
		return;
	}
	else {
		printf("Motion vector esimation between views %i (%i,%i) and %i (%i,%i) ...\n",view0->i_order,
			view0->r,view0->c, view1->i_order, view1->r, view1->c);
	}

	int nregions = static_cast<int>(view0->mv_regions.size());

	/* initialize region displacement table */
	int **region_displacements = new int*[nregions]();
	for (int iR = 0; iR < nregions; iR++) {
		region_displacements[iR] = new int[2]();
	}

	int ncomp1;

	unsigned short *im0;
	aux_read16PGMPPM(view0->path_input_ppm, view0->nc, view0->nr, ncomp1, im0);

	unsigned short *im1;
	aux_read16PGMPPM(view1->path_input_ppm, view1->nc, view1->nr, ncomp1, im1);

	int search_radius = 12; // search -search_radius:search_radius in both x,y

	float ***match_score = new float**[nregions]();
	for (int ik = 0; ik < nregions; ik++) {
		match_score[ik] = new float*[search_radius * 2 + 1]();
		for (int isr = 0; isr < search_radius * 2 + 1; isr++) {
			match_score[ik][isr] = new float[search_radius * 2 + 1]();
		}
	}

	int ***counts = new int**[nregions]();
	for (int ik = 0; ik < nregions; ik++) {
		counts[ik] = new int*[search_radius * 2 + 1]();
		for (int isr = 0; isr < search_radius * 2 + 1; isr++) {
			counts[ik][isr] = new int[search_radius * 2 + 1]();
		}
	}

	readLabelIm(view0);

	int *label_im = view0->label_im;

	int *lut_reg = new int[view0->nregions]();

	for (int ijk = 0; ijk < view0->nr*view0->nc; ijk++) {
		for (int iki = 0; iki < nregions; iki++) {
			int ik = view0->mv_regions.at(iki);
			if (*(label_im + ijk) == ik) {
				lut_reg[*(label_im + ijk)] = iki;
			}
		}
	}

	if (1) {
		unsigned short *ycbcr0 = new unsigned short[view0->nr*view0->nc * 3]();
		unsigned short *ycbcr1 = new unsigned short[view0->nr*view0->nc * 3]();


		RGB2YCbCr(im0, ycbcr0, view0->nr, view0->nc, BIT_DEPTH);
		RGB2YCbCr(im1, ycbcr1, view0->nr, view0->nc, BIT_DEPTH);

		delete[](im0);
		delete[](im1);

		im0 = ycbcr0;
		im1 = ycbcr1;
	}

	for (int ijk = 0; ijk < view0->nr*view0->nc; ijk++) {

		int iki = lut_reg[*(label_im + ijk)];

		int iy = ijk % view0->nr; //row
		int ix = (ijk - iy) / view0->nr; //col

		int iy1;
		int ix1;

		warpSubscripts_from_View0_to_View1_int(view0, view1, iy, ix, iy1, ix1);

#pragma omp parallel for
		for (int isr = -search_radius; isr <= search_radius; isr++) {
			for (int isc = -search_radius; isc <= search_radius; isc++) {

				int iy11 = iy1 + isr;
				int ix11 = ix1 + isc;

				if (iy11 >= 0 && iy11 < view0->nr && ix11 >= 0 && ix11 < view0->nc) {

					int ijk1 = ix11*view0->nr + iy11;

					for (int ic = 0; ic < 1; ic++) {
						int offc = ic*view0->nr*view0->nc;
						float err = static_cast<float>(*(im0 + ijk + offc)) - static_cast<float>(*(im1 + ijk1 + offc));
						float abse = abs(err);
						match_score[iki][isr + search_radius][isc + search_radius] += abse;
						counts[iki][isr + search_radius][isc + search_radius]++;
					}	

				}
			}
		}
	}

	std::vector< MV_REGION > mv_regions;

	/* find best matching displacement */
	for (int ik = 0; ik < nregions; ik++) {
		float lowest_score = FLT_MAX;
		for (int isr = 0; isr < search_radius * 2 + 1; isr++) {
			for (int isc = 0; isc < search_radius * 2 + 1; isc++) {

				int nk = counts[ik][isr][isc];

				if (nk > 0) {

					float clow = match_score[ik][isr][isc] / static_cast<float>( nk );

					//printf("%f\n", clow);

					if (clow < lowest_score) {
						region_displacements[ik][0] = isr - search_radius;
						region_displacements[ik][1] = isc - search_radius;
						lowest_score = clow;
					}

				}
			}
		}

		int dx = region_displacements[ik][1];
		int dy = region_displacements[ik][0];

		if (!(dy == 0 && dx == 0)) {

			MV_REGION mv_region;
			mv_region.iR = static_cast<unsigned int>( view0->mv_regions.at(ik) );
			mv_region.dx = static_cast<signed char>( dx );
			mv_region.dy = static_cast<signed char>( dy );

			mv_regions.push_back(mv_region);

		}

	}

	if (mv_regions.size() > 0) {

		std::pair< unsigned short, std::vector<MV_REGION> > tmp_mv_reg;
		tmp_mv_reg.first = static_cast<unsigned short>( view1->i_order );
		tmp_mv_reg.second = mv_regions;

		view0->mv_views.push_back(tmp_mv_reg);

		for (int iR = 0; iR < mv_regions.size(); iR++) {
			printf("iR=%i\tdy=%i\tdx=%i\n", mv_regions.at(iR).iR, mv_regions.at(iR).dy, mv_regions.at(iR).dx);
		}

	}

	printf("Regions found: %i\n", static_cast<int>(mv_regions.size()));

	for (int ik = 0; ik < nregions; ik++) {
		for (int isr = 0; isr < search_radius * 2 + 1; isr++) {
			delete[](match_score[ik][isr]);
			delete[](counts[ik][isr]);
		}
		delete[](match_score[ik]);
		delete[](counts[ik]);
	}
	delete[](match_score);
	delete[](counts);

	for (int iR = 0; iR < nregions; iR++) {
		delete[](region_displacements[iR]);
	}
	delete[](region_displacements);

	delete[](im1);
	delete[](im0);

	if (view0->label_im != NULL) {
		delete[](view0->label_im);
		view0->label_im = NULL;
	}
}