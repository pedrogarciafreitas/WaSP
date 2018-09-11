#include <vector>
#include <algorithm>

#include "motioncompensation.hh"
#include "ppm.hh"
#include "warping.hh"

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

	int search_radius = 10; // search -search_radius:search_radius in both x,y

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

	int *segp = view0->label_im;

	for (int ijk = 0; ijk < view0->nr*view0->nc; ijk++) {
		for (int iki = 0; iki < nregions; iki++) {
			int ik = view0->mv_regions.at(iki);
			if (*(segp + ijk) == ik) {
#pragma omp parallel for
				for (int isr = -search_radius; isr <= search_radius; isr++) {
					for (int isc = -search_radius; isc <= search_radius; isc++) {

						int iy = ijk % view0->nr; //row
						int ix = (ijk - iy) / view0->nr; //col

						int iy1;
						int ix1;

						warpSubscripts_from_View0_to_View1_int(view0, view1, iy + isr, ix + isc, iy1, ix1);

						if (iy1 >= 0 && iy1 < view0->nr && ix1 >= 0 && ix1 < view0->nc) {

							int ijk1 = ix1*view0->nr + iy1;

							for (int ic = 0; ic < 3; ic++) {
								int offc = ic*view0->nr*view0->nc;
								match_score[iki][isr + search_radius][isc + search_radius] += abs(((float)*(im0 + ijk + offc) - (float)*(im1 + ijk1 + offc)));
							}

							//printf("%f\n", match_score[iki][isr + search_radius][isc + search_radius]);

							counts[iki][isr + search_radius][isc + search_radius]++;

							//printf("%d\n", counts[iki][isr + search_radius][isc + search_radius]);

						}
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

					float clow = match_score[ik][isr][isc] / (float)nk;

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
		if (!(dy == 0 && dx == 0) && !(abs(dy) == search_radius && abs(dx) == search_radius)) {

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

		/*for (int iR = 0; iR < mv_regions.size(); iR++) {
			printf("iR=%i\tdy=%i\tdx=%i\n", mv_regions.at(iR).iR, mv_regions.at(iR).dy, mv_regions.at(iR).dx);
		}*/

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