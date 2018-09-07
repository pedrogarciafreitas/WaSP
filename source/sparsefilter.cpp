#include "sparsefilter.hh"
#include "fastols.hh"
#include "bitdepth.hh"

#include <cmath>
#include <cstdio>
#include <vector>
#include <algorithm>

#define NULL 0
#define MIN_REG_SZ 1024
#define NNt_reg 3
#define Ms_reg 12

void checkOutOfBounds(const int RR, const int CC, const int nr, const int nc, const int dy_0, const int dx_0, int &dy, int &dx ) {

	if (RR > nr - 1) {
		dy = dy_0 - (RR - (nr - 1));
	}
	if (RR < 0) {
		dy = dy_0 - RR;
	}
	if (CC > nc - 1) {
		dx = dx_0 - (CC - (nc - 1));
	}
	if (CC < 0) {
		dx = dx_0 - CC;
	}

}

int checkBoundaryPixels(int *label_im, int offset, const int ijk, const int NNt, const int ir, const int ic, const int dy, const int dx, const int nr, const int nc ) {
	if (*(label_im + offset) != ijk) { // if the pixel being selected is NOT part of the region

		double min_d_dy_dx = FLT_MAX;

		for (int dy1_0 = -NNt; dy1_0 <= NNt; dy1_0++) {
			for (int dx1_0 = -NNt; dx1_0 <= NNt; dx1_0++) {

				int RR1 = ir + dy + dy1_0;
				int CC1 = ic + dx + dx1_0;

				int dy1 = dy1_0;
				int dx1 = dx1_0;

				if (RR1 > nr - 1) {
					dy1 = dy1_0 - (RR1 - (nr - 1));
				}
				if (RR1 < 0) {
					dy1 = dy1_0 - RR1;
				}
				if (CC1 > nc - 1) {
					dx1 = dx1_0 - (CC1 - (nc - 1));
				}
				if (CC1 < 0) {
					dx1 = dx1_0 - CC1;
				}

				int offset2 = ir + dy + dy1 + nr*(ic + dx + dx1);

				if (*(label_im + offset2) == ijk) {
					double ddx = static_cast<double>(dx1) - static_cast<double>(dx);
					double ddy = static_cast<double>(dy1) - static_cast<double>(dy);
					double d_dydx = ddx*ddx + ddy*ddy;
					if (d_dydx < min_d_dy_dx) {
						min_d_dy_dx = d_dydx;
						offset = offset2;
					}
				}
			}
		}
	}

	return offset;
}

void applyRegionSparseFilter(view *view0) {

	//int MIN_REG_SZ = 256;

	int nregions = view0->nregions;
	int *reg_histogram = view0->reg_histogram;
	int *label_im = view0->label_im;

	int NNt = NNt_reg;
	int Ms = Ms_reg;

	int nr = view0->nr;
	int nc = view0->nc;

	unsigned short *pshort = view0->color;

	float *final_view = new float[nr*nc * 3]();

	for (int ii = 0; ii < nr*nc*3; ii++) {
		*(final_view + ii) = static_cast<float>( *(pshort + ii) );
	}

	int reg_i = 0;

	for (int ijk = 1; ijk <= nregions; ijk++) {

		if (reg_histogram[ijk] > MIN_REG_SZ) {

			std::vector< double > theta0 = view0->region_Theta.at(reg_i);

			std::vector< int > Regr0 = view0->region_Regr.at(reg_i);

			reg_i++;

			float *theta = new float[(2 * NNt + 1)*(2 * NNt + 1) + 1]();

			for (int ii = 0; ii < theta0.size(); ii++) {
				if (Regr0.at(ii) > 0) {
					theta[Regr0.at(ii) - 1] = static_cast<float>( theta0.at(ii) );
					//printf("%i\t%f\n", Regr0.at(ii), theta[Regr0.at(ii) - 1]);
				}
			}

			std::vector< unsigned int > inds;
			for (int ikj = 0; ikj < nr*nc; ikj++) {
				if (label_im[ikj] == ijk) {
					inds.push_back(ikj);
				}
			}

			for (int ee = 0; ee < inds.size(); ee++) {

				int ind = inds.at(ee);

				int ir = ind % nr; //row
				int ic = (ind - ir) / nr; //col

				int ie = 0;

				for (int icomp = 0; icomp < 3; icomp++) {
					*(final_view + ind + nr*nc*icomp) = 0;
				}

				for (int dy_0 = -NNt; dy_0 <= NNt; dy_0++) {
					for (int dx_0 = -NNt; dx_0 <= NNt; dx_0++) {

						int RR = ir + dy_0;
						int CC = ic + dx_0;

						int dy = dy_0;
						int dx = dx_0;

						checkOutOfBounds(RR, CC, nr, nc, dy_0, dx_0, dy, dx);

						int offset = ir + dy + nr*(ic + dx);
						
						offset = checkBoundaryPixels(label_im, offset, ijk, NNt, ir, ic, dy, dx, nr, nc);

						for (int icomp = 0; icomp < 3; icomp++) {
							final_view[ ind + icomp*nr*nc ] += theta[ie] * (static_cast<float>(pshort[ offset + icomp*nr*nc ]));
						}

						ie++;

					}
				}

				/* bias term */
				for (int icomp = 0; icomp < 3; icomp++) {
					final_view[ind + icomp*nr*nc] += theta[(2 * NNt + 1)*(2 * NNt + 1)];
				}

			}
		}
	}

	for (int ii = 0; ii < nr*nc * 3; ii++) {
		if (final_view[ii] < 0)
			final_view[ii] = 0;
		if (final_view[ii] > (1 << BIT_DEPTH) - 1) //(pow(2, BIT_DEPTH) - 1))
			final_view[ii] = (1 << BIT_DEPTH) - 1;// (pow(2, BIT_DEPTH) - 1);

		pshort[ii] = static_cast<unsigned short>(floor(final_view[ii] +0.5));
	}

	delete[](final_view);

}

void getRegionSparseFilter( view *view0, unsigned short *original_color_view ) {

	//int MIN_REG_SZ = 256;

	int nregions = view0->nregions;
	int *reg_histogram = view0->reg_histogram;
	int *label_im = view0->label_im;

	int nr = view0->nr;
	int nc = view0->nc;

	int NNt = NNt_reg;
	int Ms = Ms_reg;

	unsigned short *pshort = view0->color;

	for (int ijk = 1; ijk <= nregions; ijk++) {

		if ( reg_histogram[ijk] > MIN_REG_SZ ) {

			std::vector< unsigned int > inds;
			for (int ikj = 0; ikj < nr*nc; ikj++) {
				if (label_im[ikj] == ijk) {
					inds.push_back(ikj);
				}
			}

			int Npp = static_cast<int>(inds.size() * 3);
			int Npp0 = static_cast<int>(inds.size());

			int MT = (NNt * 2 + 1)*(NNt * 2 + 1) + 1; /* number of regressors */

			double *AA = new double[Npp*MT]();
			double *Yd = new double[Npp]();

			for (int ii = 0; ii < Npp; ii++) {
				*(AA + ii + (NNt * 2 + 1)*(NNt * 2 + 1)*Npp) = 1.0;
			}

			int iiu = 0;

			for (int ee = 0; ee < inds.size(); ee++) {

				int ind = inds.at(ee);

				int ir = ind % nr; //row
				int ic = (ind - ir) / nr; //col

				int ai = 0;
				for (int dy_0 = -NNt; dy_0 <= NNt; dy_0++) {
					for (int dx_0 = -NNt; dx_0 <= NNt; dx_0++) {

						int RR = ir + dy_0;
						int CC = ic + dx_0;

						int dy = dy_0;
						int dx = dx_0;

						checkOutOfBounds(RR, CC, nr, nc, dy_0, dx_0, dy, dx);

						int offset = ir + dy + nr*(ic + dx);

						/* get the desired Yd*/
						if (dy == 0 && dx == 0) {
							for (int icomp = 0; icomp < 3; icomp++) {
								*(Yd + iiu + icomp*Npp0) = ( static_cast<double>(*(original_color_view + offset + icomp*nr*nc))) / ( static_cast<double>( (1 << BIT_DEPTH) - 1) );
							}
						}

						offset = checkBoundaryPixels(label_im, offset, ijk, NNt, ir, ic, dy, dx, nr, nc);

						for (int icomp = 0; icomp < 3; icomp++) {
							/* get the regressors */
							*(AA + iiu + icomp*Npp0 + ai*Npp) = ( static_cast<double>( *(pshort + offset + icomp*nr*nc))) / ( static_cast<double>((1 << BIT_DEPTH) - 1));// (pow(2, BIT_DEPTH) - 1);
						}
						ai++;
					}
				}
				iiu++;
			}

			int *PredRegr0 = new int[MT]();
			double *PredTheta0 = new double[MT]();

			int Mtrue = FastOLS_new(&AA, &Yd, PredRegr0, PredTheta0, Ms, MT, MT, Npp);

			//printf("Region\t%i filter:\t", ijk);
			//for (int ri = 0; ri < Ms; ri++) {
			//	printf("\t%f", PredTheta0[ri]);
			//}
			//printf("\n");

			std::vector< int > PredRegr;
			std::vector< double > PredTheta;

			for (int ri = 0; ri < Ms; ri++) {
				PredRegr.push_back(PredRegr0[ri] + 1);
				PredTheta.push_back(PredTheta0[ri]);
			}

			view0->region_Regr.push_back(PredRegr);
			view0->region_Theta.push_back(PredTheta);

			if (AA != NULL) {
				delete[](AA);
			}
			if (Yd != NULL) {
				delete[](Yd);
			}

			delete[](PredRegr0);
			delete[](PredTheta0);

		}

	}
}

void applyGlobalSparseFilter(view *view0){

	unsigned char *Regr0 = view0->sparse_mask;
	int32_t *theta0 = view0->sparse_weights;

	int Ms = view0->Ms;
	int NNt = view0->NNt;
	int nr = view0->nr, nc = view0->nc;

	float *theta = new float[(2 * NNt + 1)*(2 * NNt + 1) + 1]();

	for (int ii = 0; ii < Ms; ii++){
		if (Regr0[ii] > 0){
			theta[Regr0[ii] - 1] = ((float)theta0[ii]) / (float)(1 << BIT_DEPTH_SPARSE);
			printf("%i\t%f\n", Regr0[ii], theta[Regr0[ii] - 1]);
		}
	}

	//for (int ii = 0; ii < (2 * NNt + 1)*(2 * NNt + 1) + 1; ii++) {
	//	printf("%f\n", theta[ii] );
	//}

	float *final_view = new float[nr*nc * 3]();

	unsigned short *pshort = view0->color;

	for (int ii = 0; ii < nr*nc * 3; ii++)
		final_view[ii] = pshort[ii];

	for (int rr = NNt; rr < nr - NNt; rr++){
		for (int cc = NNt; cc < nc - NNt; cc++)
		{
			for (int icomp = 0; icomp < 3; icomp++)
				final_view[rr + cc*nr + icomp*nr*nc] = 0;

			int ee = 0;

			for (int dy = -NNt; dy <= NNt; dy++){
				for (int dx = -NNt; dx <= NNt; dx++){
					for (int icomp = 0; icomp < 3; icomp++){
						final_view[rr + cc*nr + icomp*nr*nc] += theta[ee] * ((float)pshort[rr + dy + (cc + dx)*nr + icomp*nr*nc]);
					}
					ee++;
				}
			}

			/* bias term */
			for (int icomp = 0; icomp < 3; icomp++){
				final_view[rr + cc*nr + icomp*nr*nc] += theta[(2 * NNt + 1)*(2 * NNt + 1)];
			}

		}
	}


	for (int ii = 0; ii < nr*nc * 3; ii++){
		if (final_view[ii] < 0)
			final_view[ii] = 0;
		if (final_view[ii] > (1 << BIT_DEPTH) - 1) //(pow(2, BIT_DEPTH) - 1))
			final_view[ii] = (1 << BIT_DEPTH) - 1;// (pow(2, BIT_DEPTH) - 1);

		pshort[ii] = (unsigned short)(final_view[ii]);
	}

	delete[](theta);
	delete[](final_view);

}

void getGlobalSparseFilter(view *view0, unsigned short *original_color_view)
{

	int nr = view0->nr;
	int nc = view0->nc;
	int NNt = view0->NNt;
	int Ms = view0->Ms;

	unsigned char *Regr0 = new unsigned char[Ms]();
	int32_t *theta0 = new int32_t[Ms]();

	view0->sparse_weights = theta0;
	view0->sparse_mask = Regr0;

	int Npp = (nr - NNt * 2)*(nc - NNt * 2) * 3;
	int Npp0 = Npp / 3;

	int MT = (NNt * 2 + 1)*(NNt * 2 + 1) + 1; /* number of regressors */

	double *AA = new double[Npp*MT]();
	double *Yd = new double[Npp]();

	for (int ii = 0; ii < Npp; ii++)
		*(AA + ii + (NNt * 2 + 1)*(NNt * 2 + 1)*Npp) = 1.0;

	int iiu = 0;

	unsigned short *pshort = view0->color;

	//double *PHI = new double[MT*MT]();
	//double *PSI = new double[MT]();

	//int offs = 2 * NNt + 1;

	for (int ir = NNt; ir < nr - NNt; ir++){
		for (int ic = NNt; ic < nc - NNt; ic++){
			int ai = 0;
			for (int dy = -NNt; dy <= NNt; dy++){
				for (int dx = -NNt; dx <= NNt; dx++) {
					for (int icomp = 0; icomp < 3; icomp++) {

						int offset = ir + dy + nr*(ic + dx) + icomp*nr*nc;

						/* get the desired Yd*/
						if (dy == 0 && dx == 0) {
							*(Yd + iiu + icomp*Npp0) = ((double)*(original_color_view + offset)) / ( (double)(1 << BIT_DEPTH) - 1);// (pow(2, BIT_DEPTH) - 1);
						}

						/* get the regressors */
						*(AA + iiu + icomp*Npp0 + ai*Npp) = ((double)*(pshort + offset)) / ( (double)(1 << BIT_DEPTH) - 1);// (pow(2, BIT_DEPTH) - 1);

					}
					ai++;


					//for (int dy1 = -NNt; dy1 <= NNt; dy1++) {
					//	for (int dx1 = -NNt; dx1 <= NNt; dx1++) {
					//		for (int icomp = 0; icomp < 3; icomp++) {

					//			int offset_1 = ir + dy + nr*(ic + dx) + icomp*nr*nc;
					//			int offset_2 = ir + dy1 + nr*(ic + dx1) + icomp*nr*nc;

					//			double val1 = ((double)*(pshort + offset_1)) / (pow(2, BIT_DEPTH) - 1);
					//			double val2 = ((double)*(pshort + offset_2)) / (pow(2, BIT_DEPTH) - 1);

					//			double des = ((double)*(original_color_view + offset_1)) / (pow(2, BIT_DEPTH) - 1);

					//			PHI[(dy+NNt) + offs*(dx+NNt) + MT*( (dy1+NNt) + offs*(dx1+NNt))] += val1*val2;
					//			PSI[(dy+NNt) + offs*(dx+NNt)] += val1*des;

					//		}
					//	}
					//}
				}
			}
			iiu++;
		}
	}

	int *PredRegr0 = new int[MT]();
	double *PredTheta0 = new double[MT]();

	int Mtrue = FastOLS_new(&AA, &Yd, PredRegr0, PredTheta0, Ms, MT, (NNt * 2 + 1)*(NNt * 2 + 1) + 1, Npp);
	//int Mtrue = FastOLS_new(AA, Yd, PredRegr0, PredTheta0, Ms, MT, (NNt * 2 + 1)*(NNt * 2 + 1) + 1, Npp,PHI,PSI);

	if (AA != NULL) {
		delete[](AA);
	}
	if (Yd != NULL) {
		delete[](Yd);
	}

	for (int ii = 0; ii < Ms; ii++){
		*(Regr0 + ii) = ((unsigned char)*(PredRegr0 + ii) + 1);
		*(theta0 + ii) = (int32_t)floor(*(PredTheta0 + ii) * (int32_t)(1 << BIT_DEPTH_SPARSE) + 0.5 );// pow(2, BIT_DEPTH_SPARSE));
	}

	delete[](PredRegr0);
	delete[](PredTheta0);

	/* sorting of filter coeffs and sparsity mask as increasing order of regressor index, new for ES2.4 */

	std::vector<std::pair<unsigned char, int32_t>> sparsefilter;

	for (int ij = 0; ij < view0->Ms; ij++) {
		std::pair<unsigned char, int32_t> tmp_sp;
		tmp_sp.first = view0->sparse_mask[ij];
		tmp_sp.second = view0->sparse_weights[ij];
		sparsefilter.push_back(tmp_sp);
	}

	sort(sparsefilter.begin(), sparsefilter.end()); /*ascending sort based on regressor index (e.g., integer between 1 and 50 since we have 50 regressors */

	memset(view0->sparse_mask, 0x00, sizeof(unsigned char)*view0->Ms);
	memset(view0->sparse_weights, 0x00, sizeof(int32_t)*view0->Ms);

	for (int ij = 0; ij < view0->Ms; ij++) {
		view0->sparse_mask[ij] = sparsefilter.at(ij).first;
		view0->sparse_weights[ij] = sparsefilter.at(ij).second;
	}

}
