#include "sparsefilter.hh"
#include "fastols.hh"
#include "bitdepth.hh"

#include <cmath>
#include <cstdio>
#include <vector>
#include <algorithm>

#define NULL 0

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
			//printf("%i\t%f\n", Regr0[ii], theta[Regr0[ii] - 1]);
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

	int Npp = (nr - NNt * 2)*(nc - NNt * 2);
	int Npp0 = Npp / 3;

	int MT = (NNt * 2 + 1)*(NNt * 2 + 1) + 1; /* number of regressors */

	double *AA = new double[Npp*MT]();
	double *Yd = new double[Npp]();

	for (int ii = 0; ii < Npp; ii++)
		*(AA + ii + (NNt * 2 + 1)*(NNt * 2 + 1)*Npp) = 1.0 / ((double)(1 << BIT_DEPTH) - 1);

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
							*(Yd + iiu + 0*Npp0) += ((double)*(original_color_view + offset)) / ( (double)(1 << BIT_DEPTH) - 1);// (pow(2, BIT_DEPTH) - 1);
						}

						/* get the regressors */
						*(AA + iiu + 0*Npp0 + ai*Npp) += ((double)*(pshort + offset)) / ( (double)(1 << BIT_DEPTH) - 1);// (pow(2, BIT_DEPTH) - 1);

					}
					*(AA + iiu + 0 * Npp0 + ai*Npp) = *(AA + iiu + 0 * Npp0 + ai*Npp) / 3;

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
			*(Yd + iiu + 0 * Npp0) = *(Yd + iiu + 0 * Npp0) / 3;
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
