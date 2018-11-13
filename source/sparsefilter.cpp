#include "sparsefilter.hh"
#include "fastols.hh"
#include "bitdepth.hh"

#include <cmath>
#include <cstdio>
#include <vector>
#include <algorithm>

bool sortinrev(const std::pair<int, int> &a,
	const std::pair<int, int> &b)
{
	return (a.first > b.first);
}

std::vector<std::pair<int, int>> sortRegionsBySize(view *view0) {

	std::vector<std::pair<int, int>> regions_sorted_by_size_descending_order;

	for (int iR = 0; iR < view0->nregions; iR++) {

		std::pair<int, int> tmp_reg;
		tmp_reg.first = view0->reg_histogram[iR];
		tmp_reg.second = iR;

		regions_sorted_by_size_descending_order.push_back(tmp_reg);

	}

	std::sort(regions_sorted_by_size_descending_order.begin(), 
		regions_sorted_by_size_descending_order.end(), 
		sortinrev );

	return regions_sorted_by_size_descending_order;
}

void sortSparseFilter(const int Ms, unsigned char *sparse_mask, int32_t *sparse_weights) {

	/* sorting of filter coeffs (and sparsity mask) as increasing order of regressor index, new for ES2.4 */

	std::vector<std::pair<unsigned char, int32_t>> sparsefilter;

	for (int ij = 0; ij < Ms; ij++) {
		std::pair<unsigned char, int32_t> tmp_sp;
		tmp_sp.first = sparse_mask[ij];
		tmp_sp.second = sparse_weights[ij];
		sparsefilter.push_back(tmp_sp);
	}

	sort(sparsefilter.begin(), sparsefilter.end()); /*ascending sort based on regressor index (e.g., integer between 1 and 50 since we have 50 regressors */

	memset(sparse_mask, 0x00, sizeof(unsigned char)*Ms);
	memset(sparse_weights, 0x00, sizeof(int32_t)*Ms);

	for (int ij = 0; ij < Ms; ij++) {
		sparse_mask[ij] = sparsefilter.at(ij).first;
		sparse_weights[ij] = sparsefilter.at(ij).second;
	}

}

int checkOutOfBounds(const int RR, const int nr, const int dy_0) {

	/* prevents selected pixel from being outside image dimensions */

	int dy = dy_0;

	if (RR > nr - 1) {
		dy = dy_0 - (RR - (nr - 1));
	}
	if (RR < 0) {
		dy = dy_0 - RR;
	}

	return dy;
}

int checkBoundaryPixels(const int32_t *label_im, int offset, const int region_indx, const int NNt, const int ir, const int ic, const int dy, const int dx, const int nr, const int nc) {
	
	/* if selected pixel is outside the region this function selects closest pixel that is also inside the region.*/

	if (*(label_im + offset) != region_indx) { // if the pixel being selected is NOT part of the region

		double min_d_dy_dx = FLT_MAX;

		for (int dy1_0 = -NNt; dy1_0 <= NNt; dy1_0++) {
			for (int dx1_0 = -NNt; dx1_0 <= NNt; dx1_0++) {

				int RR1 = ir + dy + dy1_0;
				int CC1 = ic + dx + dx1_0;

				int dy1 = checkOutOfBounds(RR1, nr, dy1_0);
				int dx1 = checkOutOfBounds(CC1, nc, dx1_0);

				int offset2 = ir + dy + dy1 + nr*(ic + dx + dx1);

				if (*(label_im + offset2) == region_indx) {
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

unsigned short *applySparseFilterForOneRegion(view *view0, region_sparse_filter reg_sp) {

	int NNt = view0->NNt;
	int Ms = view0->Ms;

	int nr = view0->nr;
	int nc = view0->nc;

	unsigned short *pshort = view0->color;

	double *final_view = new double[nr*nc * 3]();

	unsigned short *output_image = new unsigned short[nr*nc*3]();

	memcpy(output_image, pshort, sizeof(unsigned short)*nr*nc * 3);

	for (int ii = 0; ii < nr*nc * 3; ii++) {
		*(final_view + ii) = static_cast<double>(*(pshort + ii));
	}

	std::vector< int32_t > theta0 = reg_sp.prediction_coefficients;
	std::vector< unsigned char > Regr0 = reg_sp.regressor_indices;

	int region_indx = reg_sp.region_indx;

	double *theta = new double[(2 * NNt + 1)*(2 * NNt + 1) + 1]();

	for (int ii = 0; ii < theta0.size(); ii++) {
		if (Regr0.at(ii) > 0) {
			theta[Regr0.at(ii) - 1] = static_cast<double>(theta0.at(ii)) / static_cast<double>(1 << BIT_DEPTH_SPARSE);
			//printf("%i\t%f\n", Regr0.at(ii), theta[Regr0.at(ii) - 1]);
		}
	}

	std::vector< unsigned int > inds;
	for (int ikj = 0; ikj < nr*nc; ikj++) {
		if (view0->label_im[ikj] == region_indx) {
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

				int dy = checkOutOfBounds(RR, nr, dy_0);
				int dx = checkOutOfBounds(CC, nc, dx_0);

				int offset = ir + dy + nr*(ic + dx);

				offset = checkBoundaryPixels(view0->label_im, offset, region_indx, NNt, ir, ic, dy, dx, nr, nc);

				for (int icomp = 0; icomp < 3; icomp++) {
					final_view[ind + icomp*nr*nc] += theta[ie] * static_cast<double>(pshort[offset + icomp*nr*nc]);
				}

				ie++;

			}
		}

		/* bias term */
		for (int icomp = 0; icomp < 3; icomp++) {
			final_view[ind + icomp*nr*nc] += theta[(2 * NNt + 1)*(2 * NNt + 1)];
		}

		for (int icomp = 0; icomp < 3; icomp++) {
			double tmp_val = final_view[ind + icomp*nr*nc];
			tmp_val = tmp_val < 0 ? 0 : tmp_val;
			tmp_val = tmp_val > (1 << BIT_DEPTH) - 1 ? (1 << BIT_DEPTH) - 1 : tmp_val;
			tmp_val = floor(tmp_val + 0.5);
			*(output_image + ind + icomp*nr*nc) = static_cast<unsigned short>(tmp_val);
		}

	}

	delete[](final_view);

	return output_image;
}

region_sparse_filter getSparseFilterForOneRegion(view *view0, const int region_indx) {

	int nr = view0->nr;
	int nc = view0->nc;
	int NNt = view0->NNt;
	int Ms = view0->Ms;

	std::vector< int32_t > inds;
	for (int ikj = 0; ikj < nr*nc; ikj++) {
		if (view0->label_im[ikj] == region_indx) {
			inds.push_back(ikj);
		}
	}

	int Npp = static_cast<int>(inds.size());

	int MT = (NNt * 2 + 1)*(NNt * 2 + 1) + 1; /* number of regressors */

	double *AA = new double[Npp*MT]();
	double *Yd = new double[Npp]();

	for (int ii = 0; ii < Npp; ii++) {
		*(AA + ii + (NNt * 2 + 1)*(NNt * 2 + 1)*Npp) = 1.0 / static_cast<double>((1 << BIT_DEPTH) - 1);
	}

	int iiu = 0;

	for (int ee = 0; ee < inds.size(); ee++) {

		int ind = inds.at(ee);

		for (int icomp = 0; icomp < 3; icomp++) {
			*(Yd + iiu) += static_cast<double>(*(view0->original_color_view + ind + icomp*nr*nc));
		}
		*(Yd + iiu) = *(Yd + iiu) / 3.0 / static_cast<double>((1 << BIT_DEPTH) - 1);

		int ir = ind % nr; //row
		int ic = (ind - ir) / nr; //col

		int ai = 0;
		for (int dy_0 = -NNt; dy_0 <= NNt; dy_0++) {
			for (int dx_0 = -NNt; dx_0 <= NNt; dx_0++) {

				int RR = ir + dy_0;
				int CC = ic + dx_0;

				int dy = checkOutOfBounds(RR, nr, dy_0);
				int dx = checkOutOfBounds(CC, nc, dx_0);

				int offset = ir + dy + nr*(ic + dx);

				offset = checkBoundaryPixels(view0->label_im, offset, region_indx, NNt, ir, ic, dy, dx, nr, nc);

				for (int icomp = 0; icomp < 3; icomp++) {
					*(AA + iiu + ai*Npp) += static_cast<double>(*(view0->color + offset + icomp*nr*nc));
				}
				*(AA + iiu + ai*Npp) = *(AA + iiu + ai*Npp) / 3.0 / static_cast<double>((1 << BIT_DEPTH) - 1);
				ai++;
			}
		}
		iiu++;
	}

	int *PredRegr0 = new int[MT]();
	double *PredTheta0 = new double[MT]();

	int Mtrue = FastOLS_new(&AA, &Yd, PredRegr0, PredTheta0, Ms, MT, MT, Npp);

	std::vector< unsigned char > PredRegr;
	std::vector< int32_t > PredTheta;

	for (int ri = 0; ri < Ms; ri++) {
		PredRegr.push_back(PredRegr0[ri] + 1);
		PredTheta.push_back( static_cast<int32_t>(floor( PredTheta0[ri] * static_cast<double>(1 << BIT_DEPTH_SPARSE) + 0.5 )) );
	}

	sortSparseFilter(Ms, PredRegr.data(), PredTheta.data());

	region_sparse_filter reg_sp;

	reg_sp.region_indx = region_indx;
	reg_sp.prediction_coefficients = PredTheta;
	reg_sp.regressor_indices = PredRegr;

	if (AA != nullptr) {
		delete[](AA);
	}
	if (Yd != nullptr) {
		delete[](Yd);
	}

	delete[](PredRegr0);
	delete[](PredTheta0);

	return reg_sp;

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

void getGlobalSparseFilter(view *view0)
{

	int nr = view0->nr;
	int nc = view0->nc;
	int NNt = view0->NNt;
	int Ms = view0->Ms;

	unsigned short *original_color_view = view0->original_color_view;

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

	if (AA != nullptr) {
		delete[](AA);
	}
	if (Yd != nullptr) {
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
