#include "merging.hh"
#include "medianfilter.hh"
#include "fastols.hh"
#include "bitdepth.hh"
#include "warping.hh"

#include <cstring>

void getViewMergingLSWeights_N(view *view0)
{

	/* This function puts the LS view merging weights into LSw */

	int n_references = view0->n_references;
	int nr = view0->nr;
	int nc = view0->nc;

	signed short *LScoeffs = (view0)->merge_weights;

	int MMM = 1 << n_references;// pow(2, n_references);
	
	bool *bmask = view0->bmask;

	unsigned short *seg_vp = view0->seg_vp;

	/* go through all regions, collect regressors from references */
	unsigned short **reference_view_pixels_in_classes = new unsigned short*[MMM * n_references]();

	/* also collect desired values */
	unsigned short **original_view_in_classes = new unsigned short*[MMM]();

	int *number_of_pixels_per_region = view0->number_of_pixels_per_region;

	for (int ij = 1; ij < MMM; ij++){ // region

		int NN = number_of_pixels_per_region[ij];

		original_view_in_classes[ij] = new unsigned short[NN * 3]();
		unsigned short *p3s = original_view_in_classes[ij];

		int jj = 0;

		for (int ii = 0; ii < nr*nc; ii++){
			if (seg_vp[ii] == ij){
				for (int icomp = 0; icomp < 3; icomp++)
					*(p3s + jj + icomp*NN) = *(view0->original_color_view + ii + icomp*nr*nc);
				jj++;
			}
		}

		for (int ik = 0; ik < n_references; ik++){ // reference view

			if (bmask[ij + ik * MMM])
			{

				/* allocate pixels for region */
				reference_view_pixels_in_classes[ij + ik * MMM] = new unsigned short[NN * 3]();

				unsigned short *ps = reference_view_pixels_in_classes[ij + ik * MMM];
				unsigned short *pss = view0->warped_color_views[ik];


				jj = 0;

				for (int ii = 0; ii < nr*nc; ii++){
					if (seg_vp[ii] == ij){
						for (int icomp = 0; icomp < 3; icomp++)
							*(ps + jj + icomp*NN) = *(pss + ii + icomp*nr*nc);
						jj++;
					}
				}

			}

		}

	}


	/* run fastOLS on the classes */

	double *thetas = new double[MMM * n_references]();

	for (int ij = 1; ij < MMM; ij++){

		/* form A for this class, compute A'*A (phi)
		also compute A'*y (psi), where y is the desired data from the original view */

		int M = 0;

		for (int ik = 0; ik<n_references; ik++)
		{
			if (bmask[ij + MMM * ik])
			//if (( (ij - 1) >> ik) & 1)
			{
				M++;
			}
		}

		int N = number_of_pixels_per_region[ij]; // number of rows in A

		double *AA = new double[N*M]();
		double *Yd = new double[N]();

		unsigned short *ps;

		int ikk = 0;

		for (int ik = 0; ik < n_references; ik++){
			if (bmask[ij + ik * MMM]){
				ps = reference_view_pixels_in_classes[ij + ik * MMM];
				for (int ii = 0; ii < N; ii++){
					for (int icomp = 0; icomp < 3; icomp++) {
						*(AA + ii + ikk*N) += ((double)*(ps + ii + icomp*N)) / (double)((1 << BIT_DEPTH) - 1);
					}
					*(AA + ii + ikk*N) = *(AA + ii + ikk*N) / 3;
				}
				ikk++;
			}
		}

		ps = original_view_in_classes[ij];

		for (int ii = 0; ii < N; ii++){
			for (int icomp = 0; icomp < 3; icomp++) {
				*(Yd + ii) += ((double)*(ps + ii + icomp*N)) / (double)((1 << BIT_DEPTH) - 1);
			}
			*(Yd + ii) = *(Yd + ii) / 3;
		}

		/* fastols */

		int *PredRegr0 = new int[M]();
		double *PredTheta0 = new double[M]();

		//int Mtrue = FastOLS(ATA, ATYd, YdTYd, PredRegr0, PredTheta0, M, M, M);

		int Mtrue = FastOLS_new(&AA, &Yd, PredRegr0, PredTheta0, M, M, M, N);

		if (AA != nullptr) {
			delete[](AA);
		}
		if (Yd != nullptr) {
			delete[](Yd);
		}

		/* establish the subset of reference views available for class */
		int *iks = new int[M]();
		int ee = 0;
		for (int ik = 0; ik < n_references; ik++){
			if (bmask[ij + ik * MMM]){
				*(iks + ee) = ik;
				ee++;
			}
		}

		for (int ii = 0; ii < M; ii++){
			thetas[ij + MMM * iks[PredRegr0[ii]]] = *(PredTheta0 + ii);// (signed short)floor(*(PredTheta0 + ii)*(signed short)(1 << BIT_DEPTH_MERGE) + 0.5);// pow(2, BIT_DEPTH_MERGE) + 0.5);
		}

		delete[](iks);

		delete[](PredRegr0);
		delete[](PredTheta0);
		


	}

	columnwise_collecting_of_thetas(view0, thetas);

	delete[](thetas);

	for (int ij = 0; ij < MMM; ij++){

		delete[](original_view_in_classes[ij]);

		for (int ik = 0; ik < n_references; ik++){
			delete[](reference_view_pixels_in_classes[ij + MMM * ik]);
		}

	}

	delete[](original_view_in_classes);
	delete[](reference_view_pixels_in_classes);

}

void initSegVp(view *view0) {

	int nr = view0->nr;
	int nc = view0->nc;
	int n_references = view0->n_references;

	view0->seg_vp = new unsigned short[nr*nc]();

	int MMM = 1 << view0->n_references;

	int *number_of_pixels_per_region = new int[MMM]();

	for (int ii = 0; ii < nr*nc; ii++) {

		unsigned short ci = 0;

		for (int ik = 0; ik < n_references; ik++) {
			float *pf = view0->occlusion_masks[ik];
			if (*(pf + ii) > INIT_DISPARITY_VALUE)
				ci = ci + (unsigned short)( 1 << ik );// pow(2, ik);
		}

		view0->seg_vp[ii] = ci;

		number_of_pixels_per_region[ci]++;

	}

	view0->number_of_pixels_per_region = number_of_pixels_per_region;
}

void initViewW(view *view0) {

	setBMask(view0);

	initMergingWeightArrays(view0);

	initSegVp(view0);

}

void initMergingWeightArrays(view *view0) {

	deinitMergingWeightArrays(view0);

	if (view0->merge_weights == nullptr) {
		view0->merge_weights = new signed short[ getNB(view0) / 2 ]();
	}

	if (view0->merge_weights_float == nullptr){
		view0->merge_weights_float = new float[ getNB(view0) ]();
	}

}

void deinitMergingWeightArrays(view *view0) {

	if (view0->merge_weights != nullptr) {
		delete[](view0->merge_weights);
	}

	if (view0->merge_weights_float != nullptr) {
		delete[](view0->merge_weights_float);
	}

}

void getGeomWeight(view *view0) {

	/* we don't use LS weights but something derived on geometric distance in view array*/

	float sigma = view0->sigma;

	int MMM = 1 << view0->n_references;// pow(2, (view0)->n_references);

	double *thetas = new double[getNB(view0)]();

	bool *bmask = view0->bmask;

	for (int ii = 0; ii < MMM; ii++) {
		double sumw = 0;

		for (int ij = 0; ij < view0->n_references; ij++) {
			view *view1 = view0->color_reference_views[ij];
			double vdistance = (view0->x - view1->x)*(view0->x - view1->x) + (view0->y - view1->y)*(view0->y - view1->y);
			if (bmask[ii + ij * MMM]) {
				thetas[ii + ij * MMM] = exp( -(vdistance) / (2 * sigma*sigma) );
				//printf("exp(-(%3.3f) / (2 * %3.3f*%3.3f))=%3.3f\n", vdistance, sigma, sigma, thetas[ii + ij * MMM]);
				//printf("sigma %3.3f\n", sigma);
				//printf("thetas[ii+ij*MMM]=%3.3f\n", ii, ij, MMM, thetas[ii + ij * MMM]);
				sumw += thetas[ii + ij * MMM];
			}
		}
		//printf("sumw %3.3f\n", sumw);
		for (int ij = 0; ij < view0->n_references; ij++) {
			thetas[ii + ij * MMM] = thetas[ii + ij * MMM] / sumw;
		}
	}

	columnwise_collecting_of_thetas(view0, thetas);

	delete[](thetas);

}

void columnwise_collecting_of_thetas(view *view0, double *thetas) {
	/* columnwise collecting of thetas */

	int MMM = 1 << view0->n_references;

	signed short *LScoeffs = view0->merge_weights;
	int in = 0;
	for (int ik = 0; ik < view0->n_references; ik++) {
		for (int ij = 0; ij < MMM; ij++) {
			if (view0->bmask[ij + ik * MMM]) {
				view0->merge_weights[in] = (signed short)floor(thetas[ij + MMM * ik] * (signed short)(1 << BIT_DEPTH_MERGE) + 0.5);
				in++;
			}
		}
	}
}

void mergeMedian_N(view *view0, const int ncomponents) {

	/* merge color with median */

	unsigned short *AA2 = new unsigned short[view0->nr*view0->nc*ncomponents]();

#pragma omp parallel for
	for (int ii = 0; ii < view0->nr*view0->nc; ii++) {
		for (int icomp = 0; icomp < ncomponents; icomp++) {
			std::vector< unsigned short > vals;
			for (int ik = 0; ik < view0->n_references; ik++) {

				float *pf = view0->occlusion_masks[ik];
				unsigned short *ps = view0->warped_color_views[ik];

				if (*(pf + ii) > INIT_DISPARITY_VALUE) {

					vals.push_back(*(ps + ii + icomp*view0->nr*view0->nc));

				}

			}

			if (vals.size() > 0) {
				*(AA2 + ii + icomp*view0->nr*view0->nc) = getMedian(vals);
			}
		}
	}

	memcpy(view0->color, AA2, sizeof(unsigned short)*view0->nr*view0->nc * 3);
	delete[](AA2);

}

void mergeWarped_N(view *view0, const int ncomponents)
{

	/* merge color with prediction */

	int MMM = 1 << view0->n_references;// pow(2, (view0)->n_references);

	bool *bmask = (view0)->bmask;
	
	int uu = 0;

	for (int ii = 0; ii < MMM * (view0)->n_references; ii++){
		if (bmask[ii]){
			(view0)->merge_weights_float[ii] = ((float)(view0)->merge_weights[uu++]) / (float)(1 << BIT_DEPTH_MERGE);// pow(2, BIT_DEPTH_MERGE);
			//printf("%f\t", (view0)->merge_weights_float[ii]);
		}
		else{
			(view0)->merge_weights_float[ii] = 0.0;
		}
	}

	int nr = (view0)->nr;
	int nc = (view0)->nc;
	int n_views = (view0)->n_references;
	float *LSw = (view0)->merge_weights_float;

	unsigned short *seg_vp = (view0)->seg_vp;

	float *AA1 = new float[nr*nc*ncomponents]();
	unsigned short *AA2 = new unsigned short[nr*nc*ncomponents]();

	for (int ii = 0; ii < nr*nc; ii++){

		int ci = seg_vp[ii];

		for (int ik = 0; ik < n_views; ik++){
			unsigned short *ps = view0->warped_color_views[ik];
			for (int icomp = 0; icomp < 3; icomp++){
				AA1[ii + icomp*nr*nc] = AA1[ii + icomp*nr*nc] + LSw[ci + ik * MMM] * ((float)(*(ps + ii + icomp*nr*nc)));
			}
		}

		for (int icomp = 0; icomp < 3; icomp++){
			if (AA1[ii + icomp*nr*nc] < 0)
				AA1[ii + icomp*nr*nc] = 0;
			if (AA1[ii + icomp*nr*nc] > (1 << BIT_DEPTH) - 1)//(pow(2, BIT_DEPTH) - 1))
				AA1[ii + icomp*nr*nc] = (1 << BIT_DEPTH) - 1;// (pow(2, BIT_DEPTH) - 1);

			AA2[ii + icomp*nr*nc] = (unsigned short)(floor(AA1[ii + icomp*nr*nc]));

		}
	}

	memcpy((view0)->color, AA2, sizeof(unsigned short)*nr*nc*ncomponents);

	delete[](AA1);
	delete[](AA2);

}
