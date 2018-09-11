#include "connected_components.hh"

int *get_labels(int *img, const int nr, const int nc, int &nregions, int *&reg_histogram) {

	int NRgr = 2 * nr + 4;
	int NCgr = 2 * nc + 4;

	int **gr = new int*[NRgr]();
	int **gSR = new int*[NRgr]();

	for (int ii = 0; ii < NRgr; ii++) {
		gr[ii] = new int[NCgr]();
		gSR[ii] = new int[NCgr]();
	}

	for (int ii = 0; ii < nr; ii++) {
		for (int jj = 0; jj < nc; jj++) {
			gr[3 + 2 * ii][3 + 2 * jj] = img[ii + nr*jj];
		}
	}

	get_gSR_matrix(gr, gSR, nr, nc);

	int **Cimg = new int*[nr]();
	for (int ii = 0; ii < nr; ii++) {
		Cimg[ii] = new int[nc]();
	}

	int tmpnreg = MarkRegions(gSR, Cimg, nr, nc);

	int *label_im = new int[nr*nc]();

	int ee = 0;
	for (int c = 0; c < nc; c++) {
		for (int r = 0; r < nr; r++) {
			*(label_im + ee) = Cimg[r][c];
			ee++;
		}
	}

	for (int ii = 0; ii < nr; ii++) {
		delete[](Cimg[ii]);
	}
	delete[](Cimg);

	for (int ii = 0; ii < NRgr; ii++) {
		delete[](gr[ii]);
		delete[](gSR[ii]);
	}
	delete[](gr);
	delete[](gSR);

	nregions = INT32_MIN;

	reg_histogram = new int[MAXREG]();

#pragma omp parallel for
	for (int ii = 0; ii < nr*nc; ii++) {
		nregions = *(label_im + ii) > nregions ? *(label_im + ii)+1 : nregions;

		if (*(label_im + ii) < MAXREG) {
			reg_histogram[*(label_im + ii)]++;
		}

	}

	printf("Number of regions:\t%i\n", nregions);

	return label_im;

}

void get_gSR_matrix(int **gr, int **gSR, const int nr, const int nc) {

#pragma omp parallel for
	for (int ir = 0; ir < nr; ir++) {
		for (int ic = 0; ic < nc; ic++)
		{
			if (ic > 0)
				if (gr[2 * ir + 3][2 * (ic - 1) + 3] != gr[2 * ir + 3][2 * ic + 3])
				{
					gr[2 * ir + 3][2 * ic + 2] = 1;
					gSR[ir][ic] = 1;
				}

			if (ir > 0)
				if (gr[2 * (ir - 1) + 3][2 * ic + 3] != gr[2 * ir + 3][2 * ic + 3])
				{
					gr[2 * ir + 2][2 * ic + 3] = 1;
					gSR[ir][ic] += 2;
				}
		}
	}

	for (int ii = 0; ii < nr; ii++) {
		for (int jj = 0; jj < nc; jj++) {
			gr[3 + 2 * ii][3 + 2 * jj] = 0;
		}
	}


}

int RegTreatSegment(int istart, int iend, int newregindex, int* sameregs, 
	int nsame, int ir, int* minlabels, int* usamereg, int **Cimg)
{
	int i, j, unique;
	int nusr, minusame;

	if (nsame == 0)
	{	// no need to update anything, just mark a new region
		newregindex = newregindex + 1;
		minlabels[newregindex] = newregindex;
		for (i = istart; i <= iend; i++)
			Cimg[ir][i] = newregindex;
		return newregindex;
	}

	// We have some neighbors: Find the distinct neighbours in sameregs
	nusr = 1;
	usamereg[0] = sameregs[0];
	for (i = 1; i < nsame; i++)
	{
		unique = 1;
		for (j = 0; j < nusr; j++)
			if (usamereg[j] == sameregs[i])
			{
				unique = 0;
				break;
			}
		if (unique == 1)
		{
			nusr = nusr + 1;
			usamereg[nusr - 1] = sameregs[i];
		}
	}

	// find the minimum of the representatives of sameregs
	minusame = minlabels[usamereg[0]];
	for (j = 0; j < nusr; j++) {
		if (minlabels[usamereg[j]] < minusame)
		{
			minusame = minlabels[usamereg[j]];
		}
	}

	// Update for all the minimum label
	for (j = 0; j < nusr; j++) {
		minlabels[usamereg[j]] = minusame;
	}

	for (i = istart; i <= iend; i++) {
		Cimg[ir][i] = minusame;
	}

	return newregindex;
}

int MarkRegions(int **gSR, int **Cimg, const int nr, const int nc)
{
	int i, ic, ir;
	int* minlabels;
	int* finallabels;
	int nreg, nsame, newregindex;
	int* sameregs;
	int* usamereg;
	int istart, iend;
	int change, labeli, labeli1;


	minlabels = new int[nr*nc + 10]();
	for (i = 0; i < nr*nc + 10; i++) {
		minlabels[i] = i;
	}

	finallabels = new int[nr*nc + 10]();
	sameregs = new int[nc + 10]();
	usamereg = new int[nc + 15]();   // unique values in sameregs

	newregindex = 0;

	for (ir = 0; ir<nr; ir++)
	{
		// Inspect the current line only once
		// Prepare for the first segment
		istart = 0;
		for (i = 0; i < nc + 2; i++) {
			sameregs[i] = 0;
		}
		nsame = 0;
		for (ic = 0; ic<nc; ic++)
		{
			if (gSR[ir][ic] % 2 == 1)
			{
				// There is a vertical edge at the left of ic
				// stop here for a while, treat the segment istart:(ic-1)
				iend = ic - 1;
				newregindex = RegTreatSegment(istart, iend, newregindex, sameregs, nsame, ir, minlabels, usamereg, Cimg);
				// now start a new segment
				istart = ic;
				for (i = 1; i < nc + 2; i++) {
					sameregs[i] = 0;
				}
				nsame = 0;
			}
			// if there is no hor. edge above a pixel, the pixel will have the same region index as the above pixel
			if (gSR[ir][ic] < 2)
			{
				if (ir>0)
				{

					sameregs[nsame] = Cimg[ir - 1][ic];
					nsame = nsame + 1;
				}
			}

			if (ic == (nc - 1))  // end of row
			{
				iend = (nc - 1);
				newregindex = RegTreatSegment(istart, iend, newregindex, sameregs, nsame, ir, minlabels, usamereg, Cimg);
			}
		}
	}

	// Find the minimum label representative for each label

	change = 1;
	while (change == 1)
	{
		change = 0;
		for (i = 0; i<newregindex; i++)
		{
			labeli = minlabels[i];
			labeli1 = labeli;
			while (1)
			{
				if (minlabels[labeli] < labeli)
				{
					labeli = minlabels[labeli];
				}
				else
				{
					break;
				}
			}
			if (labeli1 != labeli) {
				change = 1;
			}
			minlabels[i] = labeli;
		}
	}

	nreg = 0;
	for (i = 1; i <= newregindex; i++) {
		if (minlabels[i] == i) {
			nreg = nreg + 1;
			finallabels[i] = nreg;
		}
	}

	for (ir = 0; ir < nr; ir++) {
		for (ic = 0; ic < nc; ic++) {
			labeli = minlabels[Cimg[ir][ic]];
			Cimg[ir][ic] = finallabels[labeli];
		}
	}

	delete[](minlabels);
	delete[](finallabels);
	delete[](sameregs);
	delete[](usamereg);

	return newregindex;
}