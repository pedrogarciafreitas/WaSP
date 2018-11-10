#ifndef CONNECTEDCOMPONENTS_HH
#define CONNECTEDCOMPONENTS_HH

#include <stdint.h>
#include <cstdio>

#define MAXREG 250000

int32_t *get_labels(int *img, const int nr, const int nc, int &nregions, int *&reg_histogram);
void get_gSR_matrix(int **gr, int **gSR, const int nr, const int nc);
int RegTreatSegment(int istart, int iend, int newregindex, int* sameregs,
	int nsame, int ir, int* minlabels, int* usamereg, int **Cimg);
int MarkRegions(int **gSR, int **Cimg, const int nr, const int nc);

#endif