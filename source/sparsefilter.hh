#ifndef SPARSEFILTER_HH
#define SPARSEFILTER_HH

#include "view.hh"

void checkOutOfBounds(const int RR, const int CC, const int nr, const int nc, const int dy_0, const int dx_0, int &dy, int &dx);
int checkBoundaryPixels(int *label_im, int offset, const int ijk, const int NNt, const int ir, const int ic, const int dy, const int dx, const int nr, const int nc);

void applyRegionSparseFilter(view *view0);
void getRegionSparseFilter(view *view0, unsigned short *original_color_view);

void applyGlobalSparseFilter(view *view0);
void getGlobalSparseFilter(view *view0, unsigned short *original_color_view);

#endif