#ifndef SPARSEFILTER_HH
#define SPARSEFILTER_HH

#define REGION_SPARSE_ON true

#include "view.hh"

void sortSparseFilter(const int Ms, unsigned char *sparse_mask, int32_t *sparse_weights);

void checkOutOfBounds(const int RR, const int CC, const int nr, const int nc, const int dy_0, const int dx_0, int &dy, int &dx);
int checkBoundaryPixels(int *label_im, int offset, const int ijk, const int NNt, const int ir, const int ic, const int dy, const int dx, const int nr, const int nc);

void applyRegionSparseFilters(view *view0);
std::vector< std::pair< int, int> > validateRegionSparseFilters(view *view0, unsigned short *original_color_view);
void getRegionSparseFilter(view *view0, unsigned short *original_color_view);

void applyGlobalSparseFilter(view *view0);
void getGlobalSparseFilter(view *view0, unsigned short *original_color_view);

#endif