#ifndef SPARSEFILTER_HH
#define SPARSEFILTER_HH

#include <vector>

#include "view.hh"

void sortSparseFilter(const int Ms, unsigned char *sparse_mask, int32_t *sparse_weights);

int checkOutOfBounds(const int RR, const int nr, const int dy_0);
int checkBoundaryPixels(const int32_t *label_im, int offset, const int region_indx, const int NNt, const int ir, const int ic, const int dy, const int dx, const int nr, const int nc);

unsigned short *applySparseFilterForOneRegion(view *view0, region_sparse_filter reg_sp);
region_sparse_filter getSparseFilterForOneRegion(view *view0, const int region_indx, const int Ms, const int NNt);

bool sortinrev(const std::pair<int, int> &a, const std::pair<int, int> &b);
std::vector<std::pair<int, int>> sortRegionsBySize(view *view0);

void applyGlobalSparseFilter(view *view0);
void getGlobalSparseFilter(view *view0);

#endif