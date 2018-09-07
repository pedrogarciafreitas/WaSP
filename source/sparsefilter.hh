#ifndef SPARSEFILTER_HH
#define SPARSEFILTER_HH

#include "view.hh"

void applyRegionSparseFilter(view *view0);
void getRegionSparseFilter(view *view0, unsigned short *original_color_view);

void applyGlobalSparseFilter(view *view0);
void getGlobalSparseFilter(view *view0, unsigned short *original_color_view);

#endif