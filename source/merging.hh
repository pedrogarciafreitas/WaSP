#ifndef MERGING_HH
#define MERGING_HH

#include "view.hh"

void initSegVp(view *view0);
void initViewW(view *view0);
void initMergingWeightArrays(view *view0);
void deinitMergingWeightArrays(view *view0);

void getViewMergingLSWeights_N(view *view0);

void getGeomWeight(view *view0);
void mergeMedian_N(view *view0, const int ncomponents);
void mergeWarped_N(view *view0, const int ncomponents);

void columnwise_collecting_of_thetas(view *view0, double *thetas);

#endif