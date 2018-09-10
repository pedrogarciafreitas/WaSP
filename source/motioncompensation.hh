#ifndef MOTIONCOMPENSATION_HH
#define MOTIONCOMPENSATION_HH

#define MOTION_VECTORS true
#define MV_REG_MIN_SIZE 64
#define MV_MAX_REGS 1024

#include "view.hh"

int findMVIndex(view *view0, view* view1);
void sortRegionsBySize(view *view0);
void getMotionVectorsView0_to_View1(view *view0, view *view1);

#endif