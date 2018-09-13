#ifndef MOTIONCOMPENSATION_HH
#define MOTIONCOMPENSATION_HH

#define MOTION_VECTORS false
#define MV_REG_MIN_SIZE 256
#define MV_MAX_REGS 1024
#define MIN_RADIUS 3 

#include "view.hh"

void readLabelIm(view *view0);
int findMVIndex(view *view0, view* view1);
void sortRegionsBySize(view *view0);
void getMotionVectorsView0_to_View1(view *view0, view *view1);

#endif