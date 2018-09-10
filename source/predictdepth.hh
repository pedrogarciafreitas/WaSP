#ifndef PREDICTDEPTH_HH
#define PREDICTDEPTH_HH

#include "view.hh"

#define MEDFILT_DEPTH true
#define MEDFILT_DEPTH_SZ 3

void predictDepth(view* SAI, view *LF);


#endif