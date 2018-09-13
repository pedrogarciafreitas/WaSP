#ifndef MOTIONCOMPENSATION_HH
#define MOTIONCOMPENSATION_HH

#define MOTION_VECTORS false
#define MV_REG_MIN_SIZE 1024
#define MV_MAX_REGS 1024
#define MIN_RADIUS 6 

#include "view.hh"

template <class T>
T getMode(std::vector<T> vector)
{
	T max = *std::max_element(vector.begin(), vector.end());
	T min = *std::min_element(vector.begin(), vector.end());

	int nrange = static_cast<int>(max) - static_cast<int>(min);

	int *hh = new int[nrange+1]();

	for (int jj = 0; jj < vector.size(); jj++) {
		for (int ii = 0; ii < nrange; ii++) {

			if ((vector.at(jj) - min) == ii) {
				hh[ii]++;
			}

		}
	}
	
	T mode = 0;
	int mode_n = 0;

	for (int ii = 0; ii < nrange; ii++) {
		if (hh[ii] > mode_n) {
			mode = static_cast<T>( ii + min );
			mode_n = hh[ii];
		}
	}

	return mode;
}

void readLabelIm(view *view0);
int findMVIndex(view *view0, view* view1);
void sortRegionsBySize(view *view0);
void getMotionVectorsView0_to_View1(view *view0, view *view1);

#endif