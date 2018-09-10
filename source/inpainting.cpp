#include "inpainting.hh"
#include "medianfilter.hh"

void holefilling(unsigned short *pshort, const int ncomps, const int nr, const int nc, const unsigned short maskval)
{
	std::vector< unsigned int > holes = getHoles(pshort, nr, nc, maskval, ncomps);

	int iter = 0;

	while (holes.size() > 0) {
		holefilling_alg(holes, pshort, ncomps, nr, nc, maskval);
		holes = getHoles(pshort, nr, nc, maskval, ncomps);
		printf("holefilling iteration:\t%i\t n:%i\n", iter++, static_cast<int>(holes.size()));
	}
}

std::vector<unsigned int> getHoles(unsigned short *pshort,const int nr, const int nc, const unsigned short maskval, const int ncomps) {

	std::vector<unsigned int> holes;

	for (int ii = 0; ii < nr*nc; ii++) {

		bool is_hole = (pshort[ii] == maskval);
		if (ncomps == 3) {
			is_hole = is_hole && (pshort[ii + nr*nc] == maskval) && (pshort[ii + 2 * nr*nc] == maskval);
		}

		if (is_hole) {
			holes.push_back(ii);
		}
	}

	return holes;
}

void holefilling_alg(std::vector< unsigned int > holes,unsigned short *pshort, const int ncomps, const int nr, const int nc, const unsigned short maskval)
{

	std::vector<unsigned short> neighbours;

	int dsz = 1;

	for (int jj = 0; jj < holes.size(); jj++) {

		int ii = holes.at(jj);

		int y, x;
		y = ii%nr;
		x = (ii / nr);

		for (int icomp = 0; icomp < ncomps; icomp++) {
			neighbours.clear();
			for (int dy = -dsz; dy <= dsz; dy++) {
				for (int dx = -dsz; dx <= dsz; dx++) {
					if (!(dy == 0 && dx == 0)) {
						if ((y + dy) >= 0 && (y + dy) < nr && (x + dx) >= 0 && (x + dx) < nc) {
							bool is_hole_again = (pshort[y + dy + (x + dx)*nr] == maskval);
							if (ncomps == 3) {
								is_hole_again = is_hole_again && (pshort[y + dy + (x + dx)*nr + nr*nc] == maskval) && (pshort[y + dy + (x + dx)*nr + 2 * nr*nc] == maskval);
							}
							if (!is_hole_again) {
								neighbours.push_back(pshort[y + dy + (x + dx)*nr + icomp*nr*nc]);
							}
						}
					}
				}
			}
			if (neighbours.size() > 0) {
				pshort[ii + icomp*nr*nc] = getMedian(neighbours);
			}
		}
	}
}
