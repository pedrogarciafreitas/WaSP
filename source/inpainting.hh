#ifndef INPAINTING_HH
#define INPAINTING_HH

#include <vector>

std::vector<unsigned int> getHoles(unsigned short *pshort, const int nr, const int nc, const unsigned short maskval, const int ncomps);
void holefilling(unsigned short *pshort, const int ncomps, const int nr, const int nc, const unsigned short maskval);
void holefilling_alg(std::vector< unsigned int > holes, unsigned short *pshort, const int ncomps, const int nr, const int nc, const unsigned short maskval);

#endif