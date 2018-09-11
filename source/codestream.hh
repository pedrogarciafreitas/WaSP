#ifndef CODESTREAM_HH
#define CODESTREAM_HH

#include "view.hh"
#include "minconf.hh"
#include <iostream>

void readSparseFromBitstream(int &n_bytes_prediction, FILE *input_LF, int &NNt, int &Ms, unsigned char *&sparse_mask, int32_t *&sparse_weights, const bool regionsparse);
void writeSparseToBitstream(int &n_bytes_prediction, FILE *output_LF_file, const int Ms, unsigned char *sparse_mask, const int NNt, int32_t *sparse_weights, const bool regionsparse);

void packSparseMask(const int Ms, const unsigned char *sparse_mask, int32_t *sparse_mask_binary);
void unpackSparseMask(int32_t *sparse_mask_binary, unsigned char *sparse_mask);

void viewHeaderToCodestream(int &n_bytes_prediction, view *SAI, FILE *output_LF_file, const int yuv_transform_s);
void codestreamToViewHeader(int &n_bytes_prediction, view *SAI, FILE *input_LF, minimal_config &mconf);

#endif