#include "codestream.hh"
#include "view.hh"
#include "minconf.hh"

#include <iostream>

void readSparseFromBitstream(int &n_bytes_prediction, FILE *input_LF, int &NNt, int &Ms, unsigned char *&sparse_mask, int32_t *&sparse_weights, const bool regionsparse){

	if ( !regionsparse ) {

		unsigned char tmpNNt = 0;
		unsigned char tmpMs = 0;

		n_bytes_prediction += (int)fread(&tmpNNt, sizeof(unsigned char), 1, input_LF) * sizeof(unsigned char);
		n_bytes_prediction += (int)fread(&tmpMs, sizeof(unsigned char), 1, input_LF) * sizeof(unsigned char);

		NNt = (int)tmpNNt;
		Ms = (int)tmpMs;

	}

	sparse_weights = new int32_t[Ms]();
	sparse_mask = new unsigned char[Ms]();

	n_bytes_prediction += (int)fread(sparse_weights, sizeof(int32_t), Ms, input_LF) * sizeof(int32_t);

	int32_t sparse_mask_binary[2];

	n_bytes_prediction += (int)fread(&sparse_mask_binary[0], sizeof(int32_t), 1, input_LF) * sizeof(int32_t);
	n_bytes_prediction += (int)fread(&sparse_mask_binary[1], sizeof(int32_t), 1, input_LF) * sizeof(int32_t);

	unpackSparseMask(sparse_mask_binary, sparse_mask);
}

void writeSparseToBitstream(int &n_bytes_prediction,FILE *output_LF_file, const int Ms, unsigned char *sparse_mask, const int NNt, int32_t *sparse_weights, const bool regionsparse) {

	int32_t sparse_mask_binary[2] = { 0,0 };

	packSparseMask(Ms, sparse_mask, sparse_mask_binary);

	if (!regionsparse) {

		unsigned char tmpNNt = (unsigned char)NNt;
		unsigned char tmpMs = (unsigned char)Ms;

		n_bytes_prediction += (int)fwrite(&tmpNNt, sizeof(unsigned char), 1, output_LF_file) * sizeof(unsigned char);
		n_bytes_prediction += (int)fwrite(&tmpMs, sizeof(unsigned char), 1, output_LF_file) * sizeof(unsigned char);
	}

	n_bytes_prediction += (int)fwrite(sparse_weights, sizeof(int32_t), Ms, output_LF_file) * sizeof(int32_t);
	n_bytes_prediction += (int)fwrite(&sparse_mask_binary[0], sizeof(int32_t), 1, output_LF_file) * sizeof(int32_t);
	n_bytes_prediction += (int)fwrite(&sparse_mask_binary[1], sizeof(int32_t), 1, output_LF_file) * sizeof(int32_t);

}

void packSparseMask(const int Ms, const unsigned char *sparse_mask, int32_t *sparse_mask_binary) {

	for (int ij = 0; ij < Ms; ij++) {

		if (sparse_mask[ij] > 0) {

			if (sparse_mask[ij] <= 30) {
				sparse_mask_binary[0] = sparse_mask_binary[0] | 1 << ((int32_t)sparse_mask[ij]); // note: regressor indexing starts from 1
			}
			else {
				sparse_mask_binary[1] = sparse_mask_binary[1] | 1 << ((int32_t)sparse_mask[ij] - 31);
			}

		}

	}
}

void unpackSparseMask(int32_t *sparse_mask_binary, unsigned char *sparse_mask) {
	int ik = 0;

	for (int ij = 0; ij < 64; ij++) {
		if (ij <= 30) {
			if ((sparse_mask_binary[0] & (1 << ij))>0) {
				sparse_mask[ik] = ij; ik++;
			}
		}
		else {
			if ((sparse_mask_binary[1] & (1 << (ij - 31)))>0) {
				sparse_mask[ik] = ij; ik++;
			}
		}
	}
}

void viewHeaderToCodestream(int &n_bytes_prediction, view *SAI, FILE *output_LF_file, const int yuv_transform_s) {

	minimal_config mconf = makeMinimalConfig(SAI);

	//printf("size of minimal_config %i bytes\n", (int)sizeof(minimal_config));

	n_bytes_prediction += (int)fwrite(&mconf, sizeof(minimal_config), 1, output_LF_file) * sizeof(minimal_config);

	/* lets see what else needs to be written to bitstream */

	if (SAI->has_x_displacement) {
		n_bytes_prediction += (int)fwrite(&SAI->x, sizeof(float), 1, output_LF_file) * sizeof(float);
	}

	if (SAI->has_y_displacement) {
		n_bytes_prediction += (int)fwrite(&SAI->y, sizeof(float), 1, output_LF_file) * sizeof(float);
	}

	if (SAI->has_color_references) {
		unsigned char tmpNREF = (unsigned char)SAI->n_references;
		n_bytes_prediction += (int)fwrite(&tmpNREF, sizeof(unsigned char), 1, output_LF_file) * sizeof(unsigned char);
		for (int ij = 0; ij < SAI->n_references; ij++) {
			unsigned short nid = (unsigned short) *(SAI->references + ij);
			n_bytes_prediction += (int)fwrite(&nid, sizeof(unsigned short), 1, output_LF_file) * sizeof(unsigned short);
		}
	}

	if (SAI->has_depth_references) {
		unsigned char tmpNDREF = (unsigned char)SAI->n_depth_references;
		n_bytes_prediction += (int)fwrite(&tmpNDREF, sizeof(unsigned char), 1, output_LF_file) * sizeof(unsigned char);
		for (int ij = 0; ij < SAI->n_depth_references; ij++) {
			unsigned short nid = (unsigned short) *(SAI->depth_references + ij);
			n_bytes_prediction += (int)fwrite(&nid, sizeof(unsigned short), 1, output_LF_file) * sizeof(unsigned short);
		}
	}

	if (!SAI->use_median) {
		if (SAI->stdd < 0.001) {
			if (SAI->has_color_references) {
				/* use LS merging weights */
				n_bytes_prediction += (int)fwrite(SAI->merge_weights, sizeof(signed short), SAI->NB / 2, output_LF_file) * sizeof(signed short);
			}
		}
		else {
			/* use standard deviation */
			n_bytes_prediction += (int)fwrite(&SAI->stdd, sizeof(float), 1, output_LF_file) * sizeof(float);
		}
	}

	if ( SAI->use_global_sparse ) {

		writeSparseToBitstream(n_bytes_prediction, output_LF_file, SAI->Ms, SAI->sparse_mask, SAI->NNt, SAI->sparse_weights, false);

	}

	if (SAI->use_region_sparse) {

		unsigned int nregions_sparse = static_cast<unsigned int>( SAI->region_Regr.size() );

		n_bytes_prediction += (int)fwrite(&nregions_sparse, sizeof(unsigned int), 1, output_LF_file) * sizeof(unsigned int);

		for (unsigned int ii = 0; ii < nregions_sparse; ii++) {

			n_bytes_prediction += (int)fwrite(&SAI->valid_regions_ir[ii], sizeof(unsigned int), 1, output_LF_file) * sizeof(unsigned int);

			writeSparseToBitstream(n_bytes_prediction, output_LF_file, SAI->Ms, SAI->region_Regr.at(ii).data(), SAI->NNt, SAI->region_Theta.at(ii).data(), true);

		}

	}

	if (SAI->use_motion_vectors) {

		unsigned short n_views = static_cast<unsigned short>( SAI->mv_views.size() );

		n_bytes_prediction += (int)fwrite(&n_views, sizeof(unsigned short), 1, output_LF_file) * sizeof(unsigned short);

		for (int jj = 0; jj < n_views; jj++) {

			unsigned int n_regions = static_cast<unsigned int>( SAI->mv_views.at(jj).second.size() );

			n_bytes_prediction += (int)fwrite(&n_regions, sizeof(unsigned int), 1, output_LF_file) * sizeof(unsigned int);

			n_bytes_prediction += (int)fwrite(&SAI->mv_views.at(jj).first, sizeof(int), 1, output_LF_file) * sizeof(int);

			for (unsigned int ik = 0; ik < n_regions; ik++) {

				n_bytes_prediction += (int)fwrite(&SAI->mv_views.at(jj).second.at(ik), sizeof(MV_REGION), 1, output_LF_file) * sizeof(MV_REGION);

			}

		}

	}

	return;

}

void codestreamToViewHeader( int &n_bytes_prediction, view *SAI, FILE *input_LF, minimal_config &mconf ) {

	n_bytes_prediction += (int)fread(&mconf, sizeof(minimal_config), 1, input_LF)* sizeof(minimal_config);

	//printf("size of minimal_config %i bytes\n", (int)sizeof(minimal_config));

	setup_form_minimal_config(&mconf, SAI);

	if (SAI->has_x_displacement) {
		n_bytes_prediction += (int)fread(&SAI->x, sizeof(float), 1, input_LF) * sizeof(float);
	}

	if (SAI->has_y_displacement) {
		n_bytes_prediction += (int)fread(&SAI->y, sizeof(float), 1, input_LF) * sizeof(float);
	}


	if (SAI->has_color_references) {

		unsigned char tmpNREF = 0;

		n_bytes_prediction += (int)fread(&tmpNREF, sizeof(unsigned char), 1, input_LF) * sizeof(unsigned char);

		SAI->n_references = tmpNREF;

		SAI->references = new int[SAI->n_references]();
		for (int ij = 0; ij < SAI->n_references; ij++) {
			unsigned short nid;
			n_bytes_prediction += (int)fread(&nid, sizeof(unsigned short), 1, input_LF)* sizeof(unsigned short);
			*(SAI->references + ij) = (int)nid;
		}
	}

	if (SAI->has_depth_references) {

		unsigned char tmpNDREF = 0;

		n_bytes_prediction += (int)fread(&tmpNDREF, sizeof(unsigned char), 1, input_LF) * sizeof(unsigned char);

		SAI->n_depth_references = tmpNDREF;


		SAI->depth_references = new int[SAI->n_depth_references]();
		for (int ij = 0; ij < SAI->n_depth_references; ij++) {
			unsigned short nid;
			n_bytes_prediction += (int)fread(&nid, sizeof(unsigned short), 1, input_LF)* sizeof(unsigned short);
			*(SAI->depth_references + ij) = (int)nid;
		}
	}

	SAI->NB = (1 << SAI->n_references)*SAI->n_references;

	if (!SAI->use_median) {
		if (SAI->stdd < 0.001) {
			if (SAI->n_references > 0) {
				SAI->merge_weights = new signed short[SAI->NB / 2]();
				n_bytes_prediction += (int)fread(SAI->merge_weights, sizeof(signed short), SAI->NB / 2, input_LF) * sizeof(signed short);
			}
		}
		else {
			n_bytes_prediction += (int)fread(&SAI->stdd, sizeof(float), 1, input_LF) * sizeof(float);
		}
	}

	if ( SAI->use_global_sparse ) {

		readSparseFromBitstream(n_bytes_prediction, input_LF, SAI->NNt, SAI->Ms, SAI->sparse_mask, SAI->sparse_weights, false);
		
	}

	if (SAI->use_region_sparse) {

		unsigned int nregions_sparse = 0;
		n_bytes_prediction += (int)fread(&nregions_sparse, sizeof(unsigned int), 1, input_LF) * sizeof(unsigned int);

		for (unsigned int ii = 0; ii < nregions_sparse; ii++) {

			unsigned char *tmp_mask;
			int32_t *tmp_weights;

			int tmpIR;

			n_bytes_prediction += (int)fread(&tmpIR, sizeof(unsigned int), 1, input_LF) * sizeof(unsigned int);

			readSparseFromBitstream(n_bytes_prediction, input_LF, SAI->NNt, SAI->Ms, tmp_mask, tmp_weights, true);

			std::vector< unsigned char > vec_mask;
			std::vector< int32_t > vec_weights;

			for (int jj = 0; jj < SAI->Ms; jj++) {
				vec_mask.push_back(tmp_mask[jj]);
				vec_weights.push_back(tmp_weights[jj]);
			}

			SAI->region_Regr.push_back(vec_mask);
			SAI->region_Theta.push_back(vec_weights);
			SAI->valid_regions_ir.push_back(tmpIR);

			delete[](tmp_weights);
			delete[](tmp_mask);

		}

	}

	if (SAI->use_motion_vectors) {

		unsigned short n_views;

		n_bytes_prediction += (int)fread(&n_views, sizeof(unsigned short), 1, input_LF) * sizeof(unsigned short);

		/* Following loop reads all motion vector information to be used in forward warping this view to other views */
		
		for (int jj = 0; jj < n_views; jj++) {

			std::pair< unsigned short, std::vector<MV_REGION> > tmp_mv_view; /* For each view we store its index and a vector of regions with reg_id,dy,dx */

			unsigned int n_regions = 0;

			n_bytes_prediction += (int)fread(&n_regions, sizeof(unsigned int), 1, input_LF) * sizeof(unsigned int);

			n_bytes_prediction += (int)fread(&tmp_mv_view.first, sizeof(int), 1, input_LF) * sizeof(int);

			for (unsigned int ik = 0; ik < n_regions; ik++) {

				MV_REGION tmp_region;

				n_bytes_prediction += (int)fread(&tmp_region, sizeof(MV_REGION), 1, input_LF) * sizeof(MV_REGION);

				tmp_mv_view.second.push_back(tmp_region);

			}

			SAI->mv_views.push_back(tmp_mv_view);

		}

	}

	return;

}