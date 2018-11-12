#include <cstdio>
#include <cstdint>
#include <ctime>
#include <cstring>
#include <iostream>
#include <vector>
#include <algorithm>

#include "bitdepth.hh"
#include "merging.hh"
#include "minconf.hh"
#include "sparsefilter.hh"
#include "ycbcr.hh"
#include "view.hh"
#include "warping.hh"
#include "residualjp2.hh"
#include "ppm.hh"
#include "fileaux.hh"
#include "medianfilter.hh"
#include "psnr.hh"
#include "inpainting.hh"
#include "predictdepth.hh"
#include "codestream.hh"
#include "kmeans.hh"
#include "connected_components.hh"
#include "configLoader.hh"

#define USE_difftest_ng false

#define STD_SEARCH_LOW 10
#define STD_SEARCH_HIGH 250
#define STD_SEARCH_STEP 10

#ifndef FLT_MAX
#define FLT_MAX          3.402823466e+38F        // max value
#endif

int main(int argc, char** argv) {

	const char *input_dir = argv[1];
	const char *output_dir = argv[2];
	const char *kakadu_dir = argv[3];
	const char *config_file = argv[4];

	char kdu_compress_path[1024];
	char kdu_expand_path[1024];

	sprintf(kdu_compress_path, "%s%s", kakadu_dir, "/kdu_compress");
	sprintf(kdu_expand_path, "%s%s", kakadu_dir, "/kdu_expand");

	global_parameters *global_params = load_global_parameters(config_file);
	view *LF = load_config_and_init_LF(config_file);

	for (int ii = 0; ii < global_params->n_views_total; ii++) {
		setViewFilePaths(LF + ii, output_dir, input_dir);
	}

	/* to get effiency from multiple JP2 files, we remove parts of the files
	which are repetative over all files. For this we have a minimalistic 
	dictionary method. */
	std::vector<std::vector<unsigned char>> JP2_dict;

	/* get sub aperture image dimensions from first color image */
	int ncomp1;
	unsigned short *original_color_view = nullptr;
	aux_read16PGMPPM(LF->path_input_ppm, LF->nc, LF->nr, ncomp1, original_color_view);
	delete[](original_color_view);

	for (int ii = 0; ii < global_params->n_views_total; ii++) {
		(LF+ii)->nr = LF->nr;
		(LF+ii)->nc = LF->nc;
	}
	
	/* predict and get residual for INVERSE DEPTH at all views */
	for (int ii = 0; ii < global_params->n_views_total; ii++) {

		view *SAI = LF + ii;

		printf("Predicting inverse depth view %03d_%03d\n", SAI->c, SAI->r);

		unsigned short *original_depth_view = nullptr;

		int ncomp1;

		if (SAI->residual_rate_depth > 0) {
			SAI->depth_file_exist = aux_read16PGMPPM(SAI->path_input_pgm, SAI->nc, SAI->nr, ncomp1, original_depth_view);
		}

		SAI->depth = new unsigned short[SAI->nr*SAI->nc]();

		predictDepth(SAI, LF);

		if (SAI->residual_rate_depth > 0 && SAI->depth_file_exist) { /* residual depth if needed */

			sprintf(SAI->pgm_residual_depth_path, "%s%c%03d_%03d%s", SAI->output_dir, '/', SAI->c, SAI->r, "_depth_residual.pgm");

			sprintf(SAI->jp2_residual_depth_path_jp2, "%s%c%03d_%03d%s", SAI->output_dir, '/', SAI->c, SAI->r, "_depth_residual.jp2");

			encodeResidualJP2(SAI->nr, SAI->nc, original_depth_view, SAI->depth, SAI->pgm_residual_depth_path,
				kdu_compress_path, SAI->jp2_residual_depth_path_jp2, SAI->residual_rate_depth, 1, 0, 1);

			decodeResidualJP2(SAI->depth, kdu_expand_path, SAI->jp2_residual_depth_path_jp2, SAI->pgm_residual_depth_path, ncomp1, 0, (1 << 16) - 1, 1);

			SAI->has_depth_residual = true;
		}

		/* median filter depth */
		if (MEDFILT_DEPTH) {
			unsigned short *tmp_depth = new unsigned short[SAI->nr*SAI->nc]();
			int startt = clock();
			medfilt2D(SAI->depth, tmp_depth, 3, SAI->nr, SAI->nc);
			std::cout << "time elapsed in depth median filtering\t" << (float)((int)clock() - startt) / CLOCKS_PER_SEC << "\n";
			memcpy(SAI->depth, tmp_depth, sizeof(unsigned short)*SAI->nr*SAI->nc);
			delete[](tmp_depth);
		}

		aux_write16PGMPPM(SAI->path_out_pgm, SAI->nc, SAI->nr, 1, SAI->depth);

		delete[](original_depth_view);

		unloadInverseDepth(SAI);

	}

	/* get segmentation from inverse depth for all views*/
	for (int ii = 0; ii < global_params->n_views_total; ii++) {

		view *SAI = LF + ii;

		loadInverseDepth(SAI);

		printf("Obtaining segmentation based on inverse depth at view %03d_%03d\t", SAI->c, SAI->r);

		int *tmp_d = new int[SAI->nr*SAI->nc]();
		for (int iii = 0; iii < SAI->nr*SAI->nc; iii++) {
			*(tmp_d + iii) = static_cast<int>(*(SAI->depth + iii));
		}

		getKmeansQuantized(K_MEANS_CLUSTERS, tmp_d, SAI->nr*SAI->nc, K_MEANS_ITERATIONS); /*inplace assignment to tmp_d*/																				  

		int nregions = 0;
		int *reg_histogram = 0;
		int32_t *label_im = get_labels(tmp_d, SAI->nr, SAI->nc, nregions, reg_histogram);

		SAI->label_im = label_im;
		SAI->nregions = nregions;
		SAI->reg_histogram = reg_histogram;

		unsigned short *labels = new unsigned short[SAI->nr*SAI->nc]();
		for (int ii = 0; ii < SAI->nr*SAI->nc; ii++) {
			*(labels + ii) = static_cast<unsigned short>(  *(label_im + ii) );
		}

		char labels_file[1024];
		sprintf(labels_file, "%s%c%03d_%03d%s", SAI->output_dir, '/', SAI->c, SAI->r, "_labels.pgm");
		aux_write16PGMPPM(labels_file, SAI->nc, SAI->nr, 1, labels);

		char tmp_d_file[1024];
		sprintf(tmp_d_file, "%s%c%03d_%03d%s", SAI->output_dir, '/', SAI->c, SAI->r, "_kmeans_disparity.int32");
		FILE *tmpfile_d;
		tmpfile_d = fopen(tmp_d_file, "wb");
		fwrite(tmp_d, sizeof(int), SAI->nr*SAI->nc, tmpfile_d);
		fclose(tmpfile_d);

		delete[](tmp_d);
		
		FILE *tmpfile_im_labels;
		tmpfile_im_labels = fopen(SAI->path_label_im, "wb");
		fwrite(SAI->label_im, sizeof(int), SAI->nr*SAI->nc, tmpfile_im_labels);
		fclose(tmpfile_im_labels);

		unloadLabels(SAI);
		unloadInverseDepth(SAI);

	}

	/* warp labeling for all views with parents.*/
	for (int ii = 0; ii < global_params->n_views_total; ii++) {

		view *SAI = LF + ii;

		if (SAI->n_references > 0) {

			printf("Warping segmentation from parents (reference views) at view %03d_%03d\t", SAI->c, SAI->r);

			for (int ij = 0; ij < SAI->n_references; ij++) {

				printf("%d ", ij);

				view *ref_view = LF + SAI->references[ij];

				float y0 = ref_view->y;
				float x0 = ref_view->x;

				float y1 = SAI->y;
				float x1 = SAI->x;

				loadInverseDepth(ref_view);

				float *inverse_depth_view0 = new float[SAI->nr*SAI->nc]();
				for (int ijk = 0; ijk < SAI->nr*SAI->nc; ijk++) {
					*(inverse_depth_view0 + ijk) = static_cast<float>(ref_view->depth[ijk]);
					*(inverse_depth_view0 + ijk) = *(inverse_depth_view0 + ijk) - static_cast<float>(SAI->min_inv_d);
					*(inverse_depth_view0 + ijk) = *(inverse_depth_view0 + ijk) / static_cast<float>(1 << D_DEPTH);
				}

				unloadInverseDepth(ref_view);

				float *DM_COL = new float[SAI->nr*SAI->nc]();
				float *DM_ROW = new float[SAI->nr*SAI->nc]();

				getDisparity(y0, y1, x0, x1, inverse_depth_view0, ref_view->nr, ref_view->nc, DM_COL, DM_ROW);

				loadLabels(ref_view);

				int32_t *warpedLabels = new int32_t[SAI->nr*SAI->nc]();
				float *warped_inverse_depth = new float[SAI->nr*SAI->nc]();

				warp_0_to_1(DM_ROW, DM_COL, ref_view->label_im, SAI->nr, SAI->nc, 1, warpedLabels, warped_inverse_depth, inverse_depth_view0);

				unloadLabels(ref_view);

				writeWarpedLabelIm(SAI, ref_view, warpedLabels);

				delete[](warped_inverse_depth);
				delete[](warpedLabels);
				delete[](inverse_depth_view0);

				delete[](DM_COL);
				delete[](DM_ROW);

				

			}

			printf("\n");

		}

	}

	/* predict in default mode (VM1.0,VM1.1) and obtain residual for COLOR at all views */
	for (int ii = 0; ii < global_params->n_views_total; ii++) {

		view *SAI = LF + ii;

		printf("Predicting color view %03d_%03d\t", SAI->c, SAI->r);

		loadOriginalColor(SAI);

		SAI->color = new unsigned short[SAI->nr*SAI->nc * 3]();

		/* color prediction */
		if (SAI->n_references > 0) {

			initializeWarpingArrays(SAI);

			warp_all_references_to_current_view(SAI);

			initViewW(SAI);

			/* get LS weights */
			if (SAI->sigma < 0.0001) {
				getViewMergingLSWeights_N(SAI);
				/* merge color with prediction */
				mergeWarped_N(SAI, 3);
				/* hole filling for color*/
				holefilling(SAI->color, 3, SAI->nr, SAI->nc, 0);
			}
			else {

				int startt = clock();

				mergeMedian_N(SAI, 3);
				holefilling(SAI->color, 3, SAI->nr, SAI->nc, 0);

				//std::cout << "time elapsed in color median merging\t" << (float)((int)clock() - startt) / CLOCKS_PER_SEC << "\n";

				double psnr_med = getYCbCr_422_PSNR(SAI->color, SAI->original_color_view, SAI->nr, SAI->nc, 3, 10);

				unsigned short *tmp_m = new unsigned short[SAI->nr*SAI->nc*3]();
				memcpy(tmp_m, SAI->color, sizeof(unsigned short)*SAI->nr*SAI->nc * 3);

				

				if (global_params->STD_SEARCH) {

					float stdi = 0;
					double psnr_w = 0;

					for (SAI->sigma = 0; SAI->sigma < STD_SEARCH_HIGH; SAI->sigma += STD_SEARCH_STEP) {

						getGeomWeight(SAI);
						mergeWarped_N(SAI, 3);
						holefilling(SAI->color, 3, SAI->nr, SAI->nc, 0);

						double tpsnr = getYCbCr_422_PSNR(SAI->color, SAI->original_color_view, SAI->nr, SAI->nc, 3, BIT_DEPTH);

						if (tpsnr > psnr_w) {
							psnr_w = tpsnr;
							stdi = SAI->sigma;
						}

					}

					SAI->sigma = stdi;

					//printf("PSNR RGB median:%f\tPSNR RGB weights:%f\t SAI->std: %f\n", psnr_med, psnr_w, SAI->sigma);

				}

				getGeomWeight(SAI);
				mergeWarped_N(SAI, 3);
				holefilling(SAI->color, 3, SAI->nr, SAI->nc, 0);

				double psnr_w = getYCbCr_422_PSNR(SAI->color, SAI->original_color_view, SAI->nr, SAI->nc, 3, BIT_DEPTH);

				if (psnr_w < psnr_med) {
					delete[](SAI->color);
					SAI->color = tmp_m;
					SAI->use_median = true;
					SAI->sigma = 0.0;
				}
				else {
					delete[](tmp_m);
				}
			}

			if (SAI->seg_vp != nullptr) {
				delete[](SAI->seg_vp);
				SAI->seg_vp = nullptr;
			}

			deinitializeWarpingArrays(SAI);

			SAI->merge_psnr = getYCbCr_422_PSNR(SAI->color, SAI->original_color_view, SAI->nr, SAI->nc, 3, BIT_DEPTH);
		}


		if (SAI->NNt > 0 && SAI->Ms > 0)
		{

			SAI->use_global_sparse = true;

			int startt = clock();

			getGlobalSparseFilter(SAI);
			applyGlobalSparseFilter(SAI);

			//std::cout << "time elapsed in getGlobalSparseFilter()\t" << (float)( (int)clock() - startt ) / CLOCKS_PER_SEC << "\n";


			SAI->sparse_psnr = getYCbCr_422_PSNR(SAI->color, SAI->original_color_view, SAI->nr, SAI->nc, 3, BIT_DEPTH);
			
		}

		get_and_write_color_residual_JP2(SAI, kdu_compress_path, kdu_expand_path );

		aux_write16PGMPPM(SAI->path_out_ppm, SAI->nc, SAI->nr, 3, SAI->color);

		SAI->final_psnr = getYCbCr_422_PSNR(SAI->color, SAI->original_color_view, SAI->nr, SAI->nc, 3, BIT_DEPTH);
		
		/* to reduce memory usage */
		unloadColor(SAI);
		unloadOriginalColor(SAI);

		printf("%2.3f\t%2.3f\t%2.3f\n", SAI->merge_psnr, SAI->sparse_psnr, SAI->final_psnr);
	}

	/* write bitstream */
	char path_out_LF_data[1024];
	sprintf(path_out_LF_data, "%s%c%s", output_dir, '/', "output.LF");

	/* our bitstream starts here */
	FILE *output_LF_file;
	output_LF_file = fopen(path_out_LF_data, "wb");

	int n_bytes_prediction = 0;
	int n_bytes_residual = 0;

	int yuv_transform_s = global_params->YUV_TRANSFORM ? 1 : 0;

	n_bytes_prediction += static_cast<int>(fwrite(&global_params->n_views_total, sizeof(int), 1, output_LF_file) ) * sizeof(int);
	n_bytes_prediction += static_cast<int>(fwrite(&LF->nr, sizeof(int), 1, output_LF_file)) * sizeof(int); // needed only once per LF
	n_bytes_prediction += static_cast<int>(fwrite(&LF->nc, sizeof(int), 1, output_LF_file)) * sizeof(int); // 
	n_bytes_prediction += static_cast<int>(fwrite(&yuv_transform_s, sizeof(int), 1, output_LF_file)) * sizeof(int);
	n_bytes_prediction += static_cast<int>(fwrite(&LF->min_inv_d, sizeof(unsigned short), 1, output_LF_file)) * sizeof(unsigned short);

	for (int ii = 0; ii < global_params->n_views_total; ii++) {

		view *SAI = LF + ii;

		int n_bytes_prediction_per_view = 0;
		int n_bytes_residual_per_view = 0;

		printf("PSNR YUV (%03d_%03d): %2.3f\n", SAI->r, SAI->c, SAI->final_psnr);

		viewHeaderToCodestream(n_bytes_prediction_per_view, SAI, output_LF_file);

		if (SAI->residual_rate_color > 0) {

			if (SAI->yuv_transform && global_params->YUV_TRANSFORM) {

				int ncomp_r = SAI->has_chrominance ? 3 : 1;

				for (int icomp = 0; icomp < ncomp_r; icomp++) {

					writeResidualToDisk(SAI->ycbcr_jp2_names[icomp], output_LF_file, n_bytes_residual_per_view, JP2_dict);

				}
			}
			else {

				writeResidualToDisk(SAI->jp2_residual_path_jp2, output_LF_file, n_bytes_residual_per_view, JP2_dict);
			}
		}

		if (SAI->residual_rate_depth > 0 && SAI->depth_file_exist) {

			writeResidualToDisk(SAI->jp2_residual_depth_path_jp2, output_LF_file, n_bytes_residual_per_view, JP2_dict);

		}

		n_bytes_prediction += n_bytes_prediction_per_view;
		n_bytes_residual += n_bytes_residual_per_view;

	}

	fclose(output_LF_file);

	double num_pixels = global_params->n_views_total * LF->nr * LF->nc;

	long enc_file_size = aux_GetFileSize(path_out_LF_data);

	printf("Output: %s\nsize: %d kB\t bpp: %2.4f\n", path_out_LF_data, enc_file_size/1000, static_cast<double>( enc_file_size )/num_pixels);

	for (int ii = 0; ii < global_params->n_views_total; ii++)
	{

		cleanView(LF + ii);

	}

	delete[](LF);

	exit(0);
}