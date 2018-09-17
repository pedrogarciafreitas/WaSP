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
#include "connected_components.hh"
#include "fastols.hh"
#include "motioncompensation.hh"
#include "kmeans.hh"


#define USE_difftest_ng false

#define YUV_SEARCH_LOW 6.80f
#define YUV_SEARCH_HIGH 7.80f
#define YUV_SEARCH_STEP 0.20f
#define YUV_RATIO_DEFAULT 7.20f

#define STD_SEARCH_LOW 10
#define STD_SEARCH_HIGH 250
#define STD_SEARCH_STEP 10

#define SAVE_PARTIAL_WARPED_VIEWS false

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

	//const char *difftest_call = "C:/Local/astolap/Data/JPEG_PLENO/RIO_INPUT/ScriptJan2018/ScriptSolution/difftest_ng.exe --toycbcr --psnr ";
	//const char *difftest_call_pgm = "C:/Local/astolap/Data/JPEG_PLENO/RIO_INPUT/ScriptJan2018/ScriptSolution/difftest_ng.exe --psnr ";

	//const char *difftest_call = "D:/JPEG_VM_01/DevelopmentC/Pekka/difftest_ng.exe --toycbcr --psnr ";
	//const char *difftest_call_pgm = "D:/JPEG_VM_01/DevelopmentC/Pekka/difftest_ng.exe --psnr ";

	FILE *filept;

	filept = fopen(config_file, "rb");

	int n_views_total;
	fread(&n_views_total, sizeof(int), 1, filept); /*reading*/

	view *LF = new view[n_views_total](); /* one dimensional view vector */

	int maxR = 0, maxC = 0;

	bool YUV_TRANSFORM = false;
	bool YUV_RATIO_SEARCH = false;
	bool STD_SEARCH = false;

	const bool RESIDUAL_16BIT_bool = RESIDUAL_16BIT ? 1 : 0;

	int yuv_transform_s,yuv_ratio_search_s,std_search_s;
	fread(&yuv_transform_s, sizeof(int), 1, filept); /*reading*/
	fread(&yuv_ratio_search_s, sizeof(int), 1, filept); /*reading*/
	fread(&std_search_s, sizeof(int), 1, filept); /*reading*/

	YUV_TRANSFORM = yuv_transform_s > 0 ? true : false;
	YUV_RATIO_SEARCH = yuv_ratio_search_s > 0 ? true : false;
	STD_SEARCH = std_search_s > 0 ? true : false;

	unsigned short MINIMUM_DEPTH = 0;

	for (int ii = 0; ii < n_views_total; ii++) {

		view *SAI = LF + ii;

		initView(SAI);

		SAI->i_order = ii;

		fread(&(SAI->r), sizeof(int), 1, filept); /*reading*/
		fread(&(SAI->c), sizeof(int), 1, filept); /*reading*/

		/* find number of rows and cols */
		maxR = (SAI->r+1 > maxR) ? SAI->r + 1 : maxR;
		maxC = (SAI->c+1 > maxC) ? SAI->c + 1 : maxC;

		int xx = 0, yy = 0;

		fread(&xx, sizeof(int), 1, filept); /*reading*/
		fread(&yy, sizeof(int), 1, filept); /*reading*/

		if ( abs(xx) > 0) {
			SAI->has_x_displacement = true;
			SAI->x = float(xx) / 100000;
		}

		if ( abs(yy) > 0) {
			SAI->has_y_displacement = true;
			SAI->y = float(yy) / 100000;
		}

		int rate_color, rate_depth;

		fread(&rate_color, sizeof(int), 1, filept); /*reading*/
		fread(&rate_depth, sizeof(int), 1, filept); /*reading*/

		SAI->residual_rate_color = ((float)rate_color) / 100000;
		SAI->residual_rate_depth = ((float)rate_depth) / 100000;

		fread(&SAI->Ms, sizeof(int), 1, filept); /*reading*/
		fread(&SAI->NNt, sizeof(int), 1, filept); /*reading*/

		int stdd = 0;

		fread(&stdd, sizeof(int), 1, filept); /*reading*/
		SAI->stdd = ((float)stdd) / 100000;

		unsigned short tmpminv;

		fread(&tmpminv, sizeof(int), 1, filept); /*reading, if we have negative inverse depth,
												   for example in lenslet, we need to subtract min_inv_d
												   from the inverse depth maps*/
		if (ii == 0) {
			if (tmpminv > 0) {
				MINIMUM_DEPTH = tmpminv;
			}
		}

		if (MINIMUM_DEPTH > 0) {
			//SAI->has_min_inv_depth = true;
			SAI->min_inv_d = (int)MINIMUM_DEPTH;
		}

		fread(&(SAI->n_references), sizeof(int), 1, filept); /*reading*/

		if (SAI->n_references > 0) {

			SAI->has_color_references = true;

			SAI->references = new int[SAI->n_references]();

			fread(SAI->references, sizeof(int), SAI->n_references, filept); /*reading*/

		}

		fread(&(SAI->n_depth_references), sizeof(int), 1, filept); /*reading*/

		if (SAI->n_depth_references > 0) {

			SAI->has_depth_references = true;

			SAI->depth_references = new int[SAI->n_depth_references]();

			fread(SAI->depth_references, sizeof(int), SAI->n_depth_references, filept); /*reading*/

		}

		fread(&SAI->has_segmentation, sizeof(int), 1, filept); // new,13.06.18 /*reading*/

		sprintf(SAI->path_input_ppm, "%s%c%03d_%03d%s", input_dir, '/', SAI->c, SAI->r, ".ppm");
		sprintf(SAI->path_input_pgm, "%s%c%03d_%03d%s", input_dir, '/', SAI->c, SAI->r, ".pgm");

		sprintf(SAI->path_input_seg, "%s%c%03d_%03d%s", input_dir, '/', SAI->c, SAI->r, "_segmentation.pgm");

		sprintf(SAI->path_out_ppm, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, ".ppm");
		sprintf(SAI->path_out_pgm, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, ".pgm");

	}
	fclose(filept);

	/* 2D array format for views, useful in some cases */
	view ***LF_mat = new view**[maxR]();
	for (int ii = 0; ii < maxR; ii++) {
		LF_mat[ii] = new view*[maxC]();
	}

	for (int r = 0; r < maxR; r++) {
		for (int c = 0; c < maxC; c++) {
			LF_mat[r][c] = NULL;
		}
	}

	for (int ii = 0; ii < n_views_total; ii++) {
		view *SAI = LF + ii;
		LF_mat[SAI->r][SAI->c] = SAI;
	}

	char path_out_LF_data[1024];
	sprintf(path_out_LF_data, "%s%c%s", output_dir, '/', "output.LF");

	/* debugging codestream overhead reduction */
	char path_codestream[1024];
	sprintf(path_codestream, "%s%c%s", output_dir, '/', "LF.codestream");
	FILE *tmp_codestream;
	tmp_codestream = fopen(path_codestream, "wb");
	fclose(tmp_codestream);

	/* our bitstream starts here */
	FILE *output_LF_file;
	output_LF_file = fopen(path_out_LF_data, "wb");
	fwrite(&n_views_total, sizeof(int), 1, output_LF_file);
	fclose(output_LF_file);

	bool global_header_written = false;

	FILE *output_results_file;
	char output_results_filename[1024];
	sprintf(output_results_filename, "%s/%s", output_dir, "results.txt");

	char output_results[2048];
	int output_buffer_length = 0;

	//output_buffer_length += sprintf(output_results + output_buffer_length, "%s",
	//	"ROW\tCOL\tPSNR1\tPSNR2\tPSNR3\tPSNR4\tSTD\tYUVRATIO\tPSNR5\t\tbytes_prediction\t\tbytes_residual");
	output_results_file = fopen(output_results_filename, "w");
	//fprintf(output_results_file, "%s\n", output_results);
	fclose(output_results_file);


	/* to get effiency from multiple JP2 files, we remove parts of the files
	which are repetative over all files. For this we have a minimalistic 
	dictionary method. */

	std::vector<std::vector<unsigned char>> JP2_dict;

	for (int ii = 0; ii < n_views_total; ii++) {

		view *SAI = LF + ii;

		//char output_results[1024];
		memset(output_results, 0x00, sizeof(char) * sizeof(output_results)/sizeof(char));
		output_buffer_length = 0;
		output_buffer_length += sprintf(output_results + output_buffer_length, "%03d\t%03d", SAI->r, SAI->c);

		unsigned short *original_color_view = NULL;
		unsigned short *original_depth_view = NULL;

		int nc1, nr1, ncomp1;
		aux_read16PGMPPM(SAI->path_input_ppm, SAI->nc, SAI->nr, ncomp1, original_color_view);

		///* debug yuv */
		//unsigned short *ycbcr = new unsigned short[SAI->nr*SAI->nc * 3]();
		//RGB2YCbCr(original_color_view, ycbcr, SAI->nr, SAI->nc, 10);
		//aux_write16PGMPPM("C:/Local/astolap/Data/JPEG_PLENO/TUT-HDCA_tmp_output_lenslet/YCbCr.ppm",
		//	SAI->nc, SAI->nr, 3, ycbcr);
		//unsigned short *rgb = new unsigned short[SAI->nr*SAI->nc * 3]();
		//YCbCr2RGB(ycbcr, rgb, SAI->nr, SAI->nc, 10);
		//aux_write16PGMPPM("C:/Local/astolap/Data/JPEG_PLENO/TUT-HDCA_tmp_output_lenslet/RGB.ppm",
		//	SAI->nc, SAI->nr, 3, rgb);
		////aux_write16PGMPPM("C:/Local/astolap/Data/JPEG_PLENO/TUT-HDCA_tmp_output_lenslet/original.ppm",
		////	SAI->nc, SAI->nr, 3, original_color_view);

		//exit(0);

		bool depth_file_exist = false;
		
		if (SAI->residual_rate_depth > 0) {
			depth_file_exist = aux_read16PGMPPM(SAI->path_input_pgm, nc1, nr1, ncomp1, original_depth_view);
		}

		SAI->color = new unsigned short[SAI->nr*SAI->nc * 3]();
		SAI->depth = new unsigned short[SAI->nr*SAI->nc]();

		/*----------------DISPARITY------------------------------------------------------------*/
		/* forward warp depth */
		if (SAI->has_depth_references) {
			predictDepth(SAI, LF);
		}

		char pgm_residual_depth_path[1024];
		char jp2_residual_depth_path_jp2[1024];

		if (SAI->residual_rate_depth > 0 && depth_file_exist) { /* residual depth if needed */

			sprintf(pgm_residual_depth_path, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_depth_residual.pgm");

			sprintf(jp2_residual_depth_path_jp2, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_depth_residual.jp2");

			encodeResidualJP2(SAI->nr, SAI->nc, original_depth_view, SAI->depth, pgm_residual_depth_path,
				kdu_compress_path, jp2_residual_depth_path_jp2, SAI->residual_rate_depth, 1, 0, 1);

			decodeResidualJP2(SAI->depth, kdu_expand_path, jp2_residual_depth_path_jp2, pgm_residual_depth_path, ncomp1, 0, (1 << 16) - 1, 1);

			SAI->has_depth_residual = true;
		}

		//int QD = 150;

		int *tmp_d = new int[SAI->nr*SAI->nc]();
		for (int iii = 0; iii < SAI->nr*SAI->nc; iii++) {
			*(tmp_d + iii) = (int)(*(SAI->depth + iii));// / QD;
		}

		getKmeansQuantized(K_MEANS_CLUSTERS, tmp_d, SAI->nr*SAI->nc, K_MEANS_ITERATIONS); /*inplace assignment to tmp_d*/
																				  /*temporary test of label_im with blocks*/
		/*int *blok_im = tmp_d;
		int blocksize = 32;
		int bi = 1;
		for (int rr = 0; rr < SAI->nr; rr = rr + blocksize) {
			for (int cc = 0; cc < SAI->nc; cc = cc + blocksize) {

				for (int rr1 = 0; rr1 < blocksize; rr1++) {
					for (int cc1 = 0; cc1 < blocksize; cc1++) {

						int icc = cc + cc1;
						int irr = rr + rr1;

						if (icc >= 0 && icc < SAI->nc && irr >= 0 && irr < SAI->nr) {
							int ind = irr + (icc)*SAI->nr;
							blok_im[ind] = bi;
						}
					}
				}
				bi++;
			}
		}*/

		int nregions = 0;
		int *reg_histogram = 0;
		int *label_im = get_labels(tmp_d, SAI->nr, SAI->nc, nregions, reg_histogram);

		SAI->label_im = label_im;
		SAI->nregions = nregions;
		SAI->reg_histogram = reg_histogram;

		unsigned short *labels = new unsigned short[SAI->nr*SAI->nc]();
		for (int ii = 0; ii < SAI->nr*SAI->nc; ii++) {
			*(labels + ii) = (unsigned short) *(label_im + ii);
		}

		char labels_file[1024];
		sprintf(labels_file, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_labels.pgm");
		aux_write16PGMPPM(labels_file, SAI->nc, SAI->nr, 1, labels);

		char tmp_d_file[1024];
		sprintf(tmp_d_file, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_kmeans_disparity.int32");
		FILE *tmpfile_d;
		tmpfile_d = fopen(tmp_d_file, "wb");
		fwrite(tmp_d, sizeof(int), SAI->nr*SAI->nc, tmpfile_d);
		fclose(tmpfile_d);

		sprintf(SAI->path_label_im, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_im_labels.int32");
		FILE *tmpfile_im_labels;
		tmpfile_im_labels = fopen(SAI->path_label_im, "wb");
		fwrite(SAI->label_im, sizeof(int), SAI->nr*SAI->nc, tmpfile_im_labels);
		fclose(tmpfile_im_labels);

		if (MOTION_VECTORS) {

			sortRegionsBySize(SAI);

			//sortRegionsByDisparity(SAI);

			SAI->use_motion_vectors = true;

			/* loop through all views that use this current view (SAI) as a reference. compute motion vectors for regions */

			for (int vi = 0; vi < n_views_total; vi++) {

				for (int ef = 0; ef < (LF + vi)->n_depth_references; ef++) {
					if ((LF + vi)->depth_references[ef] == ii) {
						getMotionVectorsView0_to_View1(SAI, LF+vi);
					}
				}


				for (int ef = 0; ef < (LF + vi)->n_references; ef++) {
					if ((LF + vi)->references[ef] == ii) {
						getMotionVectorsView0_to_View1(SAI, LF+vi);
					}
				}

			}
			
		}

		delete[](labels);
		delete[](tmp_d);

		/*----------------DISPARITY ENDS------------------------------------------------------------*/

		/* color prediction */
		/* forward warp color */

		if (SAI->n_references > 0) {

			/* holds partial warped views for ii */
			unsigned short **warped_color_views = new unsigned short*[SAI->n_references]();
			unsigned short **warped_depth_views = new unsigned short*[SAI->n_references]();
			float **DispTargs = new float*[SAI->n_references]();

			for (int ij = 0; ij < SAI->n_references; ij++)
			{

				view *ref_view = LF + SAI->references[ij];

				int tmp_w, tmp_r, tmp_ncomp;

				aux_read16PGMPPM(ref_view->path_out_pgm, tmp_w, tmp_r, tmp_ncomp, ref_view->depth);
				aux_read16PGMPPM(ref_view->path_out_ppm, tmp_w, tmp_r, tmp_ncomp, ref_view->color);

				/* FORWARD warp color AND depth */
				warpView0_to_View1(ref_view, SAI, warped_color_views[ij], warped_depth_views[ij], DispTargs[ij]);

				delete[](ref_view->depth);
				delete[](ref_view->color);

				ref_view->depth = NULL;
				ref_view->color = NULL;

				char tmp_str[1024];

				if (SAVE_PARTIAL_WARPED_VIEWS) {

					sprintf(tmp_str, "%s/%03d_%03d%s%03d_%03d%s", output_dir, (ref_view)->c, (ref_view)->r, "_warped_to_", SAI->c, SAI->r, ".ppm");
					aux_write16PGMPPM(tmp_str, SAI->nc, SAI->nr, 3, warped_color_views[ij]);

					sprintf(tmp_str, "%s/%03d_%03d%s%03d_%03d%s", output_dir, (ref_view)->c, (ref_view)->r, "_warped_to_", SAI->c, SAI->r, ".pgm");
					aux_write16PGMPPM(tmp_str, SAI->nc, SAI->nr, 1, warped_depth_views[ij]);

				}

				//FILE *tmpf;
				//sprintf(tmp_str, "%s%03d_%03d%s%03d_%03d%s", output_dir, (ref_view)->c, (ref_view)->r, "_warped_to_", SAI->c, SAI->r, "_DispTarg.float");
				//tmpf = fopen(tmp_str, "wb");
				//fwrite(DispTargs[ij], sizeof(float), SAI->nr * SAI->nc, tmpf);
				//fclose(tmpf);

			}

			initViewW(SAI, DispTargs);

			/* get LS weights */
			if (SAI->stdd < 0.0001) {
				getViewMergingLSWeights_N(SAI, warped_color_views, DispTargs, original_color_view);
				/* merge color with prediction */
				mergeWarped_N(warped_color_views, DispTargs, SAI, 3);
				/* hole filling for color*/
				holefilling(SAI->color, 3, SAI->nr, SAI->nc, 0);
				if (1) {
					unsigned short *tmp_color = new unsigned short[SAI->nr*SAI->nc * 3]();
					memcpy(tmp_color, SAI->color, sizeof(unsigned short)*SAI->nr*SAI->nc * 3);

					delete[](SAI->color);
					SAI->color = new unsigned short[SAI->nr*SAI->nc * 3]();

					view *SAI_tmp = new view[1]();

					initView(SAI_tmp);

					SAI_tmp->n_references = SAI->n_references;
					SAI_tmp->references = new int[SAI->n_references]();
					memcpy(SAI_tmp->references, SAI->references, sizeof(int)*SAI->n_references);

					SAI_tmp->nr = SAI->nr;
					SAI_tmp->nc = SAI->nc;

					SAI_tmp->color = new unsigned short[SAI->nr*SAI->nc * 3]();


					/* holds partial warped views for ii */
					float **DispTargs_tmp = new float*[SAI->n_references]();
					for (int ij = 0; ij < SAI_tmp->n_references; ij++)
					{
						DispTargs_tmp[ij] = new float[SAI->nr*SAI->nc]();
					}

					for (int iR = 0; iR < SAI->nregions; iR++) {

						if (SAI->reg_histogram[iR] >= 256) {

							for (int ij = 0; ij < SAI_tmp->n_references; ij++)
							{

								memcpy(DispTargs_tmp[ij], DispTargs[ij], sizeof(float)*SAI->nr*SAI->nc);

								float *tmps = DispTargs_tmp[ij];
								for (int iuj = 0; iuj < SAI->nr*SAI->nr; iuj++) {
									if (SAI->label_im[iuj] != iR) {
										tmps[iuj] = INIT_DISPARITY_VALUE;
									}
								}

							}

							initViewW(SAI_tmp, DispTargs_tmp);

							getViewMergingLSWeights_N(SAI_tmp, warped_color_views, DispTargs_tmp, original_color_view);
							/* merge color with prediction */
							mergeWarped_N(warped_color_views, DispTargs_tmp, SAI_tmp, 3);
							/* hole filling for color*/
							//holefilling(SAI_tmp->color, 3, SAI->nr, SAI->nc, 0);

							///* clean */
							//for (int ij = 0; ij < SAI_tmp->n_references; ij++)
							//{
							//	delete[](DispTargs_tmp[ij]);
							//}
							//delete[](DispTargs_tmp);

							for (int iuj = 0; iuj < SAI->nr*SAI->nr; iuj++) {
								if (SAI->label_im[iuj] == iR) {
									SAI->color[iuj] = SAI_tmp->color[iuj];
									SAI->color[iuj + SAI->nr*SAI->nc] = SAI_tmp->color[iuj + SAI->nr*SAI->nc];
									SAI->color[iuj + SAI->nr*SAI->nc * 2] = SAI_tmp->color[iuj + SAI->nr*SAI->nc * 2];
								}
							}


							delete[](SAI_tmp->seg_vp);
							SAI_tmp->seg_vp = NULL;
							delete[](SAI_tmp->merge_weights);
							SAI_tmp->merge_weights = NULL;
							delete[](SAI_tmp->merge_weights_float);
							SAI_tmp->merge_weights_float = NULL;
							delete[](SAI_tmp->bmask);
							SAI_tmp->bmask = NULL;
							delete[](SAI_tmp->number_of_pixels_per_region);
							SAI_tmp->number_of_pixels_per_region = NULL;


							printf("\tiR %i", iR);
						}

					}

					delete[](SAI_tmp->references);
					delete[](SAI_tmp->color);
					delete[](SAI_tmp);

					/* clean */
					for (int ij = 0; ij < SAI_tmp->n_references; ij++)
					{
						delete[](DispTargs_tmp[ij]);
					}
					delete[](DispTargs_tmp);

					///* hole filling for color*/
					//holefilling(SAI->color, 3, SAI->nr, SAI->nc, 0);
					for (int ijj = 0; ijj < SAI->nr*SAI->nc; ijj++) {
						if (SAI->color[ijj] > 0 || SAI->color[ijj + SAI->nr*SAI->nc] > 0 || SAI->color[ijj + 2 * SAI->nr*SAI->nc]) {
							tmp_color[ijj] = SAI->color[ijj];
							tmp_color[ijj + SAI->nr*SAI->nc] = SAI->color[ijj + SAI->nr*SAI->nc];
							tmp_color[ijj + SAI->nr*SAI->nc * 2] = SAI->color[ijj + SAI->nr*SAI->nc * 2];
						}
					}
					memcpy(SAI->color, tmp_color, sizeof(unsigned short)*SAI->nr*SAI->nc * 3);
					delete[](tmp_color);
				}
			}
			else {
				/* get baseline with median, we then study whether weighting improves */
				/* merge color with median */
				int startt = clock();
				mergeMedian_N(warped_color_views, DispTargs, SAI, 3);
				std::cout << "time elapsed in color median merging\t" << (float)( (int)clock() - startt ) / CLOCKS_PER_SEC << "\n";
				holefilling(SAI->color, 3, SAI->nr, SAI->nc, 0);

				//double psnr_med = getYCbCr_422_PSNR(SAI->color, original_color_view, SAI->nr, SAI->nc, 3, 10);
				double psnr_med = getYCbCr_422_PSNR(SAI->color, original_color_view, SAI->nr, SAI->nc, 3, 10);

				unsigned short *tmp_m = new unsigned short[SAI->nr*SAI->nc*3]();
				memcpy(tmp_m, SAI->color, sizeof(unsigned short)*SAI->nr*SAI->nc * 3);

				double psnr_w = 0;

				if (STD_SEARCH) {

					float stdi = 0;

					for (float stds = STD_SEARCH_LOW; stds < STD_SEARCH_HIGH; stds += STD_SEARCH_STEP) {
						SAI->stdd = stds;
						/* we don't use LS weights but something derived on geometric distance in view array*/
						getGeomWeight(SAI, LF);
						/* merge color with prediction */
						mergeWarped_N(warped_color_views, DispTargs, SAI, 3);
						/* hole filling for color*/
						holefilling(SAI->color, 3, SAI->nr, SAI->nc, 0);
						double tpsnr = getYCbCr_422_PSNR(SAI->color, original_color_view, SAI->nr, SAI->nc, 3, BIT_DEPTH);
						//double tpsnr = getYCbCr_422_PSNR(SAI->color, original_color_view, SAI->nr, SAI->nc, 3, 10);

						if (tpsnr > psnr_w) {
							psnr_w = tpsnr;
							stdi = SAI->stdd;
						}

					}

					SAI->stdd = stdi;

					printf("PSNR RGB median:%f\tPSNR RGB weights:%f\t SAI->std: %f\n", psnr_med, psnr_w, SAI->stdd);

				}

				/* we don't use LS weights but something derived on geometric distance in view array*/
				getGeomWeight(SAI, LF);
				/* merge color with prediction */
				mergeWarped_N(warped_color_views, DispTargs, SAI, 3);
				/* hole filling for color*/
				holefilling(SAI->color, 3, SAI->nr, SAI->nc, 0);
				psnr_w = getYCbCr_422_PSNR(SAI->color, original_color_view, SAI->nr, SAI->nc, 3, BIT_DEPTH);

				if (psnr_w < psnr_med) {
					delete[](SAI->color);
					SAI->color = tmp_m;
					SAI->use_median = true;
					SAI->stdd = 0.0;
				}
				else {
					delete[](tmp_m);
				}
			}

			/* clean */
			for (int ij = 0; ij < SAI->n_references; ij++)
			{
				delete[](warped_color_views[ij]);
				delete[](warped_depth_views[ij]);
				delete[](DispTargs[ij]);
			}

			delete[](warped_color_views);
			delete[](warped_depth_views);
			delete[](DispTargs);

			unsigned short *tmp_med_im = new unsigned short[SAI->nr*SAI->nc * 3]();

			for (int icomp = 0; icomp < 3; icomp++) {
				medfilt2D((SAI->color) + icomp*SAI->nr*SAI->nc, tmp_med_im + icomp*SAI->nr*SAI->nc, 3, SAI->nr, SAI->nc);
			}

			double psnr_no_median_filter = getYCbCr_422_PSNR(SAI->color, original_color_view, SAI->nr, SAI->nc, 3, BIT_DEPTH);
			double psnr_median_filtered = getYCbCr_422_PSNR(tmp_med_im, original_color_view, SAI->nr, SAI->nc, 3, BIT_DEPTH);

			if (psnr_median_filtered > psnr_no_median_filter) {
				SAI->use_median_filter = true;
				memcpy(SAI->color, tmp_med_im, sizeof(unsigned short)*SAI->nr*SAI->nc * 3);
			}

			delete[](tmp_med_im);

		}


		double psnr_without_sparse = 0;

		if (SAI->n_references > 0) {

			//aux_write16PGMPPM(SAI->path_out_ppm, SAI->nc, SAI->nr, 3, SAI->color);

			//psnr_result = USE_difftest_ng ? getPSNR(NULL, SAI->path_out_ppm, SAI->path_input_ppm, difftest_call) : 0;
		
			psnr_without_sparse = getYCbCr_422_PSNR(SAI->color, original_color_view, SAI->nr, SAI->nc, 3, BIT_DEPTH);

			output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", psnr_without_sparse);

		}
		else {
			output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", 0.0);
		}

		unsigned short *colorview_temp = new unsigned short[SAI->nr*SAI->nc * 3]();
		memcpy(colorview_temp, SAI->color, sizeof(unsigned short)*SAI->nr*SAI->nc * 3);

		double psnr_with_sparse = 0;

		

		if (SAI->NNt > 0 && SAI->Ms > 0)
		{

			SAI->use_global_sparse = true;

			int startt = clock();

			getGlobalSparseFilter(SAI, original_color_view);

			std::cout << "time elapsed in getGlobalSparseFilter()\t" << (float)( (int)clock() - startt ) / CLOCKS_PER_SEC << "\n";

			applyGlobalSparseFilter(SAI);

			psnr_with_sparse = getYCbCr_422_PSNR(SAI->color, original_color_view, SAI->nr, SAI->nc, 3, BIT_DEPTH);

			//aux_write16PGMPPM(SAI->path_out_ppm, SAI->nc, SAI->nr, 3, SAI->color);
			//psnr_result = USE_difftest_ng ? getPSNR(NULL, SAI->path_out_ppm, SAI->path_input_ppm, difftest_call) : 0;
			
		}

		/* write both with and without sparse to disk */
		char w_sparse[1024], wo_sparse[1024];
		sprintf(w_sparse, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_wsparse.ppm");
		sprintf(wo_sparse, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_wosparse.ppm");

		aux_write16PGMPPM(w_sparse, SAI->nc, SAI->nr, 3, SAI->color);
		aux_write16PGMPPM(wo_sparse, SAI->nc, SAI->nr, 3, colorview_temp);

		if (SAI->use_global_sparse) { /* check validity of sparse filter */
			if ( psnr_with_sparse<psnr_without_sparse ) //<0.1
			{
				SAI->use_global_sparse = false;
				memcpy(SAI->color, colorview_temp, sizeof(unsigned short)*SAI->nr*SAI->nc * 3);
			}
		}

		if (SAI->use_global_sparse) {
			output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", psnr_with_sparse);
		}
		else {
			output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", 0.0);
		}

		///* get Y prediction error at this stage */
		//unsigned short *ycbcr_original = new unsigned short[SAI->nr*SAI->nc * 3]();
		//unsigned short *ycbcr_prediction = new unsigned short[SAI->nr*SAI->nc * 3]();

		//RGB2YCbCr(original_color_view, ycbcr_original, SAI->nr, SAI->nc, 10);
		//RGB2YCbCr(SAI->color, ycbcr_prediction, SAI->nr, SAI->nc, 10);

		//float *Y_prediction_error = new float[SAI->nr*SAI->nc]();

		//for (int jj = 0; jj < SAI->nr*SAI->nc; jj++) {
		//	Y_prediction_error[jj] = abs( static_cast<float>(ycbcr_original[jj]) -
		//		static_cast<float>(ycbcr_prediction[jj]) );
		//}

		//delete[](ycbcr_original);
		//delete[](ycbcr_prediction);

		//getKmeansQuantized(64, Y_prediction_error, SAI->nr*SAI->nc, K_MEANS_ITERATIONS); /*inplace assignment to tmp_d*/

		//float minval_err = FLT_MAX;
		//for (int jj = 0; jj < SAI->nr*SAI->nc; jj++) {
		//	minval_err = Y_prediction_error[jj] < minval_err ? Y_prediction_error[jj] : minval_err;
		//}

		//int *Y_prediction_error_int32 = new int[SAI->nr*SAI->nc]();

		//for (int jj = 0; jj < SAI->nr*SAI->nc; jj++) {
		//	Y_prediction_error_int32[jj] = static_cast<int>( floor(Y_prediction_error[jj] + minval_err + 0.5) );
		//}

		//int nregions_err = 0;
		//int *reg_histogram_err = 0;
		//int *label_im_err = get_labels(Y_prediction_error_int32, SAI->nr, SAI->nc, nregions_err, reg_histogram_err);

		//char tmplabelchar[1024];
		//sprintf(tmplabelchar, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_im_pred_err_labels.int32");
		//tmpfile_im_labels = fopen(tmplabelchar, "wb");
		//fwrite(label_im_err, sizeof(int), SAI->nr*SAI->nc, tmpfile_im_labels);
		//fclose(tmpfile_im_labels);

		//delete[](Y_prediction_error);
		//delete[](Y_prediction_error_int32);

		double psnr_with_region_sparse = 0.0;
		/* here region sparse */
		if (REGION_SPARSE_ON && ii>4) {

			SAI->Ms = 9;
			SAI->NNt = 3;

			getRegionSparseFilter(SAI, original_color_view);
			std::vector< std::pair< int, int> > valid_regions = validateRegionSparseFilters(SAI, original_color_view);

			int n_valid_regions = static_cast<int>( valid_regions.size() );

			std::vector< int > valid_regions_ir;

			std::vector< std::vector< unsigned char > > region_Regr_final;
			std::vector< std::vector< int32_t > > region_Theta_final;

			for (int ir = 0; ir < n_valid_regions; ir++) {

				valid_regions_ir.push_back(valid_regions.at(ir).second);

				region_Regr_final.push_back(SAI->region_Regr.at(valid_regions.at(ir).first));
				region_Theta_final.push_back(SAI->region_Theta.at(valid_regions.at(ir).first));

			}

			SAI->region_Regr = region_Regr_final;
			SAI->region_Theta = region_Theta_final;
			SAI->valid_regions_ir = valid_regions_ir;

			printf("size of valid regions:\t%i\n", region_Regr_final.size());

			applyRegionSparseFilters(SAI);

			psnr_with_region_sparse = getYCbCr_422_PSNR(SAI->color, original_color_view, SAI->nr, SAI->nc, 3, BIT_DEPTH);

			char w_reg_sparse[1024];
			sprintf(w_reg_sparse, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_wreg_sparse.ppm");
			aux_write16PGMPPM(w_reg_sparse, SAI->nc, SAI->nr, 3, SAI->color);

			SAI->use_region_sparse = true;
		}

		output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", psnr_with_region_sparse);
		output_buffer_length += sprintf(output_results + output_buffer_length, "\t%i", static_cast<int>( SAI->region_Regr.size()) );

		delete[](colorview_temp);

		char ppm_residual_path[1024];
		char jp2_residual_path_jp2[1024];

		char pgm_residual_Y_path[1024];
		char jp2_residual_Y_path_jp2[1024];
		char pgm_residual_Cb_path[1024];
		char jp2_residual_Cb_path_jp2[1024];
		char pgm_residual_Cr_path[1024];
		char jp2_residual_Cr_path_jp2[1024];


		char *ycbcr_pgm_names[3];
		char *ycbcr_jp2_names[3];

		float rate_a1 = 0;

		/* get color residual */
		if (SAI->residual_rate_color > 0)
		{

			/* COLOR residual here, lets try YUV */

			sprintf(pgm_residual_Y_path, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_Y_residual.pgm");
			sprintf(jp2_residual_Y_path_jp2, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_Y_residual.jp2");

			sprintf(pgm_residual_Cb_path, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_Cb_residual.pgm");
			sprintf(jp2_residual_Cb_path_jp2, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_Cb_residual.jp2");

			sprintf(pgm_residual_Cr_path, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_Cr_residual.pgm");
			sprintf(jp2_residual_Cr_path_jp2, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_Cr_residual.jp2");

			sprintf(ppm_residual_path, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_residual.ppm");
			sprintf(jp2_residual_path_jp2, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_residual.jp2");

			ycbcr_pgm_names[0] = pgm_residual_Y_path;
			ycbcr_pgm_names[1] = pgm_residual_Cb_path;
			ycbcr_pgm_names[2] = pgm_residual_Cr_path;

			ycbcr_jp2_names[0] = jp2_residual_Y_path_jp2;
			ycbcr_jp2_names[1] = jp2_residual_Cb_path_jp2;
			ycbcr_jp2_names[2] = jp2_residual_Cr_path_jp2;
			

			if (YUV_TRANSFORM) {

				int offset_v = 0;

				if (RESIDUAL_16BIT) {
					offset_v = (1<<15) - 1;
				}
				else {
					offset_v = (1<<BIT_DEPTH) - 1;
				}

				//float rate_a = 6.5 / 8.0;// 7.2 / 8.0;

				if (YUV_RATIO_SEARCH) {

					float highest_psnr = 0;

					unsigned short *tmp_im = new unsigned short[SAI->nr*SAI->nc * 3]();

					for (float rate_a = YUV_SEARCH_LOW; rate_a <= YUV_SEARCH_HIGH; rate_a += YUV_SEARCH_STEP) {

						memcpy(tmp_im, SAI->color, sizeof(unsigned short)*SAI->nr*SAI->nc * 3);

						encodeResidualJP2_YUV(SAI->nr, SAI->nc, original_color_view, tmp_im, ycbcr_pgm_names,
							kdu_compress_path, ycbcr_jp2_names, SAI->residual_rate_color, 3, offset_v, rate_a / 8.0f, RESIDUAL_16BIT_bool);

						decodeResidualJP2_YUV(tmp_im, kdu_expand_path, ycbcr_jp2_names, ycbcr_pgm_names, 3, offset_v, (1<<BIT_DEPTH) - 1, RESIDUAL_16BIT_bool);

						double psnr_result_yuv = getYCbCr_422_PSNR(tmp_im, original_color_view, SAI->nr, SAI->nc, 3, 10);

						if (psnr_result_yuv > highest_psnr) {
							highest_psnr = (float)psnr_result_yuv;
							rate_a1 = rate_a;
							printf("PSNR YUV:\t%f\tratio\t%f/%f\n", highest_psnr, rate_a1, 8.0);
						}

					}

					delete[](tmp_im);
				}
				else {
					rate_a1 = (float)YUV_RATIO_DEFAULT;
				}

				unsigned short *tmpim_w_yuv = new unsigned short[SAI->nr*SAI->nc * 3]();
				memcpy(tmpim_w_yuv, SAI->color, sizeof(unsigned short)*SAI->nr*SAI->nc * 3);

				encodeResidualJP2_YUV(SAI->nr, SAI->nc, original_color_view, SAI->color, ycbcr_pgm_names,
					kdu_compress_path, ycbcr_jp2_names, SAI->residual_rate_color, 3, offset_v, rate_a1 / (float)8.0, RESIDUAL_16BIT_bool);

				decodeResidualJP2_YUV(tmpim_w_yuv, kdu_expand_path, ycbcr_jp2_names, ycbcr_pgm_names, 3, offset_v, (1<<BIT_DEPTH) - 1, RESIDUAL_16BIT_bool);

				double psnr_result_yuv_w_trans = getYCbCr_422_PSNR(tmpim_w_yuv, original_color_view, SAI->nr, SAI->nc, 3, 10);

				/* also compete against no yuv transformation */

				offset_v = (1<<BIT_DEPTH) - 1;

				unsigned short *tmpim_wo_yuv = new unsigned short[SAI->nr*SAI->nc * 3]();
				memcpy(tmpim_wo_yuv, SAI->color, sizeof(unsigned short)*SAI->nr*SAI->nc * 3);

				encodeResidualJP2(SAI->nr, SAI->nc, original_color_view, tmpim_wo_yuv, ppm_residual_path,
					kdu_compress_path, jp2_residual_path_jp2, SAI->residual_rate_color, 3, offset_v, RESIDUAL_16BIT_bool);

				decodeResidualJP2(tmpim_wo_yuv, kdu_expand_path, jp2_residual_path_jp2, ppm_residual_path, ncomp1, offset_v, offset_v, RESIDUAL_16BIT_bool);

				double psnr_result_yuv_wo_trans = getYCbCr_422_PSNR(tmpim_wo_yuv, original_color_view, SAI->nr, SAI->nc, 3, 10);

				if (psnr_result_yuv_wo_trans > psnr_result_yuv_w_trans) {
					memcpy(SAI->color, tmpim_wo_yuv, sizeof(unsigned short)*SAI->nr*SAI->nc * 3);
					SAI->yuv_transform = false;
					rate_a1 = 0.0;
				}
				else {
					memcpy(SAI->color, tmpim_w_yuv, sizeof(unsigned short)*SAI->nr*SAI->nc * 3);
				}

				delete[](tmpim_wo_yuv);
				delete[](tmpim_w_yuv);

			}
			else {

				int offset_v = (1<<BIT_DEPTH) - 1;

				encodeResidualJP2(SAI->nr, SAI->nc, original_color_view, SAI->color, ppm_residual_path,
					kdu_compress_path, jp2_residual_path_jp2, SAI->residual_rate_color, 3, offset_v, RESIDUAL_16BIT_bool);

				decodeResidualJP2(SAI->color, kdu_expand_path, jp2_residual_path_jp2, ppm_residual_path, ncomp1, offset_v, offset_v, RESIDUAL_16BIT_bool);

			}

			SAI->has_color_residual = true;
		}

		aux_write16PGMPPM(SAI->path_out_ppm, SAI->nc, SAI->nr, 3, SAI->color);
		aux_write16PGMPPM(SAI->path_out_pgm, SAI->nc, SAI->nr, 1, SAI->depth);

		if (SAI->residual_rate_depth > 0 && depth_file_exist) 
		{
			//psnr_result = USE_difftest_ng ? getPSNR(NULL, SAI->path_out_pgm, SAI->path_input_pgm, difftest_call_pgm) : 0;
			output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", PSNR(SAI->depth, original_depth_view, SAI->nr, SAI->nc, 1) );
		}
		else {
			output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", 0.0);
		}


		//psnr_result = USE_difftest_ng ? getPSNR(NULL, SAI->path_out_ppm, SAI->path_input_ppm, difftest_call) : 0;
		output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", getYCbCr_422_PSNR(SAI->color, original_color_view, SAI->nr, SAI->nc, 3, BIT_DEPTH) );

		output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", rate_a1); 
		output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", SAI->stdd);
		

		delete[](original_color_view);
		delete[](original_depth_view);

		/* write view configuration data to bitstream */

		int n_bytes_prediction = 0, tmp_pred_bytes = 0;
		int n_bytes_residual = 0;

		output_LF_file = fopen(path_out_LF_data, "ab");

		if (!global_header_written) { //these are global
			n_bytes_prediction += (int)fwrite(&SAI->nr, sizeof(int), 1, output_LF_file) * sizeof(int); // needed only once per LF
			n_bytes_prediction += (int)fwrite(&SAI->nc, sizeof(int), 1, output_LF_file) * sizeof(int); // 
			n_bytes_prediction += (int)fwrite(&yuv_transform_s, sizeof(int), 1, output_LF_file) * sizeof(int); 
			n_bytes_prediction += (int)fwrite(&MINIMUM_DEPTH, sizeof(unsigned short), 1, output_LF_file) * sizeof(unsigned short);
			global_header_written = true;
		}

		viewHeaderToCodestream(n_bytes_prediction, SAI, output_LF_file, yuv_transform_s);

		/* debugging */
		tmp_codestream = fopen(path_codestream, "ab");
		viewHeaderToCodestream(tmp_pred_bytes, SAI, tmp_codestream, yuv_transform_s);
		fclose(tmp_codestream);

		if (SAI->residual_rate_color > 0) {

			if (SAI->yuv_transform && YUV_TRANSFORM) {
				for (int icomp = 0; icomp < 3; icomp++) {

					writeResidualToDisk(ycbcr_jp2_names[icomp], output_LF_file, n_bytes_residual, JP2_dict);

				}
			}
			else {

				writeResidualToDisk(jp2_residual_path_jp2, output_LF_file, n_bytes_residual, JP2_dict);
			}
		}

		if (SAI->residual_rate_depth > 0 && depth_file_exist) {

			writeResidualToDisk(jp2_residual_depth_path_jp2, output_LF_file, n_bytes_residual, JP2_dict);

		}

		fclose(output_LF_file);

		output_buffer_length += sprintf(output_results + output_buffer_length, "\t%i", n_bytes_prediction);
		output_buffer_length += sprintf(output_results + output_buffer_length, "\t%i", n_bytes_residual);

		output_results_file = fopen(output_results_filename, "a");
		fprintf(output_results_file, "%s\n", output_results);
		fclose(output_results_file);

		/* to reduce memory usage */
		if (SAI->color != NULL) {
			delete[](SAI->color);
			SAI->color = NULL;
		}

		if (SAI->depth != NULL) {
			delete[](SAI->depth);
			SAI->depth = NULL;
		}

		if (SAI->seg_vp != NULL) {
			delete[](SAI->seg_vp);
			SAI->seg_vp = NULL;
		}

		//if (SAI->segmentation != NULL) {
		//	delete[](SAI->segmentation);
		//	SAI->segmentation = NULL;
		//}

		if (SAI->label_im != NULL) {
			delete[](SAI->label_im);
			SAI->label_im = NULL;
		}

	}

	for (int ii = 0; ii < n_views_total; ii++)
	{

		//printf("ii=%d\n", ii);

		view *SAI = LF + ii;

		if (SAI->color != NULL)
			delete[](SAI->color);
		if (SAI->depth != NULL)
			delete[](SAI->depth);
		if (SAI->references != NULL)
			delete[](SAI->references);
		if (SAI->depth_references != NULL)
			delete[](SAI->depth_references);
		if (SAI->merge_weights != NULL)
			delete[](SAI->merge_weights);
		if (SAI->sparse_weights != NULL)
			delete[](SAI->sparse_weights);
		if (SAI->bmask != NULL)
			delete[](SAI->bmask);
		if (SAI->seg_vp != NULL)
			delete[](SAI->seg_vp);
		if (SAI->sparse_mask != NULL)
			delete[](SAI->sparse_mask);
		if (SAI->segmentation != NULL)
			delete[](SAI->segmentation);
		if (SAI->label_im != NULL)
			delete[](SAI->label_im);

	}

	
	for (int ii = 0; ii < maxR; ii++) {
		delete[](LF_mat[ii]);
	}

	delete[](LF_mat);

	delete[](LF);

	exit(0);
}