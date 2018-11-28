#include "segmentation.hh"
#include "kmeans.hh"
#include "connected_components.hh"
#include "ppm.hh"

void makeKMeansSegmentation(view *SAI, const int K_MEANS_CLUSTERS, const int K_MEANS_ITERATIONS) {

	printf("Obtaining segmentation based on inverse depth at view %03d_%03d\t", SAI->c, SAI->r);

	unloadLabels(SAI);
	loadInverseDepth(SAI);

	int *k_seg_im = getKmeansQuantized(K_MEANS_CLUSTERS, SAI->depth, SAI->nr*SAI->nc, K_MEANS_ITERATIONS); /*inplace assignment to tmp_d*/

	int nregions = 0;
	int *reg_histogram = 0;
	int32_t *label_im = get_labels(k_seg_im, SAI->nr, SAI->nc, nregions, reg_histogram);

	SAI->label_im = label_im;
	SAI->nregions = nregions;
	SAI->reg_histogram = reg_histogram;

	unsigned short *labels = new unsigned short[SAI->nr*SAI->nc]();
	for (int ii = 0; ii < SAI->nr*SAI->nc; ii++) {
		*(labels + ii) = static_cast<unsigned short>(*(label_im + ii));
	}

	char labels_file[1024];
	sprintf(labels_file, "%s%c%03d_%03d%s", SAI->output_dir, '/', SAI->c, SAI->r, "_labels.pgm");
	aux_write16PGMPPM(labels_file, SAI->nc, SAI->nr, 1, labels);
	delete[](labels);

	char tmp_d_file[1024];
	sprintf(tmp_d_file, "%s%c%03d_%03d%s", SAI->output_dir, '/', SAI->c, SAI->r, "_kmeans_disparity.int32");
	FILE *tmpfile_d;
	tmpfile_d = fopen(tmp_d_file, "wb");
	fwrite(k_seg_im, sizeof(int), SAI->nr*SAI->nc, tmpfile_d);
	fclose(tmpfile_d);
	delete[](k_seg_im);

	FILE *tmpfile_im_labels;
	tmpfile_im_labels = fopen(SAI->path_label_im, "wb");
	fwrite(SAI->label_im, sizeof(int), SAI->nr*SAI->nc, tmpfile_im_labels);
	fclose(tmpfile_im_labels);

	unloadLabels(SAI);
	unloadInverseDepth(SAI);

}