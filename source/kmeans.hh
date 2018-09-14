#ifndef KMEANS_HH
#define KMEANS_HH

#define K_MEANS_CLUSTERS 16
#define K_MEANS_ITERATIONS 16

void updateValues(const int K_clusters, double *centroids, double *dvalues, int *quantized, const int npoints);
void assignClusters( const int K_clusters, double *centroids, double *dvalues, int *quantized, const int npoints );
bool updateClusters( const int K_clusters, double *centroids, double *dvalues, int *quantized, const int npoints );
//void getKmeansQuantized( const int K_clusters, int *values, int npoints, const int iterations );

template <class T>
void getKmeansQuantized(const int K_clusters, T *values, const int npoints, const int iterations) {

	double min_value = DBL_MAX;
	double max_value = -DBL_MAX;

	double *dvalues = new double[npoints]();

	for (int ii = 0; ii < npoints; ii++) {
		min_value = values[ii] < min_value ? static_cast<double>(values[ii]) : min_value;
		max_value = values[ii] > max_value ? static_cast<double>(values[ii]) : max_value;

		dvalues[ii] = static_cast<double>(values[ii]);
	}

	double delta = (max_value - min_value) / static_cast<double>(K_clusters);

	double *centroids = new double[K_clusters]();
	for (int jj = 0; jj < K_clusters; jj++) {
		centroids[jj] = floor(min_value + jj*delta + 0.5);
	}

	int n_iter = 1;

	int *quantized = new int[npoints]();

	bool convergence = false;

	while (!convergence && n_iter <= iterations) {

		//printf("K-Means iteration:\t%i\n", n_iter);

		assignClusters(K_clusters, centroids, dvalues, quantized, npoints);
		convergence = updateClusters(K_clusters, centroids, dvalues, quantized, npoints);

		n_iter++;
	}

	updateValues(K_clusters, centroids, dvalues, quantized, npoints);

	for (int ii = 0; ii < npoints; ii++) {
		values[ii] = static_cast<T>(dvalues[ii]);
		//values[ii] = static_cast<T>(floor(dvalues[ii] + 0.5));
	}

	delete[](dvalues);
	delete[](centroids);
	delete[](quantized);
}

#endif