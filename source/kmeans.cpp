#include "kmeans.hh"
#include <vector>
#include <cmath>

void assignClusters(const int K_clusters, double *centroids, double *dvalues, int *quantized, const int npoints) {

	for (int ii = 0; ii < npoints; ii++) {
		double minval = DBL_MAX;
		int cluster_id = -1;
		for (int jj = 0; jj < K_clusters; jj++) {
			double distance = dvalues[ii] - centroids[jj];
			distance = distance*distance;
			if (distance < minval) {
				minval = distance;
				cluster_id = jj;
			}
		}
		quantized[ii] = cluster_id;
		}
}

void updateValues(const int K_clusters, double *centroids, double *dvalues, int *quantized, const int npoints) {

	for (int ii = 0; ii < npoints; ii++) {
		dvalues[ii] = centroids[quantized[ii]];
	}
}

bool updateClusters(const int K_clusters, double *centroids, double *dvalues, int *quantized, const int npoints) {
	
	bool converge = true;

	for (int jj = 0; jj < K_clusters; jj++) {
		double sum = 0;
		int n = 0;
		for (int ii = 0; ii < npoints; ii++) {
			if (quantized[ii] == jj) {
				sum += dvalues[ii];
				n++;
			}
		}
		
		if (n > 0) {
			double new_centroid = floor( (sum / static_cast<double>(n)) + 0.5 );

			if ( abs(centroids[jj] - new_centroid) > 0.1 ) {
				converge = false;
			}

			printf("Centroid id %i\told\t%f\tnew\t%f\tn: %i\n", jj, centroids[jj], new_centroid,n);

			centroids[jj] = new_centroid;
		}

	}

	return converge;
}

void getKmeansQuantized(const int K_clusters, int *values, const int npoints, const int iterations) {

	double min_value = DBL_MAX;
	double max_value = DBL_MIN;

	double *dvalues = new double[npoints]();

	for (int ii = 0; ii < npoints; ii++) {
		min_value = values[ii] < min_value ? static_cast<double>( values[ii] ) : min_value;
		max_value = values[ii] > max_value ? static_cast<double>( values[ii] ) : max_value;

		dvalues[ii] = static_cast<double>( values[ii] );
	}

	double delta = (max_value - min_value) / static_cast<double>( K_clusters );

	double *centroids = new double[K_clusters]();
	for (int jj = 0; jj < K_clusters; jj++) {
		centroids[jj] = floor( min_value + jj*delta + 0.5 );
	}

	int n_iter = 1;

	int *quantized = new int[npoints]();

	bool convergence = false;

	while ( !convergence && n_iter<iterations ) {

		printf("K-Means iteration:\t%i\n", n_iter);

		assignClusters(K_clusters, centroids, dvalues, quantized, npoints);
		convergence = updateClusters(K_clusters, centroids, dvalues, quantized, npoints);

		n_iter++;
	}

	updateValues(K_clusters, centroids, dvalues, quantized, npoints);

	for (int ii = 0; ii < npoints; ii++) {
		values[ii] = static_cast<int>( floor(dvalues[ii] + 0.5 ) );
	}
	
	delete[](dvalues);
	delete[](centroids);
	delete[](quantized);
}