#ifndef KMEANS_HH
#define KMEANS_HH

void updateValues(const int K_clusters, double *centroids, double *dvalues, int *quantized, const int npoints);
void assignClusters( const int K_clusters, double *centroids, double *dvalues, int *quantized, const int npoints );
bool updateClusters( const int K_clusters, double *centroids, double *dvalues, int *quantized, const int npoints );
void getKmeansQuantized( const int K_clusters, int *values, int npoints, const int iterations );

#endif