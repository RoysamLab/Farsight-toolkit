#ifndef KMEANS_H
#define KMEANS_H 

#include <set>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include "vnl/vnl_matrix.h"

class Kmeans
{
public:
	Kmeans();
	~Kmeans();
	void setDataToKmeans(std::vector<std::vector<double > > & points);
	void setDataToKmeans(vnl_matrix<double> points);
	void setDistasnce(char distancestring);
	void setClusterNumber(int clusternumber);
	void Clustering();
	void WriteFile(const char *filename1, const char *filename2);
	std::vector<int > getPointsToClusters();
	std::vector<std::set<int > > getClustersToPoints();
	std::vector<std::vector<int > > getClustersToPointsvector();
	std::vector<std::vector<double > > getClusterCenters();

	int num_clusters;
	int num_iterations;
	int distance_mode;	

private:
	int num_rows;
	int num_cols;

	std::vector<std::vector<double > > points;
	std::vector<std::vector<double > > centroids;
	std::vector<std::set<int > > clustersToPoints;
	std::vector<std::vector<int > > clustersToPointsvector;
	std::vector<int > pointsToClusters;

	void initialClusters();
	void initialCentroids();
	void updateCentroids();
	void updateCentroid(std::set<int > & clusterToPoints, std::vector<double > & point);
	void updateClusters(bool & move);
	void updatePointIdInClusters();
	double computeDistance(std::vector<double > & point1,std::vector<double > & point2);
};

#endif