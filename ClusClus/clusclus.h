#ifndef	CLUSCLUS_H
#define CLUSCLUS_H

#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

using namespace std;

class clusclus 
{
public:
	clusclus();
	clusclus(const char* filename);
	clusclus(double** features,int num_samples, int num_features);

	~clusclus();

	void ReadFile(const char *filename);
	void RunClusClus();
	void MergersToProgress();
	void GetMembers(int num_cluster);
	void Transpose();
	void WriteClusteringOutputToFile(const char *filename1, const char *filename2,const char *filename3, 
		 const char *filename4, const char *filename5, const char *filename6, const char *filename7);
	void WriteClusteringOutputToFile(const char *filename1);
	int  ComputeGapStatistics();
	void PrepareTreeData();
	void GetOptimalLeafOrderD();

	int       num_features;
	int       num_samples;
	int       linkmode;
	int		  num_gaps;
	int*      optimalleaforder;
	double**  features;
	double**  transposefeatures;
	double**  gap;
	double**  mergers;
	double**  treedata;

private:
//	void   NormalizeFeatures();
	void   ComputeSampleDistances();	
	double ComputeDistance(double* vector1, double* vector2);
	double MergeClusters(int num_currcluster, int *pivot1, int *pivot2);
	double UpdateClusterDistances(int num_currcluster, int pivot1, int pivot2);
	int    FindRowIndex(int col, int c);
	void   GetKids(int k, int* Tnums, int* pickedup, int* kids);

	int*     num_cluster_samples;
	int**    progress;
	int**    members;
	double*  sample_distances;
	double*  cluster_distances;
};

#endif	//end CLUSCLUS_H
