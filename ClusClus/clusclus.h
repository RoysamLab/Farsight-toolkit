#ifndef	CLUSCLUS_H
#define CLUSCLUS_H

#ifndef _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_DEPRECATE
#endif

#ifdef _OPENMP
#include "omp.h"
#endif

#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>


typedef struct ClusterTree
{
	ClusterTree()
	{
	};
	ClusterTree(int fir, int sec, double d, int par)
	{
		first = fir;
		second = sec;
		dis = d;
		parent = par;
	};
	int first;
	int second;
	double dis;
	int parent;
}ClusterTree;

class clusclus 
{
public:
	clusclus(int linkmode = 2);
	clusclus(double** features,int num_samples, int num_features, int linkmode = 2);

	~clusclus();

	void ReadFile(const char *filename);
	void Initialize(double** feature,int numsamples, int numfeatures);
	void Initialize( double** treedata, int numsamples);
	void RunClusClus();
	void Clustering();
	void MergersToProgress();
	void GetMembers(int num_cluster);
	void Transpose();
	void WriteClusteringOutputToFile(const char *filename1, const char *filename2,const char *filename3, 
		 const char *filename4, const char *filename5, const char *filename6, const char *filename7);
	void WriteClusteringOutputToFile(const char *filename1);
	int  ComputeGapStatistics();
	void PrepareTreeData();
	void GetOptimalLeafOrderD();
	void GetTreeStructure(std::vector< ClusterTree> &treeVec);

	int       num_samples;
	int       num_features;
	int		  num_gaps;
	int*      optimalleaforder;
	double**  gap;
	double**  mergers;
	double**  treedata;
	double**  features;
	double**  transposefeatures;
	int linkmode;

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
