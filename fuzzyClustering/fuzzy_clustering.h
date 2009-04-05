#ifndef __fuzzy_clustering_h
#define __fuzzy_clustering_h

#include <iostream>

using namespace std;

class FuzzyClustering
{
public:
	FuzzyClustering(int numC);
	FuzzyClustering(int numC, int numIter, int prec);
	~FuzzyClustering();
	void InitializeCenters();
	void Run();
	void RunOneStep();
	void ReadProblem(const char *filename);
	void ComputeDistances();
	void ComputeMembershipValues();
	void FindNewCenters();
	double ComputeCurrentError();
	double** GetFuzzyMemberships() { return Memberships; }
	double** ExtractTrainingSet(double C); 
	void AssignClusters();
	int* GetClusteringOutput() { return clustering_Output; }
	void WriteClusteringOutputToFile(const char *filename);
	void WriteTrainingSetToFile(const char *filename);
private:	
	void SortMembershipValues(double* vals, int* ids, int l);
	int num_Clusters;
	int num_Iterations;
	double precision;
	int num_Features;
	int num_Samples;
	double weighting_exp;
	double **Features;
	double **Centers;
	double **new_Centers;
	double **Memberships;
	double *sum_Memberships;
	double **Distances;
	double *sum_Distances;
	int **Singularities;
	double** trainingSet;
	int* clustering_Output;
	int* cluster_Distributions;
	int trainingSetSize;
};
#endif