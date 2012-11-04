/*=========================================================================
Copyright 2009 Rensselaer Polytechnic Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 
=========================================================================*/

#ifndef __fuzzy_clustering_h
#define __fuzzy_clustering_h

#include <iostream>
#include <vector>

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
	//added by Yousef on 10/13/2009
	int Validate(); 
private:	
	void SortMembershipValues(double* vals, int* ids, int l);
	//added by Yousef on 10/13/2009
	void rearrangeLabels(int* tmpLabels, int* arrangement);
	void rotate(int *c,int start, int len);
	void rotrec(int start,int lenth,int len, int *c, int** combs, int* ll);
	void NormalizeFeatures();
	//
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
	//added by Yousef on 10/13/2009
	int* labels;
};
#endif

