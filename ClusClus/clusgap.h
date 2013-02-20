#ifndef CLUSGAP_H
#define CLUSGAP_H
 
#include"clusclus.h"

class clusgap  
{
public:
	clusgap(clusclus* object, int num_trials,int numgaps);
	clusgap(clusclus* object);
	~clusgap();
	int ComputeGap();

	double** features;
	int num_samples;
	int num_features;
	int num_trials;
	int num_gaps;


private:
	void GenerateSyntheticData(int seed);
	void Stand_Devv(double *avg, double *std, double* my_vector);

	double** syntheticfeatures;
	double** dispersionmatrix ;
	class clusclus* cc;
};

#endif	//end CLUSGAP_H