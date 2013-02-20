#ifndef _ftkGTClustering_h_
#define _ftkGTClustering_h_

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_random.h>
#include <vcl_complex.h>
#include <mbl/mbl_stats_nd.h>

#include <vector>

class ftkGTClustering{

public:

	ftkGTClustering();
	~ftkGTClustering();

	void ComputeCostMatrix(); 
	void NormalizeFeatures();
	void RunReplicatorDynamics();
	void ApplySurvivalThreshold();

	void set_featureMatrix(vnl_matrix<double> feature_matrix);
	void set_repMaxIter(int max_iter);
	void set_repStoppingThr(double stopping_thr);
	void set_repSurvivalThr(double survival_thr);
	void set_repAlpha(double alpha);
	vnl_vector<double> get_repFinalPopulation();
	std::vector<bool> get_repEvolvedStrategiesFlag();

private:
	
	vnl_matrix<double> costMatrix;
	vnl_matrix<double> featureMatrix;	
	vnl_vector<double> repFinalPopulation;
	std::vector<bool> repEvolvedStrategiesFlag;

	int repMaxIter;
	double repStoppingThr;
	double repAlpha;
	double repSurvivalThr;
	std::string exitMsg;
};

#endif 
