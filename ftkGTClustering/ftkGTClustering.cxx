
#include "ftkGTClustering.h"

ftkGTClustering::ftkGTClustering(){

	this->repMaxIter = 10000; //5000;
	this->repStoppingThr = 10^-16;
	this->repAlpha = 0.01;
	this->repSurvivalThr = 0.25;
}

ftkGTClustering::~ftkGTClustering(){
}

void ftkGTClustering::ComputeCostMatrix(){

	//feature_mat.print(std::cout);
	//std::cout << "Alpha: " << this->repAlpha << std::endl;
	this->costMatrix.clear();
	this->costMatrix.set_size(this->featureMatrix.rows(), this->featureMatrix.rows());
	this->costMatrix.fill(0);
	
	for(unsigned int i = 0; i < this->costMatrix.rows(); i++){
		for(unsigned int j = 0; j < this->costMatrix.columns(); j++){						
			this->costMatrix(i, j) = vcl_exp(-1*this->repAlpha*
				vcl_sqrt((this->featureMatrix.get_row(i)-this->featureMatrix.get_row(j)).squared_magnitude()));  
		}
	}
	this->costMatrix.fill_diagonal(0);
}

void ftkGTClustering::NormalizeFeatures(){

	mbl_stats_nd stats;

	for(unsigned int i = 0; i< this->featureMatrix.rows() ; ++i)
	{
		vnl_vector<double> temp_row = this->featureMatrix.get_row(i);
		stats.obs(temp_row);	
	}

	vnl_vector<double> std_vec = stats.sd();
	vnl_vector<double> mean_vec = stats.mean();

	for(unsigned int i = 0; i < this->featureMatrix.columns() ; ++i){

		vnl_vector<double> temp_col = this->featureMatrix.get_column(i);
		if(std_vec(i) > 0)
		{	
			for(unsigned int j =0; j < temp_col.size() ; ++j)
				temp_col[j] = (temp_col[j] - mean_vec(i))/std_vec(i) ;
		}
		this->featureMatrix.set_column(i,temp_col);
	}
}

void ftkGTClustering::RunReplicatorDynamics(){

	int num_strategies = this->costMatrix.rows();
	vnl_vector<double> initial_strategies(num_strategies, (double)1.0/num_strategies);
	vnl_random random_gen; 	
	
	for(unsigned int i = 0; i < initial_strategies.size(); i++)
		initial_strategies(i) += random_gen.drand64();
	initial_strategies.normalize();
	
	vnl_vector<double> next_strategies = initial_strategies;
	vnl_vector<double> current_strategies = initial_strategies;

	int iter = 0;
	while(true){

		vnl_vector<double> Cx = this->costMatrix*current_strategies;

		double xTCx = dot_product(current_strategies, Cx);
		next_strategies = element_product(current_strategies, Cx);
		next_strategies /= xTCx;
		
		double evolution_metric = (current_strategies - next_strategies).one_norm();
		if(evolution_metric < this->repStoppingThr){
			this->exitMsg = std::string("EXIT: Evolution speed is below stoppingThreshold. ");
			break;
		}

		current_strategies = next_strategies;
		iter++;

		if(iter > this->repMaxIter){
			this->exitMsg = std::string("EXIT: Max iterations reached. ");
			break;
		}
	}
	this->repFinalPopulation = current_strategies;
}

void ftkGTClustering::ApplySurvivalThreshold(){

	this->repEvolvedStrategiesFlag.clear();
	this->repEvolvedStrategiesFlag.resize(this->repFinalPopulation.size(), false);
	double max_pop = this->repFinalPopulation.max_value();
	
	for(unsigned int i = 0; i < this->repFinalPopulation.size(); i++){
		if(this->repFinalPopulation[i] >= this->repSurvivalThr*max_pop)
			this->repEvolvedStrategiesFlag[i] = true;
	}
}

vnl_vector<double> ftkGTClustering::get_repFinalPopulation(){
	return this->repFinalPopulation;
}

void ftkGTClustering::set_featureMatrix(vnl_matrix<double> feature_matrix){
	this->featureMatrix = feature_matrix;
}

void ftkGTClustering::set_repMaxIter(int max_iter){
	this->repMaxIter = max_iter;
}

void ftkGTClustering::set_repStoppingThr(double stopping_thr){
	this->repStoppingThr = stopping_thr;
}

void ftkGTClustering::set_repAlpha(double alpha){
	this->repAlpha = alpha;
}

void ftkGTClustering::set_repSurvivalThr(double survival_thr){
	this->repSurvivalThr = survival_thr;
}

std::vector<bool> ftkGTClustering::get_repEvolvedStrategiesFlag(){
	return this->repEvolvedStrategiesFlag;
}