
#ifndef _MCLREGRESSION_H_
#define _MCLREGRESSION_H_

#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_diag_matrix.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <iostream>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/algo/vnl_determinant.h>
//#include <vnl/algo/vnl_trace.h>
#include <mbl/mbl_stats_nd.h>
#include <vtkTable.h>
#include <vtkSmartPointer.h>
#include <vtkVariant.h>
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <iostream>
#include "ftkUtils.h"

#define MAX(a,b) (((a) > (b))?(a):(b))
#define MIN(a,b) (((a) < (b))?(a):(b))

class MCLR
{
public:
	MCLR();
	~MCLR();

	struct model{
		vnl_matrix<double> w;
		//std::string method;
		vnl_vector<int> select_index;
		double sparsity_control;
		vnl_vector<double> g;
		vnl_matrix<double> FIM;
		vnl_matrix<double> CRB;
		double kstd_level;
		vnl_matrix<double> kc_ini;
	}m;


	vnl_matrix<double> x; // Data
	vnl_matrix<double> testData;
	vnl_matrix<double> trainData;
	vnl_vector<double> y; // labels (-1 for training)
	vnl_vector<double> y_ground_truth;
	std::vector< std::pair<int,int> > id_time_val;
	vnl_matrix<double> z; // used in gradient computation
	vnl_matrix<double> gradient_w; // used in gradient computation
	vnl_matrix<double> hessian;
	vnl_matrix<double> direction;
	vnl_vector<int> class_vector;

	std::vector<int> top_features;
	std::vector<double> maxInfoVector;
	vtkSmartPointer<vtkTable> test_table;

	vnl_vector<double> diff_info_3_it;//difference in g vals for last 3 iterations
	vnl_vector<double> info_3_it;	// g vals for the last three
	vnl_vector<double> std_vec;
	vnl_vector<double> mean_vec;


	double g;
	bool stop_training;
	vnl_vector<double> stop_cond;
	double delta;
	int numberOfFeatures;
	int numberOfClasses;
	int current_label;

	double max_info;
	double confidence_threshold;

public:

	void Initialize(vnl_matrix<double> data,double c,vnl_vector<double> classes, std::string str,vtkSmartPointer<vtkTable> table);
	vnl_matrix<double> act_learn_matrix;
	vnl_matrix<double> Add_Bias(vnl_matrix<double> data);

	void Get_Gradient(vnl_matrix<double> data_with_bias);
	vnl_matrix<double> Get_Hessian(vnl_matrix<double> data_with_bias,vnl_matrix<double> w);
	void Ameliorate_Hessian_Conditions();
	double Compute_Mean_Abs_Eig(vnl_symmetric_eigensystem<double> eig);
	vnl_vector<double> Column_Order_Matrix(vnl_matrix<double> mat);
	vnl_matrix<double> Reshape_Matrix(vnl_matrix<double>mat,int r,int c );
	vnl_matrix<double> Reshape_Vector(vnl_vector<double>vec,int r,int c );
	vnl_vector<double> Newton_Direction(vnl_matrix<double> hessian_matrix,vnl_vector<double> grad_vector);
	double logit_g(double alpha,vnl_matrix<double> data_with_bias);
	double logit_stepsize();
	vnl_matrix<double> Get_F_Matrix(vnl_matrix<double> data_bias,vnl_matrix<double> w_temp);
	vnl_matrix<double> Normalize_F_Sum(vnl_matrix<double> f);
	vnl_matrix<double> Test_Current_Model(vnl_matrix<double> testData);
	vnl_matrix<double> Test_Current_Model_w(vnl_matrix<double> testData, vnl_matrix<double> m_w_matrix);
	vnl_matrix<double> GetActiveLearningMatrix(){ return m.w;};
	model Get_Training_Model();
	MCLR::model Get_Temp_Training_Model(int query,int label);
	void Update_Train_Data(std::vector< std::pair<int,int> > query_label);

	void Update_trainData(std::vector< std::pair<int,int> > query_label,bool PIA);
	vnl_matrix<double> Kron(vnl_vector<double> x,vnl_vector<double> y);
	vnl_matrix<double> Kron(vnl_matrix<double> x,vnl_vector<double> y );
	FILE* FDeclare2(char *root, char *extension, char key);
	vnl_matrix <double> tableToMatrix(vtkSmartPointer<vtkTable> table,std::vector< std::pair<int,int> > id_list);
	vnl_matrix <double> tableToMatrix_w(vtkSmartPointer<vtkTable> table);
	vnl_matrix <double> Normalize_Feature_Matrix(vnl_matrix<double> feats);
	vnl_matrix <double> Normalize_Feature_Matrix_w(vnl_matrix<double> feats, vnl_vector<double> vector_1, vnl_vector<double> vector_2);
	int Active_Query();
	//std::vector<int> ALAMO();
	std::vector<int> ALAMO(int active_query);
	std::vector<int> Get_Top_Features();
	vnl_vector<double> Get_Std_Dev(){ return std_vec; };
	vnl_vector<double> Get_Mean(){ return mean_vec; };
	void Set_Number_Of_Classes(int num){ numberOfClasses = num; };
	void Set_Number_Of_Features(int num){ numberOfFeatures = num; };
	bool MyDataSortPredicate(std::pair<int, int>& lhs, std::pair<int, int>& rhs) ;
	std::vector< std::pair< std::string, vnl_vector<double> > >act_learn_model;
	std::vector< std::pair< std::string, vnl_vector<double> > > CreateActiveLearningModel(vtkSmartPointer<vtkTable> pWizard_table);
	int GetNumberOfClasses(vtkSmartPointer<vtkTable> table);
	vnl_matrix<double> Update_Temp_Train_Data(int query,vnl_matrix<double> testDataTemp );
	std::vector<int> Submodular_AL(int activeQuery,vnl_matrix<double> testDataTemp );
	vnl_matrix<double> Update_Temp_x(int query, vnl_matrix<double> tempTestData, vnl_matrix<double> xBatchMode );
	std::vector<int> Get_Feature_Order();
	vtkSmartPointer<vtkTable> Rearrange_Table(vtkSmartPointer<vtkTable> pawTable);

};
#endif
