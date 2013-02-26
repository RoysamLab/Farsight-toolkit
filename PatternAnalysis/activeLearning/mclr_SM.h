#ifndef _MCLREGRESSION_SM_H_
#define _MCLREGRESSION_SM_H_

#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_diag_matrix.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <iostream>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/algo/vnl_determinant.h>
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

#define MAX(a,b) (((a) > (b))?(a):(b))
#define MIN(a,b) (((a) < (b))?(a):(b))

class MCLR_SM 
{
public:
	MCLR_SM();	
	~MCLR_SM();

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
	vnl_matrix<double> training_data;
	vnl_matrix<double> test_data;
	vnl_matrix<double> train_data;
	vnl_vector<double> y; // labels (-1 for training)	
	vnl_vector<double> y_ground_truth;	
	std::vector< std::pair<int,int> > id_time_val;	
	vnl_matrix<double> z; // used in gradient computation	
	vnl_matrix<double> gradient_w; // used in gradient computation	 
	vnl_matrix<double> hessian;
	vnl_matrix<double> direction;  	
	vnl_vector<int> class_vector;

	std::vector<int> top_features;
	std::vector<double> max_info_vector;
	vtkSmartPointer<vtkTable> test_table;

	vnl_vector<double> diff_info_3_it;//difference in g vals for last 3 iterations
	vnl_vector<double> info_3_it;	// g vals for the last three
	vnl_vector<double> std_vec;
	vnl_vector<double> mean_vec;


	double g;	
	bool stop_training;
	vnl_vector<double> stop_cond;	
	double delta; 	
	int no_of_features;
	int no_of_classes;
	int current_label;

	double max_info;

	std::string validation;
	double confidence_threshold;

public:

	//void Initialize(vnl_matrix<double> data,double c);
	void Initialize(vnl_matrix<double> data,double c,vnl_vector<double> classes, std::string str,vtkSmartPointer<vtkTable> table );

	vnl_matrix<double> Add_Bias(vnl_matrix<double> data);

	void Get_Gradient(vnl_matrix<double> data_with_bias);
	void Get_Hessian(vnl_matrix<double> data_with_bias);
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
	vnl_matrix<double> Test_Current_Model(vnl_matrix<double> test_data);
	vnl_matrix<double> Test_Current_Model_1(vnl_matrix<double> test_data, vnl_matrix<double> m_w_matrix);
	vnl_matrix<double> GetActiveLearningMatrix(){ return m.w;};
	//vnl_matrix <double> Normalize_Feature_Matrix(vnl_matrix<double> feats);
	model Get_Training_Model();
	void Update_Train_Data(int query,int label);
	vnl_matrix<double> Kron(vnl_vector<double> x,vnl_vector<double> y);
	void Get_Label_Sample(int query);
	FILE* FDeclare2(char *root, char *extension, char key);
	vnl_matrix <double> tableToMatrix(vtkSmartPointer<vtkTable> table,std::vector< std::pair<int,int> > id_list);
	vnl_matrix <double> tableToMatrix_1(vtkSmartPointer<vtkTable> table);
	vnl_matrix <double> Normalize_Feature_Matrix(vnl_matrix<double> feats);
	vnl_matrix <double> Normalize_Feature_Matrix_1(vnl_matrix<double> feats, vnl_vector<double> vector_1, vnl_vector<double> vector_2);
	int Active_Query();
	std::vector<int> Get_Top_Features();
	vnl_vector<double> Get_Std_Dev(){ return std_vec; };
	vnl_vector<double> Get_Mean(){ return mean_vec; };
	void Set_Number_Of_Classes(int num){ no_of_classes = num; };
	void Set_Number_Of_Features(int num){ no_of_features = num; };

};
#endif
