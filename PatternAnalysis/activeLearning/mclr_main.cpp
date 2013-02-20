



#include "mclr.h"



int Num_Lines2(char *root, int mode) {
	int   row_no, column_no, zctr1 = 0, zctr2 = 0, zctr3 = 0;
	char  file_name[80], temp[80], ch;
	FILE  *fpinn;

	strcpy(file_name, root);
	if ((fpinn = fopen(file_name, "r")) == NULL) {
		printf("\n  The file %s does not exist !\n", file_name);
		printf("      Press <ENTER> to exit     \n\n");
		getchar();
		exit(1);
	}
	while (!feof(fpinn)) {
		fscanf(fpinn, "%c",&ch);
		zctr2++;
		//--- we don't want to count empty rows
		if (isspace(ch)) zctr2--;
		//
		if (ch=='\n' && zctr2 > 0){
			zctr3++;
			zctr2=0;
		}
		if (feof(fpinn) && zctr2 > 0){//--- this part required for linux
			zctr3++;
			zctr1++;
		}
	}
	rewind(fpinn);
	while (!feof(fpinn)){ //--- counts all data points
		fscanf(fpinn, "%s",&temp);
		zctr1++;
	}
	fclose(fpinn);
	zctr1     = zctr1-1;
	row_no    = zctr3;
	column_no = zctr1/row_no;
	//printf("\n\ncolumn no = %d  row_no= %d\n\n", column_no,   row_no);
	if (mode == 0) return row_no;
	if (mode == 1) return column_no;
	return -999;
}


vnl_matrix <double> Normalize_Feature_Matrix(vnl_matrix<double> feats)
{
	mbl_stats_nd stats;

	for(int i = 0; i<feats.rows() ; ++i)
	{
		vnl_vector<double> temp_row = feats.get_row(i);
		stats.obs(temp_row);	
	}

	vnl_vector<double> std_vec = stats.sd();
	vnl_vector<double> mean_vec = stats.mean();
	
	for(int i = 0; i<feats.columns()-3 ; ++i)
	{
		vnl_vector<double> temp_col = feats.get_column(i);
		if(std_vec(i) > 0)
		{	
			for(int j =0; j<temp_col.size() ; ++j)
				temp_col[j] = (temp_col[j] - mean_vec(i))/std_vec(i) ;
		}
	
		feats.set_column(i,temp_col);
	}

	return feats;
}


vnl_matrix <double> Read_From_File(char* train_fname)
{
	//Read Training File:
	//FILE * fpin = FDeclare2(train_fname, "", 'r');
	FILE *fpin;
	fpin = fopen(train_fname,"r");

	double val;
	char str [80];

	int myrows = Num_Lines2(train_fname, 0)-1;   // number of rows            
	int mycols = Num_Lines2(train_fname, 1) ;  // number of columns

	//First column contains string names
	for(int c=0; c< mycols; ++c)
	{
		//Edit -- read from file
		fscanf(fpin, "%s",str);
	}

	vnl_matrix <double> FeatsMatrix(myrows,mycols);


	for(unsigned int r=0; r < myrows ; ++r)
	{
		for(int c=0; c< mycols; ++c)
		{
			//Edit -- read from file
			fscanf(fpin, "%lf", &val);
			FeatsMatrix.put(r,c,val);
		}
	}

	fclose(fpin);
	FeatsMatrix = Normalize_Feature_Matrix(FeatsMatrix);
	return FeatsMatrix;
}



int main(int argc , char** argv)
{

	MCLR *mclr = new MCLR();
	
	vnl_matrix<double> data;
	vnl_vector<double> classes;
	double sparsity = 1;
	int active_query = 1;
	double max_info = -1e9;
	
	int max_label_examples = atoi(argv[1]);

	//Get the training data with the classes
	vnl_matrix<double> Feats = Read_From_File("C:\\ActiveLearning\\spMCLogit_active_updated_06_02_11\\alltrain_gt2.txt");
	
	mclr->Initialize(Feats,sparsity,"ground_truth");
	mclr->Get_Training_Model();

	vnl_vector<double> diff_g_3_it(3,0);//difference in g vals for last 3 iterations
	vnl_vector<double> g_3_it(3,0);	// g vals for the last three

	

	for(int iteration = 0 ; iteration< max_label_examples ; ++ iteration)
	{
		vnl_matrix<double> test_data_just_features =  mclr->test_data.transpose();
		if(mclr->validation == "ground_truth")
		test_data_just_features = test_data_just_features.get_n_rows(0,test_data_just_features.rows()-3);
		else
		test_data_just_features = test_data_just_features.get_n_rows(0,test_data_just_features.rows()-2);

		vnl_matrix<double> prob = mclr->Test_Current_Model(test_data_just_features);
		

		vnl_matrix<double> test_data_bias = mclr->Add_Bias(test_data_just_features);
		vnl_vector<double> info_vector(test_data_just_features.cols());

		//Compute Information gain
		for(int i =0; i< test_data_just_features.cols();++i )
		{
			vnl_vector<double> temp_col_prob = prob.get_column(i);
			vnl_vector<double> temp_col_data = test_data_bias.get_column(i);		
			vnl_matrix<double> q = mclr->Kron(temp_col_prob,temp_col_data);//Kronecker Product
			vnl_diag_matrix<double> identity_matrix(mclr->no_of_classes,1); // Identity Matrix;
			double infoval = (q.transpose() * mclr->m.CRB.transpose() * q).get(0,0);
			vnl_matrix<double> info_matrix(mclr->no_of_classes,mclr->no_of_classes,infoval);
			info_matrix = identity_matrix + info_matrix ;
			info_vector(i) = log(vnl_determinant(info_matrix));
			if(info_vector(i)>max_info)
			{
				active_query = i;
				max_info = info_vector(i);
			}
		}
		
		mclr->Update_Train_Data(active_query);

		mclr->Get_Training_Model();

		std::cout<< max_info <<std::endl;

		//reset the max_info for next iteration
		max_info = -1e9;
	
	}
	


		vnl_matrix<double> test_data_just_features =  mclr->test_data.transpose();
		if(mclr->validation == "ground_truth")
		test_data_just_features = test_data_just_features.get_n_rows(0,test_data_just_features.rows()-3);
		else
		test_data_just_features = test_data_just_features.get_n_rows(0,test_data_just_features.rows()-2);
		vnl_matrix<double> currprob = mclr->Test_Current_Model(test_data_just_features);
		vnl_vector<double> gt = mclr->y_ground_truth;
		

		double accuracy = 0;

		for(int i = 0; i< currprob.cols() ; ++i)
		{
			vnl_vector<double> curr_col = currprob.get_column(i);
			if(curr_col.arg_max()+1 == gt(i))
				accuracy++;
		}
		
		double p = gt.size();

		std::cout<< accuracy/p <<std::endl;
}