#include "mclr.h"

#ifdef _OPENMP
#include "omp.h"
#endif

MCLR::MCLR()
{
	current_label = -1; // keeps a track of the current label
	confidence_threshold = 0.7;
	model k;
}

MCLR::~MCLR()
{

}


struct sort_pred {
	bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right) {
		return left.first > right.first;
	}
};



// x contains the training features
// y contains the training labels
void MCLR::Initialize(vnl_matrix<double> data,double c,vnl_vector<double> classes, std::string str,vtkSmartPointer<vtkTable> table,bool PIA )
{
	int train_counter =0;
	int test_counter = 0;
	validation = str;
	//bool PIA =false;

	// test_table contains the table after the queried samples are removed from the original table
	// It is updated in the Update_Train_data 
	test_table = table; //Features->Columns ; Samples -> Rows

	for(int i = 0; i< classes.size() ; ++i)
	{
		if(classes(i) == -1)
			test_counter++;
	}

	train_data.set_size(data.rows()-test_counter,data.cols());
	test_data.set_size(test_counter,data.cols());

	test_counter = 0;
	//variable to keep track of the number of rows deleted
	int num_del_rows = 0;
	for(int i= 0;i<data.rows();++i)
	{
		if(classes(i)!=-1)
		{	
			train_data.set_row(train_counter,data.get_row(i));
			++train_counter;
			test_table->RemoveRow(i - num_del_rows);
			++num_del_rows;
		}
		else
		{
			test_data.set_row(test_counter,data.get_row(i));
			++test_counter;
		}
	}

	if(PIA)
		y.set_size(train_data.rows());

	// Transpose : Features->Rows ; Samples -> Columns
	x = train_data.transpose();// Feature matrix of training samples only ; does not contain unlabeled sample data

	int counter = 0;
	for(int i = 0; i < classes.size(); ++i)
	{
		if(classes(i)!=-1)
		{
			y(counter) = classes(i);// class value ; 
			counter++;
		}
	}

	//m.method = "direct"; // No kernel is the default
	m.sparsity_control = c; // greater the value , more is the sparsity
	m.kstd_level = 1; // kernel sigma val

	// Stop Condition
	// Used in active label selection
	stop_cond.set_size(2);
	stop_cond[0] = 1e-5;
	stop_cond[1] = 0.1;
	delta = 1e-9; // used for approximating |w| 

	no_of_features = x.rows();
	// We assume that the classes are numbered from 1 to "n" for n classes
	// so y.maxval - .minval +1 gives number of classes
	if(PIA)
	{
		no_of_classes = y.max_value()-y.min_value()+1;
		class_vector.set_size(no_of_classes);

		m.w.set_size(no_of_features+1,no_of_classes);
		m.w.fill(0);

		gradient_w.set_size(no_of_features+1,no_of_classes);
		hessian.set_size((no_of_features+1)*(no_of_classes),(no_of_features+1)*(no_of_classes)); /// Check the size again!!!

		// g :the objective function value
		g = 0;
		// Class vector
		for(int i=1;i<=no_of_classes;++i)
		{
			class_vector.put(i-1,i); 
		}

		diff_info_3_it.set_size(3);
		info_3_it.set_size(3);
	}
}


//Converts vtkTable to Vnl_Matrix
vnl_matrix <double> MCLR::tableToMatrix(vtkSmartPointer<vtkTable> table,std::vector< std::pair<int,int> > id_time)
{

	vnl_matrix <double> FeatsMatrix(table->GetNumberOfRows(),table->GetNumberOfColumns());

	//extract data from the model and get min/max values:
	for(unsigned int r=0; r < table->GetNumberOfRows() ; ++r)
	{
		for(int c=0; c< table->GetNumberOfColumns(); ++c)
		{
			FeatsMatrix.put(r,c,table->GetValue(r,c).ToDouble());
		}
	}

	// Contains the ids and time of nuclei.
	this->id_time_val = id_time;
	return FeatsMatrix;
}

vnl_matrix <double> MCLR::tableToMatrix_w(vtkSmartPointer<vtkTable> table)
{

	vnl_matrix <double> FeatsMatrix(table->GetNumberOfRows(),table->GetNumberOfColumns());
	vnl_matrix <double> FeatsMatrix_no_id(table->GetNumberOfRows(),table->GetNumberOfColumns()-1);

	//extract data from the model and get min/max values:
	for(unsigned int r=0; r < table->GetNumberOfRows() ; ++r)
	{
		for(int c=0; c< table->GetNumberOfColumns(); ++c)
		{
			FeatsMatrix.put(r,c,table->GetValue(r,c).ToDouble());
		}
	}

	return FeatsMatrix;
}


vnl_matrix <double> MCLR::Normalize_Feature_Matrix(vnl_matrix<double> feats)
{
	mbl_stats_nd stats;

	for(int i = 0; i<feats.rows() ; ++i)
	{
		vnl_vector<double> temp_row = feats.get_row(i);
		stats.obs(temp_row);	
	}

	std_vec = stats.sd();
	mean_vec = stats.mean();

	//The last column is the training column 
	for(int i = 0; i<feats.columns() ; ++i)
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

vnl_matrix <double> MCLR::Normalize_Feature_Matrix_w(vnl_matrix<double> feats, vnl_vector<double> vector_1, vnl_vector<double> vector_2)
{
	std_vec = vector_1;
	mean_vec = vector_2;


	//The last column is the training column 
	for(int i = 0; i<feats.columns() ; ++i)
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


vnl_matrix<double> MCLR::Get_F_Matrix(vnl_matrix<double> data_bias,vnl_matrix<double> w_temp)
{
	vnl_matrix<double> epow_matrix = w_temp.transpose()*data_bias;
	vnl_matrix<double> temp_f;
	temp_f.set_size(epow_matrix.rows(),epow_matrix.cols());
	for(int i=0;i<epow_matrix.rows();++i)
	{
		for(int j=0;j<epow_matrix.cols();++j)
		{
			temp_f(i,j) = exp(epow_matrix(i,j));
		}
	}
	return temp_f;
}


//f  = f./repmat(sum(f,1),[classN,1]);
vnl_matrix<double> MCLR::Normalize_F_Sum(vnl_matrix<double> f)
{
	// Matrix for normalization
	vnl_matrix<double> norm_matrix(f.rows(),f.cols());
	vnl_vector<double> norm_matrix_row;
	norm_matrix_row.set_size(f.cols());

		//repmat(sum(f,1),[classN,1]);
	//#pragma omp parallel for num_threads(16) 
	for(int i=0;i<f.cols();++i)
	{
		double sum = 0 ;  
		for(int j=0;j<no_of_classes;++j)
		{
			sum = sum + f(j,i);
		}
		norm_matrix_row(i) = sum;
	}


	for(int i=0;i<no_of_classes;++i)
	{
		norm_matrix.set_row(i,norm_matrix_row); 
	}
	
		// f  = f./repmat(sum(f,1),[classN,1]);
	//#pragma omp parallel for num_threads(16)
	for(int i=0;i<f.rows();++i)
	{
		for(int j=0;j<f.cols();++j)
		{
			f(i,j) = f(i,j)/norm_matrix(i,j);
		}
	}

	return f;
}



//	// Matrix for normalization
//	vnl_matrix<double> norm_matrix;
//	vnl_vector<double> norm_matrix_row;
//	norm_matrix_row.set_size(f.cols());
//	
//	 //repmat(sum(f,1),[classN,1]);
//
//	for(int i=0;i<f.cols();++i)
//	{
//      double sum = 0 ;  
//	  for(int j=0;j<no_of_classes;++i)
//		{
//			sum = sum + f(j,i);
//		}
//		norm_matrix_row(i) = sum;
//	}
//
//    for(int i=0;i<no_of_classes;++i)
//	{
//		norm_matrix.set_row(i,norm_matrix_row); 
//	}
//
//// f  = f./repmat(sum(f,1),[classN,1]);
//	for(int i=0;i<f.rows();++i)
//	{
//	  for(int j=0;j<f.cols();++j)
//	   {
//		  f(i,j) = f(i,j)/norm_matrix(i,j);
//	  }
//	}


void MCLR::Get_Gradient(vnl_matrix<double> data_with_bias)
{
	// f  = exp(w'*[ones(1,N);basis]); In MATLAB 
	vnl_matrix<double> f;
	f = Get_F_Matrix(data_with_bias,m.w);
	f = Normalize_F_Sum(f);

	//   grad_w = grad_w - C*(w./sqrt(w.^2+delta));
	// |w| approximation

	vnl_matrix<double> abs_w; //
	abs_w.set_size(no_of_features+1,no_of_classes);

	//#pragma omp parallel for num_threads(16) 
	for(int i=0;i<no_of_features+1;++i)
	{
		for(int j=0;j<no_of_classes;++j)
		{
			abs_w(i,j) = m.w(i,j)/ sqrt((m.w(i,j))*(m.w(i,j)) + delta);
		}
	}

	//Gradient computation
	gradient_w = data_with_bias*(z-f).transpose();	

	//Gradient with sparsity
	gradient_w = gradient_w - m.sparsity_control*abs_w;

}


vnl_matrix<double> MCLR::Add_Bias(vnl_matrix<double> data)
{
	vnl_matrix<double> data_with_bias; // Data to be modified

	vnl_vector<double> temp_vector;
	temp_vector.set_size(data.cols());
	temp_vector.fill(1);

	data_with_bias.set_size(no_of_features+1,data.cols());
	data_with_bias.set_row(0,temp_vector);
	for(int i =0; i<no_of_features ; ++i)
	{
		data_with_bias.set_row(i+1,data.get_row(i));	
	}

	return data_with_bias;
}


void MCLR::Get_Hessian(vnl_matrix<double> data_with_bias)
{	
	vnl_matrix<double> f;
	f = Get_F_Matrix(data_with_bias,m.w);
	f = Normalize_F_Sum(f);

	vnl_matrix<double> temp_hessian; // temporary Hessian Matrix
	temp_hessian.set_size((no_of_features+1)*(no_of_classes),(no_of_features+1)*(no_of_classes));

	for(int i = 1; i<=no_of_classes ;++i)
	{
		vnl_vector<double> ith_row  = f.get_row(i-1);
		
		//#pragma omp parallel for num_threads(16) 
		for(int j = 0; j< x.cols() ; ++j)
		{
			ith_row(j) = ith_row(j)*(1-ith_row(j));
		}


		vnl_diag_matrix<double> diagonal_matrix(ith_row);

		vnl_matrix<double> temp_matrix =  -data_with_bias*diagonal_matrix*data_with_bias.transpose();
		temp_hessian.update(temp_matrix,(i-1)*(no_of_features+1),(i-1)*(no_of_features+1));


		vnl_vector<int> all_class_but_i;
		all_class_but_i.set_size(class_vector.size()-1);
		int counter = 0;

		for(int k = 0; k< class_vector.size() ; ++k)
		{	
			if(class_vector(k)!=i)
			{			
				all_class_but_i.put(counter,class_vector(k));
				counter ++;
			}
		}	


		for(int k = 0; k< all_class_but_i.size() ; ++k)
		{	
			ith_row  = f.get_row(i-1);
			vnl_vector<double> kth_row  = f.get_row(all_class_but_i(k)-1);
			
		//#pragma omp parallel for num_threads(16) 
			for(int j = 0; j< x.cols() ; ++j)
			{
				//f(i,:).*(-f(k,:))
				ith_row(j) = ith_row(j)*(-kth_row(j));
			}

			vnl_diag_matrix<double> diagonal_matrix(ith_row);
			vnl_matrix<double> temp_matrix = -data_with_bias*diagonal_matrix*data_with_bias.transpose();
			temp_hessian.update(temp_matrix,(i-1)*(no_of_features+1),(all_class_but_i(k)-1)*(no_of_features+1));

		}
	}
	//hessian = hessian - diag(C*(delta./(sqrt(w(:).^2+delta)).^3));	
	//Need to vectorize w 

	vnl_vector<double> w_column_ordered = Column_Order_Matrix(m.w);

	int count = 0;
	for(int i=0;i<m.w.rows();++i)
	{
		//#pragma omp parallel for num_threads(16) 
		for(int j=0;j<m.w.cols();++j)	
		{
			//w_column_ordered(count) = (m.sparsity_control)*(delta/pow(sqrt(pow(w_column_ordered(count),2)+delta),3));
			w_column_ordered(j+i*m.w.cols()) = (m.sparsity_control)*(delta/pow(sqrt(pow(w_column_ordered(j+i*m.w.cols()),2)+delta),3));
			//count++;
		}
	}

	vnl_diag_matrix<double> diagonal_matrix(w_column_ordered);
	temp_hessian = temp_hessian - diagonal_matrix;
	//temp_hessian.extract(hessian,no_of_features+2,no_of_features+2);
	hessian = temp_hessian;

}	


void MCLR::Ameliorate_Hessian_Conditions()
{	

	//std::cout<< "Into Amel Hessian" << std::endl;
	//Compute the eigen vectors of the hessian matrix
	vnl_symmetric_eigensystem<double> eig(hessian);

	double ratio = fabs(eig.get_eigenvalue(0))/fabs(eig.get_eigenvalue((no_of_features+1)*no_of_classes)-1);
	int counter = 0;

	if(ratio>1e9)
	{
		vnl_diag_matrix<double> diag((no_of_classes)*(no_of_features+1),Compute_Mean_Abs_Eig(eig));
		hessian = hessian + diag;
		//vnl_symmetric_eigensystem<double> eig(hessian);
		//ratio = fabs(eig.get_eigenvalue(0))/fabs(eig.get_eigenvalue((no_of_features+1)*no_of_classes)-1);
	}
}


double MCLR::Compute_Mean_Abs_Eig(vnl_symmetric_eigensystem<double> eig)
{
	double sum = 0; 	
	for(int i =0; i<(no_of_features+1)*(no_of_classes);++i)
	{
		sum = sum + fabs(eig.get_eigenvalue(i));
	}	
	sum = sum/((no_of_features+1)*(no_of_classes)) ;
	return sum;
}


vnl_vector<double> MCLR::Column_Order_Matrix(vnl_matrix<double> mat)
{
	vnl_vector<double> mat_column_ordered;
	mat_column_ordered.set_size(mat.rows()*mat.cols());
	int count = 0;

	for(int j=0;j<mat.cols();++j)	
	{
		for(int i=0;i<mat.rows();++i)
		{
			mat_column_ordered(count) = mat(i,j);
			count++;
		}
	}
	return mat_column_ordered;
}



//Reshape the matrix : columns first ; Similar to MATLAB
vnl_matrix<double> MCLR::Reshape_Matrix(vnl_matrix<double>mat,int r,int c )
{
	if(mat.rows()*mat.cols() != r*c)
	{
		cout<< "Number of elements in the matrix/vector should be equal to the total number of elements in the reshaped matrix";
		getchar();
		exit(1);
	}

	vnl_matrix<double>reshaped_matrix;
	reshaped_matrix.set_size(r,c);
	int count = 0;

	for(int j=0;j<c;++j)	
	{
		for(int i=0;i<r;++i)
		{
			reshaped_matrix(i,j) = mat(count%mat.rows(),floor(static_cast<double>(count/mat.rows())));
			count++;
		}
	}
	return reshaped_matrix;
}


//Reshape the vector into a matrix : Similar to MATLAB
vnl_matrix<double> MCLR::Reshape_Vector(vnl_vector<double>vec,int r,int c )
{
	if(vec.size() != r*c)
	{
		std::cerr << "Number of elements in the matrix/vector should be equal to the total number of elements in the reshaped matrix";
		exit(1);
	}

	vnl_matrix<double>reshaped_matrix;
	reshaped_matrix.set_size(r,c);
	int count = 0;

	for(int j=0;j<c;++j)	
	{
		for(int i=0;i<r;++i)
		{
			reshaped_matrix(i,j) = vec(count);
			count++;
		}
	}
	return reshaped_matrix;
}


double MCLR::logit_stepsize()
{
	vnl_vector<double> g_temp(3,0);

	double alpha = 0;
	g_temp(0) = logit_g(alpha,Add_Bias(x));
	double a= 0;
	double	b =0  ; 
	double c = 0;
	double stepsize;
	double g_alpha;


	//find the initial three-point-pattern
	double s=0.01;  double counter=1;  

	while(1)
	{
		alpha = pow(2,counter-1) * s;
		g_alpha = logit_g(alpha,Add_Bias(x));
		//std::cout<<g_alpha<<std::endl;
		//std::cout<<m.w(0,0)<<std::endl;

		if(s<1e-30)
		{
			stepsize=-1; 
			g =g_alpha; 
			return stepsize;
		}
		else if(g_temp(1)==0 && g_alpha>=g_temp(0)) 
		{
			s = s/2;
			counter = counter-1;
		}
		else if(g_temp(1)==0 && g_alpha<g_temp(0)) 
		{
			g_temp(1) = g_alpha;
			b = alpha;
		}
		else if(g_temp(1)!=0 && g_temp(2)==0 && g_alpha<g_temp(1)) 
		{
			g_temp(1) = g_alpha;
			b = alpha;
		}
		else if(g_temp(1)!=0 && g_temp(2)==0 && g_alpha>g_temp(1)) 
		{	
			//std::cout<<g_temp(1)<<std::endl;
			g_temp(2) = g_alpha;
			c = alpha;
			break;
		}

		counter = counter+1;
		if(counter>100)
			break;
	}


	//if length(g)<3, stepsize=-1; g=0; return; end
	if(g_temp(1)==0 || g_temp(2)==0)
	{
		stepsize=-1; 
		g = 0; 
		return stepsize; 
	}

	double old_alpha = 0;

	//std::cout<<g_temp(0) <<" --- " << g_temp(1)<<"---"<< g_temp(2)<<std::endl;
	//std::cout<< a <<" --- " << b <<"---"<< c <<std::endl;

	while(1)
	{
		if(g_temp(0)*(c-b) + g_temp(1)*(a-c) + g_temp(2)*(b-a) == 0 || (g_temp(0)*(pow(c,2) - pow(b,2)) + g_temp(1)*(pow(a,2) - pow(c,2)) + g_temp(2)*(pow(b,2) - pow(a,2))) ==0)
		{
			break;
		}

		alpha = 0.5 *(g_temp(0)*(pow(c,2) - pow(b,2)) + g_temp(1)*(pow(a,2) - pow(c,2)) + g_temp(2)*(pow(b,2) - pow(a,2)))/(g_temp(0)*(c-b) + g_temp(1)*(a-c) + g_temp(2) * (b-a));

		//std::cout<<alpha<<std::endl;

		//% compute g(p+alpha*d)
		g_alpha = logit_g(alpha,Add_Bias(x));

		while(g_alpha == g_temp(1))
		{
			alpha = alpha*(1+1e-2);
			g_alpha = logit_g(alpha,Add_Bias(x));
		}

		if(alpha>b)
		{
			if(g_alpha < g_temp(1))
			{
				a = b;
				g_temp(0) = g_temp(1);
				b = alpha;
				g_temp(1) = g_alpha;
			}
			else if(g_alpha>g_temp(1))
			{
				c = alpha;
				g_temp(2) = g_alpha;
			}
		}
		else if(alpha<b)
		{	
			if(g_alpha < g_temp(1))
			{
				c = b;
				g_temp(2) = g_temp(1);
				b = alpha;
				g_temp(1) = g_alpha;
			}
			else if(g_alpha>g_temp(1))
			{
				a = alpha;
				g_temp(0) = g_alpha;
			}
		}

		//% termination criterion
		//	% the parabola is flat when (ga-2*gb+gc)<10*(-12) and therefore the search stops  
		if( fabs(c-a)<1e-3 || fabs(alpha - old_alpha)/fabs(alpha) < 1e-3 || (g_temp(0)- 2*g_temp(1) + g_temp(2) ) < 1e-12 )
		{
			break;
		}
		old_alpha = alpha;
		if(counter>100)
		{
			break;
		}
		counter = counter +1;

		//std::cout<<"Alpha==="<<alpha<<std::endl;

	}

	stepsize=alpha;
	g=-g_alpha;

	//	std::cout<<"========================================================="<<std::endl;
	//	std::cout<<"stepsize--->"<<stepsize<<std::endl;
	return stepsize;
}


double MCLR::logit_g(double alpha,vnl_matrix<double> data_with_bias)
{
	double gVal = 0;
	vnl_matrix<double> w_temp;
	w_temp = m.w + direction * alpha;

	vnl_matrix<double> f;

	f = Get_F_Matrix(data_with_bias,w_temp);


	vnl_vector<double> denominator(f.cols(),0);
	vnl_vector<double> f_vec(f.cols(),0);

	// Get the denominator
	for(int i=0;i<f.cols();++i)
	{
		vnl_vector<double> temp_col = f.get_column(i);
		denominator(i) = temp_col.sum();	
	}

	for(int i=0;i<f.cols();++i)
	{
		vnl_vector<double> temp_col = f.get_column(i);
		f_vec(i) = temp_col(y(i)-1)/denominator(i);

		if(f_vec(i)==0)
			f_vec(i) = 1e-9;
	}

	// Objective function value
	for(int i=0;i<f_vec.size();++i)
	{
		gVal = gVal+log(f_vec(i));
	}

	//g = g - C*sum(sum(sqrt(w.^2+delta)));   % consider the sparseness penalty 
	double diff_term =0;

	for(int i=0;i<no_of_features+1;++i)
	{
		for(int j=0;j<no_of_classes;++j)
		{
			diff_term += sqrt(w_temp(i,j)*w_temp(i,j)+delta);
		}
	}	

	diff_term = diff_term * m.sparsity_control;	
	gVal = diff_term-gVal;
	return gVal;
}



std::vector<int> MCLR::Get_Top_Features()
{
	//Need the indices of the top 5 features
	std::vector<int> indices(MIN(5,x.rows()));

	std::vector<std::pair<double, int> > val_idx; // value and index

	// find the maximum in each row
	for(int i = 1; i< m.w.rows() ; ++i)
	{
		vnl_vector<double> temp_row = m.w.get_row(i);
		val_idx.push_back(std::pair<double, int>(temp_row.max_value(),i-1));
	}

	// sorts by first element of the pair automatically
	std::sort(val_idx.begin(), val_idx.end());
	std::vector<std::pair<double, int> >::const_iterator itr;
	int counter = 0;


	reverse(val_idx.begin(),val_idx.end());

	for(itr = val_idx.begin(); itr != val_idx.begin()+MIN(5,x.rows()); ++itr)
	{
		indices[counter] = (*itr).second;			
		counter++;
	}

	return indices;
}




MCLR::model MCLR::Get_Training_Model()
{	
	// Set the z matrix 
	z.set_size(no_of_classes,x.cols());// Used in gradient computation
	z.fill(0);

	// Create the z matrix
	for(int i=0;i<x.cols();++i)
	{
		vnl_vector<double> temp_vector = z.get_column(i);
		temp_vector.put(y.get(i)-1,1);
		z.set_column(i,temp_vector);
	}


	vnl_vector<double> diff_g_3_it(3,0);//difference in g vals for last 3 iterations
	vnl_vector<double> g_4_it(4,0);	// g vals for the last three

	for(int i=0;i<1e0;++i)
	{	
		
	//	clock_t SomaExtraction_start_time = clock();
		//Get the gradient
		Get_Gradient(Add_Bias(x));
	//	std::cout << "Total time for Gradient is: " << (clock() - SomaExtraction_start_time) / (float) CLOCKS_PER_SEC << std::endl;
		
	//	clock_t SomaExtraction_start_time2 = clock();
		//Get the hessian
		Get_Hessian(Add_Bias(x));
	//	std::cout << "Total time for Hessian is: " << (clock() - SomaExtraction_start_time2) / (float) CLOCKS_PER_SEC << std::endl;

		//ameliorate hessian conditions
		Ameliorate_Hessian_Conditions();	


		//Get the direction 
		//		vnl_vector<double> dir = mclr->Newton_Direction(mclr->hessian,mclr->Column_Order_Matrix(mclr->gradient_w));
		//		mclr->direction = mclr->Reshape_Vector(dir,mclr->no_of_features+1,mclr->no_of_classes);
		direction = Reshape_Vector(Newton_Direction(hessian,Column_Order_Matrix(gradient_w)),no_of_features+1,no_of_classes);

		double step = logit_stepsize();		
		if(step == -1)
			step = 1e-9;

		m.w = m.w + step * direction;

		g_4_it(i%4) = g;	

		if(i>1)
			diff_g_3_it(i%3) = g_4_it(i%4) - g_4_it((i-1)%4);

		if (i>3 && (diff_g_3_it.one_norm()/3)/(g_4_it.one_norm()/4) < stop_cond(0))
			break;
	}



	m.FIM = hessian*-1;
	m.CRB = vnl_matrix_inverse<double>(hessian);
	//// quirk of vnl ...
	m.CRB = m.CRB*-1;



	//std::cout<<vnl_trace(m.CRB)<<std::endl;

	//FILE * fpin1 = FDeclare2("C:\\Users\\rkpadman\\Documents\\MATLAB\\crb.txt", "", 'w');
	//for(int abc =0;abc<hessian.rows();++abc)
	//{	
	//	vnl_vector<double> temp = hessian.get_row(abc);
	//	for(int abcd =0;abcd<hessian.cols();++abcd)
	//	{
	//		fprintf(fpin1, "%lf    ", temp[abcd]); 
	//	}
	//	fprintf(fpin1, "\n");
	//}

	//fclose(fpin1);
	top_features.clear();
	top_features = Get_Top_Features();
	return m;

}


MCLR::model MCLR::Get_Temp_Training_Model(int query,int label)
{	
	current_label = label;
	vnl_vector<double> queried_sample = test_data.get_row(query);

	vnl_matrix<double> temp_model_x = x;
	vnl_vector<double> temp_model_y = y;

	temp_model_x.set_size(temp_model_x.rows(),temp_model_x.cols()+1);
	temp_model_x.update(x,0,0);
	temp_model_x.set_column(temp_model_x.cols()-1,queried_sample);

	// Update y
	vnl_vector<double> temp_y = y;
	temp_model_y.set_size(temp_model_y.size()+1);
	temp_model_y.update(y,0);
	temp_model_y(y.size()) = label;

	// Set the z matrix 
	z.set_size(no_of_classes,temp_model_x.cols());// Used in gradient computation
	z.fill(0);

	// Create the z matrix
	for(int i=0;i<temp_model_x.cols();++i)
	{
		vnl_vector<double> temp_vector = z.get_column(i);
		temp_vector.put(temp_model_y.get(i)-1,1);
		z.set_column(i,temp_vector);
	}

	vnl_vector<double> diff_g_3_it(3,0);//difference in g vals for last 3 iterations
	vnl_vector<double> g_4_it(4,0);	// g vals for the last three

	for(int i =0;i<1e10;++i)
	{	
		//Get the gradient
		Get_Gradient(Add_Bias(temp_model_x));

		//Get the hessian
		Get_Hessian(Add_Bias(temp_model_x));

		//ameliorate hessian conditions
		Ameliorate_Hessian_Conditions();	


		//Get the direction 
		//		vnl_vector<double> dir = mclr->Newton_Direction(mclr->hessian,mclr->Column_Order_Matrix(mclr->gradient_w));
		//		mclr->direction = mclr->Reshape_Vector(dir,mclr->no_of_features+1,mclr->no_of_classes);
		direction = Reshape_Vector(Newton_Direction(hessian,Column_Order_Matrix(gradient_w)),no_of_features+1,no_of_classes);

		double step = logit_stepsize();		
		if(step == -1)
			step = 1e-9;

		m.w = m.w + step * direction;

		g_4_it(i%4) = g;	

		if(i>1)
			diff_g_3_it(i%3) = g_4_it(i%4) - g_4_it((i-1)%4);

		if (i>3 && (diff_g_3_it.one_norm()/3)/(g_4_it.one_norm()/4) < stop_cond(0))
			break;
	}



	m.FIM = hessian*-1;
	m.CRB = vnl_matrix_inverse<double>(hessian);
	// quirk of vnl ...
	m.CRB = m.CRB*-1;
	//std::cout<<vnl_trace(m.CRB)<<std::endl;

	//FILE * fpin1 = FDeclare2("C:\\Users\\rkpadman\\Documents\\MATLAB\\crb.txt", "", 'w');
	//for(int abc =0;abc<hessian.rows();++abc)
	//{	
	//	vnl_vector<double> temp = hessian.get_row(abc);
	//	for(int abcd =0;abcd<hessian.cols();++abcd)
	//	{
	//		fprintf(fpin1, "%lf    ", temp[abcd]);
	//	}
	//	fprintf(fpin1, "\n");
	//}

	//fclose(fpin1);
	top_features.clear();
	top_features = Get_Top_Features();
	return m;
}


FILE* MCLR::FDeclare2(char *root, char *extension, char key) {
	char string1[80];
	FILE *fpot;

	strcpy(string1, root);
	if (strcmp(extension,"")!=0) strcat(string1, ".");
	strcat(string1, extension);
	if ((key != 'w')&&(key!='r')) {
		printf("either read or write, wrong key %c \n", key);getchar();
		exit(1);
	}
	if (key == 'r') {
		//if (( (FILE *) fpot = fopen(string1, "r")) == NULL) {
		//changed for DOS?UNIX compatibility
		if ((fpot = fopen(string1, "r")) == NULL) {
			printf("\n");
			printf("input file %s does not exist\n", string1);
			exit(1); getchar();
		}
	}
	if (key == 'w') {
		if ((fpot = fopen(string1, "w")) == NULL) {
			printf("\n");
			printf("output file %s does not exist\n", string1);
			exit(1);getchar();
		}
	}
	return fpot;
}


//Updates test_data, x and y,id_time_val,test_table
void MCLR::Update_Train_Data(std::vector< std::pair<int,int> > query_label)
{	

	std::sort(query_label.begin(), query_label.end(), sort_pred());

	for(int i = 0; i<query_label.size();++i)
	{
		current_label = query_label.at(i).second;
		int query = query_label.at(i).first;
		vnl_vector<double> queried_sample = test_data.get_row(query);
		vnl_matrix<double> sub_test_matrix(test_data.rows()-1,test_data.cols());
		//std::cout<<list_of_ids.at(query) <<std::endl; 

		if(query!=0)
		{
			sub_test_matrix.update(test_data.get_n_rows(0,query),0,0);
			//for(int i =  ; i < test_data.rows() ; ++i)
			//	sub_test_matrix.set_row(i-1,test_data.get_row(i));
			sub_test_matrix.update(test_data.get_n_rows(query+1,test_data.rows()-query-1),query,0);	
		}
		else
			sub_test_matrix = test_data.get_n_rows(1,test_data.rows()-1);

		test_data = sub_test_matrix;
		//queried_sample = queried_sample.extract(queried_sample.size(),0);

		// if label == 0, then the user selected " I am not sure" option 
		if(current_label!=0)
		{
			vnl_matrix<double> temp_x = x;
			x.set_size(x.rows(),x.cols()+1);

			x.update(temp_x,0,0);
			x.set_column(x.cols()-1,queried_sample);

			// Update y
			vnl_vector<double> temp_y = y;
			y.set_size(y.size()+1);
			y.update(temp_y,0);
			y(temp_y.size()) = current_label;
		}

		//Update list of ids and the table
		id_time_val.erase(id_time_val.begin()+query);

		test_table->RemoveRow(query);
	}
}

void MCLR::Get_Label_Sample(int query)
{
	vnl_vector<double> temp_y = y;
	y.set_size(y.size()+1);
	y.update(temp_y,0);
	y(temp_y.size()) = y_ground_truth(query);

	if(query!=0)
	{
		vnl_vector<double> temp_y_gt = y_ground_truth;
		y_ground_truth.set_size(y_ground_truth.size()-1);
		y_ground_truth.update(temp_y_gt.extract(query,0),0);
		y_ground_truth.update(temp_y_gt.extract(temp_y_gt.size()-query-1,query+1),query);		
	}
	else
	{
		vnl_vector<double> temp_y_gt = y_ground_truth;
		y_ground_truth.set_size(y_ground_truth.size()-1);
		y_ground_truth.update(temp_y_gt.extract(temp_y_gt.size()-1,1),0);
	}
}


vnl_matrix<double> MCLR::Kron(vnl_vector<double> x,vnl_vector<double> y )
{
	vnl_matrix<double> q(x.size()*y.size(),1);
	int counter = 0;
	for(int j = 0; j < x.size() ; ++j)
		for(int k = 0; k < y.size() ; ++k)	
		{
			q(counter,0) = x(j)*y(k);
			counter++;
		}
		return q;
}




int MCLR::Active_Query()
{
	stop_training = false;

	// The max_info value and hence the algo converges if the user successively selects 
	//"I am not sure option" .. so whenever the user selects the "I am not sure option
	// delete that info value from max_info;
	if(current_label==0)
		max_info_vector.erase(max_info_vector.begin()+max_info_vector.size()-1);


	max_info = -1e9;
	int active_query = 0;
	//vnl_vector<double> diff_info_3_it(3,0);//difference in g vals for last 3 iterations
	//vnl_vector<double> g_3_it(3,0);	// g vals for the last three

	vnl_matrix<double> test_data_just_features =  test_data.transpose();
	//test_data_just_features = test_data_just_features.get_n_rows(0,test_data_just_features.rows()-1);
	vnl_matrix<double> prob = Test_Current_Model(test_data_just_features);


	vnl_matrix<double> test_data_bias = Add_Bias(test_data_just_features);
	vnl_vector<double> info_vector(test_data_just_features.cols());

	//Compute Information gain
	for(int i =0; i< test_data_just_features.cols();++i )
	{
		vnl_vector<double> temp_col_prob = prob.get_column(i);
		vnl_vector<double> temp_col_data = test_data_bias.get_column(i);		
		vnl_matrix<double> q = Kron(temp_col_prob,temp_col_data);//Kronecker Product
		vnl_diag_matrix<double> identity_matrix(no_of_classes,1); // Identity Matrix;
		double infoval = (q.transpose() * m.CRB.transpose() * q).get(0,0);
		vnl_matrix<double> info_matrix(no_of_classes,no_of_classes,infoval);
		info_matrix = identity_matrix + info_matrix ;
		info_vector(i) = log(vnl_determinant(info_matrix));
		if(info_vector(i)>max_info)
		{
			active_query = i;
			max_info = info_vector(i);
		}
	}

	max_info_vector.push_back(max_info);


	//	std::cout<<max_info_vector.size() << "--" << max_info <<std::endl;

	if(max_info_vector.size()>1)
		diff_info_3_it((max_info_vector.size()-1)%3) =  max_info_vector.at((max_info_vector.size()-1)) - max_info_vector.at((max_info_vector.size()-2));


	info_3_it((max_info_vector.size()-1)%3) = max_info;	

	if (max_info_vector.size()>=1) 
	{
		int sum = 0;
		for(int iter =0; iter < 3; iter++)
		{
			//std::cout<<fabs(diff_info_3_it(iter))/fabs(info_3_it(iter))<<std::endl;
			//std::cout<<stop_cond(1)<<"-->stop condition"<<std::endl;
			if(fabs(diff_info_3_it(iter))/fabs(info_3_it(iter)) < stop_cond(1))
			{
				sum = sum + 1;	
			}
		}
		// The algorithm is now confident about the problem
		// Training can be stopped
		if(sum==3)
			stop_training = true;
	}
	return active_query; 
}


std::vector<int> MCLR::ALAMO(int active_query)
{
	// Store the initial model 
	// will be replaced after performing margin sampling
	model k = m;

	std::vector<int> alamo_ids;
	alamo_ids.push_back(active_query);

	for(int i =0; i<no_of_classes; ++i)
	{

		// Update the data & refresh the training model and refresh the Training Dialog 		
		Get_Temp_Training_Model(active_query,i+1);
		vnl_matrix<double> curr_prob = Test_Current_Model(test_data.transpose());

		// Margin is used to calculate the next best samples 
		vnl_vector<double> margin;
		margin.set_size(curr_prob.cols());

		for(int j=0; j<curr_prob.cols();++j)
		{
			vnl_vector<double> temp_margin = curr_prob.get_column(j);
			double curr_class_val = temp_margin(temp_margin.arg_max());
			temp_margin(temp_margin.arg_max()) = -1;
			double next_class_val = temp_margin(temp_margin.arg_max()); // next probable class value
			margin(j) = 1/(	curr_class_val - next_class_val);
		}

		// This makes sure that the same id is not selected more than once!
		if(alamo_ids.end() == find(alamo_ids.begin(), alamo_ids.end(), margin.arg_max()))
		{
			alamo_ids.push_back(margin.arg_max());
		}
		else
		{
			margin(margin.arg_max()) = -1;
			alamo_ids.push_back(margin.arg_max());
		}

	}
	//replace the original model
	m = k;
	return alamo_ids;
}

vnl_vector<double> MCLR::Newton_Direction(vnl_matrix<double> hessian_matrix,vnl_vector<double> grad_vector )
{
	vnl_matrix<double> inv_hessian =  vnl_matrix_inverse<double>(hessian_matrix);
	vnl_vector<double> newton_direction = -inv_hessian*grad_vector;
	//std::cout<<newton_direction<<std::endl;
	return newton_direction;
}


vnl_matrix<double> MCLR::Test_Current_Model(vnl_matrix<double> test_data)
{	
	vnl_matrix<double> f;
	vnl_matrix<double> test_data_bias = Add_Bias(test_data);
	f = Get_F_Matrix(test_data_bias,m.w);	
	f = Normalize_F_Sum(f);
	return f;
}

vnl_matrix<double> MCLR::Test_Current_Model_w(vnl_matrix<double> test_data, vnl_matrix<double> m_w_matrix)
{	
	vnl_matrix<double> f;
	vnl_matrix<double> test_data_bias = Add_Bias(test_data);
	f = Get_F_Matrix(test_data_bias,m_w_matrix);	
	f = Normalize_F_Sum(f);
	return f;
}

//%% a deterministic method of finding the initial examples to acquire labels on
//%% no need to perform multiple Monte Carlo trials
//[dummy,idxL] = max(sum(x(:,1).^2,1));
//[dim,n] = size(x);
//for i=1:25
//    idxU = setdiff(1:n,idxL);
//    z = x(:,idxL);
//    P = eye(dim) - z*inv(z'*z)*z';
//    [dummy,idxL(end+1)] = max(sum((P*x(:,idxU)).^2,1));
//end

std::vector<std::pair<int,int> > MCLR::Plan_In_Advance(vtkSmartPointer<vtkTable> new_table, int num,std::vector< std::pair<int,int> > id_time_PIA)
{	

	vtkSmartPointer<vtkTable> table_PIA  = vtkSmartPointer<vtkTable>::New();
	table_PIA->Initialize();

	for(int col=0; col<new_table->GetNumberOfColumns(); ++col)
	{	
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName(new_table->GetColumnName(col));
		table_PIA->AddColumn(column);	
	}

	for(int row = 0; row < (int)new_table->GetNumberOfRows(); ++row)
	{	
		vtkSmartPointer<vtkVariantArray> model_data1 = vtkSmartPointer<vtkVariantArray>::New();
		for(int c =0;c<(int)table_PIA->GetNumberOfColumns();++c)
			model_data1->InsertNextValue(new_table->GetValueByName(row,table_PIA->GetColumnName(c)));
		table_PIA->InsertNextRow(model_data1);
	}

	table_PIA->RemoveColumnByName("ID");

	// class_list has all -1s
	// both class_list and id_time are dummies in this case
	// Since we are using the same functions as Active Learning

	vnl_vector<double> class_list(table_PIA->GetNumberOfRows(),-1); 			
	std::vector<std::pair<int,int> > query_label;
	query_label.resize(1);

	// No training data
	// we want to deal with all the examples belonging to the class~classval(1,2,3.....) and pick by "Planning in Advance"
	// Normalize the feature matrix

	vnl_matrix<double> Feats = Normalize_Feature_Matrix(tableToMatrix(table_PIA,id_time_PIA));
	Initialize(Feats,1,class_list,"",table_PIA,false);
	vnl_matrix<double> tempUnlabeledData = test_data.transpose();


	//zmat rows-> features and columns-> samples
	vnl_matrix<double> zMat; // refer Plan in Adavance Matlab code for interpretation
	zMat.set_size(tempUnlabeledData.rows(),num);
	//zMat.set(0);

	vnl_matrix<double> PMat; // refer Plan in Adavance Matlab code for interpretation
	PMat.set_size(tempUnlabeledData.rows(),tempUnlabeledData.rows());// test_data.cols() = dim !	
	vnl_diag_matrix<double> identity_matrix(test_data.cols(),1); // Identity Matrix;

	std::vector<std::pair<int,int> > PIA;
	PIA.resize(num);

	query_label[0].first = 0;
	query_label[0].second = 0;//dummy again
	zMat.set_column(0,tempUnlabeledData.get_column(0)); // first sample
	PIA[0] = id_time_PIA[0];


	for(int i =0;i<num-1;++i)
	{
		vnl_matrix<double> tempZMat = zMat.extract(tempUnlabeledData.rows(),i+1,0,0);
		// Update tempUnlabeledData		
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		vnl_matrix<double> sub_test_matrix(tempUnlabeledData.rows(),tempUnlabeledData.cols()-1);

		if(query_label[0].first !=0)
		{
			sub_test_matrix.update(tempUnlabeledData.get_n_columns(0,query_label[0].first),0,0);
			sub_test_matrix.update(tempUnlabeledData.get_n_columns(query_label[0].first+1,tempUnlabeledData.cols()-query_label[0].first-1),0,query_label[0].first);	
		}
		else
			sub_test_matrix = tempUnlabeledData.get_n_columns(1,tempUnlabeledData.columns()-1);

		tempUnlabeledData = sub_test_matrix;
		//id_list.erase(id_list.begin()+query_label[0].first);
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Always pick the first sample and then pick the rest
		Update_Train_Data(query_label);

		clock_t RemoveIntraSomaNodes_start_time = clock();
		vnl_matrix<double> invTempZMat = vnl_matrix_inverse<double>(tempZMat.transpose()*tempZMat);		
		PMat = identity_matrix-(tempZMat*invTempZMat*tempZMat.transpose());
		vnl_matrix<double> PxIdu = PMat*tempUnlabeledData;
		vnl_vector<double> sum_vec(tempUnlabeledData.cols(),0);

		for(int j = 0;j<PxIdu.cols(); ++j)
		{
			for(int k = 0;k< PxIdu.rows(); ++k)
			{
				sum_vec(j) = sum_vec(j)+ (PxIdu(k,j)*PxIdu(k,j));
			}
		}
		query_label[0].first = sum_vec.arg_max();
		query_label[0].second = 0;// dummy value 
		PIA[i+1] = id_time_val.at(sum_vec.arg_max()); // Will contain the ids
		zMat.set_column(i+1,tempUnlabeledData.get_column(query_label[0].first)); // first sample

	}
	return PIA;
}



int MCLR::GetNumberOfClasses(vtkSmartPointer<vtkTable> table)
{
	std::vector<std::string> prediction_names;
	//Get prediction colums, if any, for center map coloring
	prediction_names.clear();
	prediction_names = ftk::GetColumsWithString( "prediction" , table );

	int no_of_classes = -1;

	for(int row = 0; row < (int)table->GetNumberOfRows(); ++row)
	{
		if(table->GetValueByName(row,prediction_names[0].c_str())>no_of_classes)
			no_of_classes = table->GetValueByName(row,prediction_names[0].c_str()).ToInt();
	}

	return no_of_classes;
}



std::vector< std::pair< std::string, vnl_vector<double> > > MCLR::CreateActiveLearningModel(vtkSmartPointer<vtkTable> pWizard_table)
{
	act_learn_model.clear();

	vnl_vector<double> std_dev_vec;
	vnl_vector<double> mean_vec;
	act_learn_matrix = GetActiveLearningMatrix();
	std_dev_vec = Get_Std_Dev();
	mean_vec = Get_Mean();
	std::string feat_name = "bias";
	vnl_vector<double> feat_values(act_learn_matrix.get_row(0).size() + 2);
	feat_values.put(0, 0);
	feat_values.put(1, 0);
	for(int s=0; s<(int)act_learn_matrix.get_row(0).size(); ++s)
	{
		feat_values.put(s+2, act_learn_matrix.get_row(0).get(s));
	}
	std::pair< std::string, vnl_vector<double> > temp = std::make_pair(feat_name, feat_values);
	act_learn_model.push_back(temp);
	for(unsigned int col=0; col<(unsigned int)pWizard_table->GetNumberOfColumns(); ++col)
	{
		feat_name = pWizard_table->GetColumnName(col);
		feat_values.put(0, std_dev_vec.get(col));
		feat_values.put(1, mean_vec.get(col));
		for(int s=0; s<(int)act_learn_matrix.get_row(col).size(); ++s)
		{
			feat_values.put(s+2, act_learn_matrix.get_row(col+1).get(s));
		}
		std::pair< std::string, vnl_vector<double> > temp = std::make_pair(feat_name, feat_values);
		act_learn_model.push_back(temp);
	}
	return act_learn_model;
}


double MCLR::PerformPTest(vtkSmartPointer<vtkTable> featureTable,std::vector<int> ground_truth,double errorVal)
{

	vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
	column->SetName("prediction_random");
	column->SetNumberOfValues( featureTable->GetNumberOfRows() );
	featureTable->AddColumn(column);
	
	std::vector<double> error;
	double pVal =0;

	
	for(int i=0; i<30; ++i) // 25 random interations
	{	
//		clock_t SomaExtraction_start_time = clock();
		// Error for one iteration of random LOOCV
		error.push_back(LOOCV(featureTable,ground_truth));
		std::cout<<"Running random iteration " << i <<std::endl;
		if(error[i]<=errorVal)
			pVal = pVal +1;
		//std::cout << "Total time for training is: " << (clock() - SomaExtraction_start_time) / (float) CLOCKS_PER_SEC << std::endl;
		//std::cout<<"---------------------------------------------------------------------------------------------------"<<std::endl;
	}
	
	return pVal/31;
}


//Leave One Out Cross Validation
double MCLR::LOOCV(vtkSmartPointer<vtkTable> featureTable,std::vector<int> ground_truth)
{

	std::vector<int> randomclass = ground_truth;
	random_shuffle(randomclass.begin(), randomclass.end());
	double error = 0;
	vnl_vector<double> class_list(featureTable->GetNumberOfRows()); 

//#pragma omp parallel for private(test_data)
	for(int i = 0; i < (int)featureTable->GetNumberOfRows(); ++i)
	{	
		vtkSmartPointer<vtkTable> table_validation = ftk::CopyTable(featureTable);

		for(int j=0; j<table_validation->GetNumberOfRows(); ++j)
		{
			if(j==i)
			{
				table_validation->SetValueByName(j,"prediction_random",-1);
				class_list.put(j,vtkVariant(table_validation->GetValueByName(j,"prediction_random")).ToDouble());
					
			}
			else
			{
				table_validation->SetValueByName(j,"prediction_random",randomclass[j]);
				class_list.put(j,vtkVariant(table_validation->GetValueByName(j,"prediction_random")).ToDouble());
			}
		}
			

		table_validation->RemoveColumnByName("prediction_random");

		////// Final Data  to classify after the active training
		vnl_matrix<double> data_classify;
		std::vector< std::pair<int,int> > dummy_id_time;

		double sparsity = 1;
		double max_info = -1e9;
		
		// Normalize the feature matrix
		vnl_matrix<double> Feats = Normalize_Feature_Matrix(tableToMatrix(table_validation,dummy_id_time));
		Initialize(Feats,sparsity,class_list,"",table_validation,true);
		
		Get_Training_Model();

		data_classify = test_data.transpose();

		vnl_matrix<double> currprob;
		currprob = Test_Current_Model(data_classify);

	//	std::cout<<currprob.get_column(0).arg_max()+1  <<std::endl;

		if(currprob.get_column(0).arg_max()+1 != ground_truth[i] )
			error = error+1;
	}
	
	return (error/featureTable->GetNumberOfRows());
}

