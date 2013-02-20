#include "mclr.h"

#ifdef _OPENMP
#include "omp.h"
#endif

MCLR::MCLR()
{
	current_label = -1; // keeps a track of the current label
	confidence_threshold = 0.5;
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
void MCLR::Initialize(vnl_matrix<double> data,double c,vnl_vector<double> classes, std::string str,vtkSmartPointer<vtkTable> table)
{
	int train_counter =0;
	int test_counter = 0;
	
	//bool PIA =false;

	// test_table contains the table after the queried samples are removed from the original table
	// It is updated in the Update_trainData 
	test_table = table; //Features->Columns ; Samples -> Rows

	for(int i = 0; i< classes.size() ; ++i)
	{
		if(classes(i) == -1)
			test_counter++;
	}

	trainData.set_size(data.rows()-test_counter,data.cols());
	testData.set_size(test_counter,data.cols());

	test_counter = 0;
	//variable to keep track of the number of rows deleted
	int num_del_rows = 0;
	for(int i= 0;i<data.rows();++i)
	{
		if(classes(i)!=-1)
		{	
			trainData.set_row(train_counter,data.get_row(i));
			++train_counter;
			test_table->RemoveRow(i - num_del_rows);
			++num_del_rows;
		}
		else
		{
			testData.set_row(test_counter,data.get_row(i));
			++test_counter;
		}
	}

	y.set_size(trainData.rows());

	// Transpose : Features->Rows ; Samples -> Columns
	x = trainData.transpose();// Feature matrix of training samples only ; does not contain unlabeled sample data

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
	stop_cond[1] = 0.2;
	delta = 1e-9; // used for approximating |w| 

	numberOfFeatures = x.rows();
	// We assume that the classes are numbered from 1 to "n" for n classes
	// so y.maxval - .minval +1 gives number of classes

		numberOfClasses = y.max_value()-y.min_value()+1;
		class_vector.set_size(numberOfClasses);

		m.w.set_size(numberOfFeatures+1,numberOfClasses);
		m.w.fill(0);

		gradient_w.set_size(numberOfFeatures+1,numberOfClasses);
		hessian.set_size((numberOfFeatures+1)*(numberOfClasses),(numberOfFeatures+1)*(numberOfClasses)); /// Check the size again!!!

		// g :the objective function value
		g = 0;
		// Class vector
		for(int i=1;i<=numberOfClasses;++i)
		{
			class_vector.put(i-1,i); 
		}

		diff_info_3_it.set_size(3);
		info_3_it.set_size(3);
	
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

// rearrange feature matrix based on current classification and current top features
// this will reflect in the heat map at each AL iteration
vtkSmartPointer<vtkTable> MCLR::Rearrange_Table(vtkSmartPointer<vtkTable> pawTable)
{

	vnl_matrix<double> feats = tableToMatrix(pawTable,this->id_time_val);

	vnl_matrix<double> feats2 = Normalize_Feature_Matrix(tableToMatrix(pawTable, this->id_time_val));
	feats2 = feats2.transpose();

	vnl_matrix<double> rearrangeFeatsTemp(feats.rows(),feats.cols());
	vnl_matrix<double> currprob = Test_Current_Model_w(feats2, m.w);
		

	int counter = 0;
	// Re arrange the rows of the matrix according to class
	for(int i = 0; i<numberOfClasses ; ++i)
	{
		for(int j = 0; j<feats.rows() ; ++j)
		{
			vnl_vector<double> curr_col = currprob.get_column(j);

			if(curr_col.arg_max() == i && curr_col.max_value() > confidence_threshold) 
			{	
				rearrangeFeatsTemp.set_row(counter,feats.get_row(j));			
				counter++;
			}
		}
	}

	//// Class 0 examples at the end 
	//for(int i = 0; i<numberOfClasses ; ++i)
	//{
	//	for(int j = 0; j<feats.rows() ; ++j)
	//	{
	//		vnl_vector<double> curr_col = currprob.get_column(j);

	//		if(curr_col.max_value() <= confidence_threshold) 
	//		{	
	//			rearrangeFeatsTemp.set_row(counter,feats.get_row(j));			
	//			counter++;
	//		}
	//	}
	//}

	
	std::vector<int> priorityOrder = Get_Feature_Order();
	counter = 0;	
	
	vnl_matrix<double> rearrangeFeats(feats.rows(),priorityOrder.size());

	for(int j = 0; j<rearrangeFeats.cols() ; ++j)
	{
		rearrangeFeats.set_column(j,rearrangeFeatsTemp.get_column(priorityOrder[j]));
		//std::cout<< priorityOrder[j] << std::endl;
	}

	//std::cout<< " ESCAPE !!!!!!!!!! " <<std::endl;

	vtkSmartPointer<vtkTable> tableRearrange = vtkSmartPointer<vtkTable>::New();
	tableRearrange->Initialize();

	vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
	column->SetName("Dummy");
	column->SetNumberOfValues(rearrangeFeats.rows());
	tableRearrange->AddColumn(column);
	for(int row = 0; (int)row < tableRearrange->GetNumberOfRows(); ++row)  
	{
		tableRearrange->SetValue(row, 0, vtkVariant(row));
	}


	for(int i = 0; i<rearrangeFeats.cols(); ++i)
	{
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName(pawTable->GetColumnName(priorityOrder[i]));
		column->SetNumberOfValues(rearrangeFeats.rows());
		tableRearrange->AddColumn(column);
	}

	for(int row = 0; (int)row < tableRearrange->GetNumberOfRows(); ++row)  
	{
		for(int col = 1; (int)col < tableRearrange->GetNumberOfColumns(); ++col)  
		{
			tableRearrange->SetValue(row, col, vtkVariant(rearrangeFeats(row,col-1)));
		}
	}
	
	return tableRearrange;
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
		for(int j=0;j<numberOfClasses;++j)
		{
			sum = sum + f(j,i);
		}
		norm_matrix_row(i) = sum;
	}


	for(int i=0;i<numberOfClasses;++i)
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


void MCLR::Get_Gradient(vnl_matrix<double> data_with_bias)
{
	// f  = exp(w'*[ones(1,N);basis]); In MATLAB 
	vnl_matrix<double> f;
	f = Get_F_Matrix(data_with_bias,m.w);
	f = Normalize_F_Sum(f);

	//   grad_w = grad_w - C*(w./sqrt(w.^2+delta));
	// |w| approximation

	vnl_matrix<double> abs_w; //
	abs_w.set_size(numberOfFeatures+1,numberOfClasses);

	//#pragma omp parallel for num_threads(16) 
	for(int i=0;i<numberOfFeatures+1;++i)
	{
		for(int j=0;j<numberOfClasses;++j)
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

	data_with_bias.set_size(numberOfFeatures+1,data.cols());
	data_with_bias.set_row(0,temp_vector);
	for(int i =0; i<numberOfFeatures ; ++i)
	{
		data_with_bias.set_row(i+1,data.get_row(i));	
	}

	return data_with_bias;
}


vnl_matrix<double> MCLR::Get_Hessian(vnl_matrix<double> data_with_bias,vnl_matrix<double> w)
{	
	vnl_matrix<double> f;
	f = Get_F_Matrix(data_with_bias,w);
	f = Normalize_F_Sum(f);

	vnl_matrix<double> temp_hessian; // temporary Hessian Matrix
	temp_hessian.set_size((numberOfFeatures+1)*(numberOfClasses),(numberOfFeatures+1)*(numberOfClasses));

	for(int i = 1; i<=numberOfClasses ;++i)
	{
		vnl_vector<double> ith_row  = f.get_row(i-1);
		
		for(int j = 0; j< x.cols() ; ++j)
		{
			ith_row(j) = ith_row(j)*(1-ith_row(j));
		}


		vnl_diag_matrix<double> diagonal_matrix(ith_row);
		vnl_matrix<double> temp_matrix =  -data_with_bias*diagonal_matrix*data_with_bias.transpose();
		temp_hessian.update(temp_matrix,(i-1)*(numberOfFeatures+1),(i-1)*(numberOfFeatures+1));


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
			temp_hessian.update(temp_matrix,(i-1)*(numberOfFeatures+1),(all_class_but_i(k)-1)*(numberOfFeatures+1));

		}
	}
	//hessian = hessian - diag(C*(delta./(sqrt(w(:).^2+delta)).^3));	
	//Need to vectorize w 

	vnl_vector<double> w_column_ordered = Column_Order_Matrix(w);
	int count = 0;
	for(int i=0;i<w.rows();++i)
	{
		for(int j=0;j<w.cols();++j)	
		{
			w_column_ordered(j+i*w.cols()) = (m.sparsity_control)*(delta/pow(sqrt(pow(w_column_ordered(j+i*w.cols()),2)+delta),3));
			//count++;
		}
	}

	vnl_diag_matrix<double> diagonal_matrix(w_column_ordered);
	temp_hessian = temp_hessian - diagonal_matrix;
	
	return temp_hessian;
}	


void MCLR::Ameliorate_Hessian_Conditions()
{	

	//std::cout<< "Into Amel Hessian" << std::endl;
	//Compute the eigen vectors of the hessian matrix
	vnl_symmetric_eigensystem<double> eig(hessian);

	double ratio = fabs(eig.get_eigenvalue(0))/fabs(eig.get_eigenvalue((numberOfFeatures+1)*numberOfClasses)-1);
	int counter = 0;

	if(ratio>1e9)
	{
		vnl_diag_matrix<double> diag((numberOfClasses)*(numberOfFeatures+1),Compute_Mean_Abs_Eig(eig));
		hessian = hessian + diag;
		//vnl_symmetric_eigensystem<double> eig(hessian);
		//ratio = fabs(eig.get_eigenvalue(0))/fabs(eig.get_eigenvalue((numberOfFeatures+1)*numberOfClasses)-1);
	}
}


double MCLR::Compute_Mean_Abs_Eig(vnl_symmetric_eigensystem<double> eig)
{
	double sum = 0; 	
	for(int i =0; i<(numberOfFeatures+1)*(numberOfClasses);++i)
	{
		sum = sum + fabs(eig.get_eigenvalue(i));
	}	
	sum = sum/((numberOfFeatures+1)*(numberOfClasses)) ;
	return sum;
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


	}

	stepsize=alpha;
	g=-g_alpha;

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

	for(int i=0;i<numberOfFeatures+1;++i)
	{
		for(int j=0;j<numberOfClasses;++j)
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


std::vector<int> MCLR::Get_Feature_Order()
{	

	std::vector<std::pair<double, int> > val_idx; // value and index
	vnl_vector<double> normW(m.w.rows()-1);
	

	// find the maximum in each row and get them into a vector
	// This vector will be normalized in the next step.
	// Features with weight > 0.1 will be considered for the reduced heatmap
	// Remove Bias rows. hence i = 1
	
	double maxW = -1;

	for(int i = 1; i< m.w.rows() ; ++i)
	{
		vnl_vector<double> temp_row = m.w.get_row(i);
		normW(i-1) = fabs(temp_row.max_value());
		if(normW(i-1)>maxW)
			maxW = normW(i-1);
	}

	for(int i = 0; i< normW.size() ; ++i)
	{
		normW(i) = normW(i)/maxW;
	}

	//getchar();

	for(int i = 0; i< normW.size() ; ++i)
	{
		if(normW(i)>0)
		{
			val_idx.push_back(std::pair<double, int>(normW(i),i));
		}
	}

	// sorts by first element of the pair automatically
	std::sort(val_idx.begin(), val_idx.end());
	std::vector<std::pair<double, int> >::const_iterator itr;
	int counter = 0;

	reverse(val_idx.begin(),val_idx.end());

	std::vector<int> indices;

	//std::cout<<val_idx.size()<<std::endl;
	

	for(itr = val_idx.begin(); itr != val_idx.end(); ++itr)
	{	
		//std::cout<< (*itr).second <<std::endl;
		indices.push_back((*itr).second);
		counter++;
	}
	//std::cout<< " ----------------" <<std::endl;
	return indices;
}

MCLR::model MCLR::Get_Training_Model()
{	
	// Set the z matrix 
	z.set_size(numberOfClasses,x.cols());// Used in gradient computation
	z.fill(0);

	// Create the z matrix
	for(int i=0;i<x.cols();++i)
	{
		vnl_vector<double> temp_vector = z.get_column(i);
		temp_vector.put(y.get(i)-1,1);
		z.set_column(i,temp_vector);
	}


	vnl_vector<double> diff_g_3_it(3,0);//difference in g vals for last 3 iterations
	vnl_vector<double> g_3_it(3,0);	// g vals for the last three

	for(int i=0;i<1e10;++i)
	{	
		
		//Get the gradient
		Get_Gradient(Add_Bias(x));
		
		//Get the hessian
		hessian = Get_Hessian(Add_Bias(x),m.w);

		//ameliorate hessian conditions
		Ameliorate_Hessian_Conditions();	

		//Get the direction 
		//		vnl_vector<double> dir = mclr->Newton_Direction(mclr->hessian,mclr->Column_Order_Matrix(mclr->gradient_w));
		//		mclr->direction = mclr->Reshape_Vector(dir,mclr->numberOfFeatures+1,mclr->numberOfClasses);
		direction = Reshape_Vector(Newton_Direction(hessian,Column_Order_Matrix(gradient_w)),numberOfFeatures+1,numberOfClasses);

		double step = logit_stepsize();		
		if(step == -1)
			step = 1e-9;

		m.w = m.w + step * direction;

		g_3_it(i%3) = g;	

		if(i>1)
			diff_g_3_it(i%3) = g_3_it(i%3) - g_3_it((i-1)%3);

		if (i>3 && (diff_g_3_it.one_norm()/3)/(g_3_it.one_norm()/3) < stop_cond(0))
			break;
	}

	m.FIM = hessian*-1;
	m.CRB = vnl_matrix_inverse<double>(hessian);
	//// quirk of vnl ...
	m.CRB = m.CRB*-1;
	
	// Top Features selected based on sparsity of the model
	top_features.clear();
	top_features = Get_Top_Features();
	return m;
}





//Updates testData, x and y,id_time_val,test_table
void MCLR::Update_Train_Data(std::vector< std::pair<int,int> > query_label)
{	

	std::sort(query_label.begin(), query_label.end(), sort_pred());

	for(int i = 0; i<query_label.size();++i)
	{
		current_label = query_label.at(i).second;
		int query = query_label.at(i).first;
		vnl_vector<double> queried_sample = testData.get_row(query);
		vnl_matrix<double> sub_test_matrix(testData.rows()-1,testData.cols());
		//std::cout<<list_of_ids.at(query) <<std::endl; 

		if(query!=0)
		{
			sub_test_matrix.update(testData.get_n_rows(0,query),0,0);
			//for(int i =  ; i < testData.rows() ; ++i)
			//	sub_test_matrix.set_row(i-1,testData.get_row(i));
			sub_test_matrix.update(testData.get_n_rows(query+1,testData.rows()-query-1),query,0);	
		}
		else
			sub_test_matrix = testData.get_n_rows(1,testData.rows()-1);

		testData = sub_test_matrix;
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


int MCLR::Active_Query()
{
	stop_training = false;

	// The max_info value and hence the algo converges if the user successively selects 
	//"I am not sure option" .. so whenever the user selects the "I am not sure option
	// delete that info value from max_info;
	if(current_label==0)
		maxInfoVector.erase(maxInfoVector.begin()+maxInfoVector.size()-1);


	max_info = -1e9;
	int activeQuery = 0;
	//vnl_vector<double> diff_info_3_it(3,0);//difference in g vals for last 3 iterations
	//vnl_vector<double> g_3_it(3,0);	// g vals for the last three

	vnl_matrix<double> testDataTranspose =  testData.transpose();
	//testDataTranspose = testDataTranspose.get_n_rows(0,testDataTranspose.rows()-1);
	vnl_matrix<double> prob = Test_Current_Model(testDataTranspose);
	vnl_matrix<double> testDataBias = Add_Bias(testDataTranspose);
	vnl_vector<double> infoVector(testDataTranspose.cols());

	//Compute Information gain
	for(int i =0; i< testDataTranspose.cols();++i )
	{
		vnl_vector<double> tempColProb = prob.get_column(i);
		vnl_diag_matrix<double> diagProb(tempColProb);
		vnl_matrix<double> tempProb(tempColProb.size(),tempColProb.size());

		for(int p =0; p< tempColProb.size();++p)
			for(int q =0; q< tempColProb.size();++q)
				tempProb(p,q) = tempColProb(p)*tempColProb(q);

		vnl_symmetric_eigensystem<double> eig(diagProb-tempProb);
		vnl_matrix<double> DMatrix = eig.D.asMatrix();

		// Rounding really small negative numbers 
		// in D Matrix to 0 since we apply sqrt
		for(int p =0; p< tempColProb.size();++p)
		{
			for(int q =0; q< tempColProb.size();++q)
			{
				if(fabs(DMatrix(p,q))<1e-3)
					DMatrix(p,q) = 0;
			}
		}
		vnl_matrix<double> Dsqrt = DMatrix.apply(sqrt); // Square root all elements
		vnl_vector<double> tempColData = testDataBias.get_column(i);
		vnl_matrix<double> q = Kron(eig.V*Dsqrt,tempColData);//Kronecker Product
		vnl_diag_matrix<double> identity_matrix(numberOfClasses,1); // Identity Matrix;
		vnl_matrix<double> infoMatrix = identity_matrix + q.transpose() * m.CRB.transpose() * q ;
		infoVector(i) = log(vnl_determinant(infoMatrix));
		if(infoVector(i)>max_info)
		{
			activeQuery = i;
			max_info = infoVector(i);
		}
	}

	maxInfoVector.push_back(max_info);
	std::cout<<"Max Info-"  <<max_info << std::endl;

	if(maxInfoVector.size()>1)
		diff_info_3_it((maxInfoVector.size()-1)%3) =  maxInfoVector.at((maxInfoVector.size()-1)) - maxInfoVector.at((maxInfoVector.size()-2));


	info_3_it((maxInfoVector.size()-1)%3) = max_info;	

	// Minimum 3 iterations - > maxInfoVector.size()=9 
	//( assuming batch size =3)
	if (maxInfoVector.size()>=9) 
	{
		int sum = 0;
		for(int iter =0; iter < 3; iter++)
		{
			//std::cout<<fabs(diff_info_3_it(iter))/fabs(info_3_it(iter)) << "--stopcond" <<std::endl;	
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
	return activeQuery; 
}


//-----------------------------------------------------------------------------------------------------------------------------
// Batch Mode Methods : Included two methods to do batch mode : Submodularity based & Margin Sampling based
//-----------------------------------------------------------------------------------------------------------------------------

std::vector<int> MCLR::ALAMO(int activeQuery)
{
	// Store the initial model 
	// will be replaced after performing margin sampling
	model k = m;

	std::vector<int> alamo_ids;
	alamo_ids.push_back(activeQuery);

	for(int i =0; i<numberOfClasses; ++i)
	{

		// Update the data & refresh the training model and refresh the Training Dialog 		
		Get_Temp_Training_Model(activeQuery,i+1);
		vnl_matrix<double> curr_prob = Test_Current_Model(testData.transpose());

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


MCLR::model MCLR::Get_Temp_Training_Model(int query,int label)
{	
	current_label = label;
	vnl_vector<double> queried_sample = testData.get_row(query);

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
	z.set_size(numberOfClasses,temp_model_x.cols());// Used in gradient computation
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
		Get_Hessian(Add_Bias(temp_model_x),m.w);

		//ameliorate hessian conditions
		Ameliorate_Hessian_Conditions();	


		//Get the direction 
		//		vnl_vector<double> dir = mclr->Newton_Direction(mclr->hessian,mclr->Column_Order_Matrix(mclr->gradient_w));
		//		mclr->direction = mclr->Reshape_Vector(dir,mclr->numberOfFeatures+1,mclr->numberOfClasses);
		direction = Reshape_Vector(Newton_Direction(hessian,Column_Order_Matrix(gradient_w)),numberOfFeatures+1,numberOfClasses);

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
	top_features.clear();
	top_features = Get_Top_Features();
	return m;
}



std::vector<int> MCLR::Submodular_AL(int activeQuery,vnl_matrix<double> testDataTemp )
{
	
	std::vector<int> submodularALQueries;
	submodularALQueries.push_back(activeQuery);
	
	int batchSize = 0; 
	
	vnl_matrix<double> currentCRB;
	vnl_matrix<double> tempX = x; // x corresponds to the feature matrix of current labeled samples

	std::vector<int> rowIDs;

	for(int i=0;i<testDataTemp.rows();++i)
		rowIDs.push_back(i);
	
	for (int batch =0; batch<batchSize; ++batch)
	{	
		max_info = -1e9;
		tempX = Update_Temp_x(activeQuery,testDataTemp,tempX);
		testDataTemp = Update_Temp_Train_Data(activeQuery,testDataTemp);
		rowIDs.erase(rowIDs.begin()+activeQuery);
		
		vnl_matrix<double> testDataTranspose =  testDataTemp.transpose();
		//testDataTranspose = testDataTranspose.get_n_rows(0,testDataTranspose.rows()-1);
		vnl_matrix<double> prob = Test_Current_Model(testDataTranspose);
		vnl_matrix<double> testDataBias = Add_Bias(testDataTranspose);
		vnl_vector<double> infoVector(testDataTranspose.cols());
		currentCRB = Get_Hessian(Add_Bias(tempX),m.w);
		currentCRB = vnl_matrix_inverse<double>(currentCRB);
		//// quirk of vnl ...
		currentCRB = currentCRB*-1;

		//Compute Information gain	
		for(int i =0; i< testDataTranspose.cols();++i )
		{
			vnl_vector<double> tempColProb = prob.get_column(i);
			vnl_diag_matrix<double> diagProb(tempColProb);
			vnl_matrix<double> tempProb(tempColProb.size(),tempColProb.size());

			for(int p =0; p< tempColProb.size();++p)
				for(int q =0; q< tempColProb.size();++q)
					tempProb(p,q) = tempColProb(p)*tempColProb(q);

			vnl_symmetric_eigensystem<double> eig(diagProb-tempProb);
			vnl_matrix<double> DMatrix = eig.D.asMatrix();

			// Rounding really small negative numbers 
			// in D Matrix to 0 since we apply sqrt
			for(int p =0; p< tempColProb.size();++p)
			{
				for(int q =0; q< tempColProb.size();++q)
				{
					if(fabs(DMatrix(p,q))<1e-3)
						DMatrix(p,q) = 0;
				}
			}
			vnl_matrix<double> Dsqrt = DMatrix.apply(sqrt); // Square root all elements
			vnl_vector<double> tempColData = testDataBias.get_column(i);
			vnl_matrix<double> q = Kron(eig.V*Dsqrt,tempColData);//Kronecker Product
			vnl_diag_matrix<double> identity_matrix(numberOfClasses,1); // Identity Matrix;
			vnl_matrix<double> infoMatrix = identity_matrix + q.transpose() * currentCRB * q ;
			infoVector(i) = log(vnl_determinant(infoMatrix));
			
			if(infoVector(i)>max_info)
			{
				activeQuery = i;
				max_info = infoVector(i);
			}
		

		// Check for correct row ----------------------------------------
		}
		maxInfoVector.push_back(max_info);
		if(maxInfoVector.size()>1)
			diff_info_3_it((maxInfoVector.size()-1)%3) =  maxInfoVector.at((maxInfoVector.size()-1)) - maxInfoVector.at((maxInfoVector.size()-2));
	
		info_3_it((maxInfoVector.size()-1)%3) = max_info;	

		std::cout<<"Max Info-"  <<max_info << std::endl;
		submodularALQueries.push_back(rowIDs.at(activeQuery));
	}
	
	return submodularALQueries;
}



vnl_matrix<double> MCLR::Update_Temp_Train_Data(int query,vnl_matrix<double> testDataTemp )
{	

	vnl_vector<double> queried_sample = testDataTemp.get_row(query);
	vnl_matrix<double> sub_test_matrix(testDataTemp.rows()-1,testDataTemp.cols());

		if(query!=0)
		{
			sub_test_matrix.update(testDataTemp.get_n_rows(0,query),0,0);
			//for(int i =  ; i < testData.rows() ; ++i)
			//	sub_test_matrix.set_row(i-1,testData.get_row(i));
			sub_test_matrix.update(testDataTemp.get_n_rows(query+1,testDataTemp.rows()-query-1),query,0);	
		}
		else
			sub_test_matrix = testDataTemp.get_n_rows(1,testDataTemp.rows()-1);



	return sub_test_matrix;
}



vnl_matrix<double> MCLR::Update_Temp_x(int query, vnl_matrix<double> tempTestData, vnl_matrix<double> xBatchMode )
{		
		vnl_matrix<double> tempX = xBatchMode;
		tempX.set_size(xBatchMode.rows(),xBatchMode.cols()+1);
		tempX.update(xBatchMode,0,0);
		vnl_vector<double> queriedSample = tempTestData.get_row(query);
		tempX.set_column(tempX.cols()-1,queriedSample);
		return tempX;
}



//-----------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------


vnl_vector<double> MCLR::Newton_Direction(vnl_matrix<double> hessian_matrix,vnl_vector<double> grad_vector )
{
	vnl_matrix<double> inv_hessian =  vnl_matrix_inverse<double>(hessian_matrix);
	vnl_vector<double> newton_direction = -inv_hessian*grad_vector;
	//std::cout<<newton_direction<<std::endl;
	return newton_direction;
}


vnl_matrix<double> MCLR::Test_Current_Model(vnl_matrix<double> testDataLocal)
{	
	vnl_matrix<double> f;
	vnl_matrix<double> testDataBias = Add_Bias(testDataLocal);
	f = Get_F_Matrix(testDataBias,m.w);	
	f = Normalize_F_Sum(f);
	return f;
}



//-----------------------------------------------------------------------------------------------------------------------------
// Overloaded Methods
//-----------------------------------------------------------------------------------------------------------------------------


vnl_matrix<double> MCLR::Test_Current_Model_w(vnl_matrix<double> testDataLocal, vnl_matrix<double> m_w_matrix)
{	
	vnl_matrix<double> f;
	vnl_matrix<double> testDataBias = Add_Bias(testDataLocal);
	f = Get_F_Matrix(testDataBias,m_w_matrix);	
	f = Normalize_F_Sum(f);
	return f;
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

//-----------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------



int MCLR::GetNumberOfClasses(vtkSmartPointer<vtkTable> table)
{
	std::vector<std::string> prediction_names;
	//Get prediction colums, if any, for center map coloring
	prediction_names.clear();
	prediction_names = ftk::GetColumsWithString( "prediction" , table );

	int numberOfClasses = -1;

	for(int row = 0; row < (int)table->GetNumberOfRows(); ++row)
	{
		if(table->GetValueByName(row,prediction_names[0].c_str())>numberOfClasses)
			numberOfClasses = table->GetValueByName(row,prediction_names[0].c_str()).ToInt();
	}

	return numberOfClasses;
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


//-----------------------------------------------------------------------------------------------------------------------------
// C++ Methods for MATLAB functions
//-----------------------------------------------------------------------------------------------------------------------------

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


vnl_matrix<double> MCLR::Kron(vnl_matrix<double> x,vnl_vector<double> y )
{
	vnl_matrix<double> q(x.rows()*y.size(),x.cols());
	int counter;
	for(int i = 0; i < x.cols() ; ++i)
	  {
		counter = 0;
		for(int j = 0; j < x.rows() ; ++j)
		 {
			for(int k = 0; k < y.size() ; ++k)	
			{
				q(counter,i) = x(j,i)*y(k);
				counter++;
			}
		}
	}
	return q;
}

