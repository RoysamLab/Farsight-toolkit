



#include "mclr.h"


MCLR::MCLR()
{

}

MCLR::~MCLR()
{

}


// x contains the training features
// y contains the training labels
void MCLR::Initialize(vnl_matrix<double> data,double c,std::string str )
{
	
	int train_counter =0;
	int test_counter = 0;
	//vnl_vector<double> classes = data.get_column(data.cols()-2);
	validation = str;

	vnl_vector<double> classes;
	if(validation=="ground_truth")
	{
		classes = data.get_column(data.cols()-2);
	//	y_ground_truth = data.get_column(data.cols()-3);
	}

	for(int i = 0; i< classes.size() ; ++i)
	{
		if(classes(i) == -1)
			test_counter++;
	}
	
	
	train_data.set_size(data.rows()-test_counter,data.cols());
	test_data.set_size(test_counter,data.cols());
	
	if(validation=="ground_truth")
	{
		y_ground_truth.set_size(test_counter);
	}
	

	test_counter = 0;

	for(int i= 0;i<data.rows();++i)
	{
		if(data(i,data.cols()-2)!=-1)
		{	
			train_data.set_row(train_counter,data.get_row(i));
			++train_counter;
		}
		else
		{
			test_data.set_row(test_counter,data.get_row(i));
			if(validation=="ground_truth")
			{
				y_ground_truth(test_counter) = data(i,data.cols()-3);
			}
			++test_counter;
		}
	}

	// Transpose : Features->Rows ; Samples -> Columns
	x = train_data.transpose();// Feature matrix of training samples only ; does not contain unlabeled sample data
	y = train_data.get_column(train_data.cols()-2);// label value ; 


	if(validation=="ground_truth")
	// remove the the rows corresponding to id and class
	x = x.get_n_rows(0,x.rows()-3);
	else
	x = x.get_n_rows(0,x.rows()-2);	



	
	//m.method = "direct"; // No kernel is the default
	m.sparsity_control = c; // greater the value , more is the sparsity
	m.kstd_level = 1; // kernel sigma val

	// Stop Condition
	// Used in active label selection
	stop_cond.set_size(2);
	stop_cond[0] = 1e-5;
	stop_cond[1] = 1e-5;
	delta = 1e-9; // used for approximating |w| 
	
	no_of_features = x.rows();
	// We assume that the classes are numbered from 1 to "n" for n classes
	// so y.maxval - .minval +1 gives number of classes
	no_of_classes = y.max_value()-y.min_value()+1;
	class_vector.set_size(no_of_classes);

	m.w.set_size(no_of_features+1,no_of_classes);
	m.w.fill(0);

	gradient_w.set_size(no_of_features+1,no_of_classes);
	hessian.set_size((no_of_features+1)*(no_of_classes-1),(no_of_features+1)*(no_of_classes-1)); /// Check the size again!!!
	
	// g :the objective function value
	g = 0;
	// Class vector
	for(int i=1;i<=no_of_classes;++i)
	{
		class_vector.put(i-1,i); 
	}
}


//Converts vtkTable to Vnl_Matrix
vnl_matrix <double> MCLR::tableToMatrix(vtkSmartPointer<vtkTable> table)
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

//	FeatsMatrix = Normalize_Feature_Matrix(FeatsMatrix);
	return FeatsMatrix;
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

		for(int j = 0; j< x.cols() ; ++j)
		{
			ith_row(j) = ith_row(j)*(1-ith_row(j));
		}


		vnl_diag_matrix<double> diagonal_matrix(ith_row);

		vnl_matrix<double> temp_matrix =  -data_with_bias*diagonal_matrix*data_with_bias.transpose();
		temp_hessian.update(temp_matrix,(i-1)*(no_of_features+1),(i-1)*(no_of_features+1));


		//
		//for(int k = 0; k< class_vector.size() ; ++k)
		//	std::cout<<class_vector(k)<<"--------------"<<std::endl;	
			
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

		//for(int k = 0; k< class_vector.size() ; ++k)
		//	std::cout<<class_vector(k)<<"--------------"<<std::endl;
		
		for(int k = 0; k< all_class_but_i.size() ; ++k)
		{	
			ith_row  = f.get_row(i-1);
			vnl_vector<double> kth_row  = f.get_row(all_class_but_i(k)-1);
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
		  for(int j=0;j<m.w.cols();++j)	
		  {
			w_column_ordered(count) = (m.sparsity_control)*(delta/pow(sqrt(pow(w_column_ordered(count),2)+delta),3));
			count++;
		  }
		}
		
		vnl_diag_matrix<double> diagonal_matrix(w_column_ordered);
		temp_hessian = temp_hessian - diagonal_matrix;
		//temp_hessian.extract(hessian,no_of_features+2,no_of_features+2);
		hessian = temp_hessian;

		//std::cout<<"------------------------------------------------------"<<std::endl;
		//std::cout<<hessian.get_row(0)<<std::endl;
		//std::cout<<"------------------------------------------------------"<<std::endl;
		//std::cout<<hessian.get_row(115)<<std::endl;
		//std::cout<<"------------------------------------------------------"<<std::endl;

}	


void MCLR::Ameliorate_Hessian_Conditions()
{
	//Compute the eigen vectors of the covariance matrix
	vnl_symmetric_eigensystem<double> eig(hessian);
	double ratio = fabs(eig.get_eigenvalue((no_of_features+1)*2))/fabs(eig.get_eigenvalue(1));
	
	while(ratio>1e9)
	{
	   vnl_diag_matrix<double> diag((no_of_classes-1)*(no_of_features+1),Compute_Mean_Abs_Eig(eig));
	   hessian = hessian + diag;
	   vnl_symmetric_eigensystem<double> eig(hessian);
	}

}


double MCLR::Compute_Mean_Abs_Eig(vnl_symmetric_eigensystem<double> eig)
{
	double sum = 0; 	
	for(int i =0; i<(no_of_features+1)*(no_of_classes-1);++i)
	{
		sum = sum + fabs(eig.get_eigenvalue(i));
	}	
	sum = sum/((no_of_features+1)*(no_of_classes-1)) ;
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

	//vnl_matrix<double> epow_matrix = m.w.transpose()*x_with_bias;
	//f.set_size(epow_matrix.rows(),epow_matrix.cols());

	//for(int i=0;i<epow_matrix.rows();++i)
	//{
	//  for(int j=0;j<epow_matrix.cols();++j)
	//   {
	//	  f(i,j) = exp(epow_matrix(i,j));
	//  }
	//}

	f = Get_F_Matrix(data_with_bias,w_temp);

	//std::cout<<"--------------------------"<<std::endl;
	//std::cout<<f.get_n_columns(0,2)<<std::endl;
	//std::cout<<"--------------------------"<<std::endl;


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
	
	//std::cout<<m.w(0,0)<<" --- " << m.w(0,1)<<"---"<< m.w(0,2)<<std::endl;
	//std::cout<<m.w(38,0) <<" --- " << m.w(38,1) <<"---"<< m.w(38,2) <<std::endl;


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
	vnl_vector<double> g_3_it(3,0);	// g vals for the last three

	for(int i =0;i<1e10;++i)
	{	
		//Get the gradient
		Get_Gradient(Add_Bias(x));
		//Get the hessian
		Get_Hessian(Add_Bias(x));
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
		

		//std::cout<<"----------------------"<<std::endl;
		//std::cout<<hessian.get_column(0)<<std::endl;
		//std::cout<<"----------------------"<<std::endl;
		//std::cout<<hessian.get_column(1)<<std::endl;
		//std::cout<<"----------------------"<<std::endl;
		//std::cout<<hessian.get_column(2)<<std::endl;
		//std::cout<<"----------------------"<<std::endl;

		g_3_it(i%3) = g;	

		if(i>1)
			diff_g_3_it(i%3) = g_3_it(i%3) - g_3_it((i-1)%3);
    
		if (i>3 && (diff_g_3_it.one_norm()/3)/(g_3_it.one_norm()/3) < stop_cond(1))
			break;
	}


	m.FIM = hessian*-1;
	m.CRB = vnl_matrix_inverse<double>(hessian);
	

	// quirk of vnl ...
	m.CRB = m.CRB*-1;
	//std::cout<<vnl_trace(m.CRB)<<std::endl;


    FILE * fpin1 = FDeclare2("C:\\Users\\rkpadman\\Documents\\MATLAB\\crb.txt", "", 'w');

	for(int abc =0;abc<hessian.rows();++abc)
	{	
		vnl_vector<double> temp = hessian.get_row(abc);
		for(int abcd =0;abcd<hessian.cols();++abcd)
		{
			fprintf(fpin1, "%lf    ", temp[abcd]);
		}
		fprintf(fpin1, "\n");
	}


		//std::cout<<"----------------------"<<std::endl;
		//std::cout<<m.CRB.get_column(0)<<std::endl;
		//std::cout<<"----------------------"<<std::endl;
		//std::cout<<m.CRB.get_column(1)<<std::endl;
		//std::cout<<"----------------------"<<std::endl;
		//std::cout<<m.CRB.get_column(2)<<std::endl;
		//std::cout<<"----------------------"<<std::endl;
	fclose(fpin1);


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


void MCLR::Update_Train_Data(int query)
{
	vnl_vector<double> queried_sample = test_data.get_row(query);
	vnl_matrix<double> sub_test_matrix(test_data.rows()-1,test_data.cols());

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

		if(validation=="ground_truth")
			queried_sample = queried_sample.extract(queried_sample.size()-3,0);
		else
			queried_sample = queried_sample.extract(queried_sample.size()-2,0);
		

		vnl_matrix<double> temp_x = x;
		x.set_size(x.rows(),x.cols()+1);
		
		x.update(temp_x,0,0);
		x.set_column(x.cols()-1,queried_sample);
		Get_Label_Sample(query);
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


vnl_vector<double> MCLR::Newton_Direction(vnl_matrix<double> hessian_matrix,vnl_vector<double> grad_vector )
{
	vnl_matrix<double> inv_hessian =  vnl_matrix_inverse<double>(hessian_matrix);
	

	
	vnl_vector<double> newton_direction = -inv_hessian*grad_vector;
	//std::cout<<"----------------------"<<std::endl;
	//std::cout<<inv_hessian.get_row(0)<<std::endl;
	//std::cout<<"----------------------"<<std::endl;
	//newton_direction = newton_direction.normalize();
	//std::cout<<"----------------------"<<std::endl;
	//std::cout<<newton_direction<<std::endl;
	//std::cout<<"----------------------"<<std::endl;

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