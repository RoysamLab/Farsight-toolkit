



# include "agf.h"

double (* global_metric2) (double *, double *, long m) = &metric2;
///////////////////////////////////////////////////////////////////////////////////////////////////////
// CALCULATES THE PDF USING THE GAUSSIAN KERNEL
///////////////////////////////////////////////////////////////////////////////////////////////////////
double agf_calc_pdf(double **mat, long D, long n, double *vec, double var_0, 
					double Wc, agf_diag_param *diag_param) {
	
	double *d2;			//the distances (squared)
	double var_f;			//final value of the filter width (as variance)
	double tw;			//total weight
	double *weight;		//the current value for the weights
	double n1, norm;		//normalisation coeff.
	double pdf;			//final calculated value of pdf
	
	//first we calculate	all the distances:
	d2=new double[n];
	for (long i=0; i<n; i++) {
		d2[i]=(*global_metric2) (vec, mat[i], D);
	}
	
	//calculate the weights using the central "engine":
	weight=new double[n];
	diag_param->nd=agf_calc_w(d2, n, Wc, var_0, weight, var_f);
	tw=0;
	for (long i=0; i<n; i++) tw+=weight[i];
	
	//use the final filter width to normalize the pdf:
	n1=sqrt(var_f*M_PI*2);
	norm=1;
	for (long i=0; i<D; i++) norm*=n1;
	
	//norm=pow(var_f*M_PI*2, m/2.);
	//printf("var_f=%g, tw=%g, norm=%g\n", var_f, tw, norm);
	
	pdf=tw/norm/n;
	
	//set the diagnostic parameters:
	//diag_param->f=0;
	diag_param->W=tw;
	
	delete [] d2;
	delete [] weight;
	
	return pdf;

}


///////////////////////////////////////////////////////////////////////////////////////////////////////
// USED IN LIBAGF
///////////////////////////////////////////////////////////////////////////////////////////////////////
double metric2(double *v1, double *v2, long m)
{
	double d;
	double diff;
	d = -1.0;	
	d=0;
	
	for (long i=0; i<m; i++) {
		diff=v2[i]-v1[i];
		d+=diff*diff;
		// D inf instead of Cartesian
	}
	
	return  (d) ;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////
// GIVEN A SET OF DISTANCES, CALCULATES THE WEIGHTS OF AN ADAPTIVE GAUSSIAN FILTER
///////////////////////////////////////////////////////////////////////////////////////////////////////
long agf_calc_w(double *d2,     //distances squared
				long k,           //number of distances
				double Wc,         //objective total weight
				double var_0,	//initial filter width
				double *weight,	//returned weights
				double &var_f)         //returned final filter variance
{
	long nd;                      //number of squarings of the weights
	double tw_old;                 //previous value of the total weight
	double tw;                     //current value of the total weight
	
	double wtint;                  //for interpolation of final weights
	
	//calculate the weights:
	for (long i=0; i<k; i++)
	{
		weight[i]=exp(-d2[i]/var_0/2);
	}
		
	
	//repeatedly square the weights until the total is less than the threshold:
	tw=0;
	for (long i=0; i<k; i++) tw+=weight[i];
	nd=0;
	
	do {
		tw_old=tw;
		tw=0;
		for (long i=0; i<k; i++) {
			weight[i]*=weight[i];
			tw+=weight[i];
		}
		nd++;
	} while (tw > Wc && tw < tw_old);
	
	//interpolate between the current and previous weight:
	wtint=(log(tw_old)-log(Wc))/(log(tw_old)-log(tw))/2+0.5;
	
	//calculate the final filter width:
	var_f=var_0/(1 << nd)/wtint;
	
	tw=0;
	for (long i=0; i<k; i++) {
		//weight[i]=pow(weight[i], wtint);
		//instead of interpolatting, why don't we just calculate the weights directly?
		//(the differences are insignificant, but still...)
		weight[i]=exp(-d2[i]/var_f/2);
	}
	
	//return the number of iterations as a diagnostic parameter:
	return nd;
	
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// CALCULATE THE GRADIENT OF THE WEIGHTS
///////////////////////////////////////////////////////////////////////////////////////////////////////
void agf_grad_w(double **x,		//matrix of samples
				long n,			//# of dimensions
				double *xvec,		//test point
				double *w,		//weights
				double *d2,		//squared distances
				long k,			//number of weights and indices
				double var_f,		//filter width
				double **dwdx)  		//returned gradient vectors
{
	double dwcnum, dwcden;
	
	for (long j=0; j<n; j++) {
		dwcnum=0;
		dwcden=0;
		for (long i=0; i<k; i++) {
			dwcnum+=w[i]*(x[i][j]-xvec[j]);
			dwcden+=d2[i]*w[i];
		}
		for (long i=0; i<k; i++) {
			dwdx[i][j]=w[i]*(x[i][j]-xvec[j]-d2[i]*dwcnum/dwcden)/var_f;
		}
	}
	
}


//template<class type>
//void kleast(type *data, long n, long k, type *result) 
//{
//	
//	tree_lg<type> sorter;
//	for (long i=0L; i<k; i++) 
//	{
//		sorter.add(data[i]);
//	}
//	
//	
//	for (long i=k; i<n; i++) {
//		sorter.add(data[i]);
//		sorter.delete_greatest();
//	}
//	
//	sorter.decompose(result, k);
//}
//
//
//template<class type>
//tree_lg<type>::tree_lg() {
//	trunk=NULL;
//	n=0;
//}
//
//template<class type>
//tree_lg<type>::~tree_lg() {
//	delete_el(trunk);
//}
//
//
//template<class type>
//void tree_lg<type>::delete_el(tree_lg_el<type> *tel) {
//	if (tel == NULL) return;
//	delete_el(tel->right);
//	delete_el(tel->left);
//	delete tel;
//}
//
//template<class type>
//void tree_lg<type>::decompose(tree_lg_el<type> *t, type *sarr, long nd, long &iter) {
//	if (t == NULL) return;
//	if (iter >= nd) return;
//	decompose(t->left, sarr, nd, iter);
//	sarr[iter]=t->value;
//	iter++;
//	decompose(t->right, sarr, nd, iter);
//	
//}
//
//template<class type>
//long tree_lg<type>::add(type data) {
//	tree_lg_el<type> *current;
//	
//	if (trunk==NULL) {
//		trunk=new tree_lg_el<type>;
//		trunk->value=data;
//		trunk->left=NULL;
//		trunk->right=NULL;
//		n=1;
//		return n;
//	}
//	
//	current=trunk;
//	for (;;) {
//		if (data < current->value) {
//			if (current->left == NULL) {
//				current->left=new tree_lg_el<type>;
//				current=current->left;
//				break;
//			} else {
//				current=current->left;
//			}
//		} else {
//			if (current->right == NULL) {
//				current->right=new tree_lg_el<type>;
//				current=current->right;
//				break;
//			} else {
//				current=current->right;
//			}
//		}
//	}
//	
//	current->value=data;
//	current->right=NULL;
//	current->left=NULL;
//	
//	n++;
//	
//	return n;
//}
//
//template<class type>
//long tree_lg<type>::nel() {
//	return n;
//}
//
//template<class type>
//void tree_lg<type>::decompose(type *sarr, long nd) {
//	long i=0;
//	decompose(trunk, sarr, nd, i);
//}
//
//
//template<class type>
//long tree_lg<type>::delete_least() {
//	tree_lg_el<type> *current, *previous;
//	
//	if (trunk == NULL) return n;
//	
//	if (trunk->left == NULL) {
//		if (trunk->right == NULL) {
//			delete trunk;
//			trunk=NULL;
//			n=0;
//		} else {
//			current=trunk->right;
//			delete trunk;
//			trunk=current;
//		}
//	} else {
//		current=trunk;
//		while (current->left != NULL) {
//			previous=current;
//			current=current->left;
//		}
//		previous->left=current->right;
//		delete current;
//		
//	}
//	n--;
//	return n;
//}
//
//template<class type>
//long tree_lg<type>::delete_greatest() {
//	tree_lg_el<type> *current, *previous;
//	
//	if (trunk == NULL) return n;
//	
//	if (trunk->right == NULL) {
//		if (trunk->left == NULL) {
//			delete trunk;
//			trunk=NULL;
//			n=0;
//		} else {
//			current=trunk->left;
//			delete trunk;
//			trunk=current;
//		}
//	} else {
//		current=trunk;
//		while (current->right != NULL) {
//			previous=current;
//			current=current->right;
//		}
//		previous->right=current->left;
//		delete current;
//		
//	}
//	n--;
//	return n;
//}
//
//template class tree_lg<double>;
