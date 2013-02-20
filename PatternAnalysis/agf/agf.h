


//This code has been borrowed from libAGF library

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#include "math.h"


struct agf_diag_param {
		long nd;		//number of squarings
		double f;		//ratio of min. weight to max.
		double W;		//total weight  
					  };

struct agf_command_opts {
	double var_0;
	long k;
	double Wc;
	long n;
	double tol;
	double ftest;
	int normflag;
	int errflag;
	double rthresh;
	double pmin;
	char *normfile;
	int algtype;
	int metrictype;
	int jointflag;
};


double agf_calc_pdf(double **mat, long D, long n, double *vec, double var_0, double Wc, agf_diag_param *diag_param);   
long agf_calc_w(double *d2,     //distances squared
				long k,           //number of distances
				double Wc,         //objective total weight
				double var_0,	//initial filter width
				double *weight,	//returned weights
				double &var_f) ;         //returned final filter variance
double agf_calc_pdf(double **mat, long D, long n, double *vec, double var_0, 
				   long k, double Wc, agf_diag_param *diag_param);

//calculate the gradient of the weights:
 void agf_grad_w(double **x,		//matrix of samples
				long n,			//# of dimensions
				double *xvec,		//test point
				double *w,		//weights
				double *d2,		//squared distances
				long k,			//number of weights and indices
				double var_f,		//filter width
				double **dwdx);  		//returned gradient vectors
double metric2(double *v1, double *v2, long m);