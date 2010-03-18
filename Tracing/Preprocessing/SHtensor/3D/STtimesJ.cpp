

#include "mex.h"
#include "matrix.h"
#include <math.h>
#include "SphericalTensor.cpp"



void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
 	
    if(nrhs != 1) {

       return;
	} else if(nlhs>1) {
	printf("Too many output arguments\n");
    return;
	}

    int pcnt = 0;
    const mxArray *M;
    M = prhs[pcnt++];       
    const int *dims = mxGetDimensions(M);
    const int numdim = mxGetNumberOfDimensions(M);
    REAL *SHin = (REAL*) mxGetData(M);

    if (dims[0]%2 == 1)
	{
		mexPrintf("First dim. must be even!\n");
		return;
	}


    int totsiz = 1;
    for(int k = 1; k < numdim;k++)
	{
	totsiz *= dims[k];
	}

    plhs[0] = mxCreateNumericArray(numdim,dims,mxGetClassID(M),mxREAL);
    REAL *out = (REAL*) mxGetData(plhs[0]);

    int r2 = dims[0];
    for (int j = 0; j < totsiz;j++)
    {
	for (int k = 0; k < r2;k+=2)
	{
		int idx = j*r2+k;
		out[idx] = -SHin[idx+1];
		out[idx+1] = SHin[idx];
	}
    }
}




