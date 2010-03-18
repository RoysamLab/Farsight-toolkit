
#include "mex.h"
#include "matrix.h"
#include <math.h>
#include "SphericalTensor.cpp"



void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
 	
    if(nrhs != 1 && nrhs != 2  && nrhs != 3) {
	printf("\nUsage: Df = STderivDown(f,difftyp,factor)\n\n",nrhs);
	printf(" Computes Spherical Down-Derivative\n");
	printf(" Parameters:\n");
	printf("   f - 3D input image of type REAL \n");
	printf("   difftyp - central (0) / forward (-1) / backward (1)   (optional)\n");
    printf("   factor - multiply result with factor   (optional)\n\n");
    return;
	} else if(nlhs>1) {
	printf("Too many output arguments\n");
    return;
	}

    int pcnt = 0;
    const mxArray *Img;
    Img = prhs[pcnt++];       
    const int numdim = mxGetNumberOfDimensions(Img);
    const int *dims = mxGetDimensions(Img);
    int L = dims[0]/2-1;

    if (numdim != 4)
	return;



    int diff_type = 0;  
    if (nrhs >= 2)
	{
		const mxArray *typ;
		typ = prhs[pcnt++];       

		diff_type = (int) *mxGetPr(typ);
		if (diff_type < -1 || diff_type > 1)
			diff_type = 0;
	}



    REAL factor = 1.0;
    if (nrhs == 3)
	{
    		const mxArray *mxFactor = prhs[pcnt++];       
		if (mxGetClassID(mxFactor) == mxDOUBLE_CLASS)
			factor = (REAL) *((double*) mxGetData(mxFactor));
		else if (mxGetClassID(mxFactor) == mxSINGLE_CLASS)
			factor = (REAL) *((float*) mxGetData(mxFactor));
		else 
			return;

	}

    REAL *img = (REAL*) mxGetData(Img);


    int ndims[4];
    ndims[0] = 2*L; ndims[1] = dims[1]; ndims[2] = dims[2]; ndims[3] = dims[3];
    plhs[0] = mxCreateNumericArray(4,ndims,mxGetClassID(Img),mxREAL);
    REAL *res = (REAL*) mxGetData(plhs[0]);


    if (diff_type == -1)
    	STderivRealDownForward(img,&(ndims[1]),L,res,factor);
    else if (diff_type == 1)
    	STderivRealDownBackward(img,&(ndims[1]),L,res,factor);
    else if (diff_type == 0)
    	STderivRealDown(img,&(ndims[1]),L,res,factor);

}


