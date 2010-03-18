
#include "mex.h"
#include "matrix.h"
#include <math.h>
#include "SphericalTensor.cpp"


#define mxREAL_CLASS mxDOUBLE_CLASS
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
 	
    if(nrhs != 4) {
	printf("\nUsage: Df = STderivative(f)\n\n",nrhs);
	printf(" Computes xxx\n");
	printf(" Parameters:\n");
	printf("   f - 2D input image of type REAL \n");
	printf(" Return value Df contains xxx.\n\n");
    return;
	} else if(nlhs>1) {
	printf("Too many output arguments\n");
    return;
	}

    int pcnt = 0;

    const mxArray *Img1;
    Img1 = prhs[pcnt++];       
    const int *dims = mxGetDimensions(Img1);
    int L1 = dims[0]/2;
    REAL *img1 = (REAL*) mxGetData(Img1);

    const mxArray *Img2;
    Img2 = prhs[pcnt++];           
    const int *dims2 = mxGetDimensions(Img2);
    int L2 = dims2[0]/2;
    REAL *img2 = (REAL*) mxGetData(Img2);

    const mxArray *destsz= prhs[pcnt++];
    int destsizeL = 1 + (int) *mxGetPr(destsz);
	
    const mxArray *indexcg = prhs[pcnt++];       
    const int cnt = mxGetN(indexcg);

    REAL *_indexcg = (REAL*) mxGetData(indexcg);
    REAL *cg = new REAL[cnt];
    int *indextrip = new int[3*cnt];
    for (int i = 0;i < cnt;i++)
	{	
		indextrip[3*i] =  2* (int) _indexcg[4*i];
		indextrip[3*i+1] = 2* (int) _indexcg[4*i+1];
		indextrip[3*i+2] = 2* (int) _indexcg[4*i+2];
		cg[i] = _indexcg[4*i+3];
	}


    int ndims[4];
    ndims[0] = 2*destsizeL; ndims[1] = dims[1]; ndims[2] = dims[2]; ndims[3] = dims[3];
    plhs[0] = mxCreateNumericArray(4,ndims,mxGetClassID(Img1),mxREAL);
    REAL *result = (REAL*) mxGetData(plhs[0]);

    STmultiplyFourier(img1,L1,img2,L2,result,destsizeL, &(ndims[1]),  indextrip, cg, cnt);

}


