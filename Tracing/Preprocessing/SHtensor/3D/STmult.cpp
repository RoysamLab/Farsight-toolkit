
#include "mex.h"
#include "matrix.h"
#include <math.h>
#include "SphericalTensor.cpp"


template <class T> inline void swap(T& x, T& y) { T t=x; x=y; y=t; }


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
 	
    if(nrhs != 3 && nrhs != 4) {
	printf("\nUsage: Df = STmult(ST1,ST2,J,factor)\n\n",nrhs);
	printf(" Computes Spherical product\n");
	printf(" Parameters:\n");
	printf("   ST1,ST2 - Spherical Tensors to be multiplied \n");
	printf("   J - spherical rank of result\n");
    printf("   factor - result is multiplied with scalar factor (optional) \n\n");
    return;
	} else if(nlhs>1) {
	printf("Too many output arguments\n");
    return;
	}

    int pcnt = 0;

    const mxArray *Img1;
    Img1 = prhs[pcnt++];       
    const int *dims = mxGetDimensions(Img1);
    const int numdim = mxGetNumberOfDimensions(Img1); //TODO CHECK
    printf("%d\n",mxGetNumberOfDimensions(Img1));
    int L1 = dims[0]/2;
    REAL *img1 = (REAL*) mxGetData(Img1);
    

    const mxArray *Img2;
    Img2 = prhs[pcnt++];           
    const int *dims2 = mxGetDimensions(Img2);
    int L2 = dims2[0]/2;
    REAL *img2 = (REAL*) mxGetData(Img2);

    const mxArray *destsz= prhs[pcnt++];
    int destsizeL = 1 + (int) *mxGetPr(destsz);

    REAL factor = 1.0;	
    if (nrhs == 4)
	{
	const mxArray *fac= prhs[pcnt++];;
	factor = *(REAL*)mxGetData(fac);

	}
    


    int *ndims = (int *)mxCalloc(numdim,sizeof(int));
    ndims[0] = destsizeL*2; 
    int totsiz = 1;
    for(int k =1;k <numdim;k++)
    {
	    ndims[k] = dims[k];
        totsiz *= dims[k];
    }

    plhs[0] = mxCreateNumericArray(numdim,ndims,mxGetClassID(Img1),mxREAL);
    REAL *result = (REAL*) mxGetData(plhs[0]);
//    printf("I came till here\n");
    STmultiply_withCG(img1,L1,img2,L2,result,destsizeL,totsiz,factor);
}
