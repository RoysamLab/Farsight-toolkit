#ifndef MSTRINGIFY
#define MSTRINGIFY(A) A
#endif

MSTRINGIFY(

#pragma OPENCL EXTENSION cl_khr_fp64 : enable \n //note that \n is needed because of stringification because C preprocessor expects newline after macro

//This kernel only accurate up to the 70th Fibonacci number
__kernel void fibonacci (__global unsigned long* out)
{
	int iGID = get_global_id(0);
	
	double phi = (1 + sqrt(5.0)) / 2.0;
	
	out[iGID] = round(powr(phi,iGID+1) / sqrt(5.0)); 
}

);     