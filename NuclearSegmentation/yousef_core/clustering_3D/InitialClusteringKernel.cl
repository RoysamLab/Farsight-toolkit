#ifndef MSTRINGIFY
#define MSTRINGIFY(A) A
#endif

MSTRINGIFY(

#pragma OPENCL EXTENSION cl_khr_fp64 : enable \n //note that \n is needed because of stringification because C preprocessor expects newline after macro

//This kernel only accurate up to the 70th Fibonacci number
__kernel void InitialClusteringKernel (__global unsigned long* out)
{

);     