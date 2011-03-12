#ifndef MSTRINGIFY
#define MSTRINGIFY(A) A
#endif

MSTRINGIFY(

__kernel void fibonacci(__global float* out)
{
	int iGID = get_global_id(0);
	
	//float phi = (1 + native_sqrt(5)) / 2;

	//out[iGID] = (pow(phi, iGID) - pow((-1 / phi), iGID)) / native_sqrt(5);

	out[iGID] = iGID;
}

);