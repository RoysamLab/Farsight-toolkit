#ifndef MSTRINGIFY
#define MSTRINGIFY(A) A
#endif

MSTRINGIFY(

__kernel void fibonacci(__global unsigned long* out)
{
	int iGID = get_global_id(0);
	
	float phi = (1 + sqrt(5.0)) / 2.0;
	
	out[iGID] = round(native_powr(phi,iGID+1) / sqrt(5.0)); 
}

);