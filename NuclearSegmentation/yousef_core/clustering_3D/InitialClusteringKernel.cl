#ifndef MSTRINGIFY
#define MSTRINGIFY(A) A
#endif

MSTRINGIFY(

#pragma OPENCL EXTENSION cl_khr_fp64 : enable \n //note that \n is needed because of stringification because C preprocessor expects newline after macro

__kernel void InitialClusteringKernel (__global float* im_vals, __global unsigned short* local_max_vals, __global unsigned short* max_response_r, __global unsigned short* max_response_c, __global unsigned short* max_response_z , int r, int c, int z, int scale_xy, int scale_z)
{
	int iGID = get_global_id(0);
	
	int rem = ((long)iGID) % (r*c);
	int k1 = ((int)iGID-rem) / (r*c); 
	int j1 = ((long)rem) % c;
	int i1 = (rem-j1)/c;

	int min_r = (int) max((double)(0.0),(double)(i1-scale_xy));
	int min_c = (int) max((double)(0.0),(double)(j1-scale_xy));
	int min_z = (int) max((double)(0.0),(double)(k1-scale_z));
	int max_r = (int) min((double)(r-1),(double)(i1+scale_xy));
	int max_c = (int) min((double)(c-1),(double)(j1+scale_xy));                         
	int max_z = (int) min((double)(z-1),(double)(k1+scale_z));

	if(local_max_vals[(k1*r*c)+(i1*c)+j1] == 0) //if current pixel is not a seed point			
	{
		float mx = im_vals[(min_z*r*c)+(min_r*c)+min_c];//A[r1][c1][z1];
		
		max_response_r[i1 * (c * z) + j1 * z + k1] = min_r;
		max_response_r[i1 * (c * z) + j1 * z + k1] = min_c;
		max_response_r[i1 * (c * z) + j1 * z + k1] = min_z;
	    
		for(int i= min_r; i<= max_r; i++)
		{
			for(int j= min_c; j <= max_c; j++)
			{
				for(int k = min_z; k <= max_z; k++)
				{
					if(im_vals[(k*r*c)+(i*c)+j] >= mx)
					{
						mx = im_vals[(k*r*c)+(i*c)+j];//A[i][j][k];

						max_response_r[i1 * (c * z) + j1 * z + k1] = i;
						max_response_c[i1 * (c * z) + j1 * z + k1] = j;
						max_response_z[i1 * (c * z) + j1 * z + k1] = k;
					}
				}
			}
		}                                          
	}
}

);     