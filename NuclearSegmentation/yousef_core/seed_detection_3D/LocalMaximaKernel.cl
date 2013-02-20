#ifndef MSTRINGIFY
#define MSTRINGIFY(A) A
#endif

MSTRINGIFY(

#pragma OPENCL EXTENSION cl_khr_fp64 : enable \n //note that \n is needed because stringification converts all whitespace into a single space and C preprocessor expects newline after macro to delimit the macro

__kernel void LocalMaximaKernel (__global float* im_vals, __global unsigned short* out1, int r, int c, int z, double scale_xy, double scale_z)
{
	int iGID = get_global_id(0); //global index
		
	//calculate r, c, z indices as i, j, k from global index
	int rem = ((long)iGID) % (r*c);
	int k = ((int)iGID-rem) / (r*c); 
	int j = ((long)rem) % c;
	int i = (rem-j)/c;
	
	//calculate bounds
	int min_r = (int) max(0.0,i-scale_xy);
	int min_c = (int) max(0.0,j-scale_xy);
	int min_z = (int) max(0.0,k-scale_z);
	int max_r = (int)min((double)r-1,i+scale_xy);
	int max_c = (int)min((double)c-1,j+scale_xy);                         
	int max_z = (int)min((double)z-1,k+scale_z);                         
	
	//get the intensity maximum of the bounded im_vals
	float mx = im_vals[(min_z*r*c)+(min_r*c)+min_c];
    
	for(int i = min_r; i <= max_r; i++)
    {
        for(int j = min_c; j <= max_c; j++)
        {
			for(int k = min_z; k <= max_z; k++)
			{				
				if(im_vals[(k*r*c)+(i*c)+j] > mx)
					mx = im_vals[(k*r*c)+(i*c)+j];
			}
        }
    }
	
	//if the current pixel is at the maximum intensity, set it to 255 in out1 (seedImagePtr), else set it to 0
	if(im_vals[iGID] == mx)    
		out1[iGID]=255;
	else
		out1[iGID]=0;
}

);     