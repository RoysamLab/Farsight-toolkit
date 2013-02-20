#include <cuda.h>
#include <iostream>

using namespace std;

__global__ void InitialClusteringKernel_CUDA (float* im_vals, unsigned short* max_response_r, unsigned short* max_response_c, unsigned short* max_response_z , int r, int c, int z, int scale_xy, int scale_z, int offset)
{
	int iGID = blockIdx.x * blockDim.x + threadIdx.x + offset; //global index
	
	if (iGID >= r * c * z)
		return;

	int rem = ((long)iGID) % (r*c);
	int k1 = ((int)iGID-rem) / (r*c); 
	int j1 = ((long)rem) % c;
	int i1 = (rem-j1)/c;

	int min_r = (int) max((float)(0.0),(float)(i1-scale_xy));
	int min_c = (int) max((float)(0.0),(float)(j1-scale_xy));
	int min_z = (int) max((float)(0.0),(float)(k1-scale_z));
	int max_r = (int) min((float)(r-1),(float)(i1+scale_xy));
	int max_c = (int) min((float)(c-1),(float)(j1+scale_xy));                         
	int max_z = (int) min((float)(z-1),(float)(k1+scale_z));

	//If we are running on the GPU, it makes no sense to load this big array just to save computation since seed point should be the maximum anyways
	//if(local_max_vals[(k1*r*c)+(i1*c)+j1] == 0) //if current pixel is not a seed point		
	//{
		float mx = im_vals[(min_z*r*c)+(min_r*c)+min_c];//A[r1][c1][z1];
		
		//Do not access arrays in a hot loop if you can just just do all the updates at the end, much faster to use registers
		int r_temp = min_r;
		int c_temp = min_c;
		int z_temp = min_z;	    

		float im_vals_temp;
		
		for(int i= min_r; i<= max_r; i++)
		{
			for(int j= min_c; j <= max_c; j++)
			{
				for(int k = min_z; k <= max_z; k++)
				{
					im_vals_temp = im_vals[(k*r*c)+(i*c)+j];

					if( im_vals_temp >= mx)
					{
						mx = im_vals_temp;

						r_temp = i;
						c_temp = j;
						z_temp = k;
					}
				}
			}
		} 

		max_response_r[i1 * (c * z) + j1 * z + k1] = r_temp;
		max_response_c[i1 * (c * z) + j1 * z + k1] = c_temp;
		max_response_z[i1 * (c * z) + j1 * z + k1] = z_temp;
	//}
}


extern "C"
void initialClustering_CUDA (float* im_vals, unsigned short* local_max_vals, unsigned short* max_response_r, unsigned short* max_response_c, unsigned short* max_response_z , int r, int c, int z, int scale_xy, int scale_z)
{
	cout << "Entering initialClustering_CUDA" << endl;
	
	cudaError_t errorcode;

	float* dev_im_vals; 
	//unsigned short* dev_local_max_vals;
	unsigned short* dev_max_response_r;
	unsigned short* dev_max_response_c;
	unsigned short* dev_max_response_z;

	size_t free_mem, total_mem;
	cudaMemGetInfo(&free_mem, &total_mem);

	cout << free_mem / (double)(1024 * 1024) << " " << total_mem / (double)(1024 * 1024) << endl;

	cout << "Allocating " << (sizeof(*im_vals) * r * c * z)/(double)(1024*1024) << " MB of memory on GPU for im_vals" << endl;
	//cout << "Allocating " << (sizeof(*local_max_vals) * r * c * z)/(double)(1024*1024) << " MB of memory on GPU for local_max_vals" << endl;
	cout << "Allocating " << (sizeof(*max_response_r) * r * c * z)/(double)(1024*1024) << " MB of memory on GPU for max_response_r" << endl;
	cout << "Allocating " << (sizeof(*max_response_c) * r * c * z)/(double)(1024*1024) << " MB of memory on GPU for max_response_c" << endl;
	cout << "Allocating " << (sizeof(*max_response_z) * r * c * z)/(double)(1024*1024) << " MB of memory on GPU for max_response_z" << endl;
	
	//Allocate memory on device
	errorcode = cudaMalloc((void**) &dev_im_vals, r * c * z * sizeof(*im_vals));
	//errorcode = cudaMalloc((void**) &dev_local_max_vals, r * c * z * sizeof(*local_max_vals));
	errorcode = cudaMalloc((void**) &dev_max_response_r, r * c * z * sizeof(*dev_max_response_r));
	errorcode = cudaMalloc((void**) &dev_max_response_c, r * c * z * sizeof(*dev_max_response_c));
	errorcode = cudaMalloc((void**) &dev_max_response_z, r * c * z * sizeof(*dev_max_response_z));

	//cout << errorcode << endl;

	//Copy host memory contents to device contents
	cudaMemcpy(dev_im_vals, im_vals, r * c * z * sizeof(*im_vals), cudaMemcpyHostToDevice);
	//cudaMemcpy(dev_local_max_vals, local_max_vals, r * c * z * sizeof(*local_max_vals), cudaMemcpyHostToDevice);

	//prefer 48 KB L1
	CUresult drivererrorcode = cuCtxSetCacheConfig(CU_FUNC_CACHE_PREFER_L1);
	//cout << drivererrorcode << endl;

	int device;
	cudaDeviceProp device_prop;

	cudaGetDevice(&device);
	cudaGetDeviceProperties(&device_prop, device);
	
	int threadsPerBlock = device_prop.maxThreadsDim[0];
	//int threadsPerBlock = 32;
	//int numBlocks = device_prop.multiProcessorCount;
	int numBlocks = device_prop.maxGridSize[0];
	
	//Run kernel repeatedly with offset since we cannot launch too many threads at once
	for (int k = 0; k < r * c * z; k+= numBlocks * threadsPerBlock) //Run kernel on groups of pixels at a time
	{
		InitialClusteringKernel_CUDA<<< numBlocks , threadsPerBlock >>>(dev_im_vals, dev_max_response_r, dev_max_response_c, dev_max_response_z , r, c, z, scale_xy, scale_z, k);
	}
	
	//Copy device memory contents back to host memory
	cudaMemcpy(max_response_r, dev_max_response_r, r * c * z * sizeof(*max_response_r), cudaMemcpyDeviceToHost);
	cudaMemcpy(max_response_c, dev_max_response_c, r * c * z * sizeof(*max_response_c), cudaMemcpyDeviceToHost);
	cudaMemcpy(max_response_z, dev_max_response_z, r * c * z * sizeof(*max_response_z), cudaMemcpyDeviceToHost);

	cout << cudaGetErrorString(cudaGetLastError()) << endl;
	
	//Block until all precious commands are complete
	cudaThreadSynchronize();

	cudaFree(dev_im_vals);
	//cudaFree(dev_local_max_vals);
	cudaFree(max_response_r);
	cudaFree(max_response_c);
	cudaFree(max_response_z);
	
	cout << "CUDA done" << endl;
}