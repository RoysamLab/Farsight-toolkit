#include <cuda.h>
#include <iostream>

using namespace std;

__global__ void LocalMaximaKernel_CUDA(float* im_vals, unsigned short* out1, int r, int c, int z, double scale_xy, double scale_z, int offset)
{
	int iGID = blockIdx.x * blockDim.x + threadIdx.x + offset; //global index

	if (iGID >= r * c * z)
		return;
		
	//calculate r, c, z indices as i, j, k from global index
	int rem = ((long)iGID) % (r*c);
	int k = ((int)iGID-rem) / (r*c); 
	int j = ((long)rem) % c;
	int i = (rem-j)/c;
	
	//calculate bounds
	int min_r = (int) max(0.0,i-scale_xy);
	int min_c = (int) max(0.0,j-scale_xy);
	int min_z = (int) max(0.0,k-scale_z);
	int max_r = (int)min((float)r-1,i+scale_xy);
	int max_c = (int)min((float)c-1,j+scale_xy);                         
	int max_z = (int)min((float)z-1,k+scale_z);                         
	
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

extern "C"
void Detect_Local_MaximaPoints_3D_CUDA(float* im_vals, int r, int c, int z, double scale_xy, double scale_z, unsigned short* out1)
{
	cout << "Entering Detect_Local_MaximaPoints_3D_CUDA" << endl;
	
	cudaError_t errorcode;
	float* dev_im_vals; 
	unsigned short* dev_out1;

	//cout << "Allocating " << r * c * z * sizeof(*im_vals) / (double)(1024 * 1024) << " MB of memory on device" << endl;
	//Allocate memory for im_vals and out1
	errorcode = cudaMalloc((void**) &dev_im_vals, r * c * z * sizeof(*im_vals));
	//cout << errorcode << endl;
	errorcode = cudaMalloc((void**) &dev_out1, r * c * z * sizeof(*out1));
	
	//Copy im_vals content into device space
	errorcode = cudaMemcpy(dev_im_vals, im_vals, r * c * z * sizeof(*im_vals), cudaMemcpyHostToDevice);
	//cout << errorcode << endl;
	
	//Prefer 48KB L1 cache
	CUresult drivererrorcode = cuCtxSetCacheConfig(CU_FUNC_CACHE_PREFER_L1);
	//cout << drivererrorcode << endl;

	int device;
	cudaDeviceProp device_prop;

	cudaGetDevice(&device);
	cudaGetDeviceProperties(&device_prop, device);

	/*cout << device_prop.maxGridSize[0] << endl;
	cout << device_prop.maxThreadsDim[0] << endl;*/

	int threadsPerBlock = device_prop.maxThreadsDim[0];
	//int threadsPerBlock = 32;
	int numBlocks = device_prop.multiProcessorCount;
	//int numBlocks = device_prop.maxGridSize[0];

	for (int k = 0; k < r * c * z; k+= numBlocks * threadsPerBlock) //Run kernel on 16K pixels at a time
	{
		LocalMaximaKernel_CUDA<<< numBlocks , threadsPerBlock >>>(dev_im_vals, dev_out1, r, c, z, scale_xy, scale_z, k);
	}
	errorcode = cudaMemcpy(out1, dev_out1, r * c * z * sizeof(*out1), cudaMemcpyDeviceToHost);
	
	//cout << errorcode << endl;

	//Block until all precious commands are complete
	cudaThreadSynchronize();

	cudaFree(dev_im_vals);
	cudaFree(dev_out1);

	cout << cudaGetErrorString(cudaGetLastError()) << endl;

	cout << "CUDA done" << endl;
}