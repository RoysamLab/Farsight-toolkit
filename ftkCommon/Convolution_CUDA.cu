#include <cuda.h>
#include <iostream>

using namespace st;

__global__ void ConvolutionKernel_CUDA (double* kernel, double* paddedImage, double* outputImage, int outputImage_x_size, int outputImage_y_size, int outputImage_z_size, int padded_image_x_size, int padded_image_y_size, int padded_image_z_size, int kernel_size, int offset)
{
	int outImageIndex = blockIdx.x * blockDim.x + threadIdx.x + offset; //global index

	if (outImageIndex >= outputImage_x_size * outputImage_y_size * outputImage_z_size )
		return;

	double sum = 0;

	int imageIndex_x = outImageIndex / (outputImage_z_size * outputImage_y_size);
	int imageIndex_y = outImageIndex % (outputImage_z_size * outputImage_y_size) / padded_image_z_size;
	int imageIndex_z = outImageIndex % outputImage_z_size;

	for (int k = 0; k < kernel_size; k++)
	{
		for (int l = 0; l < kernel_size; l++)
		{
			for (int m = 0; m < kernel_size; m++)
			{
				sum += kernel[m + l * kernel_size + k * kernel_size * kernel_size] * paddedImage[(imageIndex_z + m) + ((imageIndex_y + l) * padded_image_z_size) + ((imageIndex_x + k)* padded_image_z_size * padded_image_y_size)];
			}
		}
	}

	outputImage[outImageIndex] = sum;
}

double*** Convolution_CUDA(double* kernel, double* paddedImage, int padded_image_x_size, int padded_image_y_size, padded_image_z_size, int kernel_size)
{	
	cout << "Entering sumOfProduct_CUDA" << endl;

	double* outputImage = (double*) malloc(outputImage_x_size * outputImage_y_size * outputImage_z_size * sizeof(double));

	int padding = kernel_size / 2;
	int padded_image_x_size = 2 * padding + outputImage_x_size;
	int padded_image_y_size = 2 * padding + outputImage_y_size;
	int padded_image_z_size = 2 * padding + outputImage_z_size;


	cudaError_t errorcode;

	double* dev_kernel; 
	double* dev_paddedImage;
	double* dev_outputImage;

	size_t free_mem, total_mem;
	cudaMemGetInfo(&free_mem, &total_mem);

	cout << free_mem / (double)(1024 * 1024) << " " << total_mem / (double)(1024 * 1024) << endl;

	cout << "Allocating " << (sizeof(*kernel) * kernel_size * kernel_size * kernel_size)/(double)(1024*1024) << " MB of memory on GPU for kernel" << endl;
	cout << "Allocating " << (sizeof(*paddedImage) * padded_image_x_size * padded_image_y_size * padded_image_z_size)/(double)(1024*1024) << " MB of memory on GPU for paddedImage" << endl;
	cout << "Allocating " << (sizeof(*outputImage) * outputImage_x_size * outputImage_y_size * outputImage_z_size)/(double)(1024*1024) << " MB of memory on GPU for outputImage" << endl;

	//allocate memory on device
	errorcode = cudaMalloc((void**) &dev_kernel, kernel_size * kernel_size * kernel_size * sizeof(*dev_kernel));
	errorcode = cudaMalloc((void**) &dev_paddedImage, padded_image_x_size * padded_image_y_size * padded_image_z_size * sizeof(*dev_paddedImage));
	errorcode = cudaMalloc((void**) &dev_outputImage, outputImage_x_size * outputImage_y_size * outputImage_z_size * sizeof(*dev_outputImage));
	
	//cout << errorcode << endl;

	//Copy host memory contents to device contents
	cudaMemcpy(dev_kernel, kernel, kernel_size * kernel_size * kernel_size * sizeof(*kernel), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_paddedImage, paddedImage, padded_image_x_size * padded_image_y_size * padded_image_z_size * sizeof(*paddedImage), cudaMemcpyHostToDevice);

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
	for (int k = 0; k < outputImage_x_size * outputImage_y_size * outputImage_z_size; k+= numBlocks * threadsPerBlock) //Run kernel on groups of pixels at a time
	{
		ConvolutionKernel_CUDA<<< numBlocks , threadsPerBlock >>>(dev_kernel, dev_paddedImage, dev_outputImage, outputImage_x_size, outputImage_y_size, outputImage_z_size, padded_image_x_size, padded_image_y_size, padded_image_z_size, kernel_size, k);
	}
	
	//Copy device memory contents back to host memory
	cudaMemcpy(outputImage, dev_outputImage, outputImage_x_size * outputImage_y_size * outputImage_z_size * sizeof(*outputImage), cudaMemcpyDeviceToHost);

	cout << cudaGetErrorString(cudaGetLastError()) << endl;
	
	//Block until all precious commands are complete
	cudaThreadSynchronize();

	cudaFree(dev_kernel);
	cudaFree(dev_paddedImage);
	cudaFree(dev_outputImage);

	//unflatten outputImage
	for (int n = 0; n < outputImage_x_size * outputImage_y_size * outputImage_z_size; n++)
	{
		int k =	n / (outputImage_z_size * outputImage_y_size);
		int l = n % (outputImage_z_size * outputImage_y_size) / outputImage_z_size;
		int m = n % outputImage_z_size;

		output3DImage[k][l][m] = outputImage[n];
	}

	//Testing by making output equal to the input
	/*for (int n = 0; n < outputImage_x_size * outputImage_y_size * outputImage_z_size; n++)
	{
		int k =	n / (padded_image_z_size * padded_image_y_size);
		int l = n % (padded_image_z_size * padded_image_y_size) / padded_image_z_size;
		int m = n % padded_image_z_size;
		 
		if (k < padding || l < padding || m < padding || k >= outputImage_x_size || l >= outputImage_y_size || m >= outputImage_z_size)
			continue;
		else
			output3DImage[k - padding][l - padding][m - padding] = paddedImage[n];
	}*/

	free(outputImage);
	
	cout << "CUDA_Convolution done" << endl;
}

double* flattenImage(double*** image, int image_x_size, int image_y_size, int image_z_size)
{
	double* flat_image = (double*) malloc(image_x_size * image_y_size * image_z_size * sizeof(double));

	for (int k = 0; k < image_x_size; k++)
		for (int l = 0; l < image_y_size; l++)
			for (int m = 0; m < image_z_size; m++)
				flattened_image[m + image_z_size * l + image_z_size * image_y_size * k] = paddedImage[k][l][m];

	/*for (int n = 0; n < image_x_size * image_y_size * image_z_size; n++)
	{
		int k =	n / (image_z_size * image_y_size);
		int l = n % (image_z_size * image_y_size) / image_z_size;
		int m = n % image_z_size;
		 
		cout << k << " " << l << " " << m << endl;
	}*/


	return flattened_image;
}

