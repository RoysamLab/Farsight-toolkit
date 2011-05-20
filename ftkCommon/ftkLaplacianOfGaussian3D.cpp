#include "ftkLaplacianOfGaussian3D.h"
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>

#ifdef _OPENMP
#include "omp.h"
#endif

#if CUDA
#include "Convolution_CUDA.cuh"
#endif

using namespace std;

//Function prototypes (private)
double*** generateKernel(float scale, int kernel_size, double* error);
double*** zeroPadImage(double*** image, int image_x_size, int image_y_size, int image_z_size, int kernel_size);
double*** convolve(double*** kernel, double*** paddedImage, int padded_image_x_size, int padded_image_y_size, int padded_image_z_size, int kernel_size);
double sumOfProduct(double*** kernel, double*** paddedImage, int padded_image_start_x, int padded_image_start_y, int padded_image_start_z, int kernel_size);
double* flattenKernel(double*** kernel, int kernel_size);
double* flattenPaddedImage(double*** paddedImage, int padded_image_x_size, int padded_image_y_size, int padded_image_z_size);
void freeKernelMem(double*** kernel, int kernel_size);
void freePaddedImageMem(double*** paddedImage, int padded_image_x_size, int padded_image_y_size);

const double PI = atan(1.0) * 4;
	
double*** runLoG(double*** image, float scale, int image_x_size, int image_y_size, int image_z_size)
{
	//left-justify output
	cout << setiosflags(ios::right);
	cout << resetiosflags(ios::left);

	//fixed-point output
	cout << setiosflags(ios::fixed);
	
	//10 decimal point precision
	cout << setprecision(10);
	
	int kernel_size;
	double*** kernel;
	double error = 0;
	
	cout << "Finding appropriate kernel size for scale " << scale << endl;
	//repeatedly generate kernels of sucessively larger sizes until our error falls below 5%
	for (int k = 2 * scale; k < 1000; k++)
	{
		kernel_size = k;
		kernel = generateKernel(scale, kernel_size, &error);
		if (error < 10.0 && error > 0)
			break;

		for (int k = 0; k < kernel_size; k++)
		{
			for (int l = 0; l < kernel_size; l++)
				free(kernel[k][l]);
			free(kernel[k]);
		}
		free(kernel);
	}
	
	cout << "Found appropriate kernel size: " << kernel_size << endl;

	double*** paddedImage = zeroPadImage(image, image_x_size, image_y_size, image_z_size, kernel_size);

	int padded_image_x_size = image_x_size + 2 * (kernel_size / 2); //note that <2 * (kernel_size / 2)> is <2 * padding> which is one less than kernel_size for odd-sized kernel
	int padded_image_y_size = image_y_size + 2 * (kernel_size / 2);
	int padded_image_z_size = image_z_size + 2 * (kernel_size / 2);

	double*** convolvedImage = convolve(kernel, paddedImage, padded_image_x_size , padded_image_y_size, padded_image_z_size, kernel_size);

	freeKernelMem(kernel, kernel_size);
	freePaddedImageMem(paddedImage, padded_image_x_size, padded_image_y_size);

	return convolvedImage;
}

void freeKernelMem(double*** kernel, int kernel_size)
{
	for (int k = 0; k < kernel_size; k++)
	{
		for (int l = 0; l < kernel_size; l++)
			free(kernel[k][l]);
		free(kernel[k]);
	}
	free(kernel);
}

void freePaddedImageMem(double*** paddedImage, int padded_image_x_size, int padded_image_y_size)
{
	for (int k = 0; k < padded_image_x_size; k++)
	{
		for (int l = 0; l < padded_image_y_size; l++)
			free(paddedImage[k][l]);
		free(paddedImage[k]);
	}
	free(paddedImage);
}

double*** generateKernel(float scale, int kernel_size, double* error)
{
	double*** kernel = (double***) malloc(kernel_size * sizeof(double**)); 
	
	for(int k = 0; k < kernel_size; k++)
	{        
		kernel[k] = (double**) malloc(kernel_size * sizeof(double*));
		for(int l = 0; l < kernel_size; l++)
		{			
			kernel[k][l] = (double *) malloc(kernel_size * sizeof(double));
			for(int m = 0; m < kernel_size; m++)
			{				
				kernel[k][l][m] = 0;
			}
		}
	}
	
	for (int k = 0; k < kernel_size; k++)
	{
		for (int l = 0; l < kernel_size; l++)
		{
			for (int m = 0; m < kernel_size; m++)
			{
				//kernel[k][l][m] = 1 / pow(2 * PI * pow(scale, 2), 1.5) * exp( -(pow(k-(kernel_size-1)/2.0, 2) + pow(l-(kernel_size-1)/2.0,2) + pow(m-(kernel_size-1)/2.0, 2))/ (2 * pow(scale, 2))); //gaussian kernel, for testing
				kernel[k][l][m] = 1 
								* (pow(k-(kernel_size-1)/2.0, 2) + pow(l-(kernel_size-1)/2.0, 2) + pow(m-(kernel_size-1)/2.0, 2) - 3 * pow(scale, 2))			// x^2 + y^2 + z^2 - 3 * sigma^2
								/ pow(scale, 2)																													// 1 / sigma^2
								* exp( -(pow(k-(kernel_size-1)/2.0, 2) + pow(l-(kernel_size-1)/2.0,2) + pow(m-(kernel_size-1)/2.0, 2))/ (2 * pow(scale, 2)));	// exponential term of gaussian
			}
		}
	}

	//cout << "Printing out kernel for scale " << scale << endl;
	
	double sumNegValues = 0;
	double sum = 0;
	for (int k = 0; k < kernel_size; k++)
	{
		for (int l = 0; l < kernel_size; l++)
		{
			for (int m = 0; m < kernel_size; m++)
			{
				sum += kernel[m][l][k];
				if (kernel[m][l][k] < 0)
					sumNegValues += kernel[m][l][k];
				//cout << setw(13) << kernel[m][l][k] << " ";
			}
			//cout << endl;
		}
		//cout << endl;
	}
	//cout << endl;

	//cout << "Missing area: " << sum << endl;
	//cout << "Error in % (Missing area / Total negative area): " << sum/sumNegValues * 100 << endl;
	*error = sum/sumNegValues * 100;
	//cout << "Peak response: " << kernel[kernel_size/2][kernel_size/2][kernel_size/2] << endl;

	/*if (sum < kernel[kernel_size/2][kernel_size/2][kernel_size/2])
		cout << "Kernel size not large enough!!" << endl;

	kernel[kernel_size/2][kernel_size/2][kernel_size/2] -= sum;*/

	return kernel;
}



//image needs to be padded so that if statements can be avoided when convolving (branching is bad, doubly so for running on a GPU)
double*** zeroPadImage(double*** image, int image_x_size, int image_y_size, int image_z_size, int kernel_size)
{
	int padding = kernel_size / 2; //Note the truncation
	
	int padded_image_x_size = 2 * padding + image_x_size;
	int padded_image_y_size = 2 * padding + image_y_size;
	int padded_image_z_size = 2 * padding + image_z_size;

	/*cout << "Printing out image" << endl;
	double sum = 0;
	
	for (int k = 0; k < image_z_size; k++)
	{
		for (int l = 0; l < image_y_size; l++)
		{
			for (int m = 0; m < image_x_size; m++)
			{
				cout << image[m][l][k] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}
	cout << endl << endl;*/
	
	
	//Allocate memory for padded image
	double*** paddedImage = (double ***) malloc(padded_image_x_size * sizeof(double**));
	for (int k = 0; k < padded_image_x_size; k++)
	{	
		paddedImage[k] = (double **) malloc(padded_image_y_size * sizeof(double *));
		for (int l = 0; l < padded_image_y_size; l++)
		{
			paddedImage[k][l] = (double *) malloc(padded_image_z_size * sizeof(double));
			for (int m = 0; m < padded_image_z_size; m++)
			{
				paddedImage[k][l][m] = 0;
			}
		}
	}

	//Copy padding and image
	for (int k = 0; k < padded_image_x_size; k++)
		for (int l = 0; l < padded_image_y_size; l++)
			for (int m = 0; m < padded_image_z_size; m++)
			{
				if (k < padding || k >= padding + image_x_size  || l < padding || l >= padding + image_y_size || m < padding || m >= padding + image_z_size)
					paddedImage[k][l][m] = 0;
				else 
					paddedImage[k][l][m] = image[k - padding][l - padding][m - padding];
			}



	/*cout << "Printing out padded image" << endl;
	for (int m = 0; m < padded_image_z_size; m++)
	{
		for (int l = 0; l < padded_image_y_size; l++)
		{
			for (int k = 0; k < padded_image_x_size; k++)
			{
				cout << paddedImage[k][l][m] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}
	cout << endl;

	cout << endl;*/

	return paddedImage;
}

double*** convolve(double*** kernel, double*** paddedImage, int padded_image_x_size, int padded_image_y_size, int padded_image_z_size, int kernel_size)
{
	int padding = kernel_size / 2;

	int outputImage_x_size = padded_image_x_size - 2 * padding;
	int outputImage_y_size = padded_image_y_size - 2 * padding;
	int outputImage_z_size = padded_image_z_size - 2 * padding;

	//allocate memory for output image
	double*** outputImage = (double*** ) malloc(outputImage_x_size * sizeof(double**));
	for(int k = 0; k < outputImage_x_size; k++)
	{        
		outputImage[k] = (double**) malloc(outputImage_y_size * sizeof(double*));
		for(int l = 0; l < outputImage_y_size; l++)
		{			
			outputImage[k][l] = (double *) malloc(outputImage_z_size * sizeof(double));
			for(int m = 0; m < outputImage_z_size; m++)
			{				
				outputImage[k][l][m] = 0;
			}
		}
	}

	//convolution
#if CUDA
	sumOfProduct_CUDA(flattenKernel(kernel, kernel_size), flattenPaddedImage(paddedImage, padded_image_x_size, padded_image_y_size, padded_image_z_size), outputImage_x_size, outputImage_y_size, outputImage_z_size, kernel_size, outputImage);
#else
	#pragma omp parallel for
	for(int k = 0; k < outputImage_x_size; k++)
	{        
		for(int l = 0; l < outputImage_y_size; l++)
		{			
			for(int m = 0; m < outputImage_z_size; m++)
			{				
				outputImage[k][l][m] = sumOfProduct(kernel, paddedImage, k, l, m, kernel_size);
			}
		}
	}
#endif

	/*cout << "Printing out convolved image" << endl;
	for (int m = 0; m < outputImage_z_size; m++)
	{
		for (int l = 0; l < outputImage_y_size; l++)
		{
			for (int k = 0; k < outputImage_x_size; k++)
			{
				cout << outputImage[k][l][m] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}
	cout << endl;

	cout << endl;*/

	return outputImage;
}

double sumOfProduct(double*** kernel, double*** paddedImage, int padded_image_start_x, int padded_image_start_y, int padded_image_start_z, int kernel_size)
{
	int padding = kernel_size / 2;
	double sum = 0;

	for (int k = 0; k < kernel_size; k++)
		for (int l = 0; l < kernel_size; l++)
			for (int m = 0; m < kernel_size; m++)
				sum += kernel[k][l][m] * paddedImage[padded_image_start_x + k][padded_image_start_y + l][padded_image_start_z + m]; //note that this is reduction

	return sum;
}

double* flattenKernel(double*** kernel, int kernel_size)
{
	double* flat_kernel = (double*) malloc(kernel_size * kernel_size * kernel_size * sizeof(double));

	for (int k = 0; k < kernel_size; k++)
		for (int l = 0; l < kernel_size; l++)
			for (int m = 0; m < kernel_size; m++)
				flat_kernel[m + kernel_size * l + kernel_size * kernel_size * k] = kernel[k][l][m];
	
	/*for (int n = 0; n < kernel_size * kernel_size * kernel_size; n++)
	{
		int k =	n / (kernel_size * kernel_size);
		int l = n % (kernel_size * kernel_size) / kernel_size;
		int m = n % kernel_size;
		 
		cout << k << " " << l << " " << m << endl;
	}*/

	return flat_kernel;
}

double* flattenPaddedImage(double*** paddedImage, int padded_image_x_size, int padded_image_y_size, int padded_image_z_size)
{
	double* flat_padded_image = (double*) malloc(padded_image_x_size * padded_image_y_size * padded_image_z_size * sizeof(double));

	for (int k = 0; k < padded_image_x_size; k++)
		for (int l = 0; l < padded_image_y_size; l++)
			for (int m = 0; m < padded_image_z_size; m++)
				flat_padded_image[m + padded_image_z_size * l + padded_image_z_size * padded_image_y_size * k] = paddedImage[k][l][m];

	/*for (int n = 0; n < padded_image_x_size * padded_image_y_size * padded_image_z_size; n++)
	{
		int k =	n / (padded_image_z_size * padded_image_y_size);
		int l = n % (padded_image_z_size * padded_image_y_size) / padded_image_z_size;
		int m = n % padded_image_z_size;
		 
		cout << k << " " << l << " " << m << endl;
	}*/


	return flat_padded_image;
}
