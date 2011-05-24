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
double* generateGaussianKernel(float scale, int kernel_size);
double*** padImage(double*** image, int image_x_size, int image_y_size, int image_z_size, int kernel_size_X, int kernel_size_Y, int kernel_size_Z);
double*** convolveGaussian(double* kernel_X, double* kernel_Y, double* kernel_Z, double*** paddedImage, int padded_image_x_size, int padded_image_y_size, int padded_image_z_size, int kernel_size_X, int kernel_size_Y, int kernel_size_Z);
double sumOfProductX(double* kernel, double*** paddedImage, int padded_image_start_x, int padded_image_start_y, int padded_image_start_z, int kernel_size);
double sumOfProductY(double* kernel, double*** paddedImage, int padded_image_start_x, int padded_image_start_y, int padded_image_start_z, int kernel_size);
double sumOfProductZ(double* kernel, double*** paddedImage, int padded_image_start_x, int padded_image_start_y, int padded_image_start_z, int kernel_size);
double*** generateLaplacianKernel();
double*** convolveLaplacian(double*** kernel, double*** image, int image_x_size, int image_y_size, int image_z_size);
double sumOfProduct(double*** kernel, double*** paddedImage, int padded_image_start_x, int padded_image_start_y, int padded_image_start_z, int kernel_size);
double*** unpadImage(double*** padded_image, int image_x_size, int image_y_size, int image_z_size, int padding_X, int padding_Y, int padding_Z);
void freeImageMem(double*** paddedImage, int padded_image_x_size, int padded_image_y_size);

const double PI = atan(1.0) * 4;
	
double*** runLoG(double*** image, float scale_X, float scale_Y, float scale_Z, int image_x_size, int image_y_size, int image_z_size)
{
	cout << "Running LoG for scale X: " << scale_X << " Y: " << scale_Y << " Z: " << scale_Z << endl;
	//left-justify output
	cout << setiosflags(ios::right);
	cout << resetiosflags(ios::left);

	//fixed-point output
	cout << setiosflags(ios::fixed);
	
	//10 decimal point precision
	cout << setprecision(10);
	
	int kernel_size_X = 7 * 2 * scale_X;
	int kernel_size_Y = 7 * 2 * scale_Y;
	int kernel_size_Z = 7 * 2 * scale_Z;

	int padding_X = kernel_size_X / 2 + 1;
	int padding_Y = kernel_size_X / 2 + 1;
	int padding_Z = kernel_size_X / 2 + 1;
	
	double* kernel_X = generateGaussianKernel(scale_X, kernel_size_X);
	double* kernel_Y = generateGaussianKernel(scale_Y, kernel_size_Y);
	double* kernel_Z = generateGaussianKernel(scale_Z, kernel_size_Z);

	double*** paddedImage = padImage(image, image_x_size, image_y_size, image_z_size, kernel_size_X + 2, kernel_size_Y + 2, kernel_size_Z + 2);

	int padded_image_x_size = 2 * padding_X + image_x_size;
	int padded_image_y_size = 2 * padding_Y + image_y_size;
	int padded_image_z_size = 2 * padding_Z + image_z_size;

	double*** smoothedPaddedImage = convolveGaussian(kernel_X, kernel_Y, kernel_Z, paddedImage, padded_image_x_size, padded_image_y_size, padded_image_z_size, kernel_size_X, kernel_size_Y, kernel_size_Z);
	freeImageMem(paddedImage, padded_image_x_size, padded_image_y_size);

	free(kernel_X);
	free(kernel_Y);
	free(kernel_Z);
	
	double*** laplacianKernel = generateLaplacianKernel();

	double*** paddedPaddedSmoothedImage = padImage(smoothedPaddedImage, padded_image_x_size, padded_image_y_size, padded_image_z_size, 3, 3, 3);
	freeImageMem(smoothedPaddedImage, padded_image_x_size, padded_image_y_size);

	double*** paddedLoGimage = convolveLaplacian(laplacianKernel, paddedPaddedSmoothedImage, padded_image_x_size + 2, padded_image_y_size + 2, padded_image_z_size + 2);
	freeImageMem(paddedPaddedSmoothedImage, padded_image_x_size + 2, padded_image_y_size + 2);
	
	double*** LoGimage = unpadImage(paddedLoGimage, image_x_size, image_y_size, image_z_size, padding_X + 1, padding_Y + 1, padding_Z + 1);
	freeImageMem(paddedLoGimage, padded_image_x_size + 2, padded_image_y_size + 2);

  
	//For testing Laplacian	
	/*double*** laplacianKernel = generateLaplacianKernel();

	double*** paddedImage = padImage(image, image_x_size, image_y_size, image_z_size, 3, 3, 3);
	double*** paddedLoGimage = convolveLaplacian(laplacianKernel, paddedImage, image_x_size + 2, image_y_size + 2, image_z_size + 2);
	double*** LoGimage = unpadImage(paddedLoGimage, image_x_size, image_y_size, image_z_size, 1, 1, 1);*/


	return LoGimage;
}

void freeImageMem(double*** paddedImage, int padded_image_x_size, int padded_image_y_size)
{
	for (int k = 0; k < padded_image_x_size; k++)
	{
		for (int l = 0; l < padded_image_y_size; l++)
			free(paddedImage[k][l]);
		free(paddedImage[k]);
	}
	free(paddedImage);
}

double* generateGaussianKernel(float scale, int kernel_size)
{	
	double* kernel = (double*) malloc(kernel_size * sizeof(double**)); 
	
	for (int k = 0; k < kernel_size; k++)
	{
		kernel[k] = 1
					* pow(scale, 2) 
					* 1 / sqrt(2 * PI * pow(scale, 2)) 
					* exp( -(pow(k-(kernel_size-1)/2.0, 2))/ (2 * pow(scale, 2))); //1-D gaussian kernel
	}
	
	/*cout << "Printing out kernel for scale: " << scale << endl;
	for (int k = 0; k < kernel_size; k++)
	{
		cout << kernel[k] << " ";
	}
	cout << endl;*/
	
	return kernel;
}

double*** generateLaplacianKernel()
{
	double*** kernel = (double ***) malloc(3 * sizeof(double **));
	for (int k = 0; k < 3; k++)
	{
		kernel[k] = (double **) malloc(3 * sizeof(double *));
		for (int l = 0; l < 3; l++)
		{
			kernel[k][l] = (double *) malloc(3 * sizeof(double));
			for (int m = 0; m < 3; m++)
			{
				kernel[k][l][m] = 0;
			}
		}
	}
	
	/* 3D Laplacian filter, see Wikipedia article on Discrete Laplace Operator */
	kernel[1][1][0] = 1;
	
	kernel[1][0][1] = 1;
	kernel[0][1][1] = 1;
	kernel[1][1][1] = -6;
	kernel[2][1][1] = 1;
	kernel[1][2][1] = 1;
	
	kernel[1][1][2] = 1;

	return kernel;
}



//image needs to be padded (pixels value clamped) so that if statements can be avoided when convolving (branching is bad, doubly so for running on a GPU)
double*** padImage(double*** image, int image_x_size, int image_y_size, int image_z_size, int kernel_size_X, int kernel_size_Y, int kernel_size_Z)
{	
	int padding_X = kernel_size_X / 2;
	int padding_Y = kernel_size_Y / 2;
	int padding_Z = kernel_size_Z / 2;
	
	int padded_image_x_size = 2 * padding_X + image_x_size;
	int padded_image_y_size = 2 * padding_Y + image_y_size;
	int padded_image_z_size = 2 * padding_Z + image_z_size;

	cout << "Padded image size: " << padded_image_x_size << "x" << padded_image_y_size << "x" << padded_image_z_size << endl;

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
	{
		for (int l = 0; l < padded_image_y_size; l++)
		{
			for (int m = 0; m < padded_image_z_size; m++)
			{
				if (k < padding_X || k >= padding_X + image_x_size  || l < padding_Y || l >= padding_Y + image_y_size || m < padding_Z || m >= padding_Z + image_z_size)
				{	
					//Clamp pixels that are edges but not "corners". Corners are defined as places where neither the seperable Gaussian kernels nor the Laplacian kernels would ever read data from when convolving. (6 cases, 2 for each direction)
					/*if (!(l < padding_Y || l >= padding_Y + image_y_size) && !(m < padding_Z || m >= padding_Z + image_z_size))			//outside X but within Y and Z
					{
						if (k < padding_X)
							paddedImage[k][l][m] = image[0][l - padding_Y][m - padding_Z];
						else if (k >= padding_X + image_x_size)
							paddedImage[k][l][m] = image[image_x_size - 1][l - padding_Y][m - padding_Z];
					}
					else if (!(k < padding_X || k >= padding_X + image_x_size) && !(m < padding_Z || m >= padding_Z + image_z_size))	//outside Y but within X and Z
					{
						if (l < padding_Y)
							paddedImage[k][l][m] = image[k - padding_X][0][m - padding_Z];
						else if (l >= padding_Y + image_z_size)
							paddedImage[k][l][m] = image[k - padding_X][image_y_size - 1][m - padding_Z];
					}
					else if (!(k < padding_X || k >= padding_X + image_x_size) && !(l < padding_Y || l >= padding_Y + image_y_size))	 //outside Z but within X and Y
					{
						if (m < padding_Z)
							paddedImage[k][l][m] = image[k - padding_X][l - padding_Y][0];
						else if (m >= padding_Z + image_z_size)
							paddedImage[k][l][m] = image[k - padding_X][l - padding_Y][image_z_size - 1];
					}
					else*/
						paddedImage[k][l][m] = 0;
				}
				else 
					paddedImage[k][l][m] = image[k - padding_X][l - padding_Y][m - padding_Z];
			}
		}
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

double*** convolveGaussian(double* kernel_X, double* kernel_Y, double* kernel_Z, double*** paddedImage, int padded_image_x_size, int padded_image_y_size, int padded_image_z_size, int kernel_size_X, int kernel_size_Y, int kernel_size_Z)
{
	cout << "Convolving Gaussian" << endl;

	int padding_X = kernel_size_X / 2 + 1;
	int padding_Y = kernel_size_Y / 2 + 1;
	int padding_Z = kernel_size_Z / 2 + 1;

	int image_x_size = padded_image_x_size - 2 * padding_X; //note that <2 * (kernel_size / 2)> is <2 * padding> which is one less than kernel_size for odd-sized kernel
	int image_y_size = padded_image_y_size - 2 * padding_Y;
	int image_z_size = padded_image_z_size - 2 * padding_Z;

	//Allocate memory for temporary padded image (to hold result after each direction of convolution)
	double*** tempPaddedImage = (double ***) malloc(padded_image_x_size * sizeof(double**));
	for (int k = 0; k < padded_image_x_size; k++)
	{	
		tempPaddedImage[k] = (double **) malloc(padded_image_y_size * sizeof(double *));
		for (int l = 0; l < padded_image_y_size; l++)
		{
			tempPaddedImage[k][l] = (double *) malloc(padded_image_z_size * sizeof(double));
			for (int m = 0; m < padded_image_z_size; m++)
			{
				tempPaddedImage[k][l][m] = paddedImage[k][l][m];
			}
		}
	}

	//convolution

	//convolving in x direction
	#pragma omp parallel for
	for(int k = 0; k < image_x_size + 2; k++)
	{        
		for(int l = 0; l < image_y_size + 2; l++)
		{			
			for(int m = 0; m < image_z_size + 2; m++)
			{				
				tempPaddedImage[k + padding_X - 1][l + padding_Y - 1][m + padding_Z - 1] = sumOfProductX(kernel_X, paddedImage, k + padding_X - kernel_size_X / 2, l + padding_Y - 1, m + padding_Z - 1, kernel_size_X);
			}
		}
	}
	

	//convolving in y-direction
	#pragma omp parallel for
	for(int k = 0; k < image_x_size + 2; k++)
	{        
		for(int l = 0; l < image_y_size + 2; l++)
		{			
			for(int m = 0; m < image_z_size + 2; m++)
			{				
				paddedImage[k + padding_X - 1][l + padding_Y - 1][m + padding_Z - 1] = sumOfProductY(kernel_Y, tempPaddedImage, k + padding_X - 1, l + padding_Y - kernel_size_Y / 2, m + padding_Z - 1, kernel_size_Y);
			}
		}
	}

	//convolving in z-direction
	#pragma omp parallel for
	for(int k = 0; k < image_x_size + 2; k++)
	{        
		for(int l = 0; l < image_y_size + 2; l++)
		{			
			for(int m = 0; m < image_z_size + 2; m++)
			{				
				tempPaddedImage[k + padding_X - 1][l + padding_Y - 1][m + padding_Z - 1] = sumOfProductZ(kernel_Z, paddedImage, k + padding_X - 1, l + padding_Y - 1, m + padding_Z - kernel_size_Z / 2, kernel_size_Z);
			}
		}
	}

	/*
	//Testing, feed input into output
	#pragma omp parallel for
	for(int k = 0; k < image_x_size; k++)
	{        
		for(int l = 0; l < image_y_size; l++)
		{			
			for(int m = 0; m < image_z_size; m++)
			{				
				image[k][l][m] = paddedImage[k+padding_X][l+padding_Y][m+padding_Z];
			}
		}
	}*/

	/*cout << "Printing out convolved image" << endl;
	for (int m = 0; m < image_z_size; m++)
	{
		for (int l = 0; l < image_y_size; l++)
		{
			for (int k = 0; k < image_x_size; k++)
			{
				cout << image[k][l][m] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}
	cout << endl;

	cout << endl;*/

	return tempPaddedImage;
}

double sumOfProductX(double* kernel, double*** paddedImage, int padded_image_start_x, int padded_image_start_y, int padded_image_start_z, int kernel_size)
{
	double sum = 0;

	for (int k = 0; k < kernel_size; k++)
		sum += kernel[k] * paddedImage[padded_image_start_x + k][padded_image_start_y][padded_image_start_z]; //note that this is reduction

	return sum;
}

double sumOfProductY(double* kernel, double*** paddedImage, int padded_image_start_x, int padded_image_start_y, int padded_image_start_z, int kernel_size)
{
	double sum = 0;

	for (int k = 0; k < kernel_size; k++)
		sum += kernel[k] * paddedImage[padded_image_start_x][padded_image_start_y + k][padded_image_start_z]; //note that this is reduction

	return sum;
}

double sumOfProductZ(double* kernel, double*** paddedImage, int padded_image_start_x, int padded_image_start_y, int padded_image_start_z, int kernel_size)
{
	double sum = 0;

	for (int k = 0; k < kernel_size; k++)
		sum += kernel[k] * paddedImage[padded_image_start_x][padded_image_start_y][padded_image_start_z + k]; //note that this is reduction

	return sum;
}


double*** convolveLaplacian(double*** kernel, double*** padded_image, int padded_image_x_size, int padded_image_y_size, int padded_image_z_size)
{
	int image_x_size = padded_image_x_size - 2; 
	int image_y_size = padded_image_y_size - 2;
	int image_z_size = padded_image_z_size - 2;
	
	double*** outputImage = (double ***) malloc(padded_image_x_size * sizeof(double**));
	for (int k = 0; k < padded_image_x_size; k++)
	{	
		outputImage[k] = (double **) malloc(padded_image_y_size * sizeof(double *));
		for (int l = 0; l < padded_image_y_size; l++)
		{
			outputImage[k][l] = (double *) malloc(padded_image_z_size * sizeof(double));
			for (int m = 0; m < padded_image_z_size; m++)
			{
				outputImage[k][l][m] = 0;
			}
		}
	}
	
	cout << "Convolving Laplacian" << endl;

	//cout << "Convolving Laplacian over: " << image_x_size << "x" << image_y_size << "x" << image_z_size << endl;
	
	#pragma omp parallel for
	for(int k = 0; k < image_x_size; k++)
	{        
		for(int l = 0; l < image_y_size; l++)
		{			
			for(int m = 0; m < image_z_size; m++)
			{				
				outputImage[k + 1][l + 1][m + 1] = sumOfProduct(kernel, padded_image, k, l, m, 3);
			}
		}
	}

	return outputImage;
}

double sumOfProduct(double*** kernel, double*** paddedImage, int padded_image_start_x, int padded_image_start_y, int padded_image_start_z, int kernel_size)
{
	double sum = 0;

	for (int k = 0; k < kernel_size; k++)
		for (int l = 0; l < kernel_size; l++)
			for (int m = 0; m < kernel_size; m++)
				sum += kernel[k][l][m] * paddedImage[padded_image_start_x + k][padded_image_start_y + l][padded_image_start_z + m]; //note that this is reduction

	return sum;
}

double*** unpadImage(double*** padded_image, int image_x_size, int image_y_size, int image_z_size, int padding_X, int padding_Y, int padding_Z)
{
	double*** outputImage = (double ***) malloc(image_x_size * sizeof(double**));
	for (int k = 0; k < image_x_size; k++)
	{	
		outputImage[k] = (double **) malloc(image_y_size * sizeof(double *));
		for (int l = 0; l < image_y_size; l++)
		{
			outputImage[k][l] = (double *) malloc(image_z_size * sizeof(double));
			for (int m = 0; m < image_z_size; m++)
			{
				outputImage[k][l][m] = 0;
			}
		}
	}

	for (int k = 0; k < image_x_size; k++)
		for (int l = 0; l < image_y_size; l++)
			for (int m = 0; m < image_z_size; m++)
				outputImage[k][l][m] = padded_image[k + padding_X][l + padding_Y][m + padding_Z];
	
	return outputImage;
}