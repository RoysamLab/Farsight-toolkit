#include "ftkLaplacianOfGaussian3D.h"

#include <limits>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>


#ifdef _OPENMP
#include "omp.h"
#endif

#define LoGSmoothnessPadding 0

const double PI = atan(1.0) * 4;

//The constructor sets up the image sizes and output formatting and generates the kernels required for each step
template <typename TPixelType>
ftkLaplacianOfGaussian3D<TPixelType>::ftkLaplacianOfGaussian3D(TPixelType*** image, float scale_X, float scale_Y, float scale_Z, unsigned int image_x_size, unsigned int image_y_size, unsigned int image_z_size)
{
	this->image = image;
	
	this->scale_X = scale_X;
	this->scale_Y = scale_Y;
	this->scale_Z = scale_Z;
	
	this->image_x_size = image_x_size;
	this->image_y_size = image_y_size;
	this->image_z_size = image_z_size;

	//setup std::cout so that is left-justified
	std::cout << std::setiosflags(std::ios::right);
	std::cout << std::resetiosflags(std::ios::left);

	//fixed-point output
	std::cout << std::setprecision(std::ios::fixed);

	//10 decimal point precision
	std::cout << std::setprecision(10);

	//Generate the kernel dimensions and the actual kernels themselves
	gaussian_kernel_size_x = 3 * 2 * scale_X;
	gaussian_kernel_size_y = 3 * 2 * scale_Y;
	gaussian_kernel_size_z = 3 * 2 * scale_Z;

	kernel_X = generateGaussianKernel(scale_X, gaussian_kernel_size_x);
	kernel_Y = generateGaussianKernel(scale_Y, gaussian_kernel_size_y);
	kernel_Z = generateGaussianKernel(scale_Z, gaussian_kernel_size_z);

	laplacian_kernel_size_x = 3;
	laplacian_kernel_size_y = 3;
	laplacian_kernel_size_z = 3;

	laplacian_kernel = generateLaplacianKernel();

	LoGImage = NULL;
}

//The main function where all of the actual LoG calculation filter starts
template <typename TPixelType>
void ftkLaplacianOfGaussian3D<TPixelType>::RunFilter()
{
	std::cout << "Running LoG filter at scale X: " << scale_X << " Y: " << scale_Y << " Z: " << scale_Z << std::endl;

	//pad image for Gaussian
	TPixelType*** paddedImage = padImage(image, image_x_size, image_y_size, image_z_size, gaussian_kernel_size_x / 2 + LoGSmoothnessPadding, gaussian_kernel_size_y / 2 + LoGSmoothnessPadding, gaussian_kernel_size_z / 2 + LoGSmoothnessPadding);
	
	unsigned int paddedGaussianImage_x_size = image_x_size + gaussian_kernel_size_x + LoGSmoothnessPadding * 2;
	unsigned int paddedGaussianImage_y_size = image_y_size + gaussian_kernel_size_y + LoGSmoothnessPadding * 2;
	unsigned int paddedGaussianImage_z_size = image_z_size + gaussian_kernel_size_z + LoGSmoothnessPadding * 2;

	//run Gaussian over padded image
	TPixelType*** paddedGaussianImage = convolveGaussian(kernel_X, kernel_Y, kernel_Z, paddedImage, paddedGaussianImage_x_size, paddedGaussianImage_y_size, paddedGaussianImage_z_size, gaussian_kernel_size_x, gaussian_kernel_size_y, gaussian_kernel_size_z);
	freeImageMem(paddedImage, paddedGaussianImage_x_size, paddedGaussianImage_y_size);

	//pad image again for Laplacian
	TPixelType*** paddedPaddedGaussianImage = padImage(paddedGaussianImage, paddedGaussianImage_x_size, paddedGaussianImage_y_size, paddedGaussianImage_z_size, laplacian_kernel_size_x / 2, laplacian_kernel_size_y / 2, laplacian_kernel_size_z / 2);
	freeImageMem(paddedGaussianImage, paddedGaussianImage_x_size, paddedGaussianImage_y_size);

	unsigned int paddedpaddedGaussianImage_x_size = paddedGaussianImage_x_size + 2 * (laplacian_kernel_size_x / 2);
	unsigned int paddedpaddedGaussianImage_y_size = paddedGaussianImage_y_size + 2 * (laplacian_kernel_size_y / 2);
	unsigned int paddedpaddedGaussianImage_z_size = paddedGaussianImage_z_size + 2 * (laplacian_kernel_size_z / 2);

	//run Laplacian over padded image
	TPixelType*** paddedPaddedLoGImage = convolveLaplacian(laplacian_kernel, paddedPaddedGaussianImage, paddedpaddedGaussianImage_x_size, paddedpaddedGaussianImage_y_size, paddedpaddedGaussianImage_z_size);
	freeImageMem(paddedPaddedGaussianImage, paddedpaddedGaussianImage_x_size, paddedpaddedGaussianImage_y_size);

	//unpad the previous 2 paddings
	LoGImage = unpadImage(paddedPaddedLoGImage, image_x_size, image_y_size, image_z_size, gaussian_kernel_size_x / 2 + LoGSmoothnessPadding/2 + laplacian_kernel_size_x / 2, gaussian_kernel_size_y / 2 + LoGSmoothnessPadding/2 + laplacian_kernel_size_y / 2, gaussian_kernel_size_z / 2 + LoGSmoothnessPadding/2 + laplacian_kernel_size_z / 2);
	freeImageMem(paddedPaddedLoGImage, paddedpaddedGaussianImage_x_size, paddedpaddedGaussianImage_x_size);

	std::cout << "Laplacian of Gaussian complete" << std::endl;
}

//Accessor to get the LoG output
template <typename TPixelType>
TPixelType*** ftkLaplacianOfGaussian3D<TPixelType>::GetOutput()
{
	/*std::cout << "Printing out output image" << std::endl;
	for (int m = 0; m < image_z_size; m++)
	{
		for (int l = 0; l < image_y_size; l++)
		{
			for (int k = 0; k < image_x_size; k++)
			{
				std::cout << LoGImage[k][l][m] << " ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	std::cout << std::endl;*/
	
	return LoGImage;
}

//function to free all the image memory (works on Laplacian kernel too)
template <typename TPixelType>
void ftkLaplacianOfGaussian3D<TPixelType>::freeImageMem(TPixelType*** image, unsigned int image_x_size, unsigned int image_y_size)
{
	for (int k = 0; k < image_x_size; k++)
	{
		for (int l = 0; l < image_y_size; l++)
			free(image[k][l]);
		free(image[k]);
	}
	free(image);

	image = NULL;
}

//Generates the Laplacian kernel
template <typename TPixelType>
double* ftkLaplacianOfGaussian3D<TPixelType>::generateGaussianKernel(float scale, unsigned int kernel_size)
{	
	//Allocate memory for Gaussian Kernel
	double* kernel = (double*) malloc(kernel_size * sizeof(double**)); 
	
	for (int k = 0; k < kernel_size; k++)
	{
		kernel[k] = 1
					//* pow(scale, 2) //scale normalization
					* 1 / sqrt(2 * PI * pow(scale, 2)) 
					* exp( -(pow(k-(kernel_size-1)/2.0, 2)) / (2 * pow(scale, 2))); //1-D gaussian kernel
	}
	
	/*std::cout << "Printing out kernel for scale: " << scale << std::endl;
	for (int k = 0; k < kernel_size; k++)
	{
		std::cout << kernel[k] << " ";
	}
	std::cout << std::endl;*/
	
	return kernel;
}

//Generates Laplacian Kernel
template <typename TPixelType>
double*** ftkLaplacianOfGaussian3D<TPixelType>::generateLaplacianKernel()
{
	//Allocate memory for Laplacian Kernel	
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
	kernel[1][1][0] = 1/6.0;
	
	kernel[1][0][1] = 1/6.0;
	kernel[0][1][1] = 1/6.0;
	kernel[1][1][1] = -1;
	kernel[2][1][1] = 1/6.0;
	kernel[1][2][1] = 1/6.0;
	
	kernel[1][1][2] = 1/6.0;

	return kernel;
}

//image needs to be padded (and optionally pixels value clamped) so that if statements can be avoided when convolving (branching is bad, doubly so for running on a GPU)
template <typename TPixelType>
TPixelType*** ftkLaplacianOfGaussian3D<TPixelType>::padImage(TPixelType*** image, unsigned int image_x_size, unsigned int image_y_size, unsigned int image_z_size, unsigned int padding_X, unsigned int padding_Y, unsigned int padding_Z)
{	
	int padded_image_x_size = 2 * padding_X + image_x_size;
	int padded_image_y_size = 2 * padding_Y + image_y_size;
	int padded_image_z_size = 2 * padding_Z + image_z_size;
	
	std::cout << "Original image size: " << image_x_size << "x" << image_y_size << "x" << image_z_size << std::endl;
	std::cout << "Padded image size: " << padded_image_x_size << "x" << padded_image_y_size << "x" << padded_image_z_size << std::endl;

	//std::cout << "Printing out image" << std::endl;
	//
	//for (int k = 0; k < image_z_size; k++)
	//{
	//	for (int l = 0; l < image_y_size; l++)
	//	{
	//		for (int m = 0; m < image_x_size; m++)
	//		{
	//			std::cout << image[m][l][k] << " ";
	//		}
	//		std::cout << std::endl;
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout << std::endl << std::endl;
	
	
	//Allocate memory for padded image
	TPixelType*** paddedImage = (TPixelType ***) malloc(padded_image_x_size * sizeof(TPixelType**));
	for (int k = 0; k < padded_image_x_size; k++)
	{	
		paddedImage[k] = (TPixelType **) malloc(padded_image_y_size * sizeof(TPixelType *)); 
		for (int l = 0; l < padded_image_y_size; l++)
		{
			paddedImage[k][l] = (TPixelType *) malloc(padded_image_z_size * sizeof(TPixelType));
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


	std::cout << "Printing out padded image" << std::endl;
	for (int m = 0; m < padded_image_z_size; m++)
	{
		for (int l = 0; l < padded_image_y_size; l++)
		{
			for (int k = 0; k < padded_image_x_size; k++)
			{
				std::cout << paddedImage[k][l][m] << " ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	std::cout << std::endl;
	
	return paddedImage;
}

//This function is where the convolution of the Gaussian takes place
template <typename TPixelType>
TPixelType*** ftkLaplacianOfGaussian3D<TPixelType>::convolveGaussian(	double* kernel_X, double* kernel_Y, double* kernel_Z, TPixelType*** paddedImage, 
																		unsigned int padded_image_x_size, unsigned int padded_image_y_size, unsigned int padded_image_z_size,  //padded_image_size should be the image size + padding for gaussian kernel + LoGSmoothenessPadding
																		unsigned int kernel_size_X, unsigned int kernel_size_Y, unsigned int kernel_size_Z)
{
	std::cout << "Convolving Gaussian" << std::endl;

	int padding_X = kernel_size_X / 2 + LoGSmoothnessPadding;
	int padding_Y = kernel_size_Y / 2 + LoGSmoothnessPadding;
	int padding_Z = kernel_size_Z / 2 + LoGSmoothnessPadding;

	//Allocate memory for temporary padded image (to hold result after each direction of convolution)
	TPixelType*** tempPaddedImage = (TPixelType ***) malloc(padded_image_x_size * sizeof(TPixelType**));
	for (int k = 0; k < padded_image_x_size; k++)
	{	
		tempPaddedImage[k] = (TPixelType **) malloc(padded_image_y_size * sizeof(TPixelType *));
		for (int l = 0; l < padded_image_y_size; l++)
		{
			tempPaddedImage[k][l] = (TPixelType *) malloc(padded_image_z_size * sizeof(TPixelType));
			for (int m = 0; m < padded_image_z_size; m++)
			{
				tempPaddedImage[k][l][m] = paddedImage[k][l][m];
			}
		}
	}

	//convolution

	//convolving in x direction

	#pragma omp parallel for
	for(int k = 0; k < image_x_size + 2 * LoGSmoothnessPadding; k++)
	{        
		for(int l = 0; l < image_y_size + 2 * LoGSmoothnessPadding; l++)
		{			
			for(int m = 0; m < image_z_size + 2 * LoGSmoothnessPadding; m++)
			{	
				tempPaddedImage[k + padding_X - LoGSmoothnessPadding][l + padding_Y - LoGSmoothnessPadding][m + padding_Z - LoGSmoothnessPadding] = sumOfProductX(kernel_X, paddedImage, k, l + padding_Y, m + padding_Z, kernel_size_X);
			}
		}
	}
	//std::cout << "Done convolving in x-direction" << std::endl;

	//convolving in y-direction
	#pragma omp parallel for
	for(int k = 0; k < image_x_size + 2 * LoGSmoothnessPadding; k++)
	{        
		for(int l = 0; l < image_y_size + 2 * LoGSmoothnessPadding; l++)
		{			
			for(int m = 0; m < image_z_size + 2 * LoGSmoothnessPadding; m++)
			{				
				paddedImage[k + padding_X - LoGSmoothnessPadding][l + padding_Y - LoGSmoothnessPadding][m + padding_Z - LoGSmoothnessPadding] = sumOfProductY(kernel_Y, tempPaddedImage, k + padding_X, l, m + padding_Z, kernel_size_Y);
			}
		}
	}

	//convolving in z-direction
	#pragma omp parallel for
	for(int k = 0; k < image_x_size + 2 * LoGSmoothnessPadding; k++)
	{        
		for(int l = 0; l < image_y_size + 2 * LoGSmoothnessPadding; l++)
		{			
			for(int m = 0; m < image_z_size + 2 * LoGSmoothnessPadding; m++)
			{				
				tempPaddedImage[k + padding_X - LoGSmoothnessPadding][l + padding_Y - LoGSmoothnessPadding][m + padding_Z - LoGSmoothnessPadding] = sumOfProductZ(kernel_Z, paddedImage, k + padding_X, l + padding_Y, m, kernel_size_Z);
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

	std::cout << "Printing out convolved image" << std::endl;
	for (int m = 0; m < padded_image_z_size; m++)
	{
		for (int l = 0; l < padded_image_y_size; l++)
		{
			for (int k = 0; k < padded_image_x_size; k++)
			{
				std::cout << tempPaddedImage[k][l][m] << " ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	std::cout << std::endl;

	return tempPaddedImage;
}

//Sum of the products in the x-direction, to be used with the Gaussian convolution
template <typename TPixelType>
TPixelType ftkLaplacianOfGaussian3D<TPixelType>::sumOfProductX(double* kernel, TPixelType*** paddedImage, unsigned int padded_image_start_x, unsigned int padded_image_start_y, unsigned int padded_image_start_z, unsigned int kernel_size)
{
	double sum = 0;

	for (int k = 0; k < kernel_size; k++)
		sum += kernel[k] * paddedImage[padded_image_start_x + k][padded_image_start_y][padded_image_start_z];

	/*if (sum / 65535.0 > 1)
		std::cout << "Gaussian X sum exceeds range" << std::endl;*/
	
	return sum;
}

//Sum of the products in the y-direction, to be used with the Gaussian convolution
template <typename TPixelType>
TPixelType ftkLaplacianOfGaussian3D<TPixelType>::sumOfProductY(double* kernel, TPixelType*** paddedImage, unsigned int padded_image_start_x, unsigned int padded_image_start_y, unsigned int padded_image_start_z, unsigned int kernel_size)
{
	double sum = 0;

	for (int k = 0; k < kernel_size; k++)
		sum += kernel[k] * paddedImage[padded_image_start_x][padded_image_start_y + k][padded_image_start_z];
	
	/*if (sum / 65535.0 > 1)
		std::cout << "Gaussian Y sum exceeds range" << std::endl;*/
	
	return sum;
}

//Sum of the products in the z-direction, to be used with the Gaussian convolution
template <typename TPixelType>
TPixelType ftkLaplacianOfGaussian3D<TPixelType>::sumOfProductZ(double* kernel, TPixelType*** paddedImage, unsigned int padded_image_start_x, unsigned int padded_image_start_y, unsigned int padded_image_start_z, unsigned int kernel_size)
{
	double sum = 0;

	for (int k = 0; k < kernel_size; k++)
	{
		if (kernel[k] * paddedImage[padded_image_start_x][padded_image_start_y][padded_image_start_z + k] / 65535.0 > 1)
		{
			std::cout << "Gaussian Z exceeds range" << std::endl;
			std::cout << "kernel[" << k << "]: " << kernel[k] << " paddedImage[" << padded_image_start_x << "][" << padded_image_start_y << "][" << padded_image_start_z << "]: " << paddedImage[padded_image_start_x][padded_image_start_y][padded_image_start_z + k] << std::endl;
		}
		sum += kernel[k] * paddedImage[padded_image_start_x][padded_image_start_y][padded_image_start_z + k];
	}
	
	/*if (sum / 65535.0 > 1)
		std::cout << "Gaussian Z sum exceeds range" << std::endl;*/
	
	return sum;
}

//Convolve with the Laplacian
template <typename TPixelType>
TPixelType*** ftkLaplacianOfGaussian3D<TPixelType>::convolveLaplacian(double*** kernel, TPixelType*** padded_image, unsigned int padded_image_x_size, unsigned int padded_image_y_size, unsigned int padded_image_z_size)
{
	int image_x_size = padded_image_x_size - 2; 
	int image_y_size = padded_image_y_size - 2;
	int image_z_size = padded_image_z_size - 2;

	//Allocate memory for the output of the Laplacian convolution	
	TPixelType*** outputImage = (TPixelType ***) malloc(padded_image_x_size * sizeof(TPixelType**));
	for (int k = 0; k < padded_image_x_size; k++)
	{	
		outputImage[k] = (TPixelType **) malloc(padded_image_y_size * sizeof(TPixelType *));
		for (int l = 0; l < padded_image_y_size; l++)
		{
			outputImage[k][l] = (TPixelType *) malloc(padded_image_z_size * sizeof(TPixelType));
			for (int m = 0; m < padded_image_z_size; m++)
			{
				outputImage[k][l][m] = 0;
				//outputImage[k][l][m] = padded_image[k][l][m]; //for bypassing Laplacian convolution
			}
		}
	}
	
	std::cout << "Convolving Laplacian" << std::endl;

	//if TPixelType can take negative values in its range, then center around 0, otherwise center around half the range (this is needed because if you used unsigned char and short, the Laplacian convolution may take negative values (which signed types cannot represent, so we shift the center up to the middle of the range of those types)
	size_t center;
	if (std::numeric_limits<TPixelType>::is_signed)
		center = 0;
	else
		center = (std::numeric_limits<TPixelType>::max() + std::numeric_limits<TPixelType>::min())/2;

	std::cout << "Center of TPixelType range is: " << center << std::endl;
	//std::cout << "Convolving Laplacian over: " << image_x_size << "x" << image_y_size << "x" << image_z_size << std::endl;
	
	#pragma omp parallel for
	for(int k = 0; k < image_x_size; k++)
	{        
		for(int l = 0; l < image_y_size; l++)
		{			
			for(int m = 0; m < image_z_size; m++)
			{				
				outputImage[k + 1][l + 1][m + 1] = sumOfProduct(kernel, padded_image, k, l, m, 3) + center;
			}
		}
	}

	return outputImage;
}

//Sum of Product over the entire 3D volume as defined by padded_image_start_x/y/z and the kernel_size
template <typename TPixelType>
double ftkLaplacianOfGaussian3D<TPixelType>::sumOfProduct(double*** kernel, TPixelType*** paddedImage, unsigned int padded_image_start_x, unsigned int padded_image_start_y, unsigned int padded_image_start_z, unsigned int kernel_size)
{
	double sum = 0;

	for (int k = 0; k < kernel_size; k++)
		for (int l = 0; l < kernel_size; l++)
			for (int m = 0; m < kernel_size; m++)
				sum += kernel[k][l][m] * paddedImage[padded_image_start_x + k][padded_image_start_y + l][padded_image_start_z + m]; //note that this is reduction

	return sum;
}

//Unpad the image by padding_X/Y/Z which represents the half the total padding in each direction
template <typename TPixelType>
TPixelType*** ftkLaplacianOfGaussian3D<TPixelType>::unpadImage(TPixelType*** padded_image, unsigned int image_x_size, unsigned int image_y_size, unsigned int image_z_size, unsigned int padding_X, unsigned int padding_Y, unsigned int padding_Z)
{
	//Allocate memory for the output of the unpadding
	TPixelType*** outputImage = (TPixelType ***) malloc(image_x_size * sizeof(TPixelType**));
	for (int k = 0; k < image_x_size; k++)
	{	
		outputImage[k] = (TPixelType **) malloc(image_y_size * sizeof(TPixelType *));
		for (int l = 0; l < image_y_size; l++)
		{
			outputImage[k][l] = (TPixelType *) malloc(image_z_size * sizeof(TPixelType));
			for (int m = 0; m < image_z_size; m++)
			{
				outputImage[k][l][m] = 0;
			}
		}
	}

	//Copy the output pixels from the input pixels
	for (int k = 0; k < image_x_size; k++)
		for (int l = 0; l < image_y_size; l++)
			for (int m = 0; m < image_z_size; m++)
				outputImage[k][l][m] = padded_image[k + padding_X][l + padding_Y][m + padding_Z];
	
	return outputImage;
}


//Destructor
template <typename TPixelType>
ftkLaplacianOfGaussian3D<TPixelType>::~ftkLaplacianOfGaussian3D()
{
	free(kernel_X);
	free(kernel_Y);
	free(kernel_Z);

	for (int k = 0; k < 3; k++)
	{
		for (int l = 0; l < 3; l++)
			free(laplacian_kernel[k][l]);
		free(laplacian_kernel[k]);
	}
	free(laplacian_kernel);

	freeImageMem(LoGImage, image_x_size, image_y_size);
}