
/*#if CUDA
#include "Convolution_CUDA.cuh"
#endif*/

template <typename TPixelType>
class ftkLaplacianOfGaussian3D
{
//Class variables
private:
	TPixelType ***image;
	
	float scale_X;
	float scale_Y;
	float scale_Z;

	unsigned int image_x_size;
	unsigned int image_y_size;
	unsigned int image_z_size;
	
	//Gaussian kernel
	double *kernel_X, *kernel_Y, *kernel_Z;
	unsigned int gaussian_kernel_size_x;
	unsigned int gaussian_kernel_size_y;
	unsigned int gaussian_kernel_size_z;

	//Laplacian kernel
	double ***laplacian_kernel;
	unsigned int laplacian_kernel_size_x;
	unsigned int laplacian_kernel_size_y;
	unsigned int laplacian_kernel_size_z;

	TPixelType ***LoGImage;

//Functions
public:	
	//image is 3-dimensional array holding the image that needs the LoG filter run on. You can specify the smoothing scale in each direction separately.
	ftkLaplacianOfGaussian3D(TPixelType*** image, float scale_X, float scale_Y, float scale_Z, unsigned int image_x_size, unsigned int image_y_size, unsigned int image_z_size);	//Construct
	void			RunFilter();																																					//Run the filter
	TPixelType***	GetOutput();																																					//Get the output, the original image is not affected

private:
	//Purposely put in private as no one should be calling default constructor/destructor
	ftkLaplacianOfGaussian3D();	
	~ftkLaplacianOfGaussian3D();

	//The private functions
	double*			generateGaussianKernel(float scale, unsigned int kernel_size);
	TPixelType***	padImage(TPixelType*** image, unsigned int image_x_size, unsigned int image_y_size, unsigned int image_z_size, unsigned int kernel_size_X, unsigned int kernel_size_Y, unsigned int kernel_size_Z);
	TPixelType***	convolveGaussian(double* kernel_X, double* kernel_Y, double* kernel_Z, TPixelType*** paddedImage, unsigned int padded_image_x_size, unsigned int padded_image_y_size, unsigned int padded_image_z_size, unsigned int kernel_size_X, unsigned int kernel_size_Y, unsigned int kernel_size_Z);
	TPixelType		sumOfProductX(double* kernel, TPixelType*** paddedImage, unsigned int padded_image_start_x, unsigned int padded_image_start_y, unsigned int padded_image_start_z, unsigned int kernel_size);
	TPixelType		sumOfProductY(double* kernel, TPixelType*** paddedImage, unsigned int padded_image_start_x, unsigned int padded_image_start_y, unsigned int padded_image_start_z, unsigned int kernel_size);
	TPixelType		sumOfProductZ(double* kernel, TPixelType*** paddedImage, unsigned int padded_image_start_x, unsigned int padded_image_start_y, unsigned int padded_image_start_z, unsigned int kernel_size);
	double***		generateLaplacianKernel();
	TPixelType***	convolveLaplacian(double*** kernel, TPixelType*** image, unsigned int padded_image_x_size, unsigned int padded_image_y_size, unsigned int padded_image_z_size);
	double			sumOfProduct(double*** kernel, TPixelType*** paddedImage, unsigned int padded_image_start_x, unsigned int padded_image_start_y, unsigned int padded_image_start_z, unsigned int kernel_size);
	TPixelType***	unpadImage(TPixelType*** padded_image, unsigned int image_x_size, unsigned int image_y_size, unsigned int image_z_size, unsigned int padding_X, unsigned int padding_Y, unsigned int padding_Z);
	void			freeImageMem(TPixelType*** image, unsigned int image_x_size, unsigned int image_y_size);
};

//Explicit instantiation of class
template class ftkLaplacianOfGaussian3D<double>;
template class ftkLaplacianOfGaussian3D<float>;
//template class ftkLaplacianOfGaussian3D<unsigned short>;	//works, but you may see weird results, since low intensity "resolution"
//template class ftkLaplacianOfGaussian3D<unsigned char>;		//works, but you may see weird results, since low intensity "resolution"

