void sumOfProduct_CUDA(double* kernel, double* paddedImage, int outputImage_x_size, int outputImage_y_size, int outputImage_z_size, int kernel_size, double*** outputImage);
double* flattenImage(double*** image, int image_x_size, int image_y_size, int image_z_size);
