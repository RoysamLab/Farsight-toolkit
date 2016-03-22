#include "Common.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace Common
{

void ReadImage3D(const std::string file_name, ImageType3D::Pointer& data_ptr){

	std::cout << "Reading input file " << file_name << std::endl;

	typedef itk::ImageFileReader<ImageType3D> ReaderType;
	//ReaderType::GlobalWarningDisplayOff();
	ReaderType::Pointer image_reader = ReaderType::New();
	image_reader->SetFileName(file_name);
	image_reader->Update();

	data_ptr = image_reader->GetOutput();
	data_ptr->Update();

	std::cout << "Input file size: " << data_ptr->GetBufferedRegion().GetSize() << std::endl;
}

void ReadImage3DUChar(const std::string file_name, RenderImageType3D::Pointer& data_ptr){

	std::cout << "Reading input file " << file_name << std::endl;

	typedef itk::ImageFileReader<RenderImageType3D> ReaderType;
	//ReaderType::GlobalWarningDisplayOff();
	ReaderType::Pointer image_reader = ReaderType::New();
	image_reader->SetFileName(file_name);
	image_reader->Update();

	data_ptr = image_reader->GetOutput();
	data_ptr->Update();

	std::cout << "Input file size: " << data_ptr->GetBufferedRegion().GetSize() << std::endl;
}

void ReadImage3DUShort(const std::string file_name, LabelImageType3D::Pointer& data_ptr){

	std::cout << "Reading input file " << file_name << std::endl;

	typedef itk::ImageFileReader<LabelImageType3D> ReaderType;
	//ReaderType::GlobalWarningDisplayOff();
	ReaderType::Pointer image_reader = ReaderType::New();
	image_reader->SetFileName(file_name);
	image_reader->Update();

	data_ptr = image_reader->GetOutput();
	data_ptr->Update();

	std::cout << "Input file size: " << data_ptr->GetBufferedRegion().GetSize() << std::endl;
}

void WriteTIFFImage3D(const std::string file_name, ImageType3D::Pointer& data_ptr){

	std::cout << "Writing output file" << file_name << std::endl;

	typedef itk::Image<unsigned char, 3> CharImageType3D;
	typedef itk::RescaleIntensityImageFilter<ImageType3D, CharImageType3D> RescalerType;
	RescalerType::Pointer rescaler = RescalerType::New();
	rescaler->SetOutputMinimum(0);
	rescaler->SetOutputMaximum(255);
	
	typedef itk::ImageFileWriter<CharImageType3D> WriterType;
	//WriterType::GlobalWarningDisplayOff();
	WriterType::Pointer image_writer = WriterType::New();
	image_writer->SetFileName(file_name);
	
	rescaler->SetInput(data_ptr);
	image_writer->SetInput(rescaler->GetOutput());
	image_writer->Update();
}

void CurvatureAnisotropicDiffusion(int& IterCount, int& conductance, ImageType3D::Pointer& data_ptr){

	itk::Size<3> radius;
	radius.Fill(1);
	
	int gaussian_variance = 2;
	typedef itk::DiscreteGaussianImageFilter<ImageType3D, ImageType3D> GaussType;
	GaussType::Pointer Gaussian_filter = GaussType::New();
	Gaussian_filter->SetVariance(gaussian_variance);
	Gaussian_filter->SetUseImageSpacingOn();
	
	
	itk::ImageRegionIterator<ImageType3D> iter_data(data_ptr, data_ptr->GetBufferedRegion());

	// Finding extrema of the data
	float max_val = -1000.00;
	float min_val = 1000.00;
	for (iter_data.GoToBegin(); !iter_data.IsAtEnd(); ++iter_data){
		max_val = (iter_data.Get() < max_val) ? max_val : iter_data.Get();
		min_val = (iter_data.Get() > min_val) ? min_val : iter_data.Get();
	}

	ImageType3D::Pointer temp = ImageType3D::New();
	temp->SetRegions(data_ptr->GetBufferedRegion());
	temp->Allocate();

	ImageType3D::Pointer data_ptr2 = ImageType3D::New();
	
	data_ptr2->SetRegions(data_ptr->GetBufferedRegion());
	data_ptr2->Allocate();
	itk::ImageRegionIterator<ImageType3D> iter_data2(data_ptr2, data_ptr2->GetBufferedRegion());
	itk::NeighborhoodIterator<ImageType3D> iter_neighbor(radius, data_ptr2, data_ptr2->GetBufferedRegion());

	//itk::Offset<3> f100 = {{-1,0,0}}, f010 = {{0,-1,0}}, f210 = {{-1,0, 1}}, f120 = {{ }}, f110 = {{ }},

	//unsigned int IterCount = 10; //why 10??
	const float K = (max_val - min_val)/conductance; //30.0; //why 30?? Is this conductance?
	std::cout << "K = " << K << std::endl;
  
	float zAspectRatio = data_ptr->GetSpacing()[2] / data_ptr->GetSpacing()[1];
	std::cout << "Z Aspect Ratio = " << zAspectRatio << std::endl;

	while (IterCount > 0)	{

		// copying data from vol to temp?  
		for (iter_data.GoToBegin(), iter_data2.GoToBegin(); !iter_data.IsAtEnd(); ++iter_data, ++iter_data2)	{
		  float val = iter_data.Get();
		  iter_data2.Set(val);
		}

		//should not work, filter is inactive..
		//WriteFloatImage3D(std::string("DiscreteGaussSmoothing.mhd"), data_ptr2);

		float dxf, dxb, dyf, dyb, dzf, dzb, dx, dy, dz;
		float gxf, gxb, gyf, gyb, gzf, gzb;
		float od1, od2;
		float propagation_gradient, speed;

		// Curvature anisotropic diffusion filter implementation - Why not use the one in ITK?  
		for (iter_data.GoToBegin(), iter_neighbor.GoToBegin(); !iter_data.IsAtEnd(); ++iter_data, ++iter_neighbor)	{
			//if (it.InBounds())	{
				  // 6 7 8    15 16 17   24 25 26
				  // 3 4 5    12 13 14   21 22 23
				  // 0 1 2    9  10 11   18 19 20
				
			dx = iter_neighbor.GetPixel(14) - iter_neighbor.GetPixel(12);
			dxf = iter_neighbor.GetPixel(14) - iter_neighbor.GetPixel(13);
			dxb = iter_neighbor.GetPixel(13) - iter_neighbor.GetPixel(12);
			dy = iter_neighbor.GetPixel(10) - iter_neighbor.GetPixel(16);
			dyf = iter_neighbor.GetPixel(16) - iter_neighbor.GetPixel(13);
			dyb = iter_neighbor.GetPixel(13) - iter_neighbor.GetPixel(10);
			dz = (iter_neighbor.GetPixel(4) - iter_neighbor.GetPixel(22))/zAspectRatio;
			dzf = (iter_neighbor.GetPixel(22) - iter_neighbor.GetPixel(13))/zAspectRatio;
			dzb = (iter_neighbor.GetPixel(13) - iter_neighbor.GetPixel(22))/zAspectRatio;

			 //x 12, 14
			od1 = iter_neighbor.GetPixel(17) - iter_neighbor.GetPixel(11);
			od2 = (iter_neighbor.GetPixel(23) - iter_neighbor.GetPixel(5))/zAspectRatio;
			gxf = vcl_sqrt(dxf*dxf + 0.25*(od1 + dy)*(od1 + dy) + 0.25*(od2 + dz)*(od2 + dz)) + 0.000001;

			od1 = iter_neighbor.GetPixel(15) - iter_neighbor.GetPixel(9);
			od2 = (iter_neighbor.GetPixel(21) - iter_neighbor.GetPixel(3))/zAspectRatio;
			gxb = vcl_sqrt(dxb*dxb + 0.25*(od1 + dy)*(od1 + dy) + 0.25*(od2 + dz)*(od2 + dz)) + 0.000001;

			//y	16, 10
			od1 = iter_neighbor.GetPixel(17) - iter_neighbor.GetPixel(15);
			od2 = (iter_neighbor.GetPixel(25) - iter_neighbor.GetPixel(7))/zAspectRatio;
			gyf = vcl_sqrt(0.25*(od1 + dx)*(od1 + dx) + dyf*dyf + 0.25*(od2 + dz)*(od2 + dz)) + 0.000001;

			od1 = iter_neighbor.GetPixel(11) - iter_neighbor.GetPixel(9);
			od2 = (iter_neighbor.GetPixel(19) - iter_neighbor.GetPixel(1))/zAspectRatio;
			gyb = vcl_sqrt(0.25*(od1 + dx)*(od1 + dx) + dyb*dyb + 0.25*(od2 + dz)*(od2 + dz)) + 0.000001;

			//z	4, 22
			od1 = iter_neighbor.GetPixel(23) - iter_neighbor.GetPixel(21);
			od2 = iter_neighbor.GetPixel(25) - iter_neighbor.GetPixel(19);
			gzf = vcl_sqrt(0.25*(od1 + dx)*(od1 + dx) + 0.25*(od2 + dy)*(od2 + dy) + dzf*dzf) + 0.000001;

			od1 = iter_neighbor.GetPixel(5) - iter_neighbor.GetPixel(3);
			od2 = iter_neighbor.GetPixel(7) - iter_neighbor.GetPixel(1);
			gzb = vcl_sqrt(0.25*(od1 + dx)*(od1 + dx) + 0.25*(od2 + dy)*(od2 + dy) + dzb*dzb) + 0.000001;

			speed = (dxf / gxf * vcl_exp( -1*gxf / K )) - (dxb / gxb * vcl_exp( -1*gxb / K ) )
				  + (dyf / gyf * vcl_exp( -1*gyf / K )) - (dyb / gyb * vcl_exp( -1*gyb / K ) )
				  + (dzf / gzf * vcl_exp( -1*gzf / K )) - (dzb / gzb * vcl_exp( -1*gzb / K ) );


			if (speed > 0)	{
				propagation_gradient = vnl_math_sqr( vnl_math_min(dxb, 0.0f) ) + vnl_math_sqr( vnl_math_max(dxf,  0.0f) ) +
					vnl_math_sqr( vnl_math_min(dyb, 0.0f) ) + vnl_math_sqr( vnl_math_max(dyf,  0.0f) ) +
					vnl_math_sqr( vnl_math_min(dzb, 0.0f) ) + vnl_math_sqr( vnl_math_max(dzf,  0.0f) );
			}
			else{
				propagation_gradient = vnl_math_sqr( vnl_math_max(dxb, 0.0f) ) + vnl_math_sqr( vnl_math_min(dxf,  0.0f) ) +
					vnl_math_sqr( vnl_math_max(dyb, 0.0f) ) + vnl_math_sqr( vnl_math_min(dyf,  0.0f) ) +
					vnl_math_sqr( vnl_math_max(dzb, 0.0f) ) + vnl_math_sqr( vnl_math_min(dzf,  0.0f) );
			}

			float val = iter_data.Get() + 0.1 * vcl_sqrt(propagation_gradient) * speed ;
			val = vnl_math_max(val, min_val);
			val = vnl_math_min(val, max_val);
				
			iter_data.Set(val);
				
			temp->SetPixel(iter_data.GetIndex(),  vcl_sqrt(propagation_gradient) * speed);
			//}
		}

		//WriteImage3D(std::string("curveUpdates.mhd"), temp);
		std::cout << "Completed Iteration " << IterCount << std::endl;
		IterCount--;
	}

	typedef itk::RescaleIntensityImageFilter<ImageType3D, ImageType3D> RescalerType;
	RescalerType::Pointer rescaler = RescalerType::New();
	rescaler->SetOutputMinimum(0.0);
	rescaler->SetOutputMaximum(1.0);
	rescaler->SetInput(data_ptr);
	rescaler->Update();
	data_ptr = rescaler->GetOutput();

}

void MedianFilter(const int radius, ImageType3D::Pointer& data_ptr){
		
	typedef itk::MedianImageFilter<ImageType3D, ImageType3D> MedianFilterType;
	MedianFilterType::Pointer median_filter = MedianFilterType::New();
	
	//radius = 1;
	itk::Size<3> kernel = {radius, radius, radius}; 
	median_filter->SetRadius(kernel);
	
	typedef itk::RescaleIntensityImageFilter<ImageType3D, ImageType3D> RescalerType;
	RescalerType::Pointer rescaler = RescalerType::New();
	rescaler->SetOutputMinimum(0.0);
	rescaler->SetOutputMaximum(1.0);
	
	median_filter->SetInput(data_ptr);
	rescaler->SetInput(median_filter->GetOutput());
	//rescaler->Update();
	
	data_ptr = rescaler->GetOutput();
	data_ptr->Update();

	std::cout << "Input file size: "<< data_ptr->GetBufferedRegion().GetSize() << std::endl;
}

void GVFDiffusion(float& smoothing_sigma, int& N_iter, const std::string& write_path, ImageType3D::Pointer& data_ptr){
	
	typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<ImageType3D, ImageType3D> SmoothingFilterType;
	SmoothingFilterType::Pointer smoothing_filter = SmoothingFilterType::New();
	smoothing_filter->SetInput(data_ptr);
	smoothing_filter->SetSigma(smoothing_sigma);
	smoothing_filter->Update();

	ImageType3D::Pointer smoothing_output = smoothing_filter->GetOutput();
	smoothing_output->Update();

	itk::ImageRegionIterator<ImageType3D> iter_data(smoothing_output, smoothing_output->GetBufferedRegion());

	// Calculate min and max of the dataset
	float max = -1000.0f, min = 1000.0f;
    for (iter_data.GoToBegin(); !iter_data.IsAtEnd(); ++iter_data ) {
		max = vnl_math_max(max, iter_data.Get());
		min = vnl_math_min(min, iter_data.Get());
	}
	std::cout << "max : " << max << " min : " << min << std::endl;

	for (iter_data.GoToBegin(); !iter_data.IsAtEnd(); ++iter_data ) {
		float d = iter_data.Get();
		d = (d - min)/(max - min);
		d = vcl_pow(d, 2.0f);
		iter_data.Set(d);
	}

	//initialize gx, gy, gz
	ImageType3D::Pointer gx = ImageType3D::New();	gx->SetRegions(data_ptr->GetBufferedRegion());	gx->Allocate();
	ImageType3D::Pointer gy = ImageType3D::New();	gy->SetRegions(data_ptr->GetBufferedRegion());	gy->Allocate();
	ImageType3D::Pointer gz = ImageType3D::New();	gz->SetRegions(data_ptr->GetBufferedRegion());	gz->Allocate();
	ImageType3D::Pointer g1x = ImageType3D::New();	g1x->SetRegions(data_ptr->GetBufferedRegion());	g1x->Allocate();
	ImageType3D::Pointer g1y = ImageType3D::New();	g1y->SetRegions(data_ptr->GetBufferedRegion());	g1y->Allocate();
	ImageType3D::Pointer g1z = ImageType3D::New();	g1z->SetRegions(data_ptr->GetBufferedRegion());	g1z->Allocate();
    
    // These are for the derivatives in x, y, and z directions?
	itk::Size<3> rad = {{1,1,1}};
	itk::ImageRegionIterator<ImageType3D> x_iter(gx, gx->GetBufferedRegion());
	itk::ImageRegionIterator<ImageType3D> y_iter(gy, gy->GetBufferedRegion());
	itk::ImageRegionIterator<ImageType3D> z_iter(gz, gz->GetBufferedRegion());

    // Calculating the derivative image?
	itk::NeighborhoodIterator<ImageType3D> iter_nhood(rad , smoothing_output, smoothing_output->GetBufferedRegion());
	
	const double EPS = 0.001;
	for (iter_nhood.GoToBegin(), x_iter.GoToBegin(), y_iter.GoToBegin(), z_iter.GoToBegin(); !iter_nhood.IsAtEnd(); ++iter_nhood, ++x_iter, ++y_iter, ++z_iter)	{
		float dx = EpsilonClip(EPS, 0.5*(iter_nhood.GetNext(0) - iter_nhood.GetPrevious(0)));
		float dy = EpsilonClip(EPS, 0.5*(iter_nhood.GetNext(1) - iter_nhood.GetPrevious(1)));
		float dz = EpsilonClip(EPS, 0.5*(iter_nhood.GetNext(2) - iter_nhood.GetPrevious(2)));

		x_iter.Set(dx);
		y_iter.Set(dy);
		z_iter.Set(dz);
	}

	itk::Size<3> data_size = smoothing_output->GetBufferedRegion().GetSize();
    
    // Looks like a Laplacian?
	const double filter[27] = {	0.0206, 0.0339, 0.0206,  0.0339, 0.0560, 0.0339,  0.0206, 0.0339, 0.0206,  \
									0.0339, 0.0560, 0.0339,  0.0560, 0.0923, 0.0560,  0.0339, 0.0560, 0.0339,  \
	                                0.0206, 0.0339, 0.0206,  0.0339, 0.0560, 0.0339,  0.0206, 0.0339, 0.0206   \
								};
	//int N_iter = 30;
	for (int iter = 0; iter < N_iter; ++iter) {

		double update = 0.0;
		for (long k = 1; k < data_size[2]-1; ++k) {
			for (long j = 1; j < data_size[1]-1; ++j) {
				for (long i = 1; i < data_size[0]-1; ++i)	{

					itk::Index<3> ndx = {i,j,k};

					double dx = double(gx->GetPixel(ndx)), dy = double(gy->GetPixel(ndx)), dz = double(gz->GetPixel(ndx));
					double g = (dx*dx) + (dy*dy) + (dz*dz) + EPS; 

                    // Application of the above filter
					double d1x = 0.0f, d1y = 0.0f, d1z = 0.0f;
					int f = 0;
					for (long kk = k-1; kk <= k+1 ; ++kk) {
						for (long jj = j-1; jj <= j+1 ; ++jj) {
							for (long ii = i-1; ii <= i+1 ; ++ii) {
								itk::Index<3> ndx1 = {ii,jj,kk};
								d1x += double(gx->GetPixel(ndx1))*filter[f];
								d1y += double(gy->GetPixel(ndx1))*filter[f];
								d1z += double(gz->GetPixel(ndx1))*filter[f];
								++f;
							}
						}
					}

					double g1 = (d1x*d1x) + (d1y*d1y) + (d1z*d1z) + EPS;

					if (g1 > g )  {
						g1x->SetPixel(ndx, float(d1x));
						g1y->SetPixel(ndx, float(d1y));
						g1z->SetPixel(ndx, float(d1z));
						update += (g1 - g);
					}
					else {
						g1x->SetPixel(ndx, float(dx));
						g1y->SetPixel(ndx, float(dy));
						g1z->SetPixel(ndx, float(dz));
						g1 = g;
					}
				}
			}
		}

		std::cout << "Iter # " << iter << " Update :" << update << std::endl;

		for (long k = 1; k < data_size[2]-1; ++k) {
			for (long j = 1; j < data_size[1]-1; ++j) {
				for (long i = 1; i < data_size[0]-1; ++i)	{
					itk::Index<3> ndx = {i,j,k};
						gx->SetPixel(ndx, g1x->GetPixel(ndx));
						gy->SetPixel(ndx, g1y->GetPixel(ndx));
						gz->SetPixel(ndx, g1z->GetPixel(ndx));
				}
			}
		}
	}
    // Write the gradient images and the vol image to disk
	std::string output_file;
	output_file = write_path + "_gx.mhd";
	WriteImage3D(output_file, gx);
	output_file = write_path + "_gy.mhd";
	WriteImage3D(output_file, gy);
	output_file = write_path + "_gz.mhd";
	WriteImage3D(output_file, gz);
	
	//output_file = write_path + ".tif";
	//WriteImage3D(output_file, data_ptr);
}

void GVFDiffusionFaster(float& smoothing_sigma, int& N_iter, const std::string& write_path, ImageType3D::Pointer& data_ptr){
	
	typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<ImageType3D, ImageType3D> SmoothingFilterType;
	SmoothingFilterType::Pointer smoothing_filter = SmoothingFilterType::New();
	smoothing_filter->SetInput(data_ptr);
	smoothing_filter->SetSigma(smoothing_sigma);
	smoothing_filter->Update();

	ImageType3D::Pointer smoothing_output = smoothing_filter->GetOutput();
	smoothing_output->Update();

	StatisticsFilterType::Pointer stats_filter = StatisticsFilterType::New();
	stats_filter->SetInput(smoothing_output);
	stats_filter->Update();
	float max = stats_filter->GetMaximum();
	float min = stats_filter->GetMinimum();
	
	SubtractImageFilter::Pointer img_subtacter = SubtractImageFilter::New();
	img_subtacter->SetInput(smoothing_output);
	img_subtacter->SetConstant2(min);
	
	MultiplyImageFilter::Pointer img_multiplier = MultiplyImageFilter::New();
	img_multiplier->SetInput(img_subtacter->GetOutput());
	img_multiplier->SetConstant2(1.0/(max - min));

	SquareImageFilter::Pointer img_squarer = SquareImageFilter::New();
	img_squarer->SetInput(img_multiplier->GetOutput());
	img_squarer->Update();

	ImageType3D::Pointer smoothing_output1 = img_squarer->GetOutput();

	
	/*itk::ImageRegionIterator<ImageType3D> iter_data(smoothing_output, smoothing_output->GetBufferedRegion());

	// Calculate min and max of the dataset
	float max = -1000.0f, min = 1000.0f;
    for (iter_data.GoToBegin(); !iter_data.IsAtEnd(); ++iter_data ) {
		max = vnl_math_max(max, iter_data.Get());
		min = vnl_math_min(min, iter_data.Get());
	}
	std::cout << "max : " << max << " min : " << min << std::endl;

	for (iter_data.GoToBegin(); !iter_data.IsAtEnd(); ++iter_data ) {
		float d = iter_data.Get();
		d = (d - min)/(max - min);
		d = vcl_pow(d, 2.0f);
		iter_data.Set(d);
	}*/

	//initialize gx, gy, gz
	ImageType3D::Pointer gx = ImageType3D::New();	gx->SetRegions(data_ptr->GetBufferedRegion());	gx->Allocate();
	ImageType3D::Pointer gy = ImageType3D::New();	gy->SetRegions(data_ptr->GetBufferedRegion());	gy->Allocate();
	ImageType3D::Pointer gz = ImageType3D::New();	gz->SetRegions(data_ptr->GetBufferedRegion());	gz->Allocate();
	ImageType3D::Pointer g1x = ImageType3D::New();	g1x->SetRegions(data_ptr->GetBufferedRegion());	g1x->Allocate();
	ImageType3D::Pointer g1y = ImageType3D::New();	g1y->SetRegions(data_ptr->GetBufferedRegion());	g1y->Allocate();
	ImageType3D::Pointer g1z = ImageType3D::New();	g1z->SetRegions(data_ptr->GetBufferedRegion());	g1z->Allocate();
    
    // These are for the derivatives in x, y, and z directions?
	itk::Size<3> rad = {{1,1,1}};
	itk::ImageRegionIterator<ImageType3D> x_iter(gx, gx->GetBufferedRegion());
	itk::ImageRegionIterator<ImageType3D> y_iter(gy, gy->GetBufferedRegion());
	itk::ImageRegionIterator<ImageType3D> z_iter(gz, gz->GetBufferedRegion());

    // Calculating the derivative image?
	itk::NeighborhoodIterator<ImageType3D> iter_nhood(rad , smoothing_output1, smoothing_output1->GetBufferedRegion());
	
	const double EPS = 0.001;
	for (iter_nhood.GoToBegin(), x_iter.GoToBegin(), y_iter.GoToBegin(), z_iter.GoToBegin(); !iter_nhood.IsAtEnd(); ++iter_nhood, ++x_iter, ++y_iter, ++z_iter)	{
		float dx = EpsilonClip(EPS, 0.5*(iter_nhood.GetNext(0) - iter_nhood.GetPrevious(0)));
		float dy = EpsilonClip(EPS, 0.5*(iter_nhood.GetNext(1) - iter_nhood.GetPrevious(1)));
		float dz = EpsilonClip(EPS, 0.5*(iter_nhood.GetNext(2) - iter_nhood.GetPrevious(2)));

		x_iter.Set(dx);
		y_iter.Set(dy);
		z_iter.Set(dz);
	}

	itk::Size<3> data_size = smoothing_output1->GetBufferedRegion().GetSize();
    
    // Looks like a Laplacian?
	const double filter[27] = {	0.0206, 0.0339, 0.0206,  0.0339, 0.0560, 0.0339,  0.0206, 0.0339, 0.0206,  \
									0.0339, 0.0560, 0.0339,  0.0560, 0.0923, 0.0560,  0.0339, 0.0560, 0.0339,  \
	                                0.0206, 0.0339, 0.0206,  0.0339, 0.0560, 0.0339,  0.0206, 0.0339, 0.0206   \
								};
	//int N_iter = 30;
	for (int iter = 0; iter < N_iter; ++iter) {

		double update = 0.0;
		for (long k = 1; k < data_size[2]-1; ++k) {
			for (long j = 1; j < data_size[1]-1; ++j) {
				for (long i = 1; i < data_size[0]-1; ++i)	{

					itk::Index<3> ndx = {i,j,k};

					double dx = double(gx->GetPixel(ndx)), dy = double(gy->GetPixel(ndx)), dz = double(gz->GetPixel(ndx));
					double g = (dx*dx) + (dy*dy) + (dz*dz) + EPS; 

                    // Application of the above filter
					double d1x = 0.0f, d1y = 0.0f, d1z = 0.0f;
					int f = 0;
					for (long kk = k-1; kk <= k+1 ; ++kk) {
						for (long jj = j-1; jj <= j+1 ; ++jj) {
							for (long ii = i-1; ii <= i+1 ; ++ii) {
								itk::Index<3> ndx1 = {ii,jj,kk};
								d1x += double(gx->GetPixel(ndx1))*filter[f];
								d1y += double(gy->GetPixel(ndx1))*filter[f];
								d1z += double(gz->GetPixel(ndx1))*filter[f];
								++f;
							}
						}
					}

					double g1 = (d1x*d1x) + (d1y*d1y) + (d1z*d1z) + EPS;

					if (g1 > g )  {
						g1x->SetPixel(ndx, float(d1x));
						g1y->SetPixel(ndx, float(d1y));
						g1z->SetPixel(ndx, float(d1z));
						update += (g1 - g);
					}
					else {
						g1x->SetPixel(ndx, float(dx));
						g1y->SetPixel(ndx, float(dy));
						g1z->SetPixel(ndx, float(dz));
						g1 = g;
					}
				}
			}
		}

		std::cout << "Iter # " << iter << " Update :" << update << std::endl;

		DuplicatorType::Pointer img_dupx = DuplicatorType::New();
		img_dupx->SetInputImage(g1x);
		img_dupx->Update();
		
		DuplicatorType::Pointer img_dupy = DuplicatorType::New();
		img_dupy->SetInputImage(g1y);
		img_dupy->Update();

		DuplicatorType::Pointer img_dupz = DuplicatorType::New();
		img_dupz->SetInputImage(g1z);
		img_dupz->Update();

		gx = img_dupx->GetOutput();
		gy = img_dupy->GetOutput();
		gz = img_dupz->GetOutput();

		/*for (long k = 1; k < data_size[2]-1; ++k) {
			for (long j = 1; j < data_size[1]-1; ++j) {
				for (long i = 1; i < data_size[0]-1; ++i)	{
					itk::Index<3> ndx = {i,j,k};
						gx->SetPixel(ndx, g1x->GetPixel(ndx));
						gy->SetPixel(ndx, g1y->GetPixel(ndx));
						gz->SetPixel(ndx, g1z->GetPixel(ndx));
				}
			}
		}*/
	}

    // Write the gradient images and the vol image to disk
	std::string output_file;
	output_file = write_path + "_gx.mhd";
	WriteImage3D(output_file, gx);
	output_file = write_path + "_gy.mhd";
	WriteImage3D(output_file, gy);
	output_file = write_path + "_gz.mhd";
	WriteImage3D(output_file, gz);
	
	output_file = write_path + "_pre.mhd";
	WriteImage3D(output_file, data_ptr);
}

void GVFDiffusion(float& smoothing_sigma, int& N_iter, ImageType3D::Pointer& data_ptr, ImageType3D::Pointer& gx, ImageType3D::Pointer& gy, ImageType3D::Pointer& gz){

	typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<ImageType3D, ImageType3D> SmoothingFilterType;
	SmoothingFilterType::Pointer smoothing_filter = SmoothingFilterType::New();
	smoothing_filter->SetInput(data_ptr);
	smoothing_filter->SetSigma(smoothing_sigma);
	smoothing_filter->Update();

	ImageType3D::Pointer smoothing_output = smoothing_filter->GetOutput();
	smoothing_output->Update();

	itk::ImageRegionIterator<ImageType3D> iter_data(smoothing_output, smoothing_output->GetBufferedRegion());

	// Calculate min and max of the dataset
	float max = -1000.0f, min = 1000.0f;
    for (iter_data.GoToBegin(); !iter_data.IsAtEnd(); ++iter_data ) {
		max = vnl_math_max(max, iter_data.Get());
		min = vnl_math_min(min, iter_data.Get());
	}
	std::cout << "max : " << max << " min : " << min << std::endl;

	for (iter_data.GoToBegin(); !iter_data.IsAtEnd(); ++iter_data ) {
		float d = iter_data.Get();
		d = (d - min)/(max - min);
		d = vcl_pow(d, 2.0f);
		iter_data.Set(d);
	}

	//initialize gx, gy, gz
	gx = ImageType3D::New();	gx->SetRegions(data_ptr->GetBufferedRegion());	gx->Allocate();
	gy = ImageType3D::New();	gy->SetRegions(data_ptr->GetBufferedRegion());	gy->Allocate();
	gz = ImageType3D::New();	gz->SetRegions(data_ptr->GetBufferedRegion());	gz->Allocate();
	ImageType3D::Pointer g1x = ImageType3D::New();	g1x->SetRegions(data_ptr->GetBufferedRegion());	g1x->Allocate();
	ImageType3D::Pointer g1y = ImageType3D::New();	g1y->SetRegions(data_ptr->GetBufferedRegion());	g1y->Allocate();
	ImageType3D::Pointer g1z = ImageType3D::New();	g1z->SetRegions(data_ptr->GetBufferedRegion());	g1z->Allocate();
    
    // These are for the derivatives in x, y, and z directions?
	itk::Size<3> rad = {{1,1,1}};
	itk::ImageRegionIterator<ImageType3D> x_iter(gx, gx->GetBufferedRegion());
	itk::ImageRegionIterator<ImageType3D> y_iter(gy, gy->GetBufferedRegion());
	itk::ImageRegionIterator<ImageType3D> z_iter(gz, gz->GetBufferedRegion());

    // Calculating the derivative image?
	itk::NeighborhoodIterator<ImageType3D> iter_nhood(rad , smoothing_output, smoothing_output->GetBufferedRegion());
	
	const double EPS = 0.001;
	for (iter_nhood.GoToBegin(), x_iter.GoToBegin(), y_iter.GoToBegin(), z_iter.GoToBegin(); !iter_nhood.IsAtEnd(); ++iter_nhood, ++x_iter, ++y_iter, ++z_iter)	{
		float dx = EpsilonClip(EPS, 0.5*(iter_nhood.GetNext(0) - iter_nhood.GetPrevious(0)));
		float dy = EpsilonClip(EPS, 0.5*(iter_nhood.GetNext(1) - iter_nhood.GetPrevious(1)));
		float dz = EpsilonClip(EPS, 0.5*(iter_nhood.GetNext(2) - iter_nhood.GetPrevious(2)));

		x_iter.Set(dx);
		y_iter.Set(dy);
		z_iter.Set(dz);
	}

	itk::Size<3> data_size = smoothing_output->GetBufferedRegion().GetSize();
    
    // How is this kernel derived? and what is it meant for?
	const double filter[27] = {	0.0206, 0.0339, 0.0206,  0.0339, 0.0560, 0.0339,  0.0206, 0.0339, 0.0206,  \
									0.0339, 0.0560, 0.0339,  0.0560, 0.0923, 0.0560,  0.0339, 0.0560, 0.0339,  \
	                                0.0206, 0.0339, 0.0206,  0.0339, 0.0560, 0.0339,  0.0206, 0.0339, 0.0206   \
								};
	//int N_iter = 30;
	for (int iter = 0; iter < N_iter; ++iter) {

		double update = 0.0;
		for (long k = 1; k < data_size[2]-1; ++k) {
			for (long j = 1; j < data_size[1]-1; ++j) {
				for (long i = 1; i < data_size[0]-1; ++i)	{

					itk::Index<3> ndx = {i,j,k};

					double dx = double(gx->GetPixel(ndx)), dy = double(gy->GetPixel(ndx)), dz = double(gz->GetPixel(ndx));
					double g = (dx*dx) + (dy*dy) + (dz*dz) + EPS; 

                    // Application of the above filter
					double d1x = 0.0f, d1y = 0.0f, d1z = 0.0f;
					int f = 0;
					for (long kk = k-1; kk <= k+1 ; ++kk) {
						for (long jj = j-1; jj <= j+1 ; ++jj) {
							for (long ii = i-1; ii <= i+1 ; ++ii) {
								itk::Index<3> ndx1 = {ii,jj,kk};
								d1x += double(gx->GetPixel(ndx1))*filter[f];
								d1y += double(gy->GetPixel(ndx1))*filter[f];
								d1z += double(gz->GetPixel(ndx1))*filter[f];
								++f;
							}
						}
					}

					double g1 = (d1x*d1x) + (d1y*d1y) + (d1z*d1z) + EPS;

					if (g1 > g )  {
						g1x->SetPixel(ndx, float(d1x));
						g1y->SetPixel(ndx, float(d1y));
						g1z->SetPixel(ndx, float(d1z));
						update += (g1 - g);
					}
					else {
						g1x->SetPixel(ndx, float(dx));
						g1y->SetPixel(ndx, float(dy));
						g1z->SetPixel(ndx, float(dz));
						g1 = g;
					}
				}
			}
		}

		std::cout << "Iter # " << iter << " Update :" << update << std::endl;
		
		#pragma omp parallel for
		for (long k = 1; k < data_size[2]-1; ++k) {
			for (long j = 1; j < data_size[1]-1; ++j) {
				for (long i = 1; i < data_size[0]-1; ++i)	{
					itk::Index<3> ndx = {i,j,k};
						gx->SetPixel(ndx, g1x->GetPixel(ndx));
						gy->SetPixel(ndx, g1y->GetPixel(ndx));
						gz->SetPixel(ndx, g1z->GetPixel(ndx));
				}
			}
		}
	}
}

float inline EpsilonClip(double EPS, float x){
	if ( vnl_math_abs(x) < EPS ) {
		x = (x < 0.0) ? -1*EPS : EPS ;
	}
	return (x);
}

double inline EpsilonClip(double EPS, double x){
	if ( vnl_math_abs(x) < EPS ) {
		x = (x < 0.0) ? -1*EPS : EPS ;
	}
	return (x);
}

void WriteImage3D(const std::string file_name, ImageType3D::Pointer& data_ptr){

	std::cout << "Writing output file "<< file_name << std::endl;

	typedef itk::ImageFileWriter<ImageType3D> WriterType;
	WriterType::GlobalWarningDisplayOff();
	WriterType::Pointer image_writer = WriterType::New();
	image_writer->SetFileName(file_name);
	image_writer->SetInput(data_ptr);
	image_writer->Update();
}

void RenderImage3D(RenderImageType3D::Pointer data_ptr){
	
	ITKToVTKConnectorType::Pointer ITK_to_VTK_connector = ITKToVTKConnectorType::New();

	ITK_to_VTK_connector->SetInput(data_ptr);
	ITK_to_VTK_connector->Update();

	vtkSmartPointer<vtkImageData> vtk_image = ITK_to_VTK_connector->GetOutput();

	// Testing vtk image
	//vtk_image->PrintSelf(std::cout, vtkIndent(0));

	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	renderer->SetBackground(1.0, 1.0, 1.0);
	
	vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();
	renderer->SetActiveCamera(camera);

	vtkSmartPointer<vtkRenderWindow> render_window = vtkSmartPointer<vtkRenderWindow>::New();
	render_window->AddRenderer(renderer);
	
	vtkSmartPointer<vtkRenderWindowInteractor> render_window_interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	render_window_interactor->SetRenderWindow(render_window);
	
	vtkSmartPointer<vtkPiecewiseFunction> opacity_transfer_function = vtkSmartPointer<vtkPiecewiseFunction>::New();
	opacity_transfer_function->AddPoint(2, 0.0);
	opacity_transfer_function->AddPoint(10, 0.1);
	

	vtkSmartPointer<vtkColorTransferFunction> color_transfer_function = vtkSmartPointer<vtkColorTransferFunction>::New();
	color_transfer_function->AddRGBPoint(0.0, 0.0, 0.0, 0.0);
	color_transfer_function->AddRGBPoint((10 * 255)/100, 1, 1, 1); //blue
	color_transfer_function->AddRGBPoint((45 * 255)/100, 0, .01, 0); //green
	color_transfer_function->AddRGBPoint((150 * 255)/100, .01, 0, 0); //red
	
	vtkSmartPointer<vtkVolumeProperty> volume_property = vtkSmartPointer<vtkVolumeProperty>::New();
	volume_property->SetColor(color_transfer_function);
	volume_property->SetScalarOpacity(opacity_transfer_function);
	volume_property->ShadeOn();
	volume_property->SetInterpolationTypeToLinear();
	
	vtkSmartPointer<vtkVolumeRayCastCompositeFunction> composite_function = vtkSmartPointer<vtkVolumeRayCastCompositeFunction>::New();
	
	//vtkSmartPointer<vtkVolumeRayCastMapper> volume_mapper = vtkSmartPointer<vtkVolumeRayCastMapper>::New();
	vtkSmartPointer<vtkFixedPointVolumeRayCastMapper> volume_mapper = vtkSmartPointer<vtkFixedPointVolumeRayCastMapper>::New();
	//vtkSmartPointer<vtkOpenGLVolumeTextureMapper3D> volume_mapper = vtkSmartPointer<vtkOpenGLVolumeTextureMapper3D>::New();
	//vtkSmartPointer<vtkGPUVolumeRayCastMapper> volume_mapper = vtkSmartPointer<vtkGPUVolumeRayCastMapper>::New();
	//vtkSmartPointer<vtkOpenGLGPUVolumeRayCastMapper> volume_mapper = vtkSmartPointer<vtkOpenGLGPUVolumeRayCastMapper>::New();
	//vtkSmartPointer<vtkSmartVolumeMapper> volume_mapper = vtkSmartPointer<vtkSmartVolumeMapper>::New();

	//volume_mapper->SetVolumeRayCastFunction(composite_function);
	#if VTK_MAJOR_VERSION <= 5
		volume_mapper->SetInput(vtk_image);
	#else
		volume_mapper->SetInputData(vtk_image);
	#endif
	volume_mapper->SetBlendModeToComposite();
	//volume_mapper->SetBlendModeToMaximumIntensity();
	//volume_mapper->SetScalarMode(1);
	//volume_mapper->SetSampleDistance(0.2);
	

	vtkSmartPointer<vtkVolume> volume = vtkSmartPointer<vtkVolume>::New();
	volume->SetMapper(volume_mapper); 
	volume->SetProperty(volume_property);
	volume->SetPosition(0, 0, 0);
	volume->SetPickable(0);
	volume->Update();
	
	render_window->Render();
	renderer->AddVolume(volume);
	renderer->ResetCamera();
	//renderer->Update();
	
	render_window_interactor->Initialize();
	render_window->Render();
	render_window_interactor->Start();
}

void RescaleDataForRendering(ImageType3D::Pointer data, RenderImageType3D::Pointer& data_to_render){

	typedef itk::RescaleIntensityImageFilter<ImageType3D, RenderImageType3D> RescalerType;
	RescalerType::Pointer rescaler = RescalerType::New();
	rescaler->SetOutputMinimum(0);
	rescaler->SetOutputMaximum(255);
	rescaler->SetInput(data);
	rescaler->Update();

	data_to_render = rescaler->GetOutput();
}

void CastImageUCharToFloat(RenderImageType3D::Pointer uchar_ptr, ImageType3D::Pointer& float_ptr){
	
	CastFilterType::Pointer cast_filter = CastFilterType::New();
	cast_filter->SetInput(uchar_ptr);
	cast_filter->Update();

	float_ptr = cast_filter->GetOutput();
}

PixelType NormalizeData(ImageType3D::Pointer inputData, ImageType3D::Pointer& normalizedInputData){

	// This filter gives wierd values at output !! Try the LabelStatisticsImageFilter
	/*StatisticsFilterType::Pointer statistics_filter = StatisticsFilterType::New();
	statistics_filter->SetInput(inputData);
	float volumeMean = statistics_filter->GetMean();
	float volumeStd = statistics_filter->GetSigma();
	float volumeMax = statistics_filter->GetMaximum();
	float volumeMin = statistics_filter->GetMinimum();
	statistics_filter->Update();*/

	normalizedInputData = ImageType3D::New();

	MinMaxCalculatorType::Pointer min_max_calculator = MinMaxCalculatorType::New();
	min_max_calculator->SetImage(inputData);
	min_max_calculator->Compute();

	PixelType volumeMax = min_max_calculator->GetMaximum();

	/*DuplicatorType::Pointer image_duplicator = DuplicatorType::New();
	image_duplicator->SetInputImage(inputData);
	image_duplicator->Update();
	normalizedInputData = image_duplicator->GetOutput();*/
	
	if(volumeMax != 1.0){

		DuplicatorType::Pointer image_duplicator2 = DuplicatorType::New();
		image_duplicator2->SetInputImage(inputData);
		image_duplicator2->Update();

		ImageType3D::Pointer divisor_image = image_duplicator2->GetOutput();
		divisor_image->FillBuffer(volumeMax);

		/*ImageType3D::Pointer img = ImageType3D::New();
		ImageType3D::IndexType start_idx; 
		start_idx[0] = 0; start_idx[1] = 0; start_idx[2] = 0;
		ImageType3D::SizeType size;
		size = inputData->getLargestPossibleRegion().GetSize();
		ImageType3D::RegionType region;
		region.SetSize(size); region.SetIndex(start);

		img->SetRegions(region);
		img->Allocate();*/

		DivideImageFilterType::Pointer divide_image_filter = DivideImageFilterType::New();
		divide_image_filter->SetInput1(inputData);
		divide_image_filter->SetInput2(divisor_image);
		divide_image_filter->Update();

		normalizedInputData = divide_image_filter->GetOutput();


		MinMaxCalculatorType::Pointer min_max_calculator2 = MinMaxCalculatorType::New();
		min_max_calculator2->SetImage(normalizedInputData);
		min_max_calculator2->Compute();

		PixelType volumeMax2 = min_max_calculator2->GetMaximum();

		/*typedef itk::Functor::Mult<ImageType3D::PixelType> MultiplicationFunctorType;
		typedef itk::UnaryFunctorImageFilter<ImageType3D, ImageType3D, MultiplicationFunctorType> MultiplicationUnaryFunctorFilterType;
		MultiplicationUnaryFunctorFilterType::Pointer multiplier = MultiplicationUnaryFunctorFilterType::New();
		
		typedef itk::BinaryFunctorImageFilter<ImageType3D, ImageType3D, ImageType3D, MultiplicationFunctorType> MultiplicationBinaryFunctorType;*/
		
		/*ShiftScaleFilterType::Pointer shift_scale_filter = ShiftScaleFilterType::New();
		shift_scale_filter->SetInput(inputData);
		shift_scale_filter->SetScale(1 / volumeMax);
		shift_scale_filter->Update();
		normalizedInputData = shift_scale_filter->GetOutput();*/

		// Too slow!!
		/*itk::ImageRegionIterator<ImageType3D> iter_input(inputData, inputData->GetBufferedRegion());
		itk::ImageRegionIterator<ImageType3D> iter_norm_input(normalizedInputData, normalizedInputData->GetBufferedRegion());
		
		while(!iter_input.IsAtEnd())
			iter_norm_input.Set(iter_input.Get() / volumeMax);*/
	}
	else
		normalizedInputData = inputData;

	return volumeMax;
}

} // end namespace Common

/*int inline GetSign(double value){
	if(value < 0)
		return -1;
	if(value > 0)
		return 1;
	return 0;
}*/
