#include "Vesselness.h"

Vesselness::Vesselness()
{
}

void Vesselness::RunMultiscaleVesselness(ImageType::Pointer image)
{
	for (double scale = 2.0; scale < 5.0; scale+=0.1)
	{
		VesselnessImageType::Pointer vesselness_image = RunVesselness(image, scale);
		std::ostringstream vesselnessFileNameStream;

		vesselnessFileNameStream << "Vesselness_" << scale << ".mhd";
		WriteVesselnessImage(vesselnessFileNameStream.str(), vesselness_image);
	}
}

Vesselness::VesselnessImageType::Pointer Vesselness::RunVesselness(ImageType::Pointer image, double scale)
{
	//Need the maximum intensity of the original image for a normalized Vesselness score
	typedef itk::MinimumMaximumImageCalculator< ImageType > ImageCalculatorFilterType;
	ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New();
	imageCalculatorFilter->SetImage(image);
	imageCalculatorFilter->ComputeMaximum();
	double max_intensity = imageCalculatorFilter->GetMaximum();
	
	//Allocate memory for a Vesselness image
	VesselnessImageType::Pointer vesselness_image = VesselnessImageType::New();
	VesselnessImageType::IndexType start;
	start.Fill(0);

	VesselnessImageType::SizeType size = image->GetLargestPossibleRegion().GetSize();

	VesselnessImageType::RegionType region(start, size);
	vesselness_image->SetRegions(region);
	vesselness_image->Allocate();
	vesselness_image->FillBuffer(0);
	
	//Run the Hessian
	Hessian *hessian_filter = new Hessian();
	HessianImageType::Pointer hessian_image = hessian_filter->RunHessian(image, scale);

	//Make iterator to read Hessian image
	itk::ImageRegionConstIterator< HessianImageType > hessian_image_iter(hessian_image, hessian_image->GetLargestPossibleRegion());
		
	//Make iterator to write Vesselness image
	itk::ImageRegionIterator< VesselnessImageType > vesselness_image_iter(vesselness_image, vesselness_image->GetLargestPossibleRegion());
	
	vesselness_image_iter.GoToBegin();
	hessian_image_iter.GoToBegin();

	while (!hessian_image_iter.IsAtEnd())
	{
		HessianTensorType tensor = hessian_image_iter.Get();

		HessianTensorType::EigenValuesArrayType ev;
		HessianTensorType::EigenVectorsMatrixType em;

		tensor.ComputeEigenAnalysis(ev, em);	//Compute Eigenvalues
		
		double vesselness_score = ComputeVesselness(ev[0], ev[1], ev[2], max_intensity);

		vesselness_image_iter.Set(vesselness_score);

		++vesselness_image_iter;
		++hessian_image_iter;
		}

	return vesselness_image;
}

double Vesselness::ComputeVesselness( HessianTensorType::EigenValuesArrayType::ValueType ev1, HessianTensorType::EigenValuesArrayType::ValueType ev2, HessianTensorType::EigenValuesArrayType::ValueType ev3, double maximum_intensity )
{
	double lambda1, lambda2, lambda3; //constrain lambda1 <= lambda2 <= lambda3

	double ev1_magnitude = std::abs(ev1);
	double ev2_magnitude = std::abs(ev2);
	double ev3_magnitude = std::abs(ev3);


	///std::cout << ev1_magnitude << " " << ev2_magnitude << " " << ev3_magnitude << std::endl;

	if (ev1_magnitude <= ev2_magnitude && ev1_magnitude <= ev3_magnitude) //ev1 is the smallest eigenvalue
	{
		lambda1 = ev1;
		if (ev2_magnitude < ev3_magnitude)
		{
			lambda2 = ev2;
			lambda3 = ev3;
		}
		else
		{
			lambda2 = ev3;
			lambda3 = ev2;
		}
	}
	else if (ev2_magnitude <= ev1_magnitude && ev2_magnitude <= ev3_magnitude) //ev2 is the smallest eigenvalue
	{
		lambda1 = ev2;
		if (ev1_magnitude < ev3_magnitude)
		{
			lambda2 = ev1;
			lambda3 = ev3;
		}
		else
		{
			lambda2 = ev3;
			lambda3 = ev1;
		}
	}
	else
	{
		lambda1 = ev3;
		if (ev1_magnitude < ev2_magnitude)
		{
			lambda2 = ev1;
			lambda3 = ev2;
		}
		else
		{
			lambda2 = ev2;
			lambda3 = ev1;
		}
	}


	double vesselness_score;
	if (lambda3 < 0 && lambda2 < 0)
	{
		lambda1 = std::abs(lambda1);
		lambda2 = std::abs(lambda2);
		lambda3 = std::abs(lambda3);

		double r_a = lambda1 / sqrt(lambda2 * lambda3);
		double r_b = lambda2 / lambda3;
		double s = sqrt(pow(lambda1, 2) + pow(lambda2, 2) + pow(lambda3, 2));

		//std::cout << r_a << " " << r_b << " " << s << std::endl;

		double alpha = 0.5 * maximum_intensity;
		double beta = 0.5 * maximum_intensity;
		double gamma = 0.25 * maximum_intensity;

		vesselness_score =			exp(-1 * pow(r_a, 2) / (2 * pow(alpha, 2))) 
			* (1 -	exp(-1 * pow(r_b, 2) / (2 * pow(beta, 2)))) 
			* (1 -	exp(-1 * pow(s, 2) / (2 * pow(gamma, 2))));

		//std::cout << vesselness_score << std::endl;
	}
	else
		vesselness_score = 0;

	return vesselness_score;
}

void Vesselness::WriteVesselnessImage(std::string filename, VesselnessImageType::Pointer image)
{
	typedef itk::ImageFileWriter< VesselnessImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(image);
	writer->SetFileName(filename);
	try
	{
		std::cout << "Writing Vesselness image: " << filename << std::endl;
		writer->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "writer Exception: " << err << std::endl;
	}
}