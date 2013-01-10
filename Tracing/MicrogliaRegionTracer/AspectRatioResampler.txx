#include "itkIdentityTransform.h"
#include "itkResampleImageFilter.h"

#include "AspectRatioResampler.h"

template< typename TImageType>
typename TImageType::Pointer AspectRatioResampler::ResampleImage(typename TImageType::Pointer image, const float aspect_ratio, SamplingType sampling_type)
{
	std::cout << "Resampling Image" << std::endl;
	typename TImageType::SizeType inputSize = image->GetLargestPossibleRegion().GetSize();
	
	typename TImageType::SizeType outputSize;
	typename TImageType::SpacingType outputSpacing;
	if (sampling_type == UpSample)
	{
		outputSize[0] = inputSize[0] * aspect_ratio;
		outputSize[1] = inputSize[1] * aspect_ratio;
		outputSize[2] = inputSize[2];
		
		outputSpacing[0] = 1/aspect_ratio;
		outputSpacing[1] = 1/aspect_ratio;
		outputSpacing[2] = 1.0;
	}
	else 
	{
		outputSize[0] = inputSize[0];
		outputSize[1] = inputSize[1];
		outputSize[2] = inputSize[2] / aspect_ratio;
		
		outputSpacing[0] = 1.0;
		outputSpacing[1] = 1.0;
		outputSpacing[2] = aspect_ratio;
	}

	typedef itk::IdentityTransform< double, 3 > TransformType;
	typedef itk::ResampleImageFilter< TImageType, TImageType > ResampleImageFilterType;
	typename ResampleImageFilterType::Pointer resample_filter = ResampleImageFilterType::New();
	resample_filter->SetInput(image);
	resample_filter->SetSize(outputSize);
	resample_filter->SetOutputSpacing(outputSpacing);
	resample_filter->SetTransform(TransformType::New());
	try
	{
		resample_filter->UpdateLargestPossibleRegion();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "ResampleImage() resample_filter exception: " << err << std::endl;
	}
	typename TImageType::Pointer resampled_image = resample_filter->GetOutput();
	resampled_image->DisconnectPipeline();

	//Set the spacing so that 1 unit in physical space is 1 pixel
	typename TImageType::SpacingType spacing;
	spacing.Fill(1.0);
	resampled_image->SetSpacing(spacing);

	return resampled_image;
}

template< typename TImageType>
typename TImageType::Pointer AspectRatioResampler::UnsampleImage(typename TImageType::Pointer image, const float aspect_ratio, SamplingType sampling_type)
{
	typename TImageType::SizeType inputSize = image->GetLargestPossibleRegion().GetSize();
	
	typename TImageType::SizeType outputSize;
	typename TImageType::SpacingType outputSpacing;
	
	if (sampling_type == UpSample)
	{
		outputSize[0] = inputSize[0] / aspect_ratio;
		outputSize[1] = inputSize[1] / aspect_ratio;
		outputSize[2] = inputSize[2];
		
		outputSpacing[0] = aspect_ratio;
		outputSpacing[1] = aspect_ratio;
		outputSpacing[2] = 1.0;
	}
	else
	{
		outputSize[0] = inputSize[0];
		outputSize[1] = inputSize[1];
		outputSize[2] = inputSize[2] * aspect_ratio;
		
		outputSpacing[0] = 1.0;
		outputSpacing[1] = 1.0;
		outputSpacing[2] = 1.0 / aspect_ratio;
	}

	typedef itk::IdentityTransform< double, 3 > TransformType;
	typedef itk::ResampleImageFilter< TImageType, TImageType > ResampleImageFilterType;
	typename ResampleImageFilterType::Pointer resample_filter = ResampleImageFilterType::New();
	resample_filter->SetInput(image);
	resample_filter->SetSize(outputSize);
	resample_filter->SetOutputSpacing(outputSpacing);
	resample_filter->SetTransform(TransformType::New());
	try
	{
		resample_filter->UpdateLargestPossibleRegion();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "UnsampleImage() resample_filter exception: " << err << std::endl;
	}

	typename TImageType::Pointer unsampled_image = resample_filter->GetOutput();
	unsampled_image->DisconnectPipeline();	//Disconnect pipeline so we don't propagate...
	
	typename TImageType::SpacingType spacing;
	spacing.Fill(1.0);
	unsampled_image->SetSpacing(spacing);

	return unsampled_image;
}
