#ifndef ASPECTRATIORESAMPLER_H
#define ASPECTRATIORESAMPLER_H

#include "itkImage.h"

class AspectRatioResampler
{
public:
	enum SamplingType
	{
		DownSample,
		UpSample
	};

public:
	template< typename TImageType > 
	static typename TImageType::Pointer ResampleImage(typename TImageType::Pointer image, const float aspect_ratio, SamplingType sampling_type);
	
	template< typename TImageType >
	static typename TImageType::Pointer UnsampleImage(typename TImageType::Pointer image, const float aspect_ratio, SamplingType sampling_type);
};

#include "AspectRatioResampler.txx"

#endif //ASPECTRATIORESAMPLER_H