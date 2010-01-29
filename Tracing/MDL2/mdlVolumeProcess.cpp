/*=========================================================================
Copyright 2009 Rensselaer Polytechnic Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 
=========================================================================*/
#if defined(_MSC_VER)
#pragma warning(disable : 4996)
#endif

#include "mdlVolumeProcess.h"

namespace mdl
{

//Constructor
VolumeProcess::VolumeProcess()
{
	m_inputImage = NULL;
	m_outputImage = NULL;
	useCAD = true;
	overwriteInput = true;
	debug = false;
}

void VolumeProcess::SetInput(ImageType::Pointer inImage)
{
	m_inputImage = inImage;
}

bool VolumeProcess::Update()
{
	if(!m_inputImage)
		return false;

	//Rescale the pixel values
    typedef itk::RescaleIntensityImageFilter<ImageType, ImageType>  RescaleFilterType;
    RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
    rescaleFilter->SetInput( m_inputImage );
	rescaleFilter->SetOutputMinimum(  0 );
    rescaleFilter->SetOutputMaximum( 255 );
	try
	{
		rescaleFilter->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		if(debug)
		{
			std::cerr << "ITK FILTER ERROR: " << err << std::endl ;
		}
		return false;
	}
	m_outputImage = rescaleFilter->GetOutput();

	if(debug)
	{
		std::cerr << "The Linear Mapping is done\n" << std::endl;
	}

	//----------------------Curvature based diffusion ---------------//
    if(useCAD) 
	{
		typedef itk::CurvatureAnisotropicDiffusionImageFilter< ImageType, ImageType > MCD_FilterType;
		MCD_FilterType::Pointer MCDFilter = MCD_FilterType::New();
    
		//Initialnization,  using the paper's optimal parameters
		const unsigned int numberOfIterations = 5;
		const double       timeStep = 0.0425;
		const double       conductance = 3;
		MCDFilter->SetNumberOfIterations(numberOfIterations);
		MCDFilter->SetTimeStep( timeStep );
		MCDFilter->SetConductanceParameter( conductance );
		MCDFilter->SetInput( m_inputImage );

		try
		{
			MCDFilter->Update();
		}
		catch( itk::ExceptionObject & err )
		{
			if(debug)
			{
				std::cerr << "ITK FILTER ERROR: " << err << std::endl ;
			}
			return false;
		}
		m_outputImage = MCDFilter->GetOutput();

		
		if(debug)
		{
			std::cerr << "The Curvature Diffusion is done\n" << std::endl;
		}
	}

    double itkThreshold = getItkOtsuThreshold(m_outputImage);
	if(debug)
	{
		std::cerr << "itk Threshold is " << itkThreshold << std::endl;
	}
    //------------------------------------------------------------------------------//

	// by xiao liang, using 3 sigma theory to estimate threshold;
	double xlThreshold = getXiaoLiangOtsuThreshold(m_outputImage);
	if(debug)
	{
		std::cerr << "Xiao Threshold = " << xlThreshold << std::endl;
	}

	//Need to compute the optimal threshold
	double threshold;
	if (itkThreshold > 0 && itkThreshold < 12)
	{
		threshold = itkThreshold;
		//threshold = itkThreshold / 3.3;
	}
	else
	{
		threshold = xlThreshold;
	}

	if(debug)	
		std::cerr << "OTSU optimal threshold " << threshold << std::endl;

	//Apply threshold
	itk::ImageRegionIterator< ImageType > itr( m_outputImage, m_outputImage->GetBufferedRegion() );
	for(itr.GoToBegin(); !itr.IsAtEnd(); ++itr)
	{
		if(itr.Get() < threshold) 
        {
          itr.Set(0);
        }
	}

	//Dilation of the image
	this->dialateImage(m_outputImage, 1, 3);



	//cout << "Done" << endl;
	return true;
}

VolumeProcess::ImageType::Pointer VolumeProcess::GetOutput()
{
	return m_outputImage;
}

double VolumeProcess::getItkOtsuThreshold(ImageType::Pointer img)
{
	double thresh = 0;

	typedef itk::OtsuThresholdImageFilter< ImageType, ImageType > OTSUFilterType;
    OTSUFilterType::Pointer OTSUFilter = OTSUFilterType::New();
	OTSUFilter->SetInput(img);
	OTSUFilter->SetNumberOfHistogramBins(m_NumberOfHistogramBins);
	try
	{
		OTSUFilter->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		if(debug)
		{
			std::cerr << "ITK FILTER ERROR: " << err << std::endl ;
		}
		return thresh;
	}

    thresh = OTSUFilter->GetThreshold();
	return thresh;
}

// This function is designed to compute the optimal threshold using OTSU method;
// this algoritm is implemented by xiao liang based on ITK's OTSU algorithm
double VolumeProcess::getXiaoLiangOtsuThreshold(ImageType::Pointer img)
{
	if(!img)
		return 0;

	double threshold = 0;
	unsigned char m_min = 255;
	unsigned char m_max = 0;
	double m_mean = 0.0;
	double m_variance = 0.0;

	ImageType::RegionType region = img->GetBufferedRegion();
	double numPix = region.GetSize(2)*region.GetSize(1)*region.GetSize(0);

	//Get min, max, and mean:
	itk::ImageRegionIterator< ImageType > itr( img, region );
	for(itr.GoToBegin(); !itr.IsAtEnd(); ++itr)
	{
		double val = itr.Get();
		if(val > m_max) m_max = val;
		if(val < m_min) m_min = val;
		m_mean += itr.Get();
	}
	m_mean  = m_mean/numPix;

	if(debug)
		std::cerr << "Max = " << (int)m_max << ", Min = " << (int)m_min << std::endl;

	//Do a sanity check
	if ( m_min >= m_max)
    {
		threshold=m_min;
		return threshold;
    }

	//Get the variance:
	for(itr.GoToBegin(); !itr.IsAtEnd(); ++itr)
	{
		double val = (double)itr.Get();
		m_variance += (val-m_mean)*(val-m_mean);
	}
	//These were not Xiao Liang's version:
	//m_variance = m_variance / numPix;
	//m_variance = sqrt(m_variance);

	threshold = m_mean - (m_variance/30); 
    // this step is only initialized a good experimental value for m_Threshold, because the 3D image
    // is sparse, there are lots of zero values; 

	//Create a histogram & init to zero
	double relativeFrequency[m_NumberOfHistogramBins];
	for ( unsigned char j = 0; j < m_NumberOfHistogramBins; j++ )
    {
       relativeFrequency[j] = 0.0;
    }

	double binMultiplier = (double)m_NumberOfHistogramBins/(double)(m_max-m_min);
	if(debug)
		std::cerr << "binMultiplier = " << binMultiplier << std::endl;

	unsigned int binNumber;
	for(itr.GoToBegin(); !itr.IsAtEnd(); ++itr)
	{
		double val = itr.Get();
		if ( val == m_min ) 
		{
			binNumber = 0;
		}
		else
		{
			binNumber = (unsigned int)(((val-m_min)*binMultiplier) - 1);
			if ( binNumber == m_NumberOfHistogramBins ) // in case of rounding errors
			{
				binNumber -= 1;
			}
		}
		relativeFrequency[binNumber] += 1.0;
	}

	// normalize the frequencies
	double totalMean = 0.0;
	for ( unsigned char j = 0; j < m_NumberOfHistogramBins; j++ )
    {
		relativeFrequency[j] /= numPix;
		totalMean += (j+1) * relativeFrequency[j];
    }

	// compute Otsu's threshold by maximizing the between-class variance
	double freqLeft = relativeFrequency[0];
	double meanLeft = 1.0;
	double meanRight = ( totalMean - freqLeft ) / ( 1.0 - freqLeft );

	double maxVarBetween = freqLeft * ( 1.0 - freqLeft ) * sqrt( meanLeft - meanRight );
	int maxBinNumber = 0;

	double freqLeftOld = freqLeft;
	double meanLeftOld = meanLeft;

	for ( unsigned char j = 1; j < m_NumberOfHistogramBins; j++ )
    {
		freqLeft += relativeFrequency[j];
		meanLeft = ( ((meanLeftOld * freqLeftOld)+(j+1)) * relativeFrequency[j] ) / freqLeft;
		if (freqLeft == 1.0)
		{
			meanRight = 0.0;
		}
		else
		{
			meanRight = ( totalMean - meanLeft * freqLeft ) / ( 1.0 - freqLeft );
		}
		
		double varBetween = freqLeft * ( 1.0 - freqLeft ) * sqrt( meanLeft - meanRight );
		if ( varBetween > maxVarBetween )
		{
			maxVarBetween = varBetween;
			maxBinNumber = j;
		}

		// cache old values
		freqLeftOld = freqLeft;
		meanLeftOld = meanLeft; 
    } 

	threshold = double( m_min + ( maxBinNumber + 1 ) / binMultiplier );
	return threshold; 
}

void VolumeProcess::dialateImage(ImageType::Pointer img, int iterations, int border )
{
	ImageType::Pointer out_img = ImageType::New();
	out_img->SetRegions( img->GetLargestPossibleRegion() );
	out_img->Allocate();
	out_img->FillBuffer(0);

	ImageType::RegionType fullRegion = img->GetBufferedRegion();
	ImageType::RegionType borderRegion;
	borderRegion.SetIndex(0, fullRegion.GetIndex(0)+border);
	borderRegion.SetIndex(1, fullRegion.GetIndex(1)+border);
	borderRegion.SetIndex(2, fullRegion.GetIndex(2)+border);
	borderRegion.SetSize(0, fullRegion.GetSize(0)-border);
	borderRegion.SetSize(1, fullRegion.GetSize(1)-border);
	borderRegion.SetSize(2, fullRegion.GetSize(2)-border);
	itk::ImageRegionIteratorWithIndex< ImageType > itr( img, borderRegion );
	while(iterations > 0)
	{
		for(itr.GoToBegin(); !itr.IsAtEnd(); ++itr)
		{
			double blockMax = itr.Get();
			for(int k=-1; k<=1; k++)
			{
				for(int j=-1; j<=1; j++)
				{
					for(int i=-1; i<=1; i++)
					{
						ImageType::IndexType index = itr.GetIndex();
						index[0] += i;
						index[1] += j;
						index[2] += k;
						ImageType::PixelType pix = img->GetPixel(index);
						if((double)pix > blockMax) 
						{
							blockMax = (double)pix;
						}
					}
				}
			}
			// Keep the peak of the original intensity
			if (blockMax == itr.Get() && blockMax != 0)
            {
				blockMax = blockMax + 1;
            }
			out_img->SetPixel(itr.GetIndex(), blockMax);
		}
		//Copy out_img back to in image for next dialation
		for(itr.GoToBegin(); !itr.IsAtEnd(); ++itr)
		{
			img->SetPixel(itr.GetIndex(), itr.Get());
		}
		
		iterations--;
	}
}


}



