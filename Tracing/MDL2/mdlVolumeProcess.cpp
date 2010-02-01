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
	debug = false;
}

void VolumeProcess::SetInput(ImageType::Pointer inImage)
{
	m_inputImage = inImage;

	//Init the output image to the input image:
	m_outputImage = ImageType::New();
	m_outputImage->SetRegions( m_inputImage->GetLargestPossibleRegion() );
	m_outputImage->Allocate();
	m_outputImage->FillBuffer(0);
	itk::ImageRegionIterator< ImageType > itr1( m_inputImage, m_inputImage->GetLargestPossibleRegion() );
	itk::ImageRegionIterator< ImageType > itr2( m_outputImage, m_outputImage->GetLargestPossibleRegion() );
	for(itr1.GoToBegin(), itr2.GoToBegin() ; !itr1.IsAtEnd(); ++itr1, ++itr2)
	{
		itr2.Set( itr1.Get() );
	}
}

ImageType::Pointer VolumeProcess::GetOutput()
{
	return m_outputImage;
}

bool VolumeProcess::RescaleIntensities(int min, int max)
{
	//Rescale the pixel values
    typedef itk::RescaleIntensityImageFilter<ImageType, ImageType>  RescaleFilterType;
    RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
    rescaleFilter->SetInput( m_outputImage );
	rescaleFilter->InPlaceOn();
	rescaleFilter->SetOutputMinimum( min );
    rescaleFilter->SetOutputMaximum( max );
	try
	{
		rescaleFilter->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "ITK FILTER ERROR: " << err << std::endl ;
		return false;
	}
	m_outputImage = rescaleFilter->GetOutput();
	if(debug)
		std::cerr << "Rescale Filter Done" << std::endl;
	return true;
}

bool VolumeProcess::RunCAD()
{
	typedef itk::CurvatureAnisotropicDiffusionImageFilter< ImageType, ImageType > CADFilterType;
	CADFilterType::Pointer cadFilter = CADFilterType::New();
    
	//Initialnization,  using the paper's optimal parameters
	const unsigned int numberOfIterations = 5;
	const double       timeStep = 0.0425;
	const double       conductance = 3;
	cadFilter->SetNumberOfIterations(numberOfIterations);
	cadFilter->SetTimeStep( timeStep );
	cadFilter->SetConductanceParameter( conductance );
	cadFilter->InPlaceOn();
	cadFilter->SetInput( m_outputImage );
	try
	{
		cadFilter->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "ITK FILTER ERROR: " << err << std::endl ;
		return false;
	}
	m_outputImage = cadFilter->GetOutput();
	if(debug)
		std::cerr << "CAD Filter Done" << std::endl;
	return true;
}

bool VolumeProcess::DialateImage(int iterations)
{
	int border = 2;

	ImageType::Pointer tempImg = ImageType::New();
	tempImg->SetRegions( m_outputImage->GetLargestPossibleRegion() );
	tempImg->Allocate();
	tempImg->FillBuffer(0);

	ImageType::RegionType fullRegion = m_outputImage->GetBufferedRegion();
	ImageType::RegionType borderRegion;
	borderRegion.SetIndex(0, fullRegion.GetIndex(0)+border);
	borderRegion.SetIndex(1, fullRegion.GetIndex(1)+border);
	borderRegion.SetIndex(2, fullRegion.GetIndex(2)+border);
	borderRegion.SetSize(0, fullRegion.GetSize(0)-(2*border));
	borderRegion.SetSize(1, fullRegion.GetSize(1)-(2*border));
	borderRegion.SetSize(2, fullRegion.GetSize(2)-(2*border));
	itk::ImageRegionIteratorWithIndex< ImageType > itr( m_outputImage, borderRegion );
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
						ImageType::PixelType pix = m_outputImage->GetPixel(index);
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
			tempImg->SetPixel(itr.GetIndex(), blockMax);
		}

		//Copy temp img back to image for next dialation
		itk::ImageRegionIterator< ImageType > itr1( tempImg, tempImg->GetLargestPossibleRegion() );
		itk::ImageRegionIterator< ImageType > itr2( m_outputImage, m_outputImage->GetLargestPossibleRegion() );
		for(itr1.GoToBegin(), itr2.GoToBegin() ; !itr1.IsAtEnd(); ++itr1, ++itr2)
		{
			itr2.Set( itr1.Get() );
		}

		iterations--;
	}
	
	if(debug)
		std::cerr << "Dialation Done" << std::endl;
	return true;
}

//This function removes all objects that are <= minObjSize from the foreground.
//The foreground remains grayscale after this filter
bool VolumeProcess::MaskSmallConnComp(int minObjSize)
{
	typedef itk::Image< unsigned short, Dimension > ShortImageType;
	typedef itk::ConnectedComponentImageFilter< ImageType, ShortImageType > CCFilterType;
	typedef itk::RelabelComponentImageFilter< ShortImageType, ShortImageType > RelabelType;

	CCFilterType::Pointer ccfilter = CCFilterType::New();
	RelabelType::Pointer relabel = RelabelType::New();
	
	ccfilter->SetInput (m_outputImage);
	ccfilter->FullyConnectedOn();

	relabel->SetInput( ccfilter->GetOutput() );
	relabel->SetMinimumObjectSize( minObjSize );
	relabel->InPlaceOn();

	try
    {
		relabel->Update();
    }
    catch( itk::ExceptionObject & excep )
    {
		std::cerr << "Relabel: exception caught !" << std::endl;
		std::cerr << excep << std::endl;
		return false;
    }
	
	unsigned short numObjects = relabel->GetNumberOfObjects();
	if(debug)
		std::cerr << "Connected components = " << numObjects << std::endl;

	ShortImageType::Pointer ccImage = relabel->GetOutput();

	//Use connected component image as a mask:
	itk::ImageRegionIterator< ShortImageType > itr1( ccImage, ccImage->GetLargestPossibleRegion() );
	itk::ImageRegionIterator< ImageType > itr2( m_outputImage, m_outputImage->GetLargestPossibleRegion() );
	for(itr1.GoToBegin(), itr2.GoToBegin() ; !itr1.IsAtEnd(); ++itr1, ++itr2)
	{
		if(itr1.Get() == 0)
		{
			itr2.Set( 0 );
		}
	}
	return true;
}

bool VolumeProcess::MaskUsingGraphCuts()
{
	ImageType::RegionType region = m_outputImage->GetBufferedRegion();
	int numRows = region.GetSize(0);
	int numColumns = region.GetSize(1);
	int numStacks = region.GetSize(2);
	long numPix = numStacks*numColumns*numRows;

	unsigned short * binImagePtr = new unsigned short[numPix];
	unsigned char * dataImagePtr = m_outputImage->GetBufferPointer();

	int ok = Cell_Binarization_3D(dataImagePtr, binImagePtr, numRows, numColumns, numStacks, 0, 1);	//Do Binarization
	if(!ok)
		return false;

	//Mask out pixels in the background:
	for(long i=0; i<numPix; ++i)
	{
		if( binImagePtr[i] == 0 )
		{
			dataImagePtr[i] = 0;
		}
	}

	delete binImagePtr;

	return true;
}

bool VolumeProcess::RunOtsuDenoising()
{
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
		//threshold = xlThreshold;
		threshold = itkThreshold / 3.3;
	}

	if(debug)	
		std::cerr << "OTSU optimal threshold " << threshold << std::endl;
	
	//Apply threshold (any thing below threshold is set to zero)
	itk::ImageRegionIterator< ImageType > itr( m_outputImage, m_outputImage->GetLargestPossibleRegion() );
	for(itr.GoToBegin(); !itr.IsAtEnd(); ++itr)
	{
		if(itr.Get() < threshold) 
        {
			itr.Set(0);
        }
	}
	if(debug)
		std::cerr << "OTSU Filter Done" << std::endl;
	return true;
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



}



