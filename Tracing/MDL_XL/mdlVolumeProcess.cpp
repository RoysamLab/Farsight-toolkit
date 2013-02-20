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
    //m_outputImage = inImage;	//THIS IS NOT CORRECT IT WILL CAUSE INPUT IMAGE TO CHANGE

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

bool VolumeProcess:: RunFillingZeroOnBouandary(int Bx,int By,int Bz)
{
	ImageType::Pointer tempImg = ImageType::New();
	tempImg->SetRegions( m_outputImage->GetLargestPossibleRegion() );
	tempImg->Allocate();
	tempImg->FillBuffer(0);

	ImageType::RegionType fullRegion = m_outputImage->GetBufferedRegion();

	int numColumns = fullRegion.GetSize(0);
	int numRows = fullRegion.GetSize(1);
	int numStacks = fullRegion.GetSize(2);
	int i, j,k; 

	itk::ImageRegionIteratorWithIndex< ImageType > itr( m_outputImage, fullRegion );

	for(itr.GoToBegin(); !itr.IsAtEnd(); ++itr)
	{
		
		ImageType::IndexType index = itr.GetIndex();
	    ImageType::PixelType pix = m_outputImage->GetPixel(index);
        i = index[0];
		j = index[1];
		k = index[2];
		if (i < Bx || i > numColumns -Bx || j < By || j > numRows-By ||k <Bz ||k > numStacks-Bz )
		    tempImg->SetPixel(itr.GetIndex(), 0);
		else 
			tempImg->SetPixel(itr.GetIndex(), pix);
     }
		//Copy temp img back to image 
	itk::ImageRegionIterator< ImageType > itr1( tempImg, tempImg->GetLargestPossibleRegion());
	itk::ImageRegionIterator< ImageType > itr2( m_outputImage, m_outputImage->GetLargestPossibleRegion());
	for(itr1.GoToBegin(), itr2.GoToBegin() ; !itr1.IsAtEnd(); ++itr1, ++itr2)
	{
	  itr2.Set( itr1.Get() );
	}

	if(debug)
		std::cerr << "RunFillingZero Done" << std::endl;

	return true;
  
}

ImageType::Pointer VolumeProcess::GetOutput()
{
	return m_outputImage;
}

bool VolumeProcess::RescaleIntensities(int min, int max)
{
	//Rescale the pixel values,//by xiao liang
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

bool VolumeProcess:: RunIntensityNormalize()
{
	//Rescale weights:
    typedef itk::NormalizeImageFilter< ImageType, FloatImageType3D> NormalizeImageFilter;
	NormalizeImageFilter::Pointer Normalize = NormalizeImageFilter::New();
	Normalize->SetInput(m_outputImage);
	Normalize->Update();
	m_outputImage = RescaleFloatToImageType(Normalize->GetOutput());
	return true;

}

ImageType::Pointer VolumeProcess::RescaleFloatToImageType(FloatImageType3D::Pointer img)
{
	//Rescale weights:
	typedef itk::RescaleIntensityImageFilter< FloatImageType3D, ImageType> RescaleType;
	RescaleType::Pointer rescale = RescaleType::New();
	rescale->SetOutputMaximum( 255 );
	rescale->SetOutputMinimum( 0 );
	rescale->SetInput( img );
	rescale->Update();
	return rescale->GetOutput();
}

bool VolumeProcess::RunGaussianSmoothing(float varX, float varY, float varZ, float maxErr)
{  //by xiao liang
   typedef itk::DiscreteGaussianImageFilter< ImageType, FloatImageType3D > GaussianFilterType;
   GaussianFilterType::Pointer GaussianFilter =  GaussianFilterType::New();
   GaussianFilter->SetInput(m_outputImage);
   //GaussianFilter->SetFilterDimensionality(3);
   GaussianFilterType::ArrayType maxErrTypeValue;
   maxErrTypeValue.Fill(maxErr);
   GaussianFilter->SetMaximumError( maxErrTypeValue );

   GaussianFilterType::ArrayType variance;
   variance[0] = varX;
   variance[1] = varY;
   variance[2] = varZ;
   GaussianFilter->SetVariance(variance);
   //GaussianFilter->SetMaximumKernelWidth(maxKernalWidth);
   try
    {
		GaussianFilter->Update();
    }
   catch( itk::ExceptionObject & err )
	{
		std::cerr << "ITK FILTER ERROR: " << err << std::endl ;
		return false;
	}
	m_outputImage = RescaleFloatToImageType(GaussianFilter->GetOutput());
	if(debug)
		std::cerr << "GaussianFilter Filter Done" << std::endl;
	return true;
}


bool VolumeProcess::RunRecursiveGaussianIIRRilter(float sigmaX, float sigmaY, float sigmaZ)
{  //by xiao liang
   typedef itk::RecursiveGaussianImageFilter< ImageType, ImageType > RecursiveGaussianFilterType;
   RecursiveGaussianFilterType::Pointer RecursiveGaussianFilterX =  RecursiveGaussianFilterType::New();
   RecursiveGaussianFilterType::Pointer RecursiveGaussianFilterY =  RecursiveGaussianFilterType::New();
   RecursiveGaussianFilterType::Pointer RecursiveGaussianFilterZ =  RecursiveGaussianFilterType::New();

   RecursiveGaussianFilterX->SetDirection(0);
   RecursiveGaussianFilterY->SetDirection(1);
   RecursiveGaussianFilterZ->SetDirection(2);
   
   RecursiveGaussianFilterX->SetOrder(RecursiveGaussianFilterType::ZeroOrder);
   RecursiveGaussianFilterY->SetOrder(RecursiveGaussianFilterType::ZeroOrder);
   RecursiveGaussianFilterZ->SetOrder(RecursiveGaussianFilterType::ZeroOrder);

   RecursiveGaussianFilterX->SetNormalizeAcrossScale(false);
   RecursiveGaussianFilterY->SetNormalizeAcrossScale(false);
   RecursiveGaussianFilterZ->SetNormalizeAcrossScale(false);

   RecursiveGaussianFilterX->SetSigma(sigmaX);
   RecursiveGaussianFilterX->SetSigma(sigmaY);
   RecursiveGaussianFilterX->SetSigma(sigmaZ);

   RecursiveGaussianFilterX->SetInput(m_outputImage);
   RecursiveGaussianFilterY->SetInput(RecursiveGaussianFilterX->GetOutput());
   RecursiveGaussianFilterZ->SetInput(RecursiveGaussianFilterY->GetOutput());
  
 
   try
    {
		RecursiveGaussianFilterZ->Update();
    }
   catch( itk::ExceptionObject & err )
	{
		std::cerr << "ITK FILTER ERROR: " << err << std::endl ;
		return false;
	}
	m_outputImage = RecursiveGaussianFilterZ->GetOutput();
	if(debug)
		std::cerr << " RecursiveGaussianFilter Done" << std::endl;
	return true;
}

bool VolumeProcess::RunCAD(unsigned int numberOfIterations,double timeStep,double conductance)
{
	typedef itk::CurvatureAnisotropicDiffusionImageFilter< ImageType, ImageType > CADFilterType;
	CADFilterType::Pointer cadFilter = CADFilterType::New();
    
	//Initialnization,  using the paper's optimal parameters
	//const unsigned int numberOfIterations = 5;
	//const double       timeStep = 0.0425;
	//const double       conductance = 3;
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
		std::cerr << "Curvature Anisotropic Diffusion Filter Done" << std::endl;
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
	int numColumns = region.GetSize(0);
	int numRows = region.GetSize(1);
	int numStacks = region.GetSize(2);
	long numPix = numStacks*numColumns*numRows;

	unsigned short * binImagePtr = new unsigned short[numPix];
	PixelType * dataImagePtr = m_outputImage->GetBufferPointer();

	//int ok = Cell_Binarization_3D(dataImagePtr, binImagePtr, numRows, numColumns, numStacks, 0, 1);	//Do Binarization
    int ok = Neuron_Binarization_3D(dataImagePtr, binImagePtr, numRows, numColumns, numStacks, 0, 1);
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

	delete [] binImagePtr;

	return true;
}

bool VolumeProcess::RunBinaryForDistanceMapUsingManualThreshold(float threshold)
{
// this function is corrsponding to the Xiaosong's original MDL code
	if(debug)
	{
		std::cerr << "We excute the Manual Thrsholding by " << threshold << std::endl;
	}
    
	//Apply threshold (any thing below threshold is set to zero)
	itk::ImageRegionIterator< ImageType > itr( m_outputImage, m_outputImage->GetLargestPossibleRegion() );
	for(itr.GoToBegin(); !itr.IsAtEnd(); ++itr)
	{
		if(itr.Get() < threshold) 
        {
			itr.Set(255);
        }
		else 
		{
			itr.Set(0);
		}
	}
	return true;
}
bool VolumeProcess::RunBinaryForDistanceMapUsingGraphCuts()
{
	ImageType::RegionType region = m_outputImage->GetBufferedRegion();
	int numColumns = region.GetSize(0);
	int numRows = region.GetSize(1);
	int numStacks = region.GetSize(2);
	long numPix = numStacks*numColumns*numRows;

	unsigned short * binImagePtr = new unsigned short[numPix];
	PixelType * dataImagePtr = m_outputImage->GetBufferPointer();

	//int ok = Cell_Binarization_3D(dataImagePtr, binImagePtr, numRows, numColumns, numStacks, 0, 1);	//Do Binarization
    int ok = Neuron_Binarization_3D(dataImagePtr, binImagePtr, numRows, numColumns, numStacks, 0, 1);
	if(!ok)
		return false;

	// background is bright and forground is dark:
	for(long i=0; i<numPix; ++i)
	{
		if( binImagePtr[i] == 0 )
		{
			dataImagePtr[i] = 255;
		}
		else 
            dataImagePtr[i] = 0;
	}

	delete [] binImagePtr;

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

//This code was written by Xiaosong Yuan and modified by Xiao Liang
bool VolumeProcess::RunAnisotropicDiffusion(int timesDiffuse, bool iso)
{
	ImageType::RegionType region = m_outputImage->GetLargestPossibleRegion();
	int sizeX = region.GetSize(0);
	int sizeY = region.GetSize(1);
	int sizeZ = region.GetSize(2);
	int sls = sizeX*sizeY;		//slice size pixels
	int sz = sls*sizeZ;			//total number of pixels
	
	//Create Gradient vector Image:
	struct mdlVector{float xd; float yd; float zd;};
	mdlVector *gradVec;
	gradVec = (mdlVector *)malloc(sz*sizeof(mdlVector));

	//Get the buffer to the input data:
	PixelType * volin = (PixelType*)malloc(sz*sizeof(PixelType));
	PixelType * input = m_outputImage->GetBufferPointer();
	memcpy(volin,input,sz*sizeof(PixelType));

	//Create buffer of output data:
	PixelType * volout = (PixelType*)malloc(sz*sizeof(PixelType));

	// define positive half kernel of derivative 
	double kernelWeight[3][3];
	kernelWeight[0][0] = 1; kernelWeight[0][1] = 2; kernelWeight[0][2] = 1;
	kernelWeight[1][0] = 2; kernelWeight[1][1] = 3; kernelWeight[1][2] = 2;
	kernelWeight[2][0] = 1; kernelWeight[2][1] = 2; kernelWeight[2][2] = 1;

	int border = 1;
	double AveGradient=0.0 ;  // by xiao liang
	double tempGradient;

	//timesDiffuse = 1;
	while (timesDiffuse > 0 )
	{
		if(debug)
			std::cerr << "Anisotropic Diffusion #" << timesDiffuse << " ..." << std::endl;

		//initial to zeros
		for (int idx=0; idx<sz; idx++)
		{
			volout[idx] = 0; 
		}

		// Compute Gradient
		for (int k = border; k < sizeZ-border; k++)
		{
			for (int j = border; j < sizeY-border; j++)
			{
				for (int i = border; i < sizeX-border; i++)
				{
					if(volin[k*sls + j*sizeX + i]>0) // by xiao
					{
						int idx = k*sls	+ j*sizeX + i;
						gradVec[idx].xd = 0;
						gradVec[idx].yd = 0;
						gradVec[idx].zd = 0;
                
						for (int d2 = -1; d2 <= 1; d2++)
						{
							for (int d1 = -1; d1 <= 1; d1++) 
							{
								int iidx1 = (k+d2) *sls + (j+d1) *sizeX + i-1;
								int iidx2 = (k+d2) *sls + (j+d1) *sizeX + i+1;
								gradVec[idx].xd = (float) (gradVec[idx].xd + kernelWeight[d2+1][d1+1]*(volin[iidx2] - volin[iidx1]));

								iidx1 = (k+d2) *sls + (j-1) *sizeX + (i+d1);
								iidx2 = (k+d2) *sls + (j+1) *sizeX + (i+d1);
								gradVec[idx].yd =(float) (gradVec[idx].yd + kernelWeight[d2+1][d1+1]*(volin[iidx2] - volin[iidx1]));

								iidx1 = (k-1) *sls + (j+d2) *sizeX + (i+d1);
								iidx2 = (k+1) *sls + (j+d2) *sizeX + (i+d1);
								gradVec[idx].zd = (float)(gradVec[idx].zd + kernelWeight[d2+1][d1+1]*(volin[iidx2] - volin[iidx1]));
							}
						}

						tempGradient = gradVec[idx].xd * gradVec[idx].xd + gradVec[idx].yd *gradVec[idx].yd+ gradVec[idx].zd * gradVec[idx].zd;
						tempGradient = sqrt(tempGradient);
						AveGradient += tempGradient;            
					} //end if
				} // end for i
			} // end for j
		} // end for k

		if(sz>0)
			AveGradient/=sz;
		else 
			AveGradient=400;

		double k_factor = AveGradient;
		double lamda = 0.02;  // 0.02

		if(!iso)	//Do anisotropic diffusion
		{
			for (int k=border; k<sizeZ-border; k++) 
			{
				for (int j=border; j<sizeY-border; j++) 
				{
					for (int i=border; i<sizeX-border; i++)
					{
						if(volin[k*sls + j*sizeX + i]>0) // by xiao
						{
							int idx = k*sls + j*sizeX + i;
							double diverg = 0;
							int iidx1 = k *sls + j *sizeX + i-1;
							int iidx2 = k *sls + j *sizeX + i+1;
							double Dxy1 = exp( -(gradVec[iidx1].xd/k_factor) * (gradVec[iidx1].xd/k_factor) );
							double Dxy2 = exp( -(gradVec[iidx2].xd/k_factor) * (gradVec[iidx2].xd/k_factor) );
							diverg = diverg + (Dxy2 * gradVec[iidx2].xd - Dxy1 * gradVec[iidx1].xd);
				          
							iidx1 = k *sls + (j-1) *sizeX + i;
							iidx2 = k *sls + (j+1) *sizeX + i;
							Dxy1 = exp( -(gradVec[iidx1].yd/k_factor) * (gradVec[iidx1].yd/k_factor) );
							Dxy2 = exp( -(gradVec[iidx2].yd/k_factor) * (gradVec[iidx2].yd/k_factor) );
							diverg = diverg + (Dxy2 * gradVec[iidx2].yd - Dxy1 * gradVec[iidx1].yd);
				          
							iidx1 = (k-1) *sls + j *sizeX + i;
							iidx2 = (k+1) *sls + j *sizeX + i;
							Dxy1 = exp( -(gradVec[iidx1].zd/k_factor) * (gradVec[iidx1].zd/k_factor) );
							Dxy2 = exp( -(gradVec[iidx2].zd/k_factor) * (gradVec[iidx2].zd/k_factor) );
							diverg = diverg + (Dxy2 * gradVec[iidx2].zd - Dxy1 * gradVec[iidx1].zd);

							double voxelUpdate = volin[idx] + lamda * diverg;
							if (voxelUpdate<0)   voxelUpdate=0;
							if (voxelUpdate>255)   voxelUpdate=255;
							volout[idx] = (int)(voxelUpdate+0.5);  // rounding to int
						} // end if
					} // end for i
				} // end for j
			} // end for k
		}
		else		//Do Isotropic diffusion
		{
			for (int k=border; k<sizeZ-border; k++) 
			{
				for (int j=border; j<sizeY-border; j++) 
				{
					for (int i=border; i<sizeX-border; i++)
					{
						if(volin[k*sls + j*sizeX + i]>0) // by xiao
						{
							int idx = k*sls + j*sizeX + i;
							double diverg = 0;
							int iidx1 = k *sls + j *sizeX + i-1;
							int iidx2 = k *sls + j *sizeX + i+1;
							diverg = diverg + (gradVec[iidx2].xd - gradVec[iidx1].xd);
							iidx1 = k *sls + (j-1) *sizeX + i;
							iidx2 = k *sls + (j+1) *sizeX + i;
							diverg = diverg + (gradVec[iidx2].yd - gradVec[iidx1].yd);
							iidx1 = (k-1) *sls + j *sizeX + i;
							iidx2 = (k+1) *sls + j *sizeX + i;
							diverg = diverg + (gradVec[iidx2].zd - gradVec[iidx1].zd);

							double voxelUpdate = volin[idx] + lamda * diverg;
							if (voxelUpdate<0)   voxelUpdate=0;
							if (voxelUpdate>255)   voxelUpdate=255;
							volout[idx] = (int)(voxelUpdate+0.5);  // rounding to int
						} //end if
					} // end for i
				} // end for j
			} // end for k
		} // end if(useIso)

		timesDiffuse--;

		// copy volout back to volin for the next dilation
		for (int idx=0; idx<sz; idx++) 
		{
			volin[idx] = volout[idx];
		}
	} // end while(timesDiffuse)

	/* Removed by Xiao Liang
	// Post-smoothing
	double blockAve;
	int ii, jj, kk;
	for (int k=border; k<sizeZ-border; k++)
		for (int j=border; j<sizeY-border; j++)
			for (int i=border; i<sizeX-border; i++)   
			{
				blockAve = 0;
				for (kk=-1; kk<=1; kk++)
					for (jj=-1; jj<=1; jj++)
						for (ii=-1; ii<=1; ii++) 
						{
							blockAve += volin[(k+kk)*sizeX*sizeY + (j+jj)*sizeX + (i+ii)];
				}
				volout[k *sizeX*sizeY + j *sizeX + i] = blockAve / 27;
	}
	*/ 

	//Set the output ITK image to the buffer
	itk::ImageRegionIterator< ImageType > itr( m_outputImage, m_outputImage->GetBufferedRegion() );
	long idx = 0;
	for(itr.GoToBegin(); !itr.IsAtEnd(); ++itr)
	{
		itr.Set( volout[idx++] );
	}

	//free memory:
	free(volin);
	free(volout);
	free(gradVec);

	if(debug)
	{
		std::cerr << "mdl Anisotropic Diffusion is Done" << std::endl;
	}

	return true;
}


bool VolumeProcess::RunManualThreshold(double threshold)
{
  // this function is corrsponding to the Xiaosong's original MDL code
	if(debug)
	{
		std::cerr << "We excute the Manual Thrsholding by " << threshold << std::endl;
	}
    
	//Apply threshold (any thing below threshold is set to zero)
	itk::ImageRegionIterator< ImageType > itr( m_outputImage, m_outputImage->GetLargestPossibleRegion() );
	for(itr.GoToBegin(); !itr.IsAtEnd(); ++itr)
	{
		if(itr.Get() < threshold) 
        {
			itr.Set(0);
        }
	}
	return true;

}

bool  VolumeProcess::distTransform(unsigned char *f, int L, int M, int N) 
{
  if(debug)
	{
		std::cerr << "We are dosing distance Transform!" << std::endl;
	}
  int i,j,k,n;
 //float *buff , df, db, d, w;
  
  //float  df, db; //d 
  long df,db;
  long *fDist,*buff;
  // modify by xl
  long idx, slsz, sz;
  
  slsz = L*M;		// slice size
  sz = slsz*N;

  fDist = new long[L*M*N];

  for (idx = 0; idx < slsz*N; idx++) {
      //if (f[idx] > 0)   fDist[idx] = 0;
	  //else fDist[idx] = 5000;
	  if (f[idx] > 0)   fDist[idx] = 5000;
	  else fDist[idx] = 0;
  }

  int maxdim = MAX(L,M);
  maxdim = MAX(maxdim,N);
  
    //buff = new float[maxdim+10];
  buff = new long[maxdim+10];

  // Using Algorithm 3 from Appendix 

  // Step 1  forward scan
  
  for (k = 0; k < N; k++)
    for (j = 0; j < M; j++)
    {
       df = L;
       for (i = 0; i < L; i++)
       {
         idx = k*slsz + j*L + i;
         if (fDist[idx] !=0)
           df = df + 1;
         else
           df = 0;
         fDist[idx] = df*df;
       }
     }
 
  //  Step 1 backward scan
  
  for (k = 0; k < N; k++)
    for (j = 0; j < M; j++)
    {
      db = L;
      for (i = L-1; i >=0; i--)
      {
        idx = k*slsz + j*L + i;
        if (fDist[idx] !=0)
          db = db + 1;
        else
          db = 0; 
        fDist[idx] = MIN(fDist[idx],db*db);
      }
    }

  // Step 2
 
  long d,w;  // add by xiao liang

  for (k = 0; k < N; k++)
    for (i = 0; i < L; i++)
    {
      for (j =0; j < M; j++)
        buff[j] = fDist[k*slsz + j*L +i];
    
      for (j = 0; j < M; j++)
      {
        d = buff[j];
        if (d != 0)
        {
          int rmax, rstart, rend;
          rmax = (int) floor(sqrt((double)d)) + 1;
          rstart = MIN(rmax, (j-1));
          rend = MIN(rmax, (M-j));
          for (n = -rstart; n < rend; n++)
          {
              if (j+n >= 0 && j+n < M)
              {
                w = buff[j+n] + n*n;
                if (w < d)  d = w;
              }
          }
        }
        idx = k*slsz + j*L +i;
        fDist[idx] = d;
      }
    }

  // Step 3
  for (j = 0; j < M; j++)
    for (i = 0; i < L; i++)
    {
      for (k =0; k < N; k++)
        buff[k] = fDist[k*slsz + j*L +i];
    
      for (k = 0; k < N; k++)
      {
        d = buff[k];
        if (d != 0)
        {
          int rmax, rstart, rend;
          rmax = (int) floor(sqrt((double)d)) + 1;
          rstart = MIN(rmax, (k-1));
          rend = MIN(rmax, (N-k));
          for (n = -rstart; n < rend; n++)
          {
              if (k+n >= 0 && k+n < N)
              {
                w = buff[k+n] + n*n;
                if (w < d)  d = w;
              }
          }
        }
        idx = k*slsz + j*L +i;
        fDist[idx] = d;
      }
    }

  for (idx = 0; idx < slsz*N; idx++) {
      fDist[idx] = (long) sqrt(float(fDist[idx]));
  }


  double dMax = 0;
  for(idx=0; idx<sz; idx++)   {  // Scale the dist result to 255
	  if (fDist[idx] > dMax)   dMax = fDist[idx];
  }
  for(idx=0; idx<sz; idx++)   {  // Scale the dist result to 255
	   // f[idx] = fDist[idx] * 255/ dMax;
	   f[idx] = (unsigned char) fDist[idx]; // by xiao (long)
  }
  
  delete []buff;
  delete []fDist;
  return true;
}


bool VolumeProcess::RunDistanceTransform(void)
{
	ImageType::RegionType region = m_outputImage->GetBufferedRegion();
	int numColumns = region.GetSize(0);
	int numRows = region.GetSize(1);
	int numStacks = region.GetSize(2);
	long numPix = numStacks*numColumns*numRows;

	PixelType * dataImagePtr = m_outputImage->GetBufferPointer();
	unsigned char * binImagePtr = new unsigned char[numPix];

	itk::ImageRegionIterator< ImageType > itr( m_outputImage, m_outputImage->GetLargestPossibleRegion() );

    long idx = 0;
	for(itr.GoToBegin(); !itr.IsAtEnd(); ++itr)
	{
		binImagePtr[idx++]=itr.Get();
	}

	bool ok= distTransform(binImagePtr, numRows, numColumns, numStacks);	//Do distance transform
	if(!ok)
		return false;

	if(debug)
	{
		std::cerr << " Distance Transform is done and we will replace image with the value of DT!" << std::endl;
	}
	//repalce the image by the distance transform
	
	idx = 0;
	for(itr.GoToBegin(); !itr.IsAtEnd(); ++itr)
	{
		itr.Set( binImagePtr[idx++] );
	}

	delete [] binImagePtr;
	return true;
}

bool VolumeProcess:: RunDanielssonDistanceMap(void)
{
    typedef itk::DanielssonDistanceMapImageFilter<ImageType, FloatImageType3D>  DT_Type;
	DT_Type::Pointer DTfilter = DT_Type::New();
	DTfilter->SetInput( m_outputImage );
	//DTfilter->GetDistanceMap();
	try
	{
		DTfilter->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "ITK FILTER ERROR: " << err << std::endl ;
		return false;
	}
	m_outputImage = RescaleFloatToImageType(DTfilter->GetDistanceMap());
	
	if(debug)
		std::cerr << "DanielssonDistanceMap Filter Done" << std::endl;
	return true;
}

bool VolumeProcess::RunVotingBinaryHoleFilling(int radiusX, int radiusY, int radiusZ)
{
    if (radiusX == 0 && radiusY ==0 && radiusZ ==0)
		return true;

	typedef itk::VotingBinaryHoleFillingImageFilter<ImageType, ImageType>  VotingBinaryHoleFillingFilterType;
	VotingBinaryHoleFillingFilterType::Pointer VotingBinaryHoleFillingFilter = VotingBinaryHoleFillingFilterType::New();
	ImageType::SizeType indexRadius;
	indexRadius[0] = radiusX; 
	indexRadius[1] = radiusY;
    indexRadius[2] = radiusZ;
	VotingBinaryHoleFillingFilter->SetRadius(indexRadius);
	VotingBinaryHoleFillingFilter->SetMajorityThreshold(2);
	VotingBinaryHoleFillingFilter->SetInput(m_outputImage);
    try
	{
		VotingBinaryHoleFillingFilter->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "ITK FILTER ERROR: " << err << std::endl ;
		return false;
	}
	m_outputImage = VotingBinaryHoleFillingFilter->GetOutput();
    if(debug)
		std::cerr << "RunVotingBinaryHoleFilling Filter Done" << std::endl;
	return true;
}


bool VolumeProcess::RunVotingBinaryIterativeHoleFilling(int radiusX, int radiusY, int radiusZ, int iterativeNumber)
{
    if (radiusX == 0 && radiusY ==0 && radiusZ ==0)
		return true;

	typedef itk::VotingBinaryIterativeHoleFillingImageFilter<ImageType>  VotingBinaryHoleFillingFilterType;
	VotingBinaryHoleFillingFilterType::Pointer VotingBinaryHoleFillingFilter = VotingBinaryHoleFillingFilterType::New();
	ImageType::SizeType indexRadius;
	indexRadius[0] = radiusX; 
	indexRadius[1] = radiusY;
    indexRadius[2] = radiusZ;
	VotingBinaryHoleFillingFilter->SetRadius(indexRadius);
	VotingBinaryHoleFillingFilter->SetMajorityThreshold(2);
	VotingBinaryHoleFillingFilter->SetInput(m_outputImage);
	VotingBinaryHoleFillingFilter->SetMaximumNumberOfIterations(iterativeNumber);
    try
	{
		VotingBinaryHoleFillingFilter->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "ITK FILTER ERROR: " << err << std::endl ;
		return false;
	}
	m_outputImage = VotingBinaryHoleFillingFilter->GetOutput();
    if(debug)
		std::cerr << "RunVotingBinaryHoleFilling Filter Done" << std::endl;
	return true;

}

bool VolumeProcess::RunBinayMedianHoleFilling(int radiusX,int radiusY, int radiusZ)
{

	if (radiusX == 0 && radiusY ==0 && radiusZ ==0)
		return true;

	typedef itk::BinaryMedianImageFilter<ImageType, ImageType>  BinaryMedianImageFilterType;
	BinaryMedianImageFilterType::Pointer BinaryMedianImageFilter = BinaryMedianImageFilterType::New();
	ImageType::SizeType indexRadius;
	indexRadius[0] = radiusX; 
	indexRadius[1] = radiusY;
    indexRadius[2] = radiusZ;
	BinaryMedianImageFilter->SetRadius(indexRadius);
	BinaryMedianImageFilter->SetInput(m_outputImage);
    try
	{
		BinaryMedianImageFilter->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "ITK FILTER ERROR: " << err << std::endl ;
		return false;
	}
	m_outputImage = BinaryMedianImageFilter->GetOutput();
    if(debug)
		std::cerr << "RunBinaryMedianImageFilter Done" << std::endl;
	return true;

}

bool VolumeProcess::NonlinearMappingSigmoidFilter(double alpha,  double beta, double min, double max)
{   // this function is a non-linear Mapping, by Xiao L.
    typedef itk::SigmoidImageFilter<ImageType, ImageType>  SigmoidImageFilterType;
    SigmoidImageFilterType::Pointer SigmoidImageFilter = SigmoidImageFilterType::New();
    SigmoidImageFilter->SetInput( m_outputImage );
	SigmoidImageFilter->SetAlpha(alpha);
	SigmoidImageFilter->SetBeta(beta);
	SigmoidImageFilter->SetOutputMinimum( min );
    SigmoidImageFilter->SetOutputMaximum( max );
	try
	{
		SigmoidImageFilter->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "ITK FILTER ERROR: " << err << std::endl ;
		return false;
	}
	m_outputImage = SigmoidImageFilter->GetOutput();
	if(debug)
		std::cerr << "Rescale Filter Done" << std::endl;
	return true;
}


}



