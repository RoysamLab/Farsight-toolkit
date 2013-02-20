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

/**
 \brief Main function for tracing in 3D volume. The input image can be
 \author $ Author: Amit Mukherjee, James Alex Tyrrell $
 \version $Revision 1.0 $
 \date May 05 2009
*/
// Copyright Farsight Toolkit 2009, Rensselaer Polytechnic institute Troy NY 12180.


#include <iostream>
#include <fstream>

#include "itkImage.h"
#include "itkRGBPixel.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkMedianImageFilter.h"
#include "itkStatisticsImageFilter.h"

#include "itkArray.h"
#include "itkCovariantVector.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkShapedNeighborhoodIterator.h"
#include "itkSymmetricEigenAnalysis.h"
#include "vnl/vnl_math.h"
#include <vnl/vnl_vector_fixed.h>

#include "TraceConfig.h"
#include "SeedContainer3D.h"
#include "Seed2Seg.h"
#include "TraceContainer3D.h"



typedef	 float PixelType ;
typedef itk::Image< PixelType, 2 >   ImageType2D;
typedef itk::Image< PixelType, 3 >   ImageType3D;
//typedef   itk::HessianRecursiveGaussianImageFilter<  ImageType3D > HessianFilterType;
//typedef   HessianFilterType::OutputImageType            HessianImageType;
//typedef   HessianImageType::PixelType                   HessianPixelType;


typedef vnl_vector_fixed<double,3> Vect3;

void ImageDenoise(ImageType3D::Pointer&, int );
void UpdateMultiscale( ImageType3D::Pointer& , ImageType3D::Pointer& );
//void GetFeature( ImageType3D::Pointer&, ImageType3D::Pointer& , const float);
void GetFeature( ImageType3D::Pointer&, ImageType3D::Pointer&);
void ImageStatistics(ImageType3D::Pointer & );
static void WriteSeedImage(ImageType2D::Pointer, SeedContainer3D::Pointer,std::string );
static void WriteStartSegmentsImage(ImageType2D::Pointer, Seed2Seg::Pointer ,std::string );


int main (int argc, char *argv[])
{

	std::cout << "SuperEllipsoid Trace3D version-0.1\tSept 11, 2007" << std::endl << std::endl;
	if (argc != 2)	{
		std::cout << "Usage: "<<argv[0] << " [InputParameterFilename.xml]" <<std::endl;
		exit(1);
	}

	std::cout << "Reading "<< argv[1] << std::endl;
	TraceConfig::Pointer m_Config = TraceConfig::New();
	if ( !m_Config->LoadParameters(argv[1]) ) {
		std::cout << "Problem encountered in reading "<< argv[1] << std::endl;
		return EXIT_FAILURE;
	}
	std::cout <<  std::endl <<  std::endl<< 
		"Parameter file parsed successfully, " <<m_Config->getNumberOfDataFiles() <<" images in the processing list!!" <<  std::endl<<  std::endl;


	for (unsigned int k = 0; k < m_Config->getNumberOfDataFiles(); k++) {
		std::cout << "Tracing " << k+1 <<" out of " << m_Config->getNumberOfDataFiles() 
			<< " image: " << m_Config->getInputFileName(k) << std::endl<< std::endl;

		ImageType3D::Pointer image3D = ImageType3D::New();
		ImageType2D::Pointer MIPimage = ImageType2D::New();
		typedef itk::ImageFileReader<ImageType3D> ReaderType;
		ReaderType::GlobalWarningDisplayOff();
		ReaderType::Pointer reader = ReaderType::New();

		reader->SetFileName(m_Config->getInputFileName(k));
		image3D = (reader->GetOutput());
		try {
			image3D->Update();
		}
		catch (itk::ExceptionObject &e)	{
			std::cout << "Exception caught in opening input image file!!! " << std::endl << e << std::endl;
			return EXIT_FAILURE;
		}
		std::cout << "Image of size " << image3D->GetBufferedRegion().GetSize() << " read successfully " << std::endl;

		ImageDenoise(image3D, m_Config->getHessianFlag());
		ImageStatistics(image3D);

		try	{
			SeedContainer3D::Pointer m_Seeds = SeedContainer3D::New();
			m_Seeds->Configure(m_Config);
			m_Seeds->Detect(image3D, MIPimage);
			WriteSeedImage(MIPimage, m_Seeds , m_Config->getOutputFileName(k));
	
			Seed2Seg::Pointer m_SS = Seed2Seg::New();
			m_SS->Configure(m_Config);
			m_SS->ComuputeStartSegments(m_Seeds , image3D, m_Config);
			m_SS->SortStartSegments();
			WriteStartSegmentsImage(MIPimage, m_SS , m_Config->getOutputFileName(k));

			TraceContainer3D::Pointer m_Tracer = TraceContainer3D::New();
			m_Tracer->Configure(m_Config);
			m_Tracer->ComputeTrace(image3D, m_SS) ;
			//m_Tracer->WriteTraceToTxtFile(m_Config->getOutputFileName(k));
			m_Tracer->WriteTraceToXMLFile(m_Config->getOutputFileName(k));

			m_Seeds = NULL;
			m_SS = NULL;
			m_Tracer = NULL;
			image3D = NULL;
			MIPimage = NULL;
		}
		catch (itk::ExceptionObject & err )	    {
			std::cout << "ExceptionObject caught !" << std::endl;
			std::cout << err << std::endl;
			return EXIT_FAILURE;
		}
	}

	return EXIT_SUCCESS;
}

//void ImageDenoise(ImageType3D::Pointer& vol)	{
void ImageDenoise(ImageType3D::Pointer& im3D, int hessianFlag)	{
	// use median filtering based denoising to get rid of impulsive noise
	std::cout << "Denoising Image..." << std::endl;
	typedef itk::MedianImageFilter<ImageType3D, ImageType3D> MedianFilterType;
	MedianFilterType::Pointer medfilt = MedianFilterType::New();
	medfilt->SetInput(im3D);
	ImageType3D::SizeType rad = {1, 1, 1};
	medfilt->SetRadius(rad);

	ImageType3D::Pointer vol = medfilt->GetOutput();
	vol->Update();

	if (hessianFlag == 1)	{
 		ImageType3D::PixelType sigmas[] = { 2.0f, 2.88f, 4.0f, 5.68f, 8.0f };
		ImageType3D::Pointer temp = ImageType3D::New();
		ImageType3D::Pointer maxF = ImageType3D::New();

		temp->SetRegions(vol->GetLargestPossibleRegion());
		temp->Allocate();
		temp->FillBuffer(0.0);
		maxF->SetRegions(vol->GetLargestPossibleRegion());
		maxF->Allocate();
		maxF->FillBuffer(0.0);

		
		typedef itk::SmoothingRecursiveGaussianImageFilter< ImageType3D , ImageType3D> GFilterType;
		GFilterType::Pointer gauss = GFilterType::New();

		for (unsigned int i = 0; i < 4; ++i)	{
			std::cout << "Performing 3D Line Filtering using Hessian at " << sigmas[i] ;
			gauss->SetInput( vol );
			gauss->SetSigma( sigmas[i] );
			gauss->SetNormalizeAcrossScale(true);
			ImageType3D::Pointer svol = gauss->GetOutput();
			svol->Update();

			std::cout << "...";
			//GetFeature( vol, temp, sigmas[i] );
			GetFeature( temp, svol );
			std::cout << "....done" << std::endl;
			UpdateMultiscale( temp, maxF);
		}
		im3D = maxF;
	}
	else {
		im3D = vol;
	}
}

void ImageStatistics(ImageType3D::Pointer & im3D)	{
	typedef itk::StatisticsImageFilter<ImageType3D>  StatisticsType;
	typedef itk::ImageRegionIterator<ImageType3D> IteratorType;

	StatisticsType::Pointer statistics = StatisticsType::New();
	statistics->SetInput(im3D );
	statistics->Update();
	float imin = statistics->GetMinimum();
	float imax = statistics->GetMaximum();
	float imean = statistics->GetMean();
	float istd = statistics->GetSigma();

	std::cout << "Input Image Statistics: min:"<< imin << " max:" << imax << " mean:"<< imean << " std:" << istd << std::endl;
	if ((imean - imin) < (imax - imean))	{
		std::cout << "Inverting intensities" << std::endl;
		IteratorType it(im3D, im3D->GetRequestedRegion());
		for (it.GoToBegin(); !it.IsAtEnd(); ++it)	{
			float d = it.Get();
			float imax2 = imean+4*istd;
			d = (d > (imax2)) ? (imax2) : d;
			d = 1 - (d - imin)/(imax2 - imin);
			it.Set(255.0*d);
		}
	}
	else {
		IteratorType it(im3D, im3D->GetRequestedRegion());
			for (it.GoToBegin(); !it.IsAtEnd(); ++it)	{
			float d = it.Get();
			float imin2 = imean-4*istd;
			d = (d < (imin2)) ? (imin2) : d;
			d = (d - imin2)/(imax - imin2);
			it.Set(255.0*d);
		}
	}
}



static void WriteSeedImage(ImageType2D::Pointer im, SeedContainer3D::Pointer seedCnt ,std::string filename)	{

  typedef itk::RGBPixel<unsigned char>   RGBPixelType;
  typedef itk::Image< RGBPixelType,  2 >    RGBImageType;
  typedef itk::ImageFileWriter< RGBImageType >  RGBWriterType;

  RGBImageType::Pointer RGBImage = RGBImageType::New();
  RGBImageType::IndexType start;
  start[0] = 0;
  start[1] = 0;

  RGBImageType::SizeType size = im->GetRequestedRegion().GetSize();;

  RGBImageType::RegionType region;
  region.SetSize( size );
  region.SetIndex( start );

  RGBImage->SetRegions( region );
  RGBImage->Allocate();


  itk::ImageRegionIterator<RGBImageType> iterRGB(RGBImage,RGBImage->GetBufferedRegion());
  itk::ImageRegionConstIterator<ImageType2D> iter(im,im->GetBufferedRegion());

  RGBPixelType newPixel;
  for (iter.GoToBegin(), iterRGB.GoToBegin(); !iter.IsAtEnd(); ++iter, ++iterRGB)	{

		newPixel.SetRed(static_cast<unsigned char>(iter.Get()));
		newPixel.SetGreen(static_cast<unsigned char>(iter.Get()));
	    newPixel.SetBlue(static_cast<unsigned char>(iter.Get()));
	  	iterRGB.Set(newPixel);

  }


  RGBPixelType pixelValue;
  RGBImageType::IndexType pixelIndex;

  std::cout << "Writing "<< seedCnt->getNumberOfSeeds() << " seeds";

  for (unsigned int i = 0; i< seedCnt->getNumberOfSeeds() ;i++){

	  Vect3 pos = seedCnt->getSeed(i)->getPosition();
	  pixelIndex[0] = static_cast<long int> (pos(0));
      pixelIndex[1] = static_cast<long int> (pos(1));

      pixelValue[2] = 0;
      pixelValue[1] = 255;
      pixelValue[0] = 0;

      RGBImage->SetPixel(pixelIndex, pixelValue);
  }

  //Printing the ith seed for verification

  RGBWriterType::Pointer writer = RGBWriterType::New();
  //writer->SetFileName( "Seed_Points.png" );
  writer->SetFileName( filename + "_Seed_Points.png" );

  writer->SetInput( RGBImage );
  try
    {
      writer->Update();
      std::cout <<"...done." << std::endl;
    }
  catch (itk::ExceptionObject & e)
    {
      std::cerr << "Exception in Writer: " << e << std::endl;
      exit(0);
    }
}

static void WriteStartSegmentsImage(ImageType2D::Pointer im, Seed2Seg::Pointer sseg ,std::string filename)	{

  typedef itk::RGBPixel<unsigned char>   RGBPixelType;
  typedef itk::Image< RGBPixelType,  2 >    RGBImageType;
  typedef itk::ImageFileWriter< RGBImageType >  RGBWriterType;

  RGBImageType::Pointer RGBImage = RGBImageType::New();
  RGBImageType::IndexType start;
  start[0] = 0;
  start[1] = 0;

  RGBImageType::SizeType size = im->GetRequestedRegion().GetSize();;

  RGBImageType::RegionType region;
  region.SetSize( size );
  region.SetIndex( start );

  RGBImage->SetRegions( region );
  RGBImage->Allocate();


  itk::ImageRegionIterator<RGBImageType> iterRGB(RGBImage,RGBImage->GetBufferedRegion());
  itk::ImageRegionConstIterator<ImageType2D> iter(im,im->GetBufferedRegion());

  RGBPixelType newPixel;
  for (iter.GoToBegin(), iterRGB.GoToBegin(); !iter.IsAtEnd(); ++iter, ++iterRGB)	{

		newPixel.SetRed(static_cast<unsigned char>(iter.Get()));
		newPixel.SetGreen(static_cast<unsigned char>(iter.Get()));
	    newPixel.SetBlue(static_cast<unsigned char>(iter.Get()));
	  	iterRGB.Set(newPixel);

  }


  RGBPixelType pixelValue;
  RGBImageType::IndexType pixelIndex;

  std::cout << "Writing "<< sseg->getNumberOfStartSegments() << " Start Segments ";

  for (unsigned int i=0; i< sseg->getNumberOfStartSegments(); i++ )	{
	
	  TVessel *seg = sseg->getStartSegment(i);
	  Vect3 pos(seg->mu);
	  pixelIndex[0] = static_cast<long int> (pos(0));
      pixelIndex[1] = static_cast<long int> (pos(1));

      pixelValue[2] = 0;
      pixelValue[1] = 0;
      pixelValue[0] = 255;

      RGBImage->SetPixel(pixelIndex, pixelValue);
  }

  //Printing the ith seed for verification

  RGBWriterType::Pointer writer = RGBWriterType::New();
  //writer->SetFileName( "Seed_Points.png" );
  writer->SetFileName( filename + "_StartSegments.png" );

  writer->SetInput( RGBImage );
  try
    {
      writer->Update();
      std::cout <<"...done." << std::endl;
    }
  catch (itk::ExceptionObject & e)
    {
      std::cerr << "Exception in Writer: " << e << std::endl;
      exit(0);
    }
}
void UpdateMultiscale( ImageType3D::Pointer& temp, ImageType3D::Pointer& maxF)	{
	itk::ImageRegionIterator<ImageType3D> itt(temp, temp->GetBufferedRegion());
	itk::ImageRegionIterator<ImageType3D> itm(maxF, maxF->GetBufferedRegion());
	for (itt.GoToBegin(), itm.GoToBegin(); !itt.IsAtEnd(); ++itt, ++itm)	{
		itm.Set(vnl_math_max(itm.Get(), itt.Get()));
	}
}


void GetFeature( ImageType3D::Pointer& temp, ImageType3D::Pointer& svol)	{
		// set the diagonal terms in neighborhood iterator
	itk::Offset<3> 
		xp =  {2 ,  0 ,   0}, 
		xn =  {-2,  0,    0},
		yp =  {0,   2,	  0},
		yn =  {0,  -2,    0}, 
		zp =  {0,   0,    2},
		zn =  {0,   0,   -2}, 
		center = {0, 0 , 0};

	itk::Size<3> rad = {{1,1,1}};
	itk::NeighborhoodIterator<ImageType3D> nit(rad , svol, svol->GetBufferedRegion());
	itk::ImageRegionIterator<ImageType3D> it(svol, svol->GetBufferedRegion());
	itk::ImageRegionIterator<ImageType3D> itt(temp, temp->GetBufferedRegion());


	// set the offsets for Hessian computation
	// 6 7 8    15 16 17   24 25 26
	// 3 4 5    12 13 14   21 22 23 
	// 0 1 2    9  10 11   18 19 20

	unsigned int 
		xy1 =  17, //{ 1 ,   1 ,	0 }, 
		xy2 =  9,  //{ -1,  -1 ,	0 },
		xy3 =  15, //{ -1,   1 ,  0 },
		xy4 =  11, //{ 1 ,  -1 ,  0 }, 

		yz1 =  25, //{ 0 ,   1 ,  1 },
		yz2 =  1,  //{ 0 ,  -1 , -1 },
		yz3 =  19, //{ 0 ,  -1 ,  1 },
		yz4 =  7,  //{ 0 ,   1 , -1 }, 

		xz1 =  23, //{ 1 ,   0 ,  1 }, 
		xz2 =  3,  //{-1 ,   0 , -1 },
		xz3 =  21, //{-1 ,   0 ,  1 },
		xz4 =  5;  //{ 1 ,   0 , -1 };

	// Bright tubular structures will have low lambda_1 and large negative
	// values of lambda_2 and lambda_3.
	typedef itk::FixedArray< double, 3 > EigenValuesArrayType;
	typedef itk::Matrix< double, 3, 3 > EigenVectorMatrixType;
	typedef itk::SymmetricSecondRankTensor<double,3> TensorType;

	itk::Size<3> sz = svol->GetBufferedRegion().GetSize();
	sz[0] = sz[0] - 3; sz[1] = sz[1] - 3; sz[2] = sz[2] - 3;

	it.GoToBegin();
	nit.GoToBegin();
	itt.GoToBegin();
	itk::Vector<float,3> sp = svol->GetSpacing();

	float alpha1 = 0.5;
	float alpha2 = 2;

	while(!nit.IsAtEnd())	{
		itk::Index<3> ndx = it.GetIndex();
		if ( (ndx[0]<2) || (ndx[1]<2) || (ndx[2]<2) || (ndx[0]>sz[0]) || (ndx[1]>sz[1]) || (ndx[2]>sz[2]) )	{
			++itt;
			++it;
			++nit;
			continue;
		}

		TensorType h;
		h.Fill(0.0);
		h[0] = svol->GetPixel( ndx + xp ) + svol->GetPixel( ndx + xn ) - 2*nit.GetPixel( 13 );
		h[3] = svol->GetPixel( ndx + yp ) + svol->GetPixel( ndx + yn ) - 2*nit.GetPixel( 13 );
		h[5] = svol->GetPixel( ndx + zp ) + svol->GetPixel( ndx + zn ) - 2*nit.GetPixel( 13 );

		float p = 0.0f;

		if ( (h[0]+h[3]+h[5]) < 0)	{
			h[1] = nit.GetPixel(xy1) + nit.GetPixel(xy2) - nit.GetPixel(xy3) - nit.GetPixel(xy4);
			h[2] = nit.GetPixel(xz1) + nit.GetPixel(xz2) - nit.GetPixel(xz3) - nit.GetPixel(xz4);
			h[4] = nit.GetPixel(yz1) + nit.GetPixel(yz2) - nit.GetPixel(yz3) - nit.GetPixel(yz4);

			EigenValuesArrayType ev;
			EigenVectorMatrixType em;
			h.ComputeEigenAnalysis (ev, em);
			
			float temp;
			if(ev[0] > ev[1])	{
				temp = ev[0];
				ev[0] = ev[1];
				ev[1] = temp;
			}
			if(ev[1] > ev[2])	{
				temp = ev[1];
				ev[1] = ev[2];
				ev[2] = temp;
			}

			//Assign vesselness
			if (ev[1] < 0)	{
				//use -ev[1] as normalization, and ev[2] as the numerator
				float norm = 0;
				if(ev[2] > 0)	{
					norm = vnl_math_sqr(ev[2]/(alpha1*ev[1]));
				}
				else	{
					norm = vnl_math_sqr(ev[2]/(alpha2*ev[1]));
				}
				p = -1*ev[1]*vcl_exp(-0.5*norm);
			}
			else	{
				p = 0.0;
			}
		}
		itt.Set(p);
		++itt;
		++it;
		++nit;
	}
}

/*
void GetFeature( ImageType3D::Pointer& temp, HessianImageType::Pointer& hessian) {

	//typedef   itk::HessianRecursiveGaussianImageFilter<  ImageType3D > HessianFilterType;
	//typedef   HessianFilterType::OutputImageType            HessianImageType;
	//typedef   HessianImageType::PixelType                   HessianPixelType;

	//HessianFilterType::Pointer m_Hessian = HessianFilterType::New();
	//m_Hessian->SetInput( vol );
	//m_Hessian->SetSigma( vstd );
	//m_Hessian->SetNormalizeAcrossScale(true);

	

	// Use smart eigen value computation to get the output image
	typedef itk::ImageRegionConstIterator<HessianImageType> HessianIteratorType;
	typedef itk::ImageRegionIterator<ImageType3D> IteratorType;
	IteratorType it(temp, temp->GetRequestedRegion());
	HessianIteratorType ith(hessian, hessian->GetRequestedRegion());

	// Bright tubular structures will have low lambda_1 and large negative
	// values of lambda_2 and lambda_3.

	//typedef itk::SymmetricSecondRankTensor< PixelType, 3 > InputMatrixType;
	typedef itk::FixedArray< PixelType, 3 > EigenValuesArrayType;
	typedef itk::Matrix< PixelType, 3, 3 > EigenVectorMatrixType;
	typedef itk::SymmetricEigenAnalysis< HessianPixelType,EigenValuesArrayType, EigenVectorMatrixType > SymmetricEigenAnalysisType;

	EigenValuesArrayType ev;
	SymmetricEigenAnalysisType symmetricEigenSystem(3);

	it.GoToBegin();
	ith.GoToBegin();


	while(!it.IsAtEnd())	{
		// find the hessian response here
		float p, Tr;
		float alpha1 = 0.5;
		float alpha2 = 2;

		Tr = ith.Get()(0,0) + ith.Get()(1,1) + ith.Get()(2,2);
		if (Tr < 0 ) {
			symmetricEigenSystem.ComputeEigenValues(ith.Get(), ev);
			// 2 lambdas are high -ve (H1, H2) and one low positive L
			//ev[0], ev[1] are H1,H2  ev[2] is L
			float temp;
			if(ev[0] > ev[1])	{
				temp = ev[0];
				ev[0] = ev[1];
				ev[1] = temp;
			}
			if(ev[1] > ev[2])	{
				temp = ev[1];
				ev[1] = ev[2];
				ev[2] = temp;
			}

			//Assign vesselness
			if (ev[1] < 0)	{
				//use -ev[1] as normalization, and ev[2] as the numerator
				float norm = 0;
				if(ev[2] > 0)	{
					norm = vnl_math_sqr(ev[2]/(alpha1*ev[1]));
				}
				else	{
					norm = vnl_math_sqr(ev[2]/(alpha2*ev[1]));
				}
				p = -1*ev[1]*vcl_exp(-0.5*norm);
			}
			else	{
				p = 0.0;
			}
		}
		else {
			p = 0.0;
		}

		it.Set(p);
		++it;
		++ith;
	}
}
*/