/////////////////////////////////////////////////////////////////////////////////
// FILE: TraceMain.cpp
// This file contains the main function.
//

#include <iostream>
#include <fstream>

//#include <float.h>
#include "itkImage.h"
#include "itkRGBPixel.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkStatisticsImageFilter.h"
#include "vnl/vnl_vector_fixed.h"

#include "TraceConfig.h"
#include "SeedContainer3D.h"
#include "Seed2Seg.h"
#include "TraceContainer3D.h"

typedef itk::Image< float, 2 >   ImageType2D;
typedef itk::Image< float, 3 >   ImageType3D;

typedef vnl_vector_fixed<double,3> Vect3;


void ImageStatistics(ImageType3D::Pointer & );
static void WriteSeedImage(ImageType2D::Pointer, SeedContainer3D::Pointer,std::string );
static void WriteSeedImage2(ImageType2D::Pointer, Seed2Seg::Pointer,std::string );


// Entry point for Trace3d executable.
//
// command line:
// Trace3d image.tif [grid_spacing:10]
//
int main (int argc, char *argv[])
{

	std::cout << "SuperEllipsoid Trace3D version-0.1\tSept 11, 2007" << std::endl << std::endl;
	if (argc < 2)	{
		std::cout << "Usage: "<<argv[0] << " [ImageFile] [Grid spacing] " <<std::endl;
		exit(1);
	}

	std::cout << "Tracing "<< argv[1] << std::endl;
	TraceConfig::Pointer m_Config = TraceConfig::New();
	m_Config->SetFileNames(argv[1]);
	if (argc>2)	{
		m_Config->SetGridSpacing(argv[2]);
	}

	ImageType3D::Pointer image3D = ImageType3D::New();
	ImageType2D::Pointer MIPimage = ImageType2D::New();

	typedef itk::ImageFileReader<ImageType3D> ReaderType;
	ReaderType::GlobalWarningDisplayOff();
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(m_Config->getInputFileName());
	image3D = (reader->GetOutput());
	image3D->Update();

	std::cout << "Image of size " << image3D->GetBufferedRegion().GetSize() << " read successfully " << std::endl;

	ImageStatistics(image3D);

	//Trace
	try	{

		//1. Detect seeds
		SeedContainer3D::Pointer m_Seeds = SeedContainer3D::New();
		m_Seeds->Configure(m_Config);
		m_Seeds->Detect(image3D, MIPimage);
		WriteSeedImage(MIPimage, m_Seeds , m_Config->getOutputFileName());

		//2. Fit and prioritize seeds
		Seed2Seg::Pointer m_SS = Seed2Seg::New();
		m_SS->ComuputeStartSegments(m_Seeds , image3D, m_Config);
		m_SS->SortStartSegments();
		WriteSeedImage2(MIPimage, m_SS , m_Config->getOutputFileName());

		//3. Trace from seeds
        TraceContainer3D::Pointer m_Tracer = TraceContainer3D::New();
        m_Tracer->ComputeTrace(image3D, m_SS) ;
		m_Tracer->WriteTraceToTxtFile(m_Config->getOutputFileName());

	}
	catch (itk::ExceptionObject & err )	    {
	    std::cout << "ExceptionObject caught !" << std::endl;
	    std::cout << err << std::endl;
	    return EXIT_FAILURE;
    }

	return EXIT_SUCCESS;
}

//Basic image statistics
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

}

//for debugging and visualization, show initial seed (actually interest) points
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

  //for (unsigned int i = 1; i< 4 ;i++) {
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


//output validated seed points, i.e., interest points that have been fitted by model
static void WriteSeedImage2(ImageType2D::Pointer im, Seed2Seg::Pointer m_SS,std::string filename)	{

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

for (unsigned int i=0; i< m_SS->getNumberOfStartSegments(); i++ )	{

	 TVessel *seg = m_SS->getStartSegment(i)->CopyToNewSegment();

	  pixelIndex[0] = static_cast<long int> (seg->mu[0]);
      pixelIndex[1] = static_cast<long int> (seg->mu[1]);

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

