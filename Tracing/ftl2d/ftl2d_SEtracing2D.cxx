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

/** @file SEtracing2D.cpp
*   @brief The main program of the tracing algorithm using Superellipsoid.
*
*   @author Amit Mukherjee
*//////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <float.h>

#include <itkImage.h>
#include <itkSmartPointer.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include "ftl2d_TraceConfig.h"
#include "ftl2d_SeedContainer2D.h"
#include "ftl2d_Tracer.h"


typedef itk::Image <float, 2> ImageType;
typedef itk::ImageFileReader< ImageType >  ReaderType;

//static void WriteSeedImage(ImageType::Pointer&, SeedContainer2D *,std::string );
static void WriteTraceImage(ImageType::Pointer&, Tracer *,std::string );
static void WriteTraceImageBinary(ImageType::Pointer&, Tracer *,std::string );


int main (int argc, char *argv[])
{
	std::cout << "SuperEllipsoid Trace2D version-0.1\tSept 11, 2007" << std::endl << std::endl;

	if (argc < 2)	{
		std::cout << "Usage: "<<argv[0] << " [parameter_file]" <<std::endl;
		exit(1);
	}
	TraceConfig * myConfig = new TraceConfig();
  	myConfig->LoadParameters( argv[1]);

  	ReaderType::Pointer reader = ReaderType::New();
  	reader->SetFileName( myConfig->getInputFileName());
  	ImageType::Pointer image = reader->GetOutput();
  	image->Update();


	try
	{
		SeedContainer2D * mySeed = new SeedContainer2D(myConfig);
		mySeed->Detect(image);
	//	mySeed->Verify(image);
		mySeed->SortSeeds();
	//	WriteSeedImage(image, mySeed, myConfig->getOutputFileName());

        Tracer * myTracer = new Tracer(myConfig);
        myTracer->Run(image, mySeed, myConfig);
	//	myTracer->ApplyRules();
	//	myTracer->ReportXML(myConfig);
	//	WriteTraceImage(image, myTracer, myConfig->getOutputFileName());
		WriteTraceImageBinary(image, myTracer, myConfig->getOutputFileName());

		delete myTracer; myTracer = NULL;
		delete mySeed; mySeed = NULL;
	}

	catch( itk::ExceptionObject & err )
	    {
	    std::cout << "ExceptionObject caught !" << std::endl;
	    std::cout << err << std::endl;
	    return EXIT_FAILURE;
    }

	delete myConfig; myConfig = NULL;

	return EXIT_SUCCESS;
}


static void WriteSeedImage(ImageType::Pointer &im, SeedContainer2D *mySeeds, std::string fname)
{

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
  std::cout<<"Displaying seeds of Image Size " << size <<std::endl;


  itk::ImageRegionIterator<RGBImageType> iterRGB(RGBImage,RGBImage->GetBufferedRegion());
  itk::ImageRegionConstIterator<ImageType> iter(im,im->GetBufferedRegion());

  RGBPixelType newPixel;
  for (iter.GoToBegin(), iterRGB.GoToBegin(); !iter.IsAtEnd(); ++iter, ++iterRGB)	{

		newPixel.SetRed(static_cast<unsigned char>(iter.Get()/2));
		newPixel.SetGreen(static_cast<unsigned char>(iter.Get()/2));
	    newPixel.SetBlue(static_cast<unsigned char>(iter.Get()/2));
	  	iterRGB.Set(newPixel);

  }


  RGBPixelType pixelValue;
  RGBImageType::IndexType pixelIndex, m;

  std::cout << "Writing "<< mySeeds->getNumberOfSeeds() << " seeds"<< std::endl;

  //for (unsigned int i = 1; i< 4 ;i++) {
  for (unsigned int i = 0; i< mySeeds->getNumberOfSeeds()/4 ;i++){

	//  mySeeds->getSeed(i)->PrintSelf();

	  pixelIndex[0] = mySeeds->getSeed(i)->getx();
      pixelIndex[1] = mySeeds->getSeed(i)->gety();

      pixelValue[2] = 0;
      pixelValue[1] = 255;
      pixelValue[0] = 0;

      RGBImage->SetPixel(pixelIndex, pixelValue);
      //std::cout << pixelIndex << " - "<<pixelValue <<std::endl;
  }

  //Printing the ith seed for verification
/*  mySeeds->getSeed(1)->PrintSelf();
  mySeeds->getSeed(10)->PrintSelf();
  mySeeds->getSeed(20)->PrintSelf();
*/
  RGBWriterType::Pointer writer = RGBWriterType::New();
  //writer->SetFileName( "Seed_Points.png" );
  writer->SetFileName( fname + "_Seed_Points.png" );

  writer->SetInput( RGBImage );
  try
    {
      writer->Update();
    }
  catch (itk::ExceptionObject & e)
    {
      std::cerr << "Exception in Writer: " << e << std::endl;
      exit(0);
    }
}


static void WriteTraceImage(ImageType::Pointer &im, Tracer *tracer, std::string fname)
{

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
  std::cout<<"Displaying seeds of Image Size " << size <<std::endl;


  itk::ImageRegionIterator<RGBImageType> iterRGB(RGBImage,RGBImage->GetBufferedRegion());
  itk::ImageRegionConstIterator<ImageType> iter(im,im->GetBufferedRegion());

  RGBPixelType newPixel;
  for (iter.GoToBegin(), iterRGB.GoToBegin(); !iter.IsAtEnd(); ++iter, ++iterRGB)	{

		newPixel.SetRed(static_cast<unsigned char>(iter.Get()/2));
		newPixel.SetGreen(static_cast<unsigned char>(iter.Get()/2));
	    newPixel.SetBlue(static_cast<unsigned char>(iter.Get()/2));
	  	iterRGB.Set(newPixel);

  }


  RGBPixelType pixelValue;
  RGBImageType::IndexType pixelIndex, m;

  std::cout << "Writing Traces"<< std::endl;

  //for (unsigned int i = 0; i< mySeeds->getNumberOfSeeds()/4 ;i++){


	pixelValue[2] = 0;
	pixelValue[1] = 255;
	pixelValue[0] = 0;

	unsigned int i =0, j=0;
	Vect2 p;
	for (i=0;i<tracer->getNumberOfVessels();i++)	{
		for (j=0;j<tracer->getVessel(i)->getNumberOfSegments();j++)	{

			p = tracer->getVessel(i)->getSegment(j)->getmu();
			//print p
			pixelIndex[0] = static_cast<long int>(p[0]);
			pixelIndex[1] = static_cast<long int>(p[1]);
			RGBImage->SetPixel(pixelIndex, pixelValue);

			//std::cout<<"Trace ("<<i<<","<<j<<") - "<<pixelIndex<<std::endl;
			//std::cin.get();

		}
	}



  RGBWriterType::Pointer writer = RGBWriterType::New();

  writer->SetFileName( fname + "_Traced_Points.png" );

  writer->SetInput( RGBImage );
  try
    {
      writer->Update();
    }
  catch (itk::ExceptionObject & e)
    {
      std::cerr << "Exception in Writer: " << e << std::endl;
      exit(0);
    }
}


static void WriteTraceImageBinary(ImageType::Pointer &im, Tracer *tracer, std::string fname)
{

  typedef itk::Image< unsigned char,  2 >    BinImageType;

  BinImageType::Pointer binImage = BinImageType::New();
  BinImageType::IndexType start;
  start.Fill(0);


  BinImageType::SizeType size = im->GetRequestedRegion().GetSize();;


  BinImageType::RegionType region;
  region.SetSize( size );
  region.SetIndex( start );

  binImage->SetRegions( region );
  binImage->CopyInformation(im);
  binImage->Allocate();
  binImage->FillBuffer(0);


  itk::ImageRegionIterator<BinImageType> iterBin(binImage,binImage->GetBufferedRegion());


  unsigned char pixelValue;
  BinImageType::IndexType pixelIndex, m;



	pixelValue = 255;

	double len = 0.0;
	Vect2 p, p1, dp;

	for (unsigned int i=0;i<tracer->getNumberOfVessels();i++)	{

		if (tracer->getVessel(i)->getNumberOfSegments() <= 2)	{
			continue;
		}

		p1 = tracer->getVessel(i)->getSegment(0)->getmu() + 0.1;
		for (unsigned int j=0;j<tracer->getVessel(i)->getNumberOfSegments();j++)	{

			p = tracer->getVessel(i)->getSegment(j)->getmu();
			dp = p - p1;
			len = dp.magnitude();
			dp = dp.normalize();
			double k = 0.0;
			while( k<len)	{				pixelIndex[0] = static_cast<long int>(p[0] + k*dp[0] + 0.5);
				pixelIndex[1] = static_cast<long int>(p[1] + k*dp[1] + 0.5);
				binImage->SetPixel(pixelIndex, pixelValue);
				k++;
			}
			p1 = p;
		}
		//std::cout<<"."<<i;
	}

  std::cout << "Writing Traces ("<< tracer->getNumberOfVessels() << ")" << std::endl ;


  typedef itk::ImageFileWriter< BinImageType >  BinWriterType;
  BinWriterType::Pointer writer = BinWriterType::New();

  writer->SetFileName( fname + "_Traced_Points.tif" );
  writer->SetInput( binImage );


  try
    {
      writer->Update();
    }
  catch (itk::ExceptionObject & e)
    {
      std::cerr << "Exception in Writer: " << e << std::endl;
      exit(0);
    }
}

