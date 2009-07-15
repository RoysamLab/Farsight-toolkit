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
#include "itkStatisticsImageFilter.h"
#include <vnl/vnl_vector_fixed.h>

#include "TraceConfig.h"
#include "SeedContainer3D.h"
#include "Seed2Seg.h"
#include "TraceContainer3D.h"

typedef itk::Image< float, 2 >   ImageType2D;
typedef itk::Image< float, 3 >   ImageType3D;

typedef vnl_vector_fixed<double,3> Vect3;


void ImageStatistics(ImageType3D::Pointer & );
static void WriteSeedImage(ImageType2D::Pointer, SeedContainer3D::Pointer,std::string );


int main (int argc, char *argv[])
{

	system("cls");
	std::cout << "SuperEllipsoid Trace3D version-0.1\tSept 11, 2007" << std::endl << std::endl;
	if (argc < 2)	{
		std::cout << "Usage: "<<argv[0] << " [ImageFile] [Grid spacing] [AspectRatio]" <<std::endl;
		exit(1);
	}

	std::cout << "Tracing "<< argv[1] << std::endl;
	TraceConfig::Pointer m_Config = TraceConfig::New();
	m_Config->SetFileNames(argv[1]);
	if (argc>2)	{
		m_Config->SetGridSpacing(argv[2]);
	}

	if (argc>3)	{
		m_Config->SetAspectRatio(argv[3]);
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



	try	{
		SeedContainer3D::Pointer m_Seeds = SeedContainer3D::New();
		m_Seeds->Configure(m_Config);
		m_Seeds->Detect(image3D, MIPimage);
		WriteSeedImage(MIPimage, m_Seeds , m_Config->getOutputFileName());

		Seed2Seg::Pointer m_SS = Seed2Seg::New();
		m_SS->ComuputeStartSegments(m_Seeds , image3D, m_Config);
		m_SS->SortStartSegments();

        TraceContainer3D::Pointer m_Tracer = TraceContainer3D::New();
        m_Tracer->ComputeTrace(image3D, m_SS) ;
//		m_Tracer->ApplyRules();
		m_Tracer->WriteTraceToTxtFile(m_Config->getOutputFileName());
		m_Tracer->WriteTraceToXMLFile(m_Config->getOutputFileName());

	}
	catch (itk::ExceptionObject & err )	    {
	    std::cout << "ExceptionObject caught !" << std::endl;
	    std::cout << err << std::endl;
	    return EXIT_FAILURE;
    }

	//delete m_Config; m_Config = NULL;

	return EXIT_SUCCESS;
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
			float imax2 = imean+2*istd;
			d = (d > (imax2)) ? (imax2) : d;
			d = 1 - (d - imin)/(imax2 - imin);
			it.Set(255.0*d);
		}
	}
	else {
		IteratorType it(im3D, im3D->GetRequestedRegion());
			for (it.GoToBegin(); !it.IsAtEnd(); ++it)	{
			float d = it.Get();
			float imin2 = imean-2*istd;
			d = (d < (imin2)) ? (imin2) : d;
			d = (d - imin2)/(imax - imin2);
			it.Set(255.0*d);
		}
	}

	//itk::ImageFileWriter<ImageType3D>::Pointer w = itk::ImageFileWriter<ImageType3D>::New();
	//w->SetFileName("ProcessedInput.mhd");
	//w->SetInput(im3D);
	//w->Update();

}

/*static void WriteStartSegments(Seed2Seg::Pointer stseg, std::string filename)	{

	std::ofstream ssfile;
	std::string file = filename + std::string("_SSeg.txt");
	ssfile.open (file.c_str());

	P(stseg->getNumberOfStartSegments())
	stseg->getStartSegment(1)->PrintSelf();
	S("Before Writing")

    for (unsigned int i=0; i< stseg->getNumberOfStartSegments(); i++ )	{
		TVessel* seg = stseg->getStartSegment(i);

		ssfile << "ID = " << seg->ID << " @ [" << seg->mu[0] << "," << seg->mu[1] << "," <<  seg->mu[2] << "]" << \
		" Foregd:" << seg->f << "  Backgd:" << seg->b << "  Lhood:" <<  seg->L << "  MAD:" <<  seg->MAD << \
		" A:<" << seg->a1 << ", " << seg->a2 << ", " <<  seg->a3 << "> Q:<" << seg->q1[0]<< ", "<< seg->q1[1] << ", "<< seg->q1[2] << ", " << seg->q1[3] << ">" << \
		" R:\t" << seg->R1[0] << " " << seg->R2[0] << " " <<  seg->R3[0] << " \t" << seg->R1[1] << " " << seg->R2[1] << " " <<  seg->R3[1] << " \t" << seg->R1[2] << " " << seg->R2[2] << " " <<  seg->R3[2] << std::endl;

		std::cout << "ID = " << seg->ID << " @ [" << seg->mu[0] << "," << seg->mu[1] << "," <<  seg->mu[2] << "]" << \
		" Foregd:" << seg->f << "  Backgd:" << seg->b << "  Lhood:" <<  seg->L << "  MAD:" <<  seg->MAD << \
		" A:<" << seg->a1 << ", " << seg->a2 << ", " <<  seg->a3 << "> Q:<" << seg->q1[0]<< ", "<< seg->q1[1] << ", "<< seg->q1[2] << ", " << seg->q1[3] << ">" << \
		" R:\t" << seg->R1[0] << " " << seg->R2[0] << " " <<  seg->R3[0] << " \t" << seg->R1[1] << " " << seg->R2[1] << " " <<  seg->R3[1] << " \t" << seg->R1[2] << " " << seg->R2[2] << " " <<  seg->R3[2] << std::endl;

	}
	ssfile.close();
}
*/

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


//static void WriteTraceImage(ImageType3D::Pointer&, Tracer *,std::string )	{
//}
