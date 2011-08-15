/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: FrontPropagationLabelImageFilter.cxx,v $
  Language:  C++
  Date:      $Date: 2010-08-17 16:14:27 -0500 (Tue, 17 Aug 2010) $
  Version:   $Revision: 2065 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "itkImage.h"
#include "itkHoleFillingFrontPropagationLabelImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

int main( int argc, char *argv[] )
{
  if( argc < 5 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputLabelImage outputLabelImage";
    std::cerr << " maximumNumberOfIterations radius";
    std::cerr << std::endl;
    return 1;
    }

  typedef   unsigned char   PixelType;
  const     unsigned int    Dimension = 3;

  typedef itk::Image< PixelType, Dimension >      ImageType;

  typedef itk::HoleFillingFrontPropagationLabelImageFilter< 
    ImageType, ImageType >    ImageFilterType;


  typedef  itk::ImageFileReader< ImageType >  ReaderType;

  typedef ImageFilterType::OutputImageType        OutputImageType;
                        
  typedef  itk::ImageFileWriter<  OutputImageType  > WriterType;

  ReaderType::Pointer  reader  = ReaderType::New();

  WriterType::Pointer writer = WriterType::New();

  reader->SetFileName( argv[1] );

  ImageFilterType::Pointer floodFilter = ImageFilterType::New();

  const unsigned int maximumNumberOfIterations = atoi( argv[3] );

  floodFilter->SetMaximumNumberOfIterations( maximumNumberOfIterations );
  floodFilter->SetMajorityThreshold( 1 );
  floodFilter->SetForegroundValue( 0 );
  floodFilter->SetBackgroundValue( 255 );
  floodFilter->InPlaceOn();

  ImageType::SizeType  ballManhattanRadius;
  ballManhattanRadius.Fill( atoi( argv[4])  );
  floodFilter->SetRadius( ballManhattanRadius );

  floodFilter->SetInput( reader->GetOutput() );

  writer->SetInput( floodFilter->GetOutput() );
 
  writer->UseCompressionOn();

  writer->SetFileName( argv[2] );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }

  std::cout << "Number of iterations used = ";
  std::cout << floodFilter->GetCurrentIterationNumber() << std::endl;

  std::cout << "Number of pixels changed = ";
  std::cout << floodFilter->GetTotalNumberOfPixelsChanged() << std::endl;

  return EXIT_SUCCESS;
}
