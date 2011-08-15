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
#include "itkDownwardFrontPropagationLabelImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

int main( int argc, char *argv[] )
{
  if( argc < 6 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputLabelImage inputImage outputLabelImage";
    std::cerr << " maximumNumberOfIterations lowerThresholdValue";
    std::cerr << std::endl;
    return 1;
    }

  typedef   double          InputPixelType;
  typedef   unsigned short  MaskPixelType;

  const     unsigned int    Dimension = 3;

  typedef itk::Image< InputPixelType, Dimension >     InputImageType;
  typedef itk::Image< MaskPixelType, Dimension >      MaskImageType;

  typedef itk::DownwardFrontPropagationLabelImageFilter< 
    MaskImageType, InputImageType, MaskImageType >    ImageFilterType;


  typedef  itk::ImageFileReader< InputImageType > InputReaderType;
  typedef  itk::ImageFileReader< MaskImageType >  MaskReaderType;

  typedef ImageFilterType::OutputImageType        OutputImageType;
                        
  typedef  itk::ImageFileWriter<  OutputImageType  > WriterType;

  InputReaderType::Pointer inputReader = InputReaderType::New();
  MaskReaderType::Pointer  maskReader  = MaskReaderType::New();

  WriterType::Pointer writer = WriterType::New();

  maskReader->SetFileName( argv[1] );
  inputReader->SetFileName( argv[2] );

  ImageFilterType::Pointer floodFilter = ImageFilterType::New();

  const unsigned int maximumNumberOfIterations = atoi( argv[4] );
  const signed short lowerThreshold = atoi( argv[5] );

  floodFilter->SetMaximumNumberOfIterations( maximumNumberOfIterations );
  floodFilter->SetMajorityThreshold( 1 );
  floodFilter->SetLowerThreshold( lowerThreshold );
  floodFilter->InPlaceOn();

  InputImageType::SizeType  ballManhattanRadius;
  ballManhattanRadius.Fill( 1 );
  floodFilter->SetRadius( ballManhattanRadius );

  floodFilter->SetInput( maskReader->GetOutput() );
  floodFilter->SetFeatureImage( inputReader->GetOutput() );

  writer->SetInput( floodFilter->GetOutput() );
 
  writer->UseCompressionOn();

  writer->SetFileName( argv[3] );

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
