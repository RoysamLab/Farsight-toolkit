/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: FrontPropagationLabelImageFilter.cxx,v $
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

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
#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

int main( int argc, char *argv[] )
{
  if( argc < 5 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputBinaryImage outputBinaryImage";
    std::cerr << " maximumNumberOfIterations radius";
    std::cerr << std::endl;
    return 1;
    }

  typedef   unsigned char   PixelType;
  const     unsigned int    Dimension = 3;

  typedef itk::Image< PixelType, Dimension >      ImageType;

  typedef itk::VotingBinaryIterativeHoleFillingImageFilter< ImageType >
    ImageFilterType;


  typedef  itk::ImageFileReader< ImageType >  ReaderType;

  typedef ImageFilterType::OutputImageType        OutputImageType;
                        
  typedef  itk::ImageFileWriter<  OutputImageType  > WriterType;

  ReaderType::Pointer  reader  = ReaderType::New();

  WriterType::Pointer writer = WriterType::New();

  reader->SetFileName( argv[1] );

  ImageFilterType::Pointer holeFiller = ImageFilterType::New();

  const unsigned int maximumNumberOfIterations = atoi( argv[3] );

  holeFiller->SetMaximumNumberOfIterations( maximumNumberOfIterations );
  holeFiller->SetMajorityThreshold( 1 );
  holeFiller->SetForegroundValue( 255 );
  holeFiller->SetBackgroundValue( 0 );

  ImageType::SizeType  ballManhattanRadius;
  ballManhattanRadius.Fill( atoi( argv[4])  );
  holeFiller->SetRadius( ballManhattanRadius );

  holeFiller->SetInput( reader->GetOutput() );
  writer->SetInput( holeFiller->GetOutput() );
 
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
  std::cout << holeFiller->GetCurrentNumberOfIterations() << std::endl;

  std::cout << "Number of pixels changed = ";
  std::cout << holeFiller->GetNumberOfPixelsChanged() << std::endl;

  return EXIT_SUCCESS;
}
