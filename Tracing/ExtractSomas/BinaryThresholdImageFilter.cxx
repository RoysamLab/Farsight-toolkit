/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: BinaryThresholdImageFilter.cxx,v $
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

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif


// Software Guide : BeginCommandLineArgs
// INPUTS:  {BrainProtonDensitySlice.png}
// OUTPUTS: {BinaryThresholdImageFilterOutput.png}
// 150 180 0 255
// Software Guide : EndCommandLineArgs


#include "itkBinaryThresholdImageFilter.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

int main( int argc, char * argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImageFile outputImageFile ";  
    std::cerr << " lowerThreshold " << std::endl;
    return EXIT_FAILURE;
    }
  

  typedef  double         InputPixelType;
  typedef  unsigned char  OutputPixelType;

  typedef itk::Image< InputPixelType,  3 >   InputImageType;
  typedef itk::Image< OutputPixelType, 3 >   OutputImageType;

  typedef itk::BinaryThresholdImageFilter<
               InputImageType, OutputImageType >  FilterType;

  typedef itk::ImageFileReader< InputImageType >  ReaderType;

  typedef itk::ImageFileWriter< OutputImageType >  WriterType;



  ReaderType::Pointer reader = ReaderType::New();
  FilterType::Pointer filter = FilterType::New();

  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( filter->GetOutput() );
  reader->SetFileName( argv[1] );



  filter->SetInput( reader->GetOutput() );



  const OutputPixelType outsideValue = 0;
  const OutputPixelType insideValue  = 255;

  filter->SetOutsideValue( outsideValue );
  filter->SetInsideValue(  insideValue  );


  const InputPixelType lowerThreshold = atoi( argv[3] );
  const InputPixelType upperThreshold = 10000.0;

  filter->SetLowerThreshold( lowerThreshold );
  filter->SetUpperThreshold( upperThreshold );


  filter->Update();


  writer->SetFileName( argv[2] );
  writer->Update();

  return EXIT_SUCCESS;
}
