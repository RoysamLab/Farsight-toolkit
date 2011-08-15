/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkSignedMaurerDistanceMapImageFilterTest.cxx,v $
  Language:  C++
  Date:      $Date: 2010-09-07 09:07:48 -0500 (Tue, 07 Sep 2010) $
  Version:   $Revision: 2079 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkDanielssonDistanceMapImageFilter.h"

int main( int argc, char * argv[] )
{
   if(argc < 3)
    {
    std::cerr << "Usage: " << argv[0] << " InputImage OutputImage\n";
    return -1;
    }

  const unsigned int      ImageDimension = 3;
  typedef unsigned char   PixelType;

  typedef itk::Image<PixelType,  ImageDimension>  ImageType;

  typedef itk::ImageFileReader<ImageType>    ReaderType;
  typedef itk::ImageFileWriter<ImageType>   WriterType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(argv[1]);
  reader->Update();  

  typedef itk::DanielssonDistanceMapImageFilter
     <ImageType, ImageType>  FilterType;

  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( reader->GetOutput() );
  filter->SetSquaredDistance( false );
  filter->SetUseImageSpacing( true );
  filter->InputIsBinaryOn();
  filter->Update();

  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( filter->GetVoronoiMap() );
  writer->SetFileName( argv[2] );
  writer->Update();

  return EXIT_SUCCESS;
};
