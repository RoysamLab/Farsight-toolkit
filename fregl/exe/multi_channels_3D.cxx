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

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRGBPixel.h"
#include "itkImage.h"
#include "vnl/vnl_math.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

#include <sstream>
#include <string>

// This program takes a list of 3D gray-scale images, and produce a 3D
// color image stack. The input file contains the information
// in the following format:
//
// image_name_#1 0/1 0/1 0/1
// image_name_#2 0/1 0/1 0/1
// ...
// 

typedef unsigned char    PixelType;
typedef itk::Image< PixelType, 3 > ImageTypeGray;
typedef itk::RGBPixel< unsigned char >   colorPixelType;
typedef itk::Image< colorPixelType, 3 > ImageType;
typedef itk::ImageFileReader< ImageTypeGray > ReaderType;
typedef itk::ImageFileWriter< ImageType >  WriterType;
typedef itk::ImageRegionConstIterator< ImageTypeGray > ConstRegionIteratorType;
typedef itk::ImageRegionIterator< ImageType > RegionIteratorType;

int main(int argc, char* argv[])
{
  if (argc<3) {
    std::cerr << "Usage: " << argv[0] << " image_list 3D_color_image";
    return EXIT_FAILURE;
  }
  
  const char *filename_list = argv[1];
  //const char *outfilename   = argv[2];

  std::ifstream in_file_str( filename_list);
  if ( !in_file_str ){
    std::cerr<<"Couldn't open "<<filename_list<<std::endl;
    exit( 0 );
  }
  
  std::string filename, line_str;
  float r, g, b;
  ImageType::Pointer image_final;
  while ( in_file_str ) {
    std::getline(in_file_str, line_str);
    if (line_str.length() == 0) continue;

    std::istringstream line_stream(line_str);
    line_stream>>filename>>r>>g>>b;

    ReaderType::Pointer reader = ReaderType::New();

    reader->SetFileName( filename );
    reader->Update();
    ImageTypeGray::Pointer image = reader->GetOutput();

    // set the 2D image if not yet initialized
    if (!image_final) {
      image_final = ImageType::New();
      image_final->SetRegions( image->GetRequestedRegion() );
      image_final->Allocate();
      image_final->FillBuffer(itk::RGBPixel<unsigned char>(itk::NumericTraits<unsigned char>::Zero));
    }

    //Set the iterator
    ConstRegionIteratorType It1( image, image->GetRequestedRegion() );
    RegionIteratorType It2( image_final, image_final->GetRequestedRegion() );

 
    // Now do the max projection on each color channel separately
    PixelType pix1;
    colorPixelType pix2, pix3;
    for (It1.GoToBegin(), It2.GoToBegin(); !It1.IsAtEnd(); ++It1, ++It2) {
      pix1 = It1.Get();
      pix2 = It2.Get();
      pix3.SetRed( vnl_math_max( int(pix1 * r), pix2.GetRed() ) );
      pix3.SetGreen( vnl_math_max( int(pix1 * g), pix2.GetGreen() ) );
      pix3.SetBlue( vnl_math_max( int(pix1 * b), pix2.GetBlue() ) );
      It2.Set( pix3 );
    }
  }

  // Output the 2D projected image
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( image_final );
  writer->Update();
  
  return 0;
}
