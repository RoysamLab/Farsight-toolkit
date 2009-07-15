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
#include "itkImageSliceConstIteratorWithIndex.h"
#include "itkImageSliceIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"

#include <sstream>
#include <string>

// This program takes a list of 3D gray-scale images, and produce a 2D
// color image which contains the maximum intensity projection of
// individual color channel. The input file contains the information
// in the following format:
//
// image_name_#1 0/1 0/1 0/1
// image_name_#2 0/1 0/1 0/1
// ...
// 

typedef unsigned char                  PixelType;
typedef itk::RGBPixel< unsigned char >   colorPixelType;
typedef itk::Image< unsigned char, 3 > ImageType3D;
typedef itk::Image< colorPixelType, 2 > ImageType2D;
typedef itk::ImageFileReader< ImageType3D > ReaderType;
typedef itk::ImageLinearIteratorWithIndex< ImageType2D > LinearIteratorType;
typedef itk::ImageSliceIteratorWithIndex< ImageType3D > SliceIteratorType1;
typedef itk::ImageFileWriter< ImageType2D >  WriterType2D;

int main(int argc, char* argv[])
{
  if (argc<3) {
    std::cerr << "Usage: " << argv[0] << " image_list 2D_outimage";
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
  ImageType2D::Pointer image2D;
  while ( in_file_str ) {
    std::getline(in_file_str, line_str);
    if (line_str.length() == 0) continue;

    std::istringstream line_stream(line_str);
    line_stream>>filename>>r>>g>>b;

    ReaderType::Pointer reader = ReaderType::New();

    reader->SetFileName( filename );
    reader->Update();
    ImageType3D::Pointer image3D = reader->GetOutput();


    ImageType2D::RegionType region;
    ImageType2D::RegionType::SizeType size;
    ImageType2D::RegionType::IndexType index;

    // set the 2D image if not yet initialized
    if (!image2D) {
      image2D = ImageType2D::New();
      ImageType3D::RegionType requestedRegion = image3D->GetRequestedRegion();
      index[ 0 ] = requestedRegion.GetIndex()[ 0 ];
      index[ 1 ] = requestedRegion.GetIndex()[ 1 ];
      size[ 0 ] = requestedRegion.GetSize()[ 0 ];
      size[ 1 ] = requestedRegion.GetSize()[ 1 ];
      region.SetSize( size );
      region.SetIndex( index );
      
      image2D->SetRegions( region );
      image2D->Allocate();
      image2D->FillBuffer(itk::RGBPixel<unsigned char>(itk::NumericTraits<unsigned char>::Zero));
      std::cout<<"Done with 2D image initialization"<<std::endl;
    }
    std::cout<<"Projecting "<<filename<<std::endl;

    //Set the iterator
    //SliceIteratorType2 inputIt( color_image, color_image->GetRequestedRegion() );
    SliceIteratorType1 output3DIt( image3D, image3D->GetRequestedRegion() );
    LinearIteratorType output2DIt( image2D, image2D->GetRequestedRegion() );

    unsigned int direction[2];
    direction[0] = 0;
    direction[1] = 1;
  
    output3DIt.SetFirstDirection( direction[1] );
    output3DIt.SetSecondDirection( direction[0] );
    output2DIt.SetDirection( 1 - direction[0] );
  
    // Now do the max projection on each color channel separately
    output3DIt.GoToBegin();
    output2DIt.GoToBegin();
    PixelType pix1;
    colorPixelType pix2, pix3;
    while( !output3DIt.IsAtEnd() ) {
      while ( !output3DIt.IsAtEndOfSlice() ) {
        while ( !output3DIt.IsAtEndOfLine() ) {
          pix1 = output3DIt.Get();
          pix2 = output2DIt.Get();
          pix3.SetRed( vnl_math_max( int(pix1 * r), pix2.GetRed() ) );
          pix3.SetGreen( vnl_math_max( int(pix1 * g), pix2.GetGreen() ) );
          pix3.SetBlue( vnl_math_max( int(pix1 * b), pix2.GetBlue() ) );
          output2DIt.Set( pix3 );
          ++output3DIt;
          ++output2DIt;
        }
        output2DIt.NextLine();
        output3DIt.NextLine();
      }
      output2DIt.GoToBegin();
      output3DIt.NextSlice();
    }
  }
  // Output the 2D projected image
  WriterType2D::Pointer writer2D = WriterType2D::New();
  writer2D->SetFileName( argv[2] );
  writer2D->SetInput( image2D );
  writer2D->Update();
  
  return 0;
}
