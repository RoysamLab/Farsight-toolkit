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

//: Executable that generate a color image out of two gray images.

#include "itkImage.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "vnl/vnl_math.h"
#include "itkRGBPixel.h"
#include <sstream>	//stringstream
#include <string>

int main(int argc, char* argv[])
{
  if (argc<4) {
    std::cerr << "Usage: " << argv[0] << " 3D_inputimage1 3D_inputImage2 OutputImage1 [mask_image]";
    return EXIT_FAILURE;
  }

  typedef itk::Image< unsigned char, 3 > GrayImageType;
  typedef itk::RGBPixel< unsigned char >   InputPixelType;
  typedef itk::Image< InputPixelType, 3 > ImageType;
  typedef itk::ImageFileReader< GrayImageType > ReaderType;
  typedef itk::ImageRegionConstIterator< GrayImageType > RegionConstIterator;
  typedef itk::ImageRegionIterator< ImageType > RegionIterator;

  std::string input;
  float frac1, frac2;

  const char *filename1 = argv[1];
  //input =  argv[2];
  std::stringstream( input ) >> frac1;
  const char *filename2 = argv[2];
  //input =  argv[3];
  std::stringstream( input ) >> frac2;
  const char *outfilename = argv[3];

  const char *maskfilename = 0;
  if (argc == 5) maskfilename = argv[4];

  ReaderType::Pointer reader1 = ReaderType::New();
  ReaderType::Pointer reader2 = ReaderType::New();
  ReaderType::Pointer reader3 = ReaderType::New();
  GrayImageType::Pointer image1, image2, mask; 

  reader1->SetFileName( filename1 );
  reader1->Update();
  reader2->SetFileName( filename2 );
  reader2->Update();
  image1 = reader1->GetOutput();
  image2 = reader2->GetOutput();
  
  if ( argc == 5) {
    reader3->SetFileName( maskfilename );
    reader3->Update();
    mask = reader3->GetOutput();
  }

  ImageType::Pointer image_out = ImageType::New();
  image_out->SetRegions(image1->GetRequestedRegion());
  image_out->Allocate();
  image_out->FillBuffer(itk::RGBPixel<unsigned char>(itk::NumericTraits<unsigned char>::Zero));

  //Set the iterator
  RegionConstIterator inputIt1( image1, image1->GetRequestedRegion() );
  RegionConstIterator inputIt2( image2, image2->GetRequestedRegion() );
  RegionIterator outputIt( image_out, image_out->GetRequestedRegion() );
  InputPixelType pix;

  if (argc == 5) { // with mask image
    RegionConstIterator inputIt3( mask, mask->GetRequestedRegion() );
    for ( inputIt1.GoToBegin(), inputIt2.GoToBegin(), inputIt3.GoToBegin(), outputIt.GoToBegin(); !inputIt1.IsAtEnd();  ++inputIt1, ++inputIt2, ++inputIt3, ++outputIt)
      {
        if (inputIt3.Get() != 0) {
          pix.SetRed( 255 );
          pix.SetGreen( 255 );
          pix.SetBlue( 255 );
        }
        else {
          pix.SetRed( inputIt1.Get() );
          pix.SetGreen( inputIt2.Get() );
          pix.SetBlue( 0 );
        }
        outputIt.Set( pix );
      }
  }
  else { //without mask image
    for ( inputIt1.GoToBegin(), inputIt2.GoToBegin(), outputIt.GoToBegin(); 
          !inputIt1.IsAtEnd();  ++inputIt1, ++inputIt2, ++outputIt)
      {
        pix.SetRed( inputIt1.Get() );
        pix.SetGreen( inputIt2.Get() );
        pix.SetBlue( 0 );
        outputIt.Set( pix );
      }
  }

  typedef itk::ImageFileWriter< ImageType >  WriterType3D;

  WriterType3D::Pointer writer3D = WriterType3D::New();
  writer3D->SetFileName( outfilename );
  writer3D->SetInput( image_out );
  writer3D->Update();

  return 0;
}
