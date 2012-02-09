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

//: Executable program to fuse two 3D image stacks to one. The input
//can be gray or rgb color. The output is always a gray image. The
//maximum intensity at each pixel is taken.
//
#include "itkImage.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRGBAPixel.h"
#include "itkTIFFImageIO.h"

#include "vnl/vnl_math.h"

//#include <Common/fsc_channel_accessor.h>
#include <fregl/fregl_util.h>

/*
ImageType::Pointer
read_image( std::string const & file_name, int channel )
{
  std::cout<<"Reading the image "<<file_name<<std::endl;

  ImageType::Pointer image;

  // Get pixel information
  itk::TIFFImageIO::Pointer io = itk::TIFFImageIO::New();
  io->SetFileName(file_name);
  io->ReadImageInformation();
  int pixel_type = (int)io->GetPixelType();
  std::cout<<"Pixel Type = "<<pixel_type<<std::endl; //1 - grayscale, 2-RGB, 3-RGBA, etc.,

  if (pixel_type == 3) { //RGBA pixel type
    typedef fsc_channel_accessor<itk::RGBAPixel<unsigned char>,3 > ChannelAccType;
    ChannelAccType channel_accessor(file_name);
    image = channel_accessor.get_channel(ChannelAccType::channel_type(channel));
  }
  else if (pixel_type == 2) { //RGA pixel type
    typedef fsc_channel_accessor<itk::RGBPixel<unsigned char>,3 > ChannelAccType;
    ChannelAccType channel_accessor(file_name);
    image = channel_accessor.get_channel(ChannelAccType::channel_type(channel));
  }
  else {// Gray image
    typedef itk::ImageFileReader< ImageType > ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( file_name );
    try {
      reader->Update();
    }
    catch(itk::ExceptionObject& e) {
      vcl_cout << e << vcl_endl;
    }
    image =  reader->GetOutput();
  }
  return image;
}
*/

int main(int argc, char* argv[])
{
	typedef unsigned short PixelType;
	typedef itk::Image< PixelType, 3 > ImageType;

	typedef itk::ImageRegionConstIterator< ImageType > RegionConstIterator;
	typedef itk::ImageRegionIterator< ImageType > RegionIterator;
	
	if (argc<4) {
    std::cerr << "Usage: " << argv[0] << " Inputimage1 InputImage2 OutputGrayImage rgb_channel_1 rgb_channel_2  ";
    return EXIT_FAILURE;
  }

  std::string filename1 = argv[1];
  std::string filename2 = argv[2];
  const char *outfilename = argv[3];
  int channel1 = 0, channel2 = 0;

  if (argc > 4) channel1 = atoi(argv[4]);
  if (argc > 5) channel2 = atoi(argv[5]);
  else channel2 = channel1;

  ImageType::Pointer image1, image2, image_out;
  image1 = fregl_util< PixelType >::fregl_util_read_image( filename1, true, channel1, false );
  image2 = fregl_util< PixelType >::fregl_util_read_image( filename2, true, channel2, false );

  // Perform the fusing here
  //
  std::cout<<"Fusing the images"<<std::endl;
  image_out = ImageType::New();
  image_out->SetRegions(image1->GetRequestedRegion());
  image_out->Allocate();

  //Set the iterator
  RegionConstIterator inputIt1( image1, image1->GetRequestedRegion() );
  RegionConstIterator inputIt2( image2, image2->GetRequestedRegion() );
  RegionIterator outputIt( image_out, image_out->GetRequestedRegion() );

  for ( inputIt1.GoToBegin(), inputIt2.GoToBegin(), outputIt.GoToBegin(); 
        !inputIt1.IsAtEnd();  ++inputIt1, ++inputIt2, ++outputIt)
      {
        outputIt.Set( vnl_math_max(inputIt1.Get(), inputIt2.Get()) );
      }

  typedef itk::ImageFileWriter< ImageType >  WriterType3D;

  WriterType3D::Pointer writer3D = WriterType3D::New();
  writer3D->SetFileName( outfilename );
  writer3D->SetInput( image_out );
  writer3D->Update();

  return 0;
}

