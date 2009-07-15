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

#include "fsc_channel_accessor.h"

#include "itkImageFileReader.h"
#include "itkImageSeriesReader.h"

//#include "itkTIFFImageIO.h"
//#include "itkPNGImageIO.h"

//: Constructor
template< typename InputPixelType, unsigned int  Dim >
fsc_channel_accessor<InputPixelType,Dim>::
fsc_channel_accessor(InputImageTypePtr image)
{
  image_ = image;
}

//: Constructor
template< typename InputPixelType, unsigned int  Dim >
fsc_channel_accessor<InputPixelType,Dim>::
fsc_channel_accessor(FileNamesContainer const& file_names)
{
  load_images( file_names );
}

//: Constructor
template< typename InputPixelType, unsigned int  Dim >
fsc_channel_accessor<InputPixelType,Dim>::
fsc_channel_accessor(std::string const& file_name)
{
  load_image( file_name );
}

//: Destructor
template< typename InputPixelType, unsigned int  Dim >
fsc_channel_accessor<InputPixelType,Dim>::
~fsc_channel_accessor()
{
}

//: Load a new set of image slices
template< typename InputPixelType, unsigned int  Dim >
void
fsc_channel_accessor<InputPixelType,Dim>::
load_images( FileNamesContainer const& file_names )
{
  // Create and set the reader
  SeriesReaderTypePtr reader = SeriesReaderType::New();

  // CT: When tested in Cygwin, I ran into problem with absolute file
  // path, while relative path works fine. I'm not sure if it is the
  // problem with Cygwin, or the ImageSeriesReader, or ...
  reader->SetFileNames(file_names);

  // Create image
  try
    {
      reader->Update();
    }
  catch (itk::ExceptionObject & e)
    {
      std::cerr << "Exception in ImageSeriesReader: " << e << std::endl;
      exit(0);
    }
    image_ = reader->GetOutput();
}

//: Load one image slice
template< typename InputPixelType, unsigned int  Dim >
void
fsc_channel_accessor<InputPixelType,Dim>::
load_image( std::string const& file_name )
{
  // Create and set the reader
  ReaderTypePtr reader = ReaderType::New();
  reader->SetFileName(file_name);

  // Create image
  try
    {
      reader->Update();
    }
  catch (itk::ExceptionObject & e)
    {
      std::cerr << "Exception in ImageSeriesReader: " << e << std::endl;
      exit(0);
    }
    image_ = reader->GetOutput();
}

//: Get a channel
template< typename InputPixelType, unsigned int  Dim >
typename fsc_channel_accessor<InputPixelType,Dim>::ImageTypePtr
fsc_channel_accessor<InputPixelType,Dim>::
get_channel( channel_type type)
{
  ImageTypePtr single_channel_image;

  switch (type) {
  case GREEN: 
    {
      ImageGreenAdaptorTypePtr g_adaptor = ImageGreenAdaptorType::New();
      g_adaptor->SetImage( image_ );
      RescalerGreenAdaptorTypePtr g_rescaler_adap = RescalerGreenAdaptorType::New();
      g_rescaler_adap->SetOutputMinimum(  0  );
      g_rescaler_adap->SetOutputMaximum( 255 );
      
      g_rescaler_adap->SetInput(g_adaptor);
      try
        {
          g_rescaler_adap->Update();
        }
      catch (itk::ExceptionObject & e)
        {
          std::cerr << "Exception in RescalerAdaptor: " << e << std::endl;
          exit(0);
        }
      single_channel_image = g_rescaler_adap->GetOutput();
    }
    break;
  case RED:
    {
      ImageRedAdaptorTypePtr r_adaptor = ImageRedAdaptorType::New();
      r_adaptor->SetImage( image_ );
      RescalerRedAdaptorTypePtr r_rescaler_adap = RescalerRedAdaptorType::New();
      r_rescaler_adap->SetOutputMinimum(  0  );
      r_rescaler_adap->SetOutputMaximum( 255 );
      
      r_rescaler_adap->SetInput(r_adaptor);
      try
        {
          r_rescaler_adap->Update();
        }
      catch (itk::ExceptionObject & e)
        {
          std::cerr << "Exception in RescalerAdaptor: " << e << std::endl;
          exit(0);
        }
      single_channel_image = r_rescaler_adap->GetOutput();
    }
    break;
  case BLUE:
     {
      ImageBlueAdaptorTypePtr b_adaptor = ImageBlueAdaptorType::New();
      b_adaptor->SetImage( image_ );
      RescalerBlueAdaptorTypePtr b_rescaler_adap = RescalerBlueAdaptorType::New();
      b_rescaler_adap->SetOutputMinimum(  0  );
      b_rescaler_adap->SetOutputMaximum( 255 );
      
      b_rescaler_adap->SetInput(b_adaptor);
      try
        {
          b_rescaler_adap->Update();
        }
      catch (itk::ExceptionObject & e)
        {
          std::cerr << "Exception in RescalerAdaptor: " << e << std::endl;
          exit(0);
        }
      single_channel_image = b_rescaler_adap->GetOutput();
    }
    break;
  default:
    {
      std::cerr <<"No such channel defined"<< std::endl;
      exit(0);
    }
  }
  single_channel_image->DisconnectPipeline();
  return single_channel_image;
}
