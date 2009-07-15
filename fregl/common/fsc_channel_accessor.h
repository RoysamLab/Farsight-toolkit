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

#ifndef fsc_channel_accessor_h
#define fsc_channel_accessor_h
//:
// \file
// \brief The class accesses to a specified color channel in the input image
// \author Charlene Tsai
// \date 06 September 2007
//

#include "itkImage.h"
#include "itkImageSeriesReader.h"
#include "itkImageAdaptor.h"
#include "itkRescaleIntensityImageFilter.h"

#include <string>
#include <vector>

namespace fsc{

  //: Helper class for blue channel accesor
  template< typename IntType>
  class BlueChannelPixelAccessor  
  {
  public:
    //typedef itk::RGBPixel<unsigned char>   InternalType;
	typedef typename   IntType::Self         InternalType;
    typedef            unsigned char    ExternalType;
    
    static ExternalType Get( const InternalType & input ) 
    {
      return static_cast<ExternalType>( input.GetBlue() );
    }
  };
  
  //: Helper class for red channel accesor
  template< typename IntType>
  class RedChannelPixelAccessor  
  {
  public:
    //typedef itk::RGBPixel<unsigned char>   InternalType;
	typedef typename   IntType::Self       InternalType;
    typedef            unsigned char    ExternalType;

    static ExternalType Get( const InternalType & input ) 
    {
      return static_cast<ExternalType>( input.GetRed() );
    }
  };

  //: Helper class for green channel accesor
  template< typename IntType>
  class GreenChannelPixelAccessor  
  {
  public:
    //typedef itk::RGBPixel<unsigned char>   InternalType;
	typedef typename    IntType::Self      InternalType;
    typedef             unsigned char    ExternalType;
    
    static ExternalType Get( const InternalType & input ) 
    {
      return static_cast<ExternalType>( input.GetGreen() );
    }
  };
} //namespace 

//: The class that allows access to the input image
template< typename InputPixelType, unsigned int  Dim >
class fsc_channel_accessor
{
public: 

  enum channel_type {RED, GREEN, BLUE};

  //static const unsigned int               Dim = 3;
  //typedef itk::RGBPixel<unsigned char>          RGBInputPixelType;
  typedef typename itk::Image< InputPixelType, Dim >   RGBInputImageType;
  typedef typename RGBInputImageType::Pointer          InputImageTypePtr;
  typedef unsigned char                                PixelType;
  typedef typename itk::Image<PixelType, Dim>          ImageType;
  typedef typename ImageType::Pointer                  ImageTypePtr;

  typedef typename itk::ImageSeriesReader< RGBInputImageType  > SeriesReaderType;
  typedef typename SeriesReaderType::Pointer           SeriesReaderTypePtr; 
  typedef typename itk::ImageFileReader< RGBInputImageType  >   ReaderType;
  typedef typename ReaderType::Pointer                ReaderTypePtr;
  typedef std::vector<std::string>                     FileNamesContainer;

  //adaptor types for channels
  typedef typename itk::ImageAdaptor< RGBInputImageType, 
                             fsc::GreenChannelPixelAccessor<InputPixelType> > ImageGreenAdaptorType;
  typedef typename ImageGreenAdaptorType::Pointer ImageGreenAdaptorTypePtr;
  typedef typename itk::ImageAdaptor< RGBInputImageType, 
                             fsc::RedChannelPixelAccessor<InputPixelType> > ImageRedAdaptorType;
  typedef typename ImageRedAdaptorType::Pointer  ImageRedAdaptorTypePtr;
  typedef typename itk::ImageAdaptor< RGBInputImageType, 
                             fsc::BlueChannelPixelAccessor<InputPixelType> > ImageBlueAdaptorType;
  typedef typename ImageBlueAdaptorType::Pointer ImageBlueAdaptorTypePtr;

  //rescaler for adaptor
  typedef typename itk::RescaleIntensityImageFilter< ImageGreenAdaptorType, 
                                            ImageType 
                                            >  RescalerGreenAdaptorType;
  typedef typename RescalerGreenAdaptorType::Pointer RescalerGreenAdaptorTypePtr;
  typedef typename itk::RescaleIntensityImageFilter< ImageBlueAdaptorType, 
                                            ImageType 
                                            >  RescalerBlueAdaptorType;
  typedef typename RescalerBlueAdaptorType::Pointer RescalerBlueAdaptorTypePtr; 
  typedef typename itk::RescaleIntensityImageFilter< ImageRedAdaptorType, 
                                            ImageType 
                                            >   RescalerRedAdaptorType;
  typedef typename RescalerRedAdaptorType::Pointer RescalerRedAdaptorTypePtr;

  //: Constructor for setting the image
  fsc_channel_accessor( InputImageTypePtr image_ );

  //: Constructor for loading a series of image slices
  fsc_channel_accessor(FileNamesContainer const& file_name);

  //: Constructor for loading one 3D image
  fsc_channel_accessor(std::string const& file_name);

  //: Destructor
  virtual ~fsc_channel_accessor();
  
  //: Load a new set of image slices
  void
  load_images( FileNamesContainer const& file_name );

  //: Load one 3D image
  void
  load_image( std::string const& file_name );

  //: Get a channel
  ImageTypePtr
  get_channel( channel_type );

private: 
  InputImageTypePtr          image_;
 
};
#endif
