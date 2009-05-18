#include "itkRGBPixel.h"

#include <Common/fsc_channel_accessor.h>
#include <Common/fsc_channel_accessor.txx>

template class fsc::BlueChannelPixelAccessor<itk::RGBPixel<unsigned char> >;
template class fsc::GreenChannelPixelAccessor<itk::RGBPixel<unsigned char> >;
template class fsc::RedChannelPixelAccessor<itk::RGBPixel<unsigned char> >;

template class fsc_channel_accessor< itk::RGBPixel<unsigned char>, 3 >;
