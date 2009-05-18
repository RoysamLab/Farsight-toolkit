#include "itkRGBAPixel.h"

#include <Common/fsc_channel_accessor.h>
#include <Common/fsc_channel_accessor.txx>

template class fsc::BlueChannelPixelAccessor<itk::RGBAPixel<unsigned char> >;
template class fsc::GreenChannelPixelAccessor<itk::RGBAPixel<unsigned char> >;
template class fsc::RedChannelPixelAccessor<itk::RGBAPixel<unsigned char> >;

template class fsc_channel_accessor< itk::RGBAPixel<unsigned char>, 3 >;
