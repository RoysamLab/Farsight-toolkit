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

#include "itkRGBAPixel.h"

#include <Common/fsc_channel_accessor.h>
#include <Common/fsc_channel_accessor.txx>

template class fsc::BlueChannelPixelAccessor<itk::RGBAPixel<unsigned char> >;
template class fsc::GreenChannelPixelAccessor<itk::RGBAPixel<unsigned char> >;
template class fsc::RedChannelPixelAccessor<itk::RGBAPixel<unsigned char> >;

template class fsc_channel_accessor< itk::RGBAPixel<unsigned char>, 3 >;
