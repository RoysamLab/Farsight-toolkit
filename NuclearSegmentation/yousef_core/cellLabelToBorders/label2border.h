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

#ifndef LABEL2BORDER_H
#define LABEL2BORDER_H

#include "itkLabelBorderImageFilter.h"
#include "itkImageFileWriter.h"

void label2border(int* img,int c, int r, int z);
void label2border2(int* img,int* bordImg, int c, int r, int z); // Added by Yousef on May 29th 2008
void label2border3(int* img,int* segImg, int c, int r, int z);// Added by Yousef on 9-13-2008

#endif

