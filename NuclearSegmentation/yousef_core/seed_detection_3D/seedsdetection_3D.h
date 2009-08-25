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

#ifndef SEEDSDETECTION_3D_H
#define SEEDSDETECTION_3D_H

int Seeds_Detection_3D( float* IM, float* IM_out, int* IM_bin, int r, int c, int z, double sigma_min, double sigma_max, double scale_xy, double scale_z,int sampl_ratio,unsigned short* bImg, int UseDistMap);

#endif
