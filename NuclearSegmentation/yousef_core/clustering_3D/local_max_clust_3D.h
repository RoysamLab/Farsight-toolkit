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

#ifndef LOCAL_MAX_CLUST_3D_H
#define LOCAL_MAX_CLUST_3D_H

#include <iostream>
#include <algorithm>

void local_max_clust_3D(float* im_vals, unsigned short* local_max_vals, unsigned short* bImg, unsigned short* out1, int r, int c, int z, int scale_xy, int scale_z);

#endif
