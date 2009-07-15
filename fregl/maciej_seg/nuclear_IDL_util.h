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

#ifndef _NUCLEAR_IDL_UTIL_H_
#define _NUCLEAR_IDL_UTIL_H_

#include "maciejSegmentation.h"

//: Read in a text file containing only classification result and the seg result
// 
//  This ascii file simply contains the class label of the cells. Each
//  line is in the format of "cell_ID class_ID". The function takes
//  two main inputs: the ascii file containing the classes and the
//  .dat file, which is the binary file contaning the masks of the
//  cells. All the results are output from the Farsight IDL program.
void nuclear_IDL_read_class( std::string const& class_file,
                             std::string const& mask_file,
                             int size_x, int size_y, int size_z, 
                             maciejSegmentation& marciejSeg );

//: Read in the seg result only
// 
//  The function takes one main input: the _seg_final.dat file, which
//  is the binary file contaning the masks of the cells. The binary
//  file is the output from the Farsight IDL program.
void nuclear_IDL_read_seg( std::string const& mask_file,
                           int size_x, int size_y, int size_z, 
                           maciejSegmentation& marciejSeg );
#endif
