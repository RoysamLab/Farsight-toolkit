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
