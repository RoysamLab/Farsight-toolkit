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

#include "nuclear_IDL_util.h"

#include <fstream>
#include <sstream>

#include "rich_cell.h"

void nuclear_IDL_read_seg(std::string const& mask_file,
                          int size_x, int size_y, int size_z, 
                          maciejSegmentation& marciejSeg)
{
  // Get the fake input parameter list to feed marciejSeg
  std::vector<std::string>        parameters;
  parameters.push_back( "1" ); //channel
  parameters.push_back("500" ); //threshold
  parameters.push_back("4"); //grid
  parameters.push_back("1"); // "median_filter_radius"
  parameters.push_back("1"); //morph_opt_radius
  parameters.push_back("1"); // "mode"
  parameters.push_back("");
  parameters.push_back("");
  marciejSeg.set_parameters( parameters );
  
  // read in the binary mask file
  int array_size = size_x*size_y*size_z;
  unsigned short * memblock = new unsigned short[array_size];
  
  std::ifstream dat_file(mask_file.c_str(), std::ios::binary);
  if (dat_file.is_open()) {
    std::cerr<<"Opened .dat file "<<mask_file<<std::endl;
    dat_file.read(reinterpret_cast<char *>(memblock), sizeof(unsigned short)*array_size);
  }
  else {
    std::cerr<<"Cannot open the .dat file "<<mask_file<<std::endl;
    delete [] memblock;
    return;
  }
  dat_file.close();

  int max_id = 0;
  for (int i = 0; i<array_size; i++) {
    if (memblock[i]>max_id) max_id = memblock[i];
  }

  vcl_vector<rich_cell::Pointer> ids(max_id+1, rich_cell::Pointer());
  for (int i = 0; i<array_size; i++) {
    if (memblock[i] < 1) continue; //0 is the background

    if (!ids[ memblock[i] ]) {
      rich_cell::Pointer cell = rich_cell::New();
      ids[ memblock[i] ] = cell;
      cell->label_ =  memblock[i];
      cell->valid_ = true;
      cell->class_type_ = 1;
    }
  }

  for (unsigned int i = 0; i<ids.size(); i++) {
    if (ids[i])  marciejSeg.add_cell( ids[i] );
  }

  marciejSeg.update_cell_pixels(memblock, size_x, size_y, size_z);

  delete [] memblock;
}

void nuclear_IDL_read_class(std::string const& class_file,
                            std::string const& mask_file,
                            int size_x, int size_y, int size_z, 
                            maciejSegmentation& marciejSeg)
{
  // Get the fake input parameter list to feed marciejSeg
  std::vector<std::string>        parameters;
  parameters.push_back( "1" ); //channel
  parameters.push_back("500" ); //threshold
  parameters.push_back("4"); //grid
  parameters.push_back("1"); // "median_filter_radius"
  parameters.push_back("1"); //morph_opt_radius
  parameters.push_back("1"); // "mode"
  parameters.push_back("");
  parameters.push_back("");
  marciejSeg.set_parameters( parameters );
  int id_count = 0;
  std::string line_str;
  std::ifstream class_info( class_file.c_str() );
  std::string spaces=" /t";
  while ( class_info && !class_info.eof() ) {
    rich_cell::Pointer cell = rich_cell::New();
    std::getline(class_info, line_str);
    if (line_str.length() == 0) continue;

    std::istringstream line_stream(line_str);
    if ( line_str.find_first_of(spaces) != std::string::npos) //each line contains "Id    class"
      line_stream>>cell->label_>>cell->class_type_;
    else { // each line contain "class" where the id is the line#
      id_count++;
      cell->label_ = id_count;
      line_stream>>cell->class_type_;
    }
    cell->valid_ = true;
    marciejSeg.add_cell( cell );
  }
  class_info.close();

  // Read in the label image and update the cell pixels
  marciejSeg.update_cell_pixels(mask_file, size_x, size_y, size_z);
}

