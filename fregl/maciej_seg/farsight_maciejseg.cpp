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

#include "farsight_maciejseg.h"
#include "xml_util.h"
#include "rich_cell.h"

#include <Common/fsc_channel_accessor.h>

void 
farsight_maciejseg:: 
set_image( std::string const& image_path,
           std::string const& image_name)
{
  farsight_object<3>::set_image( image_path, image_name);
 
  if (parameters_.size() == 0) {
    std::cerr<<"Parameters should be set first!"<<std::endl;
    exit(0);
  }

  //Read the RGB image and extract only one channel
  typedef fsc_channel_accessor<itk::RGBPixel<unsigned char>,3 > ChannelAccType;

  std::string filename = image_path + image_name;
  std::cout<<"Filename = "<<filename<<std::endl;
  ChannelAccType channel_accessor(filename);
  
  int channel = atoi(parameters_[0].c_str());
  intensity_image_ = channel_accessor.get_channel(ChannelAccType::channel_type(channel));
}

farsight_maciejseg::ImageConstPtrType
farsight_maciejseg:: 
delete_object(LocationType const& point )
{
  rich_cell::Pointer cell =  maciejseg_.get_cell_at( point );
  maciejseg_.invalidate_cell( cell );
  return this->display();
}

farsight_maciejseg::ImageConstPtrType
farsight_maciejseg:: 
merge_objects(LocationListType const& points)
{
  std::vector<rich_cell::Pointer> cells;
  for (unsigned int i = 0; i<points.size(); i++) {
    rich_cell::Pointer cell =  maciejseg_.get_cell_at( points[i] );
    cells.push_back( cell );
  }
  rich_cell::Pointer cell = maciejseg_.merge_cells( cells );

  return this->display();
}

bool
farsight_maciejseg::
run()
{
 
  if (parameters_.empty()) {
    std::cerr<<"ERROR: No input parameters!"<<std::endl;
    return false;
  }
 
  /* 
  parameters_.resize(14);
  parameters_[0] = "./";
  parameters_[1] = "ID-204_";
  parameters_[2] = "tif";
  parameters_[3] = "1";
  parameters_[4] = "6";
  parameters_[5] = "2";
  parameters_[6] = "500";
  parameters_[7] = "4";
  parameters_[8] = "1";
  parameters_[9] = "1";
  parameters_[10] = "";
  parameters_[11] = "";
  parameters_[12] = "1";
  parameters_[13] = "1";
  maciejseg_.run( parameters_ );
  */

  maciejseg_.run( intensity_image_, parameters_ );
  return true;
}

void 
farsight_maciejseg::
set_parameters( std::string const& xml_filename )
{
  xml_util_read_param( xml_filename, parameters_ );
}

void
farsight_maciejseg::
write_xml(std::string const& file_path,
          std::string const& xml_filename)
{
  xml_util_write( file_path, xml_filename, image_path_,image_name_, maciejseg_ );
}

void
farsight_maciejseg::
read_xml(std::string const& file_path,
         std::string const& xml_filename)
{
  xml_util_read(file_path, xml_filename, image_path_,image_name_, maciejseg_);
}

farsight_maciejseg::ImageConstPtrType 
farsight_maciejseg::
display()
{
  //convert the internal label image to the output image which records
  //if a pixel is a boundary or not.
  std::vector<rich_cell::Pointer> const & cells = maciejseg_.all_cells();

  if (!label_image_) { //allocate the space
    
    vnl_vector_fixed<int,3> const&  image_size = maciejseg_.image_size(); 
    label_image_ = ImageType::New();

    ImageType::IndexType start;
    start[0] = 0;
    start[1] = 0;
    start[2] = 0;
    
    ImageType::SizeType size;
    size[0] = image_size[0]; // size along X
    size[1] = image_size[1]; // size along Y
    size[2] = image_size[2]; // size along Z
    
    ImageType::RegionType region;
    region.SetSize( size );
    region.SetIndex( start );
    label_image_->SetRegions(region);
    label_image_->Allocate();
  }
  label_image_->FillBuffer(0);
  
  for (unsigned int i = 0; i<cells.size(); i++) {
    if (!cells[i]->valid_) continue;

    for (unsigned int b = 0; b<cells[i]->boundary_points_.size(); b++) {
      vnl_vector_fixed< float, 3 > const & pt =  cells[i]->boundary_points_[b];
      ImageType::IndexType pos;
      pos[0] = pt[0];
      pos[1] = pt[1];
      pos[2] = pt[2];
      label_image_->SetPixel(pos, 255);
    }
  }

  return label_image_.GetPointer();
}
