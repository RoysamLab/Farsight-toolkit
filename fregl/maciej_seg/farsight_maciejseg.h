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

#ifndef _FARSIGHT_MACIEJSEG_
#define _FARSIGHT_MACIEJSEG_

//:
// \file
// \brief  The class serve as the interface between the GUI and maciej segmentation module. 
// \author Charlene Tsai
// \date   Nov 2007

#include <Common/farsight_object.h>
#include <maciej_seg/maciejSegmentation.h>

class farsight_maciejseg: public farsight_object<3>
{
public:
  //: constructor
  farsight_maciejseg(){};
  
  //: destructor
  virtual ~farsight_maciejseg(){};

  //: Set the path and name of the input image.
  //
  //  The image is loaded in this function. Please make sure the
  //  parameters are already set, since it needs to know which color
  //  channel to load.
  virtual void set_image( std::string const& image_path,
                          std::string const& image_name);

  //: Delete and object
  virtual ImageConstPtrType delete_object(LocationType const& point);

  //: Merge objects
  virtual ImageConstPtrType merge_objects(LocationListType const& points );

  //: Trigger the computation
  virtual bool run();

  //: Set the input parameters from an xml file
  virtual void set_parameters( std::string const& xml_filename );

  //: Read in the result from an xml file
  virtual void read_xml(std::string const& file_path,
                        std::string const& filename);

  //: Output the result to an xml file
  virtual void write_xml(std::string const& file_path,
                         std::string const& xml_filename);

  //: Return point to label image
  virtual ImageConstPtrType display();


protected:
  maciejSegmentation      maciejseg_;
  ImagePtrType            intensity_image_;
};

#endif
