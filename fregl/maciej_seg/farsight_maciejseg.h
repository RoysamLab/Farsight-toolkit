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
