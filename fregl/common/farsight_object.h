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

#ifndef _FARSIGHT_OBJECT_
#define _FARSIGHT_OBJECT_

//:
// \file
// \brief  The class is the base class which defines the interface between the GUI and the computation module. 
// \author Charlene Tsai
// \date   Nov 2007

#include "itkImage.h"
//#include "itkLightObject.h"
//#include "itkObjectFactory.h"

#include <string>
#include <vector>
#include "vnl/vnl_vector_fixed.h"

//: The class is templated over the dimension of the output image. 
template<unsigned dim>
class farsight_object 
{
public:
  /*
  typedef  farsight_object                  Self;
  typedef  itk::SmartPointer< Self >        Pointer;
  typedef  itk::SmartPointer< const Self >  ConstPointer;

  itkNewMacro(Self) 
  */

  typedef unsigned char                         PixelType;
  typedef typename itk::Image< PixelType, dim > ImageType;
  typedef typename ImageType::Pointer           ImagePtrType;
  typedef typename ImageType::ConstPointer      ImageConstPtrType;
  typedef float                                 CoordType;
  typedef vnl_vector<CoordType>                 LocationType;
  typedef std::vector< LocationType >           LocationListType;


  //: constructor
  //
  //  Does nothing
  farsight_object(){};

  virtual ~farsight_object(){};

  //: A very generic function which provide the module the mouse location
  virtual ImageConstPtrType mouse_click(LocationType const& point);

  //: Set the path and name of the input image.
  virtual void set_image( std::string const& image_path,
                          std::string const& image_name);
  //: Get image_name
  virtual std::string const& get_image_name();

  //: Get image_path
  virtual std::string const& get_image_path();

  //: Add an object
  virtual ImageConstPtrType add_object(LocationListType const& points);
  virtual ImageConstPtrType add_object(LocationType const& point);

  //: Delete and object
  //
  // The object which contains the given point is removed/invalidated
  virtual ImageConstPtrType delete_object(LocationType const& point);
  virtual ImageConstPtrType delete_object(LocationListType const& point);

  //: Merge objects
  //
  //  Objects which contains the given set of points are merged
  virtual ImageConstPtrType merge_objects(LocationListType const& points );

  //: Split an object
  //
  // Objects which contains the given set of points are divided
  virtual ImageConstPtrType split_object(LocationType const& points );
  virtual ImageConstPtrType split_object(LocationListType const& points );

  //: Return the pointer to the label image for display
  virtual ImageConstPtrType display() ;

  //: Trigger the computation
  virtual bool run() = 0;

  //: Set the input parameters
  virtual void set_parameters( std::vector< std::string > const& parameters);

  //: Set the input parameters from an xml file
  virtual void set_parameters( std::string const& xml_filename );

  //: Generate the default input parameters
  virtual void set_default_parameters();

  //: Output the result (with the parameters) to an xml file
  virtual void write_xml(std::string const& file_path,
                         std::string const& xml_filename) = 0;

  //: Read in the result (with the parameters) from an xml file. 
  virtual void read_xml(std::string const& file_path,
                        std::string const& xml_filename) = 0;

protected:
  std::vector< std::string > parameters_;
  ImagePtrType               label_image_;
  std::string                image_name_;
  std::string                image_path_;
};

#endif
