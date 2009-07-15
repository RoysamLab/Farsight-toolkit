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

#ifndef _MACIEJSEGMENTATION_H_
#define _MACIEJSEGMENTATION_H_

#include <vector>
#include <string>

#include "filter.h"
#include "rich_cell.h"
//#include "itkLightObject.h"

#include <vcl_memory.h>
#include <rsdl/rsdl_bins.h>
#include <vbl/vbl_ref_count.h>
#include <vbl/vbl_smart_ptr.h>

class CastPixelAccessor 
{
public:
  typedef itk::RGBPixel<unsigned char>   InternalType;
  typedef itk::RGBPixel<float>           ExternalType;
  
  static void Set(InternalType& output, const ExternalType & input )
  {
    output = static_cast<InternalType>( input );
  }
  static ExternalType Get( const InternalType & input ) 
  {
    return static_cast<ExternalType>( input );
  }
};

class maciejSegmentation: public vbl_ref_count
{
public:
  typedef vbl_smart_ptr< maciejSegmentation >  Pointer;
  typedef itk::RGBPixel<unsigned char>      pixelType;
  typedef itk::Image< pixelType,  3 >       InputLabelImageType;
  typedef itk::Image< unsigned short,3 >    InternalLabelImageType;
  typedef itk::Image< unsigned char,3 >     InternalIntensityImageType;

  typedef itk::ConstantBoundaryCondition<InternalLabelImageType> boundaryConditionType;
  typedef itk::ConstNeighborhoodIterator<InternalLabelImageType, boundaryConditionType > NeighborhoodIteratorType;

  maciejSegmentation();
  ~maciejSegmentation(){};

  //: Execute the segmentation task 
  void run(std::vector<std::string> const& in_parameters,
           filter::mInputImageFileType::Pointer & proj_img,
           filter::mInputImageFileType::Pointer & seg_result);

  void run(std::vector<std::string> const& in_parameters);

  //: The function takes an itkImage<unsigned char,3> as an input.
  //
  //  The parameter list contains:
  //  parameter[0] = channel (0-red, 1-green, 2-blue)
  //  parameter[1] = pT (default to 500)
  //  parameter[2] = g (default to 4)
  //  parameter[3] = rM (default to 1)
  //  parameter[4] = rG (default to 1)
  //  parameter[5] = mode (default to 1 for autotrain)
  //  parameter[6] = training file for glia (can be ommited)
  //  parameter[7] = training file for neuron (can be ommited)
  void run(InternalIntensityImageType::Pointer image,
           std::vector<std::string> const& in_parameters);

  //: Return the cell which contains the pixel given by pos
  rich_cell::Pointer get_cell_at(vnl_vector_fixed<float,3> pos) const;

  //: Invalidate the cell and update the label_image
  void invalidate_cell( rich_cell::Pointer cell );

  //: Merge two cells and create a new replacement
  rich_cell::Pointer merge_cells( rich_cell::Pointer cell1, 
                                  rich_cell::Pointer cell2 );
  rich_cell::Pointer merge_cells( std::vector<rich_cell::Pointer> cell1s );

  //: Return a const reference of the vector of cells. 
  std::vector<rich_cell::Pointer> const & all_cells() const; 

  //: Return a const reference of the vector of cells. 
  void shift_cells( vnl_vector_fixed<float,3> shift ); 

  //: Return the image size of the segmented image 
  vnl_vector_fixed<int,3>  image_size() const; 

  //: Return parameters
  std::vector<std::string> const& parameters() const;

  //: Add one cell at a time
  void add_cell( rich_cell::Pointer new_cell );

  //: Set parameters
  void set_parameters(std::vector<std::string> const& parameters);

  //: Given the label_image recorded in RGB format, pixels of the cells are updated.
  //
  //  Cells are only considered valid, if the label_image and the cell
  //  pixels are synchronized.
  void update_cell_pixels(InputLabelImageType::Pointer label_image);

  //: Given the label_image recorded in binary file, pixels of the cell are updated.
  void update_cell_pixels(std::string const & label_image_filename, int size_x, int size_y, int size_z);

   //: Given the label_image recorded an 1D array(in the order of ZXY), pixels of the cell are updated.
  void update_cell_pixels(unsigned short * memblock, int size_x, int size_y, int size_z);

  //: Update the cell pixels only for one cell
  //
  //  This function updates all_points_, boundary_points_, and
  //  interior_points_
  void update_cell_pixels_for_one(rich_cell::Pointer cell);

  //: Dump the label image out as an RGB image
  void output_cell_pixels(std::string const & filename) const;

  //: Return the label image
  InternalLabelImageType::ConstPointer label_image() const;

  /************ functions for queries *************************/

  //: construct the 3D binning structure
  //
  //  The bin is constructed from objects of the same class. When
  //  class==0 it implies no classification, and all objects are
  //  considered. If dup_removed is set to true, cells marked as dup
  //  are not considered.
  void construct_bins(double bin_size = 20, int class_label=0, bool dup_removed = false);

  //:  Return the k nearest cells based on Euclidean distance.
  //
  //   k is default to 1, which is equivalent to nearest_features
  void k_nearest_cells( std::vector<rich_cell::Pointer> & results, 
                        vnl_vector_fixed<float,3> const& center, 
                        unsigned int k=1 ) const;

  //: Return all cells within a given Euclidean distance
  void cells_within_radius( std::vector<rich_cell::Pointer> & results, 
                            vnl_vector_fixed<float,3> const& center,
                            double radius ) const;

  //: Create the label image for valid cells only
  void create_label_image(vnl_vector_fixed<int,3> const& image_size);

private:
  void convert_cells_to_rich_cells(std::vector<cell*> const &cells,
                                   vnl_vector_fixed<int,3> const& image_size);


  //: Update the label_image only for one cell using all_points_
  void update_label_image( rich_cell::Pointer cell );

  //: Update the features of the given cell
  void update_features( rich_cell::Pointer cell );

  std::string                          image_id_;
  std::vector<rich_cell::Pointer>      rich_cells_;
  //InternalIntensityImageType::Pointer  intensity_image_;
  InternalLabelImageType::Pointer      label_image_; //recording the cell ID that a pixel belongs to
  std::vector<std::string>             parameters_;
  bool                                 cells_valid_;
  int                                  cell_count_;

  // data members for 3d binning
  typedef rsdl_bins<3,float,rich_cell::Pointer> bin_type;
  vcl_auto_ptr< bin_type > bins_;
  bool bins_set_;

};

#endif
