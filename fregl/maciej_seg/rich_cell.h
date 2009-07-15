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

/** @file rich_cell.h
*   @brief rich_cell class
*   This is the class that represents the cells in the image.
*
*   @author Chia-Ling Tsai (2007/09/28)
*/

#ifndef __RICH_CELL_H_
#define __RICH_CELL_H_

#include <string>
#include <vector>
#include "vnl/vnl_vector_fixed.h"
#include "itkLightObject.h"
#include "itkObjectFactory.h"

#include "cell.h"
#include "itkImageRegion.h"
#include "itkRGBPixel.h"
#include "itkImage.h"

#include "itkGrayscaleErodeImageFilter.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h" 
#include "itkMedianImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkGrayscaleMorphologicalOpeningImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkGradientMagnitudeImageFilter.h"

//#include "itkObject.h"

class edit_record
{
public:
  enum status_type {ADDED, DELETED, MERGED, SPLIT, DUP};
  
  status_type      status_;
  std::vector<int> replacements_; //needed for merge and split
  std::string      date_;
  std::string      name_; 
};

class rich_cell: public itk::LightObject
{
public:
  
  typedef  rich_cell                    Self;
  typedef  itk::SmartPointer< Self >         Pointer;
  typedef  itk::SmartPointer< const Self >   ConstPointer;

  itkNewMacro(Self) 

  typedef unsigned char                   PixelType;
  typedef itk::Image< PixelType, 3 >      ImageType;
  typedef itk::Image< float, 3 >          FloatImageType;
  typedef itk::Image< unsigned short, 3 > LabelImageType;
  
  typedef vnl_vector_fixed< float, 3 >  PointType;
  typedef std::vector<PointType>      PointVectorType;
  typedef itk::ImageRegion< 3 >       RegionType;

  //median filter
typedef itk::MedianImageFilter< ImageType,ImageType >  MedianFilterType;

//structuring element for morph operations
typedef itk::BinaryBallStructuringElement<PixelType, 3> StructuringElementType;

//morphological opening filter
typedef itk::GrayscaleMorphologicalOpeningImageFilter< ImageType, ImageType, 
                                                       StructuringElementType > OpenFilterType;

typedef itk::GrayscaleErodeImageFilter< ImageType, ImageType,
                                        StructuringElementType >  ErodeFilterType;

typedef itk::GrayscaleDilateImageFilter< ImageType, ImageType, 
                                         StructuringElementType >  DilateFilterType;

typedef itk::SubtractImageFilter< ImageType, ImageType, ImageType  >  
SubFilterType;

typedef itk::GradientMagnitudeImageFilter<ImageType, FloatImageType > GradMagFilterType;

//iterators
typedef itk::ImageRegionConstIterator< ImageType > ConstRegionIteratorType;
typedef itk::ImageRegionConstIterator< LabelImageType > LabelConstRegionIteratorType;
typedef itk::ImageRegionIterator< ImageType > RegionIteratorType;

  rich_cell();
  virtual ~rich_cell(){}
  
  int volume_;
  double ave_intensity_;	//avg intensity
  PointType center_;	        //center coordinates
  int label_;		        //label
  double texture_;		//texture
  double eccentricity_;		//eccencricity
  double average_radius_;	//average radius
  double nearest_nbr_dist_;     //Nearest neighbor distance
  int nearest_nbr_label_;       //ID of nearest neighbor
  double neuronal_signal_;
  int class_type_;		//class of cell
  double score_;		//score
  double ave_nnbr_dist_;

  /* variables not stored to the output file
  double bound_ints_ratio_;		//boundary intensity ratio
  double ave_bound_gradient_;		//avg boundary gradient
  double convexity_;		//convexity
  double shape_factor_;		//shape factor  
  double bending_energy_;	//bending energy
  double vol_grad_;	//avg grad
  double radius_variance_;	//radius variance
  double percent_nbr_;	//percent pix nex to nbrs

  // More attribute from Gang's version
  double eng_intensity_dist_;     //energy of intensity distribution
  double entropy_intensity_dist_; //entropy of intensity distribution
  double intensity_variation_; //intensity variation
  double num_poles_on_surface_; //number of poles on surface
  double skew_intensity_dist_; //skew of intensity distribution
  double surface_area_; //surface area
  double orientation_;


  double avg_bound_intensity_;	//avg bound intensity
  int major_axis_;		//major axis
  int minor_axis_;		//minor axis
  double vol_grad_var_;	//grad variance
  int num_bound_pixel_;	//# of boundary pix
  double ave_inside_gradient_;		//avg inside gradient
  int total_intensity_;		//total intensity
  double var_rad_;	
  int num_shared_pixel_;	//# of shered pix
  double orientation_;		//orientation
  std::vector<int> neighbors_;
  */

  //bool isMerged_;	//is the cell merged
  //bool isTrain_;	//is the cell a trainig cell
  //double boundGradVar_;	//boudnary gradient variance
  //int depth_;	//plane depth
  //int startPlane_;	//starting plane
  //double bendEng_;	//bending energy

  std::vector<edit_record> edit_history_;
  bool valid_;
  bool dup_; //duplicated if the same nucleus is segmented in another image

  PointVectorType       all_points_;
  PointVectorType       boundary_points_;
  PointVectorType       interior_points_;
  RegionType            bounding_box_;
  
public:
  //: Compute the volume
  void update_volume();

  //: Compute the eccentricity
  //
  // Since a cell always appears shorter in the z-dimension, the minor
  // and major axes are computed only in the projected x-y plane. The
  // major and minor axis are obtained from the eigen-value
  // decomposition of the scatter matrix of the points. The eccentricity
  // is computed as sqrt( 1-(b^2/a^2) ), where a is the max eigen-value
  // and b is the min eigen-value. If the shape is close to circular,
  // the value is 0.
  void update_eccentricity_and_radius();
  
  //: Compute average signal intensity
  //
  //  \param dist_interior Starting distance (<=0) is the interior distance away
  //  from the boundary.  
  //  \param dist_exterior Ending distance (>=0) is the exterior distance away from the boundary. (If dist_exterior > dist_interior, the entire segmented area is considered)
  //
  float compute_average_intensity(ImageType::Pointer image, 
                                  LabelImageType::ConstPointer label_image,
                                  int dist_interior, int dist_exterior);
  
  //: Compute texture
  //
  //  The texture is computed as average interior gradient magnitude of
  //  un-normalized intensities
  void update_texture(ImageType::Pointer image);

  //: Compute the distance to the nearest neighbor
  //
  //  The distance is shortest distance between the boundary of the
  //  cell to the bounary of the nearest cell. The distance is
  //  approximated as the distance between the centers minus the 2
  //  average radii
  void update_nearest_nbr(std::vector<rich_cell::Pointer> const& neighbor);

protected:
  rich_cell(const Self&);         // purposely declared and not defined
  void operator=(const Self&);    // purposely declared and not defined 
};

#endif
