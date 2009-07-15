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

#ifndef fnsl3d_filter_h_
#define fnsl3d_filter_h_
//:
// \file
// \brief A class containing functions for preprocessing filters
// \author Charlene Tsai
// \date 07 September 2007
//

#include "itkImage.h"
#include "itkGrayscaleErodeImageFilter.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h" 
#include "itkMedianImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkGrayscaleMorphologicalOpeningImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"

//: A class containing functions for preprocessing filters
//
//  The code is almost a direct porting of the work done by Maciej
//  Wojton, since the module will soon be replaced by better
//  segmentation algorithms to be developed.
class fnsl3d_filter
{
public:

  /************* typedefs *************************/
  static const unsigned int Dim = 3;
  typedef unsigned char                         InputPixelType;
  typedef itk::Image<InputPixelType, Dim>       InputImageType;
  typedef float                                 PixelType;
  typedef itk::Image<PixelType, Dim>            ImageType;

  //medium filter
  typedef itk::MedianImageFilter<
  ImageType, ImageType >  MedianFilterType;
  
  //structuring element for open filter
  typedef itk::BinaryBallStructuringElement< 
    PixelType, Dim >        StructuringElementType;
  
  //morphological opening filter
  typedef itk::GrayscaleMorphologicalOpeningImageFilter<
    ImageType, ImageType, 
    StructuringElementType > OpenFilterType;
  
  //filter types for erode, dialate, open and subtract filters 
  typedef itk::GrayscaleErodeImageFilter<
    ImageType, ImageType,
    StructuringElementType >  ErodeFilterType;
  
  typedef itk::GrayscaleDilateImageFilter<
    ImageType, ImageType, 
    StructuringElementType >  DilateFilterType;
  
  typedef itk::SubtractImageFilter<
    ImageType,ImageType,
    ImageType  >              SubFilterType;
  
  //region of interest fitler
  typedef itk::RegionOfInterestImageFilter< ImageType, 
                                            ImageType > RegionFilterType;

  //binary threshold filter
  typedef itk::BinaryThresholdImageFilter< ImageType, 
                                           ImageType >  ThreshFilterType;

  //create distance map filter
  typedef itk::DanielssonDistanceMapImageFilter< ImageType, 
                                                 ImageType >  DanielFilterType;

  //create gaussian filter
  typedef itk::DiscreteGaussianImageFilter< ImageType,
                                            ImageType> GaussFilterType;

  //declare iterator types
  typedef itk::ImageRegionConstIteratorWithIndex< ImageType > ConstIndexIteratorType;
  typedef itk::ImageRegionIteratorWithIndex< ImageType > IndexIteratorType;
  typedef itk::ImageRegionConstIterator< ImageType > ConstIteratorType;
  typedef itk::ImageRegionIterator< ImageType > IteratorType;
  
  /******************* functions ****************************/
  //: Constructor
  fnsl3d_filter();
  
  //: Destructor
  virtual ~fnsl3d_filter();
  
  //: The main function to trigger all the (fixed set of) processes
  void run_filter(InputImageType::Pointer img, int rad_med = 1, 
                  int rad_grad = 1, int gridX = 1, int gridY = 1, 
                  double paramT = 0.4, double sigma = 0.5);

  //: Get the intensity image
  ImageType::Pointer get_intensity(){ return imgpIntenisty_; }

  //: Get the filtered image
  ImageType::Pointer get_filt(){ return imgpFilt_; }
  
  //: Get gradient image
  ImageType::Pointer get_gradient(){ return imgpGradient_; }

  //: Get image size
  int* get_size(void);

  // For debugging
  void save_image(ImageType::Pointer img, std::string file);

 private:
  ImageType::Pointer median_open_image(ImageType::Pointer img,int radius);

  void get_min_max_mean_std(ImageType::Pointer image, double &max, 
                            double &min, double &mean, double &std, 
                            bool std_flag);

  ImageType::Pointer morph_gradient(ImageType::Pointer img,int radius);

  ImageType::Pointer threshold_distance_map(ImageType::Pointer img,int 
                                            hor,int ver,int dep, double param);
  ImageType::Pointer combine_invert_distance(ImageType::Pointer imgDaniel, 
                                             ImageType::Pointer imgMorphGrad,
                                             double sigma);
  
private:

  ImageType::RegionType                    region_;
  int                                      image_size_[3];

  //std::vector to hold miIndex of background pixels
  std::vector<ImageType::IndexType> vIndex_;

  //temp images
  ImageType::Pointer imgpIntenisty_;
  ImageType::Pointer imgpGradient_;
  ImageType::Pointer imgpFilt_;
};
#endif
