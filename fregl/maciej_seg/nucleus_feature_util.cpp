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

#include "nucleus_feature_util.h"
#include "vnl/vnl_math.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_matrix.h"
#include "vnl/algo/vnl_svd.h"

#include "itkGrayscaleErodeImageFilter.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h" 
#include "itkMedianImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkGrayscaleMorphologicalOpeningImageFilter.h"

#include "rich_cell.h"

static double compute_eccentricity( rich_cell::PointVectorType points);
static void compute_morph_gradient( IntensityImageType::Pointer image,
                                    rich_cell::Pointer cell);

/*
static void compute_texture( IntensityImageType::Pointer image,
                             rich_cell::Pointer cell);
*/

void
compute_volume_related_features( rich_cell::Pointer cell )
{
  //update volume
  cell->volume_ = cell->all_points_.size();

  //update center
  rich_cell::PointType center(0.0);
  for (unsigned int i = 0; i<cell->all_points_.size(); i++) 
    center += cell->all_points_[i];
  cell->center_ = center/(float)cell->all_points_.size();

  //update eccentricity
  double eccentricity = compute_eccentricity(cell->all_points_);
  cell->eccentricity_ = eccentricity;

  //update convexity. It is not computed for now
  cell->convexity_ = 0;


  //update the shape factor: F=B^3/(36*pi*V*V) where B is the number
  //of boundary voxels and V the volumn voxels
  int B = cell->boundary_points_.size();
  int V = cell->all_points_.size();
  cell->shape_factor_ = B*B*B/(36*vnl_math::pi*V*V);

  //update bending_energy. It is not computed for now
  cell->bending_energy_ = 0;

  //update the variance of the radii. A radius is computed as the
  //distance between the center and a boundary point. This is not
  //exactly how Maciej computed this feature, but it should be how it
  //is computed mathematically.
  double mean_radius = 0, rad_var = 0;
  for (unsigned int i = 0; i<cell->boundary_points_.size(); i++) {
    mean_radius += (cell->boundary_points_[i] - cell->center_).magnitude();
  } 
  mean_radius /= cell->boundary_points_.size();
  for (unsigned int i = 0; i<cell->boundary_points_.size(); i++) {
    double leng = (cell->boundary_points_[i] - cell->center_).magnitude();
    rad_var += (leng - mean_radius) * (leng - mean_radius);
  } 
  cell->radius_variance_ = rad_var/cell->boundary_points_.size();
  
}

void 
compute_intensity_related_features( IntensityImageType::Pointer image, 
                                    rich_cell::Pointer cell)
{
  //update ave_intensity_
  double intensity = 0;
  for (unsigned int i = 0; i<cell->all_points_.size(); i++) {
    IntensityImageType::IndexType ind;
    ind[0] = cell->all_points_[i][0];
    ind[1] = cell->all_points_[i][1];
    ind[2] = cell->all_points_[i][2];
    intensity += image->GetPixel(ind);
  }
  cell->ave_intensity_ = intensity/(double)cell->all_points_.size();
  //upate bound_ints_ratio_
  intensity = 0;
  for (unsigned int i = 0; i<cell->boundary_points_.size(); i++) {
    IntensityImageType::IndexType ind;
    ind[0] = cell->all_points_[i][0];
    ind[1] = cell->all_points_[i][1];
    ind[2] = cell->all_points_[i][2];
    intensity += image->GetPixel(ind);
  }
  double ave_bdry_int = intensity/(double)cell->boundary_points_.size();
  cell->bound_ints_ratio_ = cell->ave_intensity_/ave_bdry_int;

  //update ave_bound_gradient_ and vol_grad_
  compute_morph_gradient( image, cell );

  //update texture_. Not computed for now, since the computation is
  //very much like gradient operation

  //compute_texture( image, cell );

}

void 
compute_neighbor_related_features( LabelImageType::Pointer mask, 
                                   rich_cell::Pointer cell)
{
  //update percent_nbr_

  //update neighbors_
}

// Eccentricity is computed differently from the way Maciej did. Since
// a cell always appears shorter in the z-dimension, the minor and
// major axes are computed only in the projected x-y plane. The
// semi-major and semi-minor are obtained from the eigen-value
// decomposition of the scatter matrix of the points.
static double
compute_eccentricity( rich_cell::PointVectorType points)
{
  //Step1: project the points to the x-y plane
  std::vector<vnl_vector_fixed<double, 2> > points_2d;
  std::vector<bool> checked(points.size(), false);
  
  for (unsigned int i = 0; i<points.size(); i++) {
    if (checked[i] = false) {
      //record the x-y position
      vnl_vector_fixed<double, 2> pt(points[i][0],points[i][1]);
      points_2d.push_back( pt );

      //scan through the list to see if other points in the remaining
      //list contains the same x-y position
      for (unsigned int j = 0; j<i; j++) {
        if (points[j][0] == points[i][0] && points[j][1] == points[i][1]) 
          checked[j] = true;
      }
    }
  }

  //Step2: compute the scatter matrix
  //
  // compute the center
  vnl_vector_fixed<double, 2> center(0.0,0.0);
  for (unsigned int i = 0; i<points_2d.size(); i++) {
    center += points_2d[i];
  }
  center /= (double)points_2d.size();

  vnl_matrix<double> cov_matrix(2,2,0.0);
  for (unsigned int i = 0; i<points_2d.size(); i++) {
    cov_matrix += outer_product(points_2d[i]-center, points_2d[i]-center);
  }
  cov_matrix /= (double)points_2d.size();

  //perform eigen-value decomposition to get the semi-major (a) and minor (b) 
  vnl_svd<double> svd_from( cov_matrix );
  
  //eccentricity=sqrt( 1-(b^2/a^2) ). If the shape is close to
  //circular, the value is 0. 
  return vcl_sqrt( 1- (svd_from.W(1)*svd_from.W(1))/(svd_from.W(0)*svd_from.W(0)) );
}

// The image is filtered using median filtering and opened before the
// computation of the morphological gradient
static
void 
compute_morph_gradient( IntensityImageType::Pointer image,
                        rich_cell::Pointer cell)
{
  image->SetRequestedRegion( cell->bounding_box_ );
  
  int radius = 1; //for both median filtering and structuring element
  const unsigned int Dim = 3;

  //median filter
  typedef itk::MedianImageFilter<
  IntensityImageType, InternalImageType >  MedianFilterType;

  //structuring element for morph operations
  typedef itk::BinaryBallStructuringElement< 
    InternalPixelType, Dim >        StructuringElementType;
  
  //morphological opening filter
  typedef itk::GrayscaleMorphologicalOpeningImageFilter<
    InternalImageType, InternalImageType, 
    StructuringElementType > OpenFilterType;
  
  typedef itk::GrayscaleErodeImageFilter<
    InternalImageType, InternalImageType,
    StructuringElementType >  ErodeFilterType;
  
  typedef itk::GrayscaleDilateImageFilter<
    InternalImageType, InternalImageType, 
    StructuringElementType >  DilateFilterType;
  
  typedef itk::SubtractImageFilter<
    InternalImageType,InternalImageType,
    InternalImageType  >              SubFilterType;

  //set neighborhood
  IntensityImageType::SizeType miIndexRadius;
  
  miIndexRadius[0] = radius; // radius along x
  miIndexRadius[1] = radius; // radius along y
  miIndexRadius[2] = 0;	// radius along z
  MedianFilterType::Pointer f_med = MedianFilterType::New();
  
  //set radius and input
  f_med->SetRadius( miIndexRadius );
  f_med->SetInput( image );
  
  OpenFilterType::Pointer f_open = OpenFilterType::New();
  ErodeFilterType::Pointer  f_erode  = ErodeFilterType::New();
  DilateFilterType::Pointer f_dilate = DilateFilterType::New();
  SubFilterType::Pointer f_sub = SubFilterType::New();
  
  StructuringElementType  structuringElement;
  structuringElement.SetRadius( radius );
  structuringElement.CreateStructuringElement();
  
  f_open->SetKernel( structuringElement );
  f_erode->SetKernel(  structuringElement );
  f_dilate->SetKernel( structuringElement );

  //connect open to medium filter
  f_open->SetInput( f_med->GetOutput() );
  f_erode->SetInput(  f_open->GetOutput() );
  f_dilate->SetInput( f_open->GetOutput() ); 

  //set inputs of substract filter for morphological gradient
  //(dilation-erosion=morphological gradient)
  f_sub->SetInput1( f_dilate->GetOutput() );
  f_sub->SetInput2( f_erode->GetOutput() );

  try
    {
      f_sub->Update();
    }
  catch (itk::ExceptionObject & e)
    {
      std::cerr << "Exception in SubFilter: " << e << std::endl;
      exit(0);
    }

  // Now we can compute the gradient values!
  
  InternalImageType::Pointer morph_image = f_sub->GetOutput();
  double vol_grad = 0;
  IntensityImageType::IndexType ind;
  for (unsigned int i = 0; i<cell->all_points_.size(); i++) {
    ind[0] = cell->interior_points_[i][0];
    ind[1] = cell->interior_points_[i][1];
    ind[2] = cell->interior_points_[i][2];
    vol_grad += morph_image->GetPixel( ind );
  }
  cell->vol_grad_ = vol_grad/(double)cell->all_points_.size();

  double bdry_grad = 0;
  for (unsigned int i = 0; i<cell->boundary_points_.size(); i++) {
    ind[0] = cell->interior_points_[i][0];
    ind[1] = cell->interior_points_[i][1];
    ind[2] = cell->interior_points_[i][2];
    bdry_grad += morph_image->GetPixel( ind );
  }
  cell->ave_bound_gradient_ = vol_grad/(double)cell->boundary_points_.size();
}

/*
static void 
compute_bounary_shared( LabelImageType::Pointer image,
                        rich_cell::Pointer cell)
{
  typedef itk::ConstantBoundaryCondition<InternalLabelImageType> boundaryConditionType;
  typedef itk::ConstNeighborhoodIterator<InternalLabelImageType, boundaryConditionType > NeighborhoodIteratorType;

  // The offsets for the neighboring pixels
  NeighborhoodIteratorType::OffsetType offset1 = {-1, 0, 0};
  NeighborhoodIteratorType::OffsetType offset2 = { 1, 0, 0};
  NeighborhoodIteratorType::OffsetType offset3 = { 0,-1, 0};
  NeighborhoodIteratorType::OffsetType offset4 = { 0, 1, 0};
  
  NeighborhoodIteratorType::RadiusType radius;
  radius.Fill(1);
  image->SetRequestedRegion(cell->bounding_box_);
  NeighborhoodIteratorType it( radius, image,
                               image->GetRequestedRegion());
  double nbrpix = 0;
  LabelImageType::Indextype ind;
  for (unsigned int i = 0; i<cell->boundary_points_(); i++) {
    
      if ( it.GetCenterPixel() == cell->label_ ) {
        texture += vcl_sqr(
      }
  }
  }
*/
