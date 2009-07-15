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

#include "rich_cell.h"

#include <algorithm>

rich_cell::
rich_cell()
  : valid_(false), dup_(false)
{
  volume_ = 0;
  ave_intensity_ = 0;	    //avg intensity
  label_ = 0;		    //label
  texture_ = 0;		    //texture
  eccentricity_ = 0;	    //eccencricity
  average_radius_ = 0;
  neuronal_signal_ = 0;
  class_type_ = 0;	    //class of cell
  nearest_nbr_dist_ = 0;
  score_ = 0;		    //score
  nearest_nbr_label_ = 0;
  ave_nnbr_dist_ = 0;

  /*
  convexity_ = 0;	    //convexity
  shape_factor_ = 0;	    //shape factor
  bending_energy_ = 0;	    //bending energy
  vol_grad_ = 0;	    //avg grad
  radius_variance_ = 0;	    //radius variance
  percent_nbr_ = 0;	    //percent pix nex to nbrs
  ave_nbr_dist_ = 0;           //Average n-nearest neighbor distance 
  ave_bound_gradient_ = 0;  //avg boundary gradient
  bound_ints_ratio_ = 0;    //boundary intensity ratio
  */

  /*
  // More attribute from Gang's version
  eng_intensity_dist_ = 0;     //energy of intensity distribution
  entropy_intensity_dist_ = 0; //entropy of intensity distribution
  intensity_variation_ = 0;    //intensity variation
  num_poles_on_surface_ = 0;   //number of poles on surface
  skew_intensity_dist_ = 0;    //skew of intensity distribution
  surface_area_ = 0;           //surface area
  orientation_ = 0;
  */
}

void
rich_cell::
update_volume()
{
  volume_ =  all_points_.size();
}

void 
rich_cell::
update_eccentricity_and_radius()
{
  //Step1: project the points to the x-y plane
  std::vector<vnl_vector_fixed<double, 2> > points_2d;
  std::vector<bool> checked(all_points_.size(), false);
  
  for (unsigned int i = 0; i<all_points_.size(); i++) {
    if (!checked[i]) {
      //record the x-y position
      vnl_vector_fixed<double, 2> pt(all_points_[i][0],all_points_[i][1]);
      points_2d.push_back( pt );
      checked[i] = true;

      //scan through the list to see if other points in the remaining
      //list contains the same x-y position
      for (unsigned int j = i+1; j<all_points_.size(); j++) {
        if (all_points_[j][0] == all_points_[i][0] 
            && all_points_[j][1] == all_points_[i][1]) 
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
  std::cout<<"points_2d.size = "<<points_2d.size()<<std::endl;
  for (unsigned int i = 0; i<points_2d.size(); i++) {
    cov_matrix += outer_product(points_2d[i]-center, points_2d[i]-center);
  }
  cov_matrix /= (double)points_2d.size();

  //perform eigen-value decomposition to get the semi-major (a) and minor (b) 
  vnl_svd<double> svd_from( cov_matrix );
  
  //eccentricity=sqrt( 1-(b^2/a^2) ). If the shape is close to
  //circular, the value is 0. 
  //eccentricity_ =vcl_sqrt( 1- (svd_from.W(1)*svd_from.W(1))/(svd_from.W(0)*svd_from.W(0)) );
  eccentricity_ =vcl_sqrt( 1- vnl_math_abs(svd_from.W(1)/svd_from.W(0)) );
  float long_axis_mag = vcl_sqrt(vnl_math_abs(svd_from.W(0)));
  float short_axis_mag = vcl_sqrt(vnl_math_abs(svd_from.W(1)));
  average_radius_ = 0.5*(long_axis_mag + short_axis_mag);
}

float
rich_cell::
compute_average_intensity( ImageType::Pointer intensity_image, LabelImageType::ConstPointer label_image, int dist_interior, int dist_exterior)
{
  float sum_interior = 0;
  float sum_exterior = 0;
  int count_interior = 0;
  int count_exterior = 0;

  if (dist_exterior < dist_interior) { // the entire segmented area is taken
    for (unsigned int b = 0; b<all_points_.size(); b++) {
      vnl_vector_fixed< float, 3 > const & pt = all_points_[b];
      ImageType::IndexType pos;
      pos[0] = pt[0];
      pos[1] = pt[1];
      pos[2] = pt[2];
      sum_interior += intensity_image->GetPixel(pos);
    }

    return sum_interior/all_points_.size();
  }
   
  RegionType region = bounding_box_;
  if (dist_interior < 0) { //erode the mask 
    // Generate a mask image of the cell region. Erode the region by
    // r_interior
    RegionType::SizeType size = region.GetSize();
    RegionType::IndexType start={{0,0,0}};
    ImageType::Pointer cropped_mask = ImageType::New(); 
    RegionType mask_region;
    mask_region.SetIndex( start );
    mask_region.SetSize( size );
    cropped_mask->SetRegions( mask_region );
    cropped_mask->Allocate();
    cropped_mask->FillBuffer(0);
    LabelConstRegionIteratorType it1( label_image, region);
    RegionIteratorType it2( cropped_mask, mask_region );
    for (it1.GoToBegin(), it2.GoToBegin(); !it1.IsAtEnd(); ++it1, ++it2) {
      if (it1.Get() == label_)
        it2.Set( 255 );
    }
      
    ImageType::Pointer eroded_mask;
    ErodeFilterType::Pointer f_erode = ErodeFilterType::New();
    SubFilterType::Pointer f_sub = SubFilterType::New();
    StructuringElementType  structuringElement;
    structuringElement.SetRadius( -dist_interior );
    structuringElement.CreateStructuringElement();
    f_erode->SetKernel( structuringElement );
    f_erode->SetInput(cropped_mask);
    f_sub->SetInput1( cropped_mask  );
    f_sub->SetInput2( f_erode->GetOutput() );
    try {
      f_sub->Update();
    }
    catch (itk::ExceptionObject & e) {
      std::cerr << "Exception in SubFilter: " << e << std::endl;
      exit(0);
    }
    eroded_mask = f_sub->GetOutput();
      
    // Sum the signal in the eroded region only
    ConstRegionIteratorType it3( eroded_mask, mask_region );
    ConstRegionIteratorType it4( intensity_image, region);
    for (it3.GoToBegin(), it4.GoToBegin(); !it3.IsAtEnd(); ++it1, ++it3, ++it4) {
      if (it3.Get() > 0) {
        sum_interior += it4.Get();
        count_interior ++;
      }
    }
  }
  if (dist_exterior > 0) { //dilate the mask
    // enlarge the bounding box by r on each side.
    RegionType::SizeType image_size = intensity_image->GetLargestPossibleRegion().GetSize();
    RegionType::SizeType size = region.GetSize();
    RegionType::IndexType start = region.GetIndex();
    RegionType::IndexType end;
    end[0] = vnl_math_min(start[0]+size[0]+dist_exterior, image_size[0]);
    end[1] = vnl_math_min(start[1]+size[1]+dist_exterior, image_size[1]);
    end[2] = vnl_math_min(start[2]+size[2]+dist_exterior, image_size[2]);
    start[0] = vnl_math_max(int(start[0]-dist_exterior), 0);
    start[1] = vnl_math_max(int(start[1]-dist_exterior), 0);
    start[2] = vnl_math_max(int(start[2]-dist_exterior), 0);
    
    size[0] = end[0] - start[0];
    size[1] = end[1] - start[1];
    size[2] = end[2] - start[2];
    region.SetSize( size );
    region.SetIndex( start );
    
    // Generate a mask image of the region just found. Dilate the
    // region defined by the segmentation by r. 
    ImageType::Pointer cropped_mask = ImageType::New(); 
    RegionType mask_region;
    start[0] = start[1] = start[2] = 0;
    mask_region.SetIndex( start );
    mask_region.SetSize( size );
    cropped_mask->SetRegions( mask_region );
    cropped_mask->Allocate();
    cropped_mask->FillBuffer(0);
    LabelConstRegionIteratorType it1( label_image, region);
    RegionIteratorType it2( cropped_mask, mask_region );
    for (it1.GoToBegin(), it2.GoToBegin(); !it1.IsAtEnd(); ++it1, ++it2) {
      if (it1.Get() == label_)
        it2.Set( 255 );
    }
    ImageType::Pointer dilated_mask;
    DilateFilterType::Pointer f_dilate = DilateFilterType::New();
    SubFilterType::Pointer f_sub = SubFilterType::New();
    StructuringElementType  structuringElement;
    structuringElement.SetRadius( dist_exterior );
    structuringElement.CreateStructuringElement();
    f_dilate->SetKernel( structuringElement );
    f_dilate->SetInput(cropped_mask);
    f_sub->SetInput1( f_dilate->GetOutput() );
    f_sub->SetInput2( cropped_mask );
    
    try {
      f_sub->Update();
    }
    catch (itk::ExceptionObject & e) {
      std::cerr << "Exception in SubFilter: " << e << std::endl;
      exit(0);
    }
    dilated_mask = f_sub->GetOutput();
    
    // Sum the signal in the dilated region only
    ConstRegionIteratorType it3( dilated_mask, mask_region );
    ConstRegionIteratorType it4( intensity_image, region);
    for (it3.GoToBegin(), it4.GoToBegin(); !it3.IsAtEnd(); ++it1, ++it3, ++it4) {
      if (it3.Get() > 0) {
        sum_exterior += it4.Get();
        count_exterior ++;
      }
    }
  }
  
  // average the interior and exterior signals
  return (sum_interior+sum_exterior)/float(count_interior+count_exterior);
}

void 
rich_cell::
update_texture( ImageType::Pointer image )
{
  // Compute the gradient magnitude
  GradMagFilterType::Pointer gmFilter = GradMagFilterType::New();
  gmFilter->SetInput( image );
  try {
    gmFilter->Update();
  }
  catch (itk::ExceptionObject & e) {
      std::cerr << "Exception in SubFilter: " << e << std::endl;
      exit(0);
  }
  FloatImageType::Pointer gmImage = gmFilter->GetOutput();

  // Compute the texture as the average gradient of the first portion
  // of the sorted gradient magnitudes.

  std::cout<<"class = "<<class_type_<<std::endl;
  double interior_grad = 0;
  std::vector< float > grad_magnitudes;
  grad_magnitudes.reserve(interior_points_.size());
  FloatImageType::IndexType ind;
  for (unsigned int i = 0; i<interior_points_.size(); i++) {
    ind[0] = interior_points_[i][0];
    ind[1] = interior_points_[i][1];
    ind[2] = interior_points_[i][2];
    //interior_grad += gmImage->GetPixel( ind );
    //std::cout<<interior_points_[i]<<", grad_mag = "<<gmImage->GetPixel( ind )<<std::endl;
    grad_magnitudes.push_back(gmImage->GetPixel( ind ));
  }
  
  /*
  std::vector<float>::iterator loc = grad_magnitudes.begin()+interior_points_.size()/2;
  std::nth_element(grad_magnitudes.begin(), loc, grad_magnitudes.end());
  //texture_ = interior_grad/(double)interior_points_.size();
  */
  std::vector<float>::iterator loc = grad_magnitudes.begin()+grad_magnitudes.size()/2;
  std::partial_sort(grad_magnitudes.begin(), loc, grad_magnitudes.end());
  for (std::vector<float>::iterator itr=grad_magnitudes.begin(); itr<loc; itr++) {
    interior_grad += *itr;
  }
  texture_ = interior_grad/grad_magnitudes.size()*2;
  //std::cout<<"texture = "<<texture_<<std::endl;

  /*
  image->SetRequestedRegion( bounding_box_ );
  
  int radius = 1; //for both median filtering and structuring element

  //set neighborhood
  ImageType::SizeType miIndexRadius;
  
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
  ImageType::Pointer morph_image = f_sub->GetOutput();
  double interior_grad = 0;
  ImageType::IndexType ind;
  for (unsigned int i = 0; i<interior_points_.size(); i++) {
    ind[0] = interior_points_[i][0];
    ind[1] = interior_points_[i][1];
    ind[2] = interior_points_[i][2];
    interior_grad += morph_image->GetPixel( ind );
  }
  texture_ = interior_grad/(double)interior_points_.size();
  */

  /*
  double bdry_grad = 0;
  for (unsigned int i = 0; i<cell->boundary_points_.size(); i++) {
    ind[0] = cell->interior_points_[i][0];
    ind[1] = cell->interior_points_[i][1];
    ind[2] = cell->interior_points_[i][2];
    bdry_grad += morph_image->GetPixel( ind );
  }
  cell->ave_bound_gradient_ = vol_grad/(double)cell->boundary_points_.size();
  */
}

void 
rich_cell::
update_nearest_nbr(std::vector<rich_cell::Pointer> const& neighbors)
{
  for (unsigned i = 0; i<neighbors.size(); i++) {
    if (neighbors[i]->label_ == this->label_) continue;

    assert( !this->average_radius_ && !neighbors[i]->average_radius_);
    float center_dist = (neighbors[i]->center_ - this->center_).magnitude();
    this->nearest_nbr_dist_ = center_dist - this->average_radius_ - neighbors[i]->average_radius_;
    this->nearest_nbr_label_ = neighbors[i]->label_;
  }
}
