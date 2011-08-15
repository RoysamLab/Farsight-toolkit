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

#ifndef coarse2fine_register_h_
#define coarse2fine_register_h_

#include "itkImage.h"

#include <rgrl/rgrl_object.h>
#include <rgrl/rgrl_view_sptr.h>
#include <rgrl/rgrl_feature_sptr.h>
#include <rgrl/rgrl_feature_set_sptr.h>
#include <rgrl/rgrl_converge_status_sptr.h>
#include <rgrl/rgrl_transformation_sptr.h>

#include <rrl/rrl_pair_reg_params_sptr.h>
#include <rrl/rrl_view_based_reg_sptr.h>
#include <rrl/rrl_reg_info_sptr.h>

#include <rkpl/rkpl_features_for_registration.h>
#include <rkpl/rkpl_keypoint.h>

class fregl_coarse2fine_register: public rgrl_object {
public:
  typedef float                             InternalPixelType;
  typedef itk::Image<InternalPixelType, 3 > InternalImageType;
  typedef rkpl_features_for_registration<3> FeatureType;
  typedef vcl_vector< rkpl_keypoint<3> >    keypoint_vector;
  typedef vcl_vector< rgrl_feature_sptr >   feature_vector;
  typedef vcl_vector< rgrl_feature_set_sptr > feature_set_vector;

  fregl_coarse2fine_register( rrl_pair_reg_params_sptr reg_params, 
                              InternalImageType::Pointer from_image, 
                              InternalImageType::Pointer to_image);
  ~fregl_coarse2fine_register(){}

  bool run(rgrl_transformation_sptr  init_xform);

  // Defines type-related functions
  rgrl_type_macro( rgrl_transformation, rgrl_object );
private:
  void rkpl_corners_to_rgrl_features( keypoint_vector const& rkpl_corners,
                                      feature_vector & features );
  void rkpl_faces_to_rgrl_features( keypoint_vector const& rkpl_faces,
                                    feature_vector & features );

  //: Concatenate features
  void make_vectors_of_feature_sets( FeatureType::Pointer features,
                                     feature_set_vector & concat_matchable_corners,
                                     feature_set_vector & concat_driving_corners,
                                     feature_set_vector & concat_matchable_faces,
                                     feature_set_vector & concat_driving_faces,
                                     double min_scale,
                                     double max_scale,
                                     double upper_min_feature_scale);

  
private: 
  // Matches
  feature_set_vector concat_moving_matchable_corners_;
  feature_set_vector concat_moving_driving_corners_;
  feature_set_vector concat_moving_matchable_faces_;
  feature_set_vector concat_moving_driving_faces_;
  feature_set_vector concat_fixed_matchable_corners_;
  feature_set_vector concat_fixed_driving_corners_;
  feature_set_vector concat_fixed_matchable_faces_;
  feature_set_vector concat_fixed_driving_faces_;

  // registration engine
  rrl_pair_reg_params_sptr params_; 
  rrl_view_based_reg_sptr registration_ptr_;

  // reg result
  rrl_reg_info_sptr final_reg_result_;

  // current best view and status
  rgrl_view_sptr             init_view_;
  rgrl_view_sptr             current_view_;
  rgrl_converge_status_sptr  current_status_;
};
#endif
