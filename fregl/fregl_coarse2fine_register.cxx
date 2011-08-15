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

#include "fregl_coarse2fine_register.h"

#include <rrl/rrl_convergence_on_weighted_nas.h>
#include <rrl/rrl_pair_reg_params.h>
#include <rrl/rrl_view_generator_sptr.h> 
#include <rrl/rrl_model_selection_aicc.h>
#include <rrl/rrl_model_selection_criterion_sptr.h>
#include <rrl/rrl_view_gen_models.h>
#include <rrl/rrl_convergence_on_weighted_nas.h>
#include <rrl/rrl_view_based_reg_sym.h>

#include <rgrl/rgrl_feature_set_bins.h>
#include <rgrl/rgrl_matcher_k_nearest_pick_one.h>
#include <rgrl/rgrl_weighter_m_est.h>
#include <rgrl/rgrl_scale_est_closest.h>
#include <rgrl/rgrl_scale_est_all_weights.h>
#include <rgrl/rgrl_est_affine.h>
#include <rgrl/rgrl_data_manager.h>
#include <rgrl/rgrl_data_manager_sptr.h>
#include <rgrl/rgrl_scale.h>
#include <rgrl/rgrl_convergence_tester_sptr.h>
#include <rgrl/rgrl_feature_face_pt.h>
#include <rgrl/rgrl_feature_point.h>
#include <rgrl/rgrl_mask.h>

#include <rrel/rrel_m_est_obj.h>
#include <rrel/rrel_muset_obj.h>
#include <rrel/rrel_tukey_obj.h>

fregl_coarse2fine_register::
fregl_coarse2fine_register( rrl_pair_reg_params_sptr reg_params, 
                            InternalImageType::Pointer from_image, 
                            InternalImageType::Pointer to_image)
  :params_(reg_params)
{
  // Step0: Create features
  const bool generate_at_all_scales = false; //don't generate at all scale
  FeatureType::Pointer from_features_3d =  new FeatureType( from_image, generate_at_all_scales );
  from_features_3d->set_max_num_driving_faces(reg_params->max_driving_faces_);
  from_features_3d->generate_at_all_scales();

  FeatureType::Pointer to_features_3d =  new FeatureType( to_image, generate_at_all_scales );
  to_features_3d->set_max_num_driving_faces(reg_params->max_driving_faces_);
  to_features_3d->generate_at_all_scales();
  
  
  // Step1: Concatenate all features
  //
  make_vectors_of_feature_sets( from_features_3d,
                                concat_moving_matchable_corners_,
                                concat_moving_driving_corners_, 
                                concat_moving_matchable_faces_,
                                concat_moving_driving_faces_,
                                params_->feature_min_scale_, 
                                params_->feature_max_scale_, 
                                params_->c2f_upper_min_feature_scale_);
  make_vectors_of_feature_sets( to_features_3d,
                                concat_fixed_matchable_corners_,
                                concat_fixed_driving_corners_, 
                                concat_fixed_matchable_faces_,
                                concat_fixed_driving_faces_,
                                params_->feature_min_scale_, 
                                params_->feature_max_scale_, 
                                params_->c2f_upper_min_feature_scale_ );

  // if the size are not equivalent
  if( concat_moving_driving_faces_.size() != concat_fixed_driving_faces_.size() ) {
    
    const unsigned new_size = vnl_math_min( concat_moving_driving_faces_.size(),
                                            concat_fixed_driving_faces_.size() );
      
    concat_moving_matchable_corners_.resize( new_size);
    concat_moving_driving_corners_.resize( new_size); 
    concat_moving_matchable_faces_.resize( new_size);
    concat_moving_driving_faces_.resize( new_size);
    concat_fixed_matchable_corners_.resize( new_size);
    concat_fixed_driving_corners_.resize( new_size);
    concat_fixed_matchable_faces_.resize( new_size);
    concat_fixed_driving_faces_.resize( new_size);
  }

  //  Feature sets should have already been created.
  assert( concat_fixed_driving_faces_.size() > 0 );
}

bool
fregl_coarse2fine_register::
run( rgrl_transformation_sptr  init_xform )
{ 
  // step 1: Create a matcher
  //
  vcl_vector<rgrl_matcher_sptr> matchers;
  const int num_nearest_to_test = 3;
  const double distance_threshold_for_match = 30;
  matchers.resize( 4 );
  for( unsigned i=0; i<4; ++i ) {
    matchers[i] = new rgrl_matcher_k_nearest_pick_one( num_nearest_to_test, distance_threshold_for_match ); 
  }

  // Step 2: Create a robust weighter
  //
  const double tukey_parameter = 4;
  vcl_auto_ptr< rrel_m_est_obj > est_obj( new rrel_tukey_obj(tukey_parameter) );
  
  //  The parameter "use_precomputed_signature_wgt" is set to be
  //  true because the weight based on scale and (for face points)
  //  orientation will have been computed during the matching
  //  process.
  
  const bool use_signature_error = false;
  const bool use_precomputed_signature_wgt = true;
  rgrl_weighter_sptr wgter = new rgrl_weighter_m_est( est_obj, use_signature_error, use_precomputed_signature_wgt ) ;

  // Step 3:  Create a scale estimator.  
  //
  // muse and unwgted_scale_est are used in the first iteration. Do
  // NOT use sk refinement
  const bool max_lookup_table_size = 0;   //  Prevents costly and unnecessary preliminary generation of a table
  const bool use_sk_refinement_in_muse = false;
  vcl_auto_ptr<rrel_objective> obj( new rrel_muset_obj( max_lookup_table_size,
                                                        use_sk_refinement_in_muse ) );
  const bool estimate_signature_covar = false;
  rgrl_scale_estimator_unwgted_sptr unwgted_scale_est = 
    new rgrl_scale_est_closest( obj, estimate_signature_covar ); 
  unwgted_scale_est->set_debug_flag(2);
  unwgted_scale_est->set_warning( false );  // remove "unstable estimates" warning
  rgrl_scale_estimator_wgted_sptr wgted_scale_est = 
    new rgrl_scale_est_all_weights( estimate_signature_covar );  
  wgted_scale_est->set_debug_flag(2);

  // Step 4: Create the estimator. There is only one and it is affine.
  //
  rgrl_estimator_sptr estimator = new rgrl_est_affine(3);

  // Step 5: Create a simple view generator that includes model
  //  selection but not region growing.  Don't need an initializer,
  //  just an initial view. In fact, we don't even need model
  //  seletion, since there is only one model in our
  //  estimation. However, this seems like the least expensive model
  //  generator. Should switch to feature-based approach later.
  bool use_rho_function = true;  // in forming the error objective function
  rrl_model_selection_criterion_sptr model_sel = new rrl_model_selection_aicc(use_rho_function);
  rrl_view_generator_sptr c2f_view_generator = new rrl_view_gen_models( model_sel );
  
  
  // Step 6: Create a data manager.
  //
  bool use_multi_stage = true;
  rgrl_data_manager_sptr data_manager = new rgrl_data_manager( use_multi_stage );

  // Step 7: Determine the max number of resolutions and set up the
  // data manager. The problem is that two images could have different
  // number of scales available.  Take the minimum!
  const unsigned num_resolutions 
    = vnl_math_min( concat_moving_driving_faces_.size(), concat_fixed_driving_faces_.size() );
    
  for ( unsigned int res=0; res<num_resolutions; ++res ) {
    if( !params_->use_faces_only_ ) {
      data_manager->add_data( res, concat_moving_driving_corners_[res],
                              concat_fixed_matchable_corners_[res],
                              matchers[0], wgter, unwgted_scale_est, wgted_scale_est );
      // backward matching
      if( !params_->only_forward_matching_ )
        data_manager->add_data( res, concat_moving_matchable_corners_[res],
                                concat_fixed_driving_corners_[res],
                                matchers[1], wgter, unwgted_scale_est, wgted_scale_est );
    }                                

    if( !params_->use_corners_only_ ) {
      data_manager->add_data( res, concat_moving_driving_faces_[res],
                              concat_fixed_matchable_faces_[res],
                              matchers[2], wgter, unwgted_scale_est, wgted_scale_est ); 
      // backward matching
      if( !params_->only_forward_matching_ )
        data_manager->add_data( res, concat_moving_matchable_faces_[res],
                                concat_fixed_driving_faces_[res],
                                matchers[3], wgter, unwgted_scale_est, wgted_scale_est ); 
    }
      
    data_manager->add_estimator( res, estimator );
    data_manager->set_dimension_increase_for_next_stage( res, 1.0 );   //  All features are in coords of original res
  }
  vcl_cout << "rrl_pair_reg::coarse_to_fine_registration:  num_resolutions = "
           << num_resolutions << vcl_endl;

  //  Step 8: Set the initial scales based on what has been passed in.
  //  Need a scale for each feature type and direction.
  rgrl_scale_sptr initial_scale = new rgrl_scale();
  double initial_error_scale = 1.5*params_->c2f_upper_min_feature_scale_;
  initial_scale->set_geometric_scale( initial_error_scale );

  //  Step 9. Create a convergence tester.  Use Gary's.  The only thing
  //  there we really don't need is the test on region growing, but as
  //  long as we ensure that the regions indeed do not change we're ok. 
  rgrl_convergence_tester_sptr conv_test ;
  const bool force_check_all_init = false;
  conv_test= new rrl_convergence_on_weighted_nas( force_check_all_init,
                                                  params_->error_threshold_upper_, 
                                                  params_->error_threshold_lower_ );
  conv_test->set_rel_tol( params_->relative_error_change_for_convergence_ );

  //  2.10: Create a registration object and set limits on scale estimation.
  rrl_view_based_reg_sptr c2f_registration_ptr;

  if( !params_->only_forward_matching_ )  {
    c2f_registration_ptr 
       = new rrl_view_based_reg_sym( data_manager, 
                                     c2f_view_generator, 
                                     conv_test );
  } else {
    c2f_registration_ptr 
      = new rrl_view_based_registration( data_manager, 
                                         c2f_view_generator, 
                                         conv_test );
  }
  c2f_registration_ptr->set_debug_flag( 2 );

  c2f_registration_ptr->set_expected_min_geometric_scale(0.5);
  c2f_registration_ptr->set_expected_max_geometric_scale(10.0);
  c2f_registration_ptr->set_max_icp_iter( 8 );

  // 2.11: Create the initial view
  rgrl_mask_box moving_box = concat_moving_matchable_faces_[0]->bounding_box();
  rgrl_mask_sptr moving_box_sptr = new rgrl_mask_box( moving_box );
  rgrl_mask_box fixed_box = concat_fixed_matchable_faces_[0]->bounding_box();
  rgrl_mask_sptr fixed_box_sptr = new rgrl_mask_box( fixed_box );
      
  init_view_ = new rgrl_view( moving_box_sptr, 
                              fixed_box_sptr,
                              moving_box,
                              moving_box,
                              estimator,
                              init_xform,
                              0,
                              init_xform->inverse_transform() );   
  
  // Final step: run the refinement!
  // coarse-to-fine refinement is done with the reg engine. 
  // This for loop ensures that when reg fails at a coarse level,
  // it can still try a finer level where there are more features. 
  for( unsigned res=num_resolutions-1; res>0&&!c2f_registration_ptr->has_best_view(); --res ) {
    // make a copy of the view (the view may be rrl_view_sym type)
    rgrl_view_sptr initial_view = init_view_->self_copy();

    // set resolution to the coarsest one. 
    // No need to scale as everything is in image coordinate. 
    initial_view->set_resolution( res );
    
    vcl_cout << "Coarse-to-fine registration is starting at level=(" << res << ")..." << vcl_endl;
    c2f_registration_ptr->run( initial_view, initial_scale );
    vcl_cout << "Finished coarse-to-fine registration..." << vcl_endl; 
  }      

  // If a best_view is available, an accurate alignment has been
  // produced.
  if ( !c2f_registration_ptr->has_best_view() ) {
    vcl_cout << "Leaving rrl_pair_reg::coarse_to_fine_registration, without convergence."
             << vcl_endl;
      return false;
    }
  else {
    //final_view = reg_sptr->final_view();
    
    //Store the result!

    return true;
  }
  
}


// This code is from rrl_gdbicp_info::rkpl_faces_to_rgrl_features. So,
// this code cannot be distriubted.
//
void
fregl_coarse2fine_register::
rkpl_faces_to_rgrl_features( vcl_vector< rkpl_keypoint<3> > const& rkpl_faces,
                             feature_vector & features )
{
  vnl_vector<double> loc(3), dir(3);

  for( unsigned int i=0; i<rkpl_faces.size(); ++i )
    {
      loc[0] = rkpl_faces[i].physical_location[0];
      loc[1] = rkpl_faces[i].physical_location[1];
      loc[2] = rkpl_faces[i].physical_location[2];
      dir[0] = rkpl_faces[i].direction[0];
      dir[1] = rkpl_faces[i].direction[1];
      dir[2] = rkpl_faces[i].direction[2];
      rgrl_feature_sptr fptr = new rgrl_feature_face_pt( loc, dir );
      fptr->set_scale( rkpl_faces[i].smoothing_scale );
      features.push_back( fptr );
    }
}

// This code is from rrl_gdbicp_info::rkpl_corners_to_rgrl_features. So,
// this code cannot be distriubted.
//
void
fregl_coarse2fine_register::
rkpl_corners_to_rgrl_features( vcl_vector< rkpl_keypoint<3> > const& rkpl_corners,
                               feature_vector & features )
{
  vnl_vector<double> loc(3);

  for( unsigned int i=0; i<rkpl_corners.size(); ++i )
    {
      loc[0] = rkpl_corners[i].physical_location[0];
      loc[1] = rkpl_corners[i].physical_location[1];
      loc[2] = rkpl_corners[i].physical_location[2];
      rgrl_feature_sptr fptr = new rgrl_feature_point( loc );
      fptr->set_scale( rkpl_corners[i].smoothing_scale );
      features.push_back( fptr );
    }
}

// This code is from
// rrl_gdbicp_info::make_vectors_of_feature_sets. So, this code cannot
// be distriubted.
//
void 
fregl_coarse2fine_register::
make_vectors_of_feature_sets(FeatureType::Pointer features,
                             feature_set_vector & concat_matchable_corners,
                             feature_set_vector & concat_driving_corners,
                             feature_set_vector & concat_matchable_faces,
                             feature_set_vector & concat_driving_faces,
                             double min_scale,
                             double max_scale,
                             double upper_min_feature_scale)
{
  feature_vector as_rgrl_features;

  // check on features available
  if( !features->num_scales() ) {
    vcl_cerr<<"Warning: Empty features set."<<vcl_endl;
    return;
  }

  // Here we consider the entire set of features.
  int min_index = 0;
  int max_index = features->num_scales()-1;

  int upper_min_index = min_index;
  while ( upper_min_index < max_index &&
          upper_min_feature_scale > 0 &&
          features->scale(upper_min_index) <= upper_min_feature_scale )
    ++ upper_min_index;
  if( upper_min_index > min_index )
    --upper_min_index;

  const int bin_size_for_matchable_faces = 10;
  const int bin_size_for_driving_faces = 15;
  const int bin_size_for_matchable_corners = 20;
  const int bin_size_for_driving_corners = 25;

  //  Matchable corners
  as_rgrl_features.reserve(5000);
  concat_matchable_corners.resize( upper_min_index - min_index + 1 );
  for ( int i=max_index; i>=min_index; --i ) {
    rkpl_corners_to_rgrl_features( features->matchable_corners(i), 
                                   as_rgrl_features );
    if ( i <= upper_min_index ) {
      rgrl_feature_set_sptr fset 
        = new rgrl_feature_set_bins<3>( as_rgrl_features, 
                                        bin_size_for_matchable_corners,
                                        vcl_string("MATCHABLE:CORNER") );
      concat_matchable_corners[i-min_index] = fset;
      
      vcl_cout<<"Created a matchable corner feature set at index " << i-min_index
              << ", and size " << as_rgrl_features.size() << vcl_endl;

      DebugMacro( 3, "Created a matchable corner feature set at index " << i-min_index
                  << ", and size " << as_rgrl_features.size() << vcl_endl );
      
    }
  }
  
  //  Driving corners
  as_rgrl_features.clear();
  as_rgrl_features.reserve(5000);
  concat_driving_corners.resize( upper_min_index - min_index + 1 );
  for ( int i=max_index; i>=min_index; --i ) {
    rkpl_corners_to_rgrl_features( features->driving_corners(i), 
                                   as_rgrl_features );
    if ( i <= upper_min_index ) {
      rgrl_feature_set_sptr fset 
        = new rgrl_feature_set_bins<3>( as_rgrl_features, 
                                        bin_size_for_driving_corners,
                                        vcl_string("DRIVING:CORNER") );
      concat_driving_corners[i-min_index] = fset;
     
      vcl_cout<<"Created a driving corner feature set at index " << i-min_index
              << ", and size " << as_rgrl_features.size() << vcl_endl;

      DebugMacro(3, "Created a driving corner feature set at index " << i-min_index
                 << ", and size " << as_rgrl_features.size() << vcl_endl );
      
    }
  }

  //  Matchable faces
  as_rgrl_features.clear();
  as_rgrl_features.reserve(5000);
  concat_matchable_faces.resize( upper_min_index - min_index + 1 );
  for ( int i=max_index; i>=min_index; --i ) {
    rkpl_faces_to_rgrl_features( features->matchable_faces(i), 
                                 as_rgrl_features );
    if ( i <= upper_min_index ) {
      rgrl_feature_set_sptr fset
        = new rgrl_feature_set_bins<3>( as_rgrl_features, 
                                        bin_size_for_matchable_faces,
                                        vcl_string("MATCHABLE:FACE") );
      concat_matchable_faces[i-min_index] = fset;
     
      vcl_cout<< "Created a matchable faces feature set at index " 
              << i-min_index
              << ", and size " << as_rgrl_features.size() 
              << vcl_endl ;

      DebugMacro( 3, "Created a matchable faces feature set at index " 
                  << i-min_index
                  << ", and size " << as_rgrl_features.size() 
                  << vcl_endl );
    }
  }

  //  Driving faces
  as_rgrl_features.clear();
  as_rgrl_features.reserve(5000);
  concat_driving_faces.resize( upper_min_index - min_index + 1 );
  for ( int i=max_index; i>=min_index; --i ) {
    rkpl_faces_to_rgrl_features( features->driving_faces(i), 
                                 as_rgrl_features );
    if ( i <= upper_min_index ) {
      rgrl_feature_set_sptr fset 
        = new rgrl_feature_set_bins<3>( as_rgrl_features, 
                                        bin_size_for_driving_faces,
                                        vcl_string("DRIVING:FACE") );
      concat_driving_faces[i-min_index] = fset;
      
      vcl_cout<< "Created a driving faces feature set at index " 
              << i-min_index
              << ", and size " << as_rgrl_features.size() 
              << vcl_endl;

      DebugMacro( 3, "Created a driving faces feature set at index " 
                  << i-min_index
                  << ", and size " << as_rgrl_features.size() 
                  << vcl_endl );
    }
  }
}
