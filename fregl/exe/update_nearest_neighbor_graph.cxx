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

//: Update the ave_nbr_dist_ by computing the nearest n neighbors of the same class
//  The input is a file containing the xml file of the joint
//  transformations and the xml file containing the result sets. The
//  changes will be written back to the original xml segmentation
//  files.
//
//
//  update_nearest_neighbor_graph xml_joint_transforms xml_result_sets
// 
//  where
//   xml_joint_transforms   Name of the xml file containing the joint transformations
//   xml_result_sets        Name of the xml file containing the result sets
//

#include <vul/vul_arg.h>

#include <fregl/fregl_result_record.h>
#include <fregl/fregl_joint_register.h>
#include <fregl/fregl_space_transformer.h>

#include <maciej_seg/maciejSegmentation.h>
#include <maciej_seg/xml_util.h>

int 
main( int argc, char* argv[] )
{
  vul_arg< vcl_string > arg_xml_xform  ( 0, "Name of the xml file containing the joint transformations" );
  vul_arg< vcl_string > arg_xml_result ( 0, "Name of the xml file containing the result sets" );
  vul_arg< vcl_string > arg_anchor ("-anchor", "Name of the anchor image. If not set, all images are treated as anchor in turn");
  vul_arg< double > arg_k ("-k", "The parameter k for k-Nearest-Neighbor", 6);
  vul_arg< double > arg_x ("-x", "The pixel distance in x",1);
  vul_arg< double > arg_y ("-y", "The pixel distance in y",1);
  vul_arg< double > arg_z ("-z", "The in-between slice distance",1);

  vul_arg_parse( argc, argv );

  // Read in the xml_xform file
  fregl_joint_register::Pointer joint_register = new fregl_joint_register( arg_xml_xform() );
  std::vector<fregl_result_record::Pointer> result_records;
  result_record_read_xml(arg_xml_result(), result_records);

  // Set the space transformer, which allows only the anchor image and
  // the adjacent images to be operated on
  //
  fregl_space_transformer space_transformer(joint_register);
  bool overlap_only = true;
  std::vector<std::string> const & image_names = joint_register->image_names();

  // Take each image as the anchor in turn. Load the xml files of the
  //anchor and the overlapping images. Compute the ave_nbr_dist_ and
  //write the xml files back to disk. The neighbors are of a nucleus
  //have to be in the same class as the nucleus. Duplicated nuclei are
  //not considered in the nearest-neighbor queries
  //
  for (unsigned int img_ind = 0; img_ind<image_names.size(); img_ind++) {
    std::string anchor_name = image_names[img_ind];
    if (arg_anchor.set() && anchor_name!=arg_anchor()) continue;
    
    std::cout<<"anchor image = "<<anchor_name<<std::endl;
    bool in_anchor = true;
    space_transformer.set_anchor( anchor_name, in_anchor, overlap_only );

    // load the xml files of the anchor and the overlapping images
    std::vector<std::string> overlapping_image_names = space_transformer.image_names();
    std::vector<maciejSegmentation::Pointer> nuclear_seg_results(overlapping_image_names.size());
    int anchor_index;
    
    std::vector<std::string> xml_nuclear_names;
    std::vector<std::string> xml_nuclear_image_paths;
    std::vector<std::string> xml_nuclear_image_names;
    std::string  image_path;
    std::string  image_name;
    for (unsigned int o_img_ind = 0; o_img_ind<overlapping_image_names.size(); o_img_ind++) {
      if (overlapping_image_names[o_img_ind] == anchor_name) anchor_index = o_img_ind;
      for (unsigned int k = 0; k<result_records.size(); k++) {  
        if (result_records[k]->registration_image() == overlapping_image_names[o_img_ind]) {
          std::cout<<"read xml for "<<result_records[k]->registration_image()<<std::endl;
          // assuming everything to be in the current directory
          nuclear_seg_results[o_img_ind] = new maciejSegmentation();
          xml_util_read( "./", result_records[k]->nuclear_xml(),
                         image_path ,image_name ,*nuclear_seg_results[o_img_ind] );
          xml_nuclear_names.push_back(result_records[k]->nuclear_xml());
          xml_nuclear_image_paths.push_back( image_path );
          xml_nuclear_image_names.push_back( image_name );
          break;
        }
      }
    }

    // Repeat the same process for each class. First, scan the cell
    // list of the anchor image to find the max class number.
    //
    maciejSegmentation::Pointer anchor_nuclear_seg = nuclear_seg_results[anchor_index];
    std::vector<rich_cell::Pointer> const & cells_in_anchor = anchor_nuclear_seg->all_cells();
        
    int max_class_num = 0; 
    for (unsigned int i = 0; i<cells_in_anchor.size(); i++) {
      if (cells_in_anchor[i]->class_type_ > max_class_num) 
        max_class_num = cells_in_anchor[i]->class_type_;
    }

    for (int ci=1; ci<=max_class_num; ci++) {
      // Build the binning structure for the anchor and adjacent
      // images to enable the k-NN queries.
      double bin_size = 20;
      bool dup_removed = true;
      int class_label = ci;
      for (unsigned int j = 0; j<overlapping_image_names.size(); j++) {
        nuclear_seg_results[j]->construct_bins(bin_size, class_label, dup_removed);
      }

      // Compute the k-NN. For each "non-dup" nucleus in the anchor
      // image, compute the k-NN from the anchor image, and determine
      // the distance to the farest neighbor. Using the radius to
      // query nuclei in the overlapping image, if neighbors are found
      // in the queries, determine the nearest k neighbor distances
      // from all the neighbors returned by the queries.
      //
      for (std::vector<rich_cell::Pointer>::const_iterator ci = cells_in_anchor.begin(); ci != cells_in_anchor.end(); ci++ ) {

        // Skip the computation for invalid, wrong-class and
        // duplicated cells
        if (!(*ci)->valid_ || (*ci)->class_type_ != class_label || (*ci)->dup_)
          continue;

        // Compute the k-NN from the anchor image, and compute the max
        // distance from the k neighbors
        std::vector<rich_cell::Pointer>  nearby_cells_in_anchor;
        nuclear_seg_results[anchor_index]->k_nearest_cells( nearby_cells_in_anchor, (*ci)->center_, arg_k()+1); //+1 because itsself is also a NN
        double R = 0;
        vnl_vector_fixed< float, 3 > diff_vec, scaled_diff_vec; 
        std::vector<double> dist_list;
        std::vector< double > k_nn_dists;
        double dist, pixel_dist;;
        k_nn_dists.reserve( arg_k()+1 );
        for (unsigned int nni = 0; nni<nearby_cells_in_anchor.size(); nni++) {
          diff_vec = (*ci)->center_ - nearby_cells_in_anchor[nni]->center_;
          pixel_dist = diff_vec.magnitude();
          diff_vec[0] *= arg_x();
          diff_vec[1] *= arg_y();
          diff_vec[2] *= arg_z();
          dist = diff_vec.magnitude();
          dist_list.push_back( dist );
          k_nn_dists.push_back( dist );
          if (pixel_dist > R) R = pixel_dist;
        }

        // Compute the neighbors in the adjacent image within the
        // radius, R
        std::vector<rich_cell::Pointer> nearby_cells_neighbors;
        nearby_cells_neighbors.reserve(10);
        for (unsigned int o_img_ind = 0; o_img_ind<overlapping_image_names.size(); o_img_ind++) {
          if (o_img_ind == anchor_index) continue; 

          nearby_cells_neighbors.clear();
          vnl_vector_fixed< float, 3 > xformed_loc;
          if (!space_transformer.in_image((*ci)->center_, o_img_ind, xformed_loc)) continue;

          nuclear_seg_results[o_img_ind]->cells_within_radius(nearby_cells_neighbors, xformed_loc, R);
          for (unsigned int nni = 0; nni<nearby_cells_neighbors.size(); nni++) {
            diff_vec = xformed_loc - nearby_cells_neighbors[nni]->center_;
            diff_vec[0] *= arg_x();
            diff_vec[1] *= arg_y();
            diff_vec[2] *= arg_z();
            dist = diff_vec.magnitude();
            dist_list.push_back( dist );
          }
        }
        
        // Now find the smallest k+1 distances from the
        // nearby_cells_in_anchor and nearby_cells_neighbors. the "+1"
        // is to take care of the nucleus itself, since it is also
        // returned as one of the k-NN in the anchor image.
        double max_dist;
        int max_dist_ind;
        bool max_replaced = true;
        for (unsigned int nni = 0; nni < dist_list.size(); nni++) {
          //replace the max
          if (max_replaced) {
            max_dist = 0;
            for (unsigned int k_nni = 0; k_nni<k_nn_dists.size(); k_nni++) {
              if (max_dist > k_nn_dists[k_nni]) {
                max_dist =  k_nn_dists[k_nni];
                max_dist_ind = k_nni;
              }
            }
          }
          
          if (dist_list[nni] < max_dist) {
            max_replaced = true;
            k_nn_dists[ max_dist_ind ] = dist_list[nni];
          }
          else max_replaced = false;
        }

        // Compute the averate distance and update the ave_nbr_dist_
        double ave = 0;
        for (unsigned int k_nni = 0; k_nni<k_nn_dists.size(); k_nni++) {
          ave += k_nn_dists[ k_nni ];
        }
        (*ci)->ave_nnbr_dist_ = ave/k_nn_dists.size();

      } //for all cells in the anchor
    } //for all classes

    // Output the xml file of the anchor image
    std::cout<<"Update the xml file "<<xml_nuclear_names[anchor_index]<<" ..."<<std::endl;
    xml_util_write("./", xml_nuclear_names[anchor_index], 
                   xml_nuclear_image_paths[anchor_index],
                   xml_nuclear_image_names[anchor_index], 
                   *anchor_nuclear_seg );
    /*
    for (unsigned int j = 0; j < nuclear_seg_results.size(); j++) {
      xml_util_write("./", xml_nuclear_names[j], xml_nuclear_image_paths[j],
                     xml_nuclear_image_names[j], *(nuclear_seg_results[j]) );
    }
    */
   } //for each image as the anchor

   return 0;
}
