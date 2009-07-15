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

//: Compute the distance to the top of the anchor image.
//
//  The input is the file containing the xml file of the joint
//  transformations and the xml file containing the result sets
//  involved in the query with the anchor being the first image.

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
  vul_arg< vcl_string > arg_xml_result ( 0, "Name of the xml file containing the result records involved in the query." );
  vul_arg< int > arg_class (0, "The class", 0);
  vul_arg< vcl_string > arg_outfile("-output","Name of the output file","query_distance.txt"); 
  vul_arg< double > arg_x ("-x", "The pixel distance in x",1);
  vul_arg< double > arg_y ("-y", "The pixel distance in y",1);
  vul_arg< double > arg_z ("-z", "The in-between slice distance",1);
  vul_arg_parse( argc, argv ); 

  // Read in the xml_xform file
  fregl_joint_register::Pointer joint_register = new fregl_joint_register( arg_xml_xform() );
  /*
  std::vector<maciejSegmentation::Pointer> nuclear_seg_results;
  std::ifstream infile(arg_xml_result());
  std::string line, xml_file;
  std::string  image_path;
  std::string  image_name;
  while (infile) {
     std::getline(infile, line);
     if (line.length() == 0) continue;

     std::istringstream line_stream(line);
     line_stream>>xml_file;
     maciejSegmentation::Pointer nuclear_seg_result = new maciejSegmentation();
     xml_util_read( "./", xml_file, image_path ,image_name ,
                    *nuclear_seg_result );
     nuclear_seg_results.push_back(nuclear_seg_result);  
  }
  */

  std::vector<fregl_result_record::Pointer> result_records;
  result_record_read_xml(arg_xml_result(), result_records);

  // Set the space transformer, which allows only the anchor image and
  // the adjacent images to be operated on
  //
  fregl_space_transformer space_transformer(joint_register);
  std::string anchor_name = result_records[0]->registration_image();
  space_transformer.set_anchor( anchor_name );
  std::cout<<"anchor image = "<<anchor_name<<std::endl;
  std::vector<std::string> image_names = space_transformer.image_names();
  std::ofstream outfile(arg_outfile().c_str());

  // For each result record listed, transform the point to the anchor
  // space and record the x distance.
  for (unsigned int i = 0; i<result_records.size(); i++) {
    maciejSegmentation::Pointer nuclear_seg_result = new maciejSegmentation();

    // find the transform index of the from_image
    int from_index = -1;
    for (unsigned int img_ind = 0; img_ind<image_names.size(); img_ind++) {
      if ( image_names[img_ind] == result_records[i]->registration_image() ) {
        from_index = img_ind;
        std::string image_path, image_name;
        xml_util_read( "./", result_records[i]->nuclear_xml(), 
                       image_path ,image_name ,*nuclear_seg_result );
        break;
      }
      assert (from_index > -1 );
    }

    //transform all cells from result_records[i] to anchor space
    std::vector<rich_cell::Pointer> const & cells = nuclear_seg_result->all_cells();
    for (std::vector<rich_cell::Pointer>::const_iterator ci = cells.begin(); ci != cells.end(); ci++ ) {
      if (!(*ci)->valid_ || (*ci)->class_type_ != arg_class() || (*ci)->dup_) 
        continue;

      vnl_vector_fixed< float, 3 > xformed_loc;
      space_transformer.in_anchor((*ci)->center_, from_index, xformed_loc);
      outfile<<xformed_loc[0]*arg_x()<<"\t"<<xformed_loc[1]*arg_y()<<"\t"<<xformed_loc[2]*arg_z()<<"\t"<<(*ci)->ave_nnbr_dist_<<"\t"<<(*ci)->volume_<<std::endl;
    }
  }
  outfile.close();
  return 0;
}
