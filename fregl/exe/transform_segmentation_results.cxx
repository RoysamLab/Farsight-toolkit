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

//: Generate text files containing transformed points in the global space
//
//  The input includes the xml file of the joint transformations and
//  the xml file containing the result sets. The output will be
//  written to a directory with the name
//  "transformed_to_anchor_name". Only surace points and trace points
//  are transformed. For the nuclear results, each class has its own
//  output file.
//
//  Usage:
//
//  tranform_segmentation_results xml_joint_transforms xml_result_sets
//  anchor_image_name
//
//  where
//   xml_joint_transforms   Name of the xml file containing the joint transformations
//   xml_result_sets        Name of the xml file containing the result sets
//   anchor_image_name      Name of the anchor image shown in xml_joint_transforms
//

#include <fregl/fregl_joint_register.h>
#include <fregl/fregl_space_transformer.h>
#include <fregl/fregl_result_record.h>
#include <maciej_seg/maciejSegmentation.h>
#include <maciej_seg/xml_util.h>

#include <vnl/vnl_vector_fixed.h>

#include <vector>
#include <string>

int 
main( int argc, char* argv[] )
{
  if (argc < 4) {
    std::cerr<<"Usage: "<<argv[0]<<" xml_joint_transforms xml_result_sets anchor_image_name"<<std::endl;
    return 0;
  }
  std::string anchor = argv[3];
  std::string path = std::string("transformed_to_")+anchor;
  std::string cmd1 = "mkdir "+path;
  std::system(cmd1.c_str());
  path = path + std::string("/");
  
  // Get the space transformer ready
  //
  fregl_joint_register::Pointer joint_register = new fregl_joint_register( argv[1] );
fregl_space_transformer space_transformer(joint_register);
  bool overlap_only = false;
  bool in_anchor = false;
  space_transformer.set_anchor( anchor, in_anchor, overlap_only );
  fregl_space_transformer::PointType origin = space_transformer.origin();

  // Transform the results
  //
  std::vector<fregl_result_record::Pointer> result_records;
  result_record_read_xml(argv[2], result_records);
  std::vector<std::string> const& image_names = space_transformer.image_names();

  std::string  image_path;
  std::string  image_name;
  std::string output_file;
  std::vector< vnl_vector_fixed< float, 3 > > points;
  vnl_vector_fixed< float, 3 > point, xformed_pt, one(1,1,1);

  for (unsigned int img_ind = 0; img_ind<image_names.size(); img_ind++) {

    int record_index;
    for (unsigned int k = 0; k<result_records.size(); k++) {
      if (result_records[k]->registration_image() == image_names[img_ind]) {
        record_index = k;
        break;
      }
    }

    // Transform the nuclear boundary points *************************
    //
    if (result_records[record_index]->nuclear_xml().length() > 0) {
      points.clear();
      std::cout<<"Transforming "<<result_records[record_index]->nuclear_xml()<<" ..."<<std::endl;
      maciejSegmentation::Pointer nuclear_result = new maciejSegmentation();
      xml_util_read( "./", result_records[record_index]->nuclear_xml(),
                     image_path ,image_name ,*nuclear_result );
      
      // count the number of classes
      std::vector<rich_cell::Pointer> const& cells = nuclear_result->all_cells();
      int class_count = 0;
      for (unsigned i = 0; i<cells.size(); i++) {
        if (class_count < cells[i]->class_type_) class_count = cells[i]->class_type_;
      }
      // For each class generate one text file
      for (int ci = 1; ci<= class_count; ci++) {
        char class_num[5] ;
        std::sprintf( class_num, "%d", ci);
        output_file = path+result_records[record_index]->nuclear_xml()+std::string("_class")+ class_num+std::string(".txt");
        std::ofstream outfile;
        outfile.open(output_file.c_str());
        
        for (int i = 0; i<cells.size(); i++) {
          if (cells[i]->class_type_ != ci || !cells[i]->valid_ || cells[i]->dup_) continue;
          
          const std::vector< vnl_vector_fixed< float, 3 > >& points = cells[i]->boundary_points_;
          for (int pi = 0; pi<points.size(); pi++) {
            // transform the point and output in [z,y,x] format
            space_transformer.in_anchor( points[pi], img_ind, xformed_pt);
            outfile << xformed_pt[2]-origin[2]<<"\t"<<xformed_pt[1]-origin[1]<<"\t"<<xformed_pt[0]-origin[0]<<"\t0\t0\n";
          }
          
        }
        outfile.close();
      }
    }

    // Transform the vessel surface points ************************
    //
    points.clear();
    std::string vessel_result = result_records[record_index]->vessel_xml();
    if ( vessel_result.length() > 0 ) {
      std::ifstream infile;
      std::ofstream outfile;
      std::cout<<"Transforming "<<vessel_result<<" ..."<<std::endl;
      infile.open( vessel_result.c_str(), std::ifstream::in );
      if ( !infile ){
        std::cerr<<"Couldn't open "<<vessel_result<<std::endl;
        exit( 0 );
      }
      output_file = path+std::string("transformed_")+vessel_result;
      outfile.open(output_file.c_str());
      
      std::string line_str;
      while ( infile ) {
        std::getline( infile, line_str );
        if (line_str.length() == 0) continue;
        
        std::istringstream line_stream(line_str);
        line_stream>> point[2] >> point[1] >> point[0];
        //We have to subtract [1,1,1] from the position, since Arun starts
        //the index from [1,1,1], instead of [0,0,0]
        point = point-one;
        space_transformer.in_anchor( point, img_ind, xformed_pt);
        outfile << xformed_pt[2]-origin[2]<<"\t"<<xformed_pt[1]-origin[1]<<"\t"<<xformed_pt[0]-origin[0]<<"\t0\t0\n";
      }
      infile.close();
      outfile.close();
    }

    // Transform the astrocyte points ******************************
    //
    points.clear();
    std::string astrocyte_result = result_records[record_index]->astrocyte_xml();
    if ( astrocyte_result.length() > 0 ) {
      std::ifstream infile;
      std::ofstream outfile;
      std::cout<<"Transforming "<<astrocyte_result<<" ..."<<std::endl;
      infile.open( astrocyte_result.c_str(), std::ifstream::in );
      if ( !infile ){
        std::cerr<<"Couldn't open "<<astrocyte_result<<std::endl;
        exit( 0 );
      }
      output_file = path+std::string("transformed_")+astrocyte_result;
      outfile.open(output_file.c_str());
      
      std::string line_str;
      while ( infile ) {
        std::getline( infile, line_str );
        if (line_str.length() == 0) continue;
        
        std::istringstream line_stream(line_str);
        line_stream>> point[0] >> point[1] >> point[2];
        if (point[0] < 0) { //delineator
          outfile << point[0]<<"\t"<<point[1]<<"\t"<<point[2]<<std::endl;
        }
        else {
          space_transformer.in_anchor( point, img_ind, xformed_pt);
          outfile << xformed_pt[0]-origin[0]<<"\t"<<xformed_pt[1]-origin[1]<<"\t"<<xformed_pt[2]-origin[2]<<"\t0\t0\n";}
      }
      infile.close();
      outfile.close();
    }

    // Transform the microglia trace points ************************
    //
    points.clear();
    std::string microglia_result = result_records[record_index]->microglia_xml();
    if ( microglia_result.length() > 0 ) {
      std::ifstream infile;
      std::ofstream outfile;
      std::cout<<"Transforming "<<microglia_result<<" ..."<<std::endl;
      infile.open( microglia_result.c_str(), std::ifstream::in );
      if ( !infile ){
        std::cerr<<"Couldn't open "<<microglia_result<<std::endl;
        exit( 0 );
      }
      output_file = path+std::string("transformed_")+microglia_result;
      outfile.open(output_file.c_str());
      
      std::string line_str;
      while ( infile ) {
        std::getline( infile, line_str );
        if (line_str.length() == 0) continue;
        
        std::istringstream line_stream(line_str);
        line_stream>> point[0] >> point[1] >> point[2];
        if (point[0] < 0) { //delineator
          outfile << point[0]<<"\t"<<point[1]<<"\t"<<point[2]<<std::endl;
        }
        else {
          space_transformer.in_anchor( point, img_ind, xformed_pt);
          outfile << xformed_pt[0]-origin[0]<<"\t"<<xformed_pt[1]-origin[1]<<"\t"<<xformed_pt[2]-origin[2]<<"\t0\t0\n";}
      }
      infile.close();
      outfile.close();
    }
  }

  return 0;
}

