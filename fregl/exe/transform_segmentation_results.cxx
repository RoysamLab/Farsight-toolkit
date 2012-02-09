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
//  The input includes the xml file of the transformations, the name
//  of the reference space, and the list of (from_image_name
//  feature_file_name) pairs. The pair for the reference image should
//  be included if its features are part of the global space too. Each
//  feature_file starts with the name of the image for segmentation as
//  the first line. Following the image name are the features for
//  segmented objects, one line per object. The format is as follows
//  and only the first 4 fields are essential:
//
//  image_name_for_segmentation
//  ID x y z feature1 feature2 ...
//  ID x y z feature1 feature2 ...
//  ...
//
//  The output file contains the following information:
//
//  image_name_for_segmentation FROM from_image_name TO reference_image_name
//  ID xformed_x xformed_y xformed_z feature1 feature2 ...
//  ID xformed_x xformed_y xformed_z feature1 feature2 ...
//  ...
//
//  Only the coordinates are transformed. The features are ignored and written
//  directly to the transformed output file.
//
//  Usage:
//
//  tranform_segmentation_results xml_joint_transforms reference_image
//  feature_file_list 
//
//  where
//   xml_transforms         Name of the xml file containing the transformations
//   reference_image        Name of the reference image name found in the xml_transforms
//   feature_file_list      List of (from_image_name feature_file_name) pairs.
//
//  Optional arguments:
//
//  -output_prefix          The prefix for the output file. If not given, the default is "transformed_"  
//

#include <fregl/fregl_joint_register.h>
#include <fregl/fregl_space_transformer.h>

#include <vnl/vnl_vector_fixed.h>
#include <vul/vul_file.h>
#include <vul/vul_arg.h>

#include <vector>
#include <string>
#include <fstream>
int 
main( int argc, char* argv[] )
{
	typedef unsigned short InputPixelType;
	// 1. Get the input arguments
	vul_arg< vcl_string > arg_file_xforms  ( 0, "The xml file containing transformations." );
	vul_arg< vcl_string > arg_file_to    ( 0, "The reference image name in the transformation file." );
	vul_arg< vcl_string > arg_feature_list    ( 0, "List of (from_image_name feature_file_name) pairs." );
	vul_arg< vcl_string > arg_prefix    ( "-output_prefix", " The prefix for the output file. If not given, the default is transformed_", "transformed_");

	vul_arg_parse( argc, argv );

	// 2. Get the space transformer ready
	//
	fregl_joint_register< InputPixelType >::Pointer joint_register = new fregl_joint_register< InputPixelType >( arg_file_xforms() );
	fregl_space_transformer< InputPixelType > space_transformer(joint_register);
	bool overlap_only = false;
	bool in_anchor = false;
	space_transformer.set_anchor( arg_file_to(), in_anchor, overlap_only );
	fregl_space_transformer< InputPixelType >::PointType origin = space_transformer.origin();
	std::vector<std::string> const& image_names = space_transformer.image_names();
	// 3. Read in the pairs in the feature_file_list and process each
	// feature file
	std::ifstream inList;
	inList.open(arg_feature_list().c_str());
	if ( !inList ){
		std::cerr<<"Couldn't open "<<arg_feature_list()<<std::endl;
		exit( 0 );
	}

	std::string line_str;
	std::string from_image_name, feature_file_name;
	while ( inList ) {
		std::getline( inList, line_str );
		if (line_str.length() == 0) continue;

		std::istringstream line_stream(line_str);
		line_stream>> from_image_name >> feature_file_name;

		// Locate the transformation
		int from_img_ind = 0;
		for (unsigned int img_ind = 0; img_ind<image_names.size(); img_ind++) {
			if (from_image_name == image_names[img_ind]) {
				from_img_ind = img_ind;
				break;
			}
		}

		// Process each feature file
		//
		// Prepare the input file
		std::ifstream inFeatures;
		inFeatures.open(feature_file_name.c_str());
		if ( !inFeatures ){
			std::cerr<<"Couldn't open "<<feature_file_name<<std::endl;
			exit( 0 );
		}
		std::string result_file_name;
		inFeatures >> result_file_name;
		std::cout<<"Transforming "<<feature_file_name<<std::endl;

		//Prepare the output file
		std::string output_file = arg_prefix()+feature_file_name;
		std::ofstream outFeatures;
		outFeatures.open(output_file.c_str());
		outFeatures<<result_file_name<<" FROM "<<from_image_name<<" TO "<<arg_file_to()<<std::endl;

		// Parce each line and transform the location to the global space.
		int id;
		vnl_vector_fixed< float, 3 > point, xformed_pt;
		std::string remaining_str = "hello";

		while (inFeatures) {
			std::getline( inFeatures, line_str );
			if (line_str.length() == 0) continue;

			// Get the first 4 essential components and the remaining string
			std::istringstream line_stream1(line_str);
			line_stream1>> id >> point[0] >> point[1] >> point[2];
			remaining_str = line_str.substr((int)line_stream1.tellg()+1);

			// Transform the point, shift the point by the origin, and write
			// to the output file output
			space_transformer.in_anchor( point, from_img_ind, xformed_pt);
			outFeatures<<id<<" "<<xformed_pt[0]-origin[0]<<" "<<xformed_pt[1]-origin[1]<<" "<<xformed_pt[2]-origin[2]<<" "<<remaining_str<<std::endl;

		}
		inFeatures.close();
	}
	inList.close();

	return 0;
}

