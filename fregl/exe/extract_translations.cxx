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

//: Generate a text file containing all the translations to the top-left corner of the montage given the anchor image.
//
//  Usage:
//
//  extract_translations xml_joint_transforms reference_image outputfile
//
//  where
//   xml_transforms         Name of the xml file containing the transformations
//   reference_image        Name of the reference image name found in the xml_transforms
//   outputfile             Name of the ascii file containg the translations

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
	typedef unsigned short						InputPixelType;
	typedef itk::AffineTransform< double, 3>	TransformType;

	// 1. Get the input arguments
	vul_arg< vcl_string > arg_file_xforms  ( 0, "The xml file containing transformations." );
	vul_arg< vcl_string > arg_file_to    ( 0, "The reference image name in the transformation file." );
	vul_arg< vcl_string > arg_output    ( 0,"Name of the ascii file containg the translations");
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
	std::ofstream output;
	output.open(arg_output().c_str());

	for (unsigned int img_ind = 0; img_ind<image_names.size(); img_ind++) {
		TransformType::Pointer xform = joint_register->get_transform(image_names[img_ind], arg_file_to().c_str());
		if ( xform ) { //there might be no xform for the image pair
			TransformType::ParametersType params = xform->GetParameters();
			output<<image_names[img_ind]<<" "<<params[9]-origin[0]<<" "<<params[10]-origin[1]<<" "<<params[11]-origin[2]<<std::endl;
		}
	}

	output.close();
	return 0;
}

