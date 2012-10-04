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

//: Executable program to register a set of 3D image with prior pairwise transforms
//
//  The input is a file containing the xml files of the prior pairwise
//  results. The output is an xml file containing one record for every
//  image pair.
//
//   register_pair input_file 
//
//
// where 
//  input_file      Name of the file containing the xml files of the prior pairwise results 
//
// Optional arguments:
//  -output         Name of the output xml file. 
//  -error          Upper bound of the error to be considered as correct pairwise alignment

#include <vul/vul_arg.h>
#include "Templates/register_joint_template.h"

int
	main(  int argc, char* argv[] )
{
	vul_arg< vcl_string > arg_in_file     ( 0, "A file containing filenames of xml files, each containing a pairwise transformation." );
	vul_arg< vcl_string > arg_xml_file    ( "-output", "Output xml filename","joint_transforms.xml" );
	vul_arg< double > arg_multiplier    ( "-multiplier", "The multiplier for the error scale. 4 is a good value.",0);
	vul_arg< double > arg_error_bound    ( "-error_bound", "The upper bound for the accepted error in the range of [0,1]. The default is 1 (all pairs are accepted)",1);
	vul_arg< vcl_string > arg_roi_file ("-roi", "Text file containing the list of image names in the ROI");
	vul_arg< bool > arg_debug ("-debug","Dump out the statistics to a temp file which has the same name as the xml file with the suffix debug.txt",false);

	vul_arg_parse( argc, argv );

	int retcode = register_joint_template<unsigned short>(arg_in_file, arg_xml_file, arg_multiplier, arg_error_bound, arg_roi_file, arg_debug);

	return retcode;
}