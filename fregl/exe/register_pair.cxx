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

//: Executable program to register a pair of 3D images 
// 
//  The outcome is an xml file of a registration record. The images
//  can be gray, rgb color or rgba color. The input images are
//  assumed TIF images. The useage is
//
//   register_pair from_image_name to_image_name
//
//
// where 
//  from_image_name      name of the from_image
//  to_image_name        name of the to_image
// 
// optional arguments:
//  -channel        The channel to be extracted.
//  -background     Intensities below -background are ignored.
//  -smooth         Gaussian smoothing is required
//  -prior          xml file containing the prior transformation for initialization 

#include "Templates/register_pair_template.h"
#include <vul/vul_arg.h>

int
main(  int argc, char* argv[] )
{
	vul_arg< vcl_string > arg_file_from  ( 0, "From image file" );
	vul_arg< vcl_string > arg_file_to    ( 0, "To image file" );
	vul_arg< int >       channel         ("-channel", "The color channel (0-red, 1-green, 2-blue), or the image channel if the original image is a lsm image.",0);
	vul_arg< float >     background      ("-bg", "threshold value for the background", 30);
	vul_arg< double >      smooth        ("-smooth", "If smoothing is performed", 0.5);
	//vul_arg< int >       streaming       ("-streaming","Apply streaming with the number of sub-images when smooth is set", 1);
	//vul_arg< bool >     exhaustive      ("-exhaustive", "Running exhaustive search", false);
	vul_arg< vcl_string > gdbicp ("-gdbicp","The path where the gdbicp exe can be found. The exe of gdbicp can be found at http://www.vision.cs.rpi.edu/gdbicp","");
	vul_arg< int >         slices         ("-slices","Number of slices to register. Needed when registering large image stacks on PC.",100);
	vul_arg<bool>          remove_2d      ("-remove_2d", "Remove the intermediate 2d stuff", false);
	vul_arg<bool>          scaling_arg    ("-scaling","Substantial scaling is expected.", false);
	vul_arg< vcl_string >  prior_arg    ("-prior","xml file containing the initial transformation");

	vul_arg_parse( argc, argv );

	int retcode = register_pair_template<unsigned char>(arg_file_from, arg_file_to, channel, background, smooth, gdbicp, slices, remove_2d, scaling_arg, prior_arg);

	return retcode;
}
