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

//: Executable program to mosaic a set of images with given transformations
//
//  The input is an xml file containing either one pairwise
//  registration result or a set of transformations from joint
//  registration. The images can be gray, rgb color or rgba color. The
//  input images are assumed TIF images. The output montage is a gray
//  image. The usage is
//  
//  mosaic_images xml_file 
//
//  where
//    xml_file      Name of the xml_file containing the xforms
//    anchor_image  Name of the anchor image
//
//  Optional arguments"
//    -path         The path of the image files.
//    -old_str      The old substr in the image names to be replaced
//    -new_str      The replacement of the old substr
//    -output       The output image name.

#ifdef _OPENMP
#include "omp.h"
#endif

#include <iostream>
using std::cerr;
using std::endl;

#include <vul/vul_arg.h>
#include <vul/vul_file.h>

#include "Templates/mosaic_images_template.h"


int
main(int argc, char* argv[]) {
    vul_arg< vcl_string > arg_xml_file(0, "A xml file containing transformations");
    vul_arg< vcl_string > arg_anchor(0, "Anchor image name");
    vul_arg< int > arg_channel("-channel", "The color channel (0-red, 1-green, 2-blue), or the image channel if the original image is a lsm image.", 0);
    vul_arg< vcl_string > arg_img_path("-path", "The path of the image files.", ".");
    vul_arg< vcl_string > arg_old_str("-old_str", "The old substr in the image names to be replaced");
    vul_arg< vcl_string > arg_new_str("-new_str", "The new substr in the image names");
    vul_arg< bool > arg_3d("-3d", "Generate a 3D image as well", false);
    vul_arg< vcl_string > arg_outfile("-output", "The name of the output directory for the stack slices.");
    vul_arg< bool > arg_in_anchor("-in_anchor", "The final space is set to the anchor image", false);
    vul_arg< bool > arg_overlap("-overlap_only", "Only consider images that overlap the anchor image", false);
    vul_arg< bool > arg_nn("-nn", "Use Nearest-Neighbor interpolation", false);
	vul_arg< bool > arg_normalize("-normalize","Normalize the intensity of each tile",false);
	vul_arg< vcl_string > arg_background("-background", "Background image name", "");
	vul_arg< double > arg_sigma("-sigma","The Gaussian blur param for normalization", 50);
	vul_arg< double > arg_median("-median","The scale median for normalization",1000);
    vul_arg< int > arg_blending("-blending", "0: max (default), 1: even weighted, 2: photopleaching weighted (the fanciest).", 0);
    vul_arg< bool > arg_denoise("-denoise", "Making an attempt to remove noise of high frequencies", false);
	vul_arg< bool > arg_write_proj2d("-debug", "Write 2d projection of the intermediate mosaic", false);

    vul_arg_parse(argc, argv);

    int retcode = mosaic_images_template<unsigned char>(arg_xml_file, arg_anchor, arg_channel, arg_img_path, arg_old_str, arg_new_str, arg_3d, arg_outfile, arg_in_anchor, arg_overlap, arg_nn, arg_normalize, arg_background, arg_sigma, arg_median, arg_blending, arg_denoise, arg_write_proj2d);
    return retcode;
}
