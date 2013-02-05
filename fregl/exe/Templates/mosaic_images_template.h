#ifndef MOSAIC_IMAGES_TEMPLATE_H
#define MOSAIC_IMAGES_TEMPLATE_H

#include <iostream>

#include <vul/vul_arg.h>
#include <vul/vul_file.h>

#include <fregl/fregl_joint_register.h>
#include <fregl/fregl_space_transformer.h>
#include <fregl/fregl_util.h>
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileWriter.h"
#include "itkImageSliceIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkImageSeriesWriter.h"
#include "itkNumericSeriesFileNames.h"
#include "itkMultiThreader.h"

template <typename TPixel>
int
mosaic_images_template(
	vul_arg< vcl_string > arg_xml_file,
	vul_arg< vcl_string > arg_anchor,
	vul_arg< int > arg_channel,
	vul_arg< vcl_string > arg_img_path,
	vul_arg< vcl_string > arg_old_str,
	vul_arg< vcl_string > arg_new_str,
	vul_arg< bool > arg_3d,
	vul_arg< vcl_string > arg_outfile,
	vul_arg< bool > arg_in_anchor,
	vul_arg< bool > arg_overlap,
	vul_arg< bool > arg_nn,
	vul_arg< bool > arg_normalize, 
	vul_arg< vcl_string > arg_background,
	vul_arg< double > arg_sigma, 
	vul_arg< double > arg_median,
	vul_arg< int > arg_blending,
	vul_arg< bool > arg_denoise,
	vul_arg< bool > arg_write_proj2d
	);

#endif