#ifndef REGISTER_JOINT_TEMPLATE_H
#define REGISTER_JOINT_TEMPLATE_H

#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <vul/vul_arg.h>

#include "itkAffineTransform.h"
typedef itk::AffineTransform< double, 3>   TransformType; //temporary

#include <fregl/fregl_joint_register.h>

template <typename TPixel>
int
register_joint_template(
	vul_arg< vcl_string > arg_in_file,
	vul_arg< vcl_string > arg_xml_file,
	vul_arg< double > arg_multiplier,
	vul_arg< double > arg_error_bound,
	vul_arg< vcl_string > arg_roi_file,
	vul_arg< bool > arg_debug
	);

#endif