#ifndef REGISTER_PAIR_TEMPLATE_H
#define REGISTER_PAIR_TEMPLATE_H

#include <iostream>

#include "itkImageFileReader.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkRGBAPixel.h"
#include "itkRGBPixel.h"
#include "itkImage.h"
#include "itkTIFFImageIO.h"
#include "itkImageFileWriter.h"

#include <vul/vul_file.h>
#include <vul/vul_arg.h>
#include <vul/vul_timer.h>

#include <fregl/fregl_pairwise_register.h>
#include <fregl/fregl_joint_register.h>
#include <fregl/fregl_reg_record.h>
#include <fregl/fregl_util.h>

template <typename TPixel>
int
	register_pair_template(
	vul_arg< vcl_string >	arg_file_from,
	vul_arg< vcl_string >	arg_file_to,
	vul_arg< int >			channel,
	vul_arg< float >		background,
	vul_arg< double >       smooth,
	vul_arg< vcl_string >	gdbicp,
	vul_arg< int >          slices,
	vul_arg<bool>           remove_2d,
	vul_arg<bool>           scaling_arg,
	vul_arg< vcl_string >   prior_arg
);

template <typename TImageType3D, typename TInternalImageType>
typename TImageType3D::Pointer
	smooth_image(typename TImageType3D::Pointer image, int num_sub_images );

#endif