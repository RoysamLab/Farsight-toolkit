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


#ifndef _fregl_pairwise_register_h_
#define _fregl_pairwise_register_h_

#include <vcl_string.h>

#include "itkImage.h"
#include "itkCommand.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkAffineTransform.h"
#include <rgrl/rgrl_transformation_sptr.h>
#include "fregl/fregl_util.h"

// For debugging purpose
class CommandIteration : public itk::Command {
public:
	typedef CommandIteration  Self;
	typedef itk::Command         SuperClass;
	typedef itk::SmartPointer< Self >    Pointer;
	itkNewMacro( Self );
protected:
	CommandIteration() {};

public:
	typedef  itk::RegularStepGradientDescentOptimizer   OptimizerType;
	typedef  const OptimizerType                      * OptimizerPointer;

	void Execute(  itk::Object * caller, const itk::EventObject & event ) 
	{
		this-> Execute(  (const itk::Object *) caller,  event  );
	}

	void Execute(  const itk::Object * caller, const itk::EventObject & event )
	{
		OptimizerPointer optimizer = dynamic_cast< OptimizerPointer >( caller );
		if( typeid( event  ) == typeid( itk::IterationEvent ) )  
		{
			std::cout << optimizer->GetCurrentIteration() << " : ";
			std::cout << optimizer->GetValue() << " : ";
			std::cout << optimizer->GetCurrentPosition() << std::endl;
		}   
	}
};

template < class TPixel >
class fregl_pairwise_register
{
public:

	typedef TPixel								InputPixelType;
	typedef itk::Image< InputPixelType, 3 >		InputImageType;
	typedef typename InputImageType::Pointer	InputImageTypePointer;
	
	typedef float								InternalPixelType;
	typedef itk::Image<InternalPixelType, 3 >	InternalImageType;
	
	typedef itk::Image<InputPixelType, 2 >		ImageType2D;
	typedef typename ImageType2D::Pointer		ImageType2DPointer;

	typedef itk::AffineTransform< double, 3>   TransformType;

	//: constructor
	fregl_pairwise_register(InputImageTypePointer from_image, 
		InputImageTypePointer to_image, 
		std::string from_image_filename,
		std::string to_image_filename,
		float background = 30);

	//: Destructor
	~fregl_pairwise_register();

	//: Trigger the registration
	bool run(double& obj_value, const vcl_string & gdbicp_exe_path = "", bool scaling = false);

	//: Trigger the registration with an initial 3D transformation
	bool run(TransformType::Pointer prior_xform, double& obj_value);

	//: Trigger the registration with an x-y displacement
	bool run(double init_x, double init_y, double& obj_value);

	//: Set the from_image
	void set_from_image(InputImageTypePointer from_image);

	//: Set the to_image
	void set_to_image(InputImageTypePointer to_image);

	//: Set the flag for exhaustive search in z direction
	void set_exhausive_search(bool exhaustive);

	//: Set the stack size of the to_image if the stack is too large to run
	//
	//  The stack is extracted from the center portion of the image
	void set_stack_size(int size);

	//: Unset the stack size, so the entire image stack is registered
	void unset_stack_size();

	//: Get the transform
	TransformType::Pointer transform();

	//: Set the variance for smoothing
	void set_smoothing(double variance);

private: 

	//: Compute the z shift 
	//
	//  This function also set the region of interest to the overlap
	//  area.
	double compute_z_shift( rgrl_transformation_sptr fw_xform, 
		InputImageTypePointer from_image_3d,
		InputImageTypePointer to_image_3d,
		ImageType2DPointer from_image_2d, 
		ImageType2DPointer to_image_2d,
		float bg );

	//: Crop the image to the region of interest
	//
	//  This step is to reduce memory consumption, since it is too
	//  expensive to have image type of float.
	InternalImageType::Pointer crop_image(InputImageTypePointer image);

	//: Check the validity of the 2D xform
	//
	//  The affine components of the 2D xform should be more or less
	//  identity.
	bool valid_2d_xform( rgrl_transformation_sptr xform2d, bool scaling );

	//: Convert the transform from the rgrl format to itk format
	void convert_rgrl_to_itk_xform( rgrl_transformation_sptr affine2d,
		double z_displacement,
		TransformType::ParametersType& parameters);

	//: Read the xform from the transformation file dumped by gdbicp
	rgrl_transformation_sptr read_2d_xform( vcl_istream& reg_info );

	//: A utility function to remove the white spaces at the beginning and end of a string.
	vcl_string string_trim (const vcl_string& source, const vcl_string& pat);

private:
	std::string from_image_filename;
	std::string to_image_filename;
	InputImageTypePointer from_image_; 
	InputImageTypePointer to_image_;
	TransformType::Pointer  transform_;
	float background_;    //threshold value to be considered as background
	bool exhaustive_;
	bool stack_size_set_;
	int stack_size_;
	double smoothing_;
};

#endif
