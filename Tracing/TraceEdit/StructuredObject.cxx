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
#include "StructuredObject.h"

StructuredObject::StructuredObject()
{}
StructuredObject::~StructuredObject()
{}
void StructuredObject::circleKernel(ImageType::Pointer input, ImageType::Pointer &mask, int centerVoxel[], double radius, double azimuth, double elevation)//need to validate
{
	/*!
	 * @author Audrey Cheong
	 * @param input input itk image
	 * @param mask circle mask (foreground is white, background is black)
	 * @param centerVoxel specify the center of rotation
	 * @param radius radius of the circle
	 * @param azimuth angle of rotation
	 * @param elevation angle of rotation
	 */
	ImageType::RegionType region = input->GetLargestPossibleRegion();

	mask->SetRegions(region);
	mask->Allocate();

	SpatialObjectToImageFilterType::Pointer imageFilter = SpatialObjectToImageFilterType::New();
	
	ImageType::SizeType outputSize;
	outputSize[0] = region.GetSize()[0];
	outputSize[1] = region.GetSize()[1];
	outputSize[2] = region.GetSize()[2];
	imageFilter->SetSize( outputSize );

	EllipseType::Pointer circle = EllipseType::New();

	//define radius and position of the circle, and size of image
	EllipseType::ArrayType radiusArray;
	radiusArray[0] = radius;
	radiusArray[1] = radius;
	radiusArray[2] = 1.0;
	circle->SetRadius(radiusArray);

	EllipseType::TransformType::OffsetType offset; //position
	offset[0] = centerVoxel[0];
	offset[1] = centerVoxel[1];
	offset[2] = centerVoxel[2];
	circle->GetObjectToParentTransform()->SetOffset(offset);

	EllipseType::TransformType::CenterType center; //center of rotation
	center[0] = centerVoxel[0];
	center[1] = centerVoxel[1];
	center[2] = centerVoxel[2];
	circle->GetObjectToParentTransform()->SetCenter( center );

	//rotate
	EllipseType::TransformType::OutputVectorType axis;
	axis[0] = 0;
	axis[1] = 1;
	axis[2] = 0;
	circle->GetObjectToParentTransform()->Rotate3D(axis,elevation,false);
	axis[0] = 0;
	axis[1] = 0;
	axis[2] = 1;
	circle->GetObjectToParentTransform()->Rotate3D(axis,azimuth,false);
	circle->ComputeObjectToWorldTransform(); //set position and rotation

	imageFilter->SetInput(circle);
	circle->SetDefaultInsideValue(255);
	circle->SetDefaultOutsideValue(0);
	imageFilter->SetUseObjectValue(true);
	imageFilter->SetOutsideValue(0);
	imageFilter->Update();

	//WriterType::Pointer writer = WriterType::New();
	//writer->SetFileName("circle.tif");
	//writer->SetInput( imageFilter->GetOutput() );

	//std::cout<< "Circle radius: " << radius << std::endl;
	
	//try {
	//	imageFilter->Update();
	//	writer->Update();
	//}
	//catch (itk::ExceptionObject & excp )
 //   {
	//	std::cerr << excp << std::endl;
 //   }

	mask = imageFilter->GetOutput();
	mask->SetOrigin( input->GetOrigin() );
}