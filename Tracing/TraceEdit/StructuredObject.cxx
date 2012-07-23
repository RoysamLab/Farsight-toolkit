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
void StructuredObject::sphereKernel(ImageType::Pointer input, ImageType::Pointer &mask, int radius)
{
	ImageType::RegionType region = input->GetLargestPossibleRegion();

	mask->SetRegions(region);
	mask->Allocate();

	EllipseType::Pointer sphere = EllipseType::New();
	SpatialObjectToImageFilterType::Pointer imageFilter = SpatialObjectToImageFilterType::New();

	EllipseType::ArrayType radiusArray;
	radiusArray[0] = radius;
	radiusArray[1] = radius;
	radiusArray[2] = radius;
	sphere->SetRadius( radiusArray );

	//Position the sphere
	EllipseTransformType::Pointer ellipseTransform = EllipseTransformType::New();
	EllipseTransformType::OutputVectorType translation;
	translation[0] = radius;
	translation[1] = radius;
	translation[2] = radius;

	ellipseTransform->SetIdentity();
	ellipseTransform->Translate( translation, false );
	sphere->SetObjectToParentTransform( ellipseTransform );
	
	ImageType::SizeType outputSize;
	outputSize[0] = radius*2+1;
	outputSize[1] = radius*2+1;
	outputSize[2] = radius*2+1;
	imageFilter->SetSize( outputSize );

	//std::cout << "Mask size: " << outputSize[0] << " " <<outputSize[1] << " " << outputSize[2] <<std::endl;

	imageFilter->SetInput(sphere);
	sphere->SetDefaultInsideValue(255);
	sphere->SetDefaultOutsideValue(0);
	imageFilter->SetUseObjectValue(true);
	imageFilter->SetOutsideValue(0);
	imageFilter->Update();

	//WriterType::Pointer writer = WriterType::New();
	//writer->SetFileName("sphere.tif");
	//writer->SetInput(imageFilter->GetOutput());

	//std::cout<< "Sphere radius: " << radius << std::endl;
	//
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
void StructuredObject::taperedCylinderMask(ImageType::Pointer image, ImageType::Pointer &mask, int firstBit[], double radius1, int secondBit[], double radius2)
{//not yet functional
	//check if length/distanceXYZ is greater than 1 otherwise use circle kernel

	//try creating multiple circle images and combining them using AddImageFilter

	ImageType::RegionType region = image->GetLargestPossibleRegion();

	mask->SetRegions(region);
	mask->Allocate();

	ImageType::SizeType regionSize = region.GetSize();

	//add to this image
	EllipseType::Pointer blank = EllipseType::New();
	SpatialObjectToImageFilterType::Pointer blankFilter = SpatialObjectToImageFilterType::New();
	blankFilter->SetSize( regionSize );
	//blankFilter->SetInput(blank);
	//circle->SetDefaultInsideValue(0);
	//circle->SetDefaultOutsideValue(0);
	blankFilter->SetUseObjectValue(true);
	blankFilter->SetOutsideValue(0);
	blankFilter->Update();


	//double cur_radius = radius1;
	//double length = sqrt(pow(distanceXYZ[0],2)+pow(distanceXYZ[1],2)+pow(distanceXYZ[2],2));
	//std::cout << "Length: " << length << std::endl;
	//int start_z = regionSize[2]/2-length/2;
	//int end_z = regionSize[2]/2+length/2;
	//double plane_size = (regionSize[0]+1)*(regionSize[1]+1);

	////number of voxels along the length
	//int num_of_voxels = distanceXYZ[0];
	//for (int i = 1; i < 2; i++)
	//{
	//	if (distanceXYZ[i] > num_of_voxels)
	//		num_of_voxels = distanceXYZ[i];
	//}

	//if (num_of_voxels == 0)
	//{
	//	num_of_voxels = 1;
	//}

	//double radius_increment = (radius2-radius1)/num_of_voxels;

	//double centerPt[3];
	//centerPt[0] = firstBit[0];
	//centerPt[1] = firstBit[1];
	//centerPt[2] = firstBit[2];

	///////////
	////make a circle for each voxel along the length
	//for (int i = 0; i < num_of_voxels; i++)
	//{
	//	int radius = (int) floor(cur_radius+0.5);
	//	SpatialObjectToImageFilterType::Pointer imageFilter = SpatialObjectToImageFilterType::New();
	//	
	//	//ImageType::SizeType outputSize;
	//	//outputSize[0] = radius*2+1;
	//	//outputSize[1] = radius*2+1;
	//	//outputSize[2] = radius*2+1;
	//	imageFilter->SetSize( regionSize );

	//	EllipseType::Pointer circle = EllipseType::New();

	//	//define radius and position of the circle, and size of image
	//	EllipseType::ArrayType radiusArray;
	//	radiusArray[0] = radius;
	//	radiusArray[1] = radius;
	//	radiusArray[2] = 1.0;
	//	circle->SetRadius(radiusArray);

	//	EllipseType::TransformType::OffsetType offset;
	//	offset[0] = centerPt[0];
	//	offset[1] = centerPt[1];
	//	offset[2] = centerPt[2];
	//	circle->GetObjectToParentTransform()->SetOffset(offset);

	//	EllipseType::TransformType::CenterType center;
	//	center[0] = centerPt[0];
	//	center[1] = centerPt[1];
	//	center[2] = centerPt[2];
	//	circle->GetObjectToParentTransform()->SetCenter( center );

	//	//rotate
	//	EllipseType::TransformType::OutputVectorType axis;
	//	axis[0] = 0;
	//	axis[1] = 1;
	//	axis[2] = 0;
	//	circle->GetObjectToParentTransform()->Rotate3D(axis,elevation,false);
	//	axis[0] = 0;
	//	axis[1] = 0;
	//	axis[2] = 1;
	//	circle->GetObjectToParentTransform()->Rotate3D(axis,azimuth,false);
	//	circle->ComputeObjectToWorldTransform();

	//	imageFilter->SetInput(circle);
	//	circle->SetDefaultInsideValue(255);
	//	circle->SetDefaultOutsideValue(0);
	//	imageFilter->SetUseObjectValue(true);
	//	imageFilter->SetOutsideValue(0);
	//	imageFilter->Update();

	//	//ImageType::Pointer slantCircle = imageFilter->GetOutput();

	//	WriterType::Pointer writer = WriterType::New();
	//	writer->SetFileName("circle.tif");
	//	writer->SetInput( imageFilter->GetOutput() );

	//	std::cout<< "Circle radius: " << radius << std::endl;
	//	
	//	//try {
	//	//	imageFilter->Update();
	//	//	writer->Update();
	//	//}
	//	//catch (itk::ExceptionObject & excp )
	// //   {
	//	//	std::cerr << excp << std::endl;
	// //   }


	//}
	//	mask = imageFilter->GetOutput();
	//	mask->SetOrigin( input->GetOrigin() );
	//	cur_radius += radius_increment;
	//////


	//itk::ImageRegionIterator< ImageType > imageIterator(mask,region);
	//imageIterator.GoToBegin();
	////draw tapered cylinder along z axis
	//while(!imageIterator.IsAtEnd())
	//{
	//	double distanceX = fabs(imageIterator.GetIndex()[0] - referencePt[0]);
	//	double distanceY = fabs(imageIterator.GetIndex()[1] - referencePt[1]);
	//	double distanceR = sqrt( pow(distanceX,2)+pow(distanceY,2) );

	//	if (imageIterator.GetIndex()[2] < start_z || imageIterator.GetIndex()[2] > end_z)
	//	{
	//		imageIterator.Set(0);
	//	}
	//	else
	//	{
	//		if(distanceR > cur_radius)
	//		{
	//			imageIterator.Set(0);
	//		}
	//		else
	//		{
	//			imageIterator.Set(255);
	//		}

	//		++count;
	//		if (count == plane_size) //next radius
	//		{
	//			count = 0;
	//			cur_radius += radius_increment;
	//			referencePt[2]++;
	//		}
	//	}
	//	++imageIterator;
	//}
	//
	////rotate tapered cylinder
	////Euler3D
	//TransformType::Pointer transform = TransformType::New();
	//
	//TransformType::InputPointType centerOfRotation;
	//centerOfRotation[0] = referencePt[0];
	//centerOfRotation[1] = referencePt[1];
	//centerOfRotation[2] = regionSize[2]/2;
	//transform->SetCenter(centerOfRotation);
	//
	//double angleX = atan2(distanceXYZ[2],distanceXYZ[1]);
	//double angleY = atan2(distanceXYZ[2],distanceXYZ[0]);
	//double angleZ = atan2(distanceXYZ[1],distanceXYZ[0]);
	////transform->SetRotation(angleX,angleY,angleZ);
	//transform->SetRotation(0,0,0);

	//transform->Print( std::cout , 3 );

	//ResampleFilterType::Pointer resampler = ResampleFilterType::New();
	//resampler->SetInput( mask );
	//resampler->SetReferenceImage( mask );
	//resampler->UseReferenceImageOn( );
	//resampler->SetTransform( transform );

	//For validating accuracy
	//WriterType::Pointer writer1 = WriterType::New();
	//writer1->SetFileName( "Input.tif" );
	//writer1->SetInput( image );

	WriterType::Pointer writer2 = WriterType::New();
	writer2->SetFileName( "Mask_Output.tif" );
	writer2->SetInput( blankFilter->GetOutput() );
	
	try
	{
		//resampler->Update();
		//writer1->Update();
		writer2->Update();
		//mask = resampler->GetOutput();
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cerr << excp << std::endl;
	}

	////FixedCenterAffineTransform
	//TransformType::Pointer transform = TransformType::New();

	//TransformType::InputPointType centerOfRotation;
	//centerOfRotation[0] = referencePt[0];
	//centerOfRotation[1] = referencePt[1];
	//centerOfRotation[2] = regionSize[2]/2;
	//transform->SetCenterOfRotationComponent(centerOfRotation);

	//double vector[3];
	//for (int i = 0; i<3; i++)
	//{
	//	vector[i] = distanceXYZ[i]/length;
	//}
	//TransformType::MatrixType rotationMatrix;
	//double angleX = atan2(distanceXYZ[2],distanceXYZ[1]);
	//double angleY = atan2(distanceXYZ[2],distanceXYZ[0]);
	//double angleZ = atan2(distanceXYZ[1],distanceXYZ[0]);
	//rotationMatrix[0][0]=cos(angleX)+pow(vector[0],2)*(1-cos(angleX));
	//rotationMatrix[0][1]=vector[0]*vector[1]*(1-cos(angleX))-vector[2]*sin(angleX);
	//rotationMatrix[0][2]=vector[0]*vector[2]*(1-cos(angleX))+vector[1]*sin(angleX);
	//rotationMatrix[1][0]=vector[1]*vector[0]*(1-cos(angleY))+vector[2]*sin(angleY);
	//rotationMatrix[1][1]=cos(angleY)+pow(vector[1],2)*(1-cos(angleY));
	//rotationMatrix[1][2]=vector[1]*vector[2]*(1-cos(angleY))-vector[0]*sin(angleY);
	//rotationMatrix[2][0]=vector[2]*vector[0]*(1-cos(angleZ))-vector[1]*sin(angleZ);
	//rotationMatrix[2][1]=vector[2]*vector[1]*(1-cos(angleZ))+vector[0]*sin(angleZ);
	//rotationMatrix[2][2]=cos(angleZ)+pow(vector[2],2)*(1-cos(angleZ));

	//rotationMatrix[0]=cos(angleY)*cos(angleX);
	//rotationMatrix[1]=sin(angleZ)*sin(angleY)*cos(angleX)-cos(angleZ)*sin(angleX);
	//rotationMatrix[2]=cos(angleZ)*sin(angleY)*cos(angleX)+sin(angleZ)*sin(angleX);
	//rotationMatrix[3]=cos(angleY)*sin(angleX);
	//rotationMatrix[4]=sin(angleZ)*sin(angleY)*sin(angleX)+cos(angleZ)*cos(angleX);
	//rotationMatrix[5]=cos(angleZ)*sin(angleY)*sin(angleX)-sin(angleZ)*cos(angleX);
	//rotationMatrix[6]= -1*sin(angleY);
	//rotationMatrix[7]=sin(angleZ)*cos(angleY);
	//rotationMatrix[8]=cos(angleZ)*cos(angleY);
	//transform->SetMatrixComponent(rotationMatrix);

	//mask->SetTransform( transform );

	//ImageType::ArrayType radiusArray;
	//radiusArray[0] = radius;
	//radiusArray[1] = radius;

	//ImageType::TransformType::OffsetType offset;
	//offset[0] = radius;
	//offset[1] = radius;
	//offset[2] = radius;

	//QuickView viewer;
	//viewer.AddImage(image);
	//viewer.AddImage(mask);
	//viewer.Visualize();
}