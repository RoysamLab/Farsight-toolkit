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
#ifndef STRUCTUREDOBJECT_H
#define STRUCTUREDOBJECT_H

#define PI 3.14159265

#include "itkAddImageFilter.h"
#include "itkEuler3DTransform.h"
#include "itkEllipseSpatialObject.h"
#include "itkImageRegionConstIterator.h"
#include "itkMaskImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkSpatialObjectToImageFilter.h"

#include "ImageActors.h"
#include "itkImageFileWriter.h"

typedef itk::Euler3DTransform< double > TransformType;
typedef itk::Image< ImageActorPixelType, Dimension >   ImageType;
typedef itk::EllipseSpatialObject< Dimension > EllipseType;
typedef EllipseType::TransformType EllipseTransformType;
typedef itk::SpatialObjectToImageFilter< EllipseType, ImageType >	SpatialObjectToImageFilterType;
typedef itk::ResampleImageFilter< ImageType, ImageType > ResampleFilterType;
typedef itk::ImageRegionConstIterator< ImageType > ConstIteratorType;
typedef itk::ImageFileWriter< ImageType > WriterType;
typedef itk::AddImageFilter< ImageType, ImageType > AddImageFilterType;

class StructuredObject
{
public:
	StructuredObject();
	~StructuredObject();
		
	void circleKernel(ImageType::Pointer input, ImageType::Pointer &mask, int centerVoxel[], double radius, double azimuth, double elevation);
	
private:

};
#endif