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
#ifndef IMAGEACTORS_H_
#define IMAGEACTORS_H_

#include "vtkSmartPointer.h"
#include "vtkContourFilter.h"
#include "vtkActor.h"
#include "vtkImageData.h"
#include "vtkImageToStructuredPoints.h"
#include "vtkLODActor.h"
#include "vtkOpenGLVolumeTextureMapper3D.h"
#include "vtkPolyDataMapper.h"
#include "vtkVolume.h"
#include "vtkVolumeProperty.h"
#include "vtkProperty.h"

#include "itkImageFileReader.h"
#include "itkImageToVTKImageFilter.h"

#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <QtGui>



const unsigned int Dimension = 3;
typedef unsigned char  PixelType;
typedef itk::Image< PixelType, Dimension >   ImageType;
typedef itk::ImageFileReader< ImageType >    ReaderType;
typedef itk::ImageToVTKImageFilter<ImageType> ConnectorType;
struct imageFileHandle
{
	std::string tag;
	std::string filename;
	vtkSmartPointer<vtkImageData> ImageData;

};
class  ImageRenderActors
{
public:
	ImageRenderActors();
	//ImageRenderActors(std::string ImageSource);
	int loadImage(std::string ImageSource);
	vtkSmartPointer<vtkActor> ContourActor(vtkImageData *image);
	//vtkSmartPointer<vtkImageData> ImageReader(std::string ImageSource);
private:
	std::vector<imageFileHandle*> LoadedImages;

};
#endif