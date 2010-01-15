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
#include "ImageActors.h"

ImageRenderActors::ImageRenderActors()
{
	this->LoadedImages.clear();
}
int ImageRenderActors::loadImage(std::string ImageSource)
{
	if (ImageSource.empty())
	{
		return -1;
	}
	imageFileHandle *newImage;
	newImage->filename = ImageSource;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( ImageSource );
	//Test opening and reading the input file
	try
	{
		reader->Update();
	}
	catch( itk::ExceptionObject & exp )
	{
		std::cerr << "Exception thrown while reading the input file " << std::endl;
		std::cerr << exp << std::endl;
		//return EXIT_FAILURE;
	}
	ConnectorType::Pointer connector= ConnectorType::New();
	connector->SetInput( reader->GetOutput() );
	newImage->ImageData = connector->GetOutput();
	this->LoadedImages.push_back(newImage);
	return true;
}
vtkSmartPointer<vtkActor> ImageRenderActors::ContourActor(vtkImageData *image)
{
	vtkSmartPointer<vtkContourFilter> ContourFilter = vtkSmartPointer<vtkContourFilter>::New();
	ContourFilter->SetInput(image);
	ContourFilter->SetValue(0,10);
	ContourFilter->Update();
	vtkSmartPointer<vtkPolyDataMapper> ContourMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	ContourMapper->SetInput(ContourFilter->GetOutput());
	vtkSmartPointer<vtkActor>ContourActor = vtkSmartPointer<vtkActor>::New();
	ContourActor->SetMapper(ContourMapper);
	ContourActor->GetProperty()->SetOpacity(.5);
	ContourActor->GetProperty()->SetColor(0.5,0.5,0.5);
	ContourActor->SetPickable(0);
	return ContourActor;
}