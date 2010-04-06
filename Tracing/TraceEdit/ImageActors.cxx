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
	//define the points of %90 color
	this->r = 90.0;
	this->g = 45.0;
	this->b = 20.0;
	this->brightness = 150;
	//opacity values
	this->opacity1 = 50;
	this->opacity2 = 100;
	this->opacityTransferFunction = vtkSmartPointer<vtkPiecewiseFunction>::New();
	this->syncOpacityTransfetFunction();
	this->colorTransferFunction = vtkSmartPointer<vtkColorTransferFunction>::New();
	this->syncColorTransfetFunction();
}
int ImageRenderActors::loadImage(std::string ImageSource, std::string tag)
{
	if (ImageSource.empty())
	{
		return -1;
	}
	imageFileHandle *newImage= new imageFileHandle;
	newImage->filename = ImageSource;
	newImage->tag = tag;
	//newImage->colorTransferFunction = 0;
	newImage->ContourActor = 0;
	newImage->ContourFilter = 0;
	newImage->ContourMapper = 0;
	newImage->ImageData = 0;
	//newImage->opacityTransferFunction = 0;
	newImage->volume = 0;
	newImage->volumeMapper = 0;
	newImage->volumeProperty = 0;
	newImage->reader = ReaderType::New();
	newImage->reader->SetFileName( ImageSource );
	//Test opening and reading the input file
	try
	{
		newImage->reader->Update();
	}
	catch( itk::ExceptionObject & exp )
	{
		std::cerr << "Exception thrown while reading the input file " << std::endl;
		std::cerr << exp << std::endl;
		//return EXIT_FAILURE;
	}
	newImage->connector= ConnectorType::New();
	newImage->connector->SetInput( newImage->reader->GetOutput() );
	newImage->ImageData = newImage->connector->GetOutput();
	this->LoadedImages.push_back(newImage);
	return (int) (this->LoadedImages.size() -1);
}
int ImageRenderActors::loadImage(std::string ImageSource, std::string tag, double x, double y, double z)
{
	if (ImageSource.empty())
	{
		return -1;
	}
	imageFileHandle *newImage= new imageFileHandle;
	newImage->filename = ImageSource;
	newImage->tag = tag;
	//newImage->colorTransferFunction = 0;
	newImage->ContourActor = 0;
	newImage->ContourFilter = 0;
	newImage->ContourMapper = 0;
	newImage->ImageData = 0;
	//newImage->opacityTransferFunction = 0;
	newImage->volume = 0;
	newImage->volumeMapper = 0;
	newImage->volumeProperty = 0;
	newImage->reader = ReaderType::New();
	newImage->reader->SetFileName( ImageSource );
	newImage->x = x;
	newImage->y = y;
	newImage->z = z;
	//Test opening and reading the input file
	try
	{
		newImage->reader->Update();
	}
	catch( itk::ExceptionObject & exp )
	{
		std::cerr << "Exception thrown while reading the input file " << std::endl;
		std::cerr << exp << std::endl;
		//return EXIT_FAILURE;
	}
	newImage->connector= ConnectorType::New();
	newImage->connector->SetInput( newImage->reader->GetOutput() );
	newImage->ImageData = newImage->connector->GetOutput();
	this->LoadedImages.push_back(newImage);
	return (int) (this->LoadedImages.size() -1);
}
void ImageRenderActors::ShiftImage(int i, double x, double y, double z)
{
	this->LoadedImages[i]->x = x;
	this->LoadedImages[i]->y = y;
	this->LoadedImages[i]->z = z;
}
vtkSmartPointer<vtkActor> ImageRenderActors::ContourActor(int i)
{
	if (i == -1)
	{
		i = int (this->LoadedImages.size() - 1);
	}
	this->LoadedImages[i]->ContourFilter = vtkSmartPointer<vtkContourFilter>::New();
	this->LoadedImages[i]->ContourFilter->SetInput(this->LoadedImages[i]->ImageData);
	this->LoadedImages[i]->ContourFilter->SetValue(0,10);
	this->LoadedImages[i]->ContourFilter->Update();
	this->LoadedImages[i]->ContourMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	this->LoadedImages[i]->ContourMapper->SetInput(this->LoadedImages[i]->ContourFilter->GetOutput());
	this->LoadedImages[i]->ContourActor = vtkSmartPointer<vtkActor>::New();
	this->LoadedImages[i]->ContourActor->SetMapper(this->LoadedImages[i]->ContourMapper);
	this->LoadedImages[i]->ContourActor->GetProperty()->SetOpacity(.5);
	this->LoadedImages[i]->ContourActor->GetProperty()->SetColor(0.5,0.5,0.5);
	this->LoadedImages[i]->ContourActor->SetPosition(this->LoadedImages[i]->x, 
		this->LoadedImages[i]->y,this->LoadedImages[i]->z);
	this->LoadedImages[i]->ContourActor->SetPickable(0);
	return this->LoadedImages[i]->ContourActor;
}
vtkSmartPointer<vtkVolume> ImageRenderActors::RayCastVolume(int i)
{
	if (i == -1)
	{
		i = int (this->LoadedImages.size() - 1);
	}
	/*this->LoadedImages[i]->opacityTransferFunction = vtkSmartPointer<vtkPiecewiseFunction>::New();
	this->LoadedImages[i]->opacityTransferFunction->AddPoint(2,0.0);
	this->LoadedImages[i]->opacityTransferFunction->AddPoint(50,0.1);
	this->LoadedImages[i]->opacityTransferFunction->AddPoint(100,0.5);
	this->LoadedImages[i]->colorTransferFunction = vtkSmartPointer<vtkColorTransferFunction>::New();
	this->LoadedImages[i]->colorTransferFunction->AddRGBPoint(0.0, 0.0, 0.0, 0.0);
	this->LoadedImages[i]->colorTransferFunction->AddRGBPoint(40.0,0,0,.9);
	this->LoadedImages[i]->colorTransferFunction->AddRGBPoint(90.0,0,.9,0);
	this->LoadedImages[i]->colorTransferFunction->AddRGBPoint(150.0,.9,0,0);*/
	this->LoadedImages[i]->volumeProperty = vtkSmartPointer<vtkVolumeProperty>::New();
	this->LoadedImages[i]->volumeProperty->SetColor(this->colorTransferFunction);
	this->LoadedImages[i]->volumeProperty->SetScalarOpacity(this->opacityTransferFunction);
	this->LoadedImages[i]->volumeProperty->SetInterpolationTypeToLinear();
	this->LoadedImages[i]->volumeMapper = vtkSmartPointer<vtkOpenGLVolumeTextureMapper3D>::New();
	this->LoadedImages[i]->volumeMapper->SetSampleDistance(0.5);
	this->LoadedImages[i]->volumeMapper->SetInput(this->LoadedImages[i]->ImageData);
	this->LoadedImages[i]->volume = vtkSmartPointer<vtkVolume>::New();
	this->LoadedImages[i]->volume->SetMapper(this->LoadedImages[i]->volumeMapper);
	this->LoadedImages[i]->volume->SetProperty(this->LoadedImages[i]->volumeProperty);
	this->LoadedImages[i]->volume->SetPosition(this->LoadedImages[i]->x, 
		this->LoadedImages[i]->y,this->LoadedImages[i]->z);
	this->LoadedImages[i]->volume->SetPickable(0);
	return this->LoadedImages[i]->volume;
}
std::vector<std::string> ImageRenderActors::GetImageList()
{
	this->ImageList.clear();
	for (unsigned int i = 0; i< this->LoadedImages.size(); i++)
	{
		this->ImageList.push_back( this->LoadedImages[i]->filename);
	}
	return this->ImageList;
}
bool ImageRenderActors::isRayCast(int i)
{
	if (i == -1)
	{
		i = int (this->LoadedImages.size() - 1);
	}
	if (this->LoadedImages[i]->tag.compare("Image")==0)
	{
		return true;
	}
	else
	{
		return false;
	}
}
double ImageRenderActors::pointData(int i, int x, int y, int z)
{
	if (i == -1)
	{
		i = int (this->LoadedImages.size() - 1);
	}
	return this->LoadedImages[i]->ImageData->GetScalarComponentAsDouble(x,y,z,0);
}
std::vector<double> ImageRenderActors::GetImageSize(int i)
{
	if (i == -1)
	{
		i = int (this->LoadedImages.size() - 1);
	}
	std::vector<double> imgSize;
	imgSize.push_back((double) this->LoadedImages[i]->volume->GetMaxXBound());
	imgSize.push_back((double) this->LoadedImages[i]->volume->GetMaxYBound());
	imgSize.push_back((double) this->LoadedImages[i]->volume->GetMaxZBound());
	return imgSize;
}
vtkSmartPointer<vtkImageData> ImageRenderActors::GetImageData(int i)
{
	if (i == -1)
	{
		i = int (this->LoadedImages.size() - 1);
	}
	return this->LoadedImages[i]->ImageData;
}
std::vector<double> ImageRenderActors::getColorValues()
{
	std::vector<double> rgb;
	rgb.push_back(this->r);
	rgb.push_back(this->g);
	rgb.push_back(this->b);
	return rgb;
}
void ImageRenderActors::setColorValues(double r, double g, double b)
{
	this->r = r;
	this->g = g;
	this->b = b;
	this->syncColorTransfetFunction();
}
void ImageRenderActors::setColorValues(int i, double value)
{
	if (i==1)
	{
		this->r = value;
	}
	else if(i ==2)
	{
		this->g = value;
	}
	else
	{
		this->b = value;
	}
	this->syncColorTransfetFunction();
}
void ImageRenderActors::setBrightness(int value)
{
	this->brightness = (double)value;
	this->syncColorTransfetFunction();
}
int ImageRenderActors::getBrightness()
{
	return (int) this->brightness;
}
void ImageRenderActors::setOpacity(int value)
{
	this->opacity1 = (double ) value;
	this->syncOpacityTransfetFunction();
}
int ImageRenderActors::getOpacity()
{
	return (int) this->opacity1;
}
void ImageRenderActors::syncColorTransfetFunction()
{
	this->colorTransferFunction->RemoveAllPoints();
	this->colorTransferFunction->AddRGBPoint(0.0, 0.0, 0.0, 0.0);
	this->colorTransferFunction->AddRGBPoint((this->b*this->brightness)/100, 0, 0, .5);//blue
	this->colorTransferFunction->AddRGBPoint((this->g*this->brightness)/100, 0, .5, 0);//green
	this->colorTransferFunction->AddRGBPoint((this->r*this->brightness)/100, .5, 0, 0);//red
	for (unsigned int i = 0; i< this->LoadedImages.size(); i++)
	{
		this->LoadedImages[i]->volume->Update();
	}
}
void ImageRenderActors::syncOpacityTransfetFunction()
{
	this->opacityTransferFunction->RemoveAllPoints();
	this->opacityTransferFunction->AddPoint(2,0.0);
	this->opacityTransferFunction->AddPoint(this->opacity1,0.3);
	//this->opacityTransferFunction->AddPoint(this->opacity2,0.5);
	for (unsigned int i = 0; i< this->LoadedImages.size(); i++)
	{
		this->LoadedImages[i]->volume->Update();
	}
}
