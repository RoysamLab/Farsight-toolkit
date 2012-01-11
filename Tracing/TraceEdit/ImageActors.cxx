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
	this->opacity1Value = .1;
	this->opacity2 = 255;
	this->opacity2Value = 1;
	this->colorValue = 0;
	this->somaColorValue = 0;
	this->sliceBrightness = 500;
	this->RaycastSampleDist = .2;
	this->opacityTransferFunction = vtkSmartPointer<vtkPiecewiseFunction>::New();
	this->syncOpacityTransferFunction();
	this->colorTransferFunction = vtkSmartPointer<vtkColorTransferFunction>::New();
	this->syncColorTransferFunction();
	this->TotalImageSize.clear();
	for (int i = 0; i <6; i++)
	{
		this->TotalImageSize.push_back(0);
	}
}

ImageRenderActors::~ImageRenderActors()
{
  while(!this->LoadedImages.empty())
    {
    delete this->LoadedImages.back();
    this->LoadedImages.pop_back();
    }
}

int ImageRenderActors::loadImage(std::string ImageSource, std::string tag)
{
	return this->loadImage(ImageSource,tag, 0,0,0);
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
	newImage->renderStatus = false;
	newImage->ren2d = false;
	newImage->sliceCreated = false;
	//newImage->sliceActor = 0;
	newImage->imageResliceMapper = 0;
	newImage->imageProperty = 0;
	newImage->imageSlice = 0;
	//newImage->colorTransferFunction = 0;  
	newImage->ContourActor = 0;
	newImage->ContourFilter = 0;
	newImage->ContourMapper = 0;
	newImage->ImageData = 0;
	//newImage->opacityTransferFunction = 0;
	newImage->volume = 0;
	newImage->volumeMapper = 0;
	#ifdef USE_GPUREN
	{
		newImage->volumeMapperGPU = 0;
	}
#endif
	newImage->volumeProperty = 0;
	newImage->reader = ReaderType::New();
	newImage->reader->SetNumberOfThreads(16);
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
	double bounds[6];
	newImage->Rescale = IntensityRescaleType::New();
	newImage->Rescale->SetInput( newImage->reader->GetOutput() );
	newImage->connector = ConnectorType::New();
	newImage->connector->SetInput( newImage->Rescale->GetOutput() );
	newImage->ImageData = newImage->connector->GetOutput();
	newImage->ImageData->Update();
	newImage->ImageData->GetBounds(bounds);
	newImage->projectionConnector = ConnectorType::New();
	newImage->projectionConnector->SetInput( newImage->reader->GetOutput() );

	this->setImageBounds(bounds);
	this->LoadedImages.push_back(newImage);
	return (int) (this->LoadedImages.size() -1);
}
void ImageRenderActors::ShiftImage(int i, double x, double y, double z)
{
	this->LoadedImages[i]->x = x;
	this->LoadedImages[i]->y = y;
	this->LoadedImages[i]->z = z;
}
void ImageRenderActors::ShiftImage(int i, std::vector<double> shift)
{
	this->LoadedImages[i]->x = shift[0];
	this->LoadedImages[i]->y = shift[1];
	this->LoadedImages[i]->z = shift[2];
}
std::vector<double> ImageRenderActors::GetShiftImage(int i)
{
	std::vector<double> shift;
	shift.push_back(this->LoadedImages[i]->x);
	shift.push_back(this->LoadedImages[i]->y);
	shift.push_back(this->LoadedImages[i]->z);
	return shift;
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

vtkSmartPointer<vtkActor> ImageRenderActors::GetContourActor(int i)
{
	if (i == -1)
	{
		i = int (this->LoadedImages.size() - 1);
	}
	return this->LoadedImages[i]->ContourActor;
}

#ifdef USE_GPUREN
void ImageRenderActors::RaycastVolumeMapperGPU(int i)
{
	double max_memory = (5.0 * 1024 * 1024 * 1024) / LoadedImages.size();
	this->LoadedImages[i]->volumeMapperGPU = vtkSmartPointer<vtkOpenGLGPUVolumeRayCastMapper>::New();
	this->LoadedImages[i]->volumeMapperGPU->SetInput(this->LoadedImages[i]->ImageData);
	this->LoadedImages[i]->volumeMapperGPU->SetMaxMemoryInBytes(std::min(max_memory, 1.9 * 1024 * 1024 * 1024));
	this->LoadedImages[i]->volumeMapperGPU->SetMaxMemoryFraction(1.0);
	this->LoadedImages[i]->volumeMapperGPU->SetSampleDistance((float)this->RaycastSampleDist);
	this->LoadedImages[i]->volumeMapperGPU->SetBlendModeToComposite();
		
	std::cout << "Maximum GPU Memory: " << this->LoadedImages[i]->volumeMapperGPU->GetMaxMemoryInBytes() / (1024 * 1024.0) << " MB" << std::endl;
	std::cout << "Maximum GPU Usage Fraction: " << this->LoadedImages[i]->volumeMapperGPU->GetMaxMemoryFraction() << std::endl;
}
#else
void ImageRenderActors::TextureVolumeMapper(int i)
{
	this->LoadedImages[i]->volumeMapper = vtkSmartPointer<vtkOpenGLVolumeTextureMapper3D>::New();
	this->LoadedImages[i]->volumeMapper->SetSampleDistance((float)this->RaycastSampleDist);
	this->LoadedImages[i]->volumeMapper->SetInput(this->LoadedImages[i]->ImageData);
}
#endif
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
	this->LoadedImages[i]->volume = vtkSmartPointer<vtkVolume>::New();
	if (this->LoadedImages[i]->tag.compare("Soma")==0)
	{
		//std::cout << "setSomaColor" << std::endl;
		this->setSomaColor(this->somaColorValue);
	}
	else
	{
		this->LoadedImages[i]->volumeProperty->SetColor(this->colorTransferFunction);
	}//Image
	this->LoadedImages[i]->volumeProperty->SetScalarOpacity(this->opacityTransferFunction);
	this->LoadedImages[i]->volumeProperty->SetInterpolationTypeToLinear();
#ifdef USE_GPUREN
		//double max_memory = (5.0 * 1024 * 1024 * 1024) / LoadedImages.size();
		//this->LoadedImages[i]->volumeMapperGPU = vtkSmartPointer<vtkOpenGLGPUVolumeRayCastMapper>::New();
		//this->LoadedImages[i]->volumeMapperGPU->SetInput(this->LoadedImages[i]->ImageData);
		//this->LoadedImages[i]->volumeMapperGPU->SetMaxMemoryInBytes(std::min(max_memory, 1.9 * 1024 * 1024 * 1024));
		//this->LoadedImages[i]->volumeMapperGPU->SetMaxMemoryFraction(1.0);
		//this->LoadedImages[i]->volumeMapperGPU->SetSampleDistance((float)this->RaycastSampleDist);
		//this->LoadedImages[i]->volumeMapperGPU->SetBlendModeToComposite();
		RaycastVolumeMapperGPU(i);
		this->LoadedImages[i]->volume->SetMapper(this->LoadedImages[i]->volumeMapperGPU);
		//std::cout << "Maximum GPU Memory: " << this->LoadedImages[i]->volumeMapperGPU->GetMaxMemoryInBytes() / (1024 * 1024.0) << " MB" << std::endl;
		//std::cout << "Maximum GPU Usage Fraction: " << this->LoadedImages[i]->volumeMapperGPU->GetMaxMemoryFraction() << std::endl;
#else
		//this->LoadedImages[i]->volumeMapper = vtkSmartPointer<vtkOpenGLVolumeTextureMapper3D>::New();
		//this->LoadedImages[i]->volumeMapper->SetSampleDistance((float)this->RaycastSampleDist);
		//this->LoadedImages[i]->volumeMapper->SetInput(this->LoadedImages[i]->ImageData);
		TextureVolumeMapper(i);
		this->LoadedImages[i]->volume->SetMapper(this->LoadedImages[i]->volumeMapper);
#endif
	this->LoadedImages[i]->volume->SetProperty(this->LoadedImages[i]->volumeProperty);
	this->LoadedImages[i]->volume->SetPosition(this->LoadedImages[i]->x, 
		this->LoadedImages[i]->y,this->LoadedImages[i]->z);
	this->LoadedImages[i]->volume->SetPickable(0);
	return this->LoadedImages[i]->volume;
}
//vtkSmartPointer<vtkVolume> ImageRenderActors::RayCastSomaVolume(int i)
//{
//	if (i == -1)
//	{
//		i = int (this->LoadedImages.size() - 1);
//	}
//	this->LoadedImages[i]->volumeProperty = vtkSmartPointer<vtkVolumeProperty>::New();
//	this->LoadedImages[i]->somaVolume = vtkSmartPointer<vtkVolume>::New();
//	if (this->LoadedImages[i]->tag.compare("Soma")==0)
//	{
//		//std::cout << "setSomaColor" << std::endl;
//		this->setSomaColor(this->somaColorValue);
//	}
//	else
//	{
//		this->LoadedImages[i]->volumeProperty->SetColor(this->colorTransferFunction);
//	}//Image
//	this->LoadedImages[i]->volumeProperty->SetScalarOpacity(this->opacityTransferFunction);
//	this->LoadedImages[i]->volumeProperty->SetInterpolationTypeToLinear();
//#ifdef USE_GPUREN
//	{
//		RaycastVolumeMapperGPU(i);
//		this->LoadedImages[i]->somaVolume->SetMapper(this->LoadedImages[i]->volumeMapperGPU);
//	}
//#else
//	{
//		TextureVolumeMapper(i);
//		this->LoadedImages[i]->somaVolume->SetMapper(this->LoadedImages[i]->volumeMapper);
//	}
//#endif
//	this->LoadedImages[i]->somaVolume->SetProperty(this->LoadedImages[i]->volumeProperty);
//	this->LoadedImages[i]->somaVolume->SetPosition(this->LoadedImages[i]->x,this->LoadedImages[i]->y,this->LoadedImages[i]->z);
//	this->LoadedImages[i]->somaVolume->SetPickable(0);
//	return this->LoadedImages[i]->somaVolume;
//}
void ImageRenderActors::setSomaColor(double colorValue) //works when raycast mode button triggered
{
	//std::cout << "setSomaColor" << std::endl;
	double r, g ,b;
	this->somaColorValue = colorValue;
	if (colorValue < 0.5)
	{
		r = 1-(colorValue/0.5);
		g = colorValue/0.5;
		b = 0;
	}
	else
	{
		r = 0;
		g = 1-(colorValue-0.5)/0.5;
		b = (colorValue-0.5)/0.5;
	}
	vtkSmartPointer<vtkColorTransferFunction> colorTransferFunctionSoma = vtkColorTransferFunction::New();
	colorTransferFunctionSoma->AddRGBPoint((this->b*this->brightness)/100, r, g, b); // what's this->b and others?

	for (unsigned int i = 0; i< this->LoadedImages.size(); i++)
	{
		if (this->LoadedImages[i]->tag.compare("Soma")==0)
		{			
			if (this->LoadedImages[i]->volume != 0)
			{
				this->LoadedImages[i]->volumeProperty->SetColor(colorTransferFunctionSoma);
				this->LoadedImages[i]->volume->Update();
			}
		}
	}
}

vtkSmartPointer<vtkVolume> ImageRenderActors::GetRayCastVolume(int i)
{
	if (i == -1)
	{
		i = int (this->LoadedImages.size() - 1);
	}
	return this->LoadedImages[i]->volume;
}
//vtkSmartPointer<vtkVolume> ImageRenderActors::GetRayCastSomaVolume(int i)
//{
//	if (i == -1)
//	{
//		i = int (this->LoadedImages.size() - 1);
//	}
//	return this->LoadedImages[i]->somaVolume;
//}
bool ImageRenderActors::getRenderStatus(int i)
{
	if (i == -1)
	{
		i = int (this->LoadedImages.size() - 1);
	}
	return this->LoadedImages[i]->renderStatus;
}
bool ImageRenderActors::is2D(int i)
{
	if (i == -1)
	{
		i = int (this->LoadedImages.size() - 1);
	}
	return this->LoadedImages[i]->ren2d;
}
void ImageRenderActors::setIs2D(int i, bool Set2D)
{
	if (i == -1)
	{
		i = int (this->LoadedImages.size() - 1);
	}
	this->LoadedImages[i]->ren2d = Set2D;
}
void ImageRenderActors::setRenderStatus(int i, bool setStatus)
{
	if (i == -1)
	{
		i = int (this->LoadedImages.size() - 1);
	}
	this->LoadedImages[i]->renderStatus = setStatus;
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
	this->syncColorTransferFunction();
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
	this->syncColorTransferFunction();
}
void ImageRenderActors::setBrightness(int value)
{
	this->brightness = (double)value;
	this->syncColorTransferFunction();
}
int ImageRenderActors::getBrightness()
{
	return (int) this->brightness;
}
void ImageRenderActors::setOpacity(int value)
{
	this->opacity1 = (double ) value;
	this->syncOpacityTransferFunction();
}
int ImageRenderActors::getOpacity()
{
	return (int) this->opacity1;
}
void ImageRenderActors::setOpacityValue(double opacity)
{
	this->opacity1Value = opacity;
	this->syncOpacityTransferFunction();
}
double ImageRenderActors::getOpacityValue()
{
	return this->opacity1Value;
}
void ImageRenderActors::setColorValues(int value)
{
	this->colorValue = value;
	this->syncColorTransferFunction();
}
void ImageRenderActors::syncColorTransferFunction()
{
	this->colorTransferFunction->RemoveAllPoints();
	this->colorTransferFunction->AddRGBPoint(0.0, 0.0, 0.0, 0.0);

	ColorType colorPoint1, colorPoint2, colorPoint3;
	
	//Make the 3 points to be added to the color transformer. 
	//These 3 RGB points describe 3 curves that converts how to tranform the original RGB values (which are grayscale) to their new values.
	switch (colorValue)
	{
		case 1: //red
			colorPoint1.red = 0.11;		colorPoint1.green = 0.0;		colorPoint1.blue = 0.0; 
			colorPoint2.red = 0.59;		colorPoint2.green = 0.0;		colorPoint2.blue = 0.0;
			colorPoint3.red = 0.3;		colorPoint3.green = 0.0;		colorPoint3.blue = 0.0;
			break;
		case 2: //green
			colorPoint1.red = 0.0;		colorPoint1.green = 0.11;		colorPoint1.blue = 0.0;
			colorPoint2.red = 0.0;		colorPoint2.green = 0.59;		colorPoint2.blue = 0.0;
			colorPoint3.red = 0.0;		colorPoint3.green = 0.3;		colorPoint3.blue = 0.0;
			break;
		case 3: //blue
			colorPoint1.red = 0.0;		colorPoint1.green = 0.0;		colorPoint1.blue = 0.11;
			colorPoint2.red = 0.0;		colorPoint2.green = 0.0;		colorPoint2.blue = 0.59;
			colorPoint3.red = 0.0;		colorPoint3.green = 0.0;		colorPoint3.blue = 0.3;
			break;
		case 4: //gray
			colorPoint1.red = 0.75;		colorPoint1.green = 0.75;		colorPoint1.blue = 0.75;
			colorPoint2.red = 0.75;		colorPoint2.green = 0.75;		colorPoint2.blue = 0.75;
			colorPoint3.red = 0.75;		colorPoint3.green = 0.75;		colorPoint3.blue = 0.75;
			break;
		default:
			colorPoint1.red = 0.0;		colorPoint1.green = 0.0;		colorPoint1.blue = 0.5;
			colorPoint2.red = 0.0;		colorPoint2.green = 0.5;		colorPoint2.blue = 0.0;
			colorPoint3.red = 0.5;		colorPoint3.green = 0.0;		colorPoint3.blue = 0.0;
	}
	
	//Add those 3 points into the color transformer
	this->colorTransferFunction->AddRGBPoint((this->b*this->brightness)/100, colorPoint1.red, colorPoint1.green, colorPoint1.blue);
	this->colorTransferFunction->AddRGBPoint((this->g*this->brightness)/100, colorPoint2.red, colorPoint2.green, colorPoint2.blue);
	this->colorTransferFunction->AddRGBPoint((this->r*this->brightness)/100, colorPoint3.red, colorPoint3.green, colorPoint3.blue);

	/*double blueChannelcolor[3] = {0.0,0.0,0.0};
	double greenChannelcolor[3] = {0.0,0.0,0.0};
	double redChannelcolor[3] = {0.0,0.0,0.0};*/

	//switch (colorValue)
	//{
	//	case 1: //red
	//		blueChannelcolor[0] = .11;
	//		greenChannelcolor[0] = .59;
	//		redChannelcolor[0] = .3;
	//		break;
	//	case 2: //green
	//		blueChannelcolor[1] = .11;
	//		greenChannelcolor[1] = .59;
	//		redChannelcolor[1] = .3;
	//		break;
	//	case 3: //blue
	//		blueChannelcolor[2] = .11;
	//		greenChannelcolor[2] = .59;
	//		redChannelcolor[2] = .3;
	//		break;
	//	case 4: //gray
	//		blueChannelcolor[0] = blueChannelcolor[1] = blueChannelcolor[2] = .75;
	//		greenChannelcolor[0] = greenChannelcolor[1] = greenChannelcolor[2] = .75;
	//		redChannelcolor[0] = redChannelcolor[1] = redChannelcolor[2] = .75;
	//		break;
	//	default: 
	//		blueChannelcolor[2] = .5;
	//		greenChannelcolor[1] = .5;
	//		redChannelcolor[0] = .5;
	//		break; 
	//}
	//this->colorTransferFunction->AddRGBPoint((this->b*this->brightness)/100, blueChannelcolor[0],blueChannelcolor[1],blueChannelcolor[2]);//blue
	//this->colorTransferFunction->AddRGBPoint((this->g*this->brightness)/100, greenChannelcolor[0],greenChannelcolor[1],greenChannelcolor[2]);//green
	//this->colorTransferFunction->AddRGBPoint((this->r*this->brightness)/100, redChannelcolor[0], redChannelcolor[1], redChannelcolor[2]);//red
	for (unsigned int i = 0; i< this->LoadedImages.size(); i++)
	{
		if (this->LoadedImages[i]->volume != 0)
			this->LoadedImages[i]->volume->Update();
	}
}
void ImageRenderActors::syncOpacityTransferFunction()
{
	this->opacityTransferFunction->RemoveAllPoints();
	this->opacityTransferFunction->AddPoint(2,0.0);
	this->opacityTransferFunction->AddPoint(this->opacity1,this->opacity1Value);
	//this->opacityTransferFunction->AddPoint(this->opacity2,this->opacity2Value);
	for (unsigned int i = 0; i< this->LoadedImages.size(); i++)
	{
    if(this->LoadedImages[i]->volume)
		{
		this->LoadedImages[i]->volume->Update();
		}
	}
}
void ImageRenderActors::CreateImageResliceMapper(int i)
{
	this->LoadedImages[i]->imageResliceMapper = vtkImageResliceMapper::New();
	this->LoadedImages[i]->imageResliceMapper->SetInput(this->LoadedImages[i]->ImageData);
	//this->LoadedImages[i]->imageResliceMapper->SliceFacesCameraOn();
	this->LoadedImages[i]->imageResliceMapper->SliceAtFocalPointOn();
	//this->LoadedImages[i]->imageResliceMapper->SetSlicePlane(...);
	this->LoadedImages[i]->imageResliceMapper->SetSlabThickness(9);
	this->LoadedImages[i]->imageResliceMapper->SetSlabTypeToMax();
	this->LoadedImages[i]->imageResliceMapper->ResampleToScreenPixelsOff();
}
void ImageRenderActors::CreateImageProperty(int i)
{
	this->LoadedImages[i]->imageProperty = vtkImageProperty::New();
	this->LoadedImages[i]->imageProperty->SetColorWindow(2000); //set range of brightness
	this->LoadedImages[i]->imageProperty->SetColorLevel(this->sliceBrightness); //level of brightness within the range defined by colorwindow
	this->LoadedImages[i]->imageProperty->SetInterpolationTypeToLinear();
}
void ImageRenderActors::SetImageSliceWindowLevel(int value)
{
	//std::cout << "Slice Brightness: " << value << std::endl;
	this->sliceBrightness = value;
	for (unsigned int i = 0; i< this->LoadedImages.size(); i++)
	{
		this->LoadedImages[i]->imageProperty->SetColorLevel(sliceBrightness);
		this->LoadedImages[i]->imageSlice->Update();
	}
}
void ImageRenderActors::CreateImageSlice(int i)
{
	if (i == -1)
	{
		i = int (this->LoadedImages.size() - 1);
	}
	//this->LoadedImages[i]->ImageData->GetExtent(sliceBounds);

	this->LoadedImages[i]->imageSlice = vtkSmartPointer<vtkImageSlice>::New();
	CreateImageResliceMapper(i);
	CreateImageProperty(i);

	this->LoadedImages[i]->imageSlice->SetMapper(this->LoadedImages[i]->imageResliceMapper);
	this->LoadedImages[i]->imageSlice->SetProperty(this->LoadedImages[i]->imageProperty);
	this->LoadedImages[i]->imageSlice->SetPosition(this->LoadedImages[i]->x,this->LoadedImages[i]->y,this->LoadedImages[i]->z);
	this->LoadedImages[i]->imageSlice->SetPickable(0);
	this->LoadedImages[i]->sliceCreated = true;

	this->LoadedImages[i]->imageSlice->GetBounds(sliceBounds);
}
ImageSlicePointerType ImageRenderActors::GetImageSlice(int i)
{
	if (i == -1)
	{
		i = int (this->LoadedImages.size() - 1);
	}
	if (!LoadedImages[i]->sliceCreated)
	{
		CreateImageSlice(i);
		//std::cout << "Slice is created." << std::endl;
	}
	else
	{
		//std::cout << "Slice is already created." << std::endl;
	}
	return this->LoadedImages[i]->imageSlice;
}
void ImageRenderActors::SetSliceThickness(int numofslices)
{
	for (unsigned int i = 0; i< this->LoadedImages.size(); i++)
	{
		this->LoadedImages[i]->imageResliceMapper->SetSlabThickness(numofslices);
		this->LoadedImages[i]->imageSlice->Update();
	}
}
void ImageRenderActors::SetSlicePlane(int slicePlane)
{
	vtkPlane * xyzPlane = vtkPlane::New();
	xyzPlane->SetOrigin(0,0,0);

	switch (slicePlane) 
	{
		case 0: xyzPlane->SetNormal(0,0,1); break; //xy plane
		case 1: xyzPlane->SetNormal(0,1,0); break; //xz plane
		case 2: xyzPlane->SetNormal(1,0,0); break; //yz plane
		default: std::cerr << "View3D::rotateImage cannot handle axis = " << slicePlane << "." << std::endl;
	}

	for (unsigned int i = 0; i< this->LoadedImages.size(); i++)
	{
		this->LoadedImages[i]->imageResliceMapper->SetSlicePlane(xyzPlane);
		this->LoadedImages[i]->imageResliceMapper->SliceAtFocalPointOn();
		this->LoadedImages[i]->imageSlice->SetMapper(this->LoadedImages[i]->imageResliceMapper);
	}
}
std::vector<int> ImageRenderActors::MinCurrentMaxSlices(int i)
{
	if (i == -1)
	{
		i = int (this->LoadedImages.size() - 1);
	}
	std::vector<int> slices;
	//slices.push_back(this->LoadedImages[i]->sliceActor->GetSliceNumberMin()); // These need numbers from display extent
	//slices.push_back(this->LoadedImages[i]->sliceActor->GetZSlice());
	//slices.push_back(this->LoadedImages[i]->sliceActor->GetSliceNumberMax());
	
	slices.push_back(0); // Audrey test
	slices.push_back(0);
	slices.push_back(this->sliceBounds[5]);
	
	return slices;
}
//This function is turned off in View3D and not connected.
void ImageRenderActors::SetSliceNumber(int i, int num) // i is number of images, num is the chosen z slice
{
	if (i == -1)
	{
		i = int (this->LoadedImages.size() - 1);
	}
	if (!this->LoadedImages[i]->ren2d)
	{
		return;
	}
	std::vector<int> slices;
	slices = this->MinCurrentMaxSlices(i);
	//std::cout << num << " min " << slices[0] << " current " << slices[1] << " max " << slices[2] << std::endl;

	//if (slices[2] < num)
	//{
	//	this->LoadedImages[i]->sliceActor->SetZSlice(slices[2]);
	//}
	//else if (slices[0] > num)
	//{
	//	this->LoadedImages[i]->sliceActor->SetZSlice(slices[0]);
	//}
	//else
	//{
	//	std::cout << "Extent: " << LoadedImages[i]->ImageData->GetExtent()[0] << " " << LoadedImages[i]->ImageData->GetExtent()[1] <<  " " << LoadedImages[i]->ImageData->GetExtent()[2] << " " <<  LoadedImages[i]->ImageData->GetExtent()[3] << " " <<  LoadedImages[i]->ImageData->GetExtent()[4]<< " " <<  LoadedImages[i]->ImageData->GetExtent()[5] << std::endl;
	//	this->LoadedImages[i]->sliceActor->SetZSlice(num);
	//}
	//this->LoadedImages[i]->imageSlice->set
}
vtkSmartPointer<vtkImageActor> ImageRenderActors::createProjection(int i, int method, int projection_dim)
{
	//std::cout << "projection_dim: " << projection_dim << std::endl;
	this->LoadedImages[i]->ProjectionActor = vtkSmartPointer<vtkImageActor>::New();
	if (method == 2)
	{
		this->LoadedImages[i]->MinProjection = MinProjectionType::New();
		this->LoadedImages[i]->MinProjection->SetProjectionDimension(projection_dim);
		this->LoadedImages[i]->MinProjection->SetNumberOfThreads(16);
		this->LoadedImages[i]->MinProjection->SetInput(this->LoadedImages[i]->Rescale->GetOutput());
		//this->LoadedImages[i]->projectionConnector->SetInput(this->LoadedImages[i]->MinProjection->GetOutput());
				
		this->LoadedImages[i]->MinProjection->Update();
		ImageType2D::Pointer Min_proj_2D_image = LoadedImages[i]->MinProjection->GetOutput();
		
		//Make destination image		
		ImageType::Pointer Min_proj_3D_image = ImageType::New();
		ImageType::RegionType region;
		ImageType::IndexType start;
		start[0] = 0;
		start[1] = 0;
		start[2] = 0;

		ImageType2D::SizeType input_size = LoadedImages[i]->MinProjection->GetOutput()->GetLargestPossibleRegion().GetSize();
		ImageType::SizeType output_size;

		output_size[0] = input_size[0];
		output_size[1] = input_size[1];
		output_size[2] = 1;

		region.SetSize(output_size);
		region.SetIndex(start);

		Min_proj_3D_image->SetRegions(region);
		Min_proj_3D_image->Allocate();

		itk::ImageRegionIterator<ImageType2D>	Min_proj_2D_image_iterator(Min_proj_2D_image, Min_proj_2D_image->GetRequestedRegion());
		itk::ImageRegionIterator<ImageType>		Min_proj_3D_image_iterator(Min_proj_3D_image, Min_proj_3D_image->GetRequestedRegion());

		for (Min_proj_2D_image_iterator = Min_proj_2D_image_iterator.Begin(); !Min_proj_2D_image_iterator.IsAtEnd(); ++Min_proj_2D_image_iterator, ++Min_proj_3D_image_iterator)
			Min_proj_3D_image_iterator.Set(Min_proj_2D_image_iterator.Get());

		//Permute image axes
		LoadedImages[i]->itkPermute = itkPermuteFilterType::New();
		unsigned int permuteAxes[3];
		minXBound = 0; maxXBound = 0;
		minYBound = 0; maxYBound = 0;
		minZBound = 0; maxZBound = 0;
		switch(projection_dim)
		{
			case 0:
				permuteAxes[0] = 2;
				permuteAxes[1] = 1;
				permuteAxes[2] = 0;
				minXBound = 0; maxXBound = 0;
				minYBound = 0; maxYBound = output_size[1]-1;
				minZBound = 0; maxZBound = output_size[0]-1;
				break;
			case 1:
				permuteAxes[0] = 0;
				permuteAxes[1] = 2;
				permuteAxes[2] = 1;
				minXBound = 0; maxXBound = output_size[0]-1;
				minYBound = 0; maxYBound = 0;
				minZBound = 0; maxZBound = output_size[1]-1;
				break;
			case 2:
				permuteAxes[0] = 0;
				permuteAxes[1] = 1;
				permuteAxes[2] = 2;
				minXBound = 0; maxXBound = output_size[0]-1;
				minYBound = 0; maxYBound = output_size[1]-1;
				minZBound = 0; maxZBound = 0;
				break;
			default:
				permuteAxes[0] = 0;
				permuteAxes[1] = 1;
				permuteAxes[2] = 2;
				break;
		}
		LoadedImages[i]->itkPermute->SetOrder(permuteAxes);
		LoadedImages[i]->itkPermute->SetInput(Min_proj_3D_image);

		this->LoadedImages[i]->projectionConnector->SetInput(LoadedImages[i]->itkPermute->GetOutput());
	}
	else if (method == 1)
	{
		this->LoadedImages[i]->MeanProjection = MeanProjectionType::New();
		this->LoadedImages[i]->MeanProjection->SetProjectionDimension(projection_dim);
		this->LoadedImages[i]->MeanProjection->SetNumberOfThreads(16);
		this->LoadedImages[i]->MeanProjection->SetInput(this->LoadedImages[i]->Rescale->GetOutput());

		this->LoadedImages[i]->MeanProjection->Update();
		ImageType2D::Pointer mean_proj_2D_image = LoadedImages[i]->MeanProjection->GetOutput();
		
		//Make destination image		
		ImageType::Pointer mean_proj_3D_image = ImageType::New();
		ImageType::RegionType region;
		ImageType::IndexType start;
		start[0] = 0;
		start[1] = 0;
		start[2] = 0;

		ImageType2D::SizeType input_size = LoadedImages[i]->MeanProjection->GetOutput()->GetLargestPossibleRegion().GetSize();
		ImageType::SizeType output_size;

		output_size[0] = input_size[0];
		output_size[1] = input_size[1];
		output_size[2] = 1;

		region.SetSize(output_size);
		region.SetIndex(start);

		mean_proj_3D_image->SetRegions(region);
		mean_proj_3D_image->Allocate();

		itk::ImageRegionIterator<ImageType2D>	mean_proj_2D_image_iterator(mean_proj_2D_image, mean_proj_2D_image->GetRequestedRegion());
		itk::ImageRegionIterator<ImageType>		mean_proj_3D_image_iterator(mean_proj_3D_image, mean_proj_3D_image->GetRequestedRegion());

		for (mean_proj_2D_image_iterator = mean_proj_2D_image_iterator.Begin(); !mean_proj_2D_image_iterator.IsAtEnd(); ++mean_proj_2D_image_iterator, ++mean_proj_3D_image_iterator)
			mean_proj_3D_image_iterator.Set(mean_proj_2D_image_iterator.Get());

		//Permute image axes
		LoadedImages[i]->itkPermute = itkPermuteFilterType::New();
		unsigned int permuteAxes[3];
		minXBound = 0; maxXBound = 0;
		minYBound = 0; maxYBound = 0;
		minZBound = 0; maxZBound = 0;
		switch(projection_dim)
		{
			case 0:
				permuteAxes[0] = 2;
				permuteAxes[1] = 1;
				permuteAxes[2] = 0;
				minXBound = 0; maxXBound = 0;
				minYBound = 0; maxYBound = output_size[1]-1;
				minZBound = 0; maxZBound = output_size[0]-1;
				break;
			case 1:
				permuteAxes[0] = 0;
				permuteAxes[1] = 2;
				permuteAxes[2] = 1;
				minXBound = 0; maxXBound = output_size[0]-1;
				minYBound = 0; maxYBound = 0;
				minZBound = 0; maxZBound = output_size[1]-1;
				break;
			case 2:
				permuteAxes[0] = 0;
				permuteAxes[1] = 1;
				permuteAxes[2] = 2;
				minXBound = 0; maxXBound = output_size[0]-1;
				minYBound = 0; maxYBound = output_size[1]-1;
				minZBound = 0; maxZBound = 0;
				break;
			default:
				permuteAxes[0] = 0;
				permuteAxes[1] = 1;
				permuteAxes[2] = 2;
				break;
		}
		LoadedImages[i]->itkPermute->SetOrder(permuteAxes);
		LoadedImages[i]->itkPermute->SetInput(mean_proj_3D_image);

		this->LoadedImages[i]->projectionConnector->SetInput(LoadedImages[i]->itkPermute->GetOutput());
	}
	else // 0 is maximum projection
	{
		this->LoadedImages[i]->MaxProjection = MaxProjectionType::New();
		this->LoadedImages[i]->MaxProjection->SetProjectionDimension(projection_dim);
		this->LoadedImages[i]->MaxProjection->SetNumberOfThreads(16);
		this->LoadedImages[i]->MaxProjection->SetInput(this->LoadedImages[i]->Rescale->GetOutput());
		
		/*LoadedImages[i]->Rescale->Update();
		ImageType::Pointer max_proj_input_image = LoadedImages[i]->Rescale->GetOutput();
		std::cout << "Max Projection Input Image Size: " << max_proj_input_image->GetLargestPossibleRegion().GetSize() << std::endl;*/
		
		this->LoadedImages[i]->MaxProjection->Update();
		ImageType2D::Pointer max_proj_2D_image = LoadedImages[i]->MaxProjection->GetOutput();
		//std::cout << "Max Projection Output Image Size: " << max_proj_2D_image->GetLargestPossibleRegion().GetSize() << std::endl;
		
		//Make destination image		
		ImageType::Pointer max_proj_3D_image = ImageType::New();
		ImageType::RegionType region;
		ImageType::IndexType start;
		start[0] = 0;
		start[1] = 0;
		start[2] = 0;

		ImageType2D::SizeType input_size = LoadedImages[i]->MaxProjection->GetOutput()->GetLargestPossibleRegion().GetSize();
		ImageType::SizeType output_size;

		output_size[0] = input_size[0];
		output_size[1] = input_size[1];
		output_size[2] = 1;

		region.SetSize(output_size);
		region.SetIndex(start);

		max_proj_3D_image->SetRegions(region);
		max_proj_3D_image->Allocate();

		itk::ImageRegionIterator<ImageType2D>	max_proj_2D_image_iterator(max_proj_2D_image, max_proj_2D_image->GetRequestedRegion());
		itk::ImageRegionIterator<ImageType>		max_proj_3D_image_iterator(max_proj_3D_image, max_proj_3D_image->GetRequestedRegion());

		for (max_proj_2D_image_iterator = max_proj_2D_image_iterator.Begin(); !max_proj_2D_image_iterator.IsAtEnd(); ++max_proj_2D_image_iterator, ++max_proj_3D_image_iterator)
			max_proj_3D_image_iterator.Set(max_proj_2D_image_iterator.Get());

		//Permute image axes
		LoadedImages[i]->itkPermute = itkPermuteFilterType::New();
		unsigned int permuteAxes[3];
		minXBound = 0; maxXBound = 0;
		minYBound = 0; maxYBound = 0;
		minZBound = 0; maxZBound = 0;
		switch(projection_dim)
		{
			case 0:
				permuteAxes[0] = 2;
				permuteAxes[1] = 1;
				permuteAxes[2] = 0;
				minXBound = 0; maxXBound = 0;
				minYBound = 0; maxYBound = output_size[1]-1;
				minZBound = 0; maxZBound = output_size[0]-1;
				break;
			case 1:
				permuteAxes[0] = 0;
				permuteAxes[1] = 2;
				permuteAxes[2] = 1;
				minXBound = 0; maxXBound = output_size[0]-1;
				minYBound = 0; maxYBound = 0;
				minZBound = 0; maxZBound = output_size[1]-1;
				break;
			case 2:
				permuteAxes[0] = 0;
				permuteAxes[1] = 1;
				permuteAxes[2] = 2;
				minXBound = 0; maxXBound = output_size[0]-1;
				minYBound = 0; maxYBound = output_size[1]-1;
				minZBound = 0; maxZBound = 0;
				break;
			default:
				permuteAxes[0] = 0;
				permuteAxes[1] = 1;
				permuteAxes[2] = 2;
				break;
		}
		LoadedImages[i]->itkPermute->SetOrder(permuteAxes);
		LoadedImages[i]->itkPermute->SetInput(max_proj_3D_image);

		this->LoadedImages[i]->projectionConnector->SetInput(LoadedImages[i]->itkPermute->GetOutput());
		//this->LoadedImages[i]->projectionConnector->SetInput(max_proj_3D_image);
		//std::cout << "3D image size" << max_proj_3D_image->GetLargestPossibleRegion().GetSize() << std::endl;
		//max_proj_3D_image->GetLargestPossibleRegion().GetSize()[0]
		//std::cout << "ItkPermute size" << LoadedImages[i]->itkPermute << std::endl;
		
		//vtkSmartPointer<vtkImageEllipsoidSource > source = vtkSmartPointer<vtkImageEllipsoidSource >::New();
		//source->SetWholeExtent(0, 20, 0, 20, 0, 0);
		//source->SetCenter(10,10,0);
		//source->SetRadius(2,5,0);
		//source->Update();
		
		//LoadedImages[i]->PermuteFilter = PermuteFilterType::New();
		//LoadedImages[i]->PermuteFilter->SetInputConnection(source->GetOutputPort());

		////LoadedImages[i]->PermuteFilter->SetInputConnection(source->GetOutputPort());
		//LoadedImages[i]->PermuteFilter->SetFilteredAxes(0,2,1);
		//try
		//{
		//	LoadedImages[i]->PermuteFilter->Update();
		//}
		//catch (itk::ExceptionObject &err)
		//{
		//	std::cout << err << std::endl;
		//	return NULL;
		//}	
		////double center[3];


		//static double xyplane[16] = {
		//		1, 0, 0, 0,
		//		0, 1, 0, 0,
		//		0, 0, 1, 0,
		//		0, 0, 0, 1 };
		//static double yzplane[16] = {
		//		0, 0,-1, 0,
		//		1, 0, 0, 0,
		//		0,-1, 0, 0,
		//		0, 0, 0, 1 };

		//// Set slice orientation
		//LoadedImages[i]->resliceAxes = vtkSmartPointer<vtkMatrix4x4>::New();
		//LoadedImages[i]->resliceAxes->DeepCopy(xyplane);
		//// Set the point through which to slice
		////LoadedImages[i]->resliceAxes->SetElement(0, 3, center[0]);
		////LoadedImages[i]->resliceAxes->SetElement(1, 3, center[1]);
		////LoadedImages[i]->resliceAxes->SetElement(2, 3, center[2]);
		//LoadedImages[i]->resliceAxes->SetElement(0, 3, output_size[0]/2);
		//LoadedImages[i]->resliceAxes->SetElement(1, 3, output_size[1]/2);
		//LoadedImages[i]->resliceAxes->SetElement(2, 3, 0);

		//// Extract slice in the desired orientation
		//LoadedImages[i]->reslice = vtkSmartPointer<vtkImageReslice>::New();
		//this->LoadedImages[i]->PermuteFilter->Update();
		//LoadedImages[i]->reslice->SetInput(this->LoadedImages[i]->PermuteFilter->GetOutput());
		//LoadedImages[i]->reslice->SetResliceAxes(LoadedImages[i]->resliceAxes);
		//LoadedImages[i]->reslice->Update();
		//int* extent_output = LoadedImages[i]->PermuteFilter->GetOutputExtent();
		//std::cout << extent_output[0] << std::endl;
	}

	//this->LoadedImages[i]->ProjectionActor->SetInput(this->LoadedImages[i]->reslice->GetOutput());
	this->LoadedImages[i]->ProjectionActor->SetInput(this->LoadedImages[i]->projectionConnector->GetOutput());
	this->LoadedImages[i]->ProjectionActor->SetPosition(this->LoadedImages[i]->x, this->LoadedImages[i]->y, this->LoadedImages[i]->z);
	this->LoadedImages[i]->ProjectionActor->SetPickable(0);
	this->LoadedImages[i]->ProjectionActor->SetDisplayExtent(minXBound,maxXBound,minYBound,maxYBound,minZBound,maxZBound);

	return this->LoadedImages[i]->ProjectionActor;
}
vtkSmartPointer<vtkImageActor> ImageRenderActors::GetProjectionImage(int i)
{
	if (i == -1)
	{
		i = int (this->LoadedImages.size() - 1);
	}
	return this->LoadedImages[i]->ProjectionActor;
}
void ImageRenderActors::setImageBounds(double bounds[])
{ 
	if (this->TotalImageSize[0] > bounds[0])
	{
		this->TotalImageSize[0] = bounds[0];
	}
	if (this->TotalImageSize[1] < bounds[1])
	{
		this->TotalImageSize[1] = bounds[1];
	}
	if (this->TotalImageSize[2] > bounds[2])
	{
		this->TotalImageSize[2] = bounds[2];
	}
	if (this->TotalImageSize[3] < bounds[3])
	{
		this->TotalImageSize[3] = bounds[3];
	}
	if (this->TotalImageSize[4] > bounds[4])
	{
		this->TotalImageSize[4] = bounds[4];
	}
	if (this->TotalImageSize[5] < bounds[5])
	{
		this->TotalImageSize[5] = bounds[5];
	}
}
void ImageRenderActors::getImageBounds(double bounds[])
{
	for (int i = 0; i < 6; i++)
	{
		bounds[i] = this->TotalImageSize.at(i);
	}
}
//Xmin,Xmax,Ymin,Ymax,Zmin,Zmax
double* ImageRenderActors::getSliceBounds()
{
	return sliceBounds;
}
