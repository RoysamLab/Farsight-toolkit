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

#include "vtkActor.h"
#include "vtkContourFilter.h"
#include "vtkColorTransferFunction.h"
#include "vtkImageData.h"
#include "vtkImageToStructuredPoints.h"
#include "vtkImageActor.h"
#include "vtkLODActor.h"
#include "vtkOpenGLVolumeTextureMapper3D.h"
#include "itkMaximumProjectionImageFilter.h"
#include "itkMinimumProjectionImageFilter.h"
#include "itkMeanProjectionImageFilter.h"
#ifdef USE_GPUREN
#include <vtkOpenGLGPUVolumeRayCastMapper.h>
#endif
#include "vtkPiecewiseFunction.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkSmartPointer.h"
#include "vtkVolume.h"
#include "vtkVolumeProperty.h"

#include "itkImageFileReader.h"
#include "itkImageToVTKImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkExtractImageFilter.h"
#include "vtkImagePlaneWidget.h"
#include "vtkImageReslice.h"
#include "itkPermuteAxesImageFilter.h"

//slicer
#include "vtkImageResliceMapper.h"
#include "vtkImageMapper3D.h"
#include "vtkImageSlice.h"
#include "vtkImageProperty.h"
#include "vtkPlane.h"

#include <stdio.h>
#include <string>
#include <vector>
//#include <set>
//#include <QtGui>


const unsigned int Dimension = 3;

/** Standard class typedefs. */

typedef unsigned char  ImageActorPixelType;
typedef itk::Image< ImageActorPixelType, Dimension >   ImageType;
typedef itk::Image< ImageActorPixelType, 2 >   ImageType2D;
typedef itk::ImageFileReader< ImageType >    ReaderType;
typedef itk::ImageToVTKImageFilter<ImageType> ConnectorType;
typedef itk::MaximumProjectionImageFilter < ImageType, ImageType2D> MaxProjectionType;
typedef itk::MinimumProjectionImageFilter < ImageType, ImageType2D> MinProjectionType;
typedef itk::MeanProjectionImageFilter < ImageType, ImageType2D> MeanProjectionType;
typedef itk::RescaleIntensityImageFilter< ImageType, ImageType> IntensityRescaleType;
typedef vtkSmartPointer<vtkImageActor> ImageActorPointerType;
//typedef vtkSmartPointer<vtkImagePermute> PermuteFilterType;
//typedef vtkSmartPointer<vtkImageResliceMapper> ImageResliceMapper;
typedef itk::PermuteAxesImageFilter<ImageType> itkPermuteFilterType;
typedef vtkSmartPointer<vtkImageSlice> ImageSlicePointerType;

/** \struct imageFileHandle
 * \brief Basic data structure for the ImageActors class that holds different data elements used to handle the image
 *
 * The imageFileHandle data structure is the base for the ImageActors class. It stores the different parameters that
 * are needed to handle an image file. This data structure contains parameters to change the render status for the
 * image and specify the number of dimensions it will be rendered in (2D or 3D) amongst others.
 */

struct imageFileHandle
{
	/** A tag for the image. */
	std::string tag;
	
	
	std::string filename;
	bool renderStatus;
	bool ren2d;
	bool sliceCreated;
	vtkSmartPointer<vtkImageData> ImageData;
	ReaderType::Pointer reader;
	ConnectorType::Pointer connector;
	ConnectorType::Pointer projectionConnector;
	MaxProjectionType::Pointer MaxProjection;
	MeanProjectionType::Pointer MeanProjection;
	MinProjectionType::Pointer MinProjection;
	IntensityRescaleType::Pointer Rescale;
	itkPermuteFilterType::Pointer itkPermute;
	double x,y,z;
//!Contour Filter pointers
	vtkSmartPointer<vtkContourFilter> ContourFilter;
	vtkSmartPointer<vtkPolyDataMapper> ContourMapper;
	vtkSmartPointer<vtkActor> ContourActor;
//Raycast pointers
	vtkSmartPointer<vtkVolumeProperty> volumeProperty;
	vtkSmartPointer<vtkOpenGLVolumeTextureMapper3D> volumeMapper;
//image slicer
	vtkSmartPointer<vtkImageActor> ProjectionActor;
	vtkSmartPointer<vtkImageResliceMapper> imageResliceMapper;
	vtkSmartPointer<vtkImageSlice> imageSlice;
	vtkSmartPointer<vtkImageProperty> imageProperty;
#ifdef USE_GPUREN
	vtkSmartPointer<vtkOpenGLGPUVolumeRayCastMapper> volumeMapperGPU;
#endif
	vtkSmartPointer<vtkVolume> volume;
	vtkSmartPointer<vtkVolume> somaVolume;
};
class  ImageRenderActors
{
public:
	ImageRenderActors();
	~ImageRenderActors();
	int loadImage(std::string ImageSource, std::string tag);
	int loadImage(std::string ImageSource, std::string tag, double x, double y, double z);
//render actors
	vtkSmartPointer<vtkActor> ContourActor(int i);
	void setSomaColor(double colorValue);
	vtkSmartPointer<vtkActor> GetContourActor(int i);
	void CreateImageResliceMapper(int i);
	void CreateImageProperty(int i);
	void SetImageSliceWindowLevel(int value);
	void CreateImageSlice(int i);
	ImageSlicePointerType GetImageSlice(int i);
	void SetSliceThickness(int numofslices);
	void SetSlicePlane(int slicePlane);
	vtkSmartPointer<vtkImageActor> createProjection(int i, int method, int projection_dim);
	vtkSmartPointer<vtkImageActor> GetProjectionImage(int i);
	void RaycastVolumeMapperGPU(int i);
	void TextureVolumeMapper(int i);
	vtkSmartPointer<vtkVolume> RayCastVolume(int i);
	vtkSmartPointer<vtkVolume> GetRayCastVolume(int i);
	bool getRenderStatus(int i);
	void setRenderStatus(int i, bool setStatus);
//file information
	std::vector<std::string> GetImageList();
	std::string FileNameOf(int i){ return this->LoadedImages[i]->filename;};
	unsigned int NumberOfImages() {return (unsigned int)this->LoadedImages.size();};
	bool isRayCast(int i);
	bool is2D(int i);
	void setIs2D(int i, bool Set2D);
	void ShiftImage(int i, double x, double y, double z);
	void ShiftImage(int i, std::vector<double> shift);
	std::vector<double> GetShiftImage(int i);
	double pointData(int i, int x, int y, int z);
	std::vector<double> GetImageSize(int i);
	std::vector<int> MinCurrentMaxSlices(int i);
	void SetSliceNumber(int i, int num);
	vtkSmartPointer<vtkImageData> GetImageData(int i);
	std::vector<double> getColorValues();
	void setColorValues(double r, double g, double b);
	void setColorValues(int i, double value);
	void setBrightness(int value);
	int getBrightness();
	void setOpacity(int value);//this is the threshold
	int getOpacity();
	void setOpacityMax(int value);//this is the threshold
	int getOpacityMax();
	void setOpacityValue(double opacity);
	double getOpacityValue();
	void setOpacityValueMax(double opacity);
	double getOpacityValueMax();
	void setColorValues(int value);
	void getImageBounds(double bounds[]);
	void setImageBounds(double bounds[]);
	double* getSliceBounds();
	void setSomaBrightness(int value);
	int getSomaBrightness();
	void setSomaOpacity(int value);//this is the threshold
	int getSomaOpacity();
	void setSomaOpacityMax(int value);//this is the threshold
	int getSomaOpacityMax();
	void setSomaOpacityValue(double value);
	double getSomaOpacityValue();

	struct ColorType
	{
		double red;
		double green;
		double blue;
	};
private:
	bool useGPURendering;
	void syncColorTransferFunction();
	void syncOpacityTransferFunction();
	void syncSomaColorTransferFunction();
	void syncSomaOpacityTransferFunction();
	vtkSmartPointer<vtkPiecewiseFunction> opacityTransferFunction;
	vtkSmartPointer<vtkPiecewiseFunction> opacityTransferFunctionSoma;
	vtkSmartPointer<vtkColorTransferFunction> colorTransferFunction;
	vtkSmartPointer<vtkColorTransferFunction> colorTransferFunctionSoma;
	std::vector<imageFileHandle*> LoadedImages;
	std::vector<std::string> ImageList;
	std::vector<double> TotalImageSize;
	double r,g,b, opacity1, opacity2, opacity1Value, opacity2Value, brightness, RaycastSampleDist;
	double somaColorValue, somaBrightness;
	double somaOpacity, somaOpacityValue;
	int colorValue, sliceBrightness;
	int minXBound, maxXBound, minYBound, maxYBound, minZBound, maxZBound;
	double sliceBounds[6];

};
#endif
