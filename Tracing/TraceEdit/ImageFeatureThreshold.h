 #ifndef IMAGEFEATURETHRESHOLD_H_
 #define IMAGEFEATURETHRESHOLD_H_

#include <QtGui>
#include <QApplication>
#include <QThread>
#include "itkImageFileReader.h"
#include "itkBinaryThresholdImageFilter.h"
#include <itkRegionOfInterestImageFilter.h>
#include <itkCastImageFilter.h>
#include "itkImageFileWriter.h"


#include <iostream>
#include "time.h"

#include "vtkRenderer.h" 
#include <QVTKWidget.h>
#include "vtkSmartPointer.h"
#include "vtkSmartVolumeMapper.h"
#include "vtkVolume.h"
#include "vtkVolumeProperty.h"








#include <vtkSmartPointer.h>
#include <vtkCubeSource.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkActor2D.h>
#include <vtkPolyDataMapper2D.h>
#include <vtkProperty2D.h>
#include "itkImageToVTKImageFilter.h"
#include "vtkActor.h"
#include "vtkContourFilter.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkColorTransferFunction.h"



#include "vtkPiecewiseFunction.h"






#include "itkImage.h"
#include "itkCovariantVector.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkGradientImageFilter.h"
 
#include <itkImageToVTKImageFilter.h>
 
#include "vtkPointData.h"
#include "vtkFloatArray.h"
#include "vtkImageData.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderWindow.h"
#include "vtkSmartPointer.h"
#include "vtkImageActor.h"
#include "vtkActor.h"
#include "vtkInteractorStyleImage.h"
#include "vtkRenderer.h"
#include "vtkGlyph3DMapper.h"
#include "vtkArrowSource.h"

#include "ImageActors.h"
#include <QVTKWidget.h>
#include <QMainWindow>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkImageViewer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleImage.h>
#include <vtkRenderer.h>
#include <vtkJPEGReader.h>
#include "vtkTable.h"
#include <QVTKWidget.h>


//struct imageFeatureThreshold
//{
//	LabelImageType3D::Pointer somaImage;
//	ImageType3D::Pointer inputImage;
//	vtkSmartPointer<vtkTable> features;
//	vtkSmartPointer< vtkTable > SWCTable;
//	int costThreshold;
//	float intensity_threshold;
//	float contrast_threshold;
//	int cost_threshold;
//	float debris_threshold;
//	int offshoot;
//	int device;
//	std::string userInput; // yes/no based on whether to select the trace for computation
//};
 class ImageFeatureThreshold 
 {
     
 
	 typedef float PixelType;
	 typedef itk::Image< PixelType, 3 >  ImageType3D;
	 typedef itk::Image< unsigned int, 3 > LabelImageType3D;
	 typedef itk::Image< unsigned char, 3 >   UnsignedCharImageType;
	 typedef itk::RegionOfInterestImageFilter< UnsignedCharImageType, UnsignedCharImageType> RegionOfInterestFilter;

 public:

	ImageFeatureThreshold();
	~ImageFeatureThreshold();
	void saveInputImage(ImageType3D::Pointer &image){_InputImage=image;};
	void saveSomaImage(LabelImageType3D::Pointer &image){_SomaImage=image;};
	void saveFeatures(vtkSmartPointer< vtkTable > features){_Features=features;};
	void saveImageFeature(double featureValue){_ImageFeatures.push_back(featureValue);};
	void saveTraces(vtkSmartPointer< vtkTable > SWCTable){_SWCTable=SWCTable;};
	void saveCostThreshold(int costThreshold){_costThreshold =costThreshold;};
	void saveintensity_threshold(float intensity_threshold){_intensity_threshold=intensity_threshold;};
	void saveContrastThreshold(float contrast_threshold){_contrast_threshold=contrast_threshold;};
	void saveDebrisThreshold(float debris_threshold){_debris_threshold=debris_threshold;};
	void saveoffshoot(int offshoot){_offshoot=offshoot;};
	void savedevice(int device){_device=device;};
	void saveUserInput(std::string userInput){_userInput=userInput;};
	void setSWCFilePath(std::string filePath){_filePath = filePath;};
	void setDirPath(std::string dirLocation){_dirLocation = dirLocation;};
	void setCropImageFilePath(std::string filePath){_cropImageFilePath = filePath;};
	void setCropSomaImageFilePath(std::string filePath){_cropsomaImageFilePath = filePath;};
	void setGVFTracer(bool gvfTracing){_gvfTracing = gvfTracing;};
	void setxoff(std::string ssx_off){_ssx_off = ssx_off;};
	void setyoff(std::string ssy_off){_ssy_off = ssy_off;};
	void setzoff(std::string ssz_off){_ssz_off = ssz_off;};
	ImageType3D::Pointer getInputImage(){return _InputImage;};
	LabelImageType3D::Pointer getSomaImage(){return _SomaImage;};
	void getFeatures(){/*return _Features;*/};
	vtkSmartPointer< vtkTable > getTraces(){return _SWCTable;};
	int getCostThreshold(){return _costThreshold; };
	float getintensity_threshold(){return _intensity_threshold;};
	float getContrastThreshold(){return _contrast_threshold;};
	float getDebrisThreshold(){return _debris_threshold;};
	int getoffshoot(){return _offshoot;};
	int getdevice(){return _device;};
	bool isGVFTracer(){return _gvfTracing;};
	std::string getUserInput(){return _userInput;};
	std::string getSWCFilePath(){return _filePath;};
	std::string getDir(){return _dirLocation;};
	std::string getCropImageFilePath(){return _cropImageFilePath;};
	std::string getCropSomaImageFIlePath(){return _cropsomaImageFilePath;};
	std::string getxoff(){return _ssx_off;};
	std::string getyoff(){return _ssy_off ;};
	std::string getzoff(){return _ssz_off ;};




 protected:
 
 private:
	/*std::vector<imageFeatureThreshold*> computedImageFeatureThreshold;*/
	LabelImageType3D::Pointer _SomaImage;
	ImageType3D::Pointer _InputImage;
	vtkSmartPointer<vtkTable> _Features;
	std::vector<double> _ImageFeatures;//MeanIntenstiy, Intensity SD
	vtkSmartPointer< vtkTable > _SWCTable;
	int _costThreshold;
	float _intensity_threshold;
	float _contrast_threshold;
	float _debris_threshold;
	int _offshoot;
	int _device;
	bool _gvfTracing;
	std::string _filePath;
	std::string _userInput;
	std::string _dirLocation;
	std::string _cropImageFilePath;
	std::string _cropsomaImageFilePath;
	std::string _ssx_off;
	std::string _ssy_off;
	std::string _ssz_off;
 };



 #endif