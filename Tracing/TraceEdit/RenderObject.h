 #ifndef RENDEROBJECT_H_
 #define RENDEROBJECT_H_

#include <QtGui>
#include <QApplication>
#include <QThread>
#include "itkImageFileReader.h"
#include "itkBinaryThresholdImageFilter.h"
#include <itkRegionOfInterestImageFilter.h>
#include <itkCastImageFilter.h>
#include "itkImageFileWriter.h"
#include <itkLogImageFilter.h>

#include <iostream>
#include "time.h"

#include "vtkRenderer.h" 
#include <QVTKWidget.h>
#include "vtkSmartPointer.h"
#include "vtkSmartVolumeMapper.h"
#include "vtkVolume.h"
#include "vtkVolumeProperty.h"


#include <boost/lexical_cast.hpp>





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
//#include <itkMultiplyByConstantImageFilter.h>





#include "itkImage.h"
#include "itkCovariantVector.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkGradientImageFilter.h"
 
#include <itkImageToVTKImageFilter.h>
#include <itkPathToChainCodePathFilter.h>
#include <itkImageToPathFilter.h>
#include <itkGaborImageSource.h>
#include "itkConvolutionImageFilter.h"
#include "itkPolyLineParametricPath.h"

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
//#include "ImageViewer.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkOtsuThresholdImageFilter.h"
#include <itkSubtractImageFilter.h>
#include "itkRescaleIntensityImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkAbsImageFilter.h"
#include "ImageFeatureThreshold.h"
#include "Tracing/MultipleNeuronTracer/MultipleNeuronTracer.h"



class RenderObject : public QObject
 {
     Q_OBJECT

	 typedef float PixelType;
	 typedef itk::Image< PixelType, 3 >  ImageType3D;
	 typedef itk::Image< unsigned int, 3 > LabelImageType3D;
	 typedef itk::Image< unsigned char, 3 >   UnsignedCharImageType;
	 typedef itk::RescaleIntensityImageFilter<ImageType3D, ImageType3D> RescalerType;
	 typedef itk::MedianImageFilter<ImageType3D, ImageType3D> MedianFilterType;

	 typedef itk::RegionOfInterestImageFilter< UnsignedCharImageType, UnsignedCharImageType> RegionOfInterestFilter;
	 typedef itk::StatisticsImageFilter<ImageType3D> StatisticsFilterType;
	 //typedef itk::ImageToPathFilter< UnsignedCharImageType, ChainCodePath> ImageToPathFilter;
	 typedef itk::ImageToVTKImageFilter<UnsignedCharImageType> ConnectorType;
	 // to Implement Chain Code
	 typedef itk::Index<3> IndexType;
	 typedef itk::LogImageFilter<ImageType3D,ImageType3D> LogImageFilter;
	 typedef itk::InvertIntensityImageFilter <ImageType3D> InvertIntensityImageFilterType;
	 typedef itk::SubtractImageFilter <ImageType3D, ImageType3D> SubtractImageFilterType;

 public:

	RenderObject(QObject *parent = 0);
	~RenderObject();


	LabelImageType3D::Pointer somaImage;
	ImageType3D::Pointer inputImage;
	ImageType3D::Pointer inputImage_disp;
	ImageType3D::Pointer objectnessImage;
	ImageType3D::Pointer sampledImage;
	ImageType3D::Pointer sampledObjectnessImage;
	ImageType3D::Pointer _PaddedCurvImage;
		float minIntensityThreshold;
	float maxIntensityThreshold;
	float maxcount;
	float minContrastThreshold;
	float maxContrastThreshold;
	
	std::vector<IndexType> SWCPoints;
	std::string imageFile;
	std::string somaFile;
	std::string seedsFile;
	int cost_threshold;
	float intensity_threshold;
	float contrast_threshold;
	std::vector<float> intensityThreshold;
	std::vector<float> contrastThreshold;

	
	float debris_threshold;
	int offshoot;
	int device;

	float widthx;
	float widthy;
	float widthz;

	int startx;
	int starty;
	int startz;
	void ParseFileName(char* fullName);
	void sampleTraces(std::string thresholdVal);
	void writeImage(ImageType3D::Pointer &image, std::string fileName);
	void getSWCPoints(std::string filePath, std::vector<IndexType> &swcPoints);
	void getSiftFeatures(ImageType3D::Pointer &image);
	void GetFeature(float);
	void LoadCurvImage_1(ImageType3D::Pointer &image,unsigned int pad);
	void saveImageStatsFeatures();
	std::vector<float> loadThresholds(float minVal,float maxVal,float count);
	std::vector<float> computeFeatures(ImageType3D::Pointer &image);
	void saveImageFeatures(ImageFeatureThreshold* featureThreshold){_featureThreshold=featureThreshold;};
	void WriteCenterTrace(vtkSmartPointer< vtkTable > , int , int , int , std::string );
	std::vector< itk::Index<3> > getSomaTable( std::vector< itk::Index<3> > centroid_list, int x, int y, int z,int xTile, int yTile, int zTile,int xSize,int ySize, int zSize);
	std::vector< itk::Index<3> > getCentroidList();
signals:
	void finished();
	void finished(ImageFeatureThreshold* featureThreshold);
	
 public slots:
     void runMNTTracer();
	 void displayImage();
	
 protected:
 
 private:
	 ImageFeatureThreshold* _featureThreshold;
 };




 #endif