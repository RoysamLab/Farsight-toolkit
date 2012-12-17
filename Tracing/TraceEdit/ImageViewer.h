 #ifndef IMAGEVIEWER_H_
 #define IMAGEVIEWER_H_

#include <QtGui>
#include <QApplication>
#include <QThread>
#include "itkImageFileReader.h"
#include "itkBinaryThresholdImageFilter.h"
#include <itkRegionOfInterestImageFilter.h>
#include <itkCastImageFilter.h>
#include "itkImageFileWriter.h"
#include "Tracing/MultipleNeuronTracer/MultipleNeuronTracer.h"
#include "CellTraceModel.h"
#include "CellTrace.h"
#include <iostream>
#include "time.h"


#include "vtkRenderer.h" 
#include <QVTKWidget.h>
#include "vtkSmartPointer.h"
#include "vtkSmartVolumeMapper.h"
#include "vtkVolume.h"
#include "vtkVolumeProperty.h"
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>






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
 #include "itkSpatialObjectToImageFilter.h"
#include "itkEllipseSpatialObject.h"
#include <itkImageToVTKImageFilter.h>
#include <itkChainCodePath.h>

 
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

#include "TraceBit.h"
#include "TraceGap.h"
#include "TraceLine.h"
#include "TraceObject.h"
#include "TraceModel.h"
#include "ImageFeatureThreshold.h"
#include "itkStatisticsImageFilter.h"

#include <itkPathToChainCodePathFilter.h>
#include "itkPolyLineParametricPath.h"


#include <qapplication.h>
#include <qdir.h>
#include <qevent.h>
#include <qfileinfo.h>
#include <qimage.h>





#if QT_VERSION < 0x040000
#define modifiers state
#define AltModifier AltButton
#define setWindowTitle setCaption
#endif

#if QT_VERSION < 0x030000
#define flush flushX
#endif


class ImageViewer : public QMainWindow
 {
     Q_OBJECT

	 typedef float PixelType;
	 typedef itk::Image< PixelType, 3 >  ImageType3D;
	 typedef itk::Image< unsigned int, 3 > LabelImageType3D;
	 typedef itk::Image< unsigned char, 3 >   UnsignedCharImageType;
	 typedef itk::StatisticsImageFilter<ImageType3D> StatisticsImageFilterType;


	 typedef itk::RegionOfInterestImageFilter< UnsignedCharImageType, UnsignedCharImageType> RegionOfInterestFilter;
	 typedef itk::ImageToVTKImageFilter<UnsignedCharImageType> ConnectorType;

 public:

	ImageViewer(QWidget *parent = 0);
	~ImageViewer();
	void setImageFeatureThreshold(ImageFeatureThreshold* featureThreshold);
	void showImage();
	void showTraces();
	void showTracesXML();
	void UpdateLineActor();
	void UpdateBranchActor();
	void createSphere();
	void Create3DSyntheticImagesToy();
	void saveScreenShot(std::string);
	void saveTraceFeatures(std::string);
	void calculatePathLength();
	QStringList findFiles(const QString&);
	
 signals:
	void acceptImage(ImageFeatureThreshold* featureThreshold);

public slots:
	void saveImage();
	void closeWindow();
	void RayCastOpacityChanged(int value);     
	void RayCastOpacityValueChanged(double value);
	void RayCastBrightnessChanged(int value);
	void RayCastColorValueChanged(int value);

	 
	
 protected:
 
 private:
	QVTKWidget *QVTK;
	QAction *accept;
	QAction *reject;
	QToolBar *EditsToolBar;
	QWidget *CentralWidget;
	QLabel *costThresholdLabel, *contrastThresholdLabel, *debrisThresholdLabel, *deviceLabel,*offShootLabel,*intensityThresholdLabel,*noOfTraceLinesLabel;
	
	
	
	
	QSpinBox * OpacitySpin, * BrightnessSpin, * SomaOpacitySpin, * SomaBrightnessSpin;
	QDoubleSpinBox * OpacityValueSpin, * SomaOpacityValueSpin, * SomaColorSpin;
	QSlider *OpacitySlider, *BrightnessSlider, *SomaOpacitySlider, *SomaBrightnessSlider;
	QToolBar *BranchToolBar, *RaycastBar, *SomaBar, *SlicerBar;
	QComboBox *ColorProfileCombo;

	vtkSmartPointer<vtkRenderer> Renderer;
	ImageRenderActors *ImageActors;
	ImageFeatureThreshold* _featureThreshold;
	
	vtkSmartPointer<vtkPolyData> poly_line_data;
	vtkSmartPointer<vtkActor> LineActor;
	vtkSmartPointer<vtkPolyDataMapper> LineMapper;

	vtkSmartPointer<vtkPolyData> poly;
	vtkSmartPointer<vtkPolyDataMapper> polymap;
	vtkSmartPointer<vtkActor> BranchActor;
	TraceObject* tobj;
	TraceModel *TreeModel;
	
	
     
 };


 #endif