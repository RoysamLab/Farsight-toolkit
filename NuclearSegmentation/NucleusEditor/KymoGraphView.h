#ifndef _TRACKING_KYMO_VIEW_H
#define _TRACKING_KYMO_VIEW_H

#include "Trace.h"
#include "vtkImageViewer2.h"
#include "vtkImageActor.h"
#include "QMainWindow.h"
#include "QSlider.h"
#include "QLabel.h"
#include "vtkSmartPointer.h"
#include "QVTKWidget.h"
#include "QStatusBar.h"
#include "vtkImageData.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkImageMapToColors.h"
#include "vtkCallbackCommand.h"
#include "vtkPointPicker.h"
#include "vtkCellPicker.h"
#include "vtkRegularPolygonSource.h"
#include "vtkCubeSource.h"
#include "vtkInteractorStyleImage.h"
#include "vtkAppendPolyData.h"
#include "vtkImageShiftScale.h"
#include "vtkSliderWidget.h"
#include "vtkSliderRepresentation2D.h"
//#include "Trace.h"
#include "vtkCamera.h"

#include "QGridLayout.h"
#include "QDockWidget.h"
#include <QtGui>
#include <QObject>
#include <QtGui/QDialog>

//stl includes
#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>



//standard c++ includes
#include <stdio.h>
#include <stdlib.h>

//ITK includes
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageSliceConstIteratorWithIndex.h>
#include <itkImageLinearIteratorWithIndex.h>
#include <itkConstNeighborhoodIterator.h>
#include "itkImageToVTKImageFilter.h"
#include "itkLabelGeometryImageFilter.h"
#include "ftkLabelImageToFeatures.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkScalarImageTextureCalculator.h"

//VTK includes
#include "vtkPolyData.h"
#include "vtkContourFilter.h"
#include "vtkActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataNormals.h"
#include "vtkProperty.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkTextActor.h"
#include "vtkTextProperty.h"
#include "vtkCellArray.h"
#include "vtkPoints.h"
#include <vtkOpenGLVolumeTextureMapper2D.h>
#include <vtkOpenGLVolumeTextureMapper3D.h>
#include <vtkVolumeRayCastMapper.h>
#include <vtkVolumeRayCastCompositeFunction.h>
#include <vtkVolumeRayCastMIPFunction.h>
#include <vtkFixedPointVolumeRayCastMapper.h>
#include <vtkVolume.h>
#include <vtkImageImport.h>
#include <vtkVolumeProperty.h>
#include <vtkPiecewiseFunction.h>
#include <vtkColorTransferFunction.h>

//Macros
#define MAX(a,b) (((a) > (b))?(a):(b))
#define MIN(a,b) (((a) < (b))?(a):(b))
#define CLAMP(a,min,max) (((a)<(min))?(min):(((a)>(max))?(max):(a)))
#define DEBUG1 
#define DEBUG2 printf
#define DEBUG3
#define PROGRESS printf


typedef unsigned char InputPixelType;
typedef unsigned char OutputPixelType;

typedef itk::Image<InputPixelType,3> InputImageType;
typedef itk::Image<InputPixelType,2> Input2DImageType;

typedef itk::Image<short,3> LabelImageType;
typedef itk::Image<short,2> Label2DImageType;

typedef itk::Vector<unsigned char, 3> VectorPixelType;
typedef itk::Image<VectorPixelType, 3> ColorImageType;
typedef itk::Image<VectorPixelType, 2> Color2DImageType;

typedef itk::ImageRegionConstIterator<LabelImageType> ConstLabelIteratorType;
typedef itk::ImageRegionIterator<LabelImageType> LabelIteratorType;
typedef itk::ImageRegionIterator<InputImageType> IteratorType;


typedef itk::ImageLinearIteratorWithIndex< Color2DImageType > LinearColorIteratorType;
typedef itk::ImageSliceConstIteratorWithIndex< ColorImageType > SliceColorIteratorType;

typedef itk::ImageLinearIteratorWithIndex< Input2DImageType > LinearIteratorType;
typedef itk::ImageSliceConstIteratorWithIndex< InputImageType > SliceIteratorType;;

typedef itk::ImageRegionConstIterator<Input2DImageType> Const2DIteratorType;
typedef itk::ImageRegionIterator<Input2DImageType> twoDIteratorType;

typedef itk::ImageToVTKImageFilter<InputImageType> ConnectorType;
typedef itk::ImageToVTKImageFilter<Input2DImageType> Connector2DType;

class vtkSlider2DKymoCallbackBrightness : public vtkCommand
{
public:
  static vtkSlider2DKymoCallbackBrightness *New() 
    { return new vtkSlider2DKymoCallbackBrightness; }
  virtual void Execute(vtkObject *caller, unsigned long, void*)
    {
      vtkSliderWidget *sliderWidget = 
        reinterpret_cast<vtkSliderWidget*>(caller);
      double value = static_cast<vtkSliderRepresentation *>(sliderWidget->GetRepresentation())->GetValue();

	  	vtkSmartPointer<vtkPiecewiseFunction> opacityTransferFunction = vtkSmartPointer<vtkPiecewiseFunction>::New();
		opacityTransferFunction->AddPoint(2,0.0);
		//opacityTransferFunction->AddPoint(2+256*(1-value),0.2);
		opacityTransferFunction->AddPoint(50,value);
		volume->GetProperty()->SetScalarOpacity(opacityTransferFunction);
    }
  vtkSlider2DKymoCallbackBrightness() {

  }

vtkVolume* volume;
};

class vtkSlider2DKymoCallbackContrast : public vtkCommand
{
public:
  static vtkSlider2DKymoCallbackContrast *New() 
    { return new vtkSlider2DKymoCallbackContrast; }
  virtual void Execute(vtkObject *caller, unsigned long, void*)
    {
      vtkSliderWidget *sliderWidget = 
        reinterpret_cast<vtkSliderWidget*>(caller);
      double value = static_cast<vtkSliderRepresentation *>(sliderWidget->GetRepresentation())->GetValue();

		vtkColorTransferFunction *colorTransferFunction = vtkColorTransferFunction::New();
		colorTransferFunction->AddRGBPoint(0.0, 0.0, 0.0, 0.0);
		colorTransferFunction->AddRGBPoint(255*(1-value),1,1,1);
		volume->GetProperty()->SetColor(colorTransferFunction);
    }
  vtkSlider2DKymoCallbackContrast() {

  }

  vtkVolume* volume;
};


class TrackingKymoView: public QObject
{
	Q_OBJECT
public:
	TrackingKymoView(ftk::Image::Pointer myImg,std::vector<std::vector<ftk::IntrinsicFeatures>> vvfeatures )
	{
		myfeatures = vvfeatures;
		my4DImg = myImg;
		m_currenttime = 0;
		m_mainwindow = new QMainWindow();
		m_mainwindow->setWindowTitle(tr("Tracking: KymoGraph"));
		m_imageview = new QVTKWidget(m_mainwindow);
		m_mainwindow->setCentralWidget(m_imageview);

		statusLabel = new QLabel(QObject::tr(" Ready"));
		m_mainwindow->statusBar()->addWidget(statusLabel,1);
		
		m_vtkrenderer = vtkSmartPointer<vtkRenderer>::New();
		m_vtkrenderer->BackingStoreOn();
		//m_imageview->setAutomaticImageCacheEnabled(true);
		m_imageview->GetRenderWindow()->AddRenderer(m_vtkrenderer);
		m_imageview->GetRenderWindow()->Render();
		GenerateImages();
		AddSliders();
		GenerateTracks();
		m_vtkrenderer->AddActor(getTrackPoints(m_tobj->CollectTraceBits()));
		m_trackpoly = m_tobj->GetVTKPolyData();
		m_trackmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		m_trackmapper->SetInput(m_trackpoly);
		m_trackactor = vtkSmartPointer<vtkActor>::New();
		m_trackactor->SetMapper(m_trackmapper);
		m_trackactor->GetProperty()->SetLineWidth(3.0);
		m_vtkrenderer->AddActor(m_trackactor);

//		connect(m_model,SIGNAL(labelsChanged(int)),SLOT(refreshImages(int)));

		m_mainwindow->show();
		m_imageview->show();
		m_imageview->setMinimumSize(my4DImg->GetImageInfo()->numColumns*1,my4DImg->GetImageInfo()->numRows*1);
		m_vtkrenderer->Render();
		SaveAnimation();
	}
	void GenerateImages();
	void GenerateTracks();
	~TrackingKymoView()
	{
		delete m_mainwindow;
		delete m_imageview;
		delete statusLabel;

	}
/*	void SetSelection(TrackSelection* select)
	{
		m_select = select;
		m_selectionactor = vtkSmartPointer<vtkActor>::New();
		m_vtkrenderer->AddActor(m_selectionactor);
		connect(m_select,SIGNAL(changed(void)),SLOT(refreshSelection(void)));
	}*/
public slots:
		void refreshSelection();
//		void refreshImages(int t);
		void refreshImages(void);
private:

	vtkSmartPointer<vtkActor> getTrackPoints(std::vector<TraceBit> vec);
	void SaveAnimation();
	void AddSliders();
//    Input2DImageType::Pointer getProjection(InputImageType::Pointer im);
	vtkSmartPointer<vtkVolume> getOneVTKVolume(vtkSmartPointer<vtkImageData> vtkim, float colors[3]);

	QMainWindow * m_mainwindow;
	QVTKWidget * m_imageview;
	QLabel * statusLabel;

//	TrackSelection* m_select;
	TraceObject * m_tobj;
	int m_currenttime;

	ftk::Image::Pointer my4DImg;
	std::vector<std::vector<ftk::IntrinsicFeatures>> myfeatures;

	vtkSmartPointer<vtkRenderer> m_vtkrenderer;
	vtkSmartPointer<vtkActor> m_selectionactor;
	vtkSmartPointer<vtkImageData> m_vtkim;
	vtkSmartPointer<vtkVolume> m_vtkvolume;
	vtkSmartPointer<vtkActor> m_trackactor;
	vtkSmartPointer<vtkPolyDataMapper> m_trackmapper;
	vtkSmartPointer<vtkPolyData> m_trackpoly;

	vtkSlider2DKymoCallbackContrast *m_callback_contrast;
	vtkSlider2DKymoCallbackBrightness *m_callback_brightness;

};

#endif