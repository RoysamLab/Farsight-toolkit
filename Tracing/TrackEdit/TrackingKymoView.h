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

#ifndef _TRACKING_KYMO_VIEW_H
#define _TRACKING_KYMO_VIEW_H

#include "TrackingModel.h"
#include "TrackSelection.h"
#include "helpers.h"
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
#include "Trace.h"

#include "QGridLayout.h"
#include "QDockWidget.h"
#include <QtGui>
#include <QObject>



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

	  	/*vtkSmartPointer<vtkPiecewiseFunction> colorTransferFunction = vtkSmartPointer<vtkPiecewiseFunction>::New();
		colorTransferFunction->AddPoint(1.0,0.0f,0.0f,0.0f);
		colorTransferFunction->AddPoint(255*(1-value),0.5f,0.5f,0.0f);*/
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
	TrackingKymoView(TrackingModel* model)
	{
		m_model = model;
		features = m_model->getFeatures();
		m_currenttime = 0;
		m_mainwindow = new QMainWindow();
		m_imageview = new QVTKWidget(m_mainwindow);
		m_mainwindow->setCentralWidget(m_imageview);
		m_mainwindow->show();

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

		connect(m_model,SIGNAL(labelsChanged(int)),SLOT(refreshImages(int)));

		m_mainwindow->show();
		m_imageview->show();
		m_imageview->setMinimumSize(m_model->GetSizeX()*1,m_model->GetSizeY()*1);
		m_vtkrenderer->Render();
	}
	void GenerateImages();
	void GenerateTracks();
	~TrackingKymoView()
	{
		delete m_mainwindow;
		delete m_imageview;
		delete statusLabel;

	}
	void SetSelection(TrackSelection* select)
	{
		m_select = select;
		m_selectionactor = vtkSmartPointer<vtkActor>::New();
		m_vtkrenderer->AddActor(m_selectionactor);
		connect(m_select,SIGNAL(changed(void)),SLOT(refreshSelection(void)));
	}
public slots:
		void refreshSelection();
		void refreshImages(int t);
private:

	vtkSmartPointer<vtkActor> getTrackPoints(std::vector<TraceBit> vec);
	void AddSliders();

	TrackingModel * m_model;

	QMainWindow * m_mainwindow;
	QVTKWidget * m_imageview;
	QLabel * statusLabel;

	TrackSelection* m_select;
	TraceObject * m_tobj;
	int m_currenttime;
	typedef ftk::LabelImageFeatures FeatureType;
	std::vector<std::vector<FeatureType>> *features;

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