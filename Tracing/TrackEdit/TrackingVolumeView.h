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

#ifndef _TRACKING_VOLUME_VIEW_H
#define _TRACKING_VOLUME_VIEW_H

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

//class QVTKWidgetFarsight: public QVTKWidget
//{
//public:
//	QVTKWidgetFarsight(QWidget *parent):QVTKWidget(parent)
//	{
//	}
//protected:
//	virtual void paintEvent(QPaintEvent * event)
//	{
//		if (!this->cachedImageCleanFlag)
//		{
//			vtkRenderWindowInteractor* iren = NULL;
//			if(this->mRenWin)
//				iren = this->mRenWin->GetInteractor();
//
//			if(!iren || !iren->GetEnabled())
//				return;
//
//			iren->Render();
//		}
//
//	//	QPainter painter(this);
//	//	painter.drawPixmap(0,0,this->cachedImage);
//		
//		// if we have a saved image, use it
//  if (this->cachedImageCleanFlag)
//    {
//    vtkUnsignedCharArray* array = vtkUnsignedCharArray::SafeDownCast(
//      this->mCachedImage->GetPointData()->GetScalars());
//    // put cached image into back buffer if we can
//    this->mRenWin->SetPixelData(0, 0, this->width()-1, this->height()-1,
//                                array, !this->mRenWin->GetDoubleBuffer());
//    // swap buffers, if double buffering
//    this->mRenWin->Frame();
//    // or should we just put it on the front buffer?
//    return;
//    }
//
//
//	}
//};


class vtkSlider2DCallbackBrightness : public vtkCommand
{
public:
  static vtkSlider2DCallbackBrightness *New() 
    { return new vtkSlider2DCallbackBrightness; }
  virtual void Execute(vtkObject *caller, unsigned long, void*)
    {
      vtkSliderWidget *sliderWidget = 
        reinterpret_cast<vtkSliderWidget*>(caller);
      double value = static_cast<vtkSliderRepresentation *>(sliderWidget->GetRepresentation())->GetValue();

	  	vtkSmartPointer<vtkPiecewiseFunction> opacityTransferFunction = vtkSmartPointer<vtkPiecewiseFunction>::New();
		opacityTransferFunction->AddPoint(2,0.0);
		//opacityTransferFunction->AddPoint(2+256*(1-value),0.2);
		opacityTransferFunction->AddPoint(50,value);
		for(int counter =0; counter < volume->size(); counter++)
		{
     (*(this->volume))[counter]->GetProperty()->SetScalarOpacity(opacityTransferFunction);
		}
    }
  vtkSlider2DCallbackBrightness() {

  }

std::vector<vtkSmartPointer<vtkVolume>>* volume;
};

class vtkSlider2DCallbackContrast : public vtkCommand
{
public:
  static vtkSlider2DCallbackContrast *New() 
    { return new vtkSlider2DCallbackContrast; }
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
		for(int counter=0; counter< volume->size(); counter++)
		{
			(*(this->volume))[counter]->GetProperty()->SetColor(colorTransferFunction);
		}
    }
  vtkSlider2DCallbackContrast() {

  }

  std::vector<vtkSmartPointer<vtkVolume>>* volume;
};

class TrackingVolumeView: public QObject
{
	Q_OBJECT
public:
	TrackingVolumeView(TrackingModel* model)
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
		
		dockwidget = new QDockWidget(m_mainwindow);
		dockwidget->setAllowedAreas(Qt::RightDockWidgetArea\
			|Qt::LeftDockWidgetArea|Qt::TopDockWidgetArea|Qt::BottomDockWidgetArea);
		m_timeslider = new QSlider(dockwidget);
		dockwidget->setWidget(m_timeslider);
		m_mainwindow->addDockWidget(Qt::RightDockWidgetArea,dockwidget);
		m_vtkrenderer = vtkSmartPointer<vtkRenderer>::New();
		m_vtkrenderer->BackingStoreOn();
		//m_imageview->setAutomaticImageCacheEnabled(true);
		m_imageview->GetRenderWindow()->AddRenderer(m_vtkrenderer);
		m_imageview->GetRenderWindow()->Render();

		m_timeslider->setMinimum(0);
		m_timeslider->setMaximum(m_model->getNumTimePoints()-1);
		m_timeslider->setSingleStep(1);
		m_timeslider->setPageStep(1);
		m_timeslider->setTickInterval(1);
		m_timeslider->setTickPosition(QSlider::TicksRight);

		AddSliders();

		GenerateImages();
		
		GenerateTracks();

		m_vtkrenderer->AddActor(getTrackPoints(m_tobj->CollectTraceBits()));

		

		m_trackpoly = m_tobj->GetVTKPolyData();
		m_trackmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		m_trackmapper->SetInput(m_trackpoly);
		m_trackactor = vtkSmartPointer<vtkActor>::New();
		m_trackactor->SetMapper(m_trackmapper);
		m_trackactor->GetProperty()->SetColor(1,0,0);
		m_vtkrenderer->AddActor(m_trackactor);

		
		connect(m_timeslider,SIGNAL(valueChanged(int)),SLOT(setTime(int)));
		connect(m_model,SIGNAL(labelsChanged(int)),SLOT(refreshImages(int)));

		m_mainwindow->show();
		m_timeslider->show();
		m_imageview->show();
		m_imageview->setMinimumSize(m_model->GetSizeX()*1,m_model->GetSizeY()*1);
	}
	void GenerateImages();
	void GenerateTracks();
	~TrackingVolumeView()
	{
		delete m_mainwindow;
		delete m_imageview;
		delete statusLabel;
		delete m_timeslider;
		delete dockwidget;
	}
	void SetSelection(TrackSelection* select)
	{
		m_select = select;
		m_selectionactor = vtkSmartPointer<vtkActor>::New();
		m_vtkrenderer->AddActor(m_selectionactor);
		connect(m_select,SIGNAL(changed(void)),SLOT(refreshSelection(void)));
	}
public slots:
		void setTime(int t);
		void refreshSelection();
		void refreshImages(int t);
private:

	vtkSmartPointer<vtkActor> TrackingVolumeView::getTrackPoints(std::vector<TraceBit> vec);
	void TrackingVolumeView::AddSliders();

	TrackingModel * m_model;

	QMainWindow * m_mainwindow;
	QVTKWidget * m_imageview;
	QSlider * m_timeslider;
	QLabel * statusLabel;
	QDockWidget * dockwidget;

	TrackSelection* m_select;
	TraceObject * m_tobj;
	int m_currenttime;
	typedef ftk::LabelImageFeatures FeatureType;
	std::vector<std::vector<FeatureType>> *features;

	vtkSmartPointer<vtkRenderer> m_vtkrenderer;
	vtkSmartPointer<vtkActor> m_selectionactor;
	std::vector<ConnectorType::Pointer> m_volconnects;
	std::vector<vtkSmartPointer<vtkImageData>> m_vtkims;
	std::vector<vtkSmartPointer<vtkVolume>> m_vtkvolumes;
	vtkSmartPointer<vtkVolume> m_currentvol;
	vtkSmartPointer<vtkActor> m_trackactor;
	vtkSmartPointer<vtkPolyDataMapper> m_trackmapper;
	vtkSmartPointer<vtkPolyData> m_trackpoly;

	vtkSlider2DCallbackContrast *m_callback_contrast;
	vtkSlider2DCallbackBrightness *m_callback_brightness;

};

#endif