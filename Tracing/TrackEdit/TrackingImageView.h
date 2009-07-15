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

#ifndef _TRACKING_IMAGE_VIEW_H
#define _TRACKING_IMAGE_VIEW_H

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

#include "QGridLayout.h"
#include "QDockWidget.h"
#include <QtGui>
#include <QObject>




class TrackingImageView: public QObject
{
	Q_OBJECT
public:
	TrackingImageView(TrackingModel *model)
	{
		DEBUG3("Entering TrackingImageView constructor\n");
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

		GenerateImages();
		
		
		m_vtkimageshiftscale = vtkSmartPointer<vtkImageShiftScale>::New();
		m_vtkimageshiftscale->SetInput(showimages[0]);
		m_vtkimageshiftscale->SetScale(3);
		m_vtkimageshiftscale->SetOutputScalarTypeToUnsignedChar();
		m_vtkimageshiftscale->ClampOverflowOn();

		m_vtkimageactor = vtkSmartPointer<vtkImageActor>::New();
		m_vtkimageactor->SetInput(m_vtkimageshiftscale->GetOutput());
		m_vtkcontouractor = showborders[0];
		
		m_vtkrenderer->AddActor(m_vtkimageactor);
		m_vtkrenderer->AddActor(m_vtkcontouractor);

		m_imageview->GetRenderWindow()->AddRenderer(m_vtkrenderer);
		m_imageview->GetRenderWindow()->Render();
		m_timeslider->setMinimum(0);
		m_timeslider->setMaximum(m_model->getNumTimePoints()-1);
		m_timeslider->setSingleStep(1);
		m_timeslider->setPageStep(1);
		m_timeslider->setTickInterval(1);
		m_timeslider->setTickPosition(QSlider::TicksRight);
		
		InitializeCallbacks();
		connect(m_timeslider,SIGNAL(valueChanged(int)),SLOT(setTime(int)));
		connect(m_model,SIGNAL(labelsChanged(int)),SLOT(refreshImages(int)));
		m_mainwindow->show();
		m_timeslider->show();
		m_imageview->show();
		m_imageview->setMinimumSize(m_model->GetSizeX()*1.5,m_model->GetSizeY()*1.5);
		DEBUG3("Exiting constructor\n");
	}
	~TrackingImageView()
	{
		delete m_mainwindow;
		delete m_imageview;
		delete statusLabel;
		delete m_timeslider;
		delete dockwidget;

		//I dont know why I'm deleting smartpointers. vtkDebugLeaks keeps complaining! and still :(
		sliderRep1->Delete();
	    sliderWidget1->Delete();
	    m_vtkimageshiftscale->Delete();
	    m_vtkimageactor->Delete();
	    m_vtkcontouractor->Delete();
	    m_vtkrenderer->Delete();
	    m_imcallback->Delete();
	    m_lccallback->Delete();
	    m_kycallback->Delete();
	    m_brcallback->Delete();
	    m_picker->Delete();
	}
	void InitializeCallbacks();
	void InitializeSliders();
	typedef itk::Image<unsigned char,3> InputImageType;
	typedef itk::Image<unsigned short,3> LabelImageType;
	void GenerateImages();

	void SetModel(TrackingModel* model){ m_model = model;}
	void SetSelection(TrackSelection* select)
	{
		m_select = select;
		m_selectionactor = vtkSmartPointer<vtkActor>::New();
		m_vtkrenderer->AddActor(m_selectionactor);
		connect(m_select,SIGNAL(changed(void)),SLOT(refreshSelection(void)));
	}
	TrackingModel* GetModel();
	void Render();
	public slots:
		void setTime(int t);
		void refreshSelection();
		void refreshImages(int t);
private:
	static void changeBrightness(vtkObject * object, unsigned long eid, void * clientdata, void * callerdata);
	static void pickPixel(vtkObject * object, unsigned long eid, void * clientdata, void * callerdata);
	static void keyPress(vtkObject * object, unsigned long eid, void *clientdata, void * callerdata);
	std::vector<vtkSmartPointer<vtkImageData>> showimages;
	std::vector<vtkSmartPointer<vtkActor>> showborders;
	std::vector<std::vector<vtkSmartPointer<vtkTextActor>>> shownumbers;
	vtkSmartPointer<vtkActor> m_selectionactor;

	TrackingModel* m_model;
	TrackSelection* m_select;
	int m_currenttime;
	typedef ftk::LabelImageFeatures FeatureType;
	std::vector<std::vector<FeatureType>> *features;

	QMainWindow * m_mainwindow;
	QVTKWidget * m_imageview;
	QSlider * m_timeslider;
	QLabel * statusLabel;
	QDockWidget * dockwidget;

	vtkSmartPointer<vtkSliderRepresentation2D> sliderRep1;
	vtkSmartPointer<vtkSliderWidget> sliderWidget1;
	vtkSmartPointer<vtkImageShiftScale> m_vtkimageshiftscale;
	vtkSmartPointer<vtkImageActor> m_vtkimageactor;
	vtkSmartPointer<vtkActor> m_vtkcontouractor;
	vtkSmartPointer<vtkRenderer> m_vtkrenderer;
	vtkSmartPointer<vtkCallbackCommand> m_imcallback;
	vtkSmartPointer<vtkCallbackCommand> m_lccallback;
	vtkSmartPointer<vtkCallbackCommand> m_kycallback;
	vtkSmartPointer<vtkCallbackCommand> m_brcallback;
	vtkSmartPointer<vtkCellPicker> m_picker;

};

#endif