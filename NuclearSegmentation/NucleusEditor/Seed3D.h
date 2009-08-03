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

#ifndef SEED_EDITOR_H
#define SEED_EDITOR_H

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

//QT Includes:
#include <QtGui/QMainWindow>
#include <QtGui/QMenuBar>
#include <QtGui/QMenu>
#include <QtGui/QMessageBox>
#include <QtGui/QFileDialog>
#include <QtCore/QFileInfo>
#include "itkObjectFactory.h"
#include "itkSmartPointer.h"
#include "itkMacro.h"
#include "itkLightObject.h"
#include "itkImageFileReader.h"
#include "itkImageToVTKImageFilter.h"
#include "itkImage.h"

#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkCamera.h"
#include "vtkProperty.h"

#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkPointPicker.h"
#include "vtkCellPicker.h"
#include "vtkCommand.h"
#include "vtkRendererCollection.h"

#include "vtkPolyDataMapper.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkConeSource.h"
#include "vtkCallbackCommand.h"
#include <stdio.h>

#include "vtkPolyDataNormals.h"
#include "vtkCleanPolyData.h"
#include "vtkImageData.h"
#include "vtkLODActor.h"

#include "vtkImageToStructuredPoints.h"
#include "vtkStructuredPointsReader.h"
#include "vtkStructuredPoints.h"
#include "vtkPiecewiseFunction.h"
#include "vtkColorTransferFunction.h"
#include "vtkVolumeProperty.h"
#include "vtkVolume.h"
#include "vtkOpenGLVolumeTextureMapper3D.h"
#include "vtkSliderWidget.h"
#include "vtkSliderRepresentation2D.h"
#include "vtkPlaybackWidget.h"
#include "vtkPlaybackRepresentation.h"
#include "vtkGlyph3D.h"
#include "vtkFloatArray.h"
#include "vtkCellPicker.h"
#include "vtkCallbackCommand.h"
#include "vtkPointData.h"
#include "vtkSphereWidget.h"
#include "vtkImageReslice.h"
#include "vtkVolume.h"
#include "vtkDataSetMapper.h"
#include "math.h"
#include "vtkCoordinate.h"
#include "vtkImplicitPlaneRepresentation.h"
#include "vtkImplicitPlaneWidget2.h"
#include "vtkPlane.h"
#include "vtkCutter.h"
#include "vtkImageToPolyDataFilter.h"
#include "vtkAppendPolyData.h"
#include "vtkGeometryFilter.h"
#include "vtkSeedRepresentation.h"
#include "vtkOutlineFilter.h"
#include "vtkSeedWidget.h"
#include "vtkImageReslice.h"
#include "vtkLookupTable.h"
#include "vtkImageMapToColors.h"
#include "vtkImageActor.h"
#include "vtkMatrix4x4.h"
#include "vtkImageMapper.h"
#include "vtkInteractorStyleImage.h"
#include "vtkPointHandleRepresentation3D.h"
#include "vtkSphereRepresentation.h"
#include "vtkHandleWidget.h"
#include "vtkInteractorStyleRubberBand2D.h"
#include "vtkPointHandleRepresentation2D.h"
#include "vtkPointHandleRepresentation3D.h"
#include "vtkPlane.h"
#include "vtkProp3DCollection.h"
#include <QObject>
#include <QtGui>
#include <QVTKWidget.h>
#include <iostream>
#include <fstream>
#include "NuclearSegmentation/ftkNuclearSegmentation.h"


class Seed3D : public QMainWindow
{
    Q_OBJECT;
public:
	Seed3D(QWidget * parent = 0, Qt::WindowFlags flags = 0);
	void GetImage(ftk::Image::Pointer data,vector<Seed> );
	private slots:
	//void loadImage(void);
	void PlaceSeed();
    void Apply();
    void AddSeed(); 
	void Check();
	void DeleteSeed();
	//void saveResult();
	void UndoDeleteSeeds();
    void SplitSeeds();
	void MergeSeeds();
	
	void DeleteObjects();
	
private:
	/*void createMenus();
	void createStatusBar();*/

	
	
	
	
	vtkSmartPointer<vtkRenderWindowInteractor> Interactor;
	vtkSmartPointer<vtkRenderWindowInteractor> Interactor1;
    vtkSmartPointer<vtkRenderWindowInteractor> Interactor2;
	vtkSliderRepresentation2D *sliderRep;
	vtkSliderRepresentation2D *sliderRep2;
	vtkSliderRepresentation2D *sliderRep3;
	vtkSliderWidget *sliderWidget;
	vtkSliderWidget *sliderWidget2;
	vtkSliderWidget *sliderWidget3;
    vtkSmartPointer<vtkCellPicker> CellPicker;
	vtkSmartPointer<vtkPointPicker> PointPicker;
	vtkSmartPointer<vtkCallbackCommand> isPicked;
	vtkSmartPointer<vtkCallbackCommand> keyPress;

    
	QMenu *fileMenu;
	QAction *loadAction;
	QAction *saveAction;
	QAction *exitAction;
	QLabel *statusLabel;
	QString lastPath;
	QString fileName;
	QString fileNameSeed;
	QWidget *browse;
    QVTKWidget *QVTK;
    QVTKWidget *QVTK1;	
    QVTKWidget *QVTK2;
	QMessageBox* msgBox;
  
	//Qt widgets on the main window
    QCheckBox *AddBox;
    QCheckBox *DeleteBox;
    QCheckBox *UndoDelBox;	
	QCheckBox *SplitBox;	
	QCheckBox *MergeBox;	
    QPushButton *PlaceButton;
    QPushButton *ApplyButton;
    vtkSmartPointer<vtkRenderer> Renderer;
    vtkSmartPointer<vtkRenderer> Renderer1;
    vtkSmartPointer<vtkRenderer> Renderer2;
    vtkSphereWidget *sphereWidget;
    vtkSmartPointer<vtkGlyph3D> Glyph;
    vtkSmartPointer<vtkGlyph3D> delglyph;
    vtkSmartPointer<vtkGlyph3D> addglyph;



	int mode;
	int stateAdd;
	int stateDelete;
	int stateMerge;
	int stateSplit;
	int stateUndoDel;
	int counter;
	int flag;
	int iRender;
	typedef ftk::Object::Point point;
	
	std::vector<point> dup_points;
	std::vector<point> MarkedPoints;
	std::vector<point> MarkedPoints2add;
	vtkFloatArray* pcoords;
	vtkFloatArray* Addpcoords;
	vtkFloatArray* delpcoords;
	vtkPoints* point1;
	vtkPoints* point2;
	vtkPoints* point3;
	vtkPoints* allpoints;
	vtkSmartPointer<vtkPolyData> polydata1;
	vtkSmartPointer<vtkPolyData> polydata2;
	vtkSmartPointer<vtkPolyData> polydata3;
	vtkSmartPointer<vtkSphereSource> Sphere;
	vtkSmartPointer<vtkPolyDataMapper> SphereMapper;
	vtkSmartPointer<vtkActor> SphereActor;
	vtkSmartPointer<vtkPolyDataMapper> DelSphereMapper;
	vtkSmartPointer<vtkActor> DelSphereActor;
	vtkSmartPointer<vtkPolyDataMapper> AddSphereMapper;
	vtkSmartPointer<vtkActor> AddSphereActor;
    vtkSmartPointer<vtkImageData> VTKim;
	vtkSmartPointer<vtkImageReslice> Reslice;
	vtkSmartPointer<vtkImageReslice> Reslice1;
	vtkSmartPointer<vtkDataSetMapper>im_mapper;
	vtkSmartPointer<vtkDataSetMapper>im_mapper1;
	vtkSmartPointer<vtkActor> imActor;
	vtkSmartPointer<vtkActor> imActor1;

	vtkSmartPointer<vtkHandleWidget> widget;
	vtkSmartPointer<vtkPointHandleRepresentation3D> handle;
	vtkSmartPointer<vtkHandleWidget> widget1;
	vtkSmartPointer<vtkPointHandleRepresentation2D> handle1;
	vtkSmartPointer<vtkHandleWidget> widget2;
	vtkSmartPointer<vtkPointHandleRepresentation2D> handle2;
	vtkSmartPointer<vtkAppendPolyData> apd;
	vtkSmartPointer<vtkPolyDataMapper> VolumeMapper;
	vtkSmartPointer<vtkActor> VolumeActor;
	vtkSmartPointer<vtkVolume> Volume;
	vtkSmartPointer<vtkPolyData> poly;
	vtkSmartPointer<vtkPolyDataMapper> polymap;

    static void PickCell(vtkObject* caller, unsigned long event, void* clientdata, void* callerdata);
    void rayCast(char*,char*);
	
    vector<point> GetSeedpts(vector<Seed> seeds);

	std::vector<point> spPoint;
	
	double* y;
	double* wp;
 };

#endif
