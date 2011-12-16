#ifndef _IMAGE_3DVIEW_H
#define _IMAGE_3DVIEW_H

// qt includes:
#include <QtGui/QMainWindow>
#include <QObject>
#include <QGridLayout>
#include <QDockWidget>
#include <QtGui>
// vtk includes:
#include "vtkImageViewer2.h"
#include "vtkImageActor.h"
#include "vtkSmartPointer.h"
#include "QVTKWidget.h"
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
#include "vtkCamera.h"
#include "vtkStringArray.h"
#include "vtkLabeledDataMapper.h"
#include "vtkSelectVisiblePoints.h"
#include "vtkActor2D.h"
#include "vtkIntArray.h"
#include "vtkLabelPlacementMapper.h"
#include "vtkPointSetToLabelHierarchy.h"
#include <vtkTextProperty.h>
#include <vtkImagePlaneWidget.h>
#include <vtkCallbackCommand.h>
#include <vtkActor2DCollection.h>
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
#include <vtkProperty2D.h>
#include <vtkGPUVolumeRayCastMapper.h>
#include <vtkVolumeRayCastMapper.h>
#include <vtkVolumeRayCastCompositeFunction.h>
#include <vtkVolumeRayCastMIPFunction.h>
#include <vtkFixedPointVolumeRayCastMapper.h>
#include <vtkVolume.h>
#include <vtkImageImport.h>
#include <vtkVolumeProperty.h>
#include <vtkPiecewiseFunction.h>
#include <vtkColorTransferFunction.h>
#include <vtkContourFilter.h>
#include "vtkRenderWindow.h"
#include <vtkSphereSource.h>
#include <vtkSurfaceReconstructionFilter.h>
#include <vtkReverseSense.h>
#include <vtkTextActor3D.h>
#include <vtkBoxWidget.h>
#include <vtkDataSetMapper.h>
#include <vtkPolygon.h>
#include <vtkCleanPolyData.h>
#include <vtkDelaunay3D.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkMarchingContourFilter.h>
#include <vtkContourGrid.h>
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkMarchingCubes.h"
#include "vtkPolyDataConnectivityFilter.h"
//stl includes
#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <map>

//standard c++ includes
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


// Farsight Includes:
#include "ftkGUI/LabelImageViewQT.h"
#include "ftkGUI/ObjectSelection.h"


//Macros
#define Z_SPACING 1
#define MAX(a,b) (((a) > (b))?(a):(b))
#define MIN(a,b) (((a) < (b))?(a):(b))





struct CellActors{
		vtkSmartPointer<vtkActor> boundaryActor;
		vtkSmartPointer<vtkActor> centroidActor;
};

class Image3DView: public QObject
{
	Q_OBJECT

public:
	// The label image must be passed:
	Image3DView(ftk::Image::Pointer image,ftk::Image::Pointer labimage, LabelImageViewQT * imview = NULL, ObjectSelection * sels = NULL);	// Constructor
	~Image3DView();  // Destructor
	 static void HandleKeyPress(vtkObject* caller, unsigned long event, void* clientdata, void* callerdata);

public slots:
		void RefreshTimeActors(void);
		void RefreshTrackActors(void);
		void ToggleLabelVisibility(void);
		void Toggle3DStackTrackView(void);
		void ToggleContourVisibility(void);
		void ToggleVolumeVisibility(void);
		void ToggleTrackVisibility(void);
		void ChangeOpacity(void);
		void ChangeBrightness(void);

private:
	QMainWindow * MainWindow;
	QVTKWidget * QVTKView;
	// widgets:
	QDockWidget * DockWidget;
	QSlider * opacitySlider;
	QSlider * brightnessSlider;
	QLabel * opacityLabel;
	QLabel * brightnessLabel;

	vtkSmartPointer<vtkRenderer> Renderer;
	vtkSmartPointer<vtkCallbackCommand> keyPress;
	vtkSmartPointer<vtkRenderWindow > RenderWindow;
	vtkSmartPointer<vtkRenderWindowInteractor> Interactor;

	// Functions:
	void SetupRenderWindow(void);
	void Render3DStack(int currentT);
	void CreateActors(void);
	void CreateBoundaryActors(int currentT);
	void CreateCentroidActor(int currentT);
	void CreateLabelActor(int currentT);
	void CreateInteractorStyle(void);
	void CreateVolumes(void);
	void GenerateTracks(void);
	void CreateTracks(void);
	void CreateTrackPoints(void);
	void CreateBoundingBox(void);
	void GetBoudaryActorColor(int id, unsigned char color[3]);



	// Flags:
	bool labelsVisible;
	bool stacksVisible;
	bool volumesVisible;
	bool countoursVisible;
	bool tracksVisible;
	// Data:	
	LabelImageViewQT * ImageView;
	ObjectSelection * Selection;
	ftk::Image::Pointer LabelImageData;
	ftk::Image::Pointer ImageData;
	ftk::Image::PtrMode mode;
	TraceObject * TraceData;
	vtkSmartPointer<vtkPolyData> TrackPoly;
	vtkSmartPointer<vtkActor> bitsActor;
	vtkSmartPointer<vtkActor> tracksActor;
	vtkSmartPointer<vtkActor> boxActor;

	std::vector<std::map<int, ftk::Object::Point> > CenterMap;
	std::vector<std::map<int, ftk::Object::Box> > bBoxMap;
	std::vector<std::map<int, CellActors> > CellActorsMap;
	std::vector<vtkSmartPointer<vtkActor2D> > ImageLabelsVector;
	std::vector<vtkSmartPointer<vtkVolume> > ImageVolumes;





};
#endif
