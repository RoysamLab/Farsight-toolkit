#ifndef _IMAGE_3DVIEW_H
#define _IMAGE_3DVIEW_H

#include <QtGui/QMainWindow>
#include <QObject>


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




struct CellActors{
		vtkSmartPointer<vtkActor> boundaryActor;
		vtkSmartPointer<vtkActor> centroidActor;
};

class Image3DView: public QObject
{
	Q_OBJECT

public:
	// The label image must be passed:
	Image3DView(ftk::Image::Pointer image,std::vector<std::map<int, ftk::Object::Point> > CenterMapVector,LabelImageViewQT * imview = NULL, ObjectSelection * sels = NULL);	// Constructor
	~Image3DView();  // Destructor
	 static void HandleKeyPress(vtkObject* caller, unsigned long event, void* clientdata, void* callerdata);

public slots:
		void RefreshTimeActors(void);
		void RefreshTrackActors(void);
		void ToggleLabelVisibility(void);
		void Toggle3DStackTrackView(void);

private:
	// Functions:
	void SetupRenderWindow(void);
	void Render3DStack(int currentT);
	void CreateActors(void);
	void CreateBoundaryActors(int currentT);
	void CreateCentroidActor(int currentT);
	void CreateLabelActor(int currentT);
	void CreateInteractorStyle(void);


	// Flags:
	bool labelsVisible;
	bool stacksVisible;
	// Data:	
	LabelImageViewQT * ImageView;
	ObjectSelection * Selection;
	ftk::Image::Pointer LabelImageData;
	ftk::Image::PtrMode mode;
	vtkSmartPointer<vtkImageData> VTKImage;

	std::vector<std::map<int, ftk::Object::Point> > CenterMap;
	std::vector<std::map<int, ftk::Object::Box> > bBoxMap;
	std::vector<std::map<int, CellActors> > CellActorsMap;
	std::vector<vtkSmartPointer<vtkActor2D> > ImageLabelsVector;

	QMainWindow * MainWindow;
	QVTKWidget * QVTKView;
	vtkSmartPointer<vtkRenderer> Renderer;
	vtkSmartPointer<vtkCallbackCommand> keyPress;
	vtkSmartPointer<vtkRenderWindow > RenderWindow;
	vtkSmartPointer<vtkRenderWindowInteractor> Interactor;



};
#endif