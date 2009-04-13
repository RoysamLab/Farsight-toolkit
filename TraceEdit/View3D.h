#ifndef VIEW3D_H_
#define VIEW3D_H_
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

#include "vtkLineSource.h"
#include "vtkTubeFilter.h"
#include "vtkPolyDataMapper.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkCubeSource.h"
#include "vtkCallbackCommand.h"
#include <stdio.h>

#include "vtkImageReader2.h"
#include "vtkContourFilter.h"
#include "vtkPolyDataNormals.h"
#include "vtkCleanPolyData.h"
#include "vtkImageData.h"

#include "vtkImageToStructuredPoints.h"
#include "vtkStructuredPointsReader.h"
#include "vtkStructuredPoints.h"
#include "vtkPiecewiseFunction.h"
#include "vtkColorTransferFunction.h"
#include "vtkVolumeProperty.h"
#include "vtkVolumeRayCastMapper.h"
#include "vtkVolumeRayCastCompositeFunction.h"
#include "vtkVolume.h"
#include "vtkOpenGLVolumeTextureMapper3D.h"
#include "vtkOpenGLVolumeTextureMapper2D.h"
#include "vtkSliderWidget.h"
#include "vtkSliderRepresentation2D.h"
#include "Trace.h"
#include "vtkPlaybackWidget.h"
#include "vtkPlaybackRepresentation.h"

struct compTrace{
TraceLine *Trace1;
TraceLine *Trace2;
int endPT1, endPT2;

double dist; 
};

class View3d
{
public:

  View3d();
  ~View3d();
  void Initialize(int argc, char **argv);

  //general render window variables
	vtkSmartPointer<vtkRenderer> ren;
	vtkSmartPointer<vtkRenderWindow> renWin;
	vtkSmartPointer<vtkRenderWindowInteractor> iren;
	vtkSmartPointer<vtkCamera> cam;
	void RenderWin();

	void addAct(vtkActor *Actor);	
	void AddPointsAsPoints (std::vector<TraceBit> vec);
	void AddBranchIllustrators();
	void AddVolumeSliders();
	void AddContourThresholdSliders();
	void AddPlaybackWidget(char*);
	vtkActor *lineAct;
	vtkPolyDataMapper *lineMap;
	vtkActor* LineAct();

//interactor variables and point picking
	void interact();
	static void PickCell(vtkObject* caller, unsigned long event, void* clientdata, void* callerdata);
	static void SetMode(vtkObject* caller, unsigned long event, void* clientdata, void* callerdata);
    void HighlightSelected(TraceLine* tline, int t);
	vtkCallbackCommand* isPicked;
	vtkCallbackCommand* keyPress;
	
	vtkCellPicker *cell_picker;
	std::vector<int> IDList;
	
	void deleteTrace(View3d* view,TraceLine *tline);
	void MinEndPoints(View3d* view,std::vector<TraceLine*> traceList);
	vtkSphereSource *sphere;
	vtkPolyDataMapper *sphereMap;
	vtkActor *sphereAct;
//  img reading	and contour->3d
	void readImg(char* sourceFile);
	vtkSmartPointer<vtkPolyDataMapper> volMap;
	vtkSmartPointer<vtkActor> volAct;
	TraceObject* tobj;
	vtkVolume * vol;
	vtkSmartPointer<vtkContourFilter> cfilt;
//raycast
	void rayCast(char* raySource);
	vtkSmartPointer<vtkPolyData> poly_line_data;
private:


	
	//idk 

};
#endif
