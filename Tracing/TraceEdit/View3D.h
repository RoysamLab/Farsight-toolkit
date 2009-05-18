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

#include <QObject>
#include <QtGui>
#include <QVTKWidget.h>

struct compTrace{
TraceLine *Trace1;
TraceLine *Trace2;
int endPT1, endPT2;

double dist; 
};

class View3D : public QWidget 
{
Q_OBJECT;
public:
  View3D(int argc, char **argv);
  ~View3D();
  void Initialize();
  void CreateGUIObjects();
  void CreateLayout();
  void CreateInteractorStyle();
  void CreateActors();
  void UpdateLineActor();
	void UpdateBranchActor();
  void CreateSphereActor();
  bool setTol();
	void AddPointsAsPoints (std::vector<TraceBit> vec);
	void AddVolumeSliders();
	void AddContourThresholdSliders();
	void AddPlaybackWidget(char*);
	static void PickCell(vtkObject* caller, unsigned long event, void* clientdata, void* callerdata);
	static void HandleKeyPress(vtkObject* caller, unsigned long event, void* clientdata, void* callerdata);
  void HighlightSelected(TraceLine* tline, double SelectColor);
	void DeleteTrace(TraceLine *tline);
	void MinEndPoints(std::vector<TraceLine*> traceList);
	void readImg(char* sourceFile);
	void rayCast(char* raySource);

  //todo: make these private with accessors
	vtkSmartPointer<vtkRenderer> Renderer;
	vtkSmartPointer<vtkActor> BranchActor;

public slots:
  void ListSelections();
  void ClearSelection();
  void DeleteTraces();
  void MergeTraces();
  void SplitTraces();
  void FlipTraces();
  void WriteToSWCFile();
  void ShowSettingsWindow();
  void HideSettingsWindow();
  void ApplyNewSettings();

protected:
  void closeEvent(QCloseEvent *event);
  void Rerender();

private:
	int gapTol;
	int gapMax;
	float smallLine;
	float lineWidth;
	double SelectColor;

  //VTK render window embedded in a Qt widget
  QVTKWidget *QVTK;

  //Qt widgets on the main window
  QPushButton *ListButton;
  QPushButton *ClearButton;
  QPushButton *DeleteButton;
  QPushButton *MergeButton;
  QPushButton *SplitButton;
  QPushButton *FlipButton;
  QPushButton *WriteButton;
  QPushButton *SettingsButton;

  //Qt widgets for the settings window
  QWidget *SettingsWidget;
  QLineEdit *MaxGapField;
  QLineEdit *GapToleranceField;
  QLineEdit *LineLengthField;
  QLineEdit *ColorValueField;
  QLineEdit *LineWidthField;
  QPushButton *ApplySettingsButton;
  QPushButton *CancelSettingsButton;

	//stuff for tol and selection
  //general render window variables
	vtkSmartPointer<vtkRenderWindowInteractor> Interactor;
	vtkSmartPointer<vtkActor> LineActor;
	vtkSmartPointer<vtkPolyDataMapper> LineMapper;

  //interactor variables and point picking
	vtkSmartPointer<vtkCallbackCommand> isPicked;
	vtkSmartPointer<vtkCallbackCommand> keyPress;
	
	vtkSmartPointer<vtkCellPicker> CellPicker;
	std::vector<int> IDList;
	vtkSmartPointer<vtkSphereSource> Sphere;
	vtkSmartPointer<vtkPolyDataMapper> SphereMapper;
	vtkSmartPointer<vtkActor> SphereActor;
  //  img reading	and contour->3d
	vtkSmartPointer<vtkPolyDataMapper> VolumeMapper;
	vtkSmartPointer<vtkActor> VolumeActor;
	TraceObject* tobj;
	vtkSmartPointer<vtkVolume> Volume;
	vtkSmartPointer<vtkContourFilter> ContourFilter;
  //raycast
	vtkSmartPointer<vtkPolyData> poly_line_data;
	vtkSmartPointer<vtkPolyData> poly;
	vtkSmartPointer<vtkPolyDataMapper> polymap;
};
#endif
