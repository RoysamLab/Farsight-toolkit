#ifndef Seed3D_H_
#define Seed3D_H_
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
#include "vtkPropPicker.h"
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
#include "vtkImageViewer2.h"
#include "vtkInteractorStyleImage.h"
#include "vtkRegressionTestImage.h"
//#include "ftkNuclearSegmentation.h"
#include "vtkPointHandleRepresentation3D.h"
#include "vtkSphereRepresentation.h"
#include "vtkHandleWidget.h"
#include "vtkInteractorStyleRubberBand2D.h"
#include "vtkPointHandleRepresentation2D.h"
#include "vtkPointHandleRepresentation3D.h"
#include "vtkBoundedPlanePointPlacer.h"
#include "vtkPlane.h"
#include "vtkProp3DCollection.h"
#include <QObject>
#include <QtGui>
#include <QVTKWidget.h>	

/*struct compTrace{
TraceLine *Trace1;
TraceLine *Trace2;
int endPT1, endPT2;

double dist; 
};*/

typedef struct
{
float x;
float y;
float z;
} point;




class Seed3D : public QWidget 
{
Q_OBJECT;
public:

  Seed3D(int argc, char **argv);
  ~Seed3D();
  void Initialize();
  void CreateGUIObjects();
  void CreateLayout();
  void CreateInteractorStyle();


  bool setTol();
	void AddVolumeSliders();
	static void PickCell(vtkObject* caller, unsigned long event, void* clientdata, void* callerdata);
        void rayCast(char* raySource,char* fileSource );
        std::vector<point> ReadPoints(char* pointzource);
	void CreateViewer(char* raySource);
        void loadSegmentation(char* segSource);
        void Check();
//        void FindDistance(double x[3]);

  //todo: make these private with accessors
	vtkSmartPointer<vtkRenderer> Renderer;
	vtkSmartPointer<vtkRenderer> Renderer1;
	vtkSmartPointer<vtkRenderer> Renderer2;
	

vtkSmartPointer<vtkActor> BranchActor;
std::vector<point> p;
std::vector<point> dup_points;
std::vector<point> MarkedPoints;
std::vector<point> MarkedPoints2add;
vtkFloatArray* pcoords;
vtkFloatArray* Addpcoords;
vtkFloatArray* Splitpcoords;
vtkFloatArray* delpcoords;


vtkPoints* point1;
vtkPoints* point2;
vtkPoints* point3;
vtkPoints* allpoints;
vtkPolyData* polydata1;
vtkPolyData* polydata2;
vtkPolyData* polydata3;
vtkSphereSource *sphere1; 
vtkGlyph3D *glyph1;



vtkSphereWidget *sphereWidget;
vtkLODActor *selectActor;
vtkLODActor *selectActor1;
vtkActor *actor;
vtkSmartPointer<vtkImageViewer2> ImageViewer;
int mode;
int stateAdd;
int stateDelete;
int stateMerge;
int stateSplit;
int stateUndoDel;
int counter;
int flag;

/*public slots:
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
*/

public slots:

  void PlaceSeed();
  void Apply();
  void AddSeed(); 
  void DeleteSeed(); 
  void SplitSeeds();
  void MergeSeeds();
  void UndoDeleteSeeds();

private:
	int gapTol;
	int gapMax;
	float smallLine;
	float lineWidth;
	double SelectColor;
      
  //VTK render window embedded in a Qt widget
  QVTKWidget *QVTK;
  QVTKWidget *QVTK1;	
  QVTKWidget *QVTK2;
  


  //Qt widgets on the main window
  QCheckBox *AddBox;
  //QPushButton *ClearButton;
  QCheckBox *DeleteBox;
  QCheckBox *MergeBox;
  QCheckBox *SplitBox;
  QCheckBox *UndoDelBox;	
  QPushButton *PlaceButton;
  QPushButton *ApplyButton;


  //Qt widgets for the settings window
  QWidget *SettingsWidget;
  QLineEdit *MaxGapField;
  QLineEdit *GapToleranceField;
  QLineEdit *LineLengthField;
  QLineEdit *ColorValueField;
  QLineEdit *LineWidthField;
  QPushButton *ApplySettingsButton;
  QPushButton *CancelSettingsButton;


//QMessage Box

 QMessageBox* msgBox;
 


	//stuff for tol and selection
  //general render window variables
	vtkSmartPointer<vtkRenderWindowInteractor> Interactor;
	vtkSmartPointer<vtkRenderWindowInteractor> Interactor1;
        vtkSmartPointer<vtkRenderWindowInteractor> Interactor2;
        vtkSmartPointer<vtkActor> LineActor;
	vtkSmartPointer<vtkPolyDataMapper> LineMapper;

  //interactor variables and point picking
	vtkSmartPointer<vtkCallbackCommand> isPicked;
	vtkSmartPointer<vtkCallbackCommand> keyPress;
	
	vtkSmartPointer<vtkCellPicker> CellPicker;
	vtkSmartPointer<vtkPointPicker> PointPicker;
	vtkSmartPointer<vtkPropPicker> JustPicker;
        
        std::vector<int> IDList;

	vtkSmartPointer<vtkSphereSource> Sphere;
	vtkSmartPointer<vtkPolyDataMapper> SphereMapper;
	vtkSmartPointer<vtkActor> SphereActor;
	vtkSmartPointer<vtkPolyDataMapper> DelSphereMapper;
	vtkSmartPointer<vtkActor> DelSphereActor;
	vtkSmartPointer<vtkPolyDataMapper> AddSphereMapper;
	vtkSmartPointer<vtkActor> AddSphereActor;



	vtkSmartPointer<vtkConeSource> Cone;
	vtkSmartPointer<vtkPolyDataMapper> ConeMapper;
	vtkSmartPointer<vtkActor> ConeActor;


  //  img reading	and contour->3d
	vtkSmartPointer<vtkPolyDataMapper> VolumeMapper;
	vtkSmartPointer<vtkActor> VolumeActor;
	//SeedObject* tobj;
	vtkSmartPointer<vtkVolume> Volume;
  //raycast
	vtkSmartPointer<vtkPolyData> poly_line_data;
	vtkSmartPointer<vtkPolyData> poly;
	vtkSmartPointer<vtkPolyDataMapper> polymap;

  // Seeds
        vtkSmartPointer<vtkActor2D> textactor;
        vtkSmartPointer<vtkGlyph3D> Glyph;
	vtkSmartPointer<vtkGlyph3D> delglyph;
        vtkSmartPointer<vtkGlyph3D> addglyph;
        vtkSmartPointer<vtkImageData> VTKim;
        vtkSmartPointer<vtkImageData> VTKim2;
        vtkSmartPointer<vtkImageReslice> Reslice;
        vtkSmartPointer<vtkImageReslice> Reslice1;
        vtkSmartPointer<vtkDataSetMapper>im_mapper;
        vtkSmartPointer<vtkDataSetMapper>im_mapper1;
        vtkSmartPointer<vtkActor> imActor;
        vtkSmartPointer<vtkActor> imActor1;
	double* y;

        vtkSmartPointer<vtkHandleWidget> widget;
        vtkSmartPointer<vtkPointHandleRepresentation3D> handle;


        vtkSmartPointer<vtkHandleWidget> widget1;
        vtkSmartPointer<vtkPointHandleRepresentation2D> handle1;
        vtkSmartPointer<vtkCamera> Int1Cam;
        double* wp;
	vtkSmartPointer<vtkHandleWidget> widget2;
        vtkSmartPointer<vtkPointHandleRepresentation2D> handle2;
        vtkSmartPointer<vtkAppendPolyData> apd;
	};

#endif
