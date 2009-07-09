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
#include <QAction>
#include <QtGui>
#include <QVTKWidget.h>
#include <string>
#include "ftkGUI/PlotWindow.h"
#include "ftkGUI/HistoWindow.h"

struct compTrace{
TraceLine *Trace1;
TraceLine *Trace2;
int endPT1, endPT2;
double angle;
double dist; 
double maxdist;
double length;
double smoothness;
double cost;
};

class View3D : public QMainWindow 
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
	bool CheckFileExists(const char *filename);
	
	void AddPointsAsPoints (std::vector<TraceBit> vec);
	void AddVolumeSliders();
	void AddContourThresholdSliders();
	void AddPlaybackWidget(char*);
	void readImg(std::string sourceFile);
	void rayCast(char* raySource);

	static void PickCell(vtkObject* caller, unsigned long event, void* clientdata, void* callerdata);
	static void HandleKeyPress(vtkObject* caller, unsigned long event, void* clientdata, void* callerdata);

	void HighlightSelected(TraceLine* tline, double SelectColor);
	void DeleteTrace(TraceLine *tline);
	void MinEndPoints(std::vector<TraceLine*> traceList);
	void SelectedComp();
	void ShowMergeStats();
	void traceStatistics();

	bool setTol();
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
	void SLine();
	void ToggleSomas();
	void ShowLoadSomaWindow();
	void GetSomaPath();
	void GetSomaFile();
	void HideLoadSomaWindow();
	void ShowSomaSettingsWindow();
	void HideSomaSettingsWindow();
	void ApplySomaSettings();

protected:
	void closeEvent(QCloseEvent *event);
	void Rerender();

private:
	double gapTol;
	int gapMax;
	int smallLine;
	float lineWidth;
	double SelectColor;
	std::string SomaFile;
	double somaopacity;

    //VTK render window embedded in a Qt widget
	QVTKWidget *QVTK;
	QAction *loadAction;

	//Qt widgets on the main window
    QWidget *CentralWidget;

	QPushButton *ListButton;
	QPushButton *ClearButton;
	QPushButton *DeleteButton;
	QPushButton *MergeButton;
	QPushButton *SplitButton;
	QPushButton *FlipButton;
	QPushButton *WriteButton;
	QPushButton *SettingsButton;
	QPushButton *AutomateButton;
//merge statistics
	QStandardItemModel *model;
	QItemSelectionModel *selModel;
	QTableView *table;
//plots	
	PlotWindow *plot;
	HistoWindow *histo;
//tobj statistics	
	QStandardItemModel *treeModel;
	QItemSelectionModel *TreeSelModel;
	QTableView *TreeTable;

	//QT widgets for the menu bar
	QMenu *fileMenu;
	QAction *exitAction;

	QMenu *somaMenu;
	QAction *loadSoma;
	QAction *somaSettings;
	QAction *viewSomas;

	//Qt widgets for the settings window
	QWidget *SettingsWidget;
	QLineEdit *MaxGapField;
	QLineEdit *GapToleranceField;
	QLineEdit *LineLengthField;
	QLineEdit *ColorValueField;
	QLineEdit *LineWidthField;
	QPushButton *ApplySettingsButton;
	QPushButton *CancelSettingsButton;

	//Qt Widgets for the soma file reader window
	QWidget *LoadSomaWidget;
	QLineEdit *SomaFileField;
	QPushButton *OpenSomaButton;
	QPushButton *CancelSomaButton;
	QPushButton *BrowseSomaButton;
	QRegExp regex;


	//Qt Widgets for the soma settings window
	QWidget *SomaSettingsWidget;
	QLineEdit *SomaOpacityField;
	QPushButton *ApplySomaSettingsButton;
	QPushButton *CancelSomaSettingsButton;


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

    //merge info
	std::vector<compTrace> compList;
	std::vector<compTrace> grayList;
	QString myText;	QString dtext;	QString grayText;

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
