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

#include "vtkSmartPointer.h"
//#include "UndoBuffer.h"
#include "undobuff.h"
#include <stdio.h>
#include <string>
#include <QtGui>

//forward declarations
class HistoWindow;
class PlotWindow;
class QAction;
class QVTKWidget;
class TraceBit;
class TraceGap;
class TraceLine;
class TraceObject;
class vtkActor;
class vtkCallbackCommand;
class vtkContourFilter;
class vtkCellPicker;
class vtkFloatArray;
class vtkGlyph3D;
class vtkLODActor;
class vtkObject;
class vtkPoints;
class vtkPolyData;
class vtkPolyDataMapper;
class vtkRenderer;
class vtkRenderWindowInteractor;
class  vtkSliderWidget;
class vtkSphereSource;
class vtkVolume;
//class MergeModel; 
class TraceModel;
class MergeModel;
class QTableView;
class ScatterView;
class TableWindow;


class View3D : public QMainWindow 
{
Q_OBJECT;
public:

	vtkSmartPointer<vtkPolyData> poly_line_data;
	View3D(int argc, char **argv);
	~View3D();
	void Initialize();
	void CreateGUIObjects();
	void CreateLayout();
  void CreateModelsAndViews();
	void CreateInteractorStyle();
	void CreateActors();
	void UpdateLineActor();
	void UpdateBranchActor();
	void CreateSphereActor();
	
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
	void AddChildren(TraceLine *trunk, std::vector<TraceLine*> childTraces);
	void MergeSelectedTraces();
	void setupLinkedSpace();
	void ShowMergeStats();
	void CalculateGaps();

	bool CheckFileExists(const char *filename);


	bool setTol();
	//todo: make these private with accessors
	vtkSmartPointer<vtkRenderer> Renderer;
	vtkSmartPointer<vtkActor> BranchActor;

public slots:
	void ListSelections();
	void ClearSelection();
	void DeleteTraces();
	void AddNewBranches();
	void MergeTraces();
	void SplitTraces();
	void FlipTraces();
	void SaveToFile();
	void ShowTreeData();
	void ShowSettingsWindow();
	void HideSettingsWindow();
	void ApplyNewSettings();
	void SLine();

	void ShowLoadSomaWindow();
	void GetSomaPath();
	void GetSomaFile();
	void HideLoadSomaWindow();
	void ToggleSomas();

	void ShowSomaSettingsWindow();
	void HideSomaSettingsWindow();
	void ApplySomaSettings();
	//void dataChanged(const QModelIndex &topLeft, const QModelIndex &bottomRight);
    void DereferenceGapsPlotView();
	void DereferenceTreePlotView();

	void UndoAction();
	void RedoAction();


protected slots:
	void updateSelectionHighlights();
	void updateTraceSelectionHighlights();

protected:
	void closeEvent(QCloseEvent *event);
	void Rerender();

private:

	//Declares an undoBuffer
	typedef undoBuffer<std::pair<std::string, TraceObject> > bufferType;
	bufferType *undoBuff;

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

	QAction *ListButton;
	QAction *ClearButton;
	QAction *DeleteButton;
	QAction *MergeButton;
	QAction *BranchButton;
	QAction *SplitButton;
	QAction *FlipButton;
	QAction *WriteButton;
	QAction *SettingsButton;
	QAction *AutomateButton;

	QAction *UndoButton;
	QAction *RedoButton;

  //qt model view objects
  MergeModel *GapModel;
  QTableView *MergeTableView; 
  ScatterView *MergeScatterView; 

//merge statistics	
	QTableView *GapsTableView;
	Qt::SortOrder Ascending;
//plots	
	PlotWindow *GapsPlotView;
	HistoWindow *histo;
//tobj statistics	
	TraceModel *TreeModel;
	QTableView *TreeTable;
	PlotWindow *TreePlot;
	//QT widgets for the menu bar
	QMenu *fileMenu;
	QToolBar *EditsToolBar;

	QAction *saveAction;
	QAction *exitAction;

	QMenu *somaMenu;
	QAction *loadSoma;
	QAction *loadSeed;
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
	QRegExp somaRegex;

	//Qt Widgets for the seed point file reader window

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

  //ID numbers of the selected traces
	std::vector<int> SelectedTraceIDs;

    //merge info
	std::vector<TraceGap*> candidateGaps;
	QString myText;	QString dtext;	QString grayText;
	MergeModel *MergeGaps;

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
	//vtkSmartPointer<vtkPolyData> poly_line_data;
	vtkSmartPointer<vtkPolyData> poly;
	vtkSmartPointer<vtkPolyDataMapper> polymap;

  //VTK widgets
  vtkSliderWidget *OpacitySlider;
  vtkSliderWidget *BrightnessSlider;
};
#endif
