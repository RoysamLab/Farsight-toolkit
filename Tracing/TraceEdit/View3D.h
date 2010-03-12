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
//#include "undobuff.h"
#include <stdio.h>
#include <string>
#include <set>
#include <QtGui>
#include <QDate>
#include <QTime>

//forward declarations
class HistoWindow;
class PlotWindow;
class QAction;
class QVTKWidget;
class TraceBit;
class TraceGap;
class TraceLine;
class TraceObject;
class branchPT;
class ImageRenderActors;
class ProjectManager;

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
class vtkSliderWidget;
class vtkSphereSource;
class vtkVolume;
  
class TraceModel;
class MergeModel;
class QTableView;
class ScatterView;
class TableWindow;
class ImageRenderActors;


class View3D : public QMainWindow 
{
Q_OBJECT;
public:

	vtkSmartPointer<vtkPolyData> poly_line_data;
	View3D(QWidget * parent = 0);
	View3D(TraceObject* Traces);
	~View3D();
	void Initialize();
  void LoadFiles(int argc, char **argv);
	void CreateBootLoader();
	void CreateGUIObjects();
	void CreateLayout();
	void CreateInteractorStyle();
	void CreateActors();
	void UpdateLineActor();
	void UpdateBranchActor();
	void CreateSphereActor();
	
	void AddPointsAsPoints (std::vector<TraceBit> vec);
	
	static void PickCell(vtkObject* caller, unsigned long event, void* clientdata, void* callerdata);
	static void HandleKeyPress(vtkObject* caller, unsigned long event, void* clientdata, void* callerdata);

	void HighlightSelected(TraceLine* tline, double SelectColor);
	void DeleteTrace(TraceLine *tline);
	void AddChildren(TraceLine *trunk, std::vector<TraceLine*> childTraces);
	void MergeSelectedTraces();
	void setupLinkedSpace();
	void ShowMergeStats();
	void CalculateGaps();
	void TraceBitImageIntensity(int ImgID);

	//todo: make these private with accessors
	vtkSmartPointer<vtkRenderer> Renderer;
	vtkSmartPointer<vtkActor> BranchActor;

public slots:
	void ListSelections();
	void ClearSelection();
	void SelectTrees();
	void DeleteTraces();
	void SetRoots();
	void AddNewBranches();
	void ExplodeTree();
	void MergeTraces();
	void SplitTraces();
	void FlipTraces();
	void SaveToFile();
	void ShowTreeData();
	void ShowSettingsWindow();
	void HideSettingsWindow();
	void ApplyNewSettings();
	void SLine();
	void LoadTraces();
	void LoadImageData();
	void LoadSomaFile();
	void LoadProject();
	void SetTraceType(int newType);
	void ReloadState();
	void SetImgInt();

//these are for bootloadfile
	QString getTraceFile();
	QString getImageFile();
	QString getSomaFile();
	void OkToBoot();

protected slots:
	void updateSelectionHighlights();
	void updateTraceSelectionHighlights();

protected:
	void closeEvent(QCloseEvent *event);
	void Rerender();

private:

	//for the widget the buttons must be used instead of actions 
	//the get____file functions called as renderer is not initalized
	double uMperVoxel; //0 x, 1 y, 2 z
	QSettings TraceEditSettings;
	QDockWidget * InformationDisplays, *BootDock, * settingsDock;
	QTextEdit * EditLogDisplay;
	QString UserName, LabName, ProjectName;
	QStringList TraceFiles, Image, SomaFile, tempTraceFile;
	QWidget * bootLoadFiles;
	QPushButton * BootTrace;
	QPushButton * BootSoma;
	QPushButton * BootImage;
	QPushButton * BootProject;
	QPushButton * okBoot;
	QPushButton * Reload;
	QDoubleSpinBox *scale;
	QComboBox * GetAUserName, *GetLab, *GetProject;
	QGroupBox * ImageListBox;
	QListView * ImageListView; 
	bool translateImages;
	QDate  Date;
	QTime  Time;

	//Declares an undoBuffer
	//typedef undoBuffer<std::pair<std::string, TraceObject> > bufferType;
	//bufferType *undoBuff;
//settings for display
	int SmallLineLength;
	float lineWidth;
	double SelectColor;
	int numSplit, numDeleted, numMerged;
	QLabel *SplitLabel, *DeleteLabel, *MergeLabel;
    //VTK render window embedded in a Qt widget
	QVTKWidget *QVTK;	
    QWidget *CentralWidget;
	QMenu *fileMenu;
	QMenu *ShowToolBars;
	QToolBar *EditsToolBar;	
	QToolBar *BranchToolBar;

	//Qt widgets on the main window
	QAction *saveAction;
	QAction *exitAction;
	QAction *loadTraceAction;
	QAction *loadTraceImage;
	QAction *loadSoma;
	QAction *ListButton;
	QAction *ClearButton;
	QAction *SelectTreeAction;
	QAction *DeleteButton;
	QAction *MergeButton;
	QAction *BranchButton;
	QAction *SplitButton;
	QAction *FlipButton;
	QAction *WriteButton;
	QAction *SettingsButton;
	QAction *AutomateButton;
	QAction *explodeTree;
	QAction *root;
	QAction *ImageIntensity;

	/*QAction *UndoButton;
	QAction *RedoButton;*/
    //merge info
	std::vector<TraceGap*> candidateGaps;
	QString myText;	QString dtext;	QString grayText;

	Qt::SortOrder Ascending;
	TraceModel *TreeModel;
	MergeModel *MergeGaps;
	TableWindow * FTKTable;
	PlotWindow *TreePlot;
	PlotWindow *GapsPlotView;
	TableWindow *GapsTableView;

	//Qt widgets for the settings window
	QWidget *SettingsWidget;
	QLineEdit *MaxGapField;
	QLineEdit *GapToleranceField;
	QSpinBox *LineLengthField;
	QSpinBox *ColorValueField;
	QSpinBox *LineWidthField;
	QDialogButtonBox *ApplySettingsButton;
	QComboBox *typeCombo;
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

	vtkSmartPointer<vtkSphereSource> Sphere;
	vtkSmartPointer<vtkPolyDataMapper> SphereMapper;
	vtkSmartPointer<vtkActor> SphereActor;
	
	ImageRenderActors *ImageActors;
 //  img reading and contour->3d
	TraceObject* tobj;
	/*vtkSmartPointer<vtkPolyDataMapper> VolumeMapper;
	vtkSmartPointer<vtkActor> VolumeActor;
	vtkSmartPointer<vtkVolume> Volume;
	vtkSmartPointer<vtkContourFilter> ContourFilter;*/

  //raycast
	vtkSmartPointer<vtkPolyData> poly;
	vtkSmartPointer<vtkPolyDataMapper> polymap;

  //VTK widgets
  QSlider *OpacitySlider;
  QSlider *BrightnessSlider;
};
#endif
