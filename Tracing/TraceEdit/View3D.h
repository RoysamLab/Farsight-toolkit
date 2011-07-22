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
class CellTrace;
class ImageRenderActors;
class ProjectManager;
class StatisticsToolbar;

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
class vtkImagePlaneWidget;
class vtkWindowToImageFilter;
class vtkJPEGWriter;
  
class TraceModel;
class MergeModel;
class CellTraceModel;
class QTableView;
class ScatterView;
class TableWindow;
class ImageRenderActors;
class vtkPointWidget;
class vtkOrientationMarkerWidget;
class vtkAxesActor;
class vtkQuad;
class vtkLinearExtrusionFilter;
class vtkCellLocator;

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
	void AddDebugPoints(std::vector<TraceBit> vec);// shown as big yellow cubes
	
	static void PickCell(vtkObject* caller, unsigned long event, void* clientdata, void* callerdata);
	static void HandleKeyPress(vtkObject* caller, unsigned long event, void* clientdata, void* callerdata);
	void pointer3DLocation(double pos[]);
	void saveRenderWindow(const char *filename);

	void HighlightSelected(TraceLine* tline, double SelectColor);
	void DeleteTrace(TraceLine *tline);
	void UnSafeDeleteTrace(TraceLine *tline);
	void AddChildren(TraceLine *trunk, std::vector<TraceLine*> childTraces);
	void MergeSelectedTraces();
	void setupLinkedSpace();
	void ShowMergeStats();
	void CalculateGaps();
	void FlipTree(TraceLine* thisLine);
	void TraceBitImageIntensity(int ImgID);
	void CloseTreePlots();
	void HandleHippocampalDataset();
	void smoothzrecursive(TraceLine*,int);
	std::vector<int> getHippocampalTraceIDsToDelete(int,int);
	std::vector<int> getHippocampalTraceIDsToDelete_v2(int,int);
	void DeleteEmptyLeafNodesRecursive(TraceLine*);
	bool readProject(QString projectFile);
	//todo: make these private with accessors
	vtkSmartPointer<vtkRenderer> Renderer;
	vtkSmartPointer<vtkActor> BranchActor;\
	void ShowProjectTable();

public slots:
	void choosetoRender(int row, int col);
	void changeDimension(int row, int col);
	void ListSelections();
	void ClearSelection();
	void SelectTrees();
	void DeleteTraces();
	void SetRoots();
	void AddNewBranches();
	void ExplodeTree();
	void BreakBranch();
	void MergeTraces();
	void SplitTraces();
	void FlipTraces();
	void SaveToFile();
	void SaveProjectFile();
	void SaveSelected();
	void ShowTreeData();
	void ShowCellAnalysis();
	void ShowSettingsWindow();
	void HideSettingsWindow();
	void HideCellAnalysis();
	void ApplyNewSettings();
	void AutomaticEdits();
	void ShowAutomatedEdits();

	void SLine(double d);
	void FakeSpines(double d);
	void FakeBridges(double d);
	void HalfBridges(double d);

	void LoadTraces();
	void LoadImageData();
	void LoadSomaFile();
	void LoadProject();
	void SetTraceType(int newType);
	void ReloadState();
	void SetImgInt();
	void SaveScreenShot();
	void AutoCellExport();

	void chooseInteractorStyle(int iren);
	void SetProjectionMethod(int style);
	void removeImageActors();
	void raycastToSlicer();
	void ToggleColorByTrees();
	void setSlicerZValue(int value);
	void RayCastOpacityChanged(int value);
	void RayCastOpacityValueChanged(double value);
	void RayCastBrightnessChanged(int value);
	void showPTin3D(double value);
	void getPosPTin3D();
	void setPTtoSoma();
	void setUsePointer(int i);
	void createNewTraceBit();

	void DrawROI();

	void focusOn();
	void setRenderFocus(double renderBounds[], int size);
	void FocusOnCell(CellTrace* SelectedCell);
//these are for bootloadfile
	QString getTraceFile();
	QString getImageFile();
	QString getSomaFile();
	void OkToBoot();
	void EditHelp();
	void About();
	void showStatistics(void);
	void updateStatistics(void);

	///////////////////////////////////////
	void rotateImage(int axis);
	//void rotateImageRight90();
	//void rotateImageLeft90();
	//void ElevateImageUp90();
	//void ElevateImageDown90();
	//void AzimuthImageRight90();
	//void AzimuthImageLeft90();

	void rotationOptions();

protected slots:
	void updateSelectionFromCell();
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
	QDockWidget * InformationDisplays, *BootDock, * settingsDock, *cursor3DDock, *AutomationDock, *projectFilesDock;
	QTextEdit * EditLogDisplay;
	
	QTableWidget * projectFilesTable;
	
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
	QCheckBox * Use2DSlicer;
	bool translateImages;
	QDate  Date;
	QTime  Time;
	QDockWidget *statisticsDockWidget;
	int flag;
	
	//Declares an undoBuffer
	//typedef undoBuffer<std::pair<std::string, TraceObject> > bufferType;
	//bufferType *undoBuff;
//settings for display
	int SmallLineLength, maxNumBits, maxPathLength, minDistToParent;
	float lineWidth;
	double SelectColor, backColorR, backColorG, backColorB;
	int numSplit, numDeleted, numMerged;
	QLabel *SplitLabel, *DeleteLabel, *MergeLabel, *BranchesLabel;
	QDoubleSpinBox *posX, *posY, *posZ;
    //VTK render window embedded in a Qt widget
	QVTKWidget *QVTK;	
    QWidget *CentralWidget;
	QMenu *fileMenu;
	QMenu *ShowToolBars;
	QMenu *DataViews;
	QMenu *analysisViews;
	QMenu *help;
	QAction *aboutAction;
	QToolBar *EditsToolBar, *BranchToolBar, *RacastBar, *SlicerBar;

	//Qt widgets on the main window
	QAction *saveAction;
	QAction *saveSelectedAction;
	QAction *exitAction;
	QAction *loadTraceAction;
	QAction *loadTraceImage;
	QAction *CloseAllImage;
	QAction *SetRaycastToSlicer;
	QAction *ColorByTreesAction;
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
	QAction *ScreenshotAction;
	QAction *SettingsButton;
	QPushButton *AutomateButton;
	QAction *BreakButton;
	QAction *explodeTree;
	QAction *root;
	QAction *ImageIntensity;
	QPushButton *MoveSphere;
	QPushButton *updatePT3D;
	QCheckBox *ShowPointer;
	QPushButton *setSoma;
	QPushButton *createNewBitButton;
	QPushButton *createNewROIPointButton;
	QAction *FocusAction;
	QAction *AutoCellExportAction;
	QAction *ShowPlots;
	QAction *CellAnalysis;
	QAction *showStatisticsAction;
	QAction *updateStatisticsAction;
	/////////////////////////
	QAction *rotateleft90;
	QAction *rotateright90;
	QAction *rotatecustom;
	QAction *rotatereset;
	/////////////////////////

	/*QAction *UndoButton;
	QAction *RedoButton;*/
    //merge info
	std::vector<TraceGap*> candidateGaps;
	QString myText;	QString dtext;	QString grayText;

	Qt::SortOrder Ascending;
	TraceModel *TreeModel;
	MergeModel *MergeGaps;
	CellTraceModel *CellModel;
	TableWindow * FTKTable;
	PlotWindow *TreePlot;
	PlotWindow *GapsPlotView;
	TableWindow *GapsTableView;
	TableWindow * FL_MeasureTable;
	PlotWindow *FL_MeasurePlot;
	StatisticsToolbar * statisticsToolbar;

	//Qt widgets for the settings window
	QWidget *SettingsWidget;
	QSpinBox *MaxGapField;
	QDoubleSpinBox *GapToleranceField;
	//QSpinBox *LineLengthField;
	QDoubleSpinBox *ColorValueField;
	QSpinBox *LineWidthField;
	QCheckBox * markTraceBits;
	bool renderTraceBits;
	QDoubleSpinBox *BackgroundRBox,*BackgroundGBox,*BackgroundBBox;
	QDoubleSpinBox *RollBox, *ElevationBox, *AzimuthBox;
	QDialogButtonBox *ApplySettingsButton;
	QComboBox *typeCombo, *StyleCombo, *ProjectionCombo, *RotateImageUpCombo;
	QPushButton *updateRotationButton;
	//QPushButton *RotateRight90Button, *RotateLeft90Button, *ElevateUp90Button, *ElevateDown90Button, *AzimuthRight90Button, *AzimuthLeft90Button;

	//Automation widgets
	QWidget * AutomationWidget;
	//QSpinBox *MaxSpineBit, *MaxBridgeBits, *MaxHalfBridgeBits;
	QDoubleSpinBox *MaxSpinePathLength, * MinDistanceToParent, *MaxSpineBit, *MaxBridgeBits, *MaxHalfBridgeBits, *LineLengthField;
	QRadioButton *SmallLinesButton, *FalseSpinesButton, *FalseBridgesButton, *HalfBridgesButton;
	QGroupBox *SmallLinesGroup, *FakeSpinesGroup, *FakeBridgeGroup, *HalfBridgeGroup;
	//stuff for tol and selection
    //general render window variables
	vtkSmartPointer<vtkRenderWindowInteractor> Interactor;
	vtkSmartPointer<vtkActor> LineActor;
	vtkSmartPointer<vtkActor> PointsActor;
	vtkSmartPointer<vtkPolyDataMapper> LineMapper;

	// save screenshots
	vtkSmartPointer<vtkWindowToImageFilter> WindowToImage;
	vtkSmartPointer<vtkJPEGWriter> JPEGWriter;

    //interactor variables and point picking
	QWidget * CursorActionsWidget;
	vtkSmartPointer<vtkCallbackCommand> isPicked;
	vtkSmartPointer<vtkCallbackCommand> keyPress;
	vtkSmartPointer<vtkCellPicker> CellPicker;

  //ID numbers of the selected traces
	std::vector<int> SelectedTraceIDs;

	vtkSmartPointer<vtkPointWidget> pointer3d;
	bool ShowPointer3DDefault;
	vtkSmartPointer<vtkSphereSource> Sphere;
	vtkSmartPointer<vtkPolyDataMapper> SphereMapper;
	vtkSmartPointer<vtkActor> SphereActor;
	//double pointer3DPos[3];
	std::vector<TraceLine*> stems;
	vtkSmartPointer<vtkAxesActor> axes;
	vtkSmartPointer<vtkOrientationMarkerWidget> UCSMarker;
	
	ImageRenderActors *ImageActors;
	TraceObject* tobj;
  //raycast
	vtkSmartPointer<vtkPolyData> poly;
	vtkSmartPointer<vtkPolyDataMapper> polymap;
	typedef vtkSmartPointer<vtkVolume> VolumePointerType;
	std::vector<VolumePointerType> m_volumes;

	bool viewIn2D;
	int projectionStyle;
	void createSlicerSlider();
	void setSlicerBarValues(int i);
	QSpinBox * SliceSpinBox;
	QSlider * SliceSlider;

	void createRayCastSliders();
	QSpinBox * OpacitySpin, * BrightnessSpin;
	QDoubleSpinBox * OpacityValueSpin;
	QSlider *OpacitySlider, *BrightnessSlider;

	TraceObject * debug_object;
	vtkSmartPointer<vtkActor> debug_actor;

	//double prev_rollAngle;
	//double prev_azimuthAngle;
	//double prev_elevationAngle;
	struct projection {
		double roll;
		double azimuth;
		double elevation;
	} projection_base;
};
#endif
