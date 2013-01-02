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

#ifndef TRACEVIEW3D_H_
#define TRACEVIEW3D_H_

//#include "UndoBuffer.h"
//#include "undobuff.h"
#include <stdio.h>
#include <string>
#include <set>
#include <QtGui>
#include <QInputDialog>
#include <time.h>

#include "ftkGUI/PlotWindow.h"
#include "ftkGUI/HistoWindow.h"
#include "ftkGUI/TableWindow.h"
#include "ftkGUI/FTKRenderWindow.h"
#include "ftkGUI/StatisticsToolbar.h"
#include "itkImageFileReader.h"
#include "itkImageToVTKImageFilter.h"
#include "vnl/vnl_cost_function.h"
#include "vnl/algo/vnl_conjugate_gradient.h"
#include "vnl/algo/vnl_powell.h"
#include <boost/math/special_functions.hpp>

#include <QVTKWidget.h>

#include "vtkActor.h"
#include "vtkAxesActor.h"
#include "vtkCallbackCommand.h"
#include "vtkCamera.h"
#include "vtkCellArray.h"
#include "vtkCellLocator.h"
#include "vtkCellPicker.h"
#include "vtkColorTransferFunction.h"
#include "vtkCommand.h"
#include "vtkContourFilter.h"
#include "vtkCubeSource.h"
#include "vtkFloatArray.h"
#include "vtkGlyph3D.h"
#include "vtkImageData.h"
#include "vtkImageToStructuredPoints.h"
#include "vtkInteractorStyleImage.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkInteractorStyleRubberBandZoom.h"
#include "vtkInteractorStyleImage.h"
#include "vtkImagePlaneWidget.h"
#include "vtkPNGWriter.h"
#include "vtkLegendScaleActor.h"
#include "vtkLinearExtrusionFilter.h"
#include "vtkLODActor.h"
#include "vtkObjectFactory.h"
#include "vtkOBJReader.h"
#include "vtkOrientationMarkerWidget.h"
#include "vtkPiecewiseFunction.h"
#include "vtkPlaybackRepresentation.h"
#include "vtkPlaybackWidget.h"
#include "vtkPlotGrid.h"
#include "vtkPointData.h"
#include "vtkPointPicker.h"
#include "vtkPoints.h"
#include "vtkPointWidget.h"
#include "vtkPolygon.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkRendererCollection.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkSliderRepresentation2D.h"
#include "vtkSliderWidget.h"
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkVolume.h"
#include "vtkVolumeProperty.h"
#include "vtkWindowToImageFilter.h"

#include "PatternAnalysis/activeLearning/mclr_SM.h"
#include "ftkGUI/GenericALDialog.h"
#include "ftkGUI/TrainingDialog.h"
#include "ftkGUI/PatternAnalysisWizard.h"
#ifdef USE_SPD
#include "SPD/spdtestwindow.h"
//#include "SPD/spdtestwindowForNewSelection.h"
#endif
#ifdef USE_Clusclus
#include "ClusClus/HeatmapWindow.h"
#include "ClusClus/Heatmap.h"
#endif
#include "ftkVesselTracer/ftkVesselTracer.h"
#include "branchPT.h"
#include "CellTrace.h"
#include "CellTraceModel.h"
#include "ftkCommon/ftkProjectManager.h"
#include "ftkUtils.h"
#include "ImageActors.h"
#include "MergeModel.h"

#include "SelectiveClustering.h"
#include "QvtkTableView.h"
#include "SelectionUtilities.h"

#include "cellexport.h"
#include "FeatureRelation.h"
#include "GridActors.h"
#include "screenshot.h"

#include "TraceBit.h"
#include "TraceGap.h"
#include "TraceLine.h"
#include "TraceObject.h"
#include "TraceModel.h"
#include "NodeModel.h"

#include "VolumeOfInterest.h"
  
#ifdef USE_QT_TESTING
#include "GUITester.h"
#endif

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
	void AddDebugPoints(vtkSmartPointer<vtkTable> centroidsTable);// shown as big yellow cubes
	
	static void PickCell(vtkObject* caller, unsigned long event, void* clientdata, void* callerdata);
	static void HandleKeyPress(vtkObject* caller, unsigned long event, void* clientdata, void* callerdata);
	void pointer3DLocation(double pos[]);
	void saveRenderWindow(const char *filename);
	QImage Get_AL_Snapshot(CellTrace* currentCell);
	QImage vtkImageDataToQImage(vtkImageData * imageData);

	void HighlightSelected(TraceLine* tline, double SelectColor);
	void DeleteTrace(TraceLine *tline);
	void AddChildren(TraceLine *trunk, std::vector<TraceLine*> childTraces);
	void MergeSelectedTraces();
	void setupLinkedSpace();
	void ShowMergeStats();
	void CalculateGaps();
	void FlipTree(TraceLine* thisLine);
	void TraceBitImageIntensity(int ImgID);
	void TraceBitImageIntensityWeighted(int ImgID);
	void CloseTreePlots();
	void CloseNodePlots();

	bool readProject(QString projectFile);
	//todo: make these private with accessors
	vtkSmartPointer<vtkRenderer> Renderer;
	vtkSmartPointer<vtkActor> BranchActor;
	void ShowProjectTable();
	void AssociateDebrisToNuclei( vtkSmartPointer<vtkTable> debrisTable);

public slots:
	void choosetoRender(int row, int col);
	void changeDimension(int row, int col);
	void ListSelections();
	void ClearSelection();
	void FastClearSelection();
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
	void CalculateDelaunay3D();
	void ShowDelaunay3D();
	void ShowEllipsoid();
	void ShowNodeData();
	void ShowTreeData();
	void ShowCellAnalysis();
	void ShowSettingsWindow();
	void HideSettingsWindow();
	void adjustEditorSettingsSize(bool changesize);
	void HideCellAnalysis();
	void ApplyNewSettings();
	void AutomaticEdits();
	void ShowAutomatedEdits();

	void IntensityFeature();
	void IntensityWeightedFeature();

	//tracing gui slots
	void PickTracer(int choice);
	void RunTracer();
	void StartMNTracerAmit(int costThreshold);
	void StartVesselBallTracer();
	void seedPointFileDialog();
	void somaFileDialog();

	void SLine(double d);
	void FakeSpines(double d);
	void FakeBridges(double d);
	void HalfBridges(double d);

	void openTracingDialog();

	void LoadTraces();
	void LoadImageData();
	void LoadSomaFile();
	void LoadProject();
	void SetTraceType(int newType);
	void ReloadState();
	void SetImgInt();
	void SetImgWeightInt();
	void SaveScreenShot();
	void AutoCellExport();

	void StartActiveLearning();

	void setHighlightSettings(int value);
	void chooseInteractorStyle(int iren);
	void SetProjectionMethod(int style);
	void removeImageActors();
	//void raycastToSlicer();
	void ToggleColorByTrees();
	void ToggleGridlines();
	void AdjustGridlines(int value);
	void setSlicerZValue(int value);
	void setSliceThickness(int sliceThickness);
	void setSliceWindowLevel(int value);
	void RayCastOpacityChanged(int value);
	void RayCastOpacityValueChanged(double value);
	void RayCastBrightnessChanged(int value);
	void RayCastColorValueChanged(int value);
	void chooseSomaRender(int value);
	void SomaOpacityChanged(int value);
	void SomaOpacityValueChanged(double value);
	void SomaBrightnessChanged(int value);
	void SomaColorValueChanged(double value);
	void showPTin3D(double value);
	void getPosPTin3D();
	void setPTtoSoma();
	void ArunCenterline();
	void setUsePointer(int i);
	void createNewTraceBit();

	void AddROIPoint();
	void DrawROI();
	void ReadVOI();
	void WriteVOI();
	void ToggleVOI();
	void CalculateDistanceToDevice();

	void CalculateCellToCellDistanceGraph();
	void readNucleiTable();
	void AssociateNeuronToNuclei();
	void ShowSeedPoints();

	void readDebrisTable();

	void focusOn();
	void setRenderFocus(double renderBounds[], int size);
	void FocusOnCell(CellTrace* SelectedCell);

	void SPDAnalysis();
	void AddLabel();
	void ClusclusAnalysis();
	void Highlighted_selected();
	void BiclusAnalysis();
	void Getlabel(std::set<long int> labels);
	void FeatureDistributionAnalysis();
	void KNearestNeighborAnalysis();
	double ComputeAverageDistance(std::vector< std::pair<unsigned int, double> >);
	void NeighborsWithinRadiusAnalysis();
	void SpectralCluserting();
	void selectedFeaturesClustering();

//these are for bootloadfile
	QString getTraceFile();
	QString getImageFile();
	QString getSomaFile();
	void OkToBoot();
	void EditHelp();
	void About();
	void showStatistics(void);
	void updateStatistics(void);

	void rotateImage(int axis);
	void rotationOptions();

  int runTests();

protected slots:
	void updateSelectionFromCell();
	void updateSelectionHighlights();
	void updateTraceSelectionHighlights();
	void updateNodeSelection();
	void CropBorderCells();
	void SaveComputedCellFeaturesTable();

	void setSlicerMode();
	void setProjectionMode();
	void setRaycastMode();

	void activateSaveAllButton();

	void setContourMode();
	void setRaycastSomaMode();
	void clearSettings();
	void recordTest();
	void resizeForTesting();
	void UpdateKeepSelectedTraces();

protected:
	void closeEvent(QCloseEvent *event);
	void Rerender();

private:

	//for the widget the buttons must be used instead of actions 
	//the get____file functions called as renderer is not initalized
	double uMperVoxel; //0 x, 1 y, 2 z
	QSettings TraceEditSettings;
	QDockWidget * InformationDisplays, *BootDock, * settingsDock, *cursor3DDock, *AutomationDock, *projectFilesDock, *vesselSegDock;
	QTextEdit * EditLogDisplay;
	QProgressBar * ProgressBar;
        QLabel *       ProgressDescription;
	
	QString UserName, LabName, ProjectName;
	QStringList TraceFiles, Image, SomaFile, tempTraceFile, NucleiFile;
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
	QMenu * soma_sub_menu;
	bool translateImages;
	QDate  Date;
	QTime  Time;
	QDockWidget *statisticsDockWidget;
	int flag;
	
//settings for display
	int SmallLineLength, maxNumBits, maxPathLength, minDistToParent;
	float lineWidth;
	double SelectColor, SelectTipColor, backColorR, backColorG, backColorB;
	QGroupBox * selectionSettings, * displaySettings, * rotationSettings, * BackgroundSettings, * GridlineSettings;
	int numSplit, numDeleted, numMerged;
	QLabel *SplitLabel, *DeleteLabel, *MergeLabel, *BranchesLabel;
	QDoubleSpinBox *posX, *posY, *posZ;
    //VTK render window embedded in a Qt widget
	QVTKWidget *QVTK;	
	QWidget *CentralWidget;
	QMenu *fileMenu;
	QMenu *ShowToolBars;
	QMenu *processingMenu;
	QMenu *DataViews;
	QMenu *analysisViews;
	QMenu *help;
	QDialog *tracingGui;
	QAction *aboutAction;
	QToolBar *EditsToolBar, *BranchToolBar, *RaycastBar, *SomaBar, *SlicerBar;

  #ifdef USE_QT_TESTING
  //testing GUI elements
	QMenu *testingMenu;
  QAction *recordAction;
  QAction *playAction;
  QAction *clearAction;
  QAction *resizeAction;
  #endif

	//Qt widgets on the main window
	QAction *saveAction;
	QAction *SaveComputedCellFeaturesTableAction;
	QAction *saveSelectedAction;
	QAction *saveProjectAction;
	QAction *exitAction;
	QAction *loadTraceAction;
	QAction *loadTraceImage;
	QAction *CloseAllImage;

	QAction *SetSlicer;
	QAction *SetProjection;
	QAction *SetRaycast;

	QAction *SetContour;
	QAction *SetSomaRaycast;

	QPushButton *ArunVesselTracingButton;
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
	QAction *ImageWeightedIntensity;
	QAction *GridAction;

	QPushButton *MoveSphere;
	QPushButton *updatePT3D;
	QCheckBox *ShowPointer;
	QPushButton *setSoma;
	QPushButton *createNewBitButton;

	QPushButton *createNewROIPointButton;
	QPushButton *ExtrudeROIButton;
	QPushButton *ReadBinaryVOIButton;
	QPushButton *WriteVOIButton;
	QPushButton *ToggleBinaryVOIButton;
	QPushButton *CalculateDistanceToDeviceButton;
	QPushButton *CalculateCellDistanceButton;

	QAction *LoadNucleiTable;
	QAction *AssociateCellToNucleiAction;
	QAction *LoadSeedPointsAsGlyphs;

	QAction *LoadDebrisTable;
	QAction *ConvexHullAction;

	QAction *FocusAction;
	QAction *AutoCellExportAction;
	QAction *ShowPlots;
	QAction *ShowNodePlots;
	QAction *CellAnalysis;
	QAction *StartActiveLearningAction;
	QAction *showStatisticsAction;
	QAction *updateStatisticsAction;
	QAction *SPDAction;
	QAction *SPDAnalysisAction;
	QAction *AddLabelAction;
	QAction *ClusclusAction;
	QAction *selectedAction;
	QAction *BiClusAction;
	QAction *FeatureDistributionAction;
	QAction *SpectralClusteringAction;
	QAction *KNearestAction;
	QAction *InRadiusAction;

	/*QAction *UndoButton;
	QAction *RedoButton;*/
    //merge info
	std::vector<TraceGap*> candidateGaps;
	QString myText;	QString dtext;	QString grayText;

	Qt::SortOrder Ascending;
	NodeModel *Node_Model;
	TraceModel *TreeModel;
	MergeModel *MergeGaps;
	CellTraceModel *CellModel;
	TableWindow * NodeTable;
	PlotWindow *NodePlot;
	TableWindow * FTKTable;
	PlotWindow *TreePlot;
	PlotWindow *GapsPlotView;
	TableWindow *GapsTableView;
	TableWindow * FL_MeasureTable;
	PlotWindow *FL_MeasurePlot;
	HistoWindow *FL_histo;
	FTKRenderWindow *FeatDistWin;
	StatisticsToolbar * statisticsToolbar;

	bool bKeepSelectedTraces;   /// whether keep the last rendered traces

#ifdef USE_SPD
	SPDtestWindow *SPDWin;
	//SPDWindowForNewSelection *SPDWin;
#endif

#ifdef	USE_Clusclus
	Heatmap *HeatmapWin;
	Heatmap *HeatmapWins;
	BiHeatmap *Biheatmap;
#endif

	bool renderTraceBits;
	bool projectLoadedState;
	//Qt widgets for the settings window
	QWidget *SettingsWidget;
	QToolBox *SettingsToolBox;
	QSpinBox *MaxGapField;
	QDoubleSpinBox *GapToleranceField;
	//QSpinBox *LineLengthField;
	QDoubleSpinBox *ColorValueField;
	QDoubleSpinBox *TipColor;
	QSpinBox *LineWidthField;
	QCheckBox * markTraceBits, * convexHull, * ellipsoid;
	QDoubleSpinBox *BackgroundRBox,*BackgroundGBox,*BackgroundBBox;
	QDoubleSpinBox *RollBox, *ElevationBox, *AzimuthBox;
	QDialogButtonBox *ApplySettingsButton;
	QComboBox *typeCombo, *HighlightCombo, *StyleCombo, *ProjectionCombo, *RotateImageUpCombo;
	QPushButton *updateRotationButton;
	QTabWidget * tabWidget;
	QWidget * gridLineTab;
	int gridTabIndex;
		
	QTableWidget * projectFilesTable;
	bool projectFilesTableCreated;
	//QTableWidgetItem * Item2D, * Item3D;
	
	//Automation widgets
	QWidget * AutomationWidget;
	//QSpinBox *MaxSpineBit, *MaxBridgeBits, *MaxHalfBridgeBits;
	QDoubleSpinBox *MaxSpinePathLength, * MinDistanceToParent, *MaxSpineBit, *MaxBridgeBits, *MaxHalfBridgeBits, *LineLengthField;
	QRadioButton *SmallLinesButton, *FalseSpinesButton, *FalseBridgesButton, *HalfBridgesButton;
	QGroupBox *SmallLinesGroup, *FakeSpinesGroup, *FakeBridgeGroup, *HalfBridgeGroup, *BorderCellsCroppingGroup;
	QPushButton *CropBorderCellsButton;
	//stuff for tol and selection
    //general render window variables
	vtkSmartPointer<vtkRenderWindowInteractor> Interactor;
	vtkSmartPointer<vtkActor> LineActor;
	vtkSmartPointer<vtkActor> PointsActor;
	vtkSmartPointer<vtkActor> CentroidsActor;
	vtkSmartPointer<vtkPolyDataMapper> LineMapper;

	// save screenshots
	vtkSmartPointer<vtkWindowToImageFilter> WindowToImage;
	vtkSmartPointer<vtkPNGWriter> PNGWriter;
	ScreenShotDialog * savescreenshotDialog;

    //interactor variables and point picking
	QWidget * CursorActionsWidget;
	vtkSmartPointer<vtkCallbackCommand> isPicked;
	vtkSmartPointer<vtkCallbackCommand> keyPress;
	vtkSmartPointer<vtkCellPicker> CellPicker;

	//tracing gui widgets
	QComboBox * tracerCombo;
	QGroupBox * mntBox;
	QSpinBox * mntCostThreshold;
	QPushButton * mntSeedsButton;
	QPushButton * mntSomaButton;
	QString seedsFile;
	QString somaFile;
	QGroupBox *vbtBox;
	QCheckBox *vbtPreprocess;
	QSpinBox *vbtUseVesselness;
	QFormLayout *vbtLayout;
	
	//vessel segmentation
	QWidget * vesselSegWidget;

	// Vessel segmentation using ball tracing
	ftkVesselTracer *VBT;

	//ID numbers of the selected traces
	std::vector<int> SelectedTraceIDs;
	std::map<int, int> indMapFromIndToVertex;

	// for nucleus
	vtkSmartPointer<vtkUnsignedCharArray> colors;
	vtkSmartPointer<vtkPolyDataMapper>  spheremap;
	vtkSmartPointer<vtkPolyData> point_poly;
	vtkSmartPointer<vtkGlyph3D> glyphs;

	//
	MCLR_SM *mclr;
	GenericALDialog *ALDialog;

	//
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
	GridlineActors *Gridlines;
	TraceObject* tobj;
	vtkSmartPointer<vtkTable> nucleiTable;
//!raycast Functions 
	vtkSmartPointer<vtkPolyData> poly;
	vtkSmartPointer<vtkPolyDataMapper> polymap;
	typedef vtkSmartPointer<vtkVolume> VolumePointerType;
	std::vector<VolumePointerType> m_volumes;
	vtkSmartPointer<vtkOBJReader> VolReader_OBJ;
	vtkSmartPointer<vtkLegendScaleActor> scalebar;
//! Set up 2d slicer controls
	bool viewIn2D;
	int projectionStyle;
	void createSlicerSlider();
	bool SlicerBarCreated;
	void setSlicerBarValues(int i);
	QSpinBox * SliceSpinBox, * SliceThicknessSpinBox;
	QSlider * SliceSlider, * SliceBrightnessSlider;

	void createRayCastSliders(), createSomaSliders();
	QSpinBox * OpacitySpin, * BrightnessSpin, * SomaOpacitySpin, * SomaBrightnessSpin;
	QDoubleSpinBox * OpacityValueSpin, * SomaOpacityValueSpin, * SomaColorSpin;
	QSlider *OpacitySlider, *BrightnessSlider, *SomaOpacitySlider, *SomaBrightnessSlider;
	QComboBox *ColorProfileCombo;

	struct projection {
		double roll;
		double azimuth;
		double elevation;
	} projection_base;
	int projection_axis;
//!Gridline controls
	bool gridShown;
	QComboBox * GridOrientationBox, * GridDimensionMode;
	QSpinBox * HeightSpaceBox, * WidthSpaceBox, * DepthSpaceBox, * LineWidthBox;
	QSlider * GridRSlider, * GridGSlider, * GridBSlider;
	QSlider * GridOpacitySlider;
//!Delaunay triangulation - 3D convex hull using terminal tips
	std::vector<CellTrace*> delaunayCellsSelected;
	std::vector<CellTrace*> ellipsoidCellsSelected;

//!Feature Graph
	FeatureRelation * feature;

//!ROI data objects
	std::vector<double*> ROIPoints;
	//vtkSmartPointer<vtkPolyData> ROIExtrudedpolydata;	
	vtkSmartPointer<vtkActor> ROIactor;

	VolumeOfInterest * VOIType;
	bool bshowDevice;
	/*! \enum RenderModeEnum
	* \brief Sets default the rendering style
	* Slicer shows a 2d slice at a time
	* Projection 2d orthographic max min or median intensity projections
	* Raycast 3d volumetric rendering
	*/
	enum RenderModeEnum { SLICER, PROJECTION, RAYCAST, SLICERRAYCAST};
	enum HighlightModeEnum { TREE, SEGMENT, TIP };
	bool viewContour;

	enum RenderModeEnum renderMode;
	enum HighlightModeEnum highlightMode;

	void ClearRenderer(int i);

  #ifdef USE_QT_TESTING
  //testing support
  GUITester *Tester;
  QString TestInputFile;
  QString TestBaselineImageFileName;
  #endif
  bool SaveSettingsOnExit;

	std::string imageFileName;
};

class QueryDialog : public QDialog
{
	Q_OBJECT
public:
	QueryDialog(int QueryType, QVector<QString> classes, bool diffusion_graph, QWidget *parent = 0);
	std::vector<unsigned int> parseIDs(void);
	unsigned int parseK(void);
	double parseRad(void);
	unsigned short getSourceClass(void);
	unsigned short getDestClass(void);
	bool getKMutual(void);

private:
	QLabel * idLabel;
	QLineEdit * IDs;
	QHBoxLayout *idLayout;
	QLabel * kLabel;
	QLineEdit * K;
	QHBoxLayout *kLayout;
	QLabel * radLabel;
	QLineEdit * RAD;
	QHBoxLayout *radLayout;
	QLabel * classLabel1;
	QComboBox * classCombo1;
	QHBoxLayout * classLayout1;
	QLabel * classLabel2;
	QComboBox * classCombo2;
	QHBoxLayout * classLayout2;
	QCheckBox * check;
	QPushButton * okButton;
	QHBoxLayout * bLayout;
	QVBoxLayout * layout;
};

#endif
