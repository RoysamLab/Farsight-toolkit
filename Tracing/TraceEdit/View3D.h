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
class vtkSphereSource;
class vtkVolume;

class MergeModel;
class QTableView;
class ScatterView;

struct point
{
	float x, y, z;
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
	void generateSeeds();
	std::vector<point> readSeeds(std::string seedSource);
	void rayCast(char* raySource);

	static void PickCell(vtkObject* caller, unsigned long event, void* clientdata, void* callerdata);
	static void HandleKeyPress(vtkObject* caller, unsigned long event, void* clientdata, void* callerdata);

	void HighlightSelected(TraceLine* tline, double SelectColor);
	void DeleteTrace(TraceLine *tline);
	void SelectedComp();
	void ShowMergeStats();
	void CalculateGaps();

	bool CheckFileExists(const char *filename);
	void ShowSeeds();

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

	void ShowLoadSomaWindow();
	void GetSomaPath();
	void GetSomaFile();
	void HideLoadSomaWindow();
	void ToggleSomas();

	void ShowLoadSeedWindow();
	void GetSeedPath();
	void GetSeedFile();
	void HideLoadSeedWindow();
	void ToggleSeeds();

	void ShowSomaSettingsWindow();
	void HideSomaSettingsWindow();
	void ApplySomaSettings();

protected:
	void closeEvent(QCloseEvent *event);
	void Rerender();

private:

	int smallLine;
	float lineWidth;
	double SelectColor;
	std::string SomaFile;
	std::string SeedFile;
	double somaopacity;
	std::vector<point> sdPts;

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

  //qt model view objects
  MergeModel *GapModel;
  QTableView *MergeTableView; 
  ScatterView *MergeScatterView; 

//merge statistics
	QStandardItemModel *model;
	QItemSelectionModel *selModel;
	QTableView *table;
	Qt::SortOrder Ascending;
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
	QAction *loadSeed;
	QAction *somaSettings;
	QAction *viewSomas;
	QAction *viewSeeds;


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
	QWidget *LoadSeedWidget;
	QLineEdit *SeedFileField;
	QPushButton *OpenSeedButton;
	QPushButton *CancelSeedButton;
	QPushButton *BrowseSeedButton;
	QRegExp seedRegex;

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
	std::vector<int> IDList;

    //merge info
	std::vector<TraceGap> candidateGaps;
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
	
	//Seed points
	vtkSmartPointer<vtkFloatArray> pcoords;
	vtkSmartPointer<vtkPoints> points;
	vtkSmartPointer<vtkPolyData> seedPoly;
	vtkSmartPointer<vtkGlyph3D> seedGlyph;
	vtkSmartPointer<vtkLODActor> glyphActor;
	vtkSmartPointer<vtkPolyDataMapper> glyphMapper;
	vtkSmartPointer<vtkSphereSource> sphereSource;
};
#endif
