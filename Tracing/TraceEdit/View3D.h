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
	void LoadTraces();
	void LoadImageData();
	void GetSomaFile();

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

    //VTK render window embedded in a Qt widget
	QVTKWidget *QVTK;
	QAction *loadTraceAction;
	QAction *loadTraceImage;
	QAction *saveAction;
	QAction *exitAction;
	QAction *loadSoma;
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
    //merge info
	std::vector<TraceGap*> candidateGaps;
	QString myText;	QString dtext;	QString grayText;
	
	QMenu *fileMenu;
	QToolBar *EditsToolBar;	
	Qt::SortOrder Ascending;
	TraceModel *TreeModel;
	MergeModel *MergeGaps;
	PlotWindow *GapsPlotView;
	QTableView *GapsTableView;
	//HistoWindow *histo;
	QTableView *TreeTable;
	PlotWindow *TreePlot;

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

  //ID numbers of the selected traces
	std::vector<int> SelectedTraceIDs;


	vtkSmartPointer<vtkSphereSource> Sphere;
	vtkSmartPointer<vtkPolyDataMapper> SphereMapper;
	vtkSmartPointer<vtkActor> SphereActor;

 //  img reading and contour->3d
	vtkSmartPointer<vtkPolyDataMapper> VolumeMapper;
	vtkSmartPointer<vtkActor> VolumeActor;
	TraceObject* tobj;
	vtkSmartPointer<vtkVolume> Volume;
	vtkSmartPointer<vtkContourFilter> ContourFilter;

  //raycast
	vtkSmartPointer<vtkPolyData> poly;
	vtkSmartPointer<vtkPolyDataMapper> polymap;

  //VTK widgets
  vtkSliderWidget *OpacitySlider;
  vtkSliderWidget *BrightnessSlider;
};
#endif
