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

#ifndef NUCLEUS_EDITOR_H
#define NUCLEUS_EDITOR_H

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

//QT Includes:
#include <QtGui/QMainWindow>
#include <QtGui/QMenuBar>
#include <QtGui/QMenu>
#include <QtGui/QMessageBox>
#include <QtGui/QFileDialog>
#include <QtGui/QToolBar>
#include <QtGui/QProgressBar>
#include <QtGui/QDialog>
#include <QtGui/QVBoxLayout>
#include <QtGui/QRadioButton>
#include <QtCore/QFileInfo>
#include <QtCore/QThread>

//Farsight Includes:
//#include "SegmentationCommon/ftkSegmentationResult.h"
#include "NuclearSegmentation/ftkNuclearSegmentation.h"
#include "ftkImage/ftkImage.h"
#include "SegmentationModel.h"
#include "SegmentationWindow.h"
//#include "NuclearSegmentationWizard.h"
#include "ftkGUI/TableWindow.h"
#include "ftkGUI/PlotWindow.h"
#include "ftkGUI/ImageBrowser5D.h"
#include "ftkGUI/HistoWindow.h"
#include "vtkSliderWidget.h"
#include "vtkSliderRepresentation2D.h"
#include "vtkCallbackCommand.h"
#include "Seed3D.h"
//#include "SegmentationView.h"

class ParamsFileDialog;
class MarginDialog;
class Load;
class Binarize;
class SeedDetect;
class Cluster;
class Finalize;
class Features;

class NucleusEditor : public QMainWindow
{
    Q_OBJECT;

public:
	NucleusEditor(QWidget * parent = 0, Qt::WindowFlags flags = 0);
	~NucleusEditor();

protected:
	void closeEvent(QCloseEvent *event);

private slots:
	void setMouseStatus(int,int,int);
	void loadImage(void);
	void loadResult(void);
	bool saveResult(void);
	void segmentImage(void);
	void abortSegment(void);
	void stopSegment(void);
	bool checkSaveSeg(void);
	void segment(void);
	void about(void);
	//void loadDatFile(void);
	//void closeWidget(QWidget *);
	void toggleBounds();
	void toggleIDs();
	void CreateNewPlotWindow();
	void CreateNewTableWindow();
	void CreateNewHistoWindow();
	void ShowHistogram();
	void CreateNewSegWindow();
	void view3D(); 
	//DEMO - python is necessary for the demo
	void OpenPythonWindow();
	bool BrowseForPythonExecutable();

	//For Editing Menu
	void setEditsEnabled(bool val);
	void clearSelections(void);
	void mergeCells(void);
	void deleteCells(void);
	//void splitCells(void);	//Optimally split several cells - not implemented
	void splitCell(void);		//Split single cell along current z
	void addCell(void);
	void applyExclusionMargin(void);
	void changeClass(void);

	void startSplitting(void);	//begin splitting mode - user must select seeds
	void endSplitting(void);	//end splitting mode and do the splits

signals:
    
private:
	void createMenus();
	void createStatusBar();
	void createSegmentToolBar();
	void clearModel();
	void newModel();
	
	std::vector<PlotWindow *> pltWin;
	std::vector<TableWindow *> tblWin;
	HistoWindow * hisWin;
	SegmentationWindow *segWin;

	ftk::Image::Pointer NewFTKImage(std::vector<std::string> filenames);
	ftk::Image::Pointer NewFTKImage(std::string filename);

	QMenu *fileMenu;
	QAction *loadAction;
	QAction *saveAction;
	QAction *xmlAction;
	QAction *segmentAction;
	QAction *exitAction;
	//QAction *datAction;

	QMenu *viewMenu;
	QAction *showBoundsAction;
	QAction *showIDsAction;
	QAction *newScatterAction;
	QAction *showHistoAction;
	QMenu *helpMenu;
	QAction *aboutAction;
	QAction *pythonAction;
	QAction *imageIntensityAction;
	QAction *seed3DAction;


// 3D Viewer variables
	Seed3D *Seeds;
	unsigned char segFlag;
	/*QWidget *browse;
	vtkSmartPointer<vtkRenderer> Renderer;
	QVTKWidget *QVTK;
	vtkSmartPointer<vtkRenderWindowInteractor> Interactor;
	vtkSmartPointer<vtkPointPicker>PointPicker;
	vtkSmartPointer<vtkCallbackCommand> isPicked;
	vtkSliderRepresentation2D *sliderRep;
	vtkSliderRepresentation2D *sliderRep2;
	vtkSliderWidget *sliderWidget;
	vtkSliderWidget *sliderWidget2;*/
	
	//For Editing Menu
	QMenu *editMenu;
	QAction *clearSelectAction;
	QAction *mergeAction;
	QAction *deleteAction;
	QAction *addAction;
	QMenu *splitMenu;	
	QAction *splitAction;			//for split along z direction
	QAction *splitStartAction;		//for regular seed split
	QAction *splitEndAction;
	QAction *exclusionAction;
	QAction *classAction;
	
	QLabel *statusLabel;

	ftk::NuclearSegmentation *seg;
	ftk::Image::Pointer myImg;
	SegmentationModel *currentModel;

	QString lastPath;
	QString myImgName;
	int segmentState;
	QAction * segmentAbort;
	QAction * segmentStop;
	QAction * segmentContinue;
	QLabel * segmentTaskLabel;
	QProgressBar * segmentProgress;
	QToolBar * segmentTool;
	Load * loadThread;
	Binarize * binaryThread;
	SeedDetect * seedThread;
	Cluster * clusterThread;
	Finalize * finalizeThread;
	Features * featuresThread;

	//DEMO variable
	bool ConfirmClosePython();
	QSettings *settings;
	QProcess *pythonProcess;
	QLabel *pythonLabel;
	QLabel *currentPythonLabel;
	QPushButton *browseForPythonButton;


	//bool editStatus; //false shows that we cannot apply editing to cells	
 };

class ParamsFileDialog : public QDialog
{
	Q_OBJECT
public:
	ParamsFileDialog(QString lastPth, QWidget *parent = 0);
	QString getFileName();
private slots:
	void ParamBrowse(QString);
private:
	QRadioButton *autoButton;
	QRadioButton *fileButton;
	QComboBox *fileCombo;
	QString lastPath;
	QPushButton *okButton;
};

class MarginDialog : public QDialog
{
	Q_OBJECT
public:
	MarginDialog(QWidget *parent = 0);
	int getMargin();
	int getZ();
private:
	QSpinBox * marginSpin;
	QSpinBox * zSpin;
	QPushButton *okButton;
};

class Load : public QThread
{
public:
	Load(ftk::NuclearSegmentation *seg = NULL);
	void run();

private:
	ftk::NuclearSegmentation *mySeg;
};


class Binarize : public QThread
{
public:
	Binarize(ftk::NuclearSegmentation *seg = NULL);
	void run();

private:
	ftk::NuclearSegmentation *mySeg;
};

class SeedDetect : public QThread
{
public:
	SeedDetect(ftk::NuclearSegmentation *seg = NULL);
	void run();

private:
	ftk::NuclearSegmentation *mySeg;
};

class Cluster : public QThread
{
public:
	Cluster(ftk::NuclearSegmentation *seg = NULL);
	void run();

private:
	ftk::NuclearSegmentation *mySeg;
};

class Finalize : public QThread
{
public:
	Finalize(ftk::NuclearSegmentation *seg = NULL);
	void run();

private:
	ftk::NuclearSegmentation *mySeg;
};

class Features : public QThread
{
public:
	Features(ftk::NuclearSegmentation *seg = NULL);
	void run();

private:
	ftk::NuclearSegmentation *mySeg;
};

#endif
