/* 
 * Copyright 2009 Rensselaer Polytechnic Institute
 * This program is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
 * for more details.
 * 
 * You should have received a copy of the GNU General Public License along 
 * with this program; if not, write to the Free Software Foundation, Inc., 
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

/*=========================================================================

  Program:   Farsight Biological Image Segmentation and Visualization Toolkit
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision: 0.00 $

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

#include "ProjectFilenamesDialog.h"
#include "ftkProjectProcessor.h"
#include "ftkProjectFiles.h"

//Farsight Includes:
#include "NuclearSegmentation/ftkNuclearSegmentation.h"
#include "NuclearSegmentation/CytoplasmSegmentation/CytoplasmSegmentation.h"
#include "NuclearSegmentation/Nuclear_Association/ftkNuclearAssociationRules.h"
#include "ftkCommon/ftkLabelImageToFeatures.h"
#include "ftkCommon/ftkUtils.h"
#include "ftkGUI/PatternAnalysisWizard.h"
#include "ftkGUI/TableWindow.h"
#include "ftkGUI/PlotWindow.h"
#include "ftkGUI/ImageBrowser5D.h"
#include "ftkGUI/HistoWindow.h"
#include "ftkGUI/LabelImageViewQT.h"

//VTK includes:
#include "vtkQtTableView.h"

//ITK Preprocessing includes:
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkSigmoidImageFilter.h"
#include "itkGrayscaleErodeImageFilter.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h" 
#include "itkGrayscaleMorphologicalOpeningImageFilter.h"
#include "itkGrayscaleMorphologicalClosingImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkAffineTransform.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkCurvatureAnisotropicDiffusionImageFilter.h"

class ParamsFileDialog;
class MarginDialog;
class ProcessThread;
class PreprocessParamsDialog;

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

	//Loading:
	void askLoadImage(void);
	void loadImage(QString fileName);
	void askLoadResult(void);
	void loadResult(QString fileName);
	void askLoadTable(void);
	void loadTable(QString fileName);
	void loadProject(void);
	//Processing:
	void processProject(void);
	void startProcess(void);
	void process(void);
	void abortProcess(void);
	void continueProcess(void);
	void deleteProcess(void);
	//Saving:
	bool saveProject(void);
	bool askSaveImage(void);
	bool saveImage(void);
	bool askSaveResult(void);
	bool saveResult(void);
	bool askSaveTable(void);
	bool saveTable(void);
	void createDefaultLogName(void);
	
	//Views:
	void toggleBounds();
	void CreateNewPlotWindow();
	void CreateNewTableWindow();
	void CreateNewHistoWindow();
	void updateViews();
	void viewClosing(QWidget * view);
	void closeViews();

	//For Editing Menu
	void setEditsEnabled(bool val);
	void clearSelections(void);
	void addCell(int x1, int y1, int x2, int y2, int z);//Once box is drawn we call this to add a cell
	void mergeCells(void);
	void deleteCells(void);
	void splitCellAlongZ(void);		//Split single cell along current z
	void splitCells(void);
	void splitCell(int x1, int y1, int z1, int x2, int y2, int z2);
	void applyExclusionMargin(void);
	void changeClass(void);
	void markVisited(void);

	//***************************************************
	// Preprocessing Menu
	void setPreprocessingEnabled(bool val);
	void AnisotropicDiffusion(void);
	void MedianFilter(void);
	void SigmoidFilter(void);
	void GrayscaleErode(void);
	void GrayscaleDilate(void);
	void GrayscaleOpen(void);
	void GrayscaleClose(void);
	void CurvAnisotropicDiffusion(void);
	//void Resample(void);
	//*****************************************************

	//For Tools menu
	void segmentNuclei(void);
	void startEditing(void);
	void stopEditing(void);
	void startSVM();
	void startKPLS();

	void about(void);

	void menusEnabled(bool val);

signals:
    
private:
	void createMenus();
	void createStatusBar();
	void createProcessToolBar();
	void createPreprocessingMenu();

	int requestChannel(ftk::Image::Pointer img);	//Request a channel from this image
	
	LabelImageViewQT *segView;
	std::vector<PlotWindow *> pltWin;
	std::vector<TableWindow *> tblWin;
	std::vector<HistoWindow *> hisWin;
	PatternAnalysisWizard *pWizard;

	QMenu *fileMenu;
	QAction *loadImageAction;
	QAction *loadLabelAction;
	QAction *loadTableAction;
	QAction *loadProjectAction;
	QAction *processProjectAction;
	QAction *saveProjectAction;
	QAction *saveImageAction;
	QAction *saveLabelAction;
	QAction *saveTableAction;
	QAction *saveDisplayAction;
	QAction *exitAction;

	QMenu *viewMenu;
	QAction *showBoundsAction;
	QAction *newTableAction;
	QAction *newScatterAction;
	QAction *newHistoAction;
	QAction *imageIntensityAction;

	QMenu *toolMenu;
	QAction *segmentNucleiAction;
	QAction *editNucleiAction;
	QAction *svmAction;		//Start the One-Class SVM outlier detecter
	QAction *kplsAction;	//Start the KPLS Classifier
	
	//For Editing Menu
	QMenu *editMenu;
	QAction *clearSelectAction;
	QAction *addAction;
	QAction *mergeAction;
	QAction *deleteAction;
	QAction *splitZAction;			//for split along z direction
	QAction *splitAction;			//for split along x-y direction
	QAction *exclusionAction;
	QAction *classAction;
	QAction *visitAction;			//Mark an object as visited

	QMenu *helpMenu;
	QAction *aboutAction;
	
	//************************************************************************
	//Preprocess menu
	QMenu *PreprocessMenu;
    QAction *AnisotropicAction;
	QAction *MedianAction;
	QAction *SigmoidAction;
	QAction *GSErodeAction;
	QAction *GSDilateAction;
	QAction *GSOpenAction;
	QAction *GSCloseAction;
	QAction *CurvAnisotropicAction;
	//QAction *ResampleAction;
	vector<double> filterParams;
	std::string FN;
	typedef unsigned char   InpPixelType;
	typedef itk::Image<unsigned char, 3> InpImageType;		
    typedef itk::Image<float, 3> FloatImageType;
	//*********************************************************************
	
	QLabel *statusLabel;						//Shown at bottom of main window
	QString lastPath;							//Last path that has been navigated to
	QString standardImageTypes;

	ftk::NuclearSegmentation *nucSeg;			//Used for editing a nuclear segmentation
	ftk::ProjectProcessor *pProc;				//My project processor
	ftk::Image::Pointer myImg;					//My currently visible image
	ftk::Image::Pointer labImg;					//Currently visible label image
	ObjectSelection * selection;				//object selection list
	vtkSmartPointer<vtkTable> table;			//table
	ftk::ProjectFiles projectFiles;				//files in the currently visible project
	ftk::ProjectDefinition projectDefinition;	//the project definition currently being used.
	
	//Processing toolbar and thread pointers:
	bool abortProcessFlag;
	bool continueProcessFlag;
	QAction * processAbort;
	QAction * processContinue;
	QLabel * processTaskLabel;
	QProgressBar * processProgress;
	QToolBar * processToolbar;
	ProcessThread *processThread;
};


class ParamsFileDialog : public QDialog
{
	Q_OBJECT
public:
	ParamsFileDialog(QString lastPth, QVector<QString> chs, QWidget *parent = 0);
	QString getFileName();
	int getChannelNumber();
private slots:
	void ParamBrowse(QString);
private:
	QLabel *channelLabel;
	QComboBox *channelCombo;
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

class ProcessThread : public QThread
{
public:
	ProcessThread(ftk::ProjectProcessor *proc = NULL);
	void run();
private:
	ftk::ProjectProcessor *myProc;
};


//***************************************************************************
//***************************************************************************
// FOR PRE-PROCESSING:
//***************************************************************************
class PreprocessParamsDialog : public QDialog
{
	Q_OBJECT
public:
	PreprocessParamsDialog(QString lastPth, QVector<QString> chs, unsigned char id,QWidget *parent = 0);
	QString getFileName();
	int getChannelNumber();
	vector<double> getParams(unsigned char id); 
	

private:
	QLabel *channelLabel;
	QComboBox *channelCombo;
	QLabel *QTParamLabel1;
	QLabel *QTParamLabel2;
	QLabel *QTParamLabel3;
	QLabel *QTParamLabel4;
	QLabel *QTParamLabel5;
	QLabel *QTParamLabel6;
	QLabel *QTParamLabel7;
	QLabel *QTParamLabel8;
	QLabel *QTParamLabel9;

	QLineEdit *QTParam1;
	QLineEdit *QTParam2;	
	QLineEdit *QTParam3;	
	QLineEdit *QTParam4;	
	QLineEdit *QTParam5;	
	QLineEdit *QTParam6;	
	QLineEdit *QTParam7;	
	QLineEdit *QTParam8;	
	QLineEdit *QTParam9;	
	
	QString lastPath;
	QPushButton *okButton;
	
	double paramVal1;
	double paramVal2;	
	double paramVal3;
	double paramVal4;
	double paramVal5;
	double paramVal6;
	double paramVal7;
	double paramVal8;
	double paramVal9;
		
	vector<double> parameters;
	
};

#endif
