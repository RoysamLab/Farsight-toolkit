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
#include "NuclearSegmentation/ftkNuclearSegmentation.h"
#include "NuclearSegmentation/CytoplasmSegmentation/CytoplasmSegmentation.h"
#include "NuclearSegmentation/Nuclear_Association/ftkNuclearAssociationRules.h"
#include "ftkCommon/ftkLabelImageToFeatures.h"
#include "ftkCommon/ftkUtils.h"
#include "ftkImage/ftkImage.h"
#include "ftkGUI/PatternAnalysisWizard.h"
#include "ftkGUI/TableWindow.h"
#include "ftkGUI/PlotWindow.h"
#include "ftkGUI/ImageBrowser5D.h"
#include "ftkGUI/HistoWindow.h"
#include "ftkGUI/LabelImageViewQT.h"

//VTK includes:
#include "vtkQtTableView.h"

class ParamsFileDialog;
class MarginDialog;
class NucSegThread;
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
	bool saveTable(void);
	void segmentImage(void);
	void abortSegment(void);
	void stopSegment(void);
	bool checkSaveSeg(void);
	void segment(void);
	void about(void);
	void cytoSeg(void);

	void menusEnabled(bool val);
	void setMenusForResult(bool val);
	void toggleBounds();
	void toggleIDs();
	void CreateNewPlotWindow();
	void CreateNewTableWindow();
	void CreateNewHistoWindow();
	void ShowHistogram();

	//For Editing Menu
	void setEditsEnabled(bool val);
	void clearSelections(void);
	void addCell(int x1, int y1, int x2, int y2, int z);//Once box is drawn we call this to add a cell
	void mergeCells(void);
	void deleteCells(void);
	void splitCellAlongZ(void);		//Split single cell along current z
	void splitCell(int x1, int y1, int z1, int x2, int y2, int z2);
	void applyExclusionMargin(void);
	void changeClass(void);
	void markVisited(void);

	//For Tools menu
	void startAssociations();
	void startSVM();
	void startKPLS();

	void updateViews();

signals:
    
private:
	void createMenus();
	void createStatusBar();
	void createSegmentToolBar();
	void quitNucSeg();
	bool loadXMLImage(std::string filename);
	int requestChannel(ftk::Image::Pointer img);
	
	LabelImageViewQT *segView;
	std::vector<PlotWindow *> pltWin;
	std::vector<TableWindow *> tblWin;
	HistoWindow * hisWin;
	PatternAnalysisWizard *pWizard;

	QMenu *fileMenu;
	QAction *loadAction;
	QAction *saveAction;
	QAction *saveDisplayAction;
	QAction *saveTableAction;
	QAction *xmlAction;
	QAction *exitAction;

	QMenu *viewMenu;
	QAction *showBoundsAction;
	QAction *showIDsAction;
	QAction *newScatterAction;
	QAction *showHistoAction;
	QAction *imageIntensityAction;

	QMenu *helpMenu;
	QAction *aboutAction;

	QMenu *toolMenu;
	QAction *segmentAction;
	QAction *cytoAction;
	QAction *assocAction;
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
	
	QLabel *statusLabel;			//Shown at bottom of main window

	ftk::NuclearSegmentation *nucSeg;//Used for segmentation execution, loading, and editing
	ftk::Image::Pointer myImg;		//My currently visible image
	ftk::Image::Pointer labImg;		//Currently visible label image
	int nucChannel;
	int cytChannel;
	ObjectSelection * selection;	//object selection list
	vtkSmartPointer<vtkTable> table;//table

	QString lastPath;
	
	int segmentState;
	QAction * segmentAbort;
	QAction * segmentStop;
	QAction * segmentContinue;
	QLabel * segmentTaskLabel;
	QProgressBar * segmentProgress;
	QToolBar * segmentTool;

	NucSegThread *nucsegThread;
	Features * featuresThread;
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

class NucSegThread : public QThread
{
public:
	NucSegThread(ftk::NuclearSegmentation *seg = NULL);
	void run();
	int getNext(void){ return step; };
private:
	ftk::NuclearSegmentation *mySeg;
	int step;
};

class Features : public QThread
{
public:
	Features(ftk::Image::Pointer dImg, ftk::Image::Pointer lImg, int chan, std::string pre);
	void run();
	vtkSmartPointer<vtkTable> GetTable(void){ return table; };

private:
	ftk::Image::Pointer dataImg;
	ftk::Image::Pointer labImg;
	vtkSmartPointer<vtkTable> table;
	int channel;
	std::string prefix;

};

#endif
