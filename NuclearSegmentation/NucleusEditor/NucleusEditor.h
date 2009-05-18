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
#include <QtCore/QFileInfo>


//Farsight Includes:
//#include "SegmentationCommon/ftkSegmentationResult.h"
#include "NuclearSegmentation/ftkNuclearSegmentation.h"
#include "FTKImage/ftkImage.h"
#include "SegmentationModel.h"
#include "SegmentationWindow.h"
#include "ftkGUI/TableWindow.h"
#include "ftkGUI/PlotWindow.h"

class NucleusEditor : public QMainWindow
{
    Q_OBJECT;

public:
	NucleusEditor(QWidget * parent = 0, Qt::WindowFlags flags = 0);

protected:
	void closeEvent(QCloseEvent *event);

private slots:
	void loadResult(void);
	void saveResult(void);
	void about(void);

	//void closeWidget(QWidget *);
	void CreateNewPlotWindow();
	void CreateNewTableWindow();
	SegmentationWindow* CreateNewSegmentationWindow(void);

signals:
    
private:
	void createMenus();
	void createStatusBar();
	void clearModel();
	void newModel();

	vector<PlotWindow *> pltWin;
	vector<TableWindow *> tblWin;
	SegmentationWindow *segWin;

	ftk::Image::Pointer NewFTKImage(std::vector<std::string> filenames);
	ftk::Image::Pointer NewFTKImage(std::string filename);

	QMenu *fileMenu;
	QAction *saveAction;
	QAction *xmlAction;
	QAction *exitAction;

	QMenu *viewMenu;
	QAction *newScatterAction;
	QMenu *helpMenu;
	QAction *aboutAction;

	QLabel *statusLabel;

	ftk::NuclearSegmentation *segResult;
	SegmentationModel *currentModel;

	QString lastPath;
 };


#endif