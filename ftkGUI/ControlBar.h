#ifndef CONTROLBAR_H
#define CONTROLBAR_H

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "FarsightConfig.h"

//VTK INCLUDES:
//#include <vtkQtTableView.h>
//#include <vtkQtBarChartView.h>
//#include <vtkSelectionLink.h>
//#include <vtkDataRepresentation.h>

//QT INCLUDES
#include <QtCore/QThread>
#include <QtCore/QByteArray>
#include <QtCore/QFileInfo>
#include <QtCore/QDateTime>
#include <QtCore/QProcess>
#include <QtCore/QSettings>
#include <QtGui/QMainWindow>
#include <QtGui/QWidget>
#include <QtGui/QMenuBar>
#include <QtGui/QLabel>
#include <QtGui/QStatusBar>
#include <QtGui/QVBoxLayout>
#include <QtGui/QDockWidget>
#include <QtGui/QPushButton>
#include <QtGui/QToolButton>
#include <QtGui/QButtonGroup>
#include <QtGui/QItemSelectionModel>
#include <QtGui/QMessageBox>
#include <QtGui/QDialog>
#include <QtGui/QFileDialog>

//GUI INCLUDES
//#include "ModuleWidget.h"
#include "SegmentationModel.h"
#include "SegmentationWindow.h"
#include "TableWindow.h"
#include "PlotWindow.h"
#include "ImageBrowser5D.h"
#include "Histogram.h"
//#include "TableModel.h"
//#include "TableView.h"

//OTHER LOCAL INCLUDES
//#include "SegmentationCommon/ftkSegmentationResult.h"
#include "ftkImage/ftkImage.h"
//#include "Qt/Python/QtPythonDialog.h"

//*******************************************
// ADD NEW MODULES HERE
// FOLLOW BUILD_NUCLEI EXAMPLE
//*******************************************
#ifdef BUILD_NUCLEI
	#include "NuclearSegmentation/ftkNuclearSegmentation.h"
#endif

class ControlBar : public QMainWindow
{
    Q_OBJECT;

public:
    ControlBar(const char *c);

protected:
	void closeEvent(QCloseEvent *event);

private slots:
	void loadImage(void);
	void loadImageVTK(void);
	void loadImageSeries(void);
	void loadResult(void);
	void loadMetaResult(void);
	void loadOutliers(void);
	void saveResult(void);
	void showHistogram(void);
	void showOutliers();
	void refreshViews();
	void test(void);
	void about(void);

	void closeWidget(QWidget *);
	void EditPreferences();
	void InitializePreferencesDialog();
	void CreateNewPlotWindow();
	void CreateNewTableWindow();
    bool BrowseForPythonExecutable();
	void OpenPythonWindow();
	SegmentationWindow* CreateNewSegmentationWindow();
//  void initPythonInterpretor();

//*******************************************
// ADD NEW MODULES HERE
// FOLLOW BUILD_YOUSEF EXAMPLE
//*******************************************
#ifdef BUILD_NUCLEI
	void newYousefSeg(void);
#endif

signals:
    
private:
	void createMenus();
	void createStatusBar();
	void createEditWidget();
	void startModule();
	void clearModel();
	void newModel();
	void displayHelp();
  bool ConfirmClosePython();

	QDateTime outlierModTime;

	vector<PlotWindow *> pltWin;
	vector<TableWindow *> tblWin;
	vector<SegmentationWindow *> segWin;

	ftk::Image::Pointer NewFTKImage(std::vector<std::string> filenames);
	ftk::Image::Pointer NewFTKImage(std::string filename);
	int isLoaded(std::string filename);

  QMenu *fileMenu;
	QAction *openAction;
	QAction *openActionBeta;
	QAction *openSeriesAction;
	QAction *loadOutliersAction;
	QAction *saveAction;
	QAction *xmlAction;
	QAction *metaAction;
	QAction *testAction;
	QAction *histoAction;
  QAction *exitAction;

	QMenu *editMenu;
	QAction *preferencesAction;

	QMenu *viewMenu;
	QAction *newScatterAction;
	QAction *outlierAction;
	QAction *refreshAction;
	QAction *pythonAction;

	QMenu *helpMenu;
	QAction *aboutAction;

	QMenu *segmentMenu;

  QSettings *settings;
  QProcess *pythonProcess;

  QDialog *preferencesDialog;
  QGridLayout *preferencesLayout;
  QLabel *pythonLabel;
  QLabel *currentPythonLabel;
  QPushButton *browseForPythonButton;
  QPushButton *submitPreferencesButton;
  QPushButton *cancelPreferencesButton;
//  QtPythonDialog *PythonDialog;
  char *argv0;
  QString pythonFiles;
  QString exeFiles;


//*******************************************
// ADD NEW MODULES HERE
// FOLLOW BUILD_YOUSEF EXAMPLE
//*******************************************
#ifdef BUILD_NUCLEI
	QAction *yousefAction;
#endif

	QLabel *statusLabel;

	//ModuleWidget *module;

	//ftk::SegmentationResult *segResult;
	ftk::NuclearSegmentation *segResult;
	SegmentationModel *currentModel;

	std::vector< ftk::Image::Pointer > loadedImages;
	QString lastPath;
 };

#endif
