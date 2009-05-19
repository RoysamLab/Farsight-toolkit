#ifndef MODULEWIDGET_H
#define MODULEWIDGET_H

#include "FarsightConfig.h"


//QT INCLUDES
#include <QtGui/QWidget>
#include <QtGui/QLabel>
#include <QtGui/QVBoxLayout>
#include <QtGui/QDockWidget>
#include <QtGui/QPushButton>
#include <QtGui/QButtonGroup>
#include <QtGui/QTextBrowser>
#include <QtGui/QComboBox>
#include <QtGui/QFileDialog>
#include <QtGui/QSpinBox>

//GUI INCLUDES
//#include "PlotWindow.h"
//#include "TableWindow.h"
#include "SegmentationWindow.h"
#include "SegmentationModel.h"

//OTHER LOCAL INCLUDES
#include "ftkImage/ftkImage.h"
#include "NuclearSegmentation/ftkNuclearSegmentation.h"

#include <sstream>


class ModuleWidget : public QWidget
{
    Q_OBJECT;

public:
    ModuleWidget(vector<string>);

signals:
	void closing(QWidget *widget);

protected:
	void closeEvent(QCloseEvent *event);

private slots:
	void ImageBrowse(QString);
	void ParamBrowse(QString);
	void Initialize(void);
	void runModule(int);
	void saveResults(void);
	void featToXML(void);
	void showResults(void);
	void closeWidget(QWidget *widget);

private:
	void SetupUi(void);
	void CreateNewSegmentationWindow();
	void closeChildren(void);

	//Classes I know of:
	ftk::NuclearSegmentation *segmentation;

	//My windows:
	SegmentationWindow *segWin;

	//Widgets shown in this widget
	QTextBrowser *msgBox;
	QComboBox *imageFileCombo;
	QComboBox *paramFileCombo;

	//Variables
	QString lastPath;
	vector<string> possibleFilenames;
};

#endif
