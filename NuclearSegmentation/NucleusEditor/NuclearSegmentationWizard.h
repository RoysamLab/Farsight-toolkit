#ifndef NUCLEARSEGMENTATIONWIZARD_H
#define NUCLEARSEGMENTATIONWIZARD_H

//QT INCLUDES
#include <QtGui/QWizard>
#include <QtGui/QLabel>
//#include <QtGui/QCheckBox>
#include <QtGui/QVBoxLayout>
#include <QtGui/QComboBox>
#include <QtGui/QFileDialog>
#include <QtGui/QTextBrowser>
#include <QtCore/QFileInfo>
#include <QtCore/QFile>
#include <QtCore/QByteArray>
#include <QtCore/QTextStream>


//OTHER FARSIGHT INCLUDES
//#include "ftkImage/ftkImage.h"
#include "NuclearSegmentation/ftkNuclearSegmentation.h"
#include "SegmentationWindow.h"

class NuclearSegmentationWizard : public QWizard
{
    Q_OBJECT;

public:
	enum { Page_Input, Page_Parameters, Page_Binarize, Page_Seeds, Page_Cluster, Page_Finalize, Page_Exit };

    NuclearSegmentationWizard(QWidget *parent = 0); 

signals:

protected:
	void initializePage(int id);
	void cleanupPage(int id);
	int nextId() const;
	bool initSegmentation(void);
	
private slots:
	void executeNextStep(int whichButton);

private:
	ftk::NuclearSegmentation *seg;
	int lastStep;

};

class InputPage : public QWizardPage
{
	Q_OBJECT;

public:
	InputPage(QWidget *parent = 0);
	bool isComplete() const;
private slots:
	void ImageBrowse(QString);
private:
	QLabel *topLabel;
	QComboBox *imageFileCombo;
	SegmentationWindow *sWin;
	QString lastPath;
};

class ParametersPage : public QWizardPage
{
	Q_OBJECT;

public:
	ParametersPage(QWidget *parent = 0);
	bool isComplete() const;
private slots:
	void ParamBrowse(QString);

private:
	QLabel *paramLabel;
	QComboBox *paramFileCombo;
	QTextBrowser *fileText;
	QString lastPath;
};

class BinarizePage : public QWizardPage
{
	Q_OBJECT;
public:
	BinarizePage(QWidget *parent = 0);
	bool isComplete() const;
public slots:
	void ShowImages(ftk::Image::Pointer data, ftk::Image::Pointer label);
private:
	SegmentationWindow *sWin;
	bool hasImage;

};

class SeedsPage : public QWizardPage
{
	Q_OBJECT;
public:
	SeedsPage(QWidget *parent = 0);
	bool isComplete() const;
public slots:
	void ShowImages(ftk::Image::Pointer data, ftk::Image::Pointer label);
private:
	SegmentationWindow *sWin;
	bool hasImage;
};

class ClusterPage : public QWizardPage
{
	Q_OBJECT;
public:
	ClusterPage(QWidget *parent = 0);
	bool isComplete() const;
public slots:
	void ShowImages(ftk::Image::Pointer data, ftk::Image::Pointer label);
private:
	SegmentationWindow *sWin;
	bool hasImage;
};

class FinalizePage : public QWizardPage
{
	Q_OBJECT;
public:
	FinalizePage(QWidget *parent = 0);
	bool isComplete() const;
public slots:
	void ShowImages(ftk::Image::Pointer data, ftk::Image::Pointer label);
private:
	SegmentationWindow *sWin;
	bool hasImage;
};

class ExitPage : public QWizardPage
{
	Q_OBJECT;
public:
	ExitPage(QWidget *parent = 0){};
	//bool isComplete() const;
public slots:
private:
};

#endif