#ifndef ACTIVE_VALIDATION_WIZARD_H
#define ACTIVE_VALIDATION_WIZARD_H

#include <QtGui/QAction>
#include <QtGui/QWidget>
#include <QtGui/QVBoxLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QGridLayout>
#include <QtGui/QDialog>
#include <QtGui/QColorDialog>
#include <QtGui/QLabel>
#include <QtGui/QPushButton>
#include <QtCore/QStringList>
#include <QtCore/QSignalMapper>
#include <QtGui/QTableWidget>
#include <QtGui/QTableWidgetItem>
#include <QtGui/QLineEdit>
#include <QtGui/QPlainTextEdit>
#include <QHeaderView>
#include <QRadioButton>
#include <QButtonGroup>
#include <QLCDNumber>
#include <QFileDialog>
#include <QDir>
#include <fstream>
#include <QMessageBox>

#include <vtkTable.h>
#include <vtkSmartPointer.h>
#include <vtkVariant.h>

#include <vector>
#include <algorithm>

#include "ActiveValidation.h"

class Active_Validation_Wizard : public QWidget
{
	Q_OBJECT

public:
	Active_Validation_Wizard(QWidget *parent=0);
	~Active_Validation_Wizard();

private:
    QLabel *dataFileNamedata;
    QPushButton *browseButtondata;
    QPushButton *loadButtondata;

    QLabel *dataFileNamelable;
    QPushButton *browseButtonlable;
    QPushButton *loadButtonlable;
	
	QLabel *featureNumLabel;
    QLabel *featureNum;
    QLabel *sampleNumLabel;
    QLabel *sampleNum;

    QLabel *numbinLable;
    QLineEdit *numbinBox;
    QLabel *deltaLable;
    QLineEdit *deltaBox;
	QLabel *runtimeLable;
    QLineEdit *runtimeBox;

    QPushButton *runButton;

	QString FileNamedata;
	QString FileNamelable;

private slots:
    void browsedata();
    void loaddata();
    void browselable();
    void loadlable();
	void Run();

private:
	void ReadFiledata(const char *filename);
	void ReadFilelable(const char *filename);
	int numsamp;
	int numfeat;
	std::vector<std::vector<double > > data;
	std::vector<int > label;
	std::vector<double > phat;
	std::vector<double > varphat;
	std::vector<int> numiteration;
	std::vector<int > numsampled;
	std::vector<double > varphattrue;

	void WriteFile();

private:
	double sumnumsam;
	double sumnumit;
};

#endif