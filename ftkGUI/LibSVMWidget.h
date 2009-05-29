#ifndef LIBSVMWIDGET_H
#define LIBSVMWIDGET_H

#include <PatternAnalysis/libsvm/svm.h>

#include <iostream>
#include <vector>

//#include <QtGui/QAction>
#include <QtGui/QWidget>
#include <QtGui/QAction>
//#include <QtGui/QGridLayout>
//#include <QtGui/QComboBox>
#include <QtGui/QLabel>
#include <QtGui/QPushButton>
//#include <QtGui/QCloseEvent>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QGroupBox>
#include <QtGui/QVBoxLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QScrollArea>
#include <QtGui/QDoubleSpinBox>
#include <QtCore/QAbstractItemModel>
#include <QtCore/QList>


//class PlotWindow : public QWidget
class LibSVMWidget : public QWidget
{
    Q_OBJECT;

public:
	LibSVMWidget(QAbstractItemModel *mod, QWidget *parent = 0); 

private slots:
	void go();
	void selectNone();
	void selectAll();
    
private:

	QGroupBox * initFeatureBox();
	QGroupBox * initOptionBox();

	QButtonGroup *featureGroup;

	QCheckBox *normalizeBox;
	QDoubleSpinBox *nuSpin;

	QPushButton *goButton;

	QAbstractItemModel *model;

	int columnForSVM;

 };

#endif
