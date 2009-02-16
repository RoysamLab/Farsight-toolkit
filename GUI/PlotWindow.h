#ifndef PLOTWINDOW_H
#define PLOTWINDOW_H

#include <QtGui/QAction>
#include <QtGui/QMainWindow>
#include <QtGui/QWidget>
#include <QtGui/QStatusBar>
#include <QtGui/QItemSelectionModel>
#include <QtGui/QVBoxLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QGridLayout>
#include <QtGui/QComboBox>
#include <QtGui/QLabel>
#include <QtGui/QPushButton>
#include <QtGui/QCloseEvent>

#include "ScatterView.h"

class PlotWindow : public QWidget
{
    Q_OBJECT;

public:
	PlotWindow(QWidget *parent = 0);
	PlotWindow(SegmentationModel *rModel, QWidget *parent = 0);
	void updateCombos(void); 

signals:
	void closing(QWidget *widget);

protected:
	void closeEvent(QCloseEvent *event);
	void keyPressEvent(QKeyEvent *event);

private slots:
	void comboXChange(int c);
	void comboYChange(int c);
	void updateColumnForColor();
    
private:
	ScatterView *scatter;

	QHBoxLayout *hlayout;
	QVBoxLayout *vlayout;
	QHBoxLayout *hlayoutT;
	QLabel *xlabel;
	QLabel *ylabel;
	QLabel *colorlabel;
	QPushButton *selectButton;
	QPushButton *clearButton;
	QComboBox *comboX;
	QComboBox *comboY;
	QComboBox *comboSelMode;

	SegmentationModel *resultModel;

	//void updateCombos(FTKItemModel *model);
	void setupSelectionModes(void);
	void setupUI(void);	//for initial setup

 };

#endif
