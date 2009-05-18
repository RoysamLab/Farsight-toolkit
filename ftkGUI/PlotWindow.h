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
	//PlotWindow(QWidget *parent = 0);
	//PlotWindow(SegmentationModel *rModel, QWidget *parent = 0);
	PlotWindow(QItemSelectionModel *mod, QWidget *parent = 0);
	void updateCombos(void); 

signals:
	void closing(QWidget *widget);

protected:
	void closeEvent(QCloseEvent *event);
	//void keyPressEvent(QKeyEvent *event);

private slots:
	void comboXChange(int c);
	void comboYChange(int c);
	void updateColumnForColor();
    
private:
	ScatterView *scatter;

	QHBoxLayout *hlayout;
	QVBoxLayout *vlayout;
	QLabel *xlabel;
	QLabel *ylabel;
	QComboBox *comboX;
	QComboBox *comboY;

	void setupUI(void);	//for initial setup

 };

#endif
