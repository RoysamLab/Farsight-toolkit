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
#include <QtGui/QMenuBar>
#include <QtGui/QMenu>
#include <QtGui/QActionGroup>

#include "ScatterView.h"

//class PlotWindow : public QWidget
class PlotWindow : public QMainWindow
{
    Q_OBJECT;

public:
	//PlotWindow(QWidget *parent = 0);
	//PlotWindow(SegmentationModel *rModel, QWidget *parent = 0);
	PlotWindow(QItemSelectionModel *mod, QWidget *parent = 0); 

signals:
	void closing(QWidget *widget);

protected:
	void closeEvent(QCloseEvent *event);
	//void keyPressEvent(QKeyEvent *event);

private slots:
	void xChange(QAction *action);
	void yChange(QAction *action);
	void colorChange(QAction *action);
    
private:
	QMenu *optionsMenu;
	QMenu *xMenu;
	QMenu *yMenu;
	QMenu *colorMenu;

	ScatterView *scatter;

	void setupUI(void);	//for initial setup
	void updateOptionMenus(void);
 };

#endif
