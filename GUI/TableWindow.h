#ifndef TABLEWINDOW_H
#define TABLEWINDOW_H

#include <QtGui/QAction>
#include <QtGui/QMainWindow>
#include <QtGui/QApplication>
#include <QtGui/QDesktopWidget>
#include <QtGui/QWidget>
#include <QtGui/QStatusBar>
#include <QtGui/QItemSelection>
#include <QtGui/QItemSelectionModel>
#include <QtGui/QTableView>
#include <QtGui/QVBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QCloseEvent>

#include <iostream>

#include "SegmentationModel.h"

class TableWindow : public QWidget
{
    Q_OBJECT;

public:
	TableWindow(QWidget *parent = 0);
	TableWindow(SegmentationModel *mod, QWidget *parent = 0);

	void SetModels(QItemSelectionModel *selectionModel);
	void ResizeToOptimalSize(void);

signals:
	void closing(QWidget *widget);

protected:
	void closeEvent(QCloseEvent *event);

public slots:
	void update();
    
private:
	QVBoxLayout *layout;
	QTableView *table;

	int visibleRows;

	void setup(void);
};

#endif
