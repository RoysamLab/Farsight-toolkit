#ifndef SAMPLE_EDITOR_H
#define SAMPLE_EDITOR_H

//QT Includes:
#include <QtGui/QMainWindow>
#include <QtGui/QMenuBar>
#include <QtGui/QMenu>
#include <QtGui/QMessageBox>
#include <QtGui/QFileDialog>
#include <QtGui/QStandardItemModel>
#include <QtGui/QTableView>
#include <QtCore/QFileInfo>


//Farsight Includes:
#include "ftkGUI/PlotWindow.h"
#include "ftkGUI/HistoWindow.h"

#include <vector>
#include <string>

class SampleEditor : public QMainWindow
{
    Q_OBJECT;

public:
	SampleEditor(QWidget * parent = 0, Qt::WindowFlags flags = 0);


private slots:
	void loadFile(void);

signals:
    
private:
	void createMenus();
	void createStatusBar();

	QMenu *fileMenu;
	QAction *loadAction;

	QStandardItemModel *model;
	QItemSelectionModel *selModel;

	QTableView *table;
	PlotWindow *plot;
	HistoWindow *histo;

 };


#endif