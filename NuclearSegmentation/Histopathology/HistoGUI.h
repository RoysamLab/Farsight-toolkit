/* 
 * Copyright 2009 Rensselaer Polytechnic Institute
 * This program is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
 * for more details.
 * 
 * You should have received a copy of the GNU General Public License along 
 * with this program; if not, write to the Free Software Foundation, Inc., 
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */


#ifndef HistoGUI_H
#define HistoGUI_H

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

//QT Includes:
#include <QtGui/QMainWindow>
#include <QtGui/QMenuBar>
#include <QtGui/QMenu>
#include <QtGui/QMessageBox>
#include <QtGui/QFileDialog>
#include <QtGui/QToolBar>
#include <QtGui/QProgressBar>
#include <QtGui/QDialog>
#include <QtGui/QVBoxLayout>
#include <QtGui/QRadioButton>
#include <QtCore/QFileInfo>
#include <QtCore/QThread>
#include <QtCore/QString>

//Farsight Includes:
#include "ftkHistopathology.h"
#include "ftkImage/ftkImage.h"
#include "ftkGUI/PatternAnalysisWizard.h"
#include "ftkGUI/TableWindow.h"
#include "ftkGUI/PlotWindow.h"
#include "ftkGUI/ImageBrowser5D.h"
#include "ftkGUI/HistoWindow.h"
#include "ftkGUI/LabelImageViewQT.h"

class HistoGUI : public QMainWindow
{
    Q_OBJECT;

public:
	HistoGUI(QWidget * parent = 0, Qt::WindowFlags flags = 0);
	~HistoGUI();

protected:
	void closeEvent(QCloseEvent *event);

private slots:
	void setMouseStatus(int,int,int);
	void loadImage(void);
	void loadResult(void);
	void about(void);

	void menusEnabled(bool val);
	void setMenusForResult(bool val);
	void toggleBounds();
	void toggleIDs();
	void CreateNewPlotWindow();
	void CreateNewTableWindow();
	void CreateNewHistoWindow();
	void ShowHistogram();

	void updateViews();
 
signals:
    
private:
	void createMenus();
	void createStatusBar();
	
	LabelImageViewQT *segView;
	std::vector<PlotWindow *> pltWin;
	std::vector<TableWindow *> tblWin;
	HistoWindow * hisWin;

	QMenu *fileMenu;
	QAction *loadAction;
	QAction *xmlAction;
	QAction *exitAction;

	QMenu *viewMenu;
	QAction *showBoundsAction;
	QAction *showIDsAction;
	QAction *newScatterAction;
	QAction *showHistoAction;
	QAction *imageIntensityAction;

	QMenu *helpMenu;
	QAction *aboutAction;
	
	QLabel *statusLabel;			//Shown at bottom of main window

	ftk::Histopathology *patho;		//Histopathology module
	ftk::Image::Pointer myImg;		//My currently visible image
	QString myImgName;				//Name of the currently visible image
	ObjectSelection * selection;	//object selection list
	vtkSmartPointer<vtkTable> table;//table

	QString lastPath;
 };

#endif
