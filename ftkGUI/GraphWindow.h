#ifndef GRAPHWINDOW_H
#define GRAPHWINDOW_H
 
#include <QVTKWidget.h>
#include <QtGui/QAction>
#include <QtGui/QMainWindow>
#include <QtGui/QApplication>
#include <QtGui/QDesktopWidget>
#include <QtGui/QWidget>
#include <QtGui/QStatusBar>
#include <QtGui/QMenuBar>
#include <QtGui/QPushButton>
#include <QtGui/QComboBox>
#include <QtGui/QItemSelection>
#include <QtGui/QItemSelectionModel>
#include <QtGui/QTableView>
#include <QtGui/QVBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QCloseEvent>
#include <QtGui/QDialog>
#include <QtGui/QGroupBox>
#include <QtGui/QLabel>
#include <QtGui/QCheckBox>
#include <QtGui/QButtonGroup>
#include <QtGui/QScrollArea>
#include <QtGui/QScrollBar>
#include <QtGui/QDoubleSpinBox>

#include <QtCore/QMap>
#include <QtCore/QSignalMapper>


#include "vtkTable.h"
#include <vtkTableToGraph.h>
#include <vtkViewTheme.h>
#include <vtkStringToCategory.h>
#include <vtkGraphLayout.h>
#include <vtkGraphLayoutView.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkGraphToGlyphs.h>
#include <vtkRenderer.h>
#include <vtkFast2DLayoutStrategy.h>
#include <vtkArcParallelEdgeStrategy.h>
#include <vtkPolyDataMapper.h>
#include <vtkEdgeLayout.h>
#include <vtkGraphToPolyData.h>
#include <vtkActor.h>
#include <vtkProperty.h>

#include <vtkAbstractArray.h>
#include "vtkSmartPointer.h"
#include "vtkDoubleArray.h"
#include "vtkAbstractArray.h"
#include "vtkVariantArray.h"
#include <vtkCallbackCommand.h>
#include <vtkViewUpdater.h>
#include "ObjectSelection.h"
#include <vtklookupTable.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <QApplication>
#include <QFileDialog>
#include <QFile>
#include <QCoreApplication>
#include <QTextStream>

class GraphWindow : public QMainWindow
{
    Q_OBJECT;

public:
	GraphWindow(QWidget * parent = 0);
	~GraphWindow();
	void setQtModels(QItemSelectionModel *mod);
	void setModels(vtkSmartPointer<vtkTable> table, ObjectSelection * sels = NULL);
	void SetGraphTable(vtkSmartPointer<vtkTable> table);
	void SetGraphTable(vtkSmartPointer<vtkTable> table, std::string ID1, std::string ID2);
	void SetGraphTable(vtkSmartPointer<vtkTable> table, std::string ID1, std::string ID2, std::string edgeLabel);
	void ShowGraphWindow();
	void SetSelectedIds(std::set<long int>& IDs);
	ObjectSelection * GetSelection();
	void UpdataLookupTable( std::set<long int>& IDs);

protected slots:
	static void SelectionCallbackFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData );
	void UpdateGraphView();

signals:
	void selection_Changed();

private:
	vtkSmartPointer<vtkTable> dataTable;
	std::set<long int> rootID;
	ObjectSelection * selection;

	QVTKWidget mainQTRenderWidget;
	vtkSmartPointer<vtkViewTheme> theme;
	vtkSmartPointer<vtkTableToGraph> TTG;	
	vtkSmartPointer<vtkGraphLayoutView> view;
	//SelectionAdapter * selAdapter;
	vtkSmartPointer<vtkCallbackCommand> selectionCallback;
	unsigned long observerTag;
	vtkSmartPointer<vtkLookupTable> lookupTable;
};

#endif