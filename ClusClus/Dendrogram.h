#ifndef DENDROWINDOW_H
#define DENDROWINDOW_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
 
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
#include <QApplication>
#include <QFileDialog>
#include <QFile>
#include <QCoreApplication>
#include <QTextStream>

#include "vtkTable.h"
#include <vtkSelection.h>
#include <vtkLookupTable.h>
#include <vtkMutableUndirectedGraph.h>
#include <vtkRenderedGraphRepresentation.h>
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
#include <vtkLineWidget2.h>
#include <vtkLineSource.h>
#include <vtkSelectionNode.h>
#include <vtkLineRepresentation.h>
#include <vtkIdTypeArray.h>
#include <vtkAbstractArray.h>
#include <vtkAnnotationLink.h>
#include "vtkSmartPointer.h"
#include "vtkDoubleArray.h"
#include "vtkAbstractArray.h"
#include "vtkVariantArray.h"
#include <vtkIdTypeArray.h>
#include <vtkPoints.h>
#include <vtkCallbackCommand.h>
#include <vtkFloatArray.h>
#include <vtkDataSetAttributes.h>
#include <vtkBalloonWidget.h>

#include <ClusClus/clusclus.h>
#include <ftkCommon/ftkUtils.h>
#include "ObjectSelection.h"

class Dendrogram : public QMainWindow
{
    Q_OBJECT;

public:
	Dendrogram(QWidget * parent = 0);
	~Dendrogram();
	void createDataForDendogram();
	void showGraph();
	void runClusclus();
	void setModels(vtkSmartPointer<vtkTable> table = NULL, ObjectSelection * sels = NULL, double powCof = 0.2);
	void setTreeData(int num_samples, double** treedata, int* optimalleaforder);

	double** connect_Data_Tree;
	int*    Optimal_Leaf_Order;
	std::vector<std::vector< double > > Processed_Coordinate_Data_Tree;
	int num_samples;
	
protected slots:
	void SetSelectedIds(std::set<long int>& IDs );
	static void SelectionCallbackFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData );	
	void GetSelecectedIDs();

signals:
	void selection_Changed();

private:
	void reselectIds(std::set<long int>& selectedIDs, long int id);

	vtkSmartPointer<vtkTable > table;
	vtkSmartPointer<vtkMutableUndirectedGraph> graph_Layout;
	vtkSmartPointer<vtkGraphLayoutView> graphLayoutView; 
	vtkSmartPointer<vtkCallbackCommand> selectionCallback;
	vtkSmartPointer<vtkViewTheme> theme;
	vtkSmartPointer<vtkPoints> points;
	vtkSmartPointer<vtkIntArray> vertexColors;
	vtkSmartPointer<vtkLookupTable> lookupTable;
	vtkSmartPointer<vtkIdTypeArray> v;

	vtkSmartPointer<vtkBalloonRepresentation> balloonRep1;
    vtkSmartPointer<vtkBalloonWidget> balloonWidget1;
    vtkSmartPointer<vtkBalloonRepresentation> balloonRep;
    vtkSmartPointer<vtkBalloonWidget> balloonWidget ;

	QVTKWidget mainQTRenderWidget;

	ObjectSelection * Selection;
	double powCof;
	clusclus *cc1;
	clusclus *cc2;
};

#endif
