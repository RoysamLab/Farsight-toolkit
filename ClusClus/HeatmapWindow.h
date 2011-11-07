#ifndef HEATMAPWINDOW_H
#define HEATMAPWINDOW_H

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
#include <vtkCellData.h>
#include <vtkPolyData.h>
#include <vtkPlaneSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkImageChangeInformation.h>
#include <vtkRendererCollection.h>
#include <vtkDataSetMapper.h>
#include <vtkUnstructuredGrid.h>
#include <vtkIdTypeArray.h>
#include <vtkTriangleFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkCommand.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPlaneSource.h>
#include <vtkCellPicker.h>
#include <vtkPicker.h>
#include <vtkPointPicker.h>
#include <vtkInteractorStyleImage.h>
#include <vtkProperty.h>
#include <vtkSelectionNode.h>
#include <vtkSelection.h>
#include <vtkExtractSelection.h>
#include <vtkObjectFactory.h>

#include <boost/math/distributions/normal.hpp>

#include <ftkCommon/ftkUtils.h>
#include "ObjectSelection.h"

using namespace std;

class MouseInteractorStyle;

class Heatmap : public QMainWindow
{
    Q_OBJECT;

public:
	Heatmap(QWidget * parent = 0);
	~Heatmap();
	void setDataForHeatmap(double** features, int* optimalleaforder1, int* optimalleaforder2,int num_samples, int num_features);
	void setDataForDendrograms(double** treedata1, double** treedata2);
	void creatDataForHeatmap();
	void setModels(vtkSmartPointer<vtkTable> table, ObjectSelection * sels, ObjectSelection * sels2);
	void showGraph();

	double** mapdata;
	double** connect_Data_Tree1;
	int*    Optimal_Leaf_Order1;
	vector<vector<double > > Processed_Coordinate_Data_Tree1;
	int num_samples;
	double** connect_Data_Tree2;
	int*    Optimal_Leaf_Order2;
	vector<vector<double > > Processed_Coordinate_Data_Tree2;	
	int num_features;
	friend class MouseInteractorStyle;

	ObjectSelection * Selection;
	ObjectSelection * Selection2;

protected slots:
	void SetdenSelectedIds1(std::set<long int>& IDs);
	//void SetdenSelectedIds2(std::set<long int>& IDs);
	void GetSelecectedIDs();
	//void GetSelecectedIDs2();
	static void SelectionCallbackFunction1(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData );	

private:
	QVTKWidget mainQTRenderWidget;
	vtkSmartPointer<vtkPlaneSource> aPlane;
	vtkSmartPointer<vtkFloatArray> cellData;
	vtkSmartPointer<vtkLookupTable> lookuptable;
	vtkSmartPointer<vtkPolyDataMapper> mapper;
	vtkSmartPointer<vtkActor> actor;
	vtkSmartPointer<vtkIntArray> cellColors;
	vtkSmartPointer<vtkGraphLayoutView> view;
	vtkSmartPointer<vtkViewTheme> theme;
	vtkSmartPointer<vtkIdTypeArray>ids1;
	vtkSmartPointer<vtkIdTypeArray>ids2;
	/////////////////////////////////////////////////////////
	vtkSmartPointer<vtkIdTypeArray> v1;
	vtkSmartPointer<vtkIdTypeArray> v2;
	vtkSmartPointer<vtkMutableUndirectedGraph> graph_Layout1;
	vtkSmartPointer<vtkMutableUndirectedGraph> graph_Layout2;
	vtkSmartPointer<vtkPoints> points1;
	vtkSmartPointer<vtkPoints> points2;
	vtkSmartPointer<vtkIntArray> vertexColors1;
	vtkSmartPointer<vtkIntArray> vertexColors2;
	vtkSmartPointer<vtkLookupTable> lookupTable1;
	vtkSmartPointer<vtkLookupTable> lookupTable2;
	vtkSmartPointer<vtkCallbackCommand> selectionCallback1;
	vtkSmartPointer<vtkCallbackCommand> selectionCallback2;
	
	void reselectIds1(std::set<long int>& selectedIDs, long int id);
	void reselectIds2(std::set<long int>& selectedIDs2, long int id);
	void createDataForDendogram1();
	void createDataForDendogram2();
	void showDendrogram1();
	void showDendrogram2();
	void scaleData();
	void drawPoints1();
	void drawPoints2();

};


class MouseInteractorStyle : public vtkInteractorStyleImage
{
public:
	static MouseInteractorStyle* New();
	vtkSmartPointer<vtkPolyData> Data;
	vtkIdType id1;
	vtkIdType id2;

	Heatmap *hm;

protected:
	MouseInteractorStyle();
	//~MouseInteractorStyle();

private:
	virtual void OnLeftButtonDown();
	virtual void OnLeftButtonUp();
	void computeselectedcells();
	void setselectedIds();

    vtkSmartPointer<vtkDataSetMapper> selectedMapper;
    vtkSmartPointer<vtkActor> selectedActor;
	vtkSmartPointer<vtkIdTypeArray> ids;
	vtkSmartPointer<vtkSelectionNode> selectionNode;
	vtkSmartPointer<vtkSelection> selection;
	vtkSmartPointer<vtkExtractSelection> extractSelection;
	vtkSmartPointer<vtkUnstructuredGrid> selected;

	int r1;
	int r2;
	int c1;
	int c2;
};

#endif