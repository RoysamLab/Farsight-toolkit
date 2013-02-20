#ifndef HEATMAPWINDOW_H
#define HEATMAPWINDOW_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <strstream>
#include <string>
#include <math.h>
 
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include <QVTKWidget.h>
#include <QtGui/QAction>
#include <QtGui/QMainWindow>
#include <QtGui/QApplication>
#include <QtGui/QDesktopWidget>
#include <QtGui/QWidget>
#include <QApplication>
#include <QFileDialog>
#include <QFile>
#include <QCoreApplication>
#include <QTextStream>

#include <vtkTable.h>
#include <vtkLookupTable.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>
#include <vtkViewTheme.h>
#include <vtkMutableUndirectedGraph.h>
#include <vtkGraphLayout.h>
#include <vtkGraphLayoutView.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderedGraphRepresentation.h>
#include <vtkGraphToGlyphs.h>
#include <vtkPolyDataMapper.h>
#include <vtkGraphToPolyData.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkLine.h>
#include <vtkLineWidget2.h>
#include <vtkLineSource.h>
#include <vtkLineRepresentation.h>
#include <vtkIdTypeArray.h>
#include <vtkAbstractArray.h>
#include <vtkAnnotationLink.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkVariantArray.h>
#include <vtkIdTypeArray.h>
#include <vtkPoints.h>
#include <vtkCallbackCommand.h>
#include <vtkFloatArray.h>
#include <vtkDataSetAttributes.h>
#include <vtkCellData.h>
#include <vtkPolyData.h>
#include <vtkRendererCollection.h>
#include <vtkDataSetMapper.h>
#include <vtkUnstructuredGrid.h>
#include <vtkIdTypeArray.h>
#include <vtkCommand.h>
#include <vtkCellArray.h>
#include <vtkPlaneSource.h>
#include <vtkCellPicker.h>
#include <vtkPicker.h>
#include <vtkPointPicker.h>
#include <vtkInteractorStyleImage.h>
#include <vtkExtractSelection.h>
#include <vtkObjectFactory.h>
#include <vtkStringArray.h>
#include <vtkPointSetToLabelHierarchy.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkLabelPlacementMapper.h>
#include <vtkQtLabelRenderStrategy.h>
#include <vtkPointData.h>
#include <vtkScalarBarActor.h>


#include <boost/math/distributions/normal.hpp>

#include <ftkCommon/ftkUtils.h>
#include "ClusClus/clusclus.h"
#include "ObjectSelection.h"
#include "ftkGUI/ColorMap.h"


class Heatmap : public QMainWindow
{
    Q_OBJECT;

public:
	Heatmap(QWidget * parent = 0);
	~Heatmap();
	void Initialize();
	void setDataForHeatmap(double** features, int* optimalleaforder1, int* optimalleaforder2,int num_samples, int num_features);
	void setDataForDendrograms(double** treedata1, double** treedata2 = NULL);
	void setOrders(int* optimalleaforder1,int* optimalleaforder2 = NULL);
	void setMultipleTreeData(std::vector<std::vector<std::vector<double > > > treesdata);
	void creatDataForHeatmap(double powCof);
	void setModels(vtkSmartPointer<vtkTable> table = NULL, ObjectSelection * sels = NULL, ObjectSelection * sels2 = NULL);
	void setModelsforSPD(vtkSmartPointer<vtkTable> table, ObjectSelection * sels, std::vector< int> selOrder, std::vector< int> unselOrder, std::map< int, int> *indexCluster = NULL, ObjectSelection * sels2 = NULL);
	void setModelsforSPD(vtkSmartPointer<vtkTable> table, ObjectSelection * sels, std::vector< int> sampleOrder, std::vector< int> selOrder, std::vector< int> unselOrder, std::map< int, int> *indexCluster = NULL, ObjectSelection * sels2 = NULL);
	void runClusclus();
	void runClus();
	inline void setPriority(std::vector<int> order){ priority_order = order; };
	inline void closeWindow(){ close(); };
	void runClusforSPD(std::vector< int> selOrder, std::vector< int> unselOrder);
	void runClusforSPD(std::vector< int> sampleOrder, std::vector< int> selOrder, std::vector< int> unselOrder);
	void showGraph();
	void showGraphforSPD( int selCol = 0, int unselCol = 0, bool bprogressionHeatmap = false);
	void GetSelRowCol(int &r1, int &c1, int &r2, int &c2);
	void SetSelRowCol(int r1, int c1, int r2, int c2);
	void SetInteractStyle();
	void SetSPDInteractStyle();
	void showDendrogram1();
	void showDendrogram2();	
	void showMulDendrogram1();
	void reRunClus();
	void showGraphforNe();
	void drawPointsforNe();

signals:
	void SelChanged();
	void columnToColorChanged(int value);

protected:
	virtual void closeEvent(QCloseEvent *event);

protected slots:
	void SetdenSelectedIds1(std::set<long int>& IDs, bool bfirst);
	void SetdenSelectedIdsMul(std::set<long int>& IDs);
	void GetSelecectedIDs();
	void GetSelecectedIDsForSPD();
	static void SelectionCallbackFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData );
	static void SelectionCallbackFunction1(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData );
	static void SelectionCallbackFunction2(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData );
	static void SelectionCallbackFunction3(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData );
	static void SelectionCallbackFunctionForSPD(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData );
	static void HandleKeyPress(vtkObject* caller, long unsigned eventId, void* clientData, void* callData );
	
private:
	rgb GetRGBValue(double val);
	void readmustd(double** mustd);
	void scaleData(double** mustd);
	void scaleData();
	void drawPoints1();
	void drawPoints3();
	void drawPointsForSPD();
	void drawPointsForOrderHeatmap();
	void setselectedCellIds();
	void computeselectedcells();
	void createDataForDendogram1(double powCof);
	void createDataForMulDendogram1(double powCof);
	void createDataForDendogram2(double powCof);
	void createDataForDendogram2();
	void reselectIds1(std::set<long int>& selectedIDs, long int id);
	void reselectIds2(std::set<long int>& selectedIDs2, long int id);
	void reselectIdsMul(std::set<long int>& selectedIDs, long int id);
	void reselectedwithintree(std::set<long int>& selectedIDs, long int realid, int i);
	void addDragLineforSPD(double* worldPosition);
	void selectClustersforSPD(double* worldPosition);
	void reselectClustersforSPD(std::set<long int>& selectedClusterSPD);
	void reselectIdsforSPD(long int id, std::set<long int> *clusidforSPD);
	void SetdenSelectedIdsForSPD(std::set<long int>& IDs);
	void reselectSPDIds1(std::set<long int>& selectedIDs, long int id);

public:
	int              num_samples;
	int              num_features;
	double**         mapdata;
	double**         connect_Data_Tree1;
	double**         connect_Data_Tree2;
	int*             Optimal_Leaf_Order1;
	int*             Optimal_Leaf_Order2;
	
	std::vector< std::vector< std::vector< double > > > treesdata;
	std::vector< std::vector< double > > Processed_Coordinate_Data_Tree1;
	std::vector< std::vector<double > > Processed_Coordinate_Data_Tree2;	

	ObjectSelection * Selection;
	ObjectSelection * Selection2;

private:
	QVTKWidget mainQTRenderWidget;
	vtkSmartPointer<vtkTable > table;
	vtkSmartPointer<vtkPlaneSource> aPlane;
	vtkSmartPointer<vtkFloatArray> cellData;
	vtkSmartPointer<vtkLookupTable> celllut;
	vtkSmartPointer<vtkPolyDataMapper> mapper;
	vtkSmartPointer<vtkActor> actor;
	vtkSmartPointer<vtkGraphLayoutView> view;
	vtkSmartPointer<vtkViewTheme> theme;
	vtkSmartPointer<vtkIdTypeArray>ids1;
	vtkSmartPointer<vtkIdTypeArray>ids2;

	vtkSmartPointer<vtkIdTypeArray> v;
	vtkSmartPointer<vtkMutableUndirectedGraph> graph_Layout;
	vtkSmartPointer<vtkPoints> points;
	vtkSmartPointer<vtkIntArray> vertexColors;
	vtkSmartPointer<vtkLookupTable> vetexlut;
	vtkSmartPointer<vtkCallbackCommand> selectionCallback1;
	vtkSmartPointer<vtkCallbackCommand> selectionCallback2;
	vtkSmartPointer<vtkCallbackCommand> selectionCallback3;
	vtkSmartPointer<vtkCallbackCommand> keyPress;
	vtkSmartPointer<vtkCellPicker> myCellPicker;
	vtkSmartPointer<vtkIdTypeArray> ids;

	vtkSmartPointer<vtkPoints> denpoints1;
	vtkSmartPointer<vtkCellArray> denlines1;
	vtkSmartPointer<vtkPolyData> denlinesPolyData1;
	vtkSmartPointer<vtkPolyDataMapper> denmapper1;
	vtkSmartPointer<vtkActor> denactor1;
	vtkSmartPointer<vtkUnsignedCharArray> dencolors1;

	vtkSmartPointer<vtkPoints> denpoints2;
	vtkSmartPointer<vtkCellArray> denlines2;
	vtkSmartPointer<vtkPolyData> denlinesPolyData2;
	vtkSmartPointer<vtkPolyDataMapper> denmapper2;
	vtkSmartPointer<vtkActor> denactor2;
	vtkSmartPointer<vtkUnsignedCharArray> dencolors2;

	vtkSmartPointer<vtkVariantArray> featureName;

	vtkSmartPointer<vtkLineSource> dragLineSource;
	vtkSmartPointer<vtkPolyDataMapper> dragLineMapper;
	vtkSmartPointer<vtkActor> dragLineActor;
	vtkSmartPointer<vtkScalarBarActor> scalarBar;
	vtkSmartPointer<vtkLookupTable> scalarbarLut;

	std::map<int, int> indMapFromVertexToInd;
	std::vector< std::vector<int> > indSPDMapFromIndToVertex;
	std::vector<int> indMapFromIndToVertex;
	std::vector<int> priority_order;

	int tree_num;
	std::vector<int > treespoint;

	std::map<int, int> rowMapFromOriginalToReorder;
	std::map<int, int> columnMapFromOriginalToReorder;
	std::map<int, int> rowMapForTreeData;
	std::map<int, int> columnMapForTreeData;
	std::set<long int> reselectedClusterSPD;   // for clusterSPD selection
	std::set<long int> interselectedIDs;
	std::set<long int> selectedFeatureIDs;

	int     r1;
	int     r2;
	int     c1;
	int     c2;
	int     removeActorflag;
	int     denResetflag1;
	int     denResetflag2;
	int     continueselectnum;
	bool	clusflag;
	bool    continueselect;
	bool    intersectionselect;
	bool	dragLineFlag;
	vtkIdType id1;
	vtkIdType id2;

	clusclus *cc1;
	clusclus *cc2;
};

#endif
