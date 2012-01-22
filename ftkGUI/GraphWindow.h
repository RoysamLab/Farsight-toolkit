#ifndef GRAPHWINDOW_H
#define GRAPHWINDOW_H
 
#include <QVTKWidget.h>
#include <QtGui/QMainWindow>
#include <vtkTable.h>
#include <vtkTableToGraph.h>
#include <vtkGraphLayoutView.h>
#include <vtkSmartPointer.h>
#include <vtkCallbackCommand.h>
#include <vtkLookupTable.h>
#include <vtkObject.h>
#include <vtkPoints.h>
#include "ObjectSelection.h"
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <string>
#include <set>
#include <map>
#include <vector>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include "ftkGUI/ColorMap.h"

typedef struct Point
{
	Point(){};
	Point( double ix, double iy)
	{
		x = ix;
		y = iy;
	};

	double x;
	double y;
}Point;

class GraphWindow : public QMainWindow
{
    Q_OBJECT;

public:
	GraphWindow(QWidget * parent = 0);
	~GraphWindow();
	void setModels(vtkSmartPointer<vtkTable> table, ObjectSelection * sels = NULL, std::vector<int> *indexCluster = NULL, ObjectSelection * sels2 = NULL);
	void SetGraphTable(vtkSmartPointer<vtkTable> table);
	void SetGraphTable(vtkSmartPointer<vtkTable> table, std::string ID1, std::string ID2);
	void SetGraphTable(vtkSmartPointer<vtkTable> table, std::string ID1, std::string ID2, std::string edgeLabel, std::string xCol, std::string yCol, std::string zCol);
	void SetTreeTable(vtkSmartPointer<vtkTable> table, std::string ID1, std::string ID2, std::string edgeLabel, std::set<long int>* colSels = NULL, QString filename = "");
	void ShowGraphWindow();
	ObjectSelection * GetSelection();
	void GetProgressionTreeOrder(std::vector<long int> &order);
	
protected:
	void SetSelectedIds(std::set<long int>& IDs);
	void SetSelectedIds2();
	void UpdataLookupTable( std::set<long int>& IDs);
	void CalculateCoordinates(vnl_matrix<long int>& adj_matrix, std::vector<Point>& pointList);
	void find(vnl_vector<long int>& vec, long int val, std::vector<long int>& equal, std::vector<long int>& nonequal);
	void getBackBones(vnl_matrix< long int>& shortest_hop, std::vector< long int>& branchnodes, std::vector< long int>& chains);
	void GetElementsIndexInMatrix(vnl_matrix<long int>& mat, long int rownum, long int max, vnl_matrix<double>& oldmat, vnl_matrix<double>& newmat, vnl_vector< int>& tag);
	void SortChainList( vnl_matrix<long int>& shortest_hop, std::vector<long int>& backbones, 
		std::vector< std::pair<long int, std::vector<long int> > >& chainList);
	bool IsExist(std::vector<long int>& vec, long int value);
	double Median( vnl_vector<double> vec);
	void GetOrder(long int node, std::vector<long int> &order);
	void UpdateCoordinatesByEdgeWeights(std::vector<Point>& oldPointList, vnl_matrix<double>& vertexList, std::vector<Point>& newPointList);
	void UpdateChainPointList(long int attachnode, std::vector<Point>& oldPointList, vnl_matrix<double>& vertexList, std::vector<Point>& newPointList);
	Point GetNewPointFromOldPoint( Point &oldPointFirst, Point &oldPointSecond, Point &newPointFirst, double weight);
	double GetEdgeWeight(vnl_matrix<double>& vertexList, long firstIndex, long secondIndex);
	virtual void closeEvent(QCloseEvent *event);
	void SetUserDefineProgression(long int nodeID);
	void SetProgressionStartTag(bool bstart);
	void UpdateProgressionPath();
	void GetProgressionPath(vnl_matrix<long int> &hopMat, long int startNode, long int endNode, std::vector< long int> &path);
	void ResetLookupTable(vtkSmartPointer<vtkLookupTable> lookuptable, double* color);

protected slots:
	static void SelectionCallbackFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData );
	void UpdateGraphView();

signals:
	void selection_Changed();

private:
	vtkSmartPointer<vtkTable> dataTable;
	ObjectSelection *selection;
	ObjectSelection *selection2;
	   
	std::set<long int> colSelectIDs;
	
	QVTKWidget mainQTRenderWidget;
	//vtkSmartPointer<vtkViewTheme> theme;
	vtkSmartPointer<vtkTableToGraph> TTG;	
	vtkSmartPointer<vtkGraphLayoutView> view;
	vtkSmartPointer<vtkCallbackCommand> selectionCallback;
	vtkSmartPointer<vtkLookupTable> lookupTable;
	vtkSmartPointer<vtkLookupTable> edgeLookupTable;
	unsigned long observerTag;
	vnl_vector<double> edgeWeights;
	vnl_matrix<double> vertextList;

	std::map<long int, long int> indMapFromVertexToInd;
	std::vector<long int> indMapFromIndToVertex;
	std::map<long int, long int> indMapFromVertexToClusInd;
	std::vector< std::vector<int> > indMapFromClusIndToVertex;
	std::vector< std::vector<int> > indMapFromClusIndToInd;
	std::map< std::pair< long int, long int>, long int> edgeMapping;
	QString fileName;

	vnl_matrix<long int> shortest_hop;
	std::vector< long int> backbones;
	std::vector< std::pair<long int, std::vector<long int> > > chainList;

	// User Define progression
	long int progressionStartID;
	long int progressionEndID;
	bool bProgressionStart;
	std::vector< long int> progressionPath;
};

#endif
