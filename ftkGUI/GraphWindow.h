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
#include <vtkCornerAnnotation.h>

#ifndef MYPOINT
#define MYPOINT
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
#endif

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
	void SetTreeTable(vtkSmartPointer<vtkTable> table, std::string ID1, std::string ID2, std::string edgeLabel, std::vector<double> *colorVec = NULL, 
					std::vector<double> *disVec = NULL, std::set<long int>* colSels = NULL, QString filename = "");
	void SetGraphTableToPassThrough(vtkSmartPointer<vtkTable> table, unsigned int nodesNumber, std::string ID1, std::string ID2, std::string edgeLabel, 
					std::vector<double> *colorVec = NULL, std::vector<double> *disVec = NULL, std::set<long int>* colSels = NULL, QString filename = "");
	void AdjustedLayout(vtkSmartPointer<vtkTable> table, std::string ID1, std::string ID2, std::string edgeLabel, std::vector<long int> *treeOrder = NULL, std::vector<double> *colorVec = NULL, std::vector<double> *disVec = NULL);
	void ShowGraphWindow();
	ObjectSelection * GetSelection();
	void GetTrendTreeOrder(std::vector<long int> &order);
	void ColorTreeAccordingToFeatures(vnl_vector<double> &feature, const char *featureName);
	static void GetTreeNodeBetweenDistance(vtkSmartPointer<vtkTable> table, std::string ID1, std::string ID2, std::string edgeLabel, vnl_matrix<double> &disMat);

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
	void SetUserDefineTrend(long int nodeID);
	void SetTrendStartTag(bool bstart);
	void UpdateTrendPath();
	void GetTrendPath(vnl_matrix<long int> &hopMat, long int startNode, long int endNode, std::vector< long int> &path);
	void ResetLookupTable(vtkSmartPointer<vtkLookupTable> lookuptable, double* color);
	void RestoreLookupTable();

	bool FindCycle(vnl_matrix<unsigned char> &adjacent_matrix, std::vector< std::list<int> > &cycleList, std::vector<int> &seperatePts);
	void SearchCycles(vnl_matrix<unsigned char> &adjacent_matrix, std::list<int> &cycleNode, int i, int preNode, std::set<int> &path, std::vector< std::list<int> >&cycleLongest);
	int MergeCyclesToSuperNodes(vnl_matrix<unsigned char> &adj_matrix, std::vector< std::list<int> > &cycleList, std::vector< std::vector<std::list<int> > >&superNodeList, vnl_matrix<int> &superNodeConnection);
	void MergeCircles(std::list<int> &circle1, std::list<int> &circle2, vnl_vector<int> &commonNode, std::vector< std::list<int> > &leftCycleList);
	void BreakCircles(std::list<int> &circle, vnl_vector<int> &commonNode, std::vector< std::list<int> > &lines);
	template<class T> vnl_vector<T> VectorAnd(vnl_vector< T > const &v, vnl_vector< T > const &u);
	template<class T> vnl_vector<T> VectorOr(vnl_vector< T > const &v, vnl_vector< T > const &u);
	void GetConnectedLines(vnl_vector<int> &circle1, std::list<int> &circle2, vnl_vector<int> &commonNode, std::vector< std::list<int> > &lineList);
    long int IsConnectedSuperNodeToNode(vnl_matrix<unsigned char> &adj_matrix, int node, std::list<int> &superNodePt);
    void CalculateCoordinatesForCircle(std::list<int> &circle, std::vector< std::list<int> > lineList, Point &center, double radius, std::vector<Point> &pointList);

	protected slots:
	static void SelectionCallbackFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData );
	static void HandleKeyPress(vtkObject* caller, long unsigned eventId, void* clientData, void* callData );
	void UpdateGraphView();

signals:
	void selection_Changed();

private:
	vtkSmartPointer<vtkTable> dataTable;
	vnl_vector<double> colorVector;
	vnl_vector<double> featureColorVector;
	ObjectSelection *selection;
	ObjectSelection *selection2;
	   
	std::set<long int> colSelectIDs;
	
	QVTKWidget mainQTRenderWidget;
	//vtkSmartPointer<vtkViewTheme> theme;
	vtkSmartPointer<vtkTableToGraph> TTG;	
	vtkSmartPointer<vtkGraphLayoutView> view;
	vtkSmartPointer<vtkCallbackCommand> selectionCallback;
	vtkSmartPointer<vtkCallbackCommand> keyPress;
	vtkSmartPointer<vtkLookupTable> lookupTable;
	vtkSmartPointer<vtkLookupTable> edgeLookupTable;
	vtkSmartPointer<vtkCornerAnnotation> cornAnnotation;
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
	bool bTrendStart;
	std::vector< long int> progressionPath;
};

#endif
