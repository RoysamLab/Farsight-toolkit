#ifndef GRAPHWINDOWFORNEWSELECTION_H
#define GRAPHWINDOWFORNEWSELECTION_H
 
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
#include <vtkIdTypeArray.h>
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
#include "ftkGUI/SelectiveClustering.h"
#include <vtkCornerAnnotation.h>

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

class GraphWindowForNewSelection : public QMainWindow
{
    Q_OBJECT;

public:
	GraphWindowForNewSelection(QWidget * parent = 0);
	~GraphWindowForNewSelection();
	void setModels(vtkSmartPointer<vtkTable> table, SelectiveClustering * clusterSelection = NULL);
	void SetGraphTable(vtkSmartPointer<vtkTable> table);
	void SetGraphTable(vtkSmartPointer<vtkTable> table, std::string ID1, std::string ID2);
	void SetGraphTable(vtkSmartPointer<vtkTable> table, std::string ID1, std::string ID2, std::string edgeLabel, std::string xCol, std::string yCol, std::string zCol);
	void SetTreeTable(vtkSmartPointer<vtkTable> table, std::string ID1, std::string ID2, std::string edgeLabel, std::vector<vtkIdType>* clusterSize = NULL, std::vector<double> *percentVec = NULL, std::vector<double> *disVec = NULL, QString filename = "");
	void ShowGraphWindow();
	ObjectSelection * GetSelection();
	void GetProgressionTreeOrder(std::vector<long int> &order);
	void ColorTreeAccordingToFeatures(vnl_vector<double> &feature, const char *featureName);

protected:
	void UpdataLookupTable( std::vector<vtkIdType>& IDs);
	void CalculateCoordinates(vnl_matrix<vtkIdType>& adj_matrix, std::vector<Point>& pointList);
	void find(vnl_vector<vtkIdType>& vec, vtkIdType val, std::vector<vtkIdType>& equal, std::vector<vtkIdType>& nonequal);
	void getBackBones(vnl_matrix< vtkIdType>& shortest_hop, std::vector< vtkIdType>& branchnodes, std::vector< vtkIdType>& chains);
	void GetElementsIndexInMatrix(vnl_matrix<vtkIdType>& mat, vtkIdType rownum, vtkIdType max, vnl_matrix<double>& oldmat, vnl_matrix<double>& newmat, vnl_vector< int>& tag);
	void SortChainList( vnl_matrix<vtkIdType>& shortest_hop, std::vector<vtkIdType>& backbones, 
		std::vector< std::pair<vtkIdType, std::vector<vtkIdType> > >& chainList);
	bool IsExist(std::vector<vtkIdType>& vec, vtkIdType value);
	double Median( vnl_vector<double> vec);
	void GetOrder(vtkIdType node, std::vector<long int> &order);
	void UpdateCoordinatesByEdgeWeights(std::vector<Point>& oldPointList, vnl_matrix<double>& vertexList, std::vector<Point>& newPointList);
	void UpdateChainPointList(vtkIdType attachnode, std::vector<Point>& oldPointList, vnl_matrix<double>& vertexList, std::vector<Point>& newPointList);
	Point GetNewPointFromOldPoint( Point &oldPointFirst, Point &oldPointSecond, Point &newPointFirst, double weight);
	double GetEdgeWeight(vnl_matrix<double>& vertexList, long firstIndex, long secondIndex);
	virtual void closeEvent(QCloseEvent *event);
	void SetUserDefineProgression(vtkIdType nodeID);
	void SetProgressionStartTag(bool bstart);
	void UpdateProgressionPath();
	void GetProgressionPath(vnl_matrix<vtkIdType> &hopMat, vtkIdType startNode, vtkIdType endNode, std::vector< vtkIdType> &path);
	void ResetLookupTable(vtkSmartPointer<vtkLookupTable> lookuptable, double* color);
	void RestoreLookupTable();

protected slots:
	static void SelectionCallbackFunction(vtkObject* caller, long unsigned eventId, void* clientData, void* callData );
	static void HandleKeyPress(vtkObject* caller, long unsigned eventId, void* clientData, void* callData );

signals:
	void selection_Changed();

private:
	vtkSmartPointer<vtkTable> dataTable;
	SelectiveClustering * ClusterSelections;
	vtkSmartPointer<vtkIdTypeArray> clusterIds;
	std::vector< vtkIdType> size;
	vnl_vector<double> percentVector;
	vnl_vector<double> featureColorVector;
	   
	QVTKWidget mainQTRenderWidget;
	//vtkSmartPointer<vtkViewTheme> theme;
	vtkSmartPointer<vtkTableToGraph> TTG;	
	vtkSmartPointer<vtkGraphLayoutView> view;
	vtkSmartPointer<vtkCallbackCommand> selectionCallback;
	vtkSmartPointer<vtkCallbackCommand> keyPress;
	vtkSmartPointer<vtkLookupTable> lookupTable;
	vtkSmartPointer<vtkLookupTable> edgeLookupTable;
	vtkSmartPointer<vtkCornerAnnotation> cornAnnotation;

	std::map<vtkIdType, vtkIdType> indMapFromVertexToInd;
	std::vector<vtkIdType> indVectorFromIndToVertex;
	std::map< std::pair< vtkIdType, vtkIdType>, vtkIdType> edgeMapping;
	QString fileName;

	vnl_matrix< vtkIdType> shortest_hop;
	std::vector< vtkIdType> backbones;
	std::vector< std::pair<vtkIdType, std::vector<vtkIdType> > > chainList;

	// User Define progression
	vtkIdType progressionStartID;
	vtkIdType progressionEndID;
	bool bProgressionStart;
	std::vector< vtkIdType> progressionPath;
};

#endif
