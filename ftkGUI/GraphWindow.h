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
	void setModels(vtkSmartPointer<vtkTable> table, ObjectSelection * sels = NULL);
	void SetGraphTable(vtkSmartPointer<vtkTable> table);
	void SetGraphTable(vtkSmartPointer<vtkTable> table, std::string ID1, std::string ID2);
	void SetGraphTable(vtkSmartPointer<vtkTable> table, std::string ID1, std::string ID2, std::string edgeLabel);
	void SetTreeTable(vtkSmartPointer<vtkTable> table, std::string ID1, std::string ID2, std::string edgeLabel, QString filename = "");
	void ShowGraphWindow();
	ObjectSelection * GetSelection();
	
protected:
	void SetSelectedIds(std::set<long int>& IDs);
	void UpdataLookupTable( std::set<long int>& IDs);
	void CalculateCoordinates(vnl_matrix<long int>& adj_matrix, std::vector<Point>& pointList);
	void find(vnl_vector<long int>& vec, long int val, std::vector<long int>& equal, std::vector<long int>& nonequal);
	void getBackBones(vnl_matrix< long int>& shortest_hop, std::vector< long int>& branchnodes, std::vector< long int>& chains);
	void GetElementsIndexInMatrix(vnl_matrix<long int>& mat, long int rownum, long int max, vnl_matrix<double>& oldmat, vnl_matrix<double>& newmat, vnl_vector< int>& tag);
	void SortChainList( vnl_matrix<long int>& shortest_hop, std::vector<long int>& backbones, 
		std::vector< std::pair<long int, std::vector<long int> > >& chainList);
	double Median( vnl_vector<double> vec);
	bool IsExist(std::vector<long int>& vec, long int value);
	
protected slots:
	static void SelectionCallbackFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData );
	void UpdateGraphView();

signals:
	void selection_Changed();

private:
	vtkSmartPointer<vtkTable> dataTable;
	ObjectSelection *selection;
	vtkSmartPointer<vtkPoints> points;   
	
	QVTKWidget mainQTRenderWidget;
	vtkSmartPointer<vtkViewTheme> theme;
	vtkSmartPointer<vtkTableToGraph> TTG;	
	vtkSmartPointer<vtkGraphLayoutView> view;
	vtkSmartPointer<vtkCallbackCommand> selectionCallback;
	vtkSmartPointer<vtkLookupTable> lookupTable;
	unsigned long observerTag;

	std::map<long int, long int> indMapFromVertexToInd;
	std::vector<long int> indMapFromIndToVertex;
	QString fileName;
};

#endif
