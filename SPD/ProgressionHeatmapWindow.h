#ifndef PROGRESSIONHEATMAPWINDOW_H
#define PROGRESSIONHEATMAPWINDOW_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <strstream>
#include <string>
#include <math.h>
 
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
#include <vtkLine.h>
#include <vtkBalloonRepresentation.h>
#include <vtkBalloonWidget.h>
#include <vtkSphereSource.h>
#include <vtkStringArray.h>
#include <vtkHoverWidget.h>
#include <vtkPointSetToLabelHierarchy.h>

#include <boost/math/distributions/normal.hpp>

#include <ftkCommon/ftkUtils.h>
#include "ClusClus/clusclus.h"
#include "ObjectSelection.h"
#include "ftkGUI/ColorMap.h"

/*class vtkHoverCallback : public vtkCommand
{
  public:
    static vtkHoverCallback *New()
    {
      return new vtkHoverCallback;
    }
 
    vtkHoverCallback() {}
 
    virtual void Execute(vtkObject*, unsigned long event, void *vtkNotUsed(calldata))
    {
      switch (event) 
        {
        case vtkCommand::TimerEvent:
          std::cout << "TimerEvent -> the mouse stopped moving and the widget hovered" << std::endl; 
          break;
        case vtkCommand::EndInteractionEvent:
          std::cout << "EndInteractionEvent -> the mouse started to move" << std::endl; 
          break;
        }
    }
};*/

class ProgressionHeatmap : public QMainWindow
{
    Q_OBJECT;

public:
	ProgressionHeatmap(QWidget * parent = 0);
	~ProgressionHeatmap();
	void setDataForHeatmap(double** features, int* optimalleaforder1, int* optimalleaforder2,int num_samples, int num_features);
	void setDataForSimilarMatrixHeatmap(double** features, int* optimalleaforder1, int* optimalleaforder2,int num_samples, int num_features);
	void setDataForDendrograms(double** treedata1, double** treedata2);
	void creatDataForHeatmap(double powCof);
	void creatDataForProgressionHeatmap(double powCof);
	void creatDataForSimilarMatrixHeatmap();
	void setModels(vtkSmartPointer<vtkTable> table = NULL, ObjectSelection * sels = NULL, ObjectSelection * sels2 = NULL);
	void runClusclus();
	void showGraph();
	void showSimilarMatrixGraph();
	void showSPDGraph();
	void GetSelRowCol(int &r1, int &c1, int &r2, int &c2);
	void SetSelRowCol(int r1, int c1, int r2, int c2);
	void SetInteractStyle();
	void showDendrogram1();
	void showDendrogram2();	

	int              num_samples;
	int              num_features;
	double**         mapdata;
	double**         connect_Data_Tree1;
	double**         connect_Data_Tree2;
	int*             Optimal_Leaf_Order1;
	int*             Optimal_Leaf_Order2;
	
	std::vector< std::vector< double > > Processed_Coordinate_Data_Tree1;
	std::vector< std::vector< double > > Processed_Coordinate_Data_Tree2;	

	ObjectSelection * Selection;
	ObjectSelection * Selection2;

protected:
	virtual void closeEvent(QCloseEvent *event);

signals:
	void SelChanged();

protected slots:
	void SetdenSelectedIds1(std::set<long int>& IDs, bool bfirst);
	void GetSelecectedIDs();
	static void SelectionCallbackFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData );
	static void SelectionCallbackFunction1(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData );
	static void SelectionCallbackFunction2(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData );
	static void SelectionCallbackFunction3(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData );

	static void HandleKeyPress(vtkObject* caller, long unsigned eventId, void* clientData, void* callData );

private:
	rgb GetRGBValue(double val)
	{
		int index = COLOR_MAP_SIZE * val - 1;   // when val = 1; index should be the max index
		if( index >= COLOR_MAP_SIZE)
		{
			index = COLOR_MAP_SIZE - 1;
		}
		else if( index < 0)
		{
			index = 0;
		}
		return COLORMAP[index];
	};

private:
	QVTKWidget mainQTRenderWidget;
	vtkSmartPointer<vtkTable > table;
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

	vtkSmartPointer<vtkIdTypeArray> v1;
	vtkSmartPointer<vtkMutableUndirectedGraph> graph_Layout1;
	vtkSmartPointer<vtkPoints> points1;
	vtkSmartPointer<vtkIntArray> vertexColors1;
	vtkSmartPointer<vtkLookupTable> lookupTable1;
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

	//vtkSmartPointer<vtkBalloonRepresentation> balloonRep;
 //   vtkSmartPointer<vtkBalloonWidget> balloonWidget;
 //   vtkSmartPointer<vtkBalloonRepresentation> balloonRep2;
 //   vtkSmartPointer<vtkBalloonWidget> balloonWidget2 ;
	//vtkSmartPointer<vtkHoverWidget> hoverWidget;
	//vtkSmartPointer<vtkHoverCallback> hoverCallback;

	void scaleData();
	void drawPoints1();
	void drawPoints2();
	void drawPoints3();
	void setselectedCellIds();
	void computeselectedcells();
	void createDataForDendogram1(double powCof);
	void createDataForDendogram2(double powCof);
	void reselectIds1(std::set<long int>& selectedIDs, long int id);
	void reselectIds2(std::set<long int>& selectedIDs2, long int id);

	//vtkSmartPointer<vtkTable> dataTable;
	std::map<int, int> indMapFromVertexToInd;
	std::vector<int> indMapFromIndToVertex;

	int     r1;
	int     r2;
	int     c1;
	int     c2;
	int     removeActorflag;
	int     denResetflag1;
	int     denResetflag2;
	int     continueselectnum;
	bool    continueselect;
	bool    intersectionselect;
	vtkIdType id1;
	vtkIdType id2;
	std::set<long int> interselectedIDs;

	clusclus *cc1;
	clusclus *cc2;
};

#endif
