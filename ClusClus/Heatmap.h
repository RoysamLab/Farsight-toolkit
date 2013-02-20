#ifndef HEATMAP_H
#define HEATMAP_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <strstream>
#include <string>
#include <math.h>
#include <boost/math/special_functions.hpp>

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

#include <ftkCommon/ftkUtils.h>
#include "ObjectSelection.h"
#include "ftkGUI/ColorMap.h"
#include "Biclustering.h"

struct Treestructure
{
	std::vector<std::vector<int > > treeid;
	std::vector<std::vector<double > > coordinates;	
};

class BiHeatmap : public QMainWindow
{
    Q_OBJECT;

public:
	BiHeatmap(QWidget * parent = 0);
	~BiHeatmap();
	void showTree1();
	void showTree2();
	void Initialize();
	void showHeatmap();
	void WriteFile(const char *filename1);
	void setDataForTree1(Level_id levels1);
	void setDataForTree2(Level_id levels2);
	void setDataForHeatmap(std::vector<int > & order1,std::vector<int > & order2);
	void setModels(vtkSmartPointer<vtkTable> table = NULL, ObjectSelection * sels = NULL, ObjectSelection * sels2 = NULL);
	inline void closeWindow(){ close(); };

	ObjectSelection * Selection;
	ObjectSelection * Selection2;

signals:
	void lable(std::set<long int> lables);

protected:
	virtual void closeEvent(QCloseEvent *event);

private:
	void normalize();
	void drawPoints();
	void resetTree1();
	void resetTree2();
	void updataTree1();
	void updataTree2();
	void showFeatureNames();
	void SetInteractStyle();
	void creatDataForTree1();
	void creatDataForTree2();
	void creatDataForHeatmap();	
	void setselectedCellIds();
	void computeselectedcells();
	rgb GetRGBValue(double val);
	void addDragLine1(double* worldPosition);
	void addDragLine2(double* worldPosition);
	void setSelectIds(std::set<long int>& IDs);	
	void reselectIds1(std::set<long int>& selectedIDs1, long int id);
	void reselectIds2(std::set<long int>& selectedIDs2, long int id);
	
	int r1;
	int r2;
	int c1;
	int c2;
	int    num_rows;
	int    num_cols;
	int	   num_selection_area;
	bool   treedata1;
	bool   treedata2;
	bool   dragLineFlag1;
	bool   dragLineFlag2;
	bool   removeActorflag;
	bool   localselection;
	
	
	Level_id levels1;
	Level_id levels2;
	Treestructure tree1;
	Treestructure tree2;
	std::vector<int > order1;
	std::vector<int > order2;
	std::vector<std::vector<double > > data;

	std::vector<int> indMapFromIndToVertex;
	std::map<int, int> indMapFromVertexToInd;	
	std::map<int, int> rowMapFromOriginalToReorder;
	std::map<int, int> columnMapFromOriginalToReorder;

public slots:
	void GetSelecectedIDs();

protected slots:	
	static void HandleKeyPress(vtkObject* caller, long unsigned eventId, void* clientData, void* callData );
	static void SelectionCallbackFunctionLeftButton(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData );
	static void SelectionCallbackFunctionRightButtonDown(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData );
	static void SelectionCallbackFunctionRightButtonUp(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData );
		
private:
	QVTKWidget mainQTRenderWidget;
	vtkSmartPointer<vtkViewTheme> theme;
	vtkSmartPointer<vtkGraphLayoutView> view;
	vtkSmartPointer<vtkMutableUndirectedGraph> graph_Layout;

	vtkSmartPointer<vtkCellPicker> myCellPicker;
	vtkSmartPointer<vtkCallbackCommand> keyPress;
	vtkSmartPointer<vtkCallbackCommand> selectionCallbackleft;
	vtkSmartPointer<vtkCallbackCommand> selectionCallbackrightdown;
	vtkSmartPointer<vtkCallbackCommand> selectionCallbackrightup;

	vtkSmartPointer<vtkPlaneSource> plane;
	vtkSmartPointer<vtkFloatArray> cellData;
	vtkSmartPointer<vtkLookupTable> celllut;
	vtkSmartPointer<vtkPolyDataMapper> cellmapper;
	vtkSmartPointer<vtkActor> cellactor;	
	
	vtkSmartPointer<vtkPoints> points;
	vtkSmartPointer<vtkIdTypeArray> vertex;
	vtkSmartPointer<vtkIntArray> vertexColors;
	vtkSmartPointer<vtkLookupTable> vetexlut;
	vtkSmartPointer<vtkFloatArray> pointscales1;

	vtkSmartPointer<vtkActor> treeactor1;
	vtkSmartPointer<vtkPoints> treepoints1;
	vtkSmartPointer<vtkCellArray> treelines1;	
	vtkSmartPointer<vtkPolyDataMapper> treemapper1;
	vtkSmartPointer<vtkPolyData> treelinesPolyData1;
	vtkSmartPointer<vtkUnsignedCharArray> treecolors1;

	vtkSmartPointer<vtkActor> treeactor2;
	vtkSmartPointer<vtkPoints> treepoints2;
	vtkSmartPointer<vtkCellArray> treelines2;	
	vtkSmartPointer<vtkPolyDataMapper> treemapper2;
	vtkSmartPointer<vtkPolyData> treelinesPolyData2;
	vtkSmartPointer<vtkUnsignedCharArray> treecolors2;

	vtkSmartPointer<vtkActor> dragLineActor1;
	vtkSmartPointer<vtkLineSource> dragLineSource1;
	vtkSmartPointer<vtkPolyDataMapper> dragLineMapper1;
	vtkSmartPointer<vtkActor> dragLineActor2;
	vtkSmartPointer<vtkLineSource> dragLineSource2;
	vtkSmartPointer<vtkPolyDataMapper> dragLineMapper2;

	vtkSmartPointer<vtkTable > table;
	vtkSmartPointer<vtkVariantArray> featureName;
	vtkSmartPointer<vtkIdTypeArray> cellids;
	vtkIdType rightbuttonid1;
	vtkIdType rightbuttonid2;



	vtkSmartPointer<vtkIdTypeArray>ids1;
	vtkSmartPointer<vtkIdTypeArray>ids2;
	
	int     denResetflag1;
	int     denResetflag2;
	int     continueselectnum;
	bool	clusflag;
	bool    continueselect;
	bool    intersectionselect;
	

};

#endif
