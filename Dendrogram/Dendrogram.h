#ifndef DENDROGRAM_H_
#define DENDROGRAM_H_

#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include <vtkLookupTable.h>
#include <vtkDataSetAttributes.h>
#include <vtkViewTheme.h>
#include <vtkCellPicker.h>

#include <vtkCallbackCommand.h>
#include <vtkAnnotationLink.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkSelectionNode.h>
#include <vtkIdTypeArray.h>
#include <vtkSelection.h>
#include <vtkUnsignedCharArray.h>
#include <vtkObjectFactory.h>
#include <vtkGraphWriter.h>
#include <vtkRenderedGraphRepresentation.h>
#include <vtkGlyph2D.h>
#include <vtkRegularPolygonSource.h>
#include <vtkGraphToGlyphs.h>

#include <vtkTable.h>
#include <vtkVariant.h>
#include <vtkVariantArray.h>
#include <vtkDelimitedTextReader.h>
#include <iostream>
#include <vtkTableToTreeFilter.h>
#include <vtkTreeMapView.h>
#include <vtkGroupLeafVertices.h>
#include <vtkStringToCategory.h>
#include <stdio.h>
#include <string.h>
#include <vtkGraphLayoutStrategy.h>
#include <vtkAssignCoordinatesLayoutStrategy.h>
#include <vtkSimple2DLayoutStrategy.h>
#include <vtkRandomLayoutStrategy.h>
#include <vector>
#include <time.h>
#include <string.h>
#include <vtkProperty.h>

#include <vtkLineWidget2.h>
#include <vtkLineRepresentation.h>

#include <vtkBalloonRepresentation.h>
#include <vtkBalloonWidget.h>
#include <vtkFloatArray.h>

#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkAxisActor2D.h>

#include <vtkMutableUndirectedGraph.h>
#include <vtkPoints.h>
#include <vtkGraph.h>
#include <vtkGraphLayoutView.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkLineSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
//////////////////////////////////
#include <QtGui/QMainWindow>
#include <QtGui/QWidget>
#include <QIODevice>
#include<QtGui/QWizard>
#include <QVTKWidget.h>

#include <ftkCommon/ftkUtils.h>

#include "ObjectSelection.h"

using std::vector;

typedef double *PFLOAT;
typedef PFLOAT VECTOR;
typedef PFLOAT *MATRIX;
#include "ClusClus.h"


class Dendrogram :public QMainWindow
{
	Q_OBJECT;
public :
	Dendrogram(QWidget * parent = 0, Qt::WindowFlags flags = 0);
	~Dendrogram();
	vtkSmartPointer<vtkGraphLayoutView> GetGraphLayoutView();
	vtkSmartPointer<vtkLookupTable> lookupTable;
	ObjectSelection * GetObjectSelection();
	void setSelectedIds(std::set<long int> IDs);
	int *Level_Detected;
	void setModels(vtkSmartPointer<vtkTable> tbl, ObjectSelection * sels);
	void update();
	void createDataForDendogram();
	ClusClus *Cluster;

protected slots:
	void GetSelectedLevel(int level);
	void GetSelecectedIDs();
	
signals:
	void selection_Changed(int level);

	

	
private:
	int num_row,num_col;
	int *colour_child;
	double x1;
	double x2;
	double x3;
	double y1;
	double y2;
	double y3;
	int count;
	int num_cluster;
	int clust_level;
	

	
	/////////////VTK VARIABLES///////////////
	vtkSmartPointer<vtkMutableUndirectedGraph> graph_Layout;
	vtkSmartPointer<vtkPoints> points;
	vtkSmartPointer<vtkGraphLayoutView> graphLayoutView; 
	vtkIdType *v ;	
	vtkSmartPointer<vtkTable> table;
	vtkSmartPointer<vtkCallbackCommand> selectionCallback;
	vtkSmartPointer<vtkIntArray> vertexColors;
	vtkSmartPointer<vtkViewTheme> theme;
	vtkSmartPointer<vtkBalloonRepresentation> balloonRep1;
    vtkSmartPointer<vtkBalloonWidget> balloonWidget1;
    vtkSmartPointer<vtkBalloonRepresentation> balloonRep;
    vtkSmartPointer<vtkBalloonWidget> balloonWidget ;
    
	
	vector<vector<vector<double> > > Tree3D;
	vector<vector<double> >  CharLabel;
	
 
	////////////////////////////////////////
	MATRIX Optimal_Leaf_Nodes;
	MATRIX my_data;
	MATRIX connect_Data_Tree;
	MATRIX Processed_Coordinate_Data_Tree;
	MATRIX distance_data;

	ObjectSelection * Selection;

	////////////////////functions////////////////////
	int Determine_File_Chars(char *root, int *num_data, int *numfeats) ;
	FILE *FDeclare(char *root, char *extension, char key) ;
	void VectorAllocate(VECTOR *vector, int nCols);
	void Read_Meta_Data(int num_data, int num_feats, MATRIX my_data, char *root);
	void MatrixAllocate(MATRIX *pmatrix, int nRows, int nCols,int mode);
	void AllocateCols(PFLOAT matrix[], int nRows, int nCols);
	static void SelectionCallbackFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData); ///Interactor callback
	
};
#endif





