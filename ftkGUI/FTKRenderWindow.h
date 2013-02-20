
#ifndef RENDERWINDOW_H
#define RENDERWINDOW_H

#include <QtGui/QAction>
#include <QtGui/QMainWindow>
#include <QtGui/QWidget>
#include <QtGui/QFileDialog>
#include <QtGui/QLabel>
#include <QtGui/QComboBox>
#include <QtGui/QPushButton>
#include <QtGui/QDialog>
#include <QtGui/QStatusBar>
#include <QtGui/QVBoxLayout>
#include <QtGui/QCloseEvent>
#include <QtGui/QMenuBar>
#include <QtGui/QMenu>
#include <QtGui/QActionGroup>
#include <QVTKWidget.h>

#include <vtkSmartPointer.h>
#include <vtkTable.h>
#include <vtkVariant.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#include <vtkSphereSource.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>

#include "vnl/vnl_vector.h"

#include "ObjectSelection.h"

#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#if defined __GNUC__ || defined __APPLE__
	#include <ext/hash_map>
#else
	#include <hash_map>
#endif

class DataNetwork;
struct DataNode;
struct DataLink;

#define NUM_COLORS 128
typedef struct ren_rgb
{
	ren_rgb( double rv, double gv, double bv)
	{
		r = rv;
		g = gv;
		b = bv;
	};
	double r;
	double g;
	double b;
}r_g_b;


const static r_g_b REN_COLOR_MAP[NUM_COLORS] = {r_g_b(0.0000, 0.0000, 0.5313),r_g_b(0.0000, 0.0000, 0.5625),r_g_b(0.0000, 0.0000, 0.5938),
r_g_b(0.0000, 0.0000, 0.6250),r_g_b(0.0000, 0.0000, 0.6563),r_g_b(0.0000, 0.0000, 0.6875),r_g_b(0.0000, 0.0000, 0.7188),
r_g_b(0.0000, 0.0000, 0.7500),r_g_b(0.0000, 0.0000, 0.7813),r_g_b(0.0000, 0.0000, 0.8125),r_g_b(0.0000, 0.0000, 0.8438),
r_g_b(0.0000, 0.0000, 0.8750),r_g_b(0.0000, 0.0000, 0.9063),r_g_b(0.0000, 0.0000, 0.9375),r_g_b(0.0000, 0.0000, 0.9688),
r_g_b(0.0000, 0.0000, 1.0000),r_g_b(0.0000, 0.0313, 1.0000),r_g_b(0.0000, 0.0625, 1.0000),r_g_b(0.0000, 0.0938, 1.0000),
r_g_b(0.0000, 0.1250, 1.0000),r_g_b(0.0000, 0.1563, 1.0000),r_g_b(0.0000, 0.1875, 1.0000),r_g_b(0.0000, 0.2188, 1.0000),
r_g_b(0.0000, 0.2500, 1.0000),r_g_b(0.0000, 0.2813, 1.0000),r_g_b(0.0000, 0.3125, 1.0000),r_g_b(0.0000, 0.3438, 1.0000),
r_g_b(0.0000, 0.3750, 1.0000),r_g_b(0.0000, 0.4063, 1.0000),r_g_b(0.0000, 0.4375, 1.0000),r_g_b(0.0000, 0.4688, 1.0000),
r_g_b(0.0000, 0.5000, 1.0000),r_g_b(0.0000, 0.5313, 1.0000),r_g_b(0.0000, 0.5625, 1.0000),r_g_b(0.0000, 0.5938, 1.0000),
r_g_b(0.0000, 0.6250, 1.0000),r_g_b(0.0000, 0.6563, 1.0000),r_g_b(0.0000, 0.6875, 1.0000),r_g_b(0.0000, 0.7188, 1.0000),
r_g_b(0.0000, 0.7500, 1.0000),r_g_b(0.0000, 0.7813, 1.0000),r_g_b(0.0000, 0.8125, 1.0000),r_g_b(0.0000, 0.8438, 1.0000),
r_g_b(0.0000, 0.8750, 1.0000),r_g_b(0.0000, 0.9063, 1.0000),r_g_b(0.0000, 0.9375, 1.0000),r_g_b(0.0000, 0.9688, 1.0000),
r_g_b(0.0000, 1.0000, 1.0000),r_g_b(0.0313, 1.0000, 0.9688),r_g_b(0.0625, 1.0000, 0.9375),r_g_b(0.0938, 1.0000, 0.9063),
r_g_b(0.1250, 1.0000, 0.8750),r_g_b(0.1563, 1.0000, 0.8438),r_g_b(0.1875, 1.0000, 0.8125),r_g_b(0.2188, 1.0000, 0.7813),
r_g_b(0.2500, 1.0000, 0.7500),r_g_b(0.2813, 1.0000, 0.7188),r_g_b(0.3125, 1.0000, 0.6875),r_g_b(0.3438, 1.0000, 0.6563),
r_g_b(0.3750, 1.0000, 0.6250),r_g_b(0.4063, 1.0000, 0.5938),r_g_b(0.4375, 1.0000, 0.5625),r_g_b(0.4688, 1.0000, 0.5313),
r_g_b(0.5000, 1.0000, 0.5000),r_g_b(0.5313, 1.0000, 0.4688),r_g_b(0.5625, 1.0000, 0.4375),r_g_b(0.5938, 1.0000, 0.4063),
r_g_b(0.6250, 1.0000, 0.3750),r_g_b(0.6563, 1.0000, 0.3438),r_g_b(0.6875, 1.0000, 0.3125),r_g_b(0.7188, 1.0000, 0.2813),
r_g_b(0.7500, 1.0000, 0.2500),r_g_b(0.7813, 1.0000, 0.2188),r_g_b(0.8125, 1.0000, 0.1875),r_g_b(0.8438, 1.0000, 0.1563),
r_g_b(0.8750, 1.0000, 0.1250),r_g_b(0.9063, 1.0000, 0.0938),r_g_b(0.9375, 1.0000, 0.0625),r_g_b(0.9688, 1.0000, 0.0313),
r_g_b(1.0000, 1.0000, 0.0000),r_g_b(1.0000, 0.9688, 0.0000),r_g_b(1.0000, 0.9375, 0.0000),r_g_b(1.0000, 0.9063, 0.0000),
r_g_b(1.0000, 0.8750, 0.0000),r_g_b(1.0000, 0.8438, 0.0000),r_g_b(1.0000, 0.8125, 0.0000),r_g_b(1.0000, 0.7813, 0.0000),
r_g_b(1.0000, 0.7500, 0.0000),r_g_b(1.0000, 0.7188, 0.0000),r_g_b(1.0000, 0.6875, 0.0000),r_g_b(1.0000, 0.6563, 0.0000),
r_g_b(1.0000, 0.6250, 0.0000),r_g_b(1.0000, 0.5938, 0.0000),r_g_b(1.0000, 0.5625, 0.0000),r_g_b(1.0000, 0.5313, 0.0000),
r_g_b(1.0000, 0.5000, 0.0000),r_g_b(1.0000, 0.4688, 0.0000),r_g_b(1.0000, 0.4375, 0.0000),r_g_b(1.0000, 0.4063, 0.0000),
r_g_b(1.0000, 0.3750, 0.0000),r_g_b(1.0000, 0.3438, 0.0000),r_g_b(1.0000, 0.3125, 0.0000),r_g_b(1.0000, 0.2813, 0.0000),
r_g_b(1.0000, 0.2500, 0.0000),r_g_b(1.0000, 0.2188, 0.0000),r_g_b(1.0000, 0.1875, 0.0000),r_g_b(1.0000, 0.1563, 0.0000),
r_g_b(1.0000, 0.1250, 0.0000),r_g_b(1.0000, 0.0938, 0.0000),r_g_b(1.0000, 0.0625, 0.0000),r_g_b(1.0000, 0.0313, 0.0000),
r_g_b(1.0000, 0.0000, 0.0000),r_g_b(0.9688, 0.0000, 0.0000),r_g_b(0.9375, 0.0000, 0.0000),r_g_b(0.9063, 0.0000, 0.0000),
r_g_b(0.8750, 0.0000, 0.0000),r_g_b(0.8438, 0.0000, 0.0000),r_g_b(0.8125, 0.0000, 0.0000),r_g_b(0.7813, 0.0000, 0.0000),
r_g_b(0.7500, 0.0000, 0.0000),r_g_b(0.7188, 0.0000, 0.0000),r_g_b(0.6875, 0.0000, 0.0000),r_g_b(0.6563, 0.0000, 0.0000),
r_g_b(0.6250, 0.0000, 0.0000),r_g_b(0.5938, 0.0000, 0.0000),r_g_b(0.5625, 0.0000, 0.0000),r_g_b(0.5313, 0.0000, 0.0000),
r_g_b(0.5000, 0.0000, 0.0000)};


class FTKRenderWindow : public QMainWindow
{
    Q_OBJECT;

public:
	FTKRenderWindow(QWidget *parent = 0);
	~FTKRenderWindow();
	void setModels(vtkSmartPointer<vtkTable> tbl, ObjectSelection * sels = NULL);

public slots:
	void update();

signals:
	void closing(QWidget *widget);

protected:
	void closeEvent(QCloseEvent *event);

private slots:
	void setupUI(void);	//for initial setup
	void selectFeature(QAction *action);
	void selectColorCode(QAction *action);
	void setSphereRadius();
	void createDiscreteColorMap(void);
    
private:
	void setupColorCodeMenu();
	void updateFeatureMenu();
	void updateRenderView(void);
	r_g_b GetRGBValue(double val);
	//void BuildNetwork();

	QMenu *featureMenu;
	QMenu *colorMenu;	
	QMenu *settingMenu;
	QVTKWidget *QVTK;
	QAction *radiusAction;
	DataNetwork *data_network;

	vtkSmartPointer<vtkTable> table;
	ObjectSelection * selection;
	int selectedFeature;
	int selectedColorCode;
	double radius;
	std::map< int, std::string > classColorMap;
	std::map< std::string, std::vector< double > > discreteColorMap;

	vtkSmartPointer<vtkRenderer> Renderer;
	vtkSmartPointer<vtkRenderWindowInteractor> Interactor;		

};



//############################################################################################
//############################################################################################
// A DIALOG TO GET COLOR MAPPING FOR CLASSIFICATION RESULTS
//############################################################################################

class DiscreteColorDialog : public QDialog
{
	Q_OBJECT

public:
	DiscreteColorDialog(int max_class, QVector<QString> colors, QWidget *parent = 0);
	std::map< int, std::string > GetClassColorMap(void);

private:
	QVector<QString> colors;
	std::vector< QComboBox* > classCombos;
	QPushButton * okButton;
	QHBoxLayout * bLayout;
	QVBoxLayout * layout;
};


//############################################################################################
//############################################################################################
// DATA STRUCTURES USED TO RENDER THE NETWORK
//############################################################################################

//struct DataNode
//{
//	int x,y,z;
//	int time;
//	int id;
//	std::vector< unsigned char > color;
//	std::vector< unsigned int > neighbors;	
//};
//
//
//struct DataLink
//{
//	int src_x, src_y, src_z;
//	int trg_x, trg_y, trg_z;
//};


//class DataNetwork
//{
//public:
//	DataNetwork();
//	~DataNetwork();
//
//	//void CreatePolyDataRecursive(TraceLine* , vtkSmartPointer<vtkUnsignedCharArray> , vtkSmartPointer<vtkPoints> ,vtkSmartPointer<vtkCellArray>);
//	vtkSmartPointer<vtkPolyData> GetVTKPolyData();
//	//std::vector<TraceBit> CollectTraceBits();
//	std::vector< DataNode > GetDataNodesPointer(){ return data_nodes; };
//	//std::vector< DataLink >* GetDataLinksPointer(){ return &data_links; };
//	stdext::hash_map<unsigned int, unsigned long long int> hash_node;
//	stdext::hash_map<unsigned int, unsigned long long int> hash_link;
//	vtkSmartPointer<vtkUnsignedCharArray> OriginalColors; // save the colors;
//	
//private:
//
//	std::vector<DataNode> data_nodes;	
//	std::vector<DataLink> data_links;	
//
//};

#endif
