#include <boost/config.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <vtkGraphLayoutView.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkStringArray.h>
#include <vtkVariantArray.h>
#include <vtkTable.h>
#include <vtkTableToGraph.h>
#include <vtkSmartPointer.h>
#include <vtkViewTheme.h>
#include <vtkStringToCategory.h>
#include <vtkAbstractArray.h>
#include <ftkCommon/ftkUtils.h>
#include <vtkGraphLayout.h>
#include <vtkGraphToGlyphs.h>
#include <vtkRenderer.h>
#include <vtkFast2DLayoutStrategy.h>
#include <vtkArcParallelEdgeStrategy.h>
#include <vtkPolyDataMapper.h>
#include <vtkEdgeLayout.h>
#include <vtkGraphToPolyData.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include "ftkIntrinsicFeatures.h"

#include "ftkgnt.h"

#define VTK_CREATE(type, name) vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

class FTKgraph
{
public:
	
	bool cyto_image;
	typedef ftkgnt::RAGraph GraphType;
	typedef boost::graph_traits < GraphType >::vertex_descriptor graph_vertex;
	boost::graph_traits < GraphType >::vertex_iterator vi, vi_end;
	boost::graph_traits < GraphType >::adjacency_iterator ai, ai_end;
	typedef boost::property_map <GraphType, boost::vertex_index_t >::type index_map_type;
	typedef boost::property_map <GraphType, boost::vertex_name_t >::type name_map_type;
	typedef itk::Image<unsigned char,3> IntensityImageType;
    typedef itk::Image<unsigned short,3> LabelImageType;
	typedef IntensityImageType InputImageType;
	typedef LabelImageType OutputImageType;
	typedef ftk::LabelImageToFeatures< unsigned char,  unsigned short, 3 > FeatureCalcType;
	int flag;

	index_map_type index_map;
	name_map_type name_map;

	//: constructor
	FTKgraph();
 	//: destructor
	~FTKgraph();

	vtkSmartPointer<vtkTable> AdjacencyGraph_All(InputImageType::Pointer ip_image, OutputImageType::Pointer op_image, bool CytoImage = false);
	vtkSmartPointer<vtkTable> constructGraphTable_All(InputImageType::Pointer ip_image, OutputImageType::Pointer op_image, vtkSmartPointer<vtkTable> graphtable, bool CytoImage);
	vtkSmartPointer<vtkTable> AdjacencyGraph_ID(unsigned short id, InputImageType::Pointer ip_image, OutputImageType::Pointer op_image, bool CytoImage = false); 
	vtkSmartPointer<vtkTable> constructGraphTable_ID(unsigned short id, InputImageType::Pointer ip_image, OutputImageType::Pointer op_image, vtkSmartPointer<vtkTable> graphtable, bool CytoImage);
	vtkSmartPointer<vtkTable> BuildGraphTable( GraphType g, vtkSmartPointer<vtkTable> table);
	void DisplayGraph(vtkSmartPointer<vtkTable> graphtable);
	

};