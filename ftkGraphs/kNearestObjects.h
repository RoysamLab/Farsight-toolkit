#ifndef _kNearestObjects_H_
#define _kNearestObjects_H_

//ITK includes
#include "itkVector.h"
#include "itkListSample.h"
#include "itkKdTree.h"
#include "itkKdTreeGenerator.h"
#include "itkWeightedCentroidKdTreeGenerator.h"
#include "itkEuclideanDistanceMetric.h"

//VTK includes
#include <vtkTable.h>
#include <vtkSmartPointer.h>
#include <vtkStringArray.h>
#include <vtkVariantArray.h>

//BOOST includes
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <utility>
#include <fstream>

using namespace boost;

class kNearestObjects
{
public: 
	typedef itk::Vector< float, 3 > MeasurementVectorType;
	typedef itk::Statistics::ListSample< MeasurementVectorType > SampleType;
	typedef itk::Statistics::KdTreeGenerator< SampleType > TreeGeneratorType;
	typedef itk::Statistics::KdTreeNode< MeasurementVectorType > KdTreeNodeType;
	typedef TreeGeneratorType::KdTreeType TreeType;
	typedef itk::Statistics::EuclideanDistanceMetric< MeasurementVectorType > DistanceMetricType;
	typedef property<vertex_name_t, std::string > VertexProperties;
	typedef adjacency_list <vecS, vecS, undirectedS,VertexProperties> NeighborGraph;
	typedef graph_traits<NeighborGraph>::vertex_descriptor node;
	typedef graph_traits<NeighborGraph>::edge_descriptor Edge;
	typedef property_map<NeighborGraph, vertex_name_t>::type node_name;
	graph_traits < NeighborGraph >::vertex_iterator vi, vi_end;
	graph_traits < NeighborGraph >::adjacency_iterator ai, ai_end;

	//constructor
	kNearestObjects(std::map< unsigned int, std::vector<float> > centroidMap);
	//destructor
	~kNearestObjects() {}
/* initialization functions */
	std::vector<std::vector< std::pair<unsigned int, double> > > k_nearest_neighbors_All(unsigned int k, unsigned short Class_dest, unsigned short Class_src);
	std::vector< std::vector< std::pair<unsigned int, double> > > k_nearest_neighbors_IDs(std::vector<unsigned int> IDs, unsigned int k, unsigned short Class_dest);
	std::vector< std::pair<unsigned int, double> > k_nearest_neighbors_ID(unsigned int id, unsigned int k, unsigned short Class_dest);
	
	std::vector<std::vector< std::pair<unsigned int, double> > > neighborsWithinRadius_All(double radius, unsigned short Class_dest, unsigned short Class_src);
	std::vector< std::vector< std::pair<unsigned int, double> > > neighborsWithinRadius_IDs(std::vector<unsigned int> IDs, double radius, unsigned short Class_dest);
	std::vector< std::pair<unsigned int, double> > neighborsWithinRadius_ID(unsigned int id, double radius, unsigned short Class_dest);
/* conversions vrom the vectors of pairs to useable formats*/	
	vtkSmartPointer<vtkTable> vectorsToGraphTable(std::vector< std::vector< std::pair<unsigned int, double> > > NeighborIDs);
	vtkSmartPointer<vtkTable> vectorsToGraphTable(std::vector< std::pair<unsigned int, double> > NeighborIds);
	
	NeighborGraph getNeighborGraph(std::vector< std::pair<unsigned int, double> >);
	NeighborGraph getNeighborGraph(std::vector<std::vector< std::pair<unsigned int, double> > >);
	
	vtkSmartPointer<vtkTable> graphToTable(NeighborGraph g);
	
	void setFeatureTable(vtkSmartPointer<vtkTable> table){featureTable = table; return;};

private:
	
	int allIds;
	SampleType::Pointer sample;
	TreeGeneratorType::Pointer treeGenerator;
	std::map< unsigned int, std::vector<float> > centerMap;
	std::map< unsigned int, std::vector<float> >::iterator It;
	std::map<int, MeasurementVectorType> idToCentroidMap;
	std::map<int, MeasurementVectorType>::iterator IdIt;
	TreeType::Pointer tree;
	NeighborGraph NG;
	vtkSmartPointer<vtkTable> graphtable;
	node_name nodeName;
	vtkSmartPointer<vtkTable> featureTable;
	int check;

	std::string convert2string(unsigned int id);
	int GetNodeIndex(unsigned int id);
	
	
};

#endif
