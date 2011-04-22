#ifndef _kNearestObjects_H_
#define _kNearestObjects_H_

//ITK includes
#include "itkVector.h"
#include "itkListSample.h"
#include "itkKdTree.h"
#include "itkKdTreeGenerator.h"
#include "itkWeightedCentroidKdTreeGenerator.h"
#include "itkEuclideanDistance.h"

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

using namespace boost;

class kNearestObjects
{
public: 
	typedef itk::Vector< float, 3 > MeasurementVectorType;
	typedef itk::Statistics::ListSample< MeasurementVectorType > SampleType;
	typedef itk::Statistics::KdTreeGenerator< SampleType > TreeGeneratorType;
	typedef TreeGeneratorType::KdTreeType TreeType;
	typedef itk::Statistics::EuclideanDistance< MeasurementVectorType > DistanceMetricType;
	typedef property<vertex_name_t, std::string > VertexProperties;
	typedef adjacency_list <vecS, vecS, undirectedS,VertexProperties> kNeighborGraph;
	typedef graph_traits<kNeighborGraph>::vertex_descriptor node;
	typedef graph_traits<kNeighborGraph>::edge_descriptor Edge;
	typedef property_map<kNeighborGraph, vertex_name_t>::type node_name;
	graph_traits < kNeighborGraph >::vertex_iterator vi, vi_end;
	graph_traits < kNeighborGraph >::adjacency_iterator ai, ai_end;

	//constructor
	kNearestObjects(std::map< unsigned int, std::vector<float> > centroidMap);
	//destructor
	~kNearestObjects();

	std::vector<std::vector<unsigned int>> k_nearest_neighbors_All(unsigned int k=5);
	std::vector<unsigned int> k_nearest_neighbors_ID(unsigned int id, unsigned int k=5);
	kNeighborGraph kNearestGraph(std::vector<unsigned int>);
	kNeighborGraph kNearestGraph(std::vector<std::vector<unsigned int>>);
	vtkSmartPointer<vtkTable> kNeighborTable(kNeighborGraph g);

private:
	
	int allIds;
	SampleType::Pointer sample;
	TreeGeneratorType::Pointer treeGenerator;
	std::map< unsigned int, std::vector<float> > centerMap;
	std::map< unsigned int, std::vector<float> >::iterator It;
	std::map<int, MeasurementVectorType> idToCentroidMap;
	std::map<int, MeasurementVectorType>::iterator IdIt;
	TreeType::Pointer tree;
	kNeighborGraph KNG;
	node_name nodeName;
	int check;

	std::string convert2string(unsigned int id);
	int GetNodeIndex(unsigned int id);
	
	
};

#endif