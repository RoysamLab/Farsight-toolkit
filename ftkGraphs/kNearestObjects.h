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

template <int num_dimensions>
class kNearestObjects
{
public: 
	typedef itk::Vector< double, num_dimensions > MeasurementVectorType;
	typedef itk::Statistics::ListSample< MeasurementVectorType > SampleType;
	typedef itk::Statistics::KdTreeGenerator< SampleType > TreeGeneratorType;
	typedef itk::Statistics::KdTreeNode< MeasurementVectorType > KdTreeNodeType;
	typedef typename TreeGeneratorType::KdTreeType TreeType;
	typedef itk::Statistics::EuclideanDistanceMetric< MeasurementVectorType > DistanceMetricType;
	typedef boost::property< boost::vertex_name_t, std::string > VertexProperties;
	typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS, VertexProperties> NeighborGraph;
	typedef boost::graph_traits<NeighborGraph>::vertex_descriptor node;
	typedef boost::graph_traits<NeighborGraph>::edge_descriptor Edge;
	typedef boost::property_map<NeighborGraph, boost::vertex_name_t>::type node_name;
	boost::graph_traits < NeighborGraph >::vertex_iterator vi, vi_end;
	boost::graph_traits < NeighborGraph >::adjacency_iterator ai, ai_end;

	//constructor
	kNearestObjects(std::map< unsigned int, std::vector<double> > centroidMap);
	//destructor
	~kNearestObjects() {}
/* initialization functions */
	std::vector<std::vector< std::pair<unsigned int, double> > > k_nearest_neighbors_All(unsigned int k, unsigned short Class_dest, unsigned short Class_src);
	std::vector< std::vector< std::pair<unsigned int, double> > > k_nearest_neighbors_IDs(std::vector<unsigned int> IDs, unsigned int k, unsigned short Class_dest);
	std::vector< std::pair<unsigned int, double> > k_nearest_neighbors_ID(unsigned int id, unsigned int k, unsigned short Class_dest);
	
	std::vector< std::vector< std::pair<unsigned int, double> > > neighborsWithinRadius_All(double radius, unsigned short Class_dest, unsigned short Class_src);
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
	typename SampleType::Pointer sample;
	typename TreeGeneratorType::Pointer treeGenerator;
	std::map< unsigned int, std::vector<double> > centerMap;
	std::map< unsigned int, std::vector<double> >::iterator It;
	std::map<int, MeasurementVectorType> idToCentroidMap;
	typename std::map<int, MeasurementVectorType>::iterator IdIt;
	typename TreeType::Pointer tree;
	NeighborGraph NG;
	vtkSmartPointer<vtkTable> graphtable;
	node_name nodeName;
	vtkSmartPointer<vtkTable> featureTable;
	int check;

	std::string convert2string(unsigned int id);
	int GetNodeIndex(unsigned int id);
	
	
};

#include "kNearestObjects.txx"
#endif
