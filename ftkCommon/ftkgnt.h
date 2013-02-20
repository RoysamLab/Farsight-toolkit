


#ifndef _ftkgnt_H_
#define _ftkgnt_H_

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <set>
#include <vnl/vnl_vector_fixed.h>
#include <vcl_vector.h>
#include <ftkFeatures/ftkLabelImageToFeatures.h>
#include "itkImage.h"

		struct VP {
		std::string label;
 		std::set<int> RPS;
		std::vector<double> scores;
			 };
	

class ftkgnt
{
public:

  typedef itk::Image< unsigned char, 3> InputImageType;
  typedef itk::Image< unsigned short, 3> OutputImageType;
  typedef ftk::LabelImageToFeatures< unsigned char,  unsigned short, 3 > FeatureCalcType;
  // Graph typedefs
  typedef boost::property<boost::vertex_name_t, std::string > VertexProperties;
  typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS,VertexProperties> RAGraph;
  typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS,VP> MTreeType;	
  typedef boost::graph_traits<RAGraph>::vertex_descriptor node;
  typedef boost::graph_traits<RAGraph>::edge_descriptor Edge;
  typedef boost::property_map<RAGraph, boost::vertex_name_t>::type node_name;  
  typedef boost::graph_traits<RAGraph>::adjacency_iterator AdjVertIt;
  typedef boost::property_map<MTreeType,boost::vertex_name_t>::type nodes_new; 
  typedef boost::graph_traits<MTreeType>::vertex_descriptor node_mt;
  typedef ftk::IntrinsicFeatures FeaturesType;	

    //: constructor
	ftkgnt();
 	//: destructor
	~ftkgnt();
  	
RAGraph RAG;
RAGraph BuildRAG(unsigned short id);
//RAGraph BuildRAG(unsigned short id,InputImageType::Pointer input,OutputImageType::Pointer output);
int GetNodeIndex(unsigned short id,RAGraph graph1);	
MTreeType BuildMergeTreeDcon(RAGraph R1, unsigned short id,std::vector< std::set<int> > hypothesis);
std::string convert2string(unsigned short id);
void Initialize(unsigned short id);
void runLabFilter(InputImageType::Pointer input, OutputImageType::Pointer output, bool CytoImage = false);
std::vector< std::set<int> > hypotheses;
//inline std::vector< std::set<int> > getHypothesis() {return hypotheses;};
FeatureCalcType::Pointer labelFilter;
std::vector<FeaturesType> allFeat;
std::vector< unsigned short > labelIndex;
void setFeats(std::vector<FeaturesType> allFeat,std::vector< unsigned short > labelIndex);

void setmaxVol(double vol);
void setmaxDepth(int depth);


private:

	double MAX_VOL;
	int MAX_DEPTH;
	bool cyto_image;
	
};


#endif
