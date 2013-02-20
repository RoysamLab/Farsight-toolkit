#ifndef MULTIFRAMECELLTRACKER_H
#define MULTIFRAMECELLTRACKER_H

#include "helpers.h"
//#include <time.h>
#include <math.h>
#include <iostream>                  
#include <utility>                   
#include <algorithm>
#include <map>
#include <set>
#include <limits>
// VTK includes:
#include <vtkTable.h>

#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/graphml.hpp>
#include <boost/property_map/dynamic_property_map.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/connected_components.hpp>
#include "itkRegionOfInterestImageFilter.h"
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#ifdef USE_KPLS
#include <PatternAnalysis/embrex/kpls.h>
#include <itkNormalVariateGenerator.h>
#endif 
#include <ftkCommon/ftkUtils.h>

using namespace helpers;
using namespace boost;
#if defined(_MSC_VER)
#pragma warning(disable: 4996)
#endif
//clock_t starttime,endtime,firsttime;

#define DEBUG
#ifdef DEBUG
#define SHORT(x) (strrchr(x,'\\') ? strrchr(x,'\\')+1: x)
#define _TRACE {printf("In-%s:%d:%s:\n",SHORT(__FILE__),__LINE__,__FUNCTION__);}
#define _ETRACE {printf("Entering-%s:%d:%s:\n",SHORT(__FILE__),__LINE__,__FUNCTION__);}
#define _LTRACE {printf("Leaving-%s:%d:%s:\n",SHORT(__FILE__),__LINE__,__FUNCTION__);}
#else
#define _TRACE
#define _ETRACE
#define _LTRACE
#endif

//#define TIC {starttime = clock();}
//#define TOC(x) {endtime = clock(); printf("Time for %s == %0.2lfs total == %0.2lfs\n",(x),double((endtime-starttime)*1.0/CLOCKS_PER_SEC),double((endtime-firsttime)*1.0/CLOCKS_PER_SEC));}

#define mxIsFinite(a) ((a)<1e6)
#define USE_VNL_HUNGARIAN 

#define CACHE_PREFIX "cache"
#define MAX_TIME 400
#define MAX_TAGS 4
#define MAX_LABEL 10000
#define VESSEL_CHANNEL 4 // FIXME : make it dynamic based on user input
#define PAUSE {printf("%d:>",__LINE__);scanf("%*d");}

class MultiFrameCellTracker
{
  public:
    MultiFrameCellTracker();
    ~MultiFrameCellTracker();
    void setTrackParameters(std::vector<std::pair<std::string,float> > parameters);
    void settrackresultFolders(std::vector<std::pair<std::string,std::string> > folders);
    void setTrackImages(ftk::Image::Pointer rawimage,ftk::Image::Pointer labelimage);
    void setChannelToTrack(int nucChannel){this->channel_to_track = nucChannel;};
    ftk::Image::Pointer getTrackImages(void);
    std::vector<std::vector<ftk::TrackPointFeatures> > getTrackFeatures(void);
    vtkSmartPointer<vtkTable> GetTimeFeaturesTable(void){return TimeFeaturesTable;};
    void ComputeTimeFeaturesTable(void);
    std::vector<ftk::TrackFeatures> GetTimeFeatures(void){return tfs;};

    void set_parameters_from_cmd(FeatureVariances feat_var){this->fvar = feat_var;};
    void set_inputs_from_cmd(std::vector< helpers::InputImageType::Pointer > inp_im, std::vector< helpers::LabelImageType::Pointer > lab_img);
    std::vector< helpers::LabelImageType::Pointer > get_ouput_to_cmd(){return this->output_images_;};

  private:
    struct TrackVertex
    {
      bool special;
      int t;
      int findex;
      int new_label;
      std::vector<int> vertlre;
      std::vector<int> sec_order_utility;
      bool selected_sec_order;
    };

    struct TrackEdge
    {
      int utility;
      bool coupled;		// for merge-split cells
      bool fixed;
      bool selected;		// which edge has been selected
      unsigned char type;
      std::vector<int> frontlre;
      std::vector<int> backlre;
      unsigned int index;
    };

    struct MergeCandidate
    {
      int t;
      int index1, index2;
      FeatureType f;
    };

    typedef boost::adjacency_list<boost::listS, boost::vecS, boost::bidirectionalS, TrackVertex, TrackEdge> TGraph;
    typedef std::vector< std::vector < FeatureType> > VVF;
    typedef std::vector< std::vector<TGraph::vertex_descriptor> > VVV;
    typedef std::vector< std::vector < MergeCandidate > > VVM;
    typedef std::vector< std::vector < helpers::LabelImageType::Pointer > > VVL;
    typedef std::vector< std::vector < helpers::InputImageType::Pointer > > VVR;

    bool run();
    void assignmentsuboptimal1(double *assignment, double *cost, double *distMatrixIn, int nOfRows, int nOfColumns);
    void writeGraphViz(char * filename);
    void writeGraphML (char * filename);
    void writeXGMML(char *filename);

    template <typename T>
      typename T::Pointer readImage(const char *filename)
      {
        printf("Reading %s ... \n",filename);
        typedef typename itk::ImageFileReader<T> ReaderType;
        typename ReaderType::Pointer reader = ReaderType::New();

        ReaderType::GlobalWarningDisplayOff();
        reader->SetFileName(filename);
        try
        {
          reader->Update();
        }
        catch(itk::ExceptionObject &err)
        {
          std::cerr << "ExceptionObject caught!" <<std::endl;
          std::cerr << err << std::endl;
          //return EXIT_FAILURE;
        }
        printf("Done\n");
        return reader->GetOutput();

      }

    template <typename T>
      int writeImage(typename T::Pointer im, const char* filename)
      {
        printf("Writing %s ... \n",filename);
        typedef typename itk::ImageFileWriter<T> WriterType;

        typename WriterType::Pointer writer = WriterType::New();
        writer->SetFileName(filename);
        writer->SetInput(im);
        try
        {
          writer->Update();
        }
        catch(itk::ExceptionObject &err)
        {
          std::cerr << "ExceptionObject caught!" <<std::endl;
          std::cerr << err << std::endl;
          return EXIT_FAILURE;
        }
        return EXIT_SUCCESS;
      }

    std::vector<std::map<int,int> > old_to_new;
    void set_debug_images(helpers::ColorImageType::Pointer in1,helpers::ColorImageType::Pointer in2, helpers::ColorImageType::Pointer in3)
    {
      debugimage1 = in1;
      debugimage2 = in2;
      debugimage3 = in3;
    }
    helpers::LabelImageType::Pointer getOutputAtTime(int t);
    FeatureVariances get_computed_variances()
    {
      return fvarnew;
    }
    std::string dataset_id;
    //////////////////////////////////////////////////////////////////////////////////////////////////

    //bool CompareFeaturesTime(FeatureType a, FeatureType b);
    void printFeatures(FeatureType f);
    std::vector<std::vector<bool> >  generate_all_binary_strings(int n);//Added to the class

    enum EDGE_TYPE {SPLIT,MERGE,APPEAR,DISAPPEAR,TRANSLATION};
    float overlap(float bb1[6],float bb2[6]);
    float get_boundary_dist(float x[3]);
    int compute_boundary_utility(float x[3]);
    float get_distance( float x1[3],float x2[3]);
    FeatureType get_merged_features(int, int,int);
    int add_disappear_vertices(int t);
    int add_appear_vertices(int t);
    int add_normal_edges(int tmin, int tmax);
    void populate_merge_candidates(int t);
    int add_merge_split_edges(int tmax);
    void solve_lip(void);
    void solve_higher_order(void);
    void prune(int);
    int compute_normal_utility(FeatureType f1, FeatureType f2);
    int compute_normal_utility(FeatureType f1, FeatureType f2, int counter, int counter1);
    void print_stats(void);
    int get_edge_type(TGraph::edge_descriptor);
    float compute_LRUtility(TGraph::edge_descriptor, TGraph::edge_descriptor);
    float compute_LRUtility_sum(TGraph::edge_descriptor, TGraph::edge_descriptor);
    float compute_LRUtility_product(TGraph::edge_descriptor, TGraph::edge_descriptor);
    float compute_LRUtility(FeatureType f1, FeatureType f2, FeatureType f3);
    void enforce_overlap();
    int is_overlapping(TGraph::edge_descriptor e);
    int is_alpha_overlapping(TGraph::edge_descriptor e,float alpha);
    std::pair<TGraph::edge_descriptor,bool> my_add_edge(TGraph::vertex_descriptor v1, TGraph::vertex_descriptor v2, int utility, bool coupled, bool fixed, bool selected, unsigned char type);
    void draw_line_for_edge(int num, TGraph::edge_descriptor e,helpers::VectorPixelType col1,helpers::VectorPixelType col2,int);
    void print_vertex(TGraph::vertex_descriptor v,int);
    bool is_merge_node(TGraph::vertex_descriptor vd);
    bool is_split_node(TGraph::vertex_descriptor vd);
    bool is_simple_node(TGraph::vertex_descriptor vd);
    TGraph::vertex_descriptor get_parent(TGraph::vertex_descriptor);// no error check implemented TODO
    TGraph::vertex_descriptor get_child(TGraph::vertex_descriptor);// no error check implemented TODO
    bool is_separate(MultiFrameCellTracker::TGraph::vertex_descriptor v1, MultiFrameCellTracker::TGraph::vertex_descriptor v2);
    float get_LRUtility(std::vector< TGraph::vertex_descriptor > desc);
    float get_misalignment_cost(FeatureType f1, FeatureType f2, FeatureType f3);
    void print_debug_info(void);
    void resolve_merges_and_splits(void);
    int my_connected_components(std::vector<int> &component);
    void print_all_LRUtilities(TGraph::vertex_descriptor v);
    void compute_feature_variances();
    void setData(VVF &fv, VVL &l, VVR &r);
    void summarize_tracking(ftk::Image::Pointer rawImg);	
    void createTrackFeatures(std::vector<FeatureType> summaryfvector[MAX_TIME][MAX_TAGS], std::vector<ftk::TrackFeatures> &tfs, int c,int num_t);
    void changeDataHierarchy(std::vector<ftk::TrackFeatures> vectrackfeatures);
    void convertItkImagesToftkImages(ftk::Image::Pointer labelImage,ftk::Image::Pointer dataImage,std::vector<helpers::LabelImageType::Pointer> &trackImages);
    int my_connected_components2(std::vector<int> &component);
    void compute_tracks_entropy();
    void append_track_entropy_to_tfs();
    std::map<int, std::vector<int> > ComputeEntropyUtilitiesAtTime(int t);
    void ComputeVertexEntropies(void);
    float get_LRUtility_Amin(std::vector< MultiFrameCellTracker::TGraph::vertex_descriptor > desc, int fileindex,int utilindex);

    TGraph g; 
    VVF fvector; 
    VVL limages;
    std::vector< helpers::LabelImageType::Pointer > output_images_;
    VVR rimages;
    VVV rmap;
    VVM m_cand;
    std::map<TGraph::edge_descriptor, TGraph::edge_descriptor> coupled_map;
    int UTILITY_MAX;
    FeatureVariances fvar;
    FeatureVariances fvarnew;
    int K;
    int channel_to_track;

    helpers::ColorImageType::Pointer debugimage1,debugimage2,debugimage3;
    bool edge_uniqueness(TGraph::edge_descriptor, TGraph::edge_descriptor);

    struct LREdge
    {
      TGraph::edge_descriptor front,back;
    };

    std::vector<std::map<int, std::vector<int> > > VertexUtilities;  // second order edge utilities through the vertex
    std::vector<std::map<int, float> > vertex_entropies;  // second order edge entropies through the vertex

    void writeXGMML_secondorder(char *,std::vector< LREdge >&,std::vector<int>&, IloNumArray& );

    ftk::Image::Pointer resultImages;
    std::vector<ftk::TrackFeatures> tfs;
    std::vector<std::vector<ftk::TrackPointFeatures> > timeFeatures; 	// Set a different data hierarchy for output of time features:
    vtkSmartPointer<vtkTable> TimeFeaturesTable;

    std::string numbersfile;
    std::string entropyfilename;
    std::string entropyfiledirectory;
    std::string debugprefix;
    std::string debugfilename;
    std::string debugfiledirectory;
    std::string resultfilename;
    std::string resultfiledirectory;

};
#endif




