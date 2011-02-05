#include "helpers.h"

#include  <time.h>
#include <iostream>                  
#include <utility>                   
#include <algorithm>
#include <map>
#include <set>
//#include <hash_map>
#include <limits>

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

using namespace helpers;
using namespace boost;
#if defined(_MSC_VER)
#pragma warning(disable: 4996)
#endif
clock_t start_t,end_t,first_t;

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

#define TIC {start_t = clock();}
#define TOC(x) {end_t = clock(); printf("Time for %s == %0.2lfs total == %0.2lfs\n",(x),double((end_t-start_t)*1.0/CLOCKS_PER_SEC),double((end_t-first_t)*1.0/CLOCKS_PER_SEC));}

bool file_exists(char *filename)
{
	FILE * fp = fopen(filename,"r");
	if(fp!=NULL)
	{
		fclose(fp);
		return true;
	}
	return false;
}


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

int num_files;
#define mxIsFinite(a) ((a)<1e6)

void assignmentsuboptimal1(double *assignment, double *cost, double *distMatrixIn, int nOfRows, int nOfColumns)
{
	bool infiniteValueFound, finiteValueFound, repeatSteps, allSinglyValidated, singleValidationFound;
	int n, row, col, tmpRow, tmpCol, nOfElements;
	int *nOfValidObservations, *nOfValidTracks;
	double value, minValue, *distMatrix, inf;

	inf = 1e10;

	/* make working copy of distance Matrix */
	nOfElements   = nOfRows * nOfColumns;
	distMatrix    = (double *)malloc(nOfElements * sizeof(double));
	for(n=0; n<nOfElements; n++)
		distMatrix[n] = distMatrixIn[n];

	/* initialization */
	*cost = 0;
#ifdef ONE_INDEXING
	for(row=0; row<nOfRows; row++)
		assignment[row] =  0.0;
#else
	for(row=0; row<nOfRows; row++)
		assignment[row] = -1.0;
#endif

	/* allocate memory */
	nOfValidObservations  = (int *)calloc(nOfRows,    sizeof(int));
	nOfValidTracks        = (int *)calloc(nOfColumns, sizeof(int));

	/* compute number of validations */
	infiniteValueFound = false;
	finiteValueFound  = false;
	for(row=0; row<nOfRows; row++)
		for(col=0; col<nOfColumns; col++)
			if(mxIsFinite(distMatrix[row + nOfRows*col]))
			{
				nOfValidTracks[col]       += 1;
				nOfValidObservations[row] += 1;
				finiteValueFound = true;
			}
			else
				infiniteValueFound = true;

	if(infiniteValueFound)
	{
		if(!finiteValueFound)
			return;

		repeatSteps = true;

		while(repeatSteps)
		{
			repeatSteps = false;

			/* step 1: reject assignments of multiply validated tracks to singly validated observations		 */
			for(col=0; col<nOfColumns; col++)
			{
				singleValidationFound = false;
				for(row=0; row<nOfRows; row++)
					if(mxIsFinite(distMatrix[row + nOfRows*col]) && (nOfValidObservations[row] == 1))
					{
						singleValidationFound = true;
						break;
					}

					if(singleValidationFound)
					{
						for(row=0; row<nOfRows; row++)
							if((nOfValidObservations[row] > 1) && mxIsFinite(distMatrix[row + nOfRows*col]))
							{
								distMatrix[row + nOfRows*col] = inf;
								nOfValidObservations[row] -= 1;							
								nOfValidTracks[col]       -= 1;	
								repeatSteps = true;				
							}
					}
			}

			/* step 2: reject assignments of multiply validated observations to singly validated tracks */
			if(nOfColumns > 1)			
			{	
				for(row=0; row<nOfRows; row++)
				{
					singleValidationFound = false;
					for(col=0; col<nOfColumns; col++)
						if(mxIsFinite(distMatrix[row + nOfRows*col]) && (nOfValidTracks[col] == 1))
						{
							singleValidationFound = true;
							break;
						}

						if(singleValidationFound)
						{
							for(col=0; col<nOfColumns; col++)
								if((nOfValidTracks[col] > 1) && mxIsFinite(distMatrix[row + nOfRows*col]))
								{
									distMatrix[row + nOfRows*col] = inf;
									nOfValidObservations[row] -= 1;
									nOfValidTracks[col]       -= 1;
									repeatSteps = true;								
								}
						}
				}
			}
		} /* while(repeatSteps) */

		/* for each multiply validated track that validates only with singly validated  */
		/* observations, choose the observation with minimum distance */
		for(row=0; row<nOfRows; row++)
		{
			if(nOfValidObservations[row] > 1)
			{
				allSinglyValidated = true;
				minValue = inf;
				for(col=0; col<nOfColumns; col++)
				{
					value = distMatrix[row + nOfRows*col];
					if(mxIsFinite(value))
					{
						if(nOfValidTracks[col] > 1)
						{
							allSinglyValidated = false;
							break;
						}
						else if((nOfValidTracks[col] == 1) && (value < minValue))
						{
							tmpCol   = col;
							minValue = value;
						}
					}
				}

				if(allSinglyValidated)
				{
#ifdef ONE_INDEXING
					assignment[row] = tmpCol + 1;
#else
					assignment[row] = tmpCol;
#endif
					*cost += minValue;
					for(n=0; n<nOfRows; n++)
						distMatrix[n + nOfRows*tmpCol] = inf;
					for(n=0; n<nOfColumns; n++)
						distMatrix[row + nOfRows*n] = inf;
				}
			}
		}

		/* for each multiply validated observation that validates only with singly validated  */
		/* track, choose the track with minimum distance */
		for(col=0; col<nOfColumns; col++)
		{
			if(nOfValidTracks[col] > 1)
			{
				allSinglyValidated = true;
				minValue = inf;
				for(row=0; row<nOfRows; row++)
				{
					value = distMatrix[row + nOfRows*col];
					if(mxIsFinite(value))
					{
						if(nOfValidObservations[row] > 1)
						{
							allSinglyValidated = false;
							break;
						}
						else if((nOfValidObservations[row] == 1) && (value < minValue))
						{
							tmpRow   = row;
							minValue = value;
						}
					}
				}

				if(allSinglyValidated)
				{
#ifdef ONE_INDEXING
					assignment[tmpRow] = col + 1;
#else
					assignment[tmpRow] = col;
#endif
					*cost += minValue;
					for(n=0; n<nOfRows; n++)
						distMatrix[n + nOfRows*col] = inf;
					for(n=0; n<nOfColumns; n++)
						distMatrix[tmpRow + nOfRows*n] = inf;
				}
			}
		}	
	} /* if(infiniteValueFound) */


	/* now, recursively search for the minimum element and do the assignment */
	while(true)
	{
		/* find minimum distance observation-to-track pair */
		minValue = inf;
		for(row=0; row<nOfRows; row++)
			for(col=0; col<nOfColumns; col++)
			{
				value = distMatrix[row + nOfRows*col];
				if(mxIsFinite(value) && (value < minValue))
				{
					minValue = value;
					tmpRow   = row;
					tmpCol   = col;
				}
			}

			if(mxIsFinite(minValue))
			{
#ifdef ONE_INDEXING
				assignment[tmpRow] = tmpCol+ 1;
#else
				assignment[tmpRow] = tmpCol;
#endif
				*cost += minValue;
				for(n=0; n<nOfRows; n++)
					distMatrix[n + nOfRows*tmpCol] = inf;
				for(n=0; n<nOfColumns; n++)
					distMatrix[tmpRow + nOfRows*n] = inf;			
			}
			else
				break;

	} /* while(true) */

	/* free allocated memory */
	free(nOfValidObservations);
	free(nOfValidTracks);


}

//void assignmentsuboptimal1(double *assignment, double *cost, double *distMatrixIn, int nOfRows, int nOfColumns)
#define USE_VNL_HUNGARIAN 
vcl_vector<unsigned int> getTimeAssociations(std::vector<FeaturesType> &a,std::vector<FeaturesType> &b)
{

#define DEBUG_RANK -1
	int rows =0;
	int cols =0;
	rows = a.size();
	cols = b.size();
	printf("Rows = %d Cols = %d\n",rows,cols);
	double overlap;
#ifdef USE_VNL_HUNGARIAN
	vnl_matrix<double> mat(rows,cols);
	printf("Allocated Rows = %d Cols = %d\n",mat.rows(),mat.cols());
	//int pa =0; 
	//int pb =0;
	for(int cr = 0; cr<rows; cr++)
	{

		bool enforce_overlap = false;
		for(int cc=0; cc<cols; cc++)
		{
			overlap = features_box_overlap(a[cr],b[cc]);
			if(overlap>1) // volume of overlap > 1
			{
				enforce_overlap = true;
				break;
			}
		}
		/*	if(rank==DEBUG_RANK)
		{
		if(enforce_overlap)
		printf("%d/%d: I did enforce overlap for %d\n",rank,npes,cr);
		else
		printf("%d/%d: I did not enforce overlap for %d\n",rank,npes,cr);
		}*/
		for(int cc =0; cc<cols; cc++)
		{
			mat(cr,cc) = features_diff(a[cr],b[cc],enforce_overlap);
		}
	}
	printf("About to call vnl_hungarian_algorithm\n");
	vcl_vector<unsigned int> ret = vnl_hungarian_algorithm(mat);
	printf("Returned from vnl_hungarian_algorithm\n");
	for(unsigned int counter=0; counter< ret.size(); counter++)
	{
		if(mxIsFinite(ret[counter]))
		{
			if(!mxIsFinite(mat(counter,ret[counter])))
			{
				ret[counter] = static_cast<unsigned int>(-1);
			}
		}
	}
	printf("Returning from getTimeAssociations\n");
	return ret;
#else
	double * assignment  = (double*) malloc(rows*sizeof(double));
	double * cost = (double*) malloc(sizeof(double));
	double * distMatrixIn = (double*) malloc(rows *cols*sizeof(double));

	if(assignment == NULL)
		printf("Couldn't allocate memory assignment\n");

	if(distMatrixIn == NULL)
		printf("Couldn't allocate memory for distMatrixIn\n");

	if(cost == NULL)
		printf("Couldn't allocate memory for cost\n");

	printf("%d/%d About to assign datamatrix values\n",rank,npes);
	int pa =0; 
	int pb =0;
	for(int cr = 0; cr<rows; cr++)
	{

		//if(rank==0)
		//{
		//	printf("0/124: bbox %d %d %d %d %d %d \n",a[cr].bbox.sx,a[cr].bbox.sy,a[cr].bbox.sz,a[cr].bbox.ex,a[cr].bbox.ey,a[cr].bbox.ez);
		//}
		bool enforce_overlap = false;
		for(int cc=0; cc<cols; cc++)
		{
			overlap = features_box_overlap(a[cr],b[cc]);
			if(overlap>1) // volume of overlap > 1
			{
				enforce_overlap = true;
				break;
			}
		}
		//if(rank==DEBUG_RANK)
		{
			if(enforce_overlap)
				printf("%d/%d: I did enforce overlap for %d\n",rank,npes,cr);
			else
				printf("%d/%d: I did not enforce overlap for %d\n",rank,npes,cr);
		}
		for(int cc =0; cc<cols; cc++)
		{
			distMatrixIn[cr+rows*cc] = features_diff(a[cr],b[cc],enforce_overlap);
			//if(rank==DEBUG_RANK)
			{
				//	if(cr==0)
				{
					printf("distmatrix[%d,%d]=%0.3lf\n",cr,cc,distMatrixIn[cr+rows*cc]);
				}
			}
		}

	}
	printf("%d/%d About to call assignmentsuboptimal1\n",rank,npes);
	assignmentsuboptimal1(assignment,cost,distMatrixIn,rows,cols);
	printf("%d/%d Exited assignmentsuboptimal1\n",rank,npes);
	vcl_vector<unsigned int> vec(rows);
	int assigned_count=0;
	//	if(rank==DEBUG_RANK)
	//	{
	//		for(int counter=0; counter< rows; counter++)
	//		{
	//			printf("%d/%d assignment[%d] = %0.3lf\n",rank,npes,counter,assignment[counter]);
	//		}
	//	}
	printf("%d/%d About to start assigning vec values\n",rank,npes);
	for(int cr = 0; cr<rows; cr++)
	{
		if(assignment[cr]>-0.1)
		{
			vec[cr]=static_cast<unsigned int>(assignment[cr]+0.5);
			//			if(rank==DEBUG_RANK)
			//				printf("%d/%d : assigned_nums[%d] = %d\n",rank,npes,cr,vec[cr]);
			assigned_count++;
		}
		else
		{
			vec[cr]=static_cast<unsigned int>(-1);
		}
	}
	//	if(rank==DEBUG_RANK)
	//	{
	//		printf("%d/%d: assigned_count = %d\n",rank,npes,assigned_count);
	//	}
	//free(assignment); FIXME : was giving glibc corruption errors;
	//free(cost);
	//free(distMatrixIn);
	return vec;

#endif

}

void getArrayFromStdVector(std::vector<FeaturesType> &f, FeaturesType	*&farray)
{
	farray = new FeaturesType[f.size()];
	for(unsigned int counter=0; counter<f.size(); counter++)
	{
		farray[counter]=f[counter];
	}
}

void getStdVectorFromArray(FeaturesType *farray, int n,std::vector<FeaturesType> &f)
{
	f.reserve(n);
	for(int counter=0; counter<n; counter++)
	{
		f.push_back(farray[counter]);
	}
}




//#define CACHE_PREFIX "D:/ucb dataset/output/ena/cache"
#define CACHE_PREFIX "cache"
#define MAX_TIME 200
#define MAX_TAGS 4
#define MAX_LABEL 10000
#define VESSEL_CHANNEL 4 // FIXME : make it dynamic based on user input
#define PAUSE {printf("%d:>",__LINE__);scanf("%*d");}

bool compare(FeaturesType a, FeaturesType b)
{
	return a.time<b.time;
}
void createTrackFeatures(std::vector<FeaturesType> fvector[MAX_TIME][MAX_TAGS], std::vector<ftk::TrackFeatures> &tfs, int c,int num_t)
{
	int max_track_num = 0;
	for(int t = 0; t< num_t; t++)
	{
		for(unsigned int counter=0; counter< fvector[t][c-1].size(); counter++)
		{
			max_track_num = MAX(max_track_num,fvector[t][c-1][counter].num);
		}
	}

	for(int counter=1; counter <= max_track_num; counter++)
	{
		ftk::TrackFeatures trackf;
		trackf.intrinsic_features.clear();
		for(int t = 0; t< num_t;t++)
		{
			for(unsigned int counter1 = 0; counter1 < fvector[t][c-1].size(); counter1++)
			{
				if(fvector[t][c-1][counter1].num == counter)
				{
					trackf.intrinsic_features.push_back(fvector[t][c-1][counter1]);
				}
			}
		}
		std::sort(trackf.intrinsic_features.begin(),trackf.intrinsic_features.end(),compare);
		tfs.push_back(trackf);
		//PRINTF("Added %d elements to tfs\n",counter);
	}
}






class CellTracker{

	struct TrackVertex{
		bool special;
		int t;
		int findex;
		int new_label;
		std::vector<int> vertlre;
	};

	struct TrackEdge{
		int utility;
		bool coupled;
		bool fixed;
		bool selected;
		unsigned char type;
		std::vector<int> frontlre;
		std::vector<int> backlre;
		unsigned int index;
	};


	struct MergeCandidate{
		int t;
		int index1, index2;
		FeaturesType f;
	};
public:
	typedef boost::adjacency_list<boost::listS, boost::vecS, boost::bidirectionalS, TrackVertex, TrackEdge> TGraph;
	typedef std::vector< std::vector < FeaturesType> > VVF;
	typedef std::vector< std::vector<TGraph::vertex_descriptor> > VVV;
	typedef std::vector< std::vector < MergeCandidate > > VVM;
	typedef std::vector< std::vector < LabelImageType::Pointer > > VVL;
	typedef std::vector< std::vector < InputImageType::Pointer > > VVR;

	CellTracker(VVF &fv, FeatureVariances &fvariances, VVL &l, VVR &r)
	{
		fvector  = fv;
		old_to_new.resize(fv.size());
		fvar = fvariances;
		limages = l; // copy only pointers
		rimages = r;
		UTILITY_MAX = 1<<14;
		K = 2;
	}

	void run();
	void writeGraphViz(char * filename);
	void writeGraphML (char * filename);
	void writeXGMML(char *filename);
	
	std::vector<std::map<int,int>> old_to_new;
	void set_debug_images(ColorImageType::Pointer in1,ColorImageType::Pointer in2, ColorImageType::Pointer in3)
	{
		debugimage1 = in1;
		debugimage2 = in2;
		debugimage3 = in3;
	}
	LabelImageType::Pointer getOutputAtTime(int t);
	FeatureVariances get_computed_variances()
	{
		return fvarnew;
	}
private:
	enum EDGE_TYPE {SPLIT,MERGE,APPEAR,DISAPPEAR,TRANSLATION};
	float overlap(float bb1[6],float bb2[6]);
	float get_boundary_dist(float x[3]);
	int compute_boundary_utility(float x[3]);
	float get_distance( float x1[3],float x2[3]);
	FeaturesType get_merged_features(int, int,int);
	int add_disappear_vertices(int t);
	int add_appear_vertices(int t);
	int add_normal_edges(int tmin, int tmax);
	void populate_merge_candidates(int t);
	int add_merge_split_edges(int tmax);
	void solve_lip(void);
	void solve_higher_order(void);
	void prune(int);
	int compute_normal_utility(FeaturesType f1, FeaturesType f2);
	void print_stats(void);
	int get_edge_type(TGraph::edge_descriptor);
	float compute_LRUtility(TGraph::edge_descriptor, TGraph::edge_descriptor);
	float compute_LRUtility_sum(TGraph::edge_descriptor, TGraph::edge_descriptor);
	float compute_LRUtility_product(TGraph::edge_descriptor, TGraph::edge_descriptor);
	float compute_LRUtility(FeaturesType f1, FeaturesType f2, FeaturesType f3);
	void enforce_overlap();
	int is_overlapping(TGraph::edge_descriptor e);
	int is_alpha_overlapping(TGraph::edge_descriptor e,float alpha);
	std::pair<TGraph::edge_descriptor,bool> my_add_edge(TGraph::vertex_descriptor v1, TGraph::vertex_descriptor v2, int utility, bool coupled, bool fixed, bool selected, unsigned char type);
	void draw_line_for_edge(int num, TGraph::edge_descriptor e,VectorPixelType col1,VectorPixelType col2,int);
	void print_vertex(TGraph::vertex_descriptor v,int);
	bool is_merge_node(TGraph::vertex_descriptor vd);
	bool is_split_node(TGraph::vertex_descriptor vd);
	bool is_simple_node(TGraph::vertex_descriptor vd);
	TGraph::vertex_descriptor get_parent(TGraph::vertex_descriptor);// no error check implemented TODO
	TGraph::vertex_descriptor get_child(TGraph::vertex_descriptor);// no error check implemented TODO
	bool is_separate(CellTracker::TGraph::vertex_descriptor v1, CellTracker::TGraph::vertex_descriptor v2);
	float get_LRUtility(std::vector< TGraph::vertex_descriptor > desc);
	float get_misalignment_cost(FeaturesType f1, FeaturesType f2, FeaturesType f3);
    void print_debug_info(void);
	void resolve_merges_and_splits(void);
	int my_connected_components(std::vector<int> &component);
	void print_all_LRUtilities(TGraph::vertex_descriptor v);
	void compute_feature_variances();
	
	//variables

	TGraph g; 
	VVF fvector; 
	VVL limages;
	VVR rimages;
	VVV rmap;
	VVM m_cand;
	std::map<TGraph::edge_descriptor, TGraph::edge_descriptor> coupled_map;
	int UTILITY_MAX;
	FeatureVariances fvar;
	FeatureVariances fvarnew;
	int K;

	ColorImageType::Pointer debugimage1,debugimage2,debugimage3;
	bool edge_uniqueness(TGraph::edge_descriptor, TGraph::edge_descriptor);

	struct LREdge{
		TGraph::edge_descriptor front,back;
	};
	void writeXGMML_secondorder(char *,std::vector< LREdge >&,std::vector<int>&, IloNumArray& );

};

void printFeatures(FeaturesType f)
{
	printf("\nIntrinsicFeatures:\n");
	printf("\t\n");
	printf("\t X = %d Y = %d Z = %d\n", int(f.Centroid[0]+0.5), int(f.Centroid[1]+0.5), int(f.Centroid[2]+0.5));
	printf("\t Num = %d time = %d\n", f.num, f.time);
	printf("\t Volume = %d\n\n", int(f.ScalarFeatures[FeaturesType::VOLUME]));
}

float CellTracker::overlap(float bb1[6],float bb2[6])
{
	float sx,sy,sz;
	float ex,ey,ez;
	sx = MAX(bb1[0],bb2[0]);
	sy = MAX(bb1[2],bb2[2]);
	sz = MAX(bb1[4],bb2[4]);
	ex = MIN(bb1[1],bb2[1]);
	ey = MIN(bb1[3],bb2[3]);
	ez = MIN(bb1[5],bb2[5]);

	float overlap=0;
	if((sx<ex) && (sy<ey) && (sz<ez))
	{
		overlap = (ex-sx+1)*(ey-sy+1)*(ez-sz+1);
	}

	return overlap;
}
float CellTracker::get_boundary_dist(float x[3])
{
	float minv = x[0]*fvar.spacing[0];
	minv = MIN(minv,x[1]*fvar.spacing[1]);
	minv = MIN(minv,x[2]*fvar.spacing[2]);
	minv = MIN(minv, (fvar.BoundingBox[1] - x[0])*fvar.spacing[0]);
	minv = MIN(minv, (fvar.BoundingBox[3] - x[1])*fvar.spacing[1]);
	minv = MIN(minv, (fvar.BoundingBox[5] - x[2])*fvar.spacing[2]);
	//	printf("Boundary dist = %f %f %f %f\n", minv,x[0],x[1],x[2]);
	return minv;
}

int CellTracker::compute_boundary_utility(float x[3])
{
	//return 1;

	float dist = get_boundary_dist(x);
	int retval = fvar.AD_prior*UTILITY_MAX*(exp(-dist*dist/2.0/fvar.distVariance));

	//if(get_distance(x,test) < 30 || get_distance(x,test1) <30)
	//{
	//	printf("retval = %d x = %0.0f y = %0.0f z = %0.0f\n", retval,x[0],x[1],x[2]);
	//	PAUSE;
	//}
	return retval;
}

int CellTracker::compute_normal_utility(FeaturesType f1, FeaturesType f2)
{
	float utility = 0;
	for(int counter=0; counter< FeaturesType::N; counter++)
	{
		utility += (f1.ScalarFeatures[counter]-f2.ScalarFeatures[counter])*(f1.ScalarFeatures[counter]-f2.ScalarFeatures[counter])/fvar.variances[counter];
	}

	float dist = get_distance(f1.Centroid,f2.Centroid);
	utility += (dist-fvar.distMean)*(dist-fvar.distMean)/fvar.distVariance;
	utility += (abs(f1.time-f2.time)-1)*(abs(f1.time-f2.time)-1)/fvar.timeVariance;
	float ovlap = overlap(f1.BoundingBox,f2.BoundingBox);
	ovlap = 1-(ovlap)/MIN(f1.ScalarFeatures[FeaturesType::BBOX_VOLUME],f2.ScalarFeatures[FeaturesType::BBOX_VOLUME]);
	utility += ovlap*ovlap/fvar.overlapVariance;
	utility /= 2.0;
	utility = fvar.T_prior*UTILITY_MAX*(exp(-utility));
	if(utility < 0)
	{
		printf("returning negative utility\n");
		scanf("%*d");
	}
	return utility;
}

float CellTracker::get_distance( float x1[3],float x2[3])
{
	float dist  = 0;
	dist += (x1[0]-x2[0])*(x1[0]-x2[0])*fvar.spacing[0]*fvar.spacing[0];
	dist += (x1[1]-x2[1])*(x1[1]-x2[1])*fvar.spacing[1]*fvar.spacing[1];
	dist += (x1[2]-x2[2])*(x1[2]-x2[2])*fvar.spacing[2]*fvar.spacing[2];
	dist = sqrt(dist);
	return dist;
}
int CellTracker::add_disappear_vertices(int t)
{
	_ETRACE;
	using boost::graph_traits;
	typedef graph_traits<typename TGraph>::vertex_iterator vertex_iter;
	TGraph::vertex_descriptor v;
	TGraph::edge_descriptor e;

	vertex_iter vi, vend,v_next;
	bool added;
	int ret_count = 0;
	for(boost::tie(vi,vend) = vertices(g); vi != vend; vi=v_next)
	{
		v_next = vi;
		++v_next;
		//printf("hi t = %d\n",t);
		if(g[*vi].t == t-1)
		{
			//printf("hi 1\n");
			if(g[*vi].special == 0)
			{
				//printf("hi 2\n");
				//if(get_boundary_dist(fvector[t-1][g[*vi].findex].Centroid) < 4.0*sqrt(fvar.distVariance))
				{
					v = add_vertex(g);
					g[v].special = 1;
					g[v].t = t;
					tie(e,added) = add_edge(*vi,v,g);
					if(added)
					{
						g[e].coupled = 0;
						g[e].fixed = 0;
						g[e].selected = 0;
						g[e].utility = compute_boundary_utility(fvector[t-1][g[*vi].findex].Centroid);
						ret_count ++;
					}
					float test[3] = { 182,299,6};
					float test1[3] = { 170, 285, 6};
					if( get_distance(test1,fvector[t-1][g[*vi].findex].Centroid) < 30 && t==20)
					{
						printFeatures(fvector[t-1][g[*vi].findex]);
						printf("Utility for disappearing = %d\n",g[e].utility);
					}
				}
			}
		}
	}
	_LTRACE;
	//printf(" I added %d disappear edges\n",ret_count);
	//PAUSE;
	return ret_count;
}

int CellTracker::add_appear_vertices(int t)
{
	_ETRACE;
	using boost::graph_traits;
	typedef graph_traits<typename TGraph>::vertex_iterator vertex_iter;
	TGraph::vertex_descriptor v;
	TGraph::edge_descriptor e;

	vertex_iter vi, vend;
	bool added;
	int ret_count = 0;
	for(tie(vi,vend) = vertices(g); vi != vend; ++vi)
	{
		if(g[*vi].t == t+1)
		{
			if(g[*vi].special == 0)
			{
				//printf("findex = %d fvector[%d].size() = %d\n",g[*vi].findex,t,fvector[t].size());
				//if(get_boundary_dist(fvector[t+1][g[*vi].findex].Centroid) <  4.0*sqrt(fvar.distVariance))
				{
					_TRACE;
					v = add_vertex(g);
					g[v].special = 1;
					g[v].t = t;

					tie(e,added) = add_edge(v,*vi,g);
					if(added)
					{
						g[e].coupled = 0;
						g[e].fixed = 0;
						g[e].selected = 0;
						g[e].utility = compute_boundary_utility(fvector[t+1][g[*vi].findex].Centroid);
						ret_count ++;
					}
					float test[3] = { 182,299,6};
					float test1[3] = { 170, 285, 6};
					if( get_distance(test,fvector[t+1][g[*vi].findex].Centroid) < 30 && t==20)
					{
						printFeatures(fvector[t+1][g[*vi].findex]);
						printf("Utility for appearing = %d\n",g[e].utility);
					}
					_TRACE;
				}
			}
		}
	}
	_LTRACE;
	//printf(" I added %d appear edges\n",ret_count);
	//PAUSE;
	return ret_count;
}

int CellTracker::add_normal_edges(int tmin, int tmax)
{
	_ETRACE;
	TGraph::edge_descriptor e;
	bool added = false;
	int nec = 0;
	float epsilon = 50;
	int tried1 =0,tried2 = 0 ;
	for(int counter=0; counter < fvector[tmax].size(); counter++)
	{
		// for every vertex, find the correspondences in the previous frames
		TGraph::vertex_descriptor v = rmap[tmax][counter];
		if(TGraph::null_vertex() != v) // do we have a non-null vertex? then ...
		{
			tried1++;
			for(int t = tmin; t < tmax; ++t)
			{

				for(int counter1 = 0; counter1 < fvector[t].size(); counter1++)
				{
					tried2++;
					float dist = get_distance(fvector[t][counter1].Centroid,fvector[tmax][counter].Centroid);

					if(dist < 4.0*sqrt(fvar.distVariance)*(abs(t-tmax)))
					{
					
						
						//add the edge
						//if(dist>100)
						//	printf("distance is greater than 100\n");
						tie(e,added) = add_edge(rmap[t][counter1],v,g);
						if(added)
						{
							g[e].coupled = 0;
							g[e].fixed = 0;
							g[e].selected = 0;
							float test[]={182,297,6};
							if(tmax==21 && t==19 && get_distance(fvector[tmax][counter].Centroid,test)<30)
							{
							printFeatures(fvector[t][counter1]);
							printFeatures(fvector[tmax][counter]);
							printf("g[e].utility = %d\n", g[e].utility);
							//PAUSE;
							}
							g[e].utility = compute_normal_utility(fvector[t][counter1],fvector[tmax][counter]);
							nec++;
						}

					}
				}
			}
		}
	}
	//printf("add_normal_edges: I tried %d %d\n",tried1, tried2);
	return nec;
}

FeaturesType CellTracker::get_merged_features(int t1, int i1, int i2)
{
	int t2 = t1;
	LabelImageType::Pointer p1,p2;
	InputImageType::Pointer r1,r2;
	p1 = limages[t1][i1];
	p2 = limages[t2][i2];
	r1 = rimages[t1][i1];
	r2 = rimages[t2][i2];


	LabelImageType::SizeType ls;

	int lbounds[6];

	lbounds[0] = MIN(fvector[t1][i1].BoundingBox[0],fvector[t2][i2].BoundingBox[0]);
	lbounds[2] = MIN(fvector[t1][i1].BoundingBox[2],fvector[t2][i2].BoundingBox[2]);
	lbounds[4] = MIN(fvector[t1][i1].BoundingBox[4],fvector[t2][i2].BoundingBox[4]);
	lbounds[1] = MAX(fvector[t1][i1].BoundingBox[1],fvector[t2][i2].BoundingBox[1]);
	lbounds[3] = MAX(fvector[t1][i1].BoundingBox[3],fvector[t2][i2].BoundingBox[3]);
	lbounds[5] = MAX(fvector[t1][i1].BoundingBox[5],fvector[t2][i2].BoundingBox[5]);

	ls[0] = lbounds[1]-lbounds[0]+1;
	ls[1] = lbounds[3]-lbounds[2]+1;
	ls[2] = lbounds[5]-lbounds[4]+1;

	LabelImageType::Pointer p = LabelImageType::New();
	InputImageType::Pointer r = InputImageType::New();
	LabelImageType::IndexType lindex;
	lindex.Fill(0);
	LabelImageType::RegionType lregion;
	lregion.SetIndex(lindex);
	lregion.SetSize(ls);
	p->SetRegions(lregion);
	p->Allocate();
	p->FillBuffer(0);
	r->SetRegions(lregion);
	r->Allocate();
	r->FillBuffer(0);
	LabelIteratorType liter1(p1,p1->GetLargestPossibleRegion());
	IteratorType riter1(r1,r1->GetLargestPossibleRegion());


	lindex[0] = fvector[t1][i1].BoundingBox[0]-lbounds[0];
	lindex[1] = fvector[t1][i1].BoundingBox[2]-lbounds[2];
	lindex[2] = fvector[t1][i1].BoundingBox[4]-lbounds[4];

	lregion.SetSize(p1->GetLargestPossibleRegion().GetSize());
	lregion.SetIndex(lindex);

	LabelIteratorType liter(p,lregion);
	IteratorType riter(r,lregion);
	for(liter1.GoToBegin(),riter1.GoToBegin(),liter.GoToBegin(),riter.GoToBegin();!liter1.IsAtEnd(); ++liter1,++riter1,++liter,++riter)
	{
		if(liter1.Get()==fvector[t1][i1].num)
			liter.Set(255);
		riter.Set(riter1.Get());
	}

	LabelIteratorType liter2(p2,p2->GetLargestPossibleRegion());
	IteratorType riter2(r2,r2->GetLargestPossibleRegion());

	lindex[0] = fvector[t2][i2].BoundingBox[0]-lbounds[0];
	lindex[1] = fvector[t2][i2].BoundingBox[2]-lbounds[2];
	lindex[2] = fvector[t2][i2].BoundingBox[4]-lbounds[4];
	lregion.SetIndex(lindex);
	lregion.SetSize(p2->GetLargestPossibleRegion().GetSize());

	liter = LabelIteratorType(p,lregion);
	riter = IteratorType(r,lregion);

	for(liter2.GoToBegin(),riter2.GoToBegin(),liter.GoToBegin(),riter.GoToBegin();!liter2.IsAtEnd(); ++liter2,++liter,++riter,++riter2)
	{
		if(liter2.Get()==fvector[t2][i2].num)
			liter.Set(255);
		riter.Set(riter2.Get());
	}


	std::vector<FeaturesType> f1;
	getFeatureVectorsFarsight(p,r,f1,t1,fvector[t1][i1].tag);

	FeaturesType f = f1[0];
	f.Centroid[0]+=lbounds[0];
	f.Centroid[1]+=lbounds[2];
	f.Centroid[2]+=lbounds[4];
	f.WeightedCentroid[0]+=lbounds[0];
	f.WeightedCentroid[2]+=lbounds[2];
	f.WeightedCentroid[4]+=lbounds[4];
	f.BoundingBox[0]+=lbounds[0];
	f.BoundingBox[1]+=lbounds[0];
	f.BoundingBox[2]+=lbounds[2];
	f.BoundingBox[3]+=lbounds[2];
	f.BoundingBox[4]+=lbounds[4];
	f.BoundingBox[5]+=lbounds[5];
	return f;
}
void  CellTracker::populate_merge_candidates(int t)
{
	MergeCandidate m;
	std::vector<MergeCandidate> vm;
	std::vector<MergeCandidate> nullvm;
	nullvm.clear();
	vm.clear();
	for(int counter = 0; counter < fvector[t].size(); counter++)
	{
		for(int counter1 = counter+1; counter1 < fvector[t].size(); counter1++)
		{
			//if( MIN(fvector[t][counter].ScalarFeatures[FeaturesType::BBOX_VOLUME],fvector[t][counter1].ScalarFeatures[FeaturesType::BBOX_VOLUME])/overlap(fvector[t][counter].BoundingBox,fvector[t][counter1].BoundingBox)< 4.0*sqrt(fvar.overlapVariance))
			if(get_distance(fvector[t][counter1].Centroid,fvector[t][counter].Centroid)<4.0*sqrt(fvar.distVariance))
			{
				m.t = t;
				m.index1 = counter;
				m.index2 = counter1;
				//m.f = get_merged_features(t,counter,counter1); // FIXME WRONG
				vm.push_back(m);
			}
		}
	}
	while(m_cand.size()<t)
	{
		m_cand.push_back(nullvm); //CHECK
	}
	m_cand.push_back(vm);
}


int CellTracker::add_merge_split_edges(int tmax)
{
	_TRACE
		int msec = 0;
	if(m_cand.size() < tmax)
	{
		printf("something is wrong.. plz check to make sure m_cand are populated correctly\n");
		_exit(1);
	}
	if(m_cand.size() < tmax + 1)
	{
		populate_merge_candidates(tmax);
	}

	TGraph::edge_descriptor e1,e2;
	bool added1, added2;
	for(int tcounter = MAX(tmax-K,0);tcounter <=tmax-1; tcounter++)
	{
		for(int counter=0; counter< m_cand[tcounter].size(); counter++)
		{
			for(int counter1 = 0; counter1 < fvector[tmax].size(); counter1++)
			{
				int i1 = m_cand[tcounter][counter].index1;
				int i2 = m_cand[tcounter][counter].index2;

				float centroid[3];
				for(int i = 0; i < 3; i++)
					centroid[i] = (fvector[tcounter][i1].Centroid[i] + fvector[tcounter][i2].Centroid[i])/2.0;

				if(get_distance(fvector[tcounter][i1].Centroid,fvector[tmax][counter1].Centroid)<4.0*sqrt(fvar.distVariance) && get_distance(fvector[tcounter][i2].Centroid,fvector[tmax][counter1].Centroid)<4.0*sqrt(fvar.distVariance))
				{
					FeaturesType lc1 = fvector[tcounter][i1];
					FeaturesType lc2 = fvector[tcounter][i2];

					int volume_factor = 1;

					lc1.ScalarFeatures[FeaturesType::VOLUME] = (lc1.ScalarFeatures[FeaturesType::VOLUME] + lc2.ScalarFeatures[FeaturesType::VOLUME])/volume_factor;
					lc2.ScalarFeatures[FeaturesType::VOLUME] = lc1.ScalarFeatures[FeaturesType::VOLUME];

					FeaturesType lc3 = fvector[tmax][counter1];
					lc3.ScalarFeatures[FeaturesType::VOLUME] /=volume_factor;

					int utility1 = compute_normal_utility(lc1,lc3);//fvector[tmax][counter1]);
					int utility2 = compute_normal_utility(lc2,lc3);//fvector[tmax][counter1]);

					tie(e1,added1) = add_edge(rmap[tcounter][i1],rmap[tmax][counter1],g);
					tie(e2,added2) = add_edge(rmap[tcounter][i2],rmap[tmax][counter1],g);
					if(added1&&added2)
					{
						g[e1].coupled = 1;
						g[e2].coupled = 1;
						coupled_map[e1] = e2;
						coupled_map[e2] = e1;
						g[e1].utility = (utility1+utility2)/fvar.T_prior*fvar.MS_prior;
						g[e2].utility = (utility1+utility2)/fvar.T_prior*fvar.MS_prior;
						if(utility1+utility2 == 0)
						{
							//printf("debug_for_merge :\n");
							//print_vertex(rmap[tcounter][i1],1);
							//printf("\n\n\n");
							//print_vertex(rmap[tcounter][i2],1);
							//printf("\n\n\n");
//							PAUSE;
						}
						g[e1].fixed = 0;
						g[e2].fixed = 0;
						g[e1].selected = 0;
						g[e2].selected = 0;
						msec+=2;
					}
				}
			}
		}
	}
	for(int counter=0; counter< m_cand[tmax].size(); counter++)
	{
		for(int tcounter = MAX(tmax-K,0);tcounter <=tmax-1; tcounter++)
		{		
			for(int counter1 = 0; counter1 < fvector[tcounter].size(); counter1++)
			{
				int i1 = m_cand[tmax][counter].index1;
				int i2 = m_cand[tmax][counter].index2;

				float centroid[3];
				for(int i = 0; i < 3; i++)
					centroid[i] = (fvector[tmax][i1].Centroid[i] + fvector[tmax][i2].Centroid[i])/2.0;

				if(get_distance(fvector[tmax][i1].Centroid,fvector[tcounter][counter1].Centroid)<4.0*sqrt(fvar.distVariance) && get_distance(fvector[tmax][i2].Centroid,fvector[tcounter][counter1].Centroid)<4.0*sqrt(fvar.distVariance))
				{
					FeaturesType lc1 = fvector[tmax][i1];
					FeaturesType lc2 = fvector[tmax][i2];

					int volume_factor = 1;

					lc1.ScalarFeatures[FeaturesType::VOLUME] = (lc1.ScalarFeatures[FeaturesType::VOLUME] + lc2.ScalarFeatures[FeaturesType::VOLUME])/volume_factor;
					lc2.ScalarFeatures[FeaturesType::VOLUME] = lc1.ScalarFeatures[FeaturesType::VOLUME];

					FeaturesType lc3 = fvector[tcounter][counter1];
					lc3.ScalarFeatures[FeaturesType::VOLUME] /=volume_factor;
					
					int utility1 = compute_normal_utility(lc1,lc3);//fvector[tcounter][counter1]);
					int utility2 = compute_normal_utility(lc2,lc3);//fvector[tcounter][counter1]);

					tie(e1,added1) = add_edge(rmap[tcounter][counter1],rmap[tmax][i1],g);
					tie(e2,added2) = add_edge(rmap[tcounter][counter1],rmap[tmax][i2],g);
					if(added1&&added2)
					{
						g[e1].coupled = 1;
						g[e2].coupled = 1;
						coupled_map[e1] = e2;
						coupled_map[e2] = e1;
						g[e1].utility = (utility1+utility2)/fvar.T_prior*fvar.MS_prior;
						g[e2].utility = (utility1+utility2)/fvar.T_prior*fvar.MS_prior;
						if(utility1+utility2 == 0)
						{
							if(tcounter==2 && tmax == 3)
							{
							/*printf("debug_for_merge :\n");
							print_vertex(rmap[tcounter][counter1],0);
							printf("\n\n\n");
							print_vertex(rmap[tmax][i1],0);
							printf("\n\n\n");
							print_vertex(rmap[tmax][i2],0);
							printf("\n\n\n");
							PAUSE;*/
							}
						}
						g[e1].fixed = 0;
						g[e2].fixed = 0;
						g[e1].selected = 0;
						g[e2].selected = 0;
						msec+=2;
					}
				}
			}
		}
	}
	return msec;
}
/*
struct TrackEdge{
int utility;
bool coupled;
bool fixed;
bool selected;
};
*/
int CellTracker::get_edge_type(TGraph::edge_descriptor ei)
{
	if(g[ei].coupled == 1)
	{
		TGraph::edge_descriptor ed = coupled_map[ei];
		if(source(ei,g)==source(ed,g))
		{
			g[ei].type = SPLIT;
			return SPLIT;
		}
		else if(target(ei,g) == target(ed,g))
		{
			g[ei].type = MERGE;
			return MERGE;
		}
		else
		{
			std::cerr<<"Error!! check\n";
			_exit(-1);
		}
	}
	else 
	{
		TGraph::vertex_descriptor vd1 = source(ei,g);
		TGraph::vertex_descriptor vd2 = target(ei,g);
		if(g[vd1].special == 1)
		{
			g[ei].type = APPEAR;
			return APPEAR;
		}
		if(g[vd2].special == 1)
		{
			g[ei].type = DISAPPEAR;
			return DISAPPEAR;
		}
		g[ei].type = TRANSLATION;
		return TRANSLATION;
	}
}


FeaturesType predict_feature(FeaturesType f1, FeaturesType f2)
{
	FeaturesType f3;
	f3.Centroid[0] = f2.Centroid[0] + f2.Centroid[0] - f1.Centroid[0];
	f3.Centroid[1] = f2.Centroid[1] + f2.Centroid[1] - f1.Centroid[1];
	f3.Centroid[2] = f2.Centroid[2] + f2.Centroid[2] - f1.Centroid[2];

	f3.BoundingBox[0] = 2*f2.BoundingBox[0] - f1.BoundingBox[0]; 
	f3.BoundingBox[1] = 2*f2.BoundingBox[1] - f1.BoundingBox[1]; 
	f3.BoundingBox[2] = 2*f2.BoundingBox[2] - f1.BoundingBox[2]; 
	f3.BoundingBox[3] = 2*f2.BoundingBox[3] - f1.BoundingBox[3]; 
	f3.BoundingBox[4] = 2*f2.BoundingBox[4] - f1.BoundingBox[4]; 
	f3.BoundingBox[5] = 2*f2.BoundingBox[5] - f1.BoundingBox[5]; 
	for(int counter=0; counter< FeaturesType::N; counter++)
	{
		f3.ScalarFeatures[counter] = (f2.ScalarFeatures[counter] + f1.ScalarFeatures[counter])/2.0;
	}
	f3.time = f1.time;
	return f3;
}

float CellTracker::compute_LRUtility(FeaturesType f1, FeaturesType f2, FeaturesType f3)
{
	float utility = 0;
	utility += compute_normal_utility(f1,f2);
	utility += compute_normal_utility(f2,f3);


	


	if(1)
	{
		float cencentroid[3];
		cencentroid[0] = (f1.Centroid[0] + f3.Centroid[0])/2.0;
		cencentroid[1] = (f1.Centroid[1] + f3.Centroid[1])/2.0;
		cencentroid[2] = (f1.Centroid[2] + f3.Centroid[2])/2.0;

		float dist = get_distance(f2.Centroid,cencentroid);
		utility *= exp(-dist*dist/2/fvar.distVariance);
		return utility;
	}
	else
	{
		float factor;
		float dist1 = get_distance(f1.Centroid,f2.Centroid);
		float dist2 = get_distance(f2.Centroid,f3.Centroid);
		float radius1 = pow(3.0/4.0/3.141*f1.ScalarFeatures[FeaturesType::VOLUME],1/3.0);
		float radius2 = pow(3.0/4.0/3.141*f2.ScalarFeatures[FeaturesType::VOLUME],1/3.0);
		if(dist1<radius1 || dist2 < radius2)
		{
			return utility;
		}
		else
		{


			float dir1[3],dir2[3];
			dir1[0] = f2.Centroid[0]-f1.Centroid[0];
			dir1[1] = f2.Centroid[1]-f1.Centroid[1];
			dir1[2] = f2.Centroid[2]-f1.Centroid[2];
			dir2[0] = f3.Centroid[0]-f2.Centroid[0];
			dir2[1] = f3.Centroid[1]-f2.Centroid[1];
			dir2[2] = f3.Centroid[2]-f2.Centroid[2];
			float sum1 = sqrt(dir1[0]*dir1[0]+dir1[1]*dir1[1]+dir1[2]*dir1[2]);
			float sum2 = sqrt(dir2[0]*dir2[0]+dir2[1]*dir2[1]+dir2[2]*dir2[2]);
			if (sum1 < 1e-3)
				sum1 = 1e-3;
			if (sum2 < 1e-3)
				sum2 = 1e-3;
			dir1[0] /=sum1;
			dir1[1] /=sum1;
			dir1[2] /=sum1;
			dir2[0] /=sum2;
			dir2[1] /=sum2;
			dir2[2] /=sum2;

			float proj = (1+dir1[0]*dir2[0]+dir1[1]*dir2[1]+dir1[2]*dir2[2])/2;

			float utility1 = utility*proj ;//* MIN(dist1,dist2)/MAX(dist1,dist2)*proj;

			//float test[3] = { 450,48.0,6.0};
			//if(get_distance(test, f1.Centroid) < 20 && (f1.time == 19))// || f1.time==20 || f1.time == 21 || f1.time==18))
			//{
			//	printFeatures(f1);
			//	printFeatures(f2);
			//	printFeatures(f3);
			//	printf("utility = %0.0f proj = %0.4f\n",utility,proj);
			//	//PAUSE;
			//}
			return utility1;
		}
	}

}

float CellTracker::get_misalignment_cost(FeaturesType f1, FeaturesType f2, FeaturesType f3)
{
	float cencentroid[3];
	cencentroid[0] = (f1.Centroid[0] + f3.Centroid[0])/2.0;
	cencentroid[1] = (f1.Centroid[1] + f3.Centroid[1])/2.0;
	cencentroid[2] = (f1.Centroid[2] + f3.Centroid[2])/2.0;

	float dist = get_distance(f2.Centroid,cencentroid);
	float reduction = exp(-dist*dist/2/fvar.distVariance);
	
	return reduction;
}


FeaturesType average_feature(FeaturesType f1,FeaturesType f2)
{
	FeaturesType f;
	for(int counter = 0; counter< FeaturesType::N; counter++)
	{
		f.ScalarFeatures[counter] = (f1.ScalarFeatures[counter]+f2.ScalarFeatures[counter])/2.0;
	}
	f.time = f1.time;
	f.Centroid[0] = (f1.Centroid[0] + f2.Centroid[0])/2.0;
	f.Centroid[1] = (f1.Centroid[1] + f2.Centroid[1])/2.0;
	f.Centroid[2] = (f1.Centroid[2] + f2.Centroid[2])/2.0;

	f.BoundingBox[0] = (f1.BoundingBox[0] + f2.BoundingBox[0])/2.0;
	f.BoundingBox[1] = (f1.BoundingBox[1] + f2.BoundingBox[1])/2.0;
	f.BoundingBox[2] = (f1.BoundingBox[2] + f2.BoundingBox[2])/2.0;
	f.BoundingBox[3] = (f1.BoundingBox[3] + f2.BoundingBox[3])/2.0;
	f.BoundingBox[4] = (f1.BoundingBox[4] + f2.BoundingBox[4])/2.0;
	f.BoundingBox[5] = (f1.BoundingBox[5] + f2.BoundingBox[5])/2.0;
	return f;
}
float CellTracker::compute_LRUtility(TGraph::edge_descriptor e1,TGraph::edge_descriptor e2)
{
	float utility = 0;
	utility = (g[e1].utility + g[e2].utility);
	if(get_edge_type(e1) == APPEAR || get_edge_type(e2) == DISAPPEAR)
	{
		printf("a");
		if(get_edge_type(e1)==APPEAR && get_edge_type(e2)!=DISAPPEAR)
		{
			FeaturesType f1, f2, f3;
			f2 = fvector[g[source(e2,g)].t][g[source(e2,g)].findex];
			f3 = fvector[g[target(e2,g)].t][g[target(e2,g)].findex];
			
			if(get_edge_type(e2)==SPLIT)
			{
				printf("c");
				FeaturesType f1a = predict_feature(f3,f2);
				FeaturesType f1b = predict_feature(fvector[g[target(coupled_map[e2],g)].t][g[target(coupled_map[e2],g)].findex],f2);
				f1 = average_feature(f1a,f1b);
				utility = g[e2].utility;
			}
			else
			{
				printf("d");
				f1 = predict_feature(f3, f2);
				if(get_edge_type(e2)==MERGE)
					utility = compute_normal_utility(f2,f3)/fvar.T_prior*fvar.MS_prior;
				else
					utility = compute_normal_utility(f2,f3);
			}

			float dist = get_boundary_dist(f1.Centroid)-fvar.spacing[2];
			if(dist <0)
			{
				utility += (1-exp(-dist*dist/2/fvar.distVariance))*compute_normal_utility(f1,f2)/fvar.T_prior*fvar.AD_prior;
			}
			else
			{
				utility = 10;
			}
			/*utility += compute_normal_utility(f1,f2);

			if(dist < 0)
				utility *= (1-exp(-dist*dist/2/fvar.distVariance))/10;
			else
				utility = 10;*/
			
			return utility;
			
		}
		else if(get_edge_type(e1)!=APPEAR && get_edge_type(e2)==DISAPPEAR)
		{
			printf("b");
			FeaturesType f1, f2,f3;
			f1 = fvector[g[source(e1,g)].t][g[source(e1,g)].findex];
			f2 = fvector[g[target(e1,g)].t][g[target(e1,g)].findex];

			if(get_edge_type(e1)==MERGE)
			{
				printf("c");
				FeaturesType f3a = predict_feature(f1,f2);
				FeaturesType f3b = predict_feature(fvector[g[target(coupled_map[e1],g)].t][g[target(coupled_map[e1],g)].findex],f2);
				f3 = average_feature(f3a,f3b);
				utility = g[e1].utility;
			}
			else
			{
				printf("d");
				f3 = predict_feature(f1, f2);
				if(get_edge_type(e1)==SPLIT)
					utility = compute_normal_utility(f1,f2)/fvar.T_prior*fvar.MS_prior;
				else
					utility = compute_normal_utility(f1,f2);
			}
			float dist = get_boundary_dist(f3.Centroid)-fvar.spacing[2];

			if(dist<0)
			{
					utility += (1-exp(-dist*dist/2/fvar.distVariance))*compute_normal_utility(f2,f3)/fvar.T_prior*fvar.AD_prior;
			}
			else
			{
				return 10;
			}
			/*
			utility += compute_normal_utility(f2,f3);
			if(dist < 0)
				utility *= (1-exp(-dist*dist/2/fvar.distVariance))/10;
			else
				utility = 10;
				*/
			
			return utility;
		}
		return utility;
	}
	/*FeaturesType f1,f2,f3;
	f1 = fvector[g[source(e1,g)].t][g[source(e1,g)].findex];
	f2 = fvector[g[source(e2,g)].t][g[source(e2,g)].findex];
	f3 = fvector[g[target(e2,g)].t][g[target(e2,g)].findex];*/
	FeaturesType f1a, f1b,f2, f3a,f3b;
	
	int e1_type = get_edge_type(e1);
	int e2_type = get_edge_type(e2);
	float reduction;
	if(e1_type == MERGE)
	{
		if(e2_type == SPLIT)
		{
			f1a = fvector[g[source(e1,g)].t][g[source(e1,g)].findex];
			f1b = fvector[g[source(coupled_map[e1],g)].t][g[source(coupled_map[e1],g)].findex];
			f2 = fvector[g[target(e1,g)].t][g[target(e1,g)].findex];
			f3a = fvector[g[target(e2,g)].t][g[target(e2,g)].findex];
			f3b = fvector[g[source(coupled_map[e2],g)].t][g[source(coupled_map[e2],g)].findex];
			reduction = MAX(get_misalignment_cost(f1a,f2,f3a)+get_misalignment_cost(f1b,f2,f3b),get_misalignment_cost(f1a,f2,f3b)+get_misalignment_cost(f1b,f2,f3a));
			utility = (g[e1].utility + g[e2].utility);
			utility *=reduction/2.0;
			return utility;
		}
		else
		{
			f1a = fvector[g[source(e1,g)].t][g[source(e1,g)].findex];
			f1b = fvector[g[source(coupled_map[e1],g)].t][g[source(coupled_map[e1],g)].findex];
			f2 = fvector[g[target(e1,g)].t][g[target(e1,g)].findex];
			f3a = fvector[g[target(e2,g)].t][g[target(e2,g)].findex];
			if(e2_type==MERGE)
				utility = g[e1].utility + compute_normal_utility(f2,f3a)/fvar.T_prior*fvar.MS_prior;
			else
				utility = g[e1].utility + compute_normal_utility(f2,f3a);
			reduction = get_misalignment_cost(f1a,f2,f3a)+ get_misalignment_cost(f1b,f2,f3a);
			utility *=reduction/2;
			return utility;
		}
	}
	else
	{
		if(e2_type==SPLIT)
		{
			f1a = fvector[g[source(e1,g)].t][g[source(e1,g)].findex];
			f2 = fvector[g[target(e1,g)].t][g[target(e1,g)].findex];
			f3a = fvector[g[target(e2,g)].t][g[target(e2,g)].findex];
			f3b = fvector[g[source(coupled_map[e2],g)].t][g[source(coupled_map[e2],g)].findex];
			if(e1_type==MERGE)
				utility = g[e2].utility + compute_normal_utility(f1a,f2)/fvar.T_prior*fvar.MS_prior;
			else
				utility = g[e2].utility + compute_normal_utility(f1a,f2);
			reduction = get_misalignment_cost(f1a,f2,f3a)+ get_misalignment_cost(f1a,f2,f3b);
			utility *=reduction/2;
			return utility;
		}
		else
		{
			f1a = fvector[g[source(e1,g)].t][g[source(e1,g)].findex];
			f2 = fvector[g[target(e1,g)].t][g[target(e1,g)].findex];
			f3a = fvector[g[target(e2,g)].t][g[target(e2,g)].findex];
			if(e2_type==MERGE)
			{
				//utility =  compute_normal_utility(f2,f3a)/fvar.T_prior*fvar.MS_prior;
				utility =  g[e2].utility/2;
			}
			else
			{
				utility =  compute_normal_utility(f2,f3a);
			}
			if(e1_type==SPLIT)
			{
				//utility += compute_normal_utility(f1a,f2)/fvar.T_prior*fvar.MS_prior;
				utility += g[e1].utility/2;
			}
			else
			{
				utility += compute_normal_utility(f1a,f2);
			}
			reduction = get_misalignment_cost(f1a,f2,f3a);
			utility *=reduction;
			return utility;
		}
	}
}
float CellTracker::compute_LRUtility_sum(TGraph::edge_descriptor e1,TGraph::edge_descriptor e2)
{
	float utility = 0;
	utility = -2*UTILITY_MAX;//(g[e1].utility * g[e2].utility);
	if(g[e1].utility < 0 || g[e2].utility <0)
	{
		//printf("I'm getting negative utilities here\n");
		//scanf("%*d");
	}

	
	if(get_edge_type(e1) == APPEAR || get_edge_type(e2) == DISAPPEAR)
	{
		bool special_case = false;
		if(g[source(e1,g)].t == -1 || g[target(e2,g)].t == fvar.time_last)
		{
			special_case = true;
		}
		printf("a");
		if(get_edge_type(e1)==APPEAR && get_edge_type(e2)!=DISAPPEAR)
		{
			FeaturesType f1, f2, f3;
			f2 = fvector[g[source(e2,g)].t][g[source(e2,g)].findex];
			f3 = fvector[g[target(e2,g)].t][g[target(e2,g)].findex];
			
			if(get_edge_type(e2)==SPLIT)
			{
				printf("c");
				FeaturesType f1a = predict_feature(f3,f2);
				FeaturesType f1b = predict_feature(fvector[g[target(coupled_map[e2],g)].t][g[target(coupled_map[e2],g)].findex],f2);
				f1 = average_feature(f1a,f1b);
				utility = g[e2].utility;
			}
			else
			{
				printf("d");
				f1 = predict_feature(f3, f2);
				if(get_edge_type(e2)==MERGE)
					utility = g[e2].utility/2; //compute_normal_utility(f2,f3)/fvar.T_prior*fvar.MS_prior;
				else
					utility = compute_normal_utility(f2,f3);
			}

			float dist = get_boundary_dist(f1.Centroid)-fvar.spacing[2];
			if(dist <0)
			{
				utility += (1-exp(-dist*dist/2/fvar.distVariance))*compute_normal_utility(f1,f2)/fvar.T_prior*fvar.AD_prior;
			}
			else
			{
				if(!special_case)
					utility +=-UTILITY_MAX;
				else
					utility +=0;
			}
			/*utility += compute_normal_utility(f1,f2);

			if(dist < 0)
				utility *= (1-exp(-dist*dist/2/fvar.distVariance))/10;
			else
				utility = 10;*/
			
			return utility;
			
		}
		else if(get_edge_type(e1)!=APPEAR && get_edge_type(e2)==DISAPPEAR)
		{
			printf("b");
			FeaturesType f1, f2,f3;
			f1 = fvector[g[source(e1,g)].t][g[source(e1,g)].findex];
			f2 = fvector[g[target(e1,g)].t][g[target(e1,g)].findex];

			if(get_edge_type(e1)==MERGE)
			{
				printf("c");
				FeaturesType f3a = predict_feature(f1,f2);
				FeaturesType f3b = predict_feature(fvector[g[target(coupled_map[e1],g)].t][g[target(coupled_map[e1],g)].findex],f2);
				f3 = average_feature(f3a,f3b);
				utility = g[e1].utility;
			}
			else
			{
				printf("d");
				f3 = predict_feature(f1, f2);
				if(get_edge_type(e1)==SPLIT)
					utility = g[e1].utility/2;//compute_normal_utility(f1,f2)/fvar.T_prior*fvar.MS_prior;
				else
					utility = compute_normal_utility(f1,f2);
			}
			float dist = get_boundary_dist(f3.Centroid)-fvar.spacing[2];

			if(dist<0)
			{
					utility += (1-exp(-dist*dist/2/fvar.distVariance))*compute_normal_utility(f2,f3)/fvar.T_prior*fvar.AD_prior;
			}
			else
			{
				if(!special_case)
					utility += -UTILITY_MAX;
				else
					utility +=0;
			}
			/*
			utility += compute_normal_utility(f2,f3);
			if(dist < 0)
				utility *= (1-exp(-dist*dist/2/fvar.distVariance))/10;
			else
				utility = 10;
				*/
			
			return utility;
		}
		return utility;
	}
	/*FeaturesType f1,f2,f3;
	f1 = fvector[g[source(e1,g)].t][g[source(e1,g)].findex];
	f2 = fvector[g[source(e2,g)].t][g[source(e2,g)].findex];
	f3 = fvector[g[target(e2,g)].t][g[target(e2,g)].findex];*/
	FeaturesType f1a, f1b,f2, f3a,f3b;
	
	int e1_type = get_edge_type(e1);
	int e2_type = get_edge_type(e2);
	float reduction;
	if(e1_type == MERGE)
	{
		if(e2_type == SPLIT)
		{
			f1a = fvector[g[source(e1,g)].t][g[source(e1,g)].findex];
			f1b = fvector[g[source(coupled_map[e1],g)].t][g[source(coupled_map[e1],g)].findex];
			f2 = fvector[g[target(e1,g)].t][g[target(e1,g)].findex];
			f3a = fvector[g[target(e2,g)].t][g[target(e2,g)].findex];
			f3b = fvector[g[source(coupled_map[e2],g)].t][g[source(coupled_map[e2],g)].findex];
			reduction = MAX(get_misalignment_cost(f1a,f2,f3a)+get_misalignment_cost(f1b,f2,f3b),get_misalignment_cost(f1a,f2,f3b)+get_misalignment_cost(f1b,f2,f3a));
			utility = (g[e1].utility + g[e2].utility);
			utility *=reduction/2.0;
			return utility;
		}
		else
		{
			f1a = fvector[g[source(e1,g)].t][g[source(e1,g)].findex];
			f1b = fvector[g[source(coupled_map[e1],g)].t][g[source(coupled_map[e1],g)].findex];
			f2 = fvector[g[target(e1,g)].t][g[target(e1,g)].findex];
			f3a = fvector[g[target(e2,g)].t][g[target(e2,g)].findex];
			if(e2_type==MERGE)
				utility = g[e1].utility + g[e2].utility/2;// compute_normal_utility(f2,f3a)/fvar.T_prior*fvar.MS_prior;
			else
				utility = g[e1].utility +g[e2].utility;// compute_normal_utility(f2,f3a);
			reduction = get_misalignment_cost(f1a,f2,f3a)+ get_misalignment_cost(f1b,f2,f3a);
			utility *=reduction/2;
			return utility;
		}
	}
	else
	{
		if(e2_type==SPLIT)
		{
			f1a = fvector[g[source(e1,g)].t][g[source(e1,g)].findex];
			f2 = fvector[g[target(e1,g)].t][g[target(e1,g)].findex];
			f3a = fvector[g[target(e2,g)].t][g[target(e2,g)].findex];
			f3b = fvector[g[source(coupled_map[e2],g)].t][g[source(coupled_map[e2],g)].findex];
			if(e1_type==SPLIT)
				utility = g[e2].utility + g[e1].utility/2;// compute_normal_utility(f1a,f2)/fvar.T_prior*fvar.MS_prior;
			else
				utility = g[e2].utility + g[e1].utility;// compute_normal_utility(f1a,f2);
			reduction = get_misalignment_cost(f1a,f2,f3a)+ get_misalignment_cost(f1a,f2,f3b);
			utility *=reduction/2;
			return utility;
		}
		else
		{
			f1a = fvector[g[source(e1,g)].t][g[source(e1,g)].findex];
			f2 = fvector[g[target(e1,g)].t][g[target(e1,g)].findex];
			f3a = fvector[g[target(e2,g)].t][g[target(e2,g)].findex];
			if(e2_type==MERGE)
			{
				//utility =  compute_normal_utility(f2,f3a)/fvar.T_prior*fvar.MS_prior;
				utility =  g[e2].utility/2;
			}
			else
			{
				utility =  g[e2].utility;//compute_normal_utility(f2,f3a);
			}
			if(e1_type==SPLIT)
			{
				//utility += compute_normal_utility(f1a,f2)/fvar.T_prior*fvar.MS_prior;
				utility += g[e1].utility/2;
			}
			else
			{
				utility += g[e1].utility;//compute_normal_utility(f1a,f2);
			}
			reduction = get_misalignment_cost(f1a,f2,f3a);
			utility *=reduction;
			return utility;
		}
	}
}
float CellTracker::compute_LRUtility_product(TGraph::edge_descriptor e1,TGraph::edge_descriptor e2)
{
	float utility = 0;
	utility = -4*UTILITY_MAX;//(g[e1].utility * g[e2].utility);
	if(g[e1].utility < 0 || g[e2].utility <0)
	{
		printf("I'm getting negative utilities here\n");
		scanf("%*d");
	}

	
	if(get_edge_type(e1) == APPEAR || get_edge_type(e2) == DISAPPEAR)
	{
		bool special_case = false;
		if(g[source(e1,g)].t == -1 || g[target(e2,g)].t == fvar.time_last+1)
		{
			special_case = true;
			utility = -2*UTILITY_MAX;
		}
		printf("a");
		if(get_edge_type(e1)==APPEAR && get_edge_type(e2)!=DISAPPEAR)
		{
			FeaturesType f1, f2, f3;
			f2 = fvector[g[source(e2,g)].t][g[source(e2,g)].findex];
			f3 = fvector[g[target(e2,g)].t][g[target(e2,g)].findex];
			
			if(get_edge_type(e2)==SPLIT)
			{
				printf("c");
				FeaturesType f1a = predict_feature(f3,f2);
				FeaturesType f1b = predict_feature(fvector[g[target(coupled_map[e2],g)].t][g[target(coupled_map[e2],g)].findex],f2);
				f1 = average_feature(f1a,f1b);
				utility = g[e2].utility;
			}
			else
			{
				printf("d");
				f1 = predict_feature(f3, f2);
				if(get_edge_type(e2)==MERGE)
					utility = g[e2].utility/2; //compute_normal_utility(f2,f3)/fvar.T_prior*fvar.MS_prior;
				else
					utility = g[e2].utility;
			}

			float dist = get_boundary_dist(f1.Centroid)-fvar.boundDistMean;
			if(dist <0)
			{
				utility *= (1-exp(-dist*dist/2/fvar.boundDistVariance))*compute_normal_utility(f1,f2)/fvar.T_prior*fvar.AD_prior;
			}
			else
			{
				if(!special_case)
					utility +=-UTILITY_MAX;
				else
					utility *=1;
			}
			/*utility += compute_normal_utility(f1,f2);

			if(dist < 0)
				utility *= (1-exp(-dist*dist/2/fvar.distVariance))/10;
			else
				utility = 10;*/
			
			return utility;
			
		}
		else if(get_edge_type(e1)!=APPEAR && get_edge_type(e2)==DISAPPEAR)
		{
			printf("b");
			FeaturesType f1, f2,f3;
			f1 = fvector[g[source(e1,g)].t][g[source(e1,g)].findex];
			f2 = fvector[g[target(e1,g)].t][g[target(e1,g)].findex];

			if(get_edge_type(e1)==MERGE)
			{
				printf("c");
				FeaturesType f3a = predict_feature(f1,f2);
				FeaturesType f3b = predict_feature(fvector[g[target(coupled_map[e1],g)].t][g[target(coupled_map[e1],g)].findex],f2);
				f3 = average_feature(f3a,f3b);
				utility = g[e1].utility;
			}
			else
			{
				printf("d");
				f3 = predict_feature(f1, f2);
				if(get_edge_type(e1)==SPLIT)
					utility = g[e1].utility/2;//compute_normal_utility(f1,f2)/fvar.T_prior*fvar.MS_prior;
				else
					utility = compute_normal_utility(f1,f2);
			}
			float dist = get_boundary_dist(f3.Centroid)-fvar.boundDistMean;

			if(dist<0)
			{
					utility *= (1-exp(-dist*dist/2/fvar.boundDistVariance))*compute_normal_utility(f2,f3)/fvar.T_prior*fvar.AD_prior;
			}
			else
			{
				if(!special_case)
					utility += -UTILITY_MAX;
				else
					utility *=1;
			}
			/*
			utility += compute_normal_utility(f2,f3);
			if(dist < 0)
				utility *= (1-exp(-dist*dist/2/fvar.distVariance))/10;
			else
				utility = 10;
				*/
			
			return utility;
		}
		return utility;
	}
	/*FeaturesType f1,f2,f3;
	f1 = fvector[g[source(e1,g)].t][g[source(e1,g)].findex];
	f2 = fvector[g[source(e2,g)].t][g[source(e2,g)].findex];
	f3 = fvector[g[target(e2,g)].t][g[target(e2,g)].findex];*/
	FeaturesType f1a, f1b,f2, f3a,f3b;
	
	int e1_type = get_edge_type(e1);
	int e2_type = get_edge_type(e2);
	float reduction;
	if(e1_type == MERGE)
	{
		if(e2_type == SPLIT)
		{
			f1a = fvector[g[source(e1,g)].t][g[source(e1,g)].findex];
			f1b = fvector[g[source(coupled_map[e1],g)].t][g[source(coupled_map[e1],g)].findex];
			f2 = fvector[g[target(e1,g)].t][g[target(e1,g)].findex];
			f3a = fvector[g[target(e2,g)].t][g[target(e2,g)].findex];
			f3b = fvector[g[source(coupled_map[e2],g)].t][g[source(coupled_map[e2],g)].findex];
			reduction = MAX(get_misalignment_cost(f1a,f2,f3a)+get_misalignment_cost(f1b,f2,f3b),get_misalignment_cost(f1a,f2,f3b)+get_misalignment_cost(f1b,f2,f3a));
			utility = (g[e1].utility * g[e2].utility);
			utility *=reduction/2.0;
			return utility;
		}
		else
		{
			f1a = fvector[g[source(e1,g)].t][g[source(e1,g)].findex];
			f1b = fvector[g[source(coupled_map[e1],g)].t][g[source(coupled_map[e1],g)].findex];
			f2 = fvector[g[target(e1,g)].t][g[target(e1,g)].findex];
			f3a = fvector[g[target(e2,g)].t][g[target(e2,g)].findex];
			if(e2_type==MERGE)
				utility = g[e1].utility * g[e2].utility/2;// compute_normal_utility(f2,f3a)/fvar.T_prior*fvar.MS_prior;
			else
				utility = g[e1].utility *g[e2].utility;// compute_normal_utility(f2,f3a);
			reduction = get_misalignment_cost(f1a,f2,f3a)+ get_misalignment_cost(f1b,f2,f3a);
			utility *=reduction/2;
			return utility;
		}
	}
	else
	{
		if(e2_type==SPLIT)
		{
			f1a = fvector[g[source(e1,g)].t][g[source(e1,g)].findex];
			f2 = fvector[g[target(e1,g)].t][g[target(e1,g)].findex];
			f3a = fvector[g[target(e2,g)].t][g[target(e2,g)].findex];
			f3b = fvector[g[source(coupled_map[e2],g)].t][g[source(coupled_map[e2],g)].findex];
			if(e1_type==SPLIT)
				utility = g[e2].utility * g[e1].utility/2;// compute_normal_utility(f1a,f2)/fvar.T_prior*fvar.MS_prior;
			else
				utility = g[e2].utility * g[e1].utility;// compute_normal_utility(f1a,f2);
			reduction = get_misalignment_cost(f1a,f2,f3a)+ get_misalignment_cost(f1a,f2,f3b);
			utility *=reduction/2;
			return utility;
		}
		else
		{
			f1a = fvector[g[source(e1,g)].t][g[source(e1,g)].findex];
			f2 = fvector[g[target(e1,g)].t][g[target(e1,g)].findex];
			f3a = fvector[g[target(e2,g)].t][g[target(e2,g)].findex];
			if(e2_type==MERGE)
			{
				//utility =  compute_normal_utility(f2,f3a)/fvar.T_prior*fvar.MS_prior;
				utility =  g[e2].utility/2;
			}
			else
			{
				utility =  g[e2].utility;//compute_normal_utility(f2,f3a);
			}
			if(e1_type==SPLIT)
			{
				//utility += compute_normal_utility(f1a,f2)/fvar.T_prior*fvar.MS_prior;
				utility *= g[e1].utility/2;
			}
			else
			{
				utility *= g[e1].utility;//compute_normal_utility(f1a,f2);
			}
			reduction = get_misalignment_cost(f1a,f2,f3a);
			utility *=reduction;
			return utility;
		}
	}
}

void CellTracker::draw_line_for_edge(int num, TGraph::edge_descriptor e,VectorPixelType col1,VectorPixelType col2, int shift = 0)
{
	TGraph::vertex_descriptor v1,v2;
	//printf("g[e].utility = %d\n",g[e].utility);
	//if(g[e].utility < 0)
	//{
	//	return;
	//}
	float strength = 1.0;//pow(double(g[e].utility*1.0/UTILITY_MAX),double(0.1));
	VectorPixelType color1,color2;
	color1[0] = col1[0]*strength;color1[1] = col1[1]*strength;color1[2] = col1[2]*strength;
	color2[0] = col2[0]*strength;color2[1] = col2[1]*strength;color2[2] = col2[2]*strength;
	v1 = source(e,g);
	v2 = target(e,g);

	if(1)//g[e].selected==1)
	{
		if(g[v1].special ==0 && g[v2].special ==0)
		{
			FeaturesType f1 = fvector[g[v1].t][g[v1].findex];
			FeaturesType f2 = fvector[g[v2].t][g[v2].findex];
			//if(g[v1].findex==0 && g[v1].t==4)
			//{
			//	printf("This edge is still there %d %d %d %d\n",g[v1].t,g[v1].findex,g[v2].t,g[v2].findex);
			//	scanf("%*d");
			//}
			if(get_distance(f1.Centroid,f2.Centroid)>100)
			{
				printf("source\n");
				print_vertex(v1,1);
				printf("target\n");
				print_vertex(v2,1);
				scanf("%*d");
			}
			//printf("col1 col2 %d %d %d %d %d %d\n",col1[0],col1[1],col1[2],col2[0],col2[1],col2[2]);
			//printf("color1 color2 %d %d %d %d %d %d\n",color1[0],color1[1],color1[2],color2[0],color2[1],color2[2]);
			//PAUSE;
			if(num == 1)
				drawLine(debugimage1,color1,color2,f1.Centroid[0]+shift,f1.Centroid[1],f1.time,f2.Centroid[0]+shift,f2.Centroid[1],f1.time);
			else if(num==2)
				drawLine(debugimage2,color1,color2,f1.Centroid[0]+shift,f1.Centroid[1],f1.time,f2.Centroid[0]+shift,f2.Centroid[1],f1.time);
			else
				drawLine(debugimage3,color1,color2,f1.Centroid[0]+shift,f1.Centroid[1],f1.time,f2.Centroid[0]+shift,f2.Centroid[1],f1.time);
			//drawLine(debugimage,col2,f1.Centroid[0]+shift,f1.Centroid[1],f2.time,f2.Centroid[0]+shift,f2.Centroid[1],f2.time);
		}
		else
		{
			ColorImageType::PixelType pixel;
			pixel[0] = 255;
			pixel[1] = 0;
			pixel[2] = 0;
			FeaturesType f1;
			int nshift = 0;
			if(get_edge_type(e)==APPEAR)
			{
				nshift = 2;
				pixel[0] = 0;
				pixel[1] = 255;
				pixel[2] = 0;
			}
			if(get_edge_type(e)==DISAPPEAR)
				nshift = -2;
			if(g[v1].special == 0)
			{
				//printf("I came into g[v1].special==0\n");
				//PAUSE;
				f1 = fvector[g[v1].t][g[v1].findex];
			}
			else
			{
				f1 = fvector[g[v2].t][g[v2].findex];
			}
			ColorImageType::IndexType index1,index2;
			index1[0] = MAX(MIN(f1.Centroid[0]+nshift-1,fvar.BoundingBox[1]),fvar.BoundingBox[0]);
			index1[1] = MAX(MIN(f1.Centroid[1]+nshift-1,fvar.BoundingBox[3]),fvar.BoundingBox[2]);
			index1[2] = f1.time;
			index2[0] = MAX(MIN(f1.Centroid[0]+nshift+1,fvar.BoundingBox[1]),fvar.BoundingBox[0]);
			index2[1] = MAX(MIN(f1.Centroid[1]+nshift+1,fvar.BoundingBox[3]),fvar.BoundingBox[2]);
			index2[2] = f1.time;
			for(int cox = index1[0]; cox<=index2[0]; ++cox)
				for(int coy = index1[1];coy<=index2[1]; ++coy)
					for(int coz = index1[2];coz<=index2[2]; ++coz)
					{
						//printf("-\n");
						ColorImageType::IndexType index; index[0] = cox; index[1] = coy; index[2] = coz;
						if(num==1)
							debugimage1->SetPixel(index,pixel);
						else if(num==2)
							debugimage2->SetPixel(index,pixel);
						else
							debugimage3->SetPixel(index,pixel);
					}
		}
	}
}



bool CellTracker::edge_uniqueness( TGraph::edge_descriptor e1, TGraph::edge_descriptor e2)
{
	if(g[e1].coupled == 1)
	{
		if(coupled_map[e1]==e2)
		{
			return 1;
		}
		return 0;
	}
	return 0;
}
void CellTracker::solve_higher_order()
{

	VectorPixelType col1,col2,col3,col4;
	col1[0] = 255;col1[1] = 0;col1[2] = 0;
	col2[0] = 255;col2[1] = 255;col2[2] = 0;
	col3[0] = 255;col3[1] = 255;col3[2] = 255;
	col4[0] = 255;col4[1] = 0; col4[2] = 255;



	boost::property_map<TGraph, boost::vertex_index_t>::type index;
	index = get(boost::vertex_index,g);
	//TODO - FIXME - work in progress - dont use it yet.
	IloEnv env;
	try{
		using boost::graph_traits;
		graph_traits<TGraph>::edge_iterator ei,eend;
//		std::map<TGraph::edge_descriptor,int> var_index;
		std::vector<int> utility;

		std::vector<LREdge> lredges;

		std::multimap<TGraph::edge_descriptor, int> frontmap;
		std::multimap<TGraph::edge_descriptor, int> backmap;
		std::multimap<TGraph::vertex_descriptor, int> vertmap;
		typedef std::pair<TGraph::edge_descriptor, int> EdgeMapType;
		typedef std::pair<TGraph::vertex_descriptor, int> VertexMapType;

		IloObjective obj = IloMaximize(env);
		IloRangeArray c(env);

		//find number of variables
		tie(ei,eend) = boost::edges(g);
		int ecount = 0;
		int varc = 0;

		
		printf("ecount = %d varc = %d\n",ecount,varc);
		graph_traits<TGraph>::vertex_iterator vi,vend;

		int num_v = 0;
		int in_deg_count = 0;
		int out_deg_count = 0;
		for(tie(vi,vend) = vertices(g); vi != vend; ++vi)
		{
			//printf("#");
			++num_v;
			in_deg_count += in_degree(*vi,g);
			out_deg_count += out_degree(*vi,g);

			if(in_degree(*vi,g) ==0 || out_degree(*vi,g)==0)
			{
				// no need to form any new LREdge with this vertex
				continue;
			}
			graph_traits<TGraph>::in_edge_iterator e_in,e_in_end;
			graph_traits<TGraph>::out_edge_iterator e_out,e_out_end;

			tie(e_in,e_in_end) = in_edges(*vi,g);
			tie(e_out,e_out_end) = out_edges(*vi,g);
			

			std::set<TGraph::edge_descriptor> in_unique_set;
			std::set<TGraph::edge_descriptor> out_unique_set;
			in_unique_set.clear();
			out_unique_set.clear();
			for(;e_in != e_in_end;++e_in)
			{
				if(g[*e_in].coupled == 0)
				{
					in_unique_set.insert(*e_in);
				}
				else
				{
					if(in_unique_set.count(coupled_map[*e_in])==0)
					{
						in_unique_set.insert(*e_in);
					}
				}
			}
			for(;e_out != e_out_end; ++e_out)
			{
				int type = get_edge_type(*e_out);
				if(g[*e_out].coupled == 0)
				{
					out_unique_set.insert(*e_out);
				}
				else
				{
					if(out_unique_set.count(coupled_map[*e_out])==0)
					{
						out_unique_set.insert(*e_out);
					}
				}
			}

			std::set<TGraph::edge_descriptor>::iterator i1,i2;
			for(i1 = in_unique_set.begin(); i1!= in_unique_set.end(); ++i1)
			{
				for(i2 = out_unique_set.begin(); i2!= out_unique_set.end(); ++i2)
				{
					LREdge lre;
					lre.front = *i2;
					lre.back = *i1;
					lredges.push_back(lre);
					//printf("-");
					utility.push_back(compute_LRUtility_product(lre.back,lre.front));
					//printf("|");
					g[*i2].frontlre.push_back(lredges.size()-1);
					g[*i1].backlre.push_back(lredges.size()-1);
					g[*vi].vertlre.push_back(lredges.size()-1);

					//frontmap.insert(EdgeMapType(*i2,lredges.size()-1));
					//backmap.insert(EdgeMapType(*i1,lredges.size()-1));
					//vertmap.insert(VertexMapType(*vi,lredges.size()-1));
				}
			}	
		}


		
		printf("lredges.size() = %d\n", lredges.size(),utility);
		printf("num_v = %d avg_in_degree = %0.2f avg_out_degree = %0.2f\n", num_v, in_deg_count*1.0/num_v, out_deg_count*1.0/num_v);
		//PAUSE;
		//_exit(0);
		varc = lredges.size();
		IloNumVarArray x(env,varc,0,1,ILOBOOL);
		IloNumArray numarr(env,varc);

		//IloNumExprArg expr(env);
		printf("I have changed\n");
		for(int counter=0; counter < varc; counter++)
		{
			if(utility[counter] <0)
			{
				;
				//printf("I do have <0 utility %d\n",utility[counter]);
				//scanf("%*d");
			}
			numarr[counter] = utility[counter];
			//obj.setLinearCoef(x[counter],utility[counter]);
		}
		//obj.setLinearCoef(x[counter],utility[counter]);
		obj.setLinearCoefs(x,numarr);
			//if(index[source(lredges[counter].front,g)]==638)
			//{
			//	printf("Utility for 638 is %d\n", utility[counter]);
			////	PAUSE;
			//}
			/*if(utility[counter]==0)
			{
			//	printf("Utility is zero : %d\n",counter);
				if(!(get_edge_type(lredges[counter].front)==DISAPPEAR && get_edge_type(lredges[counter].back)==APPEAR))
				{
					;
					//printf("vertex index = %d\n", index[source(lredges[counter].front,g)]);
				//	PAUSE;
				}
				
			}*/
			
			//printf("utility[%d] = %d\n",counter,utility[counter]);
		//}
		printf("step1 completed\n");
		int vcount = -1;
		for(tie(vi,vend) = vertices(g); vi != vend; ++vi)
		{
			if(g[*vi].vertlre.size()>0)
			{
				vcount++;
				c.add(IloRange(env,0.99,1.01));
				for(int colre = 0; colre< g[*vi].vertlre.size(); colre++)
				{
					c[vcount].setLinearCoef(x[g[*vi].vertlre[colre]],1);
				}
			}
			/*if(vertmap.count(*vi)!=0)
			{
				vcount++;
				c.add(IloRange(env,1,1));
				typedef std::multimap<TGraph::vertex_descriptor,int>::const_iterator Type1;
				std::pair<Type1, Type1> itrange = vertmap.equal_range(*vi);
				for(;itrange.first!=itrange.second; ++itrange.first)
				{
					Type1 it1 = itrange.first;
					c[vcount].setLinearCoef(x[it1->second],1.0);
				}
			}*/
		}
		printf("Done adding vertex constraints\n");
		//for(int counter = 0; counter < varc; counter++)
		//{
		//	x[counter].setBounds(0,1);
		//}
		printf("vertex_constraints = %d\n",vcount);
		//PAUSE;
		for(tie(ei,eend) = edges(g); ei!=eend; ++ei)
		{

			/*int type = get_edge_type(*ei);
			if(type==SPLIT || type== MERGE)
			{
				if(frontmap.count(*ei)>0 || frontmap.count(coupled_map[*ei]
			}*/
			//if(frontmap.count(*ei)>0 && backmap.count(*ei) >0)
			if(g[*ei].frontlre.size()>0 && g[*ei].backlre.size()>0)
			{
				int forward_coeff = 1;
				int backward_coeff = 1;
				int type = get_edge_type(*ei);
				if(type==SPLIT)
					backward_coeff = 2;
				if(type==MERGE)
					forward_coeff = 2;
				vcount++;
				c.add(IloRange(env,0,0));
				for(int colre = 0; colre<g[*ei].frontlre.size(); colre++)
				{
					c[vcount].setLinearCoef(x[g[*ei].frontlre[colre]],-backward_coeff);
				}
				if(type==MERGE)
				{
					TGraph::edge_descriptor ecoupled = coupled_map[*ei];
					for(int colre = 0; colre<g[ecoupled].frontlre.size(); colre++)
					{
						c[vcount].setLinearCoef(x[g[ecoupled].frontlre[colre]],-backward_coeff);
					}
				}
				for(int colre = 0; colre<g[*ei].backlre.size(); colre++)
				{
					c[vcount].setLinearCoef(x[g[*ei].backlre[colre]],forward_coeff);
				}
				if(type==SPLIT)
				{
					TGraph::edge_descriptor ecoupled = coupled_map[*ei];
					for(int colre = 0; colre<g[ecoupled].backlre.size(); colre++)
					{
						c[vcount].setLinearCoef(x[g[ecoupled].backlre[colre]],forward_coeff);
					}
				}

				/*typedef std::multimap<TGraph::edge_descriptor,int>::iterator Type1;
				std::pair<Type1,Type1> itrange = frontmap.equal_range(*ei);
				for(;itrange.first!=itrange.second; ++itrange.first)
				{
					Type1 it1 = itrange.first;
					c[vcount].setLinearCoef(x[it1->second],-backward_coeff);
				}
				if(type==MERGE)
				{
					itrange = frontmap.equal_range(coupled_map[*ei]);
					for(;itrange.first!=itrange.second; ++itrange.first)
					{
						Type1 it1 = itrange.first;
						c[vcount].setLinearCoef(x[it1->second],-backward_coeff);
					}
				}
				itrange = backmap.equal_range(*ei);
				for(;itrange.first!=itrange.second; ++itrange.first)
				{
					Type1 it1 = itrange.first;
					c[vcount].setLinearCoef(x[it1->second],forward_coeff);
				}
				if(type==SPLIT)
				{
					itrange = backmap.equal_range(coupled_map[*ei]);
					for(;itrange.first!=itrange.second; ++itrange.first)
					{
						Type1 it1 = itrange.first;
						c[vcount].setLinearCoef(x[it1->second],forward_coeff);
					}	
				}*/

			}
		}
		printf("Done adding edge constraints\n");

		printf("vcount = %d\n",vcount);
		////scanf("%*d");
		IloModel model(env);
		model.add(obj);
		model.add(c);

		//PAUSE;
		//std::cout<<model<<std::endl;
		IloCplex cplex(model);
		if(!cplex.solve())
		{
			std::cerr << " Could not solve.. error"<<std::endl;
			PAUSE;
		}
		IloNumArray vals(env);
		env.out() << "Solution status = " << cplex.getStatus() << std::endl;
		env.out() << "Solution value  = " << cplex.getObjValue() << std::endl;
		cplex.getValues(vals, x);
		/*for(int counter =0; counter < vals.getSize(); counter++)
		{
			if(vals[counter] < 0)
				printf("%0.2f ",vals[counter]);
		}*/
	//	printf("Are elements boolean? : %d\n", vals.areElementsBoolean());
	//	env.out() << "Values        = " << vals << std::endl;
		//std::cout << "Values.getSize() = "<< vals.getSize() << std::endl;
		int num_zero = 0;
		int num_one = 0;
		int num_others =0;
		writeXGMML_secondorder("C:\\Users\\arun\\Research\\Tracking\\harvard\\ch4_secondorder.xgmml",lredges,utility,vals);
		printf("varc =%d vals.getSize() = %d\n",varc,vals.getSize());
		//PAUSE;
		for(int counter=0; counter< vals.getSize(); counter++)
		{
			if(int(vals[counter]+0.5)==1)
				num_one++;
			else if(int(vals[counter]+0.5)==0)
				num_zero++;
			else
			{
				printf("vals[counter] = %f\n",vals[counter]);
				num_others++;
			}
			if(int(vals[counter]+0.5)==1)
			{
				g[lredges[counter].front].selected = 1;
				if(g[lredges[counter].front].coupled == 1)
				{
					g[coupled_map[lredges[counter].front]].selected = 1;
				}
				g[lredges[counter].back].selected = 1;
				if(g[lredges[counter].back].coupled == 1)
				{
					g[coupled_map[lredges[counter].back]].selected = 1;
				}
			}
		}
		printf("num_zero = %d num_one = %d num_others = %d\n",num_zero,num_one, num_others);
		//	//assert(g[inv_index[counter]].selected == 1);
		//	g[inv_index[counter]].selected = vals[counter];
		//	if(g[inv_index[counter]].coupled==1)
		//		g[coupled_map[inv_index[counter]]].selected = vals[counter];
		//	if(vals[counter]==1 && utility[counter] <0)
		//	{
		//		printf("wrong @ %d\n",counter);
		//	}
		//}
		//PAUSE;
		//_exit(0);
	}
	catch(IloException &e)
	{
		std::cerr << e << std::endl;
		PAUSE;
	}
	env.end();
}

std::pair<CellTracker::TGraph::edge_descriptor,bool> CellTracker::my_add_edge(TGraph::vertex_descriptor v1, TGraph::vertex_descriptor v2, int utility, bool coupled, bool fixed, bool selected, unsigned char type)
{
	printf("in my_add_edge\n");
	printf("v1-\n");print_vertex(v1,1);
	printf("v2-\n");print_vertex(v2,1);
	TGraph::edge_descriptor e;
	bool added;
	tie(e,added) = add_edge(v1,v2,g);
	if(added==0)
	{
		printf("Could not add edge...\n");
		scanf("%*d");
	}
	g[e].utility = utility;
	g[e].coupled = coupled;
	g[e].fixed = fixed;
	g[e].selected = selected;
	g[e].type = type;
	std::pair<TGraph::edge_descriptor,bool> pair1;
	pair1.first = e;
	pair1.second = added;
	return pair1;
}
void CellTracker::solve_lip()
{
	IloEnv env;
	try{
		using boost::graph_traits;
		graph_traits<TGraph>::edge_iterator ei,eend;
		std::map<TGraph::edge_descriptor,int> var_index;
		std::map<int,TGraph::edge_descriptor> inv_index;
		std::vector<int> utility;

		IloObjective obj = IloMaximize(env);
		IloRangeArray c(env);

		//find number of variables
		tie(ei,eend) = boost::edges(g);
		int ecount = 0;
		int varc = 0;
		for(;ei != eend; ++ei)
		{
			ecount++;
			if(g[*ei].fixed != 1)
			{
				if(g[*ei].selected!=1)
				{
					var_index[*ei] = varc;
					inv_index[varc] = *ei;
					g[*ei].selected = 1;
					if(g[*ei].coupled == 1)
					{
						var_index[coupled_map[*ei]] = varc;
						inv_index[varc] = coupled_map[*ei];
						g[coupled_map[*ei]].selected = 1;
					}
					utility.push_back(g[*ei].utility);
					varc++;
				}
			}
		}
		printf("ecount = %d varc = %d\n",ecount,varc);
		graph_traits<TGraph>::vertex_iterator vi,vend;

		IloBoolVarArray x(env,varc);

		for(int counter=0; counter < varc; counter++)
		{
			obj.setLinearCoef(x[counter],utility[counter]);
			//printf("utility[%d] = %d\n",counter,utility[counter]);
		}

		int vcount = -1;
		for(tie(vi,vend) = vertices(g); vi != vend; ++vi)
		{

			graph_traits<TGraph>::in_edge_iterator e_i,e_end;
			tie(e_i,e_end) = in_edges(*vi,g);
			bool once = false;
			float fixed_sum = 0;
			for(;e_i!=e_end; ++e_i)
			{
				if(g[*e_i].fixed == 0 )
				{
					if(once == false )
					{
						once = true;
						vcount++;
						c.add(IloRange(env,0,1));
					}
					c[vcount].setLinearCoef(x[var_index[*e_i]],1.0);
				}
				else
				{
					if(g[*e_i].coupled==0)
						fixed_sum += g[*e_i].selected;
					else
						fixed_sum += 0.5*g[*e_i].selected;
				}
			}
			if(once)
				c[vcount].setBounds(0,1-fixed_sum);
			graph_traits<TGraph>::out_edge_iterator e_2,e_end2;
			tie(e_2,e_end2) = out_edges(*vi,g);
			once = false;
			fixed_sum = 0;
			for(;e_2!=e_end2; ++e_2)
			{
				if(g[*e_2].fixed == 0 )
				{
					if(once == false )
					{
						once = true;
						vcount++;
						c.add(IloRange(env,0,1));
					}
					c[vcount].setLinearCoef(x[var_index[*e_2]],1.0);
				}
				else
				{
					if(g[*e_2].coupled==0)
						fixed_sum += g[*e_2].selected;
					else
						fixed_sum += 0.5*g[*e_2].selected;
				}
			}
			if(once)
				c[vcount].setBounds(0,1-fixed_sum);

		}
		printf("vcount = %d\n",vcount);
		//scanf("%*d");
		IloModel model(env);
		model.add(obj);
		model.add(c);


		//std::cout<<model<<std::endl;
		IloCplex cplex(model);
		if(!cplex.solve())
		{
			std::cerr << " Could not solve.. error"<<std::endl;
		}
		IloNumArray vals(env);
		env.out() << "Solution status = " << cplex.getStatus() << std::endl;
		env.out() << "Solution value  = " << cplex.getObjValue() << std::endl;
		cplex.getValues(vals, x);
		//env.out() << "Values        = " << vals << std::endl;
		std::cout << "Values.getSize() = "<< vals.getSize() << std::endl;
		for(int counter=0; counter< vals.getSize(); counter++)
		{
			//assert(g[inv_index[counter]].selected == 1);
			g[inv_index[counter]].selected = vals[counter];
			if(g[inv_index[counter]].coupled==1)
				g[coupled_map[inv_index[counter]]].selected = vals[counter];
			if(vals[counter]==1 && utility[counter] <0)
			{
				printf("wrong @ %d\n",counter);
			}
		}
	}
	catch(IloException &e)
	{
		std::cerr << e << std::endl;
	}
	env.end();
	//scanf("%*d");
}

void CellTracker::prune(int t)
{


	using boost::graph_traits;
	graph_traits<TGraph>::edge_iterator e_i,e_end;
	tie(e_i,e_end) = edges(g);
	for(;e_i!=e_end; ++e_i)
	{
		g[*e_i].selected = 0;
	}

}
void CellTracker::print_stats()
{
	using boost::graph_traits;

	graph_traits<TGraph>::edge_iterator e_i,e_end;

	float sne,sde,sae,sme;
	sne=sde=sae=sme=0;

	FILE*fp = fopen("C:\\Users\\arun\\Research\\Tracking\\harvard\\debug.txt","w");
	if(fp==NULL)
	{
		printf("Could not open debug.txt \n");
		_exit(-1);
	}
	tie(e_i,e_end) = edges(g);
	for(;e_i!=e_end; ++e_i)
	{
		if(g[*e_i].selected == 1)
		{
			if(g[*e_i].coupled == 1)
				sme += 0.5;
			else if(g[source(*e_i,g)].special == 1)
				sae++;
			else if(g[target(*e_i,g)].special == 1)
				sde++;
			else
				sne++;
		}
	}
	boost::graph_traits<TGraph>::vertex_iterator v_i,v_end;

	char mod;
	for(tie(v_i,v_end)=vertices(g); v_i != v_end; ++v_i)
	{
		boost::graph_traits<TGraph>::in_edge_iterator e_in,e_in_end;
		tie(e_in,e_in_end) = in_edges(*v_i,g);
		bool once = false;
		for(;e_in!=e_in_end;++e_in)
		{
			once = true;
			TrackVertex v = g[source(*e_in,g)];
			if(g[*e_in].selected == 1)
				mod = '*';
			else
				mod = ' ';
			if(v.special == 0)
				fprintf(fp,"[%d].%d%c ",v.t,fvector[v.t][v.findex].num,mod);
			else
				fprintf(fp,"A ");
		}

		if(g[*v_i].special == 0)
			fprintf(fp," -> [%d].%d ->",g[*v_i].t,fvector[g[*v_i].t][g[*v_i].findex].num);
		else if (once == false)
			fprintf(fp,"A ->");
		else
			fprintf(fp,"-> D");
		boost::graph_traits<TGraph>::out_edge_iterator e_out,e_out_end;
		tie(e_out,e_out_end) = out_edges(*v_i,g);
		for(;e_out!=e_out_end;++e_out)
		{
			TrackVertex v = g[target(*e_out,g)];
			if(g[*e_out].selected == 1)
				mod = '*';
			else
				mod = ' ';
			if(v.special == 0)
				fprintf(fp,"[%d].%d%c ",v.t,fvector[v.t][v.findex].num,mod);
			else
				fprintf(fp,"D ");
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
	printf("Normal = %d, Appear = %d, Disappear = %d, Merge/split = %d\n",int(sne+0.5),int(sae+0.5),int(sde+0.5),int(sme+0.5));



	//PAUSE;
}


int CellTracker::is_overlapping(TGraph::edge_descriptor e)
{
	//printf("Entered is_overlapping\n");
	TGraph::vertex_descriptor v1, v2;
	v1 = source(e,g);
	v2 = target(e,g);
	if(g[v1].special  == 1 || g[v2].special == 1)
		return 2;
	//printf("Not returned yet\n");
	FeaturesType f1,f2;
	//printf("C-\n");
	f1 = fvector[g[v1].t][g[v1].findex];
	//printf("1\n");
	f2 = fvector[g[v2].t][g[v2].findex];
	//printf("2\n");
	if(overlap(f1.BoundingBox, f2.BoundingBox) > 0)
	{
		return 1;
	}
	return 0;
}
int CellTracker::is_alpha_overlapping(TGraph::edge_descriptor e,float alpha) // alpha \in [0-1]
{
	//printf("Entered is_overlapping\n");
	TGraph::vertex_descriptor v1, v2;
	v1 = source(e,g);
	v2 = target(e,g);
	if(g[v1].special  == 1 || g[v2].special == 1)
		return 2;
	//printf("Not returned yet\n");
	FeaturesType f1,f2;
	//printf("C-\n");
	f1 = fvector[g[v1].t][g[v1].findex];
	//printf("1\n");
	f2 = fvector[g[v2].t][g[v2].findex];
	//printf("2\n");
	if(overlap(f1.BoundingBox, f2.BoundingBox) > alpha*MIN(f1.ScalarFeatures[FeaturesType::BBOX_VOLUME],f2.ScalarFeatures[FeaturesType::BBOX_VOLUME]))
	{
		return 1;
	}
	return 0;
}

void CellTracker::enforce_overlap()
{
	printf("Entered enforce overlap\n");
	// if A->B has overlap of bounding boxes, then all out_edges of A should have overlap and all in_edges of B should have overlap.

	TGraph::vertex_iterator vi,vend;
	float alpha = 0.30;
	for(tie(vi,vend) = vertices(g); vi!=vend; ++vi)
	{
		// in edges
		TGraph::in_edge_iterator e_in, e_in_end,e_in2,e_in_end2;
		tie(e_in,e_in_end) = in_edges(*vi,g);
		bool isovlap = false;
		//printf("Ein1\n");
		for(;e_in!=e_in_end; ++e_in)
		{
			if(is_alpha_overlapping(*e_in,alpha)==1)
			{
				isovlap = true;
				break;
			}
		}
		if(isovlap)
		{
			//printf("E_in2\n");
			tie(e_in2,e_in_end2) = in_edges(*vi,g);
			for(;e_in2!=e_in_end2; ++e_in2)
			{
				if(is_alpha_overlapping(*e_in2,alpha)==0)
				{
					
					if(g[*e_in2].coupled==1)
					{
						if(is_alpha_overlapping(coupled_map[*e_in2],alpha)==0)
						{
							g[coupled_map[*e_in2]].utility = -(UTILITY_MAX-1);
							g[*e_in2].utility = -(UTILITY_MAX-1);
						}
					}
					else
					{
						g[*e_in2].utility = -(UTILITY_MAX-1);
					}
				}
			}
			//printf("Finished E_in2\n");
			TGraph::vertex_descriptor vsource = source(*e_in,g);
			TGraph::out_edge_iterator e_out, e_out_end;
			tie(e_out,e_out_end) = out_edges(vsource,g);
			//printf("E_out\n");
			for(;e_out!=e_out_end; ++e_out)
			{
				if(is_alpha_overlapping(*e_out,alpha)==0)
				{
					
					if(g[*e_out].coupled==1)
					{
						if(is_alpha_overlapping(coupled_map[*e_out],alpha)==0)
						{
							g[coupled_map[*e_out]].utility = -(UTILITY_MAX -1);
							g[*e_out].utility = -(UTILITY_MAX-1);
						}
					}
					else
					{
						g[*e_out].utility = -(UTILITY_MAX-1);
					}
				}
			}
		}
	}

}


bool CellTracker::is_merge_node(CellTracker::TGraph::vertex_descriptor vd)
{
	bool yes = false;
	TGraph::in_edge_iterator in_e, in_end;
	tie(in_e,in_end) = in_edges(vd,g);
	int in_count = 0;
	for(;in_e!=in_end; ++in_e,++in_count);
	if(in_count == 2)
	{
		yes = true;
	}
	return yes;
}

bool CellTracker::is_split_node(CellTracker::TGraph::vertex_descriptor vd)
{

	bool yes = false;
	TGraph::out_edge_iterator out_e, out_end;
	tie(out_e,out_end) = out_edges(vd,g);
	int out_count = 0;
	for(;out_e!=out_end; ++out_e,++out_count);
	if(out_count == 2)
	{
		yes = true;
	}
	return yes;
}

bool CellTracker::is_simple_node(CellTracker::TGraph::vertex_descriptor vd)
{
	if(in_degree(vd,g) == 1 && out_degree(vd,g) == 1)
		return true;
	return false;
}

CellTracker::TGraph::vertex_descriptor CellTracker::get_parent(CellTracker::TGraph::vertex_descriptor vd)
{
	if(in_degree(vd,g)>1)
	{
		printf("Something wrong: Request for parent when there are more than one...\n");
		PAUSE;
	}
	TGraph::in_edge_iterator in_e, in_end;
	tie(in_e,in_end) = in_edges(vd,g);
	if(in_e == in_end)
		return TGraph::null_vertex();
	return source(*in_e,g);
}

CellTracker::TGraph::vertex_descriptor CellTracker::get_child(CellTracker::TGraph::vertex_descriptor vd)
{
	if(out_degree(vd,g)>1)
	{
		printf("Something wrong: Request for child when there are more than one...\n");
		PAUSE;
	}
	TGraph::out_edge_iterator out_e, out_end;
	tie(out_e,out_end) = out_edges(vd,g);
	if(out_e == out_end)
		return TGraph::null_vertex();
	return target(*out_e,g);
}
bool CellTracker::is_separate(CellTracker::TGraph::vertex_descriptor v1, CellTracker::TGraph::vertex_descriptor v2)
{
	// assuming not special vertices or sth like that

	int t1 = g[v1].t;
	int t2 = g[v2].t;
	int i1 = g[v1].findex;
	int i2 = g[v2].findex;
	LabelImageType::Pointer p1,p2;
	p1 = limages[t1][i1];
	p2 = limages[t2][i2];



	LabelImageType::SizeType ls;

	int lbounds[6];

	lbounds[0] = MIN(fvector[t1][i1].BoundingBox[0],fvector[t2][i2].BoundingBox[0]);
	lbounds[2] = MIN(fvector[t1][i1].BoundingBox[2],fvector[t2][i2].BoundingBox[2]);
	lbounds[4] = MIN(fvector[t1][i1].BoundingBox[4],fvector[t2][i2].BoundingBox[4]);
	lbounds[1] = MAX(fvector[t1][i1].BoundingBox[1],fvector[t2][i2].BoundingBox[1]);
	lbounds[3] = MAX(fvector[t1][i1].BoundingBox[3],fvector[t2][i2].BoundingBox[3]);
	lbounds[5] = MAX(fvector[t1][i1].BoundingBox[5],fvector[t2][i2].BoundingBox[5]);

	ls[0] = lbounds[1]-lbounds[0]+1;
	ls[1] = lbounds[3]-lbounds[2]+1;
	ls[2] = lbounds[5]-lbounds[4]+1;

	LabelImageType::Pointer p = LabelImageType::New();
	LabelImageType::IndexType lindex;
	lindex.Fill(0);
	LabelImageType::RegionType lregion;
	lregion.SetIndex(lindex);
	lregion.SetSize(ls);
	p->SetRegions(lregion);
	p->Allocate();
	p->FillBuffer(0);
	LabelIteratorType liter1(p1,p1->GetLargestPossibleRegion());


	lindex[0] = fvector[t1][i1].BoundingBox[0]-lbounds[0];
	lindex[1] = fvector[t1][i1].BoundingBox[2]-lbounds[2];
	lindex[2] = fvector[t1][i1].BoundingBox[4]-lbounds[4];

	lregion.SetSize(p1->GetLargestPossibleRegion().GetSize());
	lregion.SetIndex(lindex);

	LabelIteratorType liter(p,lregion);
	for(liter1.GoToBegin(),liter.GoToBegin();!liter1.IsAtEnd(); ++liter1,++liter)
	{
		if(liter1.Get()==fvector[t1][i1].num)
			liter.Set(fvector[t1][i1].num);
	}

	LabelIteratorType liter2(p2,p2->GetLargestPossibleRegion());

	lindex[0] = fvector[t2][i2].BoundingBox[0]-lbounds[0];
	lindex[1] = fvector[t2][i2].BoundingBox[2]-lbounds[2];
	lindex[2] = fvector[t2][i2].BoundingBox[4]-lbounds[4];
	lregion.SetIndex(lindex);
	lregion.SetSize(p2->GetLargestPossibleRegion().GetSize());

	liter = LabelIteratorType(p,lregion);

	for(liter2.GoToBegin(),liter.GoToBegin();!liter2.IsAtEnd(); ++liter2,++liter)
	{
		if(liter2.Get()==fvector[t2][i2].num)
			liter.Set(fvector[t2][i2].num);
	}

	typedef itk::ConstNeighborhoodIterator< LabelImageType > ConstNeighType;
	typedef itk::NeighborhoodIterator < LabelImageType > NeighborhoodType;
	typedef itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator < LabelImageType > FaceCalculatorType;

	ConstNeighType::RadiusType rad;
	float radius = 1;
	rad.Fill(radius);

	FaceCalculatorType facecalc;
	FaceCalculatorType::FaceListType facelist;


	facelist = facecalc(p,p->GetLargestPossibleRegion(),rad);

	ConstNeighType nit1,nit2;
	FaceCalculatorType::FaceListType::iterator fit;

	int found = 0;
	int val1 = fvector[t1][i1].num;
	int val2 = fvector[t2][i2].num;

	for(fit = facelist.begin(); fit!=facelist.end(); ++fit)
	{
		nit1 = ConstNeighType(rad,p,*fit);

		for(nit1.GoToBegin(); !nit1.IsAtEnd(); ++nit1)
		{
			if(nit1.GetCenterPixel()!=0)
			{
				unsigned short value = nit1.GetCenterPixel();
				unsigned short other = val1+val2-value;
				found = (nit1.GetPrevious(0)==other) + (nit1.GetNext(0)==other) + (nit1.GetPrevious(1)==other) + (nit1.GetNext(1)==other) + (nit1.GetPrevious(2)==other) + (nit1.GetNext(2)==other);
				if(found>0)
					break;
			}
		}
		if(found>0)
			break;
	}

	if(found>0)
	{
		//yay!
		return false;
	}
	return true;
	/*FeaturesType f1,f2;
	f1 = fvector[g[v1].t][g[v1].findex];
	f2 = fvector[g[v2].t][g[v2].findex];

	float min[3];
	float max[3];
	min[0] = MAX(f1.BoundingBox[0],f2.BoundingBox[0]);
	min[1] = MAX(f1.BoundingBox[2],f2.BoundingBox[2]);
	min[2] = MAX(f1.BoundingBox[4],f2.BoundingBox[4]);

	max[0] = MIN(f1.BoundingBox[1],f2.BoundingBox[1]);
	max[1] = MIN(f1.BoundingBox[3],f2.BoundingBox[3]);
	max[2] = MIN(f1.BoundingBox[5],f2.BoundingBox[5]);

	LabelImageType::Pointer lim1,lim2;
	lim1 = limages[g[v1].t][g[v1].findex];
	lim2 = limages[g[v2].t][g[v2].findex];


	LabelImageType::IndexType ind1,ind2;
	LabelImageType::SizeType size1,size2;
	LabelImageType::RegionType region1,region2;

	ind1[0] = min[0] - f1.BoundingBox[0];
	ind1[1] = min[1] - f1.BoundingBox[2];
	ind1[2] = min[2] - f1.BoundingBox[4];

	ind2[0] = min[0] - f2.BoundingBox[0];
	ind2[1] = min[1] - f2.BoundingBox[2];
	ind2[2] = min[2] - f2.BoundingBox[4];

	size1[0] = max[0] - min[0] + 1;
	size1[1] = max[1] - min[1] + 1;
	size1[2] = max[2] - min[2] + 1;

	size2 = size1;

	region1.SetIndex(ind1);region1.SetSize(size1);
	region2.SetIndex(ind2);region2.SetSize(size2);


	typedef itk::ConstNeighborhoodIterator< LabelImageType > ConstNeighType;
	typedef itk::NeighborhoodIterator < LabelImageType > NeighborhoodType;
	typedef itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator < LabelImageType > FaceCalculatorType;

	ConstNeighType::RadiusType rad;
	float radius = 1;
	rad.Fill(radius);

	FaceCalculatorType facecalc1;
	FaceCalculatorType::FaceListType facelist1;
	FaceCalculatorType facecalc2;
	FaceCalculatorType::FaceListType facelist2;

	LabelImageType::Pointer lom1,lom2;
	typedef itk::RegionOfInterestImageFilter<LabelImageType,LabelImageType> ExtractFilterType;

	region1.Print(std::cout);
	region2.Print(std::cout);
	ExtractFilterType::Pointer ext1 = ExtractFilterType::New();
	ext1->SetInput(lim1);
	ext1->SetRegionOfInterest(region1);
	ext1->Update();
	lom1 = ext1->GetOutput();

	ExtractFilterType::Pointer ext2 = ExtractFilterType::New();
	ext2->SetInput(lim2);
	ext2->SetRegionOfInterest(region2);
	ext2->Update();
	lom2 = ext2->GetOutput();

	LabelImageType::RegionType region;
	LabelImageType::IndexType ind;
	LabelImageType::SizeType size;
	ind.Fill(0);
	size = region1.GetSize();
	region.SetIndex(ind);
	region.SetSize(size);


	facelist1 = facecalc1(lom1,region,rad);
	facelist2 = facecalc2(lom2,region,rad);

	ConstNeighType nit1,nit2;
	FaceCalculatorType::FaceListType::iterator fit1,fit2;

	int found = 0;
	for(fit1 = facelist1.begin(),fit2 = facelist2.begin(); fit1!=facelist1.end(); ++fit2,++fit1)
	{
	nit1 = ConstNeighType(rad,lim1,*fit1);
	nit2 = ConstNeighType(rad,lim2,*fit2);

	for(nit1.GoToBegin(),nit2.GoToBegin(); !nit1.IsAtEnd(); ++nit1,++nit2)
	{
	if(nit1.GetCenterPixel()!=0)
	{
	found = (nit2.GetPrevious(0)!=0) + (nit2.GetNext(0)!=0) + (nit2.GetPrevious(1)!=0) + (nit2.GetNext(1)!=0) + (nit2.GetPrevious(2)!=0) + (nit2.GetNext(2)!=0);
	if(found>0)
	break;
	}
	}
	if(found>0)
	break;
	}

	if(found>0)
	{
	//yay!
	return false;
	}
	return true;*/
}

std::vector<std::vector<bool> >  generate_all_binary_strings(int n)
{
	std::vector<std::vector<bool> > out;

	int max = 1;
	int ncopy = n;
	while(ncopy--)
	{
		max <<= 1;
	}
	for(int counter = 0; counter < max; counter++)
	{
		std::vector<bool> row;
		int co = counter;
		ncopy = n;
		while(ncopy--)
		{
			row.push_back(co&1);
			co >>=1;
		}
		out.push_back(row);
	}

	/*for(int co1 = 0; co1< out.size(); co1++)
	{
	for(int co2 = 0; co2 < out[co1].size(); co2++)
	{
	printf("%d ",int(out[co1][co2]));
	}
	printf("\n");
	}*/
	return out;
}

float CellTracker::get_LRUtility(std::vector< CellTracker::TGraph::vertex_descriptor > desc)
{
	float util = 0;
	//printf("called get_LRUtility \n");

	for(int counter =0; counter < desc.size()-2; counter++)
	{
		if(g[desc[counter]].special==1 || g[desc[counter+2]].special == 1)
		{
			;
		}
		else
		{
			//printf("[1]");
			//print_vertex(desc[counter],0);
			//printf("[2]");
			//print_vertex(desc[counter+1],0);
			//printf("[3]");
			//print_vertex(desc[counter+2],0);
			FeaturesType f1,f2,f3;
			f1 = fvector[g[desc[counter]].t][g[desc[counter]].findex];
			f2 = fvector[g[desc[counter+1]].t][g[desc[counter+1]].findex];
			f3 = fvector[g[desc[counter+2]].t][g[desc[counter+2]].findex];
			util = util + compute_LRUtility(f1,f2,f3);
		}
	}
	

	return  util;
}
void CellTracker::print_debug_info(void)
{
	print_vertex(rmap[2][0],1);
	//PAUSE;
}


void CellTracker::resolve_merges_and_splits()
{
	std::vector< std::vector < boost::graph_traits<TGraph>::vertex_descriptor > > to_resolve;

	TGraph::vertex_iterator v_i,v_end, v_next;
	while(1)
	{
		to_resolve.clear();
		std::map< TGraph::vertex_descriptor, char> tr_map;
		for(tie(v_i,v_end) = vertices(g);v_i!=v_end; v_i = v_next)
		{
			v_next = v_i;
			++v_next;
			TGraph::vertex_descriptor vd;
			vd = *v_i;
			bool issplit = is_split_node(vd);
			bool ismerge = is_merge_node(vd);
			bool incomplete = false;
			if(issplit || ismerge)
			{
				if(tr_map.find(*v_i)==tr_map.end())// if already added to to_resolve, dont redo it
				{
					std::vector< TGraph::vertex_descriptor > vvd;
					TGraph::vertex_descriptor vdt = vd;
					vvd.push_back(vdt);
					tr_map[vdt] = true;
					do
					{
						if(issplit && in_degree(vdt,g) == 1)
							vdt = get_parent(vdt);
						else if( ismerge && out_degree(vdt,g) == 1)
							vdt = get_child(vdt);
						else
							break;
						if(vdt == TGraph::null_vertex())
							break;
						if(ismerge && is_merge_node(vdt))
						{
							incomplete = true;
							break;
						}
						if(issplit && is_split_node(vdt))
						{
							incomplete = true;
							break;
						}
						vvd.push_back(vdt);
						tr_map[vdt] = true;
					}while(is_simple_node(vdt));
					if(incomplete)
						continue;
					if(issplit)
						std::reverse(vvd.begin(),vvd.end());
					to_resolve.push_back(vvd);
				}
			}
		}


		for(int counter =0; counter < to_resolve.size(); counter++)
		{
			printf("START --- \n");
			for(int counter1 = 0; counter1 < to_resolve[counter].size(); counter1++)
			{
				print_vertex(to_resolve[counter][counter1],0);
			}
			printf("\nEND --- \n");
		}


		//is_separate check

		printf("is_separate tests follows:\n");
	/*	printf("Is_separate? : %d\n", int(is_separate(rmap[3][11],rmap[3][15])));
		printf("Is_separate? : %d\n", int(is_separate(rmap[3][6],rmap[3][8])));
		printf("Is_separate? : %d\n", int(is_separate(rmap[3][3],rmap[3][4])));
		printf("Is_separate? : %d\n", int(is_separate(rmap[3][10],rmap[3][17])));*/
		//PAUSE;

		bool change_made = false;
		for(int counter = 0; counter < to_resolve.size(); counter++)
		{
			
			// steps:
			// go backward from the first one 
			// go forward from the last one
			// find evidence to split
			// otherwise merge the remaining
			bool got_evidence = false;
			bool got_evidence1 = false;
			bool got_evidence2 = false;
			TGraph::vertex_descriptor vd = to_resolve[counter][0];
			std::vector < std::pair < TGraph::vertex_descriptor,TGraph::vertex_descriptor  > > before;
			std::vector < std::pair < TGraph::vertex_descriptor, TGraph::vertex_descriptor> > after;
			if(is_merge_node(vd))
			{
				TGraph::edge_descriptor e1,e2;
				TGraph::in_edge_iterator in_e, in_end;
				tie(in_e,in_end) = in_edges(vd,g);
				e1 = *in_e;
				++in_e;
				e2 = *in_e;                                                                                                                                                                                                             
				TGraph::vertex_descriptor v1,v2;
				v1 = source(e1,g);
				v2 = source(e2,g);
				printf("Node vd \n");
				print_vertex(vd,0);
				do{
					printf("Checking v1,v2...\n");
					print_vertex(v1,0);
					print_vertex(v2,0);
					//printf("Done\n");
					if(is_separate(v1,v2))
					{
						printf("YaY! I got evidence now\n");
						got_evidence = true;
						got_evidence1 = true;
						//break;
					}
					std::pair < TGraph::vertex_descriptor,TGraph::vertex_descriptor> pair1;
					pair1.first = v1;
					pair1.second = v2;
					before.push_back(pair1);
					if( is_simple_node(v1) && is_simple_node(v2))
					{
						v1 = get_parent(v1);
						v2 = get_parent(v2);
						if(g[v1].special == 1 || g[v2].special == 1)
							break;
					}
					else
					{
						printf("Node v1 is simple? %d\n", is_simple_node(v1));
						printf("Node v2 is simple? %d\n", is_simple_node(v2));
						break;
					}
				}while(1);
			}
			std::reverse(before.begin(),before.end());
			//if(!got_evidence)
			{
				vd = to_resolve[counter][to_resolve[counter].size()-1];
				if(is_split_node(vd))
				{
					TGraph::edge_descriptor e1,e2;
					TGraph::out_edge_iterator out_e, out_end;
					tie(out_e,out_end) = out_edges(vd,g);
					e1 = *out_e;
					++out_e;
					e2 = *out_e;
					TGraph::vertex_descriptor v1,v2;
					v1 = target(e1,g);
					v2 = target(e2,g);
					printf("Node vd \n");
					print_vertex(vd,0);
					do{
						printf("Checking v1,v2...\n");
						print_vertex(v1,0);
						print_vertex(v2,0);
						//printf("Done\n");
						if(is_separate(v1,v2))
						{
							printf("YaY! I got evidence now\n");
							got_evidence = true;
							got_evidence2 = true;
							//break;
						}
						std::pair < TGraph::vertex_descriptor,TGraph::vertex_descriptor> pair1;
						pair1.first = v1;
						pair1.second = v2;
						after.push_back(pair1);
						if( is_simple_node(v1) && is_simple_node(v2))
						{
							v1 = get_child(v1);
							v2 = get_child(v2);
							if(g[v1].special == 1 || g[v2].special == 1)
								break;
						}
						else
						{
							printf("Node v1 is simple? %d\n", is_simple_node(v1));
							printf("Node v2 is simple? %d\n", is_simple_node(v2));
							break;
						}
					}while(1);
				}
			}

			if(got_evidence == true)
			{
				if(to_resolve[counter].size()>3)
				{
					std::vector<TGraph::edge_descriptor> marked_for_removal;
					TGraph::edge_descriptor dummye1;
					bool dummy_added;
					if(got_evidence1 == true)
					{
						int nbefore = 2;
						int nafter = 2;
						std::vector<TGraph::vertex_descriptor> option1,option2;
						for(int counter1 = before.size()- nbefore; counter1 < before.size(); counter1++)
						{
							option1.push_back(before[counter1].first);
							option2.push_back(before[counter1].second);
						}
						for(int counter1 = 0; counter1 < MIN(nafter,to_resolve[counter].size()); counter1++)
						{
							option1.push_back(to_resolve[counter][counter1]);
							option2.push_back(to_resolve[counter][counter1]);
						}
						float option1util = get_LRUtility(option1);
						float option2util = get_LRUtility(option2);
						clear_vertex(to_resolve[counter][0],g);
						if(option1util > option2util)
						{
							//tie(e1,added) = my_add_edge(before[before.size()-1].first,spvertices[0].second,1/*dummy utility*/,0,1,1,TRANSLATION);
							tie(dummye1, dummy_added) = my_add_edge(before[before.size()-1].first,to_resolve[counter][0],1,0,1,1,TRANSLATION);
							
						}
						else
						{
							tie(dummye1, dummy_added) = my_add_edge(before[before.size()-1].second,to_resolve[counter][0],1,0,1,1,TRANSLATION);
						}
						tie(dummye1, dummy_added) = my_add_edge(to_resolve[counter][0],to_resolve[counter][1],1,0,1,1,TRANSLATION);
					}
					if(got_evidence2 == true)
					{
						int nbefore = 2;
						int nafter = 2;
						std::vector<TGraph::vertex_descriptor> option1,option2;
						
						for(int counter1 = to_resolve[counter].size()-nafter; counter1 < to_resolve[counter].size(); counter1++)
						{
							option1.push_back(to_resolve[counter][counter1]);
							option2.push_back(to_resolve[counter][counter1]);
						}
						for(int counter1 = 0; counter1 < MIN(nbefore,after.size()); counter1++)
						{
							option1.push_back(after[counter1].first);
							option2.push_back(after[counter1].second);
						}
						float option1util = get_LRUtility(option1);
						float option2util = get_LRUtility(option2);
						clear_vertex(to_resolve[counter][to_resolve[counter].size()-1],g);
						if(option1util > option2util)
						{
							//tie(e1,added) = my_add_edge(before[before.size()-1].first,spvertices[0].second,1/*dummy utility*/,0,1,1,TRANSLATION);
							tie(dummye1, dummy_added) = my_add_edge(to_resolve[counter][to_resolve[counter].size()-1],after[0].first,1,0,1,1,TRANSLATION);
							
						}
						else
						{
							tie(dummye1, dummy_added) = my_add_edge(to_resolve[counter][to_resolve[counter].size()-1],after[0].second,1,0,1,1,TRANSLATION);
						}
						tie(dummye1, dummy_added) = my_add_edge(to_resolve[counter][to_resolve[counter].size()-2],to_resolve[counter][to_resolve[counter].size()-1],1,0,1,1,TRANSLATION);
					}
					continue;
				}
				printf("started got_evidence\n");
				
				change_made = true;	
				std::vector < TGraph::vertex_descriptor > sp = to_resolve[counter];
				std::vector < std::pair < TGraph::vertex_descriptor, TGraph::vertex_descriptor > > spvertices;
				
				for(int cosp =0; cosp < sp.size(); cosp++)
				{
					TGraph::vertex_descriptor vd = sp[cosp];

					std::vector<LabelImageType::Pointer> lout;
					std::vector<InputImageType::Pointer> rout;
					std::vector<FeaturesType> fvecout;
					printf("I'm trying to split:\n");
					printf("g[vd].t = %d g[vd].findex = %d g[vd].special = %d\n",g[vd].t, g[vd].findex,g[vd].special);
				
					printFeatures(fvector[g[vd].t][g[vd].findex]);
					_TRACE;
					SplitCell(limages[g[vd].t][g[vd].findex],rimages[g[vd].t][g[vd].findex],fvector[g[vd].t][g[vd].findex],fvar,lout,rout,fvecout);
					_TRACE;

					int ttemp = g[vd].t;
					fvecout[0].time = ttemp;
					fvecout[1].time = ttemp;

					fvecout[0].num = fvector[ttemp].size()-2;
					fvecout[1].num = fvector[ttemp].size()-1;

					limages[g[vd].t].push_back(lout[0]);
					limages[g[vd].t].push_back(lout[1]);

					rimages[g[vd].t].push_back(rout[0]);
					rimages[g[vd].t].push_back(rout[1]);

					fvector[g[vd].t].push_back(fvecout[0]);
					fvector[g[vd].t].push_back(fvecout[1]);

					TGraph::vertex_descriptor v1, v2 ;
					v1 = add_vertex(g);
					v2 = add_vertex(g);

					g[v1].special = 0;
					g[v1].t = g[vd].t;
					g[v1].findex = fvector[g[vd].t].size()-2;
					
					g[v2].special = 0;
					g[v2].t = g[vd].t;
					g[v2].findex = fvector[g[vd].t].size()-1;

					std::pair < TGraph::vertex_descriptor, TGraph::vertex_descriptor> pair1;
					pair1.first = v1;
					pair1.second = v2;
					spvertices.push_back(pair1);
				}

				// I have before, spvertices and after;
				int permsize = spvertices.size()-1;
				if(after.size()!=0)
					permsize++;
				if(before.size()!=0)
					permsize++;
				std::vector<std::vector<bool>> permutations = generate_all_binary_strings(permsize);
				float utilmax = -1;
				int utilpos = -1;
				_TRACE;
				int permpc = 0;
				for(int cop = 0; cop < permutations.size(); cop++)
				{
					permpc = 0;
					std::vector<TGraph::vertex_descriptor> desc;
					if(before.size()!=0)
					{
						for(int co1 = 0; co1 < before.size(); co1++)
							desc.push_back(before[co1].first);
						if(permutations[cop][permpc++] == 0)
							desc.push_back(spvertices[0].first);
						else
							desc.push_back(spvertices[0].second);
					}
					else
					{
						desc.push_back(spvertices[0].first);
					}
					for(int co1 = 1; co1 < spvertices.size(); co1++)
					{
						if(permutations[cop][permpc++] == 0)
							desc.push_back(spvertices[co1].first);
						else
							desc.push_back(spvertices[co1].second);
					}
					if(after.size()!=0)
					{
						if(permutations[cop][permpc++] == 0)
						{
							desc.push_back(after[0].first);
							for(int co1 = 1; co1 < after.size(); co1++)
								desc.push_back(after[co1].first);
						}
						else
						{
							desc.push_back(after[0].second);
							for(int co1 = 1; co1 < after.size(); co1++)
								desc.push_back(after[co1].second);
						}
					}

					float util1 = get_LRUtility(desc);

					desc.clear();
					permpc = 0;
					if(before.size()!=0)
					{
						for(int co1 = 0; co1 < before.size(); co1++)
							desc.push_back(before[co1].second);
						if(permutations[cop][permpc++] == 0)
							desc.push_back(spvertices[0].second);
						else
							desc.push_back(spvertices[0].first);
					}
					else
					{
						desc.push_back(spvertices[0].second);
					}
					for(int co1 = 1; co1 < spvertices.size(); co1++)
					{
						if(permutations[cop][permpc++] == 0)
							desc.push_back(spvertices[co1].second);
						else
							desc.push_back(spvertices[co1].first);
					}
					if(after.size()!=0)
					{
						if(permutations[cop][permpc++] == 0)
						{
							desc.push_back(after[0].second);
							for(int co1 = 1; co1 < after.size(); co1++)
								desc.push_back(after[co1].second);
						}
						else
						{
							desc.push_back(after[0].first);
							for(int co1 = 1; co1 < after.size(); co1++)
								desc.push_back(after[co1].first);
						}
					}
					float util2 = get_LRUtility(desc);
					for(int cod1 = 0; cod1 < desc.size(); cod1++)
					{
						if(g[desc[cod1]].t == 8 && fvector[g[desc[cod1]].t][g[desc[cod1]].findex].num == 17)
						{
							printf("utils: %0.2f %0.2f\n", util1, util2);
							for(int cod = 0; cod < desc.size(); cod++)
							{
								print_vertex(desc[cod],0);
							}
							//PAUSE;
						}
					}
					float util = util1+util2;
					//	printf("utils: %0.2f %0.2f\n", util1, util2);
					//printf("util = %f\n", util);
					if(util > utilmax)
					{
						utilmax = util;
						utilpos = cop;
					}
				}
				_TRACE;
				//yay! I have utilpos. All I have to do is to connect them with the right edges. first clear all original vertices. Dont delete the vertices now.
				for(int co1 = 0; co1 < to_resolve[counter].size(); co1++)
				{
					clear_vertex(to_resolve[counter][co1],g);
				}
				TGraph::edge_descriptor e1,e2;
				bool added = false;
				permpc = 0;
				bool flipped = false;
				TGraph::vertex_descriptor prev1,prev2;
				printf("before.size() = %d after.size() = %d utilpos = %d\n",before.size(),after.size(),utilpos);
				if(before.size()!=0)
				{
					if(permutations[utilpos][permpc++] == 0)
					{
						tie(e1,added) = my_add_edge(before[before.size()-1].first,spvertices[0].first,1/*dummy utility*/,0,1,1,TRANSLATION);
						tie(e2,added) = my_add_edge(before[before.size()-1].second,spvertices[0].second,1/*dummy utility*/,0,1,1,TRANSLATION);
						prev1 = spvertices[0].first;
						prev2 = spvertices[0].second;
					}
					else
					{
						tie(e1,added) = my_add_edge(before[before.size()-1].first,spvertices[0].second,1/*dummy utility*/,0,1,1,TRANSLATION);
						tie(e2,added) = my_add_edge(before[before.size()-1].second,spvertices[0].first,1/*dummy utility*/,0,1,1,TRANSLATION);
						prev1 = spvertices[0].second;
						prev2 = spvertices[0].first;
					}
				}
				else
				{
					prev1 = spvertices[0].first;
					prev2 = spvertices[0].second;
				}
				_TRACE;

				for(int co1 = 0; co1 < spvertices.size()-1; co1++)
				{
					if(permutations[utilpos][permpc++] == 0)
					{
						//tie(e1,added) = my_add_edge(spvertices[co1].first,spvertices[co1+1].first,1/*dummy utility*/,0,1,1,TRANSLATION);
						//tie(e2,added) = my_add_edge(spvertices[co1].second,spvertices[co1+1].second,1/*dummy utility*/,0,1,1,TRANSLATION);	
						tie(e1,added) = my_add_edge(prev1,spvertices[co1+1].first,1/*dummy utility*/,0,1,1,TRANSLATION);
						tie(e2,added) = my_add_edge(prev2,spvertices[co1+1].second,1/*dummy utility*/,0,1,1,TRANSLATION);	
						prev1 = spvertices[co1+1].first;
						prev2 = spvertices[co1+1].second;
					}
					else
					{
						//tie(e1,added) = my_add_edge(spvertices[co1].first,spvertices[co1+1].second,1/*dummy utility*/,0,1,1,TRANSLATION);
						//tie(e2,added) = my_add_edge(spvertices[co1].second,spvertices[co1+1].first,1/*dummy utility*/,0,1,1,TRANSLATION);
						tie(e1,added) = my_add_edge(prev1,spvertices[co1+1].second,1/*dummy utility*/,0,1,1,TRANSLATION);
						tie(e2,added) = my_add_edge(prev2,spvertices[co1+1].first,1/*dummy utility*/,0,1,1,TRANSLATION);
						prev1 = spvertices[co1+1].second;
						prev2 = spvertices[co1+1].first;
					}
				}
				_TRACE;
			//	printf("prev1\n");print_vertex(prev1,1);
		//		printf("prev2\n");print_vertex(prev2,1);
			//	printf("after[0].first\n");print_vertex(after[0].first,1);
			//	printf("after[0].second\n");print_vertex(after[0].second,1);
		//		printf("utilpos = %d, permc = %d, permutation.size() = %d [0].size = %d\n",utilpos, permpc, permutations.size(), permutations[0].size());
				if(after.size()!=0)
				{//_TRACE;
					//for(int co_ = 0; co_ < permutations.size(); co_++)
				//	{
				//		for(int co1_ = 0; co1_ < permutations[co_].size(); co1_++)
				//		{
			//				printf("permutations[%d][%d] = %d\n",co_,co1_,int(permutations[co_][co1_]));
				//		}
			//		}
					//_TRACE;
				//	printf("permutations[0][0] = %d\n", int(permutations[0][0]));
					if(int(permutations[utilpos][permpc]) == 0)
					{
						_TRACE;
						tie(e1,added) = my_add_edge(prev1,after[0].first,1/*dummy utility*/,0,1,1,TRANSLATION);
						_TRACE;
						tie(e1,added) = my_add_edge(prev2,after[0].second,1/*dummy utility*/,0,1,1,TRANSLATION);
					}
					else
					{
						_TRACE;
						tie(e1,added) = my_add_edge(prev1,after[0].second,1/*dummy utility*/,0,1,1,TRANSLATION);
						_TRACE;
						tie(e1,added) = my_add_edge(prev2,after[0].first,1/*dummy utility*/,0,1,1,TRANSLATION);
					}
					permpc++;
				}
				_TRACE;
				if(permpc != permutations[0].size())
				{
					printf("Something is wrong with permpc\n");
					PAUSE;
				}
				_TRACE;
				///PAUSE;
			}
			else
			{
				;
			}
			//PAUSE;
		}


		//
		_TRACE;
		if(change_made == false)
			break;
	}
}

int CellTracker::my_connected_components(std::vector<int> &component)
{
	boost::property_map<TGraph, boost::vertex_index_t>::type index;
	index = get(boost::vertex_index,g);

	for(int counter = 0; counter < component.size(); counter++)
	{
		component[counter] = -1;
	}

	TGraph::vertex_iterator vi,vend;

	int curcomp = -1;
	for(tie(vi,vend)=vertices(g); vi!=vend; ++vi)
	{
		if(component[index[*vi]]==-1)
		{
			curcomp++;
			std::queue<TGraph::vertex_descriptor> q;
			q.push(*vi);
			while(!q.empty())
			{
				TGraph::vertex_descriptor top = q.front();
				q.pop();
				component[index[top]] = curcomp;
				TGraph::out_edge_iterator ei, eend;
				for(tie(ei,eend) = out_edges(top,g);ei!=eend;++ei)
				{
					if(component[index[target(*ei,g)]] == -1)
					{
						q.push(target(*ei,g));
					}
				}
				TGraph::in_edge_iterator e2,eend2;
				for(tie(e2,eend2) = in_edges(top,g);e2!=eend2; ++e2)
				{
					if(component[index[source(*e2,g)]] == -1)
					{
						q.push(source(*e2,g));
					}
				}
			}
		}
	}
}


void CellTracker::print_all_LRUtilities(TGraph::vertex_descriptor v)
{
	boost::property_map<TGraph, boost::vertex_index_t>::type index;
	index = get(boost::vertex_index,g);
	graph_traits<TGraph>::in_edge_iterator e_in,e_in_end;
	graph_traits<TGraph>::out_edge_iterator e_out,e_out_end;

	tie(e_in,e_in_end) = in_edges(v,g);
	tie(e_out,e_out_end) = out_edges(v,g);


	std::set<TGraph::edge_descriptor> in_unique_set;
	std::set<TGraph::edge_descriptor> out_unique_set;
	in_unique_set.clear();
	out_unique_set.clear();
	for(;e_in != e_in_end;++e_in)
	{
		if(g[*e_in].coupled == 0)
		{
			in_unique_set.insert(*e_in);
		}
		else
		{
			if(in_unique_set.count(coupled_map[*e_in])==0)
			{
				in_unique_set.insert(*e_in);
			}
		}
	}
	for(;e_out != e_out_end; ++e_out)
	{
		int type = get_edge_type(*e_out);
		if(g[*e_out].coupled == 0)
		{
			out_unique_set.insert(*e_out);
		}
		else
		{
			if(out_unique_set.count(coupled_map[*e_out])==0)
			{
				out_unique_set.insert(*e_out);
			}
		}
	}

	std::set<TGraph::edge_descriptor>::iterator i1,i2;
	for(i1 = in_unique_set.begin(); i1!= in_unique_set.end(); ++i1)
	{
		for(i2 = out_unique_set.begin(); i2!= out_unique_set.end(); ++i2)
		{
			TGraph::vertex_descriptor v1,v2,v3;
			v1 = source(*i1,g);
			v2 = target(*i1,g);
			v3 = target(*i2,g);
			char label[1024];
			FeaturesType f1;
			if(g[v1].special==1)
			{
				if(in_degree(v1,g)==0)
					sprintf(label,"A%d",index[v1]);
				else
					sprintf(label,"D%d",index[v1]);
			}
			else
			{
				f1 = fvector[g[v1].t][g[v1].findex];
				sprintf(label,"%d,%d",g[v1].t,f1.num);
			}
			printf("%s-",label);
			if(g[v2].special==1)
			{
				if(in_degree(v2,g)==0)
					sprintf(label,"A%d",index[v2]);
				else
					sprintf(label,"D%d",index[v2]);
			}
			else
			{
				f1 = fvector[g[v2].t][g[v2].findex];
				sprintf(label,"%d,%d",g[v2].t,f1.num);
			}
			printf("%s-",label);
			if(g[v3].special==1)
			{
				if(in_degree(v3,g)==0)
					sprintf(label,"A%d",index[v3]);
				else
					sprintf(label,"D%d",index[v3]);
			}
			else
			{
				f1 = fvector[g[v3].t][g[v3].findex];
				sprintf(label,"%d,%d",g[v3].t,f1.num);
			}
			printf("%s____[%d]__[%d]___",label,g[*i1].utility,g[*i2].utility);
			printf("%f\n",compute_LRUtility_product(*i1,*i2));
		}
	}	
}
void CellTracker::run()
{

	TGraph::vertex_descriptor vt1 = TGraph::null_vertex();
	std::vector< TGraph::vertex_descriptor > vv;
	for(int counter=0; counter< fvector.size(); counter++)
	{
		vv.clear();
		for(int counter1 = 0; counter1 < fvector[counter].size(); counter1++)
		{
			vv.push_back(vt1);
		}
		rmap.push_back(vv);
	}

	int avc = 0;
	int dvc = 0;
	int nec = 0;
	int msec = 0;

	TGraph::vertex_descriptor v;

	for(int counter=0; counter< fvector[0].size(); counter++)
	{
		v = add_vertex(g);
		g[v].special = 0;
		g[v].t = 0;
		g[v].findex = counter;
		rmap[0][counter] = v;
	}
	
	populate_merge_candidates(0);
//	avc += add_appear_vertices(-1);

	_TRACE;
	first_t = clock();
	for(int t = 1; t < fvector.size(); t++)
	{
		int tmin = MAX(0,t-K);
		for(int counter=0; counter< fvector[t].size(); counter++)
		{
			v = add_vertex(g);
			g[v].special = 0;
			g[v].t = t;
			g[v].findex = counter;
			rmap[t][counter] = v;
		}
		//printf("T = %d\n",t);
		TIC		populate_merge_candidates(t);//TOC("populate_merge_candidates");
		TIC;		nec += this->add_normal_edges(tmin,t);// TOC("add_normal_edges()");
		TIC;		msec+= this->add_merge_split_edges(t);// TOC("add_merge_split_edges()");
		TIC;		dvc += this->add_disappear_vertices(t);// TOC("add_disappear_vertices()");
		TIC;		avc += this->add_appear_vertices(t-1);// TOC("add_appear_vertices()");
		printf("total edges = %d+%d+%d+%d = %d\n",nec,dvc,avc,msec,nec+dvc+avc+msec);
		//PAUSE;
		TIC;		prune(t);//TOC("prune()");
	}
	//dvc += this->add_disappear_vertices(fvector.size());

	//enforce_overlap();
	


	VectorPixelType col1,col2,col3,col4;
	col1[0] = 255;col1[1] = 0;col1[2] = 0;
	col2[0] = 255;col2[1] = 255;col2[2] = 0;
	col3[0] = 255;col3[1] = 255;col3[2] = 255;
	col4[0] = 255;col4[1] = 0;col4[2] = 255;


	col1[0] = 255;col1[1] = 0;col1[2] = 0;
	col2[0] = 255;col2[1] = 0;col2[2] = 0;
	col3[0] = 0;col3[1] = 173;col3[2] = 255;
	col4[0] = 0;col4[1] = 173;col4[2] = 255;


	boost::graph_traits<TGraph>::edge_iterator e_i,e_next,e_end;

	/*int num_edges_removed = 0;
	for(tie(e_i,e_end) = edges(g); e_i!= e_end ; e_i = e_next)
	{
		e_next = e_i;
		++e_next;
		if(g[*e_i].utility <1 )
		{
			num_edges_removed++;
			remove_edge(*e_i,g);
		}
	}*/
	//printf("I removed %d edges\n", num_edges_removed);
	//scanf("%*d");

	for(tie(e_i,e_end) = edges(g); e_i!= e_end ; ++e_i)
	{
		g[*e_i].selected = 0;
		TGraph::vertex_descriptor vert1 = source(*e_i,g),vert2 = target(*e_i,g);
		if(g[vert1].t >= g[vert2].t)
		{
			printf("problem here:");
			print_vertex(vert1,0);
			print_vertex(vert2,0);
			//scanf("%*d");
		}
		if(g[*e_i].coupled == 0)
			draw_line_for_edge(2,*e_i,col1,col2,1);
		else
			draw_line_for_edge(2,*e_i,col3,col4,-1);

	}


/*	while(1)
	{
		int tin, indin;
		scanf("%d %d",&tin,&indin);
		if(tin==-1)
			break;
		printFeatures(fvector[tin][indin]);
	}
*/

	print_debug_info();
	solve_higher_order();
	//solve_lip();
	print_stats();
	writeXGMML("C:\\Users\\arun\\Research\\Tracking\\harvard\\ch4_new.xgmml");
	boost::property_map<TGraph, boost::vertex_index_t>::type index_v;
	index_v = get(boost::vertex_index,g);
	/*while(1)
	{
		int v1;
		cin>>v1;
		if(v1<0)
			break;
		TGraph::vertex_iterator vv1,vv2;
		for(tie(vv1,vv2) = vertices(g);vv1!=vv2; ++vv1)
		{
			if(index_v[*vv1]==v1)
			{
				print_all_LRUtilities(*vv1);
				break;
			}
		}
	}*/
	
	for(tie(e_i,e_end) = edges(g); e_i!= e_end ; e_i = e_next)
	{
		e_next = e_i;
		++e_next;
		if(g[*e_i].selected == 0)
		{
			remove_edge(*e_i,g);
		}
	}

//	compute_feature_variances();


	for(tie(e_i,e_end) = edges(g); e_i!= e_end ; ++e_i)
	{
		if(g[*e_i].coupled == 0)
			draw_line_for_edge(1,*e_i,col1,col2,1);
		else
			draw_line_for_edge(1,*e_i,col3,col4,-1);

	}
	for(tie(e_i,e_end) = edges(g); e_i!= e_end ; e_i = e_next)
	{
		e_next = e_i;
		++e_next;
		if(get_edge_type(*e_i) == APPEAR || get_edge_type(*e_i) == DISAPPEAR )
		{
			remove_edge(*e_i,g);
		}
	}
	//scanf("%*d");


	resolve_merges_and_splits();

	for(tie(e_i,e_end) = edges(g); e_i!= e_end ; ++e_i)
	{
		if(g[*e_i].coupled == 0)
			draw_line_for_edge(3,*e_i,col1,col2,1);
		else
			draw_line_for_edge(3,*e_i,col3,col4,-1);

	}

	
	/*TGraph::vertex_iterator del,del_end,del_next,del_tmp;
	tie(del,del_end) = vertices(g);
	for(; del!=del_end; del = del_next)
	{
		del_next = del;
		++del_next;
		if(in_degree(*del,g) == 0 && out_degree(*del,g)==0)
		{
			clear_vertex(*del,g);
			remove_vertex(*del,g);
		}
		tie(del_tmp,del_end) = vertices(g);
	}
	*/

	/*
	for(tie(e_i,e_end) = edges(g); e_i!=e_end; e_i = e_next)
	{
		e_next = e_i;
		++e_next;
		TGraph::edge_iterator e_temp;
		tie(e_temp,e_end) = edges(g);
		if(g[*e_i].selected == 1)
		{

			TGraph::vertex_descriptor v1,v2;
			v1 = source(*e_i,g);
			v2 = target(*e_i,g);
			if(g[v1].special == 0 && g[v1].t == 4)
			{
				float x = fvector[4][g[v1].findex].Centroid[0];
				float y = fvector[4][g[v1].findex].Centroid[1];
				if(abs(x-233) < 10 && abs(y-119)<10)
				{
					print_vertex(v1,1);
					//scanf("%*d");
				}
			}
			if(g[v1].findex == 0 && g[v1].t == 4)
			{
				printf("The edge is still there.. g[v1].special = %d \n",g[v1].special);
				//scanf("%*d");
			}
			TGraph::edge_descriptor e;
			bool added;
			tie(e,added) = add_edge(v2,v1,g);
			if(added)
			{
				g[e].coupled = g[*e_i].coupled;
				g[e].fixed = 0;
				g[e].selected = 0;
				g[e].utility = g[*e_i].utility;
			}
			if(g[*e_i].coupled==0)
			{
				if(g[v1].special == 0 && g[v2].special == 0)
				{

					FeaturesType f1 = fvector[g[v1].t][g[v1].findex];
					FeaturesType f2 = fvector[g[v2].t][g[v2].findex];
					printf("v1 = %0.0f v2 = %0.0f\n",f1.ScalarFeatures[FeaturesType::VOLUME],f2.ScalarFeatures[FeaturesType::VOLUME]);
				}
				//draw_line_for_edge(*e_i,col1,col2,-2);
			}
			else
			{
				//draw_line_for_edge(*e_i,col3,col4,2);
			}
		}
	}*/

	for(tie(e_i,e_end) = edges(g); e_i!= e_end ; e_i = e_next)
	{
		e_next = e_i;
		++e_next;
		if(g[*e_i].selected == 0)
		{
			remove_edge(*e_i,g);
		}
	}


		TGraph::vertex_iterator v_i, v_end;

		printf("Started deleting vertices\n");
	int total_vertex_removed_count = 0;
	TGraph::vertex_iterator v_next,v_temp;
	while(1)
	{
	int vertex_removed_count = 0;
	for(tie(v_i,v_end)=vertices(g);v_i!=v_end; v_i = v_next)
	{
		v_next = v_i;
		++v_next;
		
		if(g[*v_i].special == 1)//in_degree(*v_i,g) == 0 && out_degree(*v_i,g)==0)
		{
			clear_vertex(*v_i,g);
			remove_vertex(*v_i,g);
			vertex_removed_count ++;
			break;
		}
		tie(v_temp,v_end) = vertices(g);
	}
	
	//printf(" I removed %d vertices\n", vertex_removed_count);
	if(vertex_removed_count==0)
		break;

	total_vertex_removed_count++;
	}
	printf("finished deleting vertices\n");

	//PAUSE;

	/*for(tie(e_i,e_end) = edges(g); e_i!= e_end ; ++e_i)
	{
		if(g[*e_i].coupled == 0)
			draw_line_for_edge(1,*e_i,col1,col2,1);
		else
			draw_line_for_edge(1,*e_i,col3,col4,-1);

	}*/
	boost::property_map<TGraph, boost::vertex_index_t>::type index;
	index = get(boost::vertex_index,g);

	std::vector<int> component(num_vertices(g));
	int num = my_connected_components(component);

	//int num = connected_components(g,&component[0]);
	int tempnum = num;
	printf("tempnum = %d num_vertices(g) = %d\n",tempnum,num_vertices(g));
	//PAUSE;




	std::vector<int> vertex_count(num);
	std::vector<int> in_vertices_count(num);
	std::vector<int> out_vertices_count(num);
	for(int counter=0; counter < vertex_count.size(); counter++)
	{
		vertex_count[counter]=0;
		in_vertices_count[counter] = 0;
		out_vertices_count[counter] = 0;
	}

	for(tie(v_i,v_end) = vertices(g);v_i!=v_end; ++v_i)
	{
		printf("%d %d # ", index[*v_i], component[index[*v_i]]);
		vertex_count[component[index[*v_i]]]++;
		if(in_degree(*v_i,g)==0)
			in_vertices_count[component[index[*v_i]]]++;
		if(out_degree(*v_i,g)==0)
			out_vertices_count[component[index[*v_i]]]++;
	}
	//PAUSE;
	for(int counter =0; counter < num; counter++)
	{
		if(vertex_count[counter]>1)
		{
			if(!(in_vertices_count[counter]==1 && out_vertices_count[counter] == 1))
			{
				//if(in_vertices_count[counter]!=out_vertices_count[counter])
				//	printf("component+1 = %d in_count = %d out_count = %d\n",component[counter]+1,in_vertices_count[counter],out_vertices_count[counter]);
			}
		}
		else
		{
			printf("Has only one vertex in the component\n");
		}
	}
	//PAUSE;

	tie(v_i,v_end) = vertices(g);
	for(;v_i!=v_end; ++v_i)
	{
		if(vertex_count[component[index[*v_i]]]==1)
		{
			//if(g[*v_i].special ==0)
			{
				printf("here is a vertex : \n", *v_i);
				print_vertex(*v_i,1);
			}
			
		}
		else
		{
			//printf("%d ",component[index[*v_i]]);
		}
	}



	//PAUSE;
	tie(v_i,v_end) = vertices(g);

	int max1 = -1;
	for(;v_i!=v_end;++v_i)
	{
		if(g[*v_i].special == 0 )
		{
			g[*v_i].new_label = component[index[*v_i]]+1;
			old_to_new[g[*v_i].t][fvector[g[*v_i].t][g[*v_i].findex].num]=component[index[*v_i]]+1;
			max1 = MAX(max1,component[index[*v_i]]+1);
		}
	}



	printf("number of connected components = %d reduced = %d maxcomp = %d\n",num,tempnum,max1);
	//PAUSE;
}


void CellTracker::print_vertex(TGraph::vertex_descriptor v, int depth)
{

	if(g[v].special == 0)
	{
		FeaturesType f = fvector[g[v].t][g[v].findex];
		printf("fvector index = %d t = %d Centroid =(%0.0f,%0.0f,%0.0f) ",g[v].findex, g[v].t,f.Centroid[0],f.Centroid[1],f.Centroid[2]);
	}
	else
	{
		printf("special = 1");
	}
	printf(" in_degree = %d out_degree = %d\n",in_degree(v,g),out_degree(v,g));


	if(depth > 0)
	{
		TGraph::in_edge_iterator e_i,e_end;
		for(tie(e_i,e_end) = in_edges(v,g); e_i!=e_end; ++e_i)
		{
			printf("in w=%d\t",g[*e_i].utility);
			print_vertex(source(*e_i,g),depth-1);
		}
		TGraph::out_edge_iterator e_o,e_end2;
		for(tie(e_o,e_end2) = out_edges(v,g); e_o!=e_end2; ++e_o)
		{
			printf("out w=%d\t",g[*e_o].utility);
			print_vertex(target(*e_o,g),depth-1);
		}
	}


}
LabelImageType::Pointer CellTracker::getOutputAtTime(int t)
{

	LabelImageType::Pointer lim = LabelImageType::New();

	LabelImageType::SizeType lsize;
	LabelImageType::RegionType lregion;
	LabelImageType::IndexType lindex;

	lindex.Fill(0);
	lsize[0] = fvar.BoundingBox[1]-fvar.BoundingBox[0]+1;
	lsize[1] = fvar.BoundingBox[3]-fvar.BoundingBox[2]+1;
	lsize[2] = fvar.BoundingBox[5]-fvar.BoundingBox[4]+1;
	lregion.SetSize(lsize);
	lregion.SetIndex(lindex);

	lim->SetRegions(lregion);
	lim->Allocate();

	TGraph::vertex_iterator v_i,v_end;
	for(tie(v_i,v_end) = vertices(g); v_i!=v_end; ++v_i)
	{
		if(g[*v_i].special == 0 && g[*v_i].t == t)
		{
			int find = g[*v_i].findex;
			if(find > fvector[t].size()-1 || find> limages[t].size()-1)
			{
				printf("We have a problem here\n");
			}
			lsize[0] = fvector[t][find].BoundingBox[1]-fvector[t][find].BoundingBox[0]+1;
			lsize[1] = fvector[t][find].BoundingBox[3]-fvector[t][find].BoundingBox[2]+1;
			lsize[2] = fvector[t][find].BoundingBox[5]-fvector[t][find].BoundingBox[4]+1;
			lindex.Fill(0);
			lregion.SetSize(lsize);
			lregion.SetIndex(lindex);


			//limages[t][find]->GetLargestPossibleRegion().Print(std::cout);
			LabelIteratorType liter1(limages[t][find],lregion);

			lindex[0] = fvector[t][find].BoundingBox[0];
			lindex[1] = fvector[t][find].BoundingBox[2];
			lindex[2] = fvector[t][find].BoundingBox[4];
			lregion.SetIndex(lindex);

			LabelIteratorType liter2(lim,lregion);

			for(liter1.GoToBegin(),liter2.GoToBegin();!liter1.IsAtEnd();++liter1,++liter2)
			{
				if(liter1.Get()!=0)
					liter2.Set(g[*v_i].new_label);
			}

		}
	}

	return lim;
}


void CellTracker::writeGraphViz(char * filename)
{
	std::ofstream f;
	f.open(filename,std::ios::out);
	write_graphviz(f,g);
	f.close();
}

void CellTracker::writeGraphML(char * filename)
{
	std::ofstream f;
	f.open(filename,std::ios::out);
	using boost::dynamic_properties;

	dynamic_properties dp;
	dp.property("time",get(&TrackVertex::t,g));
	dp.property("special",get(&TrackVertex::special,g));
	dp.property("utility",get(&TrackEdge::utility,g));
	dp.property("coupled",get(&TrackEdge::coupled,g));
	write_graphml(f,g,dp,true);
	f.close();
}



void CellTracker::writeXGMML(char *filename)
{
	FILE *fp = fopen(filename,"w");

	fprintf(fp,"<?xml version=\"1.0\"?>\n");
	fprintf(fp,"<graph id=\"Antonio_data\" label=\"Antonio_data\" directed=\"1\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xmlns:ns1=\"http://www.w3.org/1999/xlink\" xmlns:dc=\"http://purl.org/dc/elements/1.1/\" xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" xmlns=\"http://www.cs.rpi.edu/XGMML\">\n");
	
	TGraph::vertex_iterator vi,vend;

	boost::property_map<TGraph, boost::vertex_index_t>::type index;
	index = get(boost::vertex_index,g);

	for(tie(vi,vend)=vertices(g);vi!=vend;++vi)
	{
		char label[1024];
		FeaturesType f1;
		if(g[*vi].special==1)
		{
			if(in_degree(*vi,g)==0)
				sprintf(label,"A%d",index[*vi]);
			else
				sprintf(label,"D%d",index[*vi]);
		}
		else
		{
			f1 = fvector[g[*vi].t][g[*vi].findex];
			sprintf(label,"%d,%d",g[*vi].t,f1.num);
		}
		
		fprintf(fp,"\t<node id=\"%d\" label=\"%s\"\>\n",index[*vi],label);
		if(g[*vi].special==0)
		{
			fprintf(fp,"\t\t<graphics type=\"ELLIPSE\" x=\"%d\" y=\"%d\" fill=\"#ff9999\"/>\n",f1.time*100,f1.num*100);
			fprintf(fp,"\t\t<att type=\"integer\" name=\"time\" value=\"%d\"/>\n",g[*vi].t);
			fprintf(fp,"\t\t<att type=\"integer\" name=\"num\" value=\"%d\"/>\n",f1.num);
			fprintf(fp,"\t\t<att type=\"real\" name=\"x\" value=\"%0.2f\"/>\n",f1.Centroid[0]);
			fprintf(fp,"\t\t<att type=\"real\" name=\"y\" value=\"%0.2f\"/>\n",f1.Centroid[1]);
			fprintf(fp,"\t\t<att type=\"real\" name=\"z\" value=\"%0.2f\"/>\n",f1.Centroid[2]);
			fprintf(fp,"\t\t<att type=\"integer\" name=\"volume\" value=\"%d\"/>\n",int(f1.ScalarFeatures[FeaturesType::VOLUME]));
			fprintf(fp,"\t\t<att type=\"integer\" name=\"index\" value=\"%d\"/>\n",index[*vi]);
		}
		fprintf(fp,"\t\t<att type=\"real\" name=\"special\" value=\"%d\"/>\n",g[*vi].special);
		printf("%d\n",index[*vi]);
		fprintf(fp,"\t</node>\n");
	}

	TGraph::edge_iterator ei,eend;
	int edge_index = 0;
	for(tie(ei,eend)=edges(g);ei!=eend; ++ei)
	{
		TGraph::vertex_descriptor v1 = source(*ei,g);
		TGraph::vertex_descriptor v2 = target(*ei,g);
		TGraph::out_edge_iterator oe,oend;
		bool once = false;
		std::string label ="";
		bool coupled = false;
		std::string descriptor ="";
		bool selected = false;
		for(tie(oe,oend)=out_edges(v1,g);oe!=oend;++oe)
		{
			char buff[1024];
			if(target(*oe,g)!=v2)
				continue;
			coupled = coupled | g[*oe].coupled;
			selected = selected | g[*oe].selected;
			char ch=' ';
			int edge_type = get_edge_type(*oe);
			switch(edge_type)
			{
				case APPEAR: ch='A';break;
				case DISAPPEAR: ch='D';break;
				case MERGE: ch='M';break;
				case SPLIT: ch='S';break;
				case TRANSLATION: ch='T';break;
			}
			//if(g[*oe].utility <2)
			//	continue;
			if(!once)
			{
				sprintf(buff,"%c%d",ch,g[*oe].utility);
				label = label + buff;
				sprintf(buff,"%d",*oe);
				descriptor = descriptor + buff;
				once = true;
			}
			else
			{
				sprintf(buff,",%c%d",ch,g[*oe].utility);
				label = label + buff;
				sprintf(buff,",%d",*oe);
				descriptor = descriptor + buff;
			}
		}
		if(label.size()!=0)
		{
			fprintf(fp,"\t<edge id=\"%d\" label=\"%s\" source=\"%d\" target=\"%d\">\n",++edge_index,label.c_str(),index[source(*ei,g)],index[target(*ei,g)]);
			fprintf(fp,"\t\t<att type=\"integer\" name=\"coupled\" value=\"%d\"/>\n",coupled);
			fprintf(fp,"\t\t<att type=\"integer\" name=\"utility\" value=\"%d\"/>\n",g[*ei].utility);
			fprintf(fp,"\t\t<att type=\"integer\" name=\"selected\" value=\"%d\"/>\n",selected);
			fprintf(fp,"\t\t<att type=\"integer\" name=\"type\" value=\"%d\"/>\n",g[*ei].type);
			fprintf(fp,"\t\t<att type=\"string\" name=\"label\" value=\"%s\"/>\n",label.c_str());
			//fprintf(fp,"\t\t<att type=\"string\" name=\"descriptor\" value=\"%s\"/>\n",descriptor.c_str());
			fprintf(fp,"\t</edge>\n");
		}
	}
	fprintf(fp,"</graph>\n");
	fclose(fp);
}

void CellTracker::writeXGMML_secondorder(char *filename, std::vector< CellTracker::LREdge > &lredges,std::vector<int> &utility, IloNumArray& vals )
{
	FILE *fp = fopen(filename,"w");

	fprintf(fp,"<?xml version=\"1.0\"?>\n");
	fprintf(fp,"<graph id=\"Antonio_data_secondorder\" label=\"Antonio_data_secondorder\" directed=\"1\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xmlns:ns1=\"http://www.w3.org/1999/xlink\" xmlns:dc=\"http://purl.org/dc/elements/1.1/\" xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" xmlns=\"http://www.cs.rpi.edu/XGMML\">\n");
	
	TGraph::vertex_iterator vi,vend;

	boost::property_map<TGraph, boost::vertex_index_t>::type index;
	index = get(boost::vertex_index,g);

	//for(tie(vi,vend)=vertices(g);vi!=vend;++vi)
	//{
	//	char label[1024];
	//	FeaturesType f1;
	//	if(g[*vi].special==1)
	//	{
	//		if(in_degree(*vi,g)==0)
	//			sprintf(label,"A%d",index[*vi]);
	//		else
	//			sprintf(label,"D%d",index[*vi]);
	//	}
	//	else
	//	{
	//		f1 = fvector[g[*vi].t][g[*vi].findex];
	//		sprintf(label,"%d,%d",g[*vi].t,f1.num);
	//	}
	//	
	//	fprintf(fp,"\t<node id=\"%d\" label=\"%s\"\>\n",index[*vi],label);
	//	if(g[*vi].special==0)
	//	{
	//		fprintf(fp,"\t\t<graphics type=\"ELLIPSE\" x=\"%d\" y=\"%d\" fill=\"#ff9999\"/>\n",f1.time*100,f1.num*100);
	//		fprintf(fp,"\t\t<att type=\"integer\" name=\"time\" value=\"%d\"/>\n",g[*vi].t);
	//		fprintf(fp,"\t\t<att type=\"integer\" name=\"num\" value=\"%d\"/>\n",f1.num);
	//		fprintf(fp,"\t\t<att type=\"real\" name=\"x\" value=\"%0.2f\"/>\n",f1.Centroid[0]);
	//		fprintf(fp,"\t\t<att type=\"real\" name=\"y\" value=\"%0.2f\"/>\n",f1.Centroid[1]);
	//		fprintf(fp,"\t\t<att type=\"real\" name=\"z\" value=\"%0.2f\"/>\n",f1.Centroid[2]);
	//		fprintf(fp,"\t\t<att type=\"integer\" name=\"volume\" value=\"%d\"/>\n",int(f1.ScalarFeatures[FeaturesType::VOLUME]));
	//		fprintf(fp,"\t\t<att type=\"integer\" name=\"index\" value=\"%d\"/>\n",index[*vi]);
	//	}
	//	fprintf(fp,"\t\t<att type=\"real\" name=\"special\" value=\"%d\"/>\n",g[*vi].special);
	//	printf("%d\n",index[*vi]);
	//	fprintf(fp,"\t</node>\n");
	//}

	TGraph::edge_iterator ei,eend;
	int edge_index = 1;
	for(tie(ei,eend)=edges(g);ei!=eend; ++ei)
	{
		TGraph::vertex_descriptor v1 = source(*ei,g);
		TGraph::vertex_descriptor v2 = target(*ei,g);
		bool once = false;
		std::string label ="";
		bool coupled = false;
		std::string descriptor ="";
		bool selected = false;
		int edge_type = get_edge_type(*ei);
		g[*ei].index = edge_index++;
		//if(g[*ei].utility <2)
		//	continue;
		char ch;
		switch(edge_type)
		{
				case APPEAR: ch='A';break;
				case DISAPPEAR: ch='D';break;
				case MERGE: ch='M';break;
				case SPLIT: ch='S';break;
				case TRANSLATION: ch='T';break;
		}
		label += ch;
		char buff[1024];
		switch(g[v1].special)
		{
		case 0:
			sprintf(buff,"%d,%d",g[v1].t,fvector[g[v1].t][g[v1].findex].num);
			label+=buff;
			break;
		case 1:
			label+="A";
		}
		label += "-";
		switch(g[v2].special)
		{
		case 0:
			sprintf(buff,"%d,%d",g[v2].t,fvector[g[v2].t][g[v2].findex].num);
			label+=buff;
			break;
		case 1:
			label+="D";
		}

		fprintf(fp,"\t<node id=\"%d\" label=\"%s\">\n",g[*ei].index,label.c_str());
		fprintf(fp,"\t\t<graphics type=\"RECTANGLE\" x=\"%d\" y=\"%d\" fill=\"#ff9999\"/>\n",rand()%1000,rand()%1000);
		fprintf(fp,"</node>\n");
	}


		//fprintf(fp,"\t\t<att type=\"string\" name=\"source\" value=\"%s\"/>\n",);
		//		fprintf(fp,"\t\t<att type=\"integer\" name=\"num\" value=\"%d\"/>\n",f1.num);
		//		fprintf(fp,"\t\t<att type=\"real\" name=\"x\" value=\"%0.2f\"/>\n",f1.Centroid[0]);
		//		fprintf(fp,"\t\t<att type=\"real\" name=\"y\" value=\"%0.2f\"/>\n",f1.Centroid[1]);
		//		fprintf(fp,"\t\t<att type=\"real\" name=\"z\" value=\"%0.2f\"/>\n",f1.Centroid[2]);
		//		fprintf(fp,"\t\t<att type=\"integer\" name=\"volume\" value=\"%d\"/>\n",int(f1.ScalarFeatures[FeaturesType::VOLUME]));
		//		fprintf(fp,"\t\t<att type=\"integer\" name=\"index\" value=\"%d\"/>\n",index[*vi]);
	int se_index = 0;
	for(tie(ei,eend)=edges(g);ei!=eend; ++ei)
	{
		for(int counter=0; counter < g[*ei].backlre.size(); counter++)
		{
			char buff[1024];
			sprintf(buff,"%d",utility[g[*ei].backlre[counter]]);
			fprintf(fp,"\t<edge id=\"%d\" label=\"%s\" source=\"%d\" target=\"%d\">\n",++se_index,buff,g[*ei].index,g[lredges[g[*ei].backlre[counter]].front].index);
			fprintf(fp,"\t\t<att type=\"string\" name=\"label\" value=\"%s\"/>\n",buff);
			/*fprintf(fp,"\t\t<att type=\"integer\" name=\"coupled\" value=\"%d\"/>\n",coupled);
			fprintf(fp,"\t\t<att type=\"integer\" name=\"utility\" value=\"%d\"/>\n",g[*ei].utility);*/
			fprintf(fp,"\t\t<att type=\"integer\" name=\"selected\" value=\"%d\"/>\n",int(vals[g[*ei].backlre[counter]]));
			/*fprintf(fp,"\t\t<att type=\"integer\" name=\"type\" value=\"%d\"/>\n",g[*ei].type);*/
			
			//fprintf(fp,"\t\t<att type=\"string\" name=\"descriptor\" value=\"%s\"/>\n",descriptor.c_str());
			fprintf(fp,"\t</edge>\n");
		}
	}
	fprintf(fp,"</graph>\n");
	fclose(fp);
}
void testKPLS()
{
	KPLS * kpls = new KPLS();
	kpls->SetLatentVars(5);
	kpls->SetSigma(20);

	int num_rows = 1000;
	int num_cols = 1;

	itk::Statistics::NormalVariateGenerator::Pointer nvg1, nvg2;
	nvg1 = itk::Statistics::NormalVariateGenerator::New();
	nvg2 = itk::Statistics::NormalVariateGenerator::New();

	nvg1->Initialize(1);
	nvg2->Initialize(2);




	MATRIX data = kpls->GetDataPtr(num_rows,num_cols);
	VECTOR ids = kpls->GetIDPtr();
	VECTOR training = kpls->GetTrainingPtr();


	int t = 7;

	for(int r = 0; r< num_rows; r++)
	{
		if(r < t)
		{
			data[r][0] = 1.0*rand()/RAND_MAX+2.0;
			ids[r] = r+1;
			training[r] = 1;
		}
		else if(r < 2*t)
		{
			data[r][0] = 10.0*rand()/RAND_MAX+3.0;
			ids[r] = r+1;
			training[r] = 2;
		}
		else if (r < (num_rows-2*t)/2+2*t)
		{
			data[r][0] = 1.0*rand()/RAND_MAX+2.0;
			ids[r] = r+1;
			training[r] = -1;
		}
		else
		{
			data[r][0] = 10.0*rand()/RAND_MAX+3.0;
			ids[r] = r+1;
			training[r] = -1;
		}
	}
	kpls->InitVariables();
	kpls->ScaleData();
	kpls->Train();
	kpls->Classify();

	VECTOR predictions = kpls->GetPredictions();
	int cor = 0,ncor = 0;
	for(int r = 2*t; r < num_rows; r++)
	{
		if(r < (num_rows-2*t)/2+2*t)
		{
			if(predictions[r] == 1)
				cor++;
			else
				ncor++;
		}
		else
		{
			if(predictions[r] == 2)
				cor++;
			else 
				ncor++;
		}
	}
	printf("It got cor = %d, ncor =  %d\n",cor, ncor);
	//scanf("%*d");
	exit(0);
}



void CellTracker::compute_feature_variances()
{
	TGraph::edge_iterator ei,eend;
	fvarnew =fvar;
	for(int counter = 0; counter < FeatureVariances::N; counter++)
	{
		fvarnew.variances[counter] = 0;
		fvarnew.means[counter] = 0;
	}
	fvarnew.distMean = 0;
	fvarnew.distVariance = 0;
	fvarnew.timeMean = 0;
	fvarnew.timeVariance = 0;
	fvarnew.overlapVariance = 0;
	fvarnew.overlapMean = 0;
	fvarnew.MS_prior = 0;
	fvarnew.AD_prior = 0;
	fvarnew.T_prior = 0;
	int num_samples = 0;
	for(tie(ei,eend)= edges(g); ei!=eend; ++ei)
	{
		int type = get_edge_type(*ei);
		num_samples++;
		switch(type)
		{
		case APPEAR: 
		case DISAPPEAR: fvarnew.AD_prior++;
			break;
		case SPLIT: 
		case MERGE: fvarnew.MS_prior++;
		case TRANSLATION: 
			FeaturesType f1,f2;
			f1 = fvector[g[source(*ei,g)].t][g[source(*ei,g)].findex];
			f2 = fvector[g[target(*ei,g)].t][g[target(*ei,g)].findex];
			float dist = get_distance(f1.Centroid,f2.Centroid); 
			fvarnew.distMean += dist;
			fvarnew.distVariance += dist*dist;
			fvarnew.timeMean += f2.time - f1.time;
			fvarnew.timeVariance += (f2.time - f1.time)*(f2.time - f1.time);
			float ovlap = 1-overlap(f1.BoundingBox,f2.BoundingBox)/MIN(f1.ScalarFeatures[FeaturesType::BBOX_VOLUME],f2.ScalarFeatures[FeaturesType::BBOX_VOLUME]);
			fvarnew.overlapMean += ovlap;
			fvarnew.overlapVariance +=ovlap*ovlap;
			for(int counter = 0; counter < FeatureVariances::N; counter++)
			{
				fvarnew.means[counter] += f2.ScalarFeatures[counter] - f1.ScalarFeatures[counter];
				fvarnew.variances[counter] += (f2.ScalarFeatures[counter] - f1.ScalarFeatures[counter])*(f2.ScalarFeatures[counter] - f1.ScalarFeatures[counter]);
			}
			if(type==TRANSLATION)
			{
				fvarnew.T_prior++;
			}
			break;
		}
	}
	for(int counter = 0; counter < FeatureVariances::N; counter++)
	{
		fvarnew.means[counter] /= num_samples;
		fvarnew.variances[counter] /=num_samples;
		fvarnew.variances[counter] -= fvarnew.means[counter]*fvarnew.means[counter];
	}
	fvarnew.distMean /=num_samples;
	fvarnew.distVariance /=num_samples;
	fvarnew.distVariance -= fvarnew.distMean*fvarnew.distMean;
	fvarnew.timeMean /= num_samples;
	fvarnew.timeVariance /= num_samples;
	fvarnew.timeVariance -= fvarnew.timeMean * fvarnew.timeMean;
	fvarnew.overlapMean /= num_samples;
	fvarnew.overlapVariance /= num_samples;
	fvarnew.overlapVariance -= fvarnew.overlapMean*fvarnew.overlapMean;
	fvarnew.MS_prior /= num_samples;
	fvarnew.AD_prior -=100;
	fvarnew.AD_prior /= num_samples;
	fvarnew.T_prior /= num_samples;
	float sum = fvarnew.T_prior + fvarnew.MS_prior + fvarnew.AD_prior;
	fvarnew.MS_prior /=sum;
	fvarnew.AD_prior /=sum;
	fvarnew.T_prior /=sum;
	

	printf("compute feature variances are : num_samples = %d\n", num_samples);
	for(int counter = 0; counter < FeatureVariances::N; counter++)
	{
		printf("ScalarFeatures[%d]. mean = %0.3f .variance = %0.3f\n",counter,fvarnew.means[counter],fvarnew.variances[counter]);
	}
	printf("distMean = %0.3f distVariance = %0.3f\n",fvarnew.distMean, fvarnew.distVariance);
	printf("timeMean = %0.3f timeVariance = %0.3f\n",fvarnew.timeMean, fvarnew.timeVariance);
	printf("overlapMean = %0.3f overlapVariance = %0.3f\n", fvarnew.overlapMean, fvarnew.overlapVariance);
	printf("MS_prior = %0.3f AD_prior = %0.3f T_prior = %0.3f\n", fvarnew.MS_prior, fvarnew.AD_prior, fvarnew.T_prior);
}

int optionsCreate(const char* optfile, map<std::string,std::string>& options)
{
  options.clear();
  std::ifstream fin(optfile); assert(fin.good());
  std::string name;  fin>>name;
  while(fin.good()) {
	 char cont[100];	 fin.getline(cont, 99);
	 options[name] = std::string(cont);
	 fin>>name;
  }
  fin.close();
  return 0;
}

int main(int argc1, char **argv1)
{
	FeatureVariances fvar;
	for(int counter=0; counter < FeatureVariances::N; counter++)
	{
		fvar.variances[counter] = std::numeric_limits<float>::max();
	}
	//if(argc1!=2)
	//	exit(0);

	fvar.distVariance = 50;//8.77;//50
	fvar.distMean = 5;//3.3;//5
	fvar.spacing[0] = 1;
	fvar.spacing[1] = 1;
	fvar.spacing[2] = 4;
	fvar.timeVariance = 1;//.119//1;
	fvar.overlapVariance = 1;//0.034;//1;
	fvar.overlapMean = 0;//0.2;//0;
	fvar.variances[FeatureVariances::VOLUME] = 90000;//44000;//90000;
	fvar.MS_prior = 0.4;//
	fvar.AD_prior = 0.01;
	fvar.T_prior = 1;
	fvar.boundDistMean = 4;
	fvar.boundDistVariance = 50;

	std::string debugprefix = "harvard";

	if(argc1==3)
	{
		 map<string, string> opts;  optionsCreate(argv1[2], opts);
  
		//get input data
		map<string,string>::iterator mi;
		mi = opts.find("x_spacing"); 
		if(mi!=opts.end()){ istringstream ss((*mi).second); ss>>fvar.spacing[0]; }

		mi = opts.find("y_spacing");
		if(mi!=opts.end()){ istringstream ss((*mi).second); ss>>fvar.spacing[1]; }
		
		mi = opts.find("z_spacing"); 
		if(mi!=opts.end()){ istringstream ss((*mi).second); ss>>fvar.spacing[2]; }
		
		mi = opts.find("distVariance");
		if(mi!=opts.end()){ istringstream ss((*mi).second); ss>>fvar.distVariance; }

		mi = opts.find("distMean"); 
		if(mi!=opts.end())
		{ istringstream ss((*mi).second); ss>>fvar.distMean; }

		mi = opts.find("boundaryDistVariance");
		if(mi!=opts.end()){ istringstream ss((*mi).second); ss>>fvar.boundDistVariance; }

		mi = opts.find("boundaryDistMean"); 
		if(mi!=opts.end()){ istringstream ss((*mi).second); ss>>fvar.boundDistMean; }

		mi = opts.find("timeVariance"); 
		if(mi!=opts.end()){ istringstream ss((*mi).second); ss>>fvar.timeVariance; }

		mi = opts.find("timeMean"); 
		if(mi!=opts.end()){ istringstream ss((*mi).second); ss>>fvar.timeMean; }

		mi = opts.find("overlapVariance"); 
		if(mi!=opts.end()){ istringstream ss((*mi).second); ss>>fvar.overlapVariance; }

		mi = opts.find("overlapMean"); 
		if(mi!=opts.end()){ istringstream ss((*mi).second); ss>>fvar.overlapMean; }

		mi = opts.find("volumeVariance");
		if(mi!=opts.end()){ istringstream ss((*mi).second); ss>>fvar.variances[FeatureVariances::VOLUME]; }

		mi = opts.find("merge_split_prior"); 
		if(mi!=opts.end()){ istringstream ss((*mi).second); ss>>fvar.MS_prior; }

		mi = opts.find("appear_disappear_prior"); 
		if(mi!=opts.end()){ istringstream ss((*mi).second); ss>>fvar.AD_prior; }

		mi = opts.find("translation_prior"); 
		if(mi!=opts.end()){ istringstream ss((*mi).second); ss>>fvar.T_prior; }

		mi = opts.find("debugprefix"); 
		if(mi!=opts.end()){ istringstream ss((*mi).second); ss>>debugprefix; }

		printf("fvar.spacing[0] = %f\n",fvar.spacing[0]);
		printf("fvar.spacing[1] = %f\n",fvar.spacing[1]);
		printf("fvar.spacing[2] = %f\n",fvar.spacing[2]);
		printf("fvar.distVariance = %f\n",fvar.distVariance);
		printf("fvar.distMean = %f\n",fvar.distMean);
		printf("fvar.boundDistVariance = %f\n",fvar.boundDistVariance);
		printf("fvar.boundDistMean = %f\n",fvar.boundDistMean);
		printf("fvar.timeVariance = %f\n", fvar.timeVariance);
		printf("fvar.timeMean = %f\n", fvar.timeMean);
		printf("fvar.overlapVariance = %f\n", fvar.overlapVariance);
		printf("fvar.overlapMean = %f\n", fvar.overlapMean);
		printf("fvar.volumeVariance = %f\n", fvar.variances[FeatureVariances::VOLUME]);
		printf("fvar.MS_prior = %f\n", fvar.MS_prior);
		printf("fvar.AD_prior = %f\n", fvar.AD_prior);
		printf("fvar.T_prior = %f\n", fvar.T_prior);
		printf("debugprefix = %s\n",debugprefix.c_str());
		//scanf("%*d");

	}
	else if(argc1>=4)
	{
		printf("Usage: %s tracking_filenames_file [tracking_parameters_file]",argv1[0]);
		scanf("%*d");
		return -1;
	}
	
	char buff_ar[1024];

	std::ifstream ff;
	ff.open(argv1[1],std::ios::in);
	int num_lines_count = 0;
	int argc = 1;
	char **argv = new char *[10000];
	argv[0] = new char[1];
	std::string ffbuf;
	char ffbuffer[1024];
	/*ff.getline(ffbuffer,1024);
	float spac[3];
	if(ffbuffer[0]!='\0')
	{
		sscanf(ffbuffer,"%f %f %f",spac,spac+1,spac+2);
		printf("I got spacing = %f %f %f\n", spac[0],spac[1],spac[2]);
	//	scanf("%*d");
	}*/
	while(!ff.eof())
	{	
		//ff.getline(ffbuffer,1024);
		ff.getline(ffbuffer,1024);
		//ff>>ffbuf;
		if(ff.gcount()==0)
			break;
		if(ffbuffer[strlen(ffbuffer)-1]=='\r')
			ffbuffer[strlen(ffbuffer)-1]=0;
		argv[argc] = new char [1024];
		strcpy(argv[argc],ffbuffer);
		//printf("|%s|\n",argv[argc]);
		argc++;
		
	}

	ff.close();
	printf("argc = %d",argc);
	//scanf("%*d");
//	return 0;

	
	/*//ST();
	// testKPLS();
	int num_tc = 10;
//#define BASE "C:\\Users\\arun\\Research\\Tracking\\harvard\\cache\\testTSeries-02102009-1455-624"
#define BASE "L:\\Tracking\\cache\\testTSeries-02102009-1455-624"
	int argc = num_tc*3+1+3;

	char ** argv = new char* [argc];
	int ch = 3;
	for(int counter=1; counter <argc; counter++)
	{
		argv[counter] = new char [1024];
	}
	sprintf(argv[1],BASE);
	for(int counter =1; counter<=num_tc; counter++)
	{
		sprintf(argv[1+counter],"smoothed_TSeries-02102009-1455-624_Cycle%03d_CurrentSettings_Ch%d.tif",counter,ch);
	}
	sprintf(argv[num_tc+2],BASE);
	for(int counter =1; counter<=num_tc; counter++)
	{
		sprintf(argv[counter+num_tc+2],"clabeled_TSeries-02102009-1455-624_Cycle%03d_CurrentSettings_Ch%d.tif",counter,ch);
	}
	sprintf(argv[2*num_tc+3],BASE);
	for(int counter =1; counter<=num_tc; counter++)
	{
		sprintf(argv[counter+num_tc*2+3],"labeled_tracks_TSeries-02102009-1455-624_Cycle%03d_CurrentSettings_Ch%d.tif",counter,ch);
	}

	//END DEBUG
*/
	std::vector<std::string> imfilenames,labelfilenames,outputfilenames;



	printf("Started\n");
	/*int num_files = atoi(argv[1]);
	int c;*/
	//int counter = 0;
	printf(" I got argc = %d\n",argc);
	int num_t = (argc-4)/3;

	for(int counter = 0; counter < num_t; counter++)
	{
		std::string tbuff = std::string(argv[1]) + std::string("/") + std::string(argv[2+counter]);
		imfilenames.push_back(tbuff);
	}
	for(int counter = 0; counter < num_t; counter++)
	{
		std::string tbuff = std::string(argv[2+num_t]) + std::string("/") +std::string(argv[3+num_t+counter]);
		labelfilenames.push_back(tbuff);
	}
	for(int counter = 0; counter < num_t; counter++)
	{
		std::string tbuff = std::string(argv[3+2*num_t]) + std::string("/") +std::string(argv[4+2*num_t+counter]);
		outputfilenames.push_back(tbuff);
	}


	/*for(int counter=0; counter < num_t; counter++)
	{
	printf("---t = %d %s %s %s----\n", counter, imfilenames[counter].c_str(),labelfilenames[counter].c_str(), outputfilenames[counter].c_str());
	}*/
	//scanf("%*d");
	/*
	char file_name[100][1024];
	char str[1024];
	int num_c= atoi(argv[5]);
	printf("num_files,num_t:%d %d\n",num_files,num_t);
	*/

	/*	
	FILE * fp = fopen(argv[4],"r");
	printf("filename:%s\n",argv[4]);
	c = atoi(argv[2]);*/
	//check num_files

	/*for(int i=0;i<(3*num_files);i++)
	{
	fgets(file_name[i],1024,fp);
	file_name[i][strlen(file_name[i])-1]= '\0';
	printf("files:%s\n",file_name[i]);
	if(feof(fp))
	break;
	}
	*/


	//LabelImageType::Pointer segmented[MAX_TIME][MAX_TAGS]; // FIXME
	//InputImageType::Pointer images[MAX_TIME][MAX_TAGS];
	int c=1;



	//Now track the cells alone and compute associations

	std::vector<std::vector<FeaturesType> > fvector;
	CellTracker::VVL limages;
	//FeaturesType *farray;

	LabelImageType::Pointer labeled;
	bool fexists = true;
	for(int counter=0; counter<num_t; counter++)
	{
		if(!file_exists((char*)outputfilenames[counter].c_str()))
		{
			fexists = false;
			break;
		}
	}
	if(1)//!fexists)
	{
		CellTracker::VVL limages;
		CellTracker::VVR rimages;
		LabelImageType::Pointer tempsegmented;
		InputImageType::Pointer tempimage;
		std::vector<Color2DImageType::Pointer> input,output;
		char *filename_number = "C:\\Users\\arun\\Research\\Tracking\\harvard\\numbers.bmp";
		Color2DImageType::Pointer number = readImage<Color2DImageType>(filename_number);
		fvar.time_last = num_t-1;
		for(int t =0; t<num_t; t++)
		{
			tempimage = readImage<InputImageType>(imfilenames[t].c_str());	
			tempsegmented = readImage<LabelImageType>(labelfilenames[t].c_str());
			//tempsegmented = getLargeLabels(tempsegmented,100);
			Color2DImageType::Pointer cimp = getColor2DImage(tempsegmented,2);
			std::vector<FeaturesType> f;
			std::vector<LabelImageType::Pointer> li;
			std::vector<InputImageType::Pointer> ri;
			getFeatureVectorsFarsight(tempsegmented,tempimage,f,t,c);
			for(int counter=0; counter < f.size(); counter++)
			{
				li.push_back(extract_label_image(f[counter].num,f[counter].BoundingBox,tempsegmented));
				ri.push_back(extract_raw_image(f[counter].BoundingBox,tempimage));
				annotateImage(number,cimp,f[counter].num,MAX(f[counter].Centroid[0],0),MAX(f[counter].Centroid[1],0));
				if(f[counter].ScalarFeatures[FeaturesType::VOLUME]<5)
				{
					printf("A cell has a volume of %f\n",f[counter].ScalarFeatures[FeaturesType::VOLUME]);
					printFeatures(f[counter]);
					//PAUSE;
				}
			}
			input.push_back(cimp);
			fvector.push_back(f);
			limages.push_back(li);
			rimages.push_back(ri);
		}

		ColorImageType::Pointer colin = getColorImageFromColor2DImages(input);
		//writeImage<ColorImageType>(colin,"C:\\Users\\arun\\Research\\Tracking\\harvard\\debug\\input.tif");
		std::string debugfolder = "C:\\Users\\arun\\Research\\Tracking\\harvard\\debug\\";
		std::string debugstring = debugfolder + debugprefix + "_input.tif";
		writeImage<ColorImageType>(colin,debugstring.c_str());

		InputImageType::SizeType imsize = tempimage->GetLargestPossibleRegion().GetSize();
		fvar.BoundingBox[0] = 0;
		fvar.BoundingBox[2] = 0;
		fvar.BoundingBox[4] = 0;
		fvar.BoundingBox[1] = imsize[0]-1;
		fvar.BoundingBox[3] = imsize[1]-1;
		fvar.BoundingBox[5] = imsize[2]-1;


		CellTracker ct(fvector,fvar,limages,rimages);
		ColorImageType::Pointer debugcol1 = ColorImageType::New();
		ColorImageType::Pointer debugcol2 = ColorImageType::New();
		ColorImageType::Pointer debugcol3 = ColorImageType::New();
		ColorImageType::SizeType colsize;
		ColorImageType::RegionType colregion;
		ColorImageType::IndexType colindex;
		colindex.Fill(0);
		colsize[0] = imsize[0]*1;
		colsize[1] = imsize[1]*1;
		colsize[2] = num_t;
		colregion.SetIndex(colindex);
		colregion.SetSize(colsize);
		debugcol1->SetRegions(colregion);
		debugcol1->Allocate();
		debugcol2->SetRegions(colregion);
		debugcol2->Allocate();
		debugcol3->SetRegions(colregion);
		debugcol3->Allocate();
		ct.set_debug_images(debugcol1,debugcol2,debugcol3);
		ct.run();
		printf("Rerunning with computed variances\n");
		FeatureVariances fvarnew(ct.get_computed_variances());
		fvarnew.MS_prior = 0.4;
		fvarnew.AD_prior = 0.01;
		fvarnew.T_prior = 1;
		fvarnew.timeVariance = 1;
		fvarnew.overlapVariance = 1;
		//CellTracker ct1(fvector,fvarnew,limages,rimages);
		//ct1.set_debug_images(debugcol1,debugcol2);
		//ct.writeGraphViz("C:\\Users\\arun\\Research\\Tracking\\harvard\\graphviz_testSeries.dot");
		//ct.writeGraphML("C:\\Users\\arun\\Research\\Tracking\\harvard\\graphviz_testSeries.graphml");
		//PAUSE;
		//ct1.run();
		for(int t = 0; t< num_t; t++)
		{
			printf("In final loop t = %d\n",t);
			//tempsegmented = readImage<LabelImageType>(argv[(t+1)+num_t]);
			//tempsegmented = getLargeLabels(tempsegmented,100);
			//LabelImageType::Pointer track = LabelImageType::New();
			//track->SetRegions(tempsegmented->GetLargestPossibleRegion());
			//track->Allocate();
			//track->FillBuffer(0);
			//LabelIteratorType liter1(tempsegmented,tempsegmented->GetLargestPossibleRegion());
			//LabelIteratorType liter2(track,track->GetLargestPossibleRegion());
			//liter1.GoToBegin();
			//liter2.GoToBegin();
			//std::set<int> s;
			//for(;!liter1.IsAtEnd();++liter1,++liter2)
			//{
			//	int val = liter1.Get();
			//	if(ct.old_to_new[t].count(val)!=0)
			//	{
			//		if(ct.old_to_new[t][val] !=0)
			//			liter2.Set(ct.old_to_new[t][val]);
			//		else
			//			s.insert(val);
			//	}
			//	else
			//	{
			//		s.insert(val);
			//	}
			//}
			LabelImageType::Pointer track = ct.getOutputAtTime(t);
			Color2DImageType::Pointer cimp = getColor2DImage(track,2);
			std::vector<FeaturesType> f;
			getFeatureVectorsFarsight(track,tempimage,f,t,c);
			printf("About to begin annotate loop\n");
			for(int counter=0; counter< f.size(); counter++)
			{
				std::vector<FeaturesType> conncomp = get_all_connected_components(track,f[counter]);
				for(int counter1 = 0; counter1 < conncomp.size(); counter1++)
				{
					//annotateImage(number,cimp,f[counter].num,MAX(f[counter].Centroid[0],0),MAX(f[counter].Centroid[1],0));
					annotateImage(number,cimp,f[counter].num, MAX(conncomp[counter1].Centroid[0],0),MAX(conncomp[counter1].Centroid[1],0));
				}
			}
			//PAUSE;
			printf("Finished annotate loop\n");
			output.push_back(cimp);
			//std::cout<<"could not find :" ;
			//copy(s.begin(), s.end(), std::ostream_iterator<int>(std::cout, " "));
			printf("About to call writeImage\n");
			writeImage<LabelImageType>(track,outputfilenames[t].c_str());
		}
		
		ColorImageType::Pointer colout = getColorImageFromColor2DImages(output);
		debugstring = debugfolder + debugprefix + "_debugcol1.tif";
		writeImage<ColorImageType>(debugcol1,debugstring.c_str());
		debugstring = debugfolder + debugprefix + "_debugcol2.tif";
		writeImage<ColorImageType>(debugcol2,debugstring.c_str());
		debugstring = debugfolder + debugprefix + "_debugcol3.tif";
		writeImage<ColorImageType>(debugcol3,debugstring.c_str());
		debugstring = debugfolder + debugprefix + "_output.tif";
		writeImage<ColorImageType>(colout,debugstring.c_str());
		

	}

	//  //for(int c=1; c<=2; c++)
	//	{
	//		int cur_max = 0;
	//		printf("Finished computing the features\n");
	//		unsigned int * track_labels = new unsigned int[MAX_LABEL];
	//		for(unsigned int cind =0; cind< fvector[0][c-1].size(); cind++)
	//		{
	//			track_labels[fvector[0][c-1][cind].num-1] = cind;
	//			printf("track_labels[%d] = %d\n",cind,fvector[0][c-1][cind].num);
	//		}
	//		segmented[0][c-1] = getLabelsMapped(segmented[0][c-1],fvector[0][c-1],track_labels);
	//		cur_max = fvector[0][c-1].size();
	//		delete [] track_labels;

	//		//sprintf(tracks_cache,"%s/%s",CACHE_PREFIX,file_name[(num_files*2)+(c-1)*num_t]);
	//		writeImage<LabelImageType>(segmented[0][c-1],argv[0+1+2*num_t]);

	//		printf("cur_max  = %d\n",cur_max);
	//		//cur_max is the next allowed track number
	//		for(int t = 0; t< num_t-1; t++)
	//		{
	//			vcl_vector<unsigned> indices = getTimeAssociations(fvector[t+1][c-1],fvector[t][c-1]);
	//			printf("fvector@t.size %d fvector@t+1.size %d indices.size %d\n",
	//              (int)fvector[t][c-1].size(), (int)fvector[t+1][c-1].size(),
	//              (int)indices.size());
	//			track_labels = new unsigned int[MAX_LABEL];
	//			for(unsigned int cind = 0; cind<indices.size(); cind++)
	//			{
	//				if(indices[cind] > 100000)
	//				{
	//					track_labels[fvector[t+1][c-1][cind].num-1] = cur_max++;
	//				}
	//				else
	//				{
	//					track_labels[fvector[t+1][c-1][cind].num-1] = fvector[t][c-1][indices[cind]].num-1;
	//				}
	//				//	printf("track_labels[%d] = %u\n",cind, track_labels[cind]);
	//			}

	//			segmented[t+1][c-1] = getLabelsMapped(segmented[t+1][c-1],fvector[t+1][c-1],track_labels);
	//			for(unsigned int cind =0; cind < fvector[t+1][c-1].size(); cind++)
	//			{
	//				printf("fvector[%d].num = %d\n",cind, fvector[t+1][c-1][cind].num);
	//			}
	//			//sprintf(tracks_cache,"%s/%s",CACHE_PREFIX,file_name[(num_files*2)+((c-1)*num_t)+(t+1)]);
	//			writeImage<LabelImageType>(segmented[t+1][c-1],argv[t+1+1+2*num_t]);
	//       //scanf("%*d");
	//			delete [] track_labels;
	//		}
	//		printf("Cur_max %d\n",cur_max);
	//	}
	//}


	//scanf("%*d");
	return 0;
}
