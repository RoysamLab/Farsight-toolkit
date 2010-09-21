#include "helpers.h"

#include  <time.h>
#include <iostream>                  
#include <utility>                   
#include <algorithm>
#include <map>
#include <set>
#include <limits>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
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
	};

	struct TrackEdge{
		int utility;
		bool coupled;
		bool fixed;
		bool selected;
		unsigned char type;
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
		UTILITY_MAX = 1<<16;
		K = 1;
	}

	void run();
	void writeGraphViz(char * filename);
	void writeGraphML (char * filename);
	std::vector<std::map<int,int>> old_to_new;
	void set_debug_image(ColorImageType::Pointer in)
	{
		debugimage = in;
	}
	LabelImageType::Pointer getOutputAtTime(int t);
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
	float compute_LRUtility(FeaturesType f1, FeaturesType f2, FeaturesType f3);
	void enforce_overlap();
	int is_overlapping(TGraph::edge_descriptor e);
	std::pair<TGraph::edge_descriptor,bool> my_add_edge(TGraph::vertex_descriptor v1, TGraph::vertex_descriptor v2, int utility, bool coupled, bool fixed, bool selected, unsigned char type);
	void draw_line_for_edge(TGraph::edge_descriptor e,VectorPixelType col1,VectorPixelType col2,int);
	void print_vertex(TGraph::vertex_descriptor v,int);
	bool is_merge_node(TGraph::vertex_descriptor vd);
	bool is_split_node(TGraph::vertex_descriptor vd);
	bool is_simple_node(TGraph::vertex_descriptor vd);
	TGraph::vertex_descriptor get_parent(TGraph::vertex_descriptor);// no error check implemented TODO
	TGraph::vertex_descriptor get_child(TGraph::vertex_descriptor);// no error check implemented TODO
	bool is_separate(CellTracker::TGraph::vertex_descriptor v1, CellTracker::TGraph::vertex_descriptor v2);
	float get_LRUtility(std::vector< TGraph::vertex_descriptor > desc);


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
	int K;

	ColorImageType::Pointer debugimage;


	struct LREdge{
		TGraph::edge_descriptor ef,eb;
	};

};



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
		overlap = (ex-sx)*(ey-sy)*(ez-sz);
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
	float dist = get_boundary_dist(x)+3*sqrt(fvar.distVariance);
	return UTILITY_MAX*(exp(-dist*dist/2.0/fvar.distVariance));
}

int CellTracker::compute_normal_utility(FeaturesType f1, FeaturesType f2)
{
	float utility = 0;
	for(int counter=0; counter< FeaturesType::N; counter++)
	{
		utility += (f1.ScalarFeatures[counter]-f2.ScalarFeatures[counter])*(f1.ScalarFeatures[counter]-f2.ScalarFeatures[counter])/fvar.variances[counter];
	}

	float dist = get_distance(f1.Centroid,f2.Centroid);
	utility += dist*dist/fvar.distVariance;
	utility += (abs(f1.time-f2.time)-1)*(abs(f1.time-f2.time)-1)/fvar.timeVariance;
	float ovlap = overlap(f1.BoundingBox,f2.BoundingBox);
	ovlap = 1-(ovlap)/MIN(f1.ScalarFeatures[FeaturesType::BBOX_VOLUME],f2.ScalarFeatures[FeaturesType::BBOX_VOLUME]);
	utility += ovlap*ovlap/fvar.overlapVariance;
	utility /= 2.0;
	utility = UTILITY_MAX*(exp(-utility));
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

	vertex_iter vi, vend;
	bool added;
	int ret_count = 0;
	for(tie(vi,vend) = vertices(g); vi != vend; ++vi)
	{
		//printf("hi t = %d\n",t);
		if(g[*vi].t == t-1)
		{
			//printf("hi 1\n");
			if(g[*vi].special == 0)
			{
				//printf("hi 2\n");
				if(get_boundary_dist(fvector[t-1][g[*vi].findex].Centroid) < 4.0*sqrt(fvar.distVariance))
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
				}
			}
		}
	}
	_LTRACE;
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
				if(get_boundary_dist(fvector[t+1][g[*vi].findex].Centroid) <  4.0*sqrt(fvar.distVariance))
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
					_TRACE;
				}
			}
		}
	}
	_LTRACE;
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
					if(dist < 4.0*sqrt(fvar.distVariance))
					{
						//add the edge
						if(dist>100)
							printf("distance is greater than 100\n");
						tie(e,added) = add_edge(rmap[t][counter1],v,g);
						if(added)
						{
							g[e].coupled = 0;
							g[e].fixed = 0;
							g[e].selected = 0;
							g[e].utility = compute_normal_utility(fvector[t][counter1],fvector[tmax][counter]);
							nec++;
						}
					}
				}
			}
		}
	}
	_LTRACE;
	printf("add_normal_edges: I tried %d %d\n",tried1, tried2);
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

					lc1.ScalarFeatures[FeaturesType::VOLUME] = lc1.ScalarFeatures[FeaturesType::VOLUME] + lc2.ScalarFeatures[FeaturesType::VOLUME];
					lc2.ScalarFeatures[FeaturesType::VOLUME] = lc1.ScalarFeatures[FeaturesType::VOLUME];

					int utility1 = compute_normal_utility(lc1,fvector[tmax][counter1]);
					int utility2 = compute_normal_utility(lc2,fvector[tmax][counter1]);

					tie(e1,added1) = add_edge(rmap[tcounter][i1],rmap[tmax][counter1],g);
					tie(e2,added2) = add_edge(rmap[tcounter][i2],rmap[tmax][counter1],g);
					if(added1&&added2)
					{
						g[e1].coupled = 1;
						g[e2].coupled = 1;
						coupled_map[e1] = e2;
						coupled_map[e2] = e1;
						g[e1].utility = (utility1+utility2);
						g[e2].utility = (utility1+utility2);
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

					lc1.ScalarFeatures[FeaturesType::VOLUME] = lc1.ScalarFeatures[FeaturesType::VOLUME] + lc2.ScalarFeatures[FeaturesType::VOLUME];
					lc2.ScalarFeatures[FeaturesType::VOLUME] = lc1.ScalarFeatures[FeaturesType::VOLUME];

					int utility1 = compute_normal_utility(lc1,fvector[tcounter][counter1]);
					int utility2 = compute_normal_utility(lc2,fvector[tcounter][counter1]);

					tie(e1,added1) = add_edge(rmap[tcounter][counter1],rmap[tmax][i1],g);
					tie(e2,added2) = add_edge(rmap[tcounter][counter1],rmap[tmax][i2],g);
					if(added1&&added2)
					{
						g[e1].coupled = 1;
						g[e2].coupled = 1;
						coupled_map[e1] = e2;
						coupled_map[e2] = e1;
						g[e1].utility = (utility1+utility2);
						g[e2].utility = (utility1+utility2);
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


float CellTracker::compute_LRUtility(FeaturesType f1, FeaturesType f2, FeaturesType f3)
{
	float utility = 0;
	utility += compute_normal_utility(f1,f2);
	utility += compute_normal_utility(f2,f3);


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

		float proj = (1+dir1[0]*dir2[0]+dir1[1]*dir2[1]+dir1[2]*dir2[2]);

		utility = MIN(dist1,dist2)/MAX(dist1,dist2)*proj;
		return utility;
	}

}
float CellTracker::compute_LRUtility(TGraph::edge_descriptor e1,TGraph::edge_descriptor e2)
{
	float utility = 0;
	utility = (g[e1].utility + g[e2].utility);
	if(g[e1].type == APPEAR || g[e2].type == DISAPPEAR)
	{
		return utility;
	}
	FeaturesType f1,f2,f3;
	f1 = fvector[g[source(e1,g)].t][g[source(e1,g)].findex];
	f2 = fvector[g[source(e2,g)].t][g[source(e2,g)].findex];
	f3 = fvector[g[target(e2,g)].t][g[target(e2,g)].findex];
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

		float proj = dir1[0]*dir2[0]+dir1[1]*dir2[1]+dir1[2]*dir2[2];

		utility = MIN(dist1,dist2)/MAX(dist1,dist2)*proj;
		return utility;
	}

}



void CellTracker::draw_line_for_edge(TGraph::edge_descriptor e,VectorPixelType col1,VectorPixelType col2, int shift = 0)
{
	TGraph::vertex_descriptor v1,v2;
	v1 = source(e,g);
	v2 = target(e,g);
	if(g[e].selected==1)
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
			drawLine(debugimage,col1,col2,f1.Centroid[0]+shift,f1.Centroid[1],f1.time,f2.Centroid[0]+shift,f2.Centroid[1],f1.time);
			//drawLine(debugimage,col2,f1.Centroid[0]+shift,f1.Centroid[1],f2.time,f2.Centroid[0]+shift,f2.Centroid[1],f2.time);
		}
	}
}
void CellTracker::solve_higher_order()
{

	VectorPixelType col1,col2,col3,col4;
	col1[0] = 255;col1[1] = 0;col1[2] = 0;
	col2[0] = 255;col2[1] = 255;col2[2] = 0;
	col3[0] = 255;col3[1] = 255;col3[2] = 255;
	col4[0] = 255;col4[1] = 0; col4[2] = 255;

	//TODO - FIXME - work in progress - dont use it yet.
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
					//draw a line for the edge
					if(g[*ei].utility >0)
					{
						if(g[*ei].coupled==0)
						{
							//draw_line_for_edge(*ei,col1,col2,-2);
						}
						else
						{
							//draw_line_for_edge(*ei,col3,col4,2);
						}
					}
					//
					/*var_index[*ei] = varc;
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
					*/
				}
			}
		}

		printf("ecount = %d varc = %d\n",ecount,varc);
		graph_traits<TGraph>::vertex_iterator vi,vend;



		for(tie(vi,vend) = vertices(g); vi != vend; ++vi)
		{
			if(in_degree(*vi,g) ==0 || out_degree(*vi,g)==0)
			{
				// no need to form any new LREdge with this vertex
				continue;
			}
			graph_traits<TGraph>::in_edge_iterator e_in,e_in_end;
			graph_traits<TGraph>::out_edge_iterator e_out,e_out_end;

			tie(e_in,e_in_end) = in_edges(*vi,g);
			tie(e_out,e_out_end) = out_edges(*vi,g);

			for(;e_in!=e_in_end; ++e_in)
			{
				for(; e_out != e_out_end; ++e_out)
				{
					//handle lotsa cases
					int tin, tout;
					tin = get_edge_type(*e_in);
					tout = get_edge_type(*e_out);
					int utility = compute_LRUtility(*e_in,*e_out);
					if(utility > 0)
					{
						;
					}
				}
			}
		}
		//IloBoolVarArray x(env,varc);

		//for(int counter=0; counter < varc; counter++)
		//{
		//	obj.setLinearCoef(x[counter],utility[counter]);
		//	//printf("utility[%d] = %d\n",counter,utility[counter]);
		//}

		//int vcount = -1;
		//for(tie(vi,vend) = vertices(g); vi != vend; ++vi)
		//{

		//	graph_traits<TGraph>::in_edge_iterator e_i,e_end;
		//	tie(e_i,e_end) = in_edges(*vi,g);
		//	bool once = false;
		//	float fixed_sum = 0;
		//	for(;e_i!=e_end; ++e_i)
		//	{
		//		if(g[*e_i].fixed == 0 )
		//		{
		//			if(once == false )
		//			{
		//				once = true;
		//				vcount++;
		//				c.add(IloRange(env,0,1));
		//			}
		//			c[vcount].setLinearCoef(x[var_index[*e_i]],1.0);
		//		}
		//		else
		//		{
		//			if(g[*e_i].coupled==0)
		//				fixed_sum += g[*e_i].selected;
		//			else
		//				fixed_sum += 0.5*g[*e_i].selected;
		//		}
		//	}
		//	if(once)
		//		c[vcount].setBounds(0,1-fixed_sum);
		//	graph_traits<TGraph>::out_edge_iterator e_2,e_end2;
		//	tie(e_2,e_end2) = out_edges(*vi,g);
		//	once = false;
		//	fixed_sum = 0;
		//	for(;e_2!=e_end2; ++e_2)
		//	{
		//		if(g[*e_2].fixed == 0 )
		//		{
		//			if(once == false )
		//			{
		//				once = true;
		//				vcount++;
		//				c.add(IloRange(env,0,1));
		//			}
		//			c[vcount].setLinearCoef(x[var_index[*e_2]],1.0);
		//		}
		//		else
		//		{
		//			if(g[*e_2].coupled==0)
		//				fixed_sum += g[*e_2].selected;
		//			else
		//				fixed_sum += 0.5*g[*e_2].selected;
		//		}
		//	}
		//	if(once)
		//		c[vcount].setBounds(0,1-fixed_sum);

		//}
		//printf("vcount = %d\n",vcount);
		////scanf("%*d");
		//IloModel model(env);
		//model.add(obj);
		//model.add(c);


		////std::cout<<model<<std::endl;
		//IloCplex cplex(model);
		//if(!cplex.solve())
		//{
		//	std::cerr << " Could not solve.. error"<<std::endl;
		//}
		//IloNumArray vals(env);
		//env.out() << "Solution status = " << cplex.getStatus() << std::endl;
		//env.out() << "Solution value  = " << cplex.getObjValue() << std::endl;
		//cplex.getValues(vals, x);
		////env.out() << "Values        = " << vals << std::endl;
		//std::cout << "Values.getSize() = "<< vals.getSize() << std::endl;
		//for(int counter=0; counter< vals.getSize(); counter++)
		//{
		//	//assert(g[inv_index[counter]].selected == 1);
		//	g[inv_index[counter]].selected = vals[counter];
		//	if(g[inv_index[counter]].coupled==1)
		//		g[coupled_map[inv_index[counter]]].selected = vals[counter];
		//	if(vals[counter]==1 && utility[counter] <0)
		//	{
		//		printf("wrong @ %d\n",counter);
		//	}
		//}
	}
	catch(IloException &e)
	{
		std::cerr << e << std::endl;
	}
	env.end();
}

std::pair<CellTracker::TGraph::edge_descriptor,bool> CellTracker::my_add_edge(TGraph::vertex_descriptor v1, TGraph::vertex_descriptor v2, int utility, bool coupled, bool fixed, bool selected, unsigned char type)
{
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
void CellTracker::enforce_overlap()
{
	printf("Entered enforce overlap\n");
	// if A->B has overlap of bounding boxes, then all out_edges of A should have overlap and all in_edges of B should have overlap.

	TGraph::vertex_iterator vi,vend;

	for(tie(vi,vend) = vertices(g); vi!=vend; ++vi)
	{
		// in edges
		TGraph::in_edge_iterator e_in, e_in_end,e_in2,e_in_end2;
		tie(e_in,e_in_end) = in_edges(*vi,g);
		bool isovlap = false;
		//printf("Ein1\n");
		for(;e_in!=e_in_end; ++e_in)
		{
			if(is_overlapping(*e_in)==1)
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
				if(is_overlapping(*e_in2)==0)
				{
					g[*e_in2].utility = -(UTILITY_MAX-1);
					if(g[*e_in2].coupled==1)
						g[coupled_map[*e_in2]].utility = -(UTILITY_MAX-1);
				}
			}
			//printf("Finished E_in2\n");
			TGraph::vertex_descriptor vsource = source(*e_in,g);
			TGraph::out_edge_iterator e_out, e_out_end;
			tie(e_out,e_out_end) = out_edges(vsource,g);
			//printf("E_out\n");
			for(;e_out!=e_out_end; ++e_out)
			{
				if(is_overlapping(*e_out)==0)
				{
					g[*e_out].utility = -(UTILITY_MAX-1);
					if(g[*e_out].coupled==1)
						g[coupled_map[*e_out]].utility = -(UTILITY_MAX -1);
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
		return -1;
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
		return -1;
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
	return  rand()%1000;
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
	enforce_overlap();
	solve_higher_order();
	solve_lip();
	print_stats();

	boost::graph_traits<TGraph>::edge_iterator e_i,e_next,e_end;
	for(tie(e_i,e_end) = edges(g); e_i!= e_end ; e_i = e_next)
	{
		e_next = e_i;
		++e_next;
		if(g[*e_i].selected == 0)
		{

			remove_edge(*e_i,g);
		}
	}

	for(tie(e_i,e_end) = edges(g); e_i!= e_end ; ++e_i)
	{
		TGraph::vertex_descriptor vert1 = source(*e_i,g),vert2 = target(*e_i,g);
		if(g[vert1].t >= g[vert2].t)
		{
			printf("problem here:");
			print_vertex(vert1,0);
			print_vertex(vert2,0);
			//scanf("%*d");
		}
	}



	boost::graph_traits<TGraph>::vertex_iterator v_i,v_next,v_end,v_temp;

	for(tie(v_i,v_end) = vertices(g);v_i!=v_end; ++v_i)
	{
		if(in_degree(*v_i,g) >2 || out_degree(*v_i,g) >2)
		{
			printf("problem here2:");
			print_vertex(*v_i,0);
			//scanf("%*d");
		}
	}

	std::vector< std::vector < boost::graph_traits<TGraph>::vertex_descriptor > > to_resolve;


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
						if(vdt == -1)
							break;
						vvd.push_back(vdt);
						tr_map[vdt] = true;
					}while(is_simple_node(vdt));
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
		printf("Is_separate? : %d\n", int(is_separate(rmap[3][11],rmap[3][15])));
		printf("Is_separate? : %d\n", int(is_separate(rmap[3][6],rmap[3][8])));
		printf("Is_separate? : %d\n", int(is_separate(rmap[3][3],rmap[3][4])));
		printf("Is_separate? : %d\n", int(is_separate(rmap[3][10],rmap[3][17])));
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
				change_made = true;	
				std::vector < TGraph::vertex_descriptor > sp = to_resolve[counter];
				std::vector < std::pair < TGraph::vertex_descriptor, TGraph::vertex_descriptor > > spvertices;
				for(int cosp =0; cosp < sp.size(); cosp++)
				{
					TGraph::vertex_descriptor vd = sp[cosp];

					std::vector<LabelImageType::Pointer> lout;
					std::vector<InputImageType::Pointer> rout;
					std::vector<FeaturesType> fvecout;
					_TRACE;
					SplitCell(limages[g[vd].t][g[vd].findex],rimages[g[vd].t][g[vd].findex],fvector[g[vd].t][g[vd].findex],fvar,lout,rout,fvecout);
					_TRACE;

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
				float utilpos = -1;
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

					float util = get_LRUtility(desc);
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
						desc.push_back(spvertices[0].first);
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
					util += get_LRUtility(desc);
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
				if(before.size()!=0)
				{
					if(permutations[utilpos][permpc++] == 0)
					{
						tie(e1,added) = my_add_edge(before[before.size()-1].first,spvertices[0].first,1/*dummy utility*/,0,1,1,TRANSLATION);
						tie(e2,added) = my_add_edge(before[before.size()-1].second,spvertices[0].second,1/*dummy utility*/,0,1,1,TRANSLATION);	
					}
					else
					{
						tie(e1,added) = my_add_edge(before[before.size()-1].first,spvertices[0].second,1/*dummy utility*/,0,1,1,TRANSLATION);
						tie(e2,added) = my_add_edge(before[before.size()-1].second,spvertices[0].first,1/*dummy utility*/,0,1,1,TRANSLATION);
					}
				}
				for(int co1 = 0; co1 < spvertices.size()-1; co1++)
				{
					if(permutations[utilpos][permpc++] == 0)
					{
						tie(e1,added) = my_add_edge(spvertices[co1].first,spvertices[co1+1].first,1/*dummy utility*/,0,1,1,TRANSLATION);
						tie(e2,added) = my_add_edge(spvertices[co1].second,spvertices[co1+1].second,1/*dummy utility*/,0,1,1,TRANSLATION);	
					}
					else
					{
						tie(e1,added) = my_add_edge(spvertices[co1].first,spvertices[co1+1].second,1/*dummy utility*/,0,1,1,TRANSLATION);
						tie(e2,added) = my_add_edge(spvertices[co1].second,spvertices[co1+1].first,1/*dummy utility*/,0,1,1,TRANSLATION);	
					}
				}
				if(after.size()!=0)
				{
					if(permutations[utilpos][permpc++] == 0)
					{
						tie(e1,added) = my_add_edge(spvertices[spvertices.size()-1].first,after[0].first,1/*dummy utility*/,0,1,1,TRANSLATION);
						tie(e1,added) = my_add_edge(spvertices[spvertices.size()-1].second,after[0].second,1/*dummy utility*/,0,1,1,TRANSLATION);
					}
					else
					{
						tie(e1,added) = my_add_edge(spvertices[spvertices.size()-1].first,after[0].second,1/*dummy utility*/,0,1,1,TRANSLATION);
						tie(e1,added) = my_add_edge(spvertices[spvertices.size()-1].second,after[0].first,1/*dummy utility*/,0,1,1,TRANSLATION);
					}
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
	TGraph::vertex_iterator del,del_end,del_next,del_tmp;
	tie(del,del_end) = vertices(g);
	for(; del!=del_end; del = del_next)
	{
		del_next = del;
		++del_next;
		if(in_degree(*del,g) == 0 && out_degree(*del,g)==0)
			remove_vertex(*del,g);
		tie(del_tmp,del_end) = vertices(g);
	}
	/*

	int loopc = 0;
	for(tie(v_i,v_end) = vertices(g);v_i!=v_end; v_i = v_next)
	{
	tie(v_temp,v_end) = vertices(g);
	if(v_i==v_end)
	{
	break;
	}
	loopc++;
	printf("loopc = %d\r",loopc);
	v_next = v_i;
	_TRACE;
	++v_next;
	_TRACE;
	if(1)
	{
	if(in_degree(*v_i,g)==2)
	{
	if(out_degree(*v_i,g)==2)
	{

	//case 1

	//split the cell into two 
	//remove 4 edges 
	//remove the vertex 
	//add two new entries in limages[t]
	//add two new entries in rimages[t] TODO
	//compute features and add it to fvector[t]
	//add two new vertices to the graph
	//add 4 new edges

	TGraph::vertex_descriptor vd = *v_i;

	FeaturesType ft = fvector[g[vd].t][g[vd].findex];
	if(g[vd].t == 4 && ft.num == 3)
	{
	printf("FT : num = %d X = %0.0f Y = %0.0f Z = %0.0f\n",ft.num,ft.Centroid[0],ft.Centroid[1],ft.Centroid[2]);
	}
	std::vector<LabelImageType::Pointer> lout;
	std::vector<InputImageType::Pointer> rout;
	std::vector<FeaturesType> fvecout;
	_TRACE;
	SplitCell(limages[g[vd].t][g[vd].findex],rimages[g[vd].t][g[vd].findex],fvector[g[vd].t][g[vd].findex],fvar,lout,rout,fvecout);
	_TRACE;
	limages[g[vd].t].push_back(lout[0]);
	limages[g[vd].t].push_back(lout[1]);

	rimages[g[vd].t].push_back(rout[0]);
	rimages[g[vd].t].push_back(rout[1]);

	fvector[g[vd].t].push_back(fvecout[0]);
	fvector[g[vd].t].push_back(fvecout[1]);

	int cur_t = g[vd].t;

	TGraph::in_edge_iterator it,itend;
	tie(it,itend) = in_edges(*v_i,g);

	TGraph::edge_descriptor e1,e2,e3,e4;
	e1 = *it;++it;
	e2 = *it;
	TGraph::out_edge_iterator ot,otend;
	tie(ot,otend) = out_edges(*v_i,g);
	e3 = *ot;++ot;
	e4 = *ot;
	_TRACE;
	TGraph::vertex_descriptor vin1, vin2, vout1, vout2, vcur;
	vin1 = source(e1,g);
	vin2 = source(e2,g);
	vout1 = target(e3,g);
	vout2 = target(e4,g);

	printf("At t=%d vin1=%d vin2=%d vout1=%d vout2=%d v_i.findex = %d\n",cur_t,fvector[cur_t][g[vin1].findex].num,fvector[cur_t][g[vin2].findex].num,fvector[cur_t][g[vout1].findex].num,fvector[cur_t][g[vout2].findex].num,g[*v_i].findex);
	//scanf("%*d");
	vcur = *v_i;

	print_vertex(vin1,0);
	print_vertex(vin2,0);
	print_vertex(vout1,0);
	print_vertex(vout2,0);

	//remove_edge(e1,g);
	//remove_edge(e2,g);
	//remove_edge(e3,g);
	//remove_edge(e4,g);


	fvector[cur_t][g[vcur].findex].Centroid[0]+=50;

	TGraph::vertex_descriptor v1, v2;
	v1 = add_vertex(g);
	v2 = add_vertex(g);
	g[v1].findex = fvector[cur_t].size()-2;
	g[v2].findex = fvector[cur_t].size()-1;

	g[v1].special = 0;
	g[v2].special = 0;
	g[v1].t = cur_t;
	g[v2].t = cur_t;

	printf("checking v1, v2:\n");
	print_vertex(v1,1);
	print_vertex(v2,1);
	//scanf("%*d");
	float util[4];
	util[0] = compute_LRUtility(fvector[g[vin1].t][g[vin1].findex],fvector[g[v1].t][g[v1].findex],fvector[g[vout1].t][g[vout1].findex]) + compute_LRUtility(fvector[g[vin2].t][g[vin2].findex],fvector[g[v2].t][g[v2].findex],fvector[g[vout2].t][g[vout2].findex]);
	util[1] = compute_LRUtility(fvector[g[vin1].t][g[vin1].findex],fvector[g[v1].t][g[v1].findex],fvector[g[vout2].t][g[vout2].findex]) + compute_LRUtility(fvector[g[vin2].t][g[vin2].findex],fvector[g[v2].t][g[v2].findex],fvector[g[vout1].t][g[vout1].findex]);
	util[2] = compute_LRUtility(fvector[g[vin1].t][g[vin1].findex],fvector[g[v2].t][g[v2].findex],fvector[g[vout1].t][g[vout1].findex]) + compute_LRUtility(fvector[g[vin2].t][g[vin2].findex],fvector[g[v1].t][g[v1].findex],fvector[g[vout2].t][g[vout2].findex]);
	util[3] = compute_LRUtility(fvector[g[vin1].t][g[vin1].findex],fvector[g[v2].t][g[v2].findex],fvector[g[vout2].t][g[vout2].findex]) + compute_LRUtility(fvector[g[vin2].t][g[vin2].findex],fvector[g[v1].t][g[v1].findex],fvector[g[vout1].t][g[vout1].findex]);
	int max1 = -1;
	int maxpos1 = -1;
	for(int counter= 0; counter < 4; counter++)
	{
	if(max1 < util[counter])
	{
	max1 = util[counter];
	maxpos1 = counter;
	}
	}
	bool added = false;
	if(maxpos1 !=-1)
	{
	switch(maxpos1)
	{
	case 0:
	printf("Case 0\n");
	tie(e1,added) = my_add_edge(vin1,v1,util[0]/4,0,0,1,TRANSLATION);
	tie(e2,added) = my_add_edge(vin2,v2,util[0]/4,0,0,1,TRANSLATION);
	tie(e3,added) = my_add_edge(v1,vout1,util[0]/4,0,0,1,TRANSLATION);
	tie(e4,added) = my_add_edge(v2,vout2,util[0]/4,0,0,1,TRANSLATION);
	break;
	case 1:
	printf("Case 1\n");
	tie(e1,added) = my_add_edge(vin1,v1,util[1]/4,0,0,1,TRANSLATION);
	tie(e2,added) = my_add_edge(vin2,v2,util[1]/4,0,0,1,TRANSLATION);
	tie(e3,added) = my_add_edge(v1,vout2,util[1]/4,0,0,1,TRANSLATION);
	tie(e4,added) = my_add_edge(v2,vout1,util[1]/4,0,0,1,TRANSLATION);
	break;
	case 2:
	printf("Case 2\n");
	tie(e1,added) = my_add_edge(vin2,v1,util[2]/4,0,0,1,TRANSLATION);
	tie(e2,added) = my_add_edge(vin1,v2,util[2]/4,0,0,1,TRANSLATION);
	tie(e3,added) = my_add_edge(v1,vout2,util[2]/4,0,0,1,TRANSLATION);
	tie(e4,added) = my_add_edge(v2,vout1,util[2]/4,0,0,1,TRANSLATION);
	break;
	case 3:
	printf("Case 3\n");
	tie(e1,added) = my_add_edge(vin1,v2,util[3]/4,0,0,1,TRANSLATION);
	tie(e2,added) = my_add_edge(vin2,v1,util[3]/4,0,0,1,TRANSLATION);
	tie(e3,added) = my_add_edge(v1,vout1,util[3]/4,0,0,1,TRANSLATION);
	tie(e4,added) = my_add_edge(v2,vout2,util[3]/4,0,0,1,TRANSLATION);
	break;
	default:
	printf("Cannot come to this default case in the switch. Check\n");
	break;
	}
	}
	else
	{
	printf("?? %f %f %f %f\n",util[0],util[1],util[2],util[3]);
	//scanf("%*d");
	}
	//scanf("%*d");
	printf("checking 2:\n");
	print_vertex(v1,1);
	print_vertex(v2,1);
	printf("End checking2\n");
	print_vertex(vin1,0);
	print_vertex(vin2,0);
	print_vertex(vout1,0);
	print_vertex(vout2,0);

	clear_vertex(vcur,g);
	remove_vertex(vcur,g);

	//print_vertex(vin1,0);
	//print_vertex(vin2,0);
	//print_vertex(vout1,0);
	//print_vertex(vout2,0);

	//scanf("%*d");
	break;
	}
	else if(out_degree(*v_i,g)==1)
	{
	TGraph::out_edge_iterator et,etend;
	tie(et,etend) = out_edges(*v_i,g);
	TGraph::vertex_descriptor vnext;
	vnext = target(*et,g);
	if(out_degree(vnext,g)==2)
	{
	//case 2
	}
	}
	}
	}
	}
	*/
	for(tie(v_i,v_end) = vertices(g);v_i!=v_end; ++v_i)
	{
		if(in_degree(*v_i,g) >2 || out_degree(*v_i,g) >2)
		{
			printf("problem here3:");
			print_vertex(*v_i,0);
			//scanf("%*d");
		}
	}




	VectorPixelType col1,col2,col3,col4;
	col1[0] = 255;col1[1] = 0;col1[2] = 0;
	col2[0] = 255;col2[1] = 255;col2[2] = 0;
	col3[0] = 255;col3[1] = 255;col3[2] = 255;
	col4[0] = 255;col4[1] = 0; col4[2] = 255;

	for(tie(e_i,e_end) = edges(g); e_i!=e_end; ++e_i)
	{

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
				draw_line_for_edge(*e_i,col1,col2,-2);
			}
			else
			{
				draw_line_for_edge(*e_i,col3,col4,2);
			}
		}
	}

	std::vector<int> component(num_vertices(g));
	int num = connected_components(g,&component[0]);
	int tempnum = num;

	for(tie(e_i,e_end) = edges(g); e_i!= e_end ; e_i = e_next)
	{
		e_next = e_i;
		++e_next;
		if(g[*e_i].selected == 0)
		{

			remove_edge(*e_i,g);
		}
	}



	std::vector<int> vertex_count(num);
	std::vector<int> in_vertices_count(num);
	std::vector<int> out_vertices_count(num);
	for(int counter=0; counter < vertex_count.size(); counter++)
	{
		vertex_count[counter]=0;
		in_vertices_count[counter] = 0;
		out_vertices_count[counter] = 0;
	}


	boost::property_map<TGraph, boost::vertex_index_t>::type index;
	index = get(boost::vertex_index,g);



	for(tie(v_i,v_end) = vertices(g);v_i!=v_end; ++v_i)
	{
		vertex_count[component[index[*v_i]]]++;
		if(in_degree(*v_i,g)==0)
			in_vertices_count[component[index[*v_i]]]++;
		if(out_degree(*v_i,g)==0)
			out_vertices_count[component[index[*v_i]]]++;
	}
	for(int counter =0; counter < num; counter++)
	{
		if(vertex_count[counter]>1)
		{
			if(!(in_vertices_count[counter]==1 && out_vertices_count[counter] == 1))
			{
				if(in_vertices_count[counter]!=out_vertices_count[counter])
					printf("component+1 = %d in_count = %d out_count = %d\n",component[counter]+1,in_vertices_count[counter],out_vertices_count[counter]);
			}
		}
	}

	tie(v_i,v_end) = vertices(g);
	for(;v_i!=v_end; ++v_i)
	{
		if(vertex_count[component[index[*v_i]]]==1 && g[*v_i].special == 1)
		{
			tempnum--;
		}
		else
		{
			//printf("%d ",component[index[*v_i]]);
		}
	}
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
	PAUSE;
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
			printf("in\t");
			print_vertex(source(*e_i,g),depth-1);
		}
		TGraph::out_edge_iterator e_o,e_end2;
		for(tie(e_o,e_end2) = out_edges(v,g); e_o!=e_end2; ++e_o)
		{
			printf("out\t");
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

int main()//int argc, char **argv)
{
	//ST();
	// testKPLS();
	int num_tc = 30;
#define BASE "C:\\Users\\arun\\Research\\Tracking\\harvard\\cache\\testTSeries-02102009-1455-624"
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
			}
			input.push_back(cimp);
			fvector.push_back(f);
			limages.push_back(li);
			rimages.push_back(ri);
		}

		FeatureVariances fvar;
		for(int counter=0; counter < FeatureVariances::N; counter++)
		{
			fvar.variances[counter] = std::numeric_limits<float>::max();
		}
		InputImageType::SizeType imsize = tempimage->GetLargestPossibleRegion().GetSize();
		fvar.BoundingBox[0] = 0;
		fvar.BoundingBox[2] = 0;
		fvar.BoundingBox[4] = 0;
		fvar.BoundingBox[1] = imsize[0]-1;
		fvar.BoundingBox[3] = imsize[1]-1;
		fvar.BoundingBox[5] = imsize[2]-1;
		fvar.distVariance = 50;
		fvar.spacing[0] = 0.965;
		fvar.spacing[1] = 0.965;
		fvar.spacing[2] = 4.0;
		fvar.timeVariance = 0.1;
		fvar.overlapVariance = 1;
		fvar.variances[FeatureVariances::VOLUME] = 40000;

		CellTracker ct(fvector,fvar,limages,rimages);
		ColorImageType::Pointer debugcol = ColorImageType::New();
		ColorImageType::SizeType colsize;
		ColorImageType::RegionType colregion;
		ColorImageType::IndexType colindex;
		colindex.Fill(0);
		colsize[0] = imsize[0]*2;
		colsize[1] = imsize[1]*2;
		colsize[2] = num_t;
		colregion.SetIndex(colindex);
		colregion.SetSize(colsize);
		debugcol->SetRegions(colregion);
		debugcol->Allocate();
		ct.set_debug_image(debugcol);
		ct.run();
		//ct.writeGraphViz("C:\\Users\\arun\\Research\\Tracking\\harvard\\graphviz_testSeries.dot");
		//ct.writeGraphML("C:\\Users\\arun\\Research\\Tracking\\harvard\\graphviz_testSeries.graphml");
		//PAUSE;
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
		ColorImageType::Pointer colin = getColorImageFromColor2DImages(input);
		ColorImageType::Pointer colout = getColorImageFromColor2DImages(output);
		writeImage<ColorImageType>(debugcol,"C:\\Users\\arun\\Research\\Tracking\\harvard\\debug\\debugcol.tif");
		writeImage<ColorImageType>(colin,"C:\\Users\\arun\\Research\\Tracking\\harvard\\debug\\input.tif");
		writeImage<ColorImageType>(colout,"C:\\Users\\arun\\Research\\Tracking\\harvard\\debug\\output.tif");


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
