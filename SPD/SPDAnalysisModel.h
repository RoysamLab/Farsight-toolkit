#ifndef SPDANALYSISMODEL_H
#define SPDANALYSISMODEL_H
#include <vector>
#include <vtkTable.h>
#include <vtkSmartPointer.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <map>
#include <QString>
#include "ClusClus/clusclus.h"

typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS, 
	boost::property< boost::vertex_distance_t, unsigned int>, boost::property< boost::edge_weight_t, double> > Graph;
typedef std::pair < unsigned int, unsigned int>Edge;

typedef struct Tree
{
	Tree()
	{
	};
	Tree(int fir, int sec, double c, int par)
	{
		first = fir;
		second = sec;
		cor = c;
		parent = par;
	};
	int first;
	int second;
	double cor;
	int parent;
}Tree;


enum DISTANCE_TYPE
{
	CITY_BLOCK,
	SQUARE_DISTANCE
};

class SPDAnalysisModel
{
public:
	static SPDAnalysisModel* InitInstance();
	static void DeInstance();

	void GetTableHeaders(std::vector<std::string> &headers);
	vtkSmartPointer<vtkTable> GetDataTable();
	QString GetFileName();

	bool ReadCellTraceFile(std::string fileName, bool btest);
	void ParseTraceFile(vtkSmartPointer<vtkTable> table);

	unsigned int GetSampleNum();
	unsigned int GetFeatureNum();

	void NormalizeData();
	int ClusterAgglomerate( double cor, double mer);
	void ClusterMerge( double cor, double mer);
	void HierachicalClustering(vtkSmartPointer<vtkTable> table, bool bcol = true);
	void HierachicalClustering();
	void ClusterCells(double cor);
	void GenerateMST();
	void GenerateDistanceMST();
	vtkSmartPointer<vtkTable> GetMSTTable( int MSTIndex);
	void RunEMDAnalysis();
	void GetEMDMatrixDivByMax(vnl_matrix<double> &emdMatrix);
	void GetClusClusData(clusclus& c1, clusclus& c2, double threshold, double magFactor = 0);
	vtkSmartPointer<vtkTable> GenerateProgressionTree( std::string& selectedModules);
	void GetSelectedFeatures(std::set<long int>& selectedFeatures);
	void SaveSelectedFeatureNames(QString filename, std::set<long int>& selectedFeatures);
	double GetEMDSelectedPercentage(double thres);
	double GetEMDSelectedThreshold( double per);
	void GetMatrixData(vnl_matrix<double> &mat);
	void GetClusterMapping( std::vector<int> &indToclusInd);
	void SetProgressionTag(bool bProg);
	bool GetProgressionTag();
	void GetDistanceOrder(std::vector<long int> &order);

protected:
	SPDAnalysisModel();
	~SPDAnalysisModel();
	void split( std::string& s, char delim,std::vector< std::string >* ret);
	int LineNum( const char* fileName);
	int ClusterAggFeatures( vnl_vector<unsigned int>& index, vnl_matrix<double>& mean, double cor, vnl_vector<int>& TreeIndex, int fold);
	vnl_vector<int> GetModuleSize( vnl_vector<unsigned int>& index);
	void GetCombinedMatrix( vnl_matrix<double> &datamat, vnl_vector<unsigned int>& index, unsigned int moduleId, unsigned int moduleDeleteId, vnl_matrix<double>& mat);
	void GetCombinedMatrix( vnl_matrix<double> &datamat, vnl_vector< unsigned int>& index, std::vector< unsigned int> moduleID, vnl_matrix<double>& mat);
	void GetMatrixRowMeanStd(vnl_matrix<double>& mat, vnl_vector<double>& mean, vnl_vector<double>& std);
	void StandardizeIndex(vnl_vector<unsigned int>& index);
	void EraseZeroCol(vnl_matrix<double>& mat);
	void EraseCols(vnl_matrix<double>& mat, std::vector<unsigned int> vec);
	void SubstitudeVectorElement( vnl_vector<unsigned int>& vector, unsigned int ori, unsigned int newValue);
	void DeleteMatrixColumn( vnl_matrix<double>& mat, unsigned int col);
	double CityBlockDist( vnl_matrix<double>& mat, unsigned int ind1, unsigned int ind2);
	double EuclideanBlockDist( vnl_matrix<double>& mat, unsigned int ind1, unsigned int ind2);
	int GetSingleModuleSize( vnl_vector<unsigned int>& index, unsigned int ind);
	void GetMatrixDistance( vnl_matrix<double>& data, vnl_vector<double>&distance, DISTANCE_TYPE type);
	void GetMSTMatrixDistance( vnl_vector<double>& distance, std::vector< boost::graph_traits<Graph>::vertex_descriptor>& vertex, vnl_vector<double>& MSTdistance);
	void Hist(vnl_vector<double>&distance, int num_bin, vnl_vector<double>& interval, vnl_vector<unsigned int>& histDis);
	void Hist(vnl_vector<double>&distance, vnl_vector<double>& interval, vnl_vector<unsigned int>& histDis);
	//double Dist(int *first, int *second);
	double EarthMoverDistance(vnl_vector<unsigned int>& first, vnl_vector<unsigned int>& second, vnl_matrix<double> &flowMatrix);
	bool IsExist(std::vector<unsigned int> vec, unsigned int value);
	long int GetSelectedFeatures(vnl_vector< unsigned int> & index, std::vector<unsigned int> moduleID, std::set<long int>& featureSelectedIDs);
	double VnlVecMultiply(vnl_vector<double> const &vec1, vnl_vector<double> const &vec2);

public:
	std::vector< Tree> PublicTreeData;

private:
	static SPDAnalysisModel *s_pmodel;
	std::vector< Tree> TreeData;
	vnl_vector<int> TreeIndexNew;
	//save filename
	QString filename;

	// basic data storation
	std::vector<int> CellTraceIndex;
	std::vector<std::string> FeatureNames;
	std::vector<double> DistanceToDevice;
	vnl_matrix<double> DataMatrix;			// normalized data feature for analysis
	vtkSmartPointer<vtkTable> DataTable;
	std::vector<long int> indMapFromIndToVertex;    // index mapping
	std::vector<std::string> headers;
	std::vector<int> FeatureIndex;
	
	// data for cell agglomeration
	vnl_vector<unsigned int> CellClusterIndex;
	std::vector< std::vector<int> > CellCluster;
	vnl_matrix<double> MatrixAfterCellCluster;
		
	//data for feature agglormeration
	vnl_vector<unsigned int> ClusterIndex;
	vnl_matrix<double> ModuleMean;
	vnl_vector<int> TreeIndex;     // tree node index

	//MST for each module
	bool bProgression;

	std::vector< std::vector< boost::graph_traits<Graph>::vertex_descriptor> > ModuleMST;
	std::vector< boost::graph_traits<Graph>::vertex_descriptor> DistanceMST;
	std::vector< double> MSTWeight;
	std::vector< std::vector< boost::graph_traits<Graph>::vertex_descriptor> > ModuleGraph;
	std::vector< vtkSmartPointer<vtkTable> > MSTTable;    // data to pass to the views
	std::vector< vtkSmartPointer<vtkTable> > MSTDistanceTable;    // data to pass to the views

	//data for EMD 
	vnl_matrix<double> EMDMatrix;
	std::set< long int> selectedFeatureIDs;
	vnl_vector<double> DistanceEMDVector;

	// for heatmap
	vnl_matrix<double> heatmapMatrix;
	vnl_matrix<double> heatmapMatrixNew;
	clusclus *cc1;
	clusclus *cc2;
};
#endif
