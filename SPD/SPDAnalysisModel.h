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
	SPDAnalysisModel();
	~SPDAnalysisModel();

	void GetTableHeaders(std::vector<std::string> &headers);
	vtkSmartPointer<vtkTable> GetDataTable();
	QString GetFileName();

	bool ReadCellTraceFile(std::string fileName, bool bContrast);
	void ParseTraceFile(vtkSmartPointer<vtkTable> table, bool bContrast = false);
	bool ReadRawData(std::string fileName);
	unsigned int GetSampleNum();
	unsigned int GetFeatureNum();
	unsigned int GetContrastDataSampleNum();

	int ClusterAgglomerate( double cor, double mer);
	void ClusterMerge( double cor, double mer);
	void HierachicalClustering();
	void ClusterSamples(double cor);
	int ClusterSamples( double cor, vnl_matrix<double> &mat, vnl_vector<unsigned int> &index);
								
	void GetCellClusterSize( std::vector<int> &clusterSize);
	vtkSmartPointer<vtkTable> GetDataTableAfterCellCluster();
	void GetFeatureIdbyModId(std::vector<unsigned int> &modID, std::vector<unsigned int> &featureID);

	void GetSingleLinkageClusterAverage(std::vector< std::vector< long int> > &index, vnl_matrix<double> &clusAverageMat);  // after single linkage clustering
	void GetClusterMapping( std::map< int, int> &index);
	void GetSingleLinkageClusterMapping(std::vector< std::vector< long int> > &index, std::vector<int> &newIndex);   

	void GenerateMST();
	vtkSmartPointer<vtkTable> GenerateMST( vnl_matrix<double> &mat, std::vector< unsigned int> &selFeatures, std::vector<int> &clusterNum);
    vtkSmartPointer<vtkTable> GenerateSubGraph( vnl_matrix<double> &mat, std::vector< std::vector< long int> > &clusIndex, std::vector< unsigned int> &selFeatures, std::vector<int> &clusterNum);
	void GenerateDistanceMST();
	vtkSmartPointer<vtkTable> GetMSTTable( int MSTIndex);
	void RunEMDAnalysis();
	void GetEMDMatrixDivByMax(vnl_matrix<double> &emdMatrix);
	void GetClusClusData(clusclus* c1, double threshold, std::vector< unsigned int> *disModIndex = NULL);
	vtkSmartPointer<vtkTable> GenerateProgressionTree( std::string& selectedModules);
	void GetSelectedFeatures(std::set<long int>& selectedFeatures);
	void SaveSelectedFeatureNames(QString filename, std::vector<int>& selectedFeatures);
	void SaveSelectedFeatureNames(QString filename, std::vector<unsigned int>& selectedFeatures);
	double GetEMDSelectedPercentage(double thres);
	double GetEMDSelectedThreshold( double per);
	void GetMatrixData(vnl_matrix<double> &mat);
	void GetClusterMatrixData(vnl_matrix<double> &mat);
	bool SetProgressionType(bool bProg);
	bool GetProgressionType();
	void GetDistanceOrder(std::vector<long int> &order);
	void GetClusterOrder(std::vector< std::vector< long int> > &clusIndex, std::vector<long int> &treeOrder, std::vector< int> &clusterOrder);

	void ModuleCoherenceMatchAnalysis();
	void ModuleCorrelationMatrixMatch(unsigned int kNeighbor);

	void GetClusClusDataForCorMatrix( clusclus* c1, clusclus* c2, double threshold, std::vector< unsigned int> *disModIndex = NULL);
	double GetCorMatSelectedPercentage(double thres);
	void GetCombinedDataTable(vtkSmartPointer<vtkTable> table);
	void SetMaxVertexID(int verId);
	void GetPercentage(std::vector< std::vector< long int> > &clusIndex, std::vector< double> &colorVec);
	void GetCloseToDevicePercentage( std::vector< std::vector< long int> > &clusIndex, std::vector< double> &disPer, double disThreshold);
	void GetClusterFeatureValue(std::vector< std::vector< long int> > &clusIndex, int nfeature, vnl_vector<double> &featureValue, std::string &featureName);
	vtkSmartPointer<vtkTable> GetAverModuleTable(std::vector< std::vector< long int> > &clusIndex, std::vector<long int> &TreeOrder, std::vector< double> &percentageOfSamples,
				    std::vector< double> &percentageOfNearDeviceSamples, std::vector< int> &selFeatureOrder, std::vector< int> &unselFeatureOrder, int maxId);
	void ConvertTableToMatrix(vtkSmartPointer<vtkTable> table, vnl_matrix<double> &mat, std::vector<int> &index, vnl_vector<double> &distance);
	void ConvertTableToMatrixForLayerData(vtkSmartPointer<vtkTable> table, vnl_matrix<double> &mat, std::vector<int> &index, vnl_vector<int> &clusNo);
	vtkSmartPointer<vtkTable> GetTableForHist(std::vector< int> &selFeatureOrder, std::vector< int> &unselFeatureOrder);
	int GetConnectedComponent(std::vector< unsigned int> &selFeatureID, std::vector<int> &component);
	unsigned int GetKNeighborNum();

	void BuildMSTForConnectedComponent(std::vector< unsigned int> &selFeatureID, std::vector<int> &component, int connectedNum);
	void GetComponentMinDistance(std::vector< unsigned int> selFeatureID, std::vector<int> &component, int connectedNum, vnl_matrix<double> &dis);
    bool SearchSubsetsOfFeatures(std::vector< unsigned int> &selModules);

	static double CaculatePS(unsigned int kNeighbor, unsigned int nbins, unsigned int bdevide, vnl_vector<double> vec1, vnl_vector<double> vec2);

protected:

	static void NormalizeData(vnl_matrix<double> &mat);
	void split( std::string& s, char delim,std::vector< std::string >* ret);
	int LineNum( const char* fileName);
	bool MergeMatrix( vnl_matrix<double> &firstMat, vnl_matrix<double> &secondMat, vnl_matrix<double> &combinedMat);
	void ConvertMatrixToTable(vtkSmartPointer<vtkTable> table, vnl_matrix<double> &mat, vnl_vector<double> &distance);
	void GetClusterIndexFromVnlVector(std::vector< std::vector<int> > &clusterIndex, vnl_vector<unsigned int> &index);
	void GetAverageModule( vnl_matrix<double> &mat, vnl_vector<double> &distance, std::vector< std::vector<int> > &index, vnl_matrix<double> &averageMat, vnl_vector<double> &averageDistance);
	int ClusterAggFeatures( vnl_matrix<double>& mainmatrix, vnl_vector<unsigned int>& index, vnl_matrix<double>& mean, double cor);
	vnl_vector<int> GetModuleSize( vnl_vector<unsigned int>& index);
	void GetCombinedMatrix( vnl_matrix<double> &datamat, vnl_vector<unsigned int>& index, unsigned int moduleId, unsigned int moduleDeleteId, vnl_matrix<double>& mat);
	void GetCombinedMatrix( vnl_matrix<double> &datamat, vnl_vector< unsigned int>& index, std::vector< unsigned int> moduleID, vnl_matrix<double>& mat);
	void GetCombinedMatrix( vnl_matrix<double> &datamat, int nstart, int nrow, std::vector< unsigned int> selFeatureIDs, vnl_matrix<double>& mat);
	void GetCombinedMatrixByModuleId( vnl_matrix<double> &datamat, std::vector< std::vector< unsigned int> > &featureClusterIndex, std::vector< unsigned int> &selModuleIDs, vnl_matrix<double>& mat);
	void GetCombinedInversedMatrix(vnl_matrix<double> &datamat, vnl_vector<unsigned int>& index, unsigned int moduleId, unsigned int moduleInversedId, vnl_matrix<double>& mat);
	static void GetMatrixRowMeanStd(vnl_matrix<double>& mat, vnl_vector<double>& mean, vnl_vector<double>& std);
	void StandardizeIndex(vnl_vector<unsigned int>& index);
	void EraseZeroCol(vnl_matrix<double>& mat);
	void EraseCols(vnl_matrix<double>& mat, std::vector<unsigned int> vec);
	void SubstitudeVectorElement( vnl_vector<unsigned int>& vector, unsigned int ori, unsigned int newValue);
	void DeleteMatrixColumn( vnl_matrix<double>& mat, unsigned int col);
	double CityBlockDist( vnl_matrix<double>& mat, unsigned int ind1, unsigned int ind2);
	void CityBlockDist( vnl_matrix<double>& mat, vnl_matrix<double>& matDis);
	static void EuclideanBlockDist( vnl_matrix<double>& mat, vnl_matrix<double>& matDis);
    static void EuclideanBlockDist( vnl_vector<double>& vec, vnl_matrix<double>& matDis);
	double EuclideanBlockDist( vnl_vector<double>& vec1, vnl_vector<double>& vec2);
	static void FindNearestKSample(vnl_matrix< double> &modDist, std::vector< unsigned int>& index, unsigned int kNeighbor);
	static void FindFarthestKSample(vnl_matrix< double> &modDist, std::vector< unsigned int>& index, unsigned int kNeighbor);
	static void GetKWeights(vnl_matrix< double> &modDist, std::vector< unsigned int>& index, vnl_vector<double>& weights, unsigned int kNeighbor);
	int GetConnectedComponent(std::vector< unsigned int>& index, std::vector< int>& component, unsigned int kNeighbor);
	static double EuclideanBlockDist( vnl_matrix<double>& mat, unsigned int ind1, unsigned int ind2);
	int GetSingleModuleSize( vnl_vector<unsigned int>& index, unsigned int ind);
	void GetMatrixDistance( vnl_matrix<double>& data, vnl_vector<double>&distance, DISTANCE_TYPE type);
	void GetMSTMatrixDistance( vnl_vector<double>& distance, std::vector< boost::graph_traits<Graph>::vertex_descriptor>& vertex, vnl_vector<double>& MSTdistance);
	void Hist(vnl_vector<double>&distance, int num_bin, vnl_vector<double>& interval, vnl_vector<unsigned int>& histDis);
	void Hist(vnl_vector<double>&distance, vnl_vector<double>& interval, vnl_vector<unsigned int>& histDis);
	static void Hist(vnl_vector<double>&distance, double interVal, double minVal, vnl_vector<unsigned int>& histDis, int num_bin);
	static void Hist(vnl_matrix<double>&distance, double interVal, double minVal, vnl_vector<unsigned int>& histDis, int num_bin);
	double FindMinBetweenConnectComponent(vnl_matrix<double> &dis, std::vector<int> ver1, std::vector<int> ver2);
	void GetComponentMinDistance(vnl_matrix<double> &distMat, std::vector<int> &component, int connectedNum, vnl_matrix<double> &dis);

	//double Dist(int *first, int *second);
	static double EarthMoverDistance(vnl_vector<unsigned int>& first, vnl_vector<unsigned int>& second, vnl_matrix<double> &flowMatrix, int num_bin);
	bool IsExist(std::vector<unsigned int> vec, unsigned int value);
	long int GetSelectedFeatures(vnl_vector< unsigned int> & index, std::vector<unsigned int> moduleID, std::set<long int>& featureSelectedIDs);
	double VnlVecMultiply(vnl_vector<double> const &vec1, vnl_vector<double> const &vec2);
	bool MergeTables(vtkSmartPointer<vtkTable> firstTable, vtkSmartPointer<vtkTable> secondTable, vtkSmartPointer<vtkTable> table);
	void CopyTable( vtkSmartPointer<vtkTable> oriTable, vtkSmartPointer<vtkTable> targetTable);

	double ComputeModuleDistanceAndConnection(vnl_matrix<double> &mati, vnl_matrix<double> &matj, int &rowi, int &rowj);
	void GetAverageVec(vnl_matrix<double> &mat, vnl_vector<double> &vec);
        bool IsConnected(std::multimap<int, int> &neighborGraph, std::vector<long int> &index1, std::vector<long int> &index2);

	/// for multi-level demo
	bool RunSPDforFeatureDistributionTable(std::string fileName);
	void ComputeDistributionDistance( vnl_matrix<unsigned int> &mat, vnl_matrix<double> &dismat);
	void ComputeDistributionDistance(vnl_matrix<unsigned int> &mat, vnl_vector<double> &moduleDistance);
	bool GenerateMST( vnl_matrix<double> &mat, bool bfirst);
	void RunEMDAnalysis( vnl_vector<double> &moduleDistance, int ind);

public:
	std::vector< Tree> PublicTreeData;
	std::vector< Tree> mstTreeList;

private:
	std::vector< Tree> TreeData;
	vnl_vector<int> TreeIndexNew;
	//save filename
	QString filename;

	// basic data storation, main data for analysis
	std::vector<int> indMapFromIndToVertex;    // index mapping
	int maxVertexId;
	std::map< int, int> indMapFromVertexToClus;
	vnl_vector<double> DistanceToDevice;
	vnl_vector<double> UNDistanceToDevice;

	double disCor;
	vnl_matrix<double> DataMatrix;			// normalized data feature for analysis
	vnl_matrix<double> DataMatrixAfterCellCluster;
	vtkSmartPointer<vtkTable> DataTable;
	std::vector<std::string> headers;
	std::vector<int> FeatureIndex;

	/// for combined data 
	vnl_matrix<double> MatrixAfterCellCluster;     
	vnl_matrix<double> UNMatrixAfterCellCluster; //matAfterCellCluster;

	// for constrast data
	vnl_matrix<double> ContrastDataMatrix;			// normalized data feature for analysis
	vnl_matrix<double> ContrastMatrixAfterCellCluster;
	vtkSmartPointer<vtkTable> ContrastDataTable;
	std::vector<int> indMapFromIndToContrastVertex;
	vnl_vector<double> contrastDistance;
	vnl_vector<unsigned int> constrastCellClusterIndex;
	std::vector< std::vector<int> > contrastCellCluster;
	std::map< int, int> combinedIndexMapping;
	vnl_vector<unsigned int> combinedCellClusterIndex;
	
	// data for cell agglomeration
	vtkSmartPointer<vtkTable> DataTableAfterCellCluster;  // average data after cell cluster without normalization
	
	vnl_vector<unsigned int> CellClusterIndex;
	std::vector< std::vector<int> > CellCluster;


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
	unsigned int m_kNeighbor;

	// for heatmap
	vnl_matrix<double> heatmapMatrix;
	vnl_matrix<double> heatmapMatrixNew;
	std::vector< int> clusterOrder;
	std::vector< int> moduleForSelection;

	// for spdtestwindow
	vnl_matrix<double> ModuleCompareCorMatrix;
	vnl_vector<double> DistanceCorVector;

	// for multi-level demo
	int nSampleSize;
	int nFeatureSize;
	int nBinNum;

	// for layer data
	vnl_vector<int> clusNo;
};
#endif
