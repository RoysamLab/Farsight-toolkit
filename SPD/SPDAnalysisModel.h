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
	void GenerateMST();
	vtkSmartPointer<vtkTable> GetMSTTable( int MSTIndex);
	void RunEMDAnalysis();

	void GetClusClusData(clusclus& c1, clusclus& c2, double threshold);
	vtkSmartPointer<vtkTable> GenerateProgressionTree( std::string& selectedModules);

protected:
	SPDAnalysisModel();
	~SPDAnalysisModel();
	void split( std::string& s, char delim,std::vector< std::string >* ret);
	int LineNum( const char* fileName);
	int ClusterAggFeatures( vnl_vector<unsigned int>& index, vnl_matrix<double>& mean, double cor, int fold);
	vnl_vector<int> GetModuleSize( vnl_vector<unsigned int>& index);
	void GetCombinedMatrix( vnl_vector<unsigned int>& index, unsigned int moduleId, unsigned int moduleDeleteId, vnl_matrix<double>& mat);
	void GetCombinedMatrix( vnl_vector< unsigned int>& index, std::vector< unsigned int> moduleID, vnl_matrix<double>& mat);
	void GetMatrixRowMeanStd(vnl_matrix<double>& mat, vnl_vector<double>& mean, vnl_vector<double>& std);
	void StandardizeIndex(vnl_vector<unsigned int>& index);
	void EraseZeroCol(vnl_matrix<double>& mat);
	void SubstitudeVectorElement( vnl_vector<unsigned int>& vector, unsigned int ori, unsigned int newValue);
	void DeleteMatrixColumn( vnl_matrix<double>& mat, unsigned int col);
	double CityBlockDist( vnl_matrix<double>& mat, unsigned int ind1, unsigned int ind2);
	double EuclideanBlockDist( vnl_matrix<double>& mat, unsigned int ind1, unsigned int ind2);
	int GetSingleModuleSize( vnl_vector<unsigned int>& index, unsigned int ind);
	void GetMatrixDistance( vnl_matrix<double>& data, vnl_vector<double>&distance, DISTANCE_TYPE type);
	void GetMSTMatrixDistance( vnl_vector<double>& distance, std::vector< boost::graph_traits<Graph>::vertex_descriptor>& vertex, vnl_vector<double>& MSTdistance);
	vnl_vector<unsigned int> Hist(vnl_vector<double>&distance, int num_bin, vnl_vector<double>& interval);
	vnl_vector<unsigned int> Hist(vnl_vector<double>&distance, vnl_vector<double>& interval);
	//double Dist(int *first, int *second);
	double EarthMoverDistance(vnl_vector<unsigned int>& first, vnl_vector<unsigned int>& second, vnl_matrix<double> &flowMatrix);
	bool IsExist(std::vector<unsigned int> vec, unsigned int value);

private:
	static SPDAnalysisModel *s_pmodel;

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

	//data for agglormeration
	vnl_vector<unsigned int> ClusterIndex;
	vnl_matrix<double> ModuleMean;
	vnl_vector<int> ModuleSize;

	//MST for each module
	std::vector< std::vector< boost::graph_traits<Graph>::vertex_descriptor> > ModuleMST;
	std::vector< double> MSTWeight;
	std::vector< std::vector< boost::graph_traits<Graph>::vertex_descriptor> > ModuleGraph;
	std::vector< vtkSmartPointer<vtkTable> > MSTTable;    // data to pass to the views

	//data for EMD 
	vnl_matrix<double> EMDMatrix;

	// for heatmap
	vnl_matrix<double> heatmapMatrix;
	clusclus *cc1;
	clusclus *cc2;
};
#endif
