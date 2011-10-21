#include "SPDAnalysisModel.h"
#include <fstream>
#include <string>
#include <vtkVariant.h>
#include <math.h>
#include <mbl/mbl_stats_nd.h>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include "ftkUtils.h"
#include "transportSimplex.h"

#define NUM_BIN 20

SPDAnalysisModel* SPDAnalysisModel::s_pmodel = NULL;

SPDAnalysisModel* SPDAnalysisModel::InitInstance()
{
	if( s_pmodel == NULL)
	{
		s_pmodel = new SPDAnalysisModel;

	}
	return s_pmodel;
}

void SPDAnalysisModel::DeInstance()
{
	if( s_pmodel != NULL)
	{
		delete s_pmodel;
	}
}

SPDAnalysisModel::SPDAnalysisModel()
{
	DataTable = vtkSmartPointer<vtkTable>::New();
	headers.push_back("node1");
	headers.push_back("node2");
	headers.push_back("weight");
}

SPDAnalysisModel::~SPDAnalysisModel()
{

}

vtkSmartPointer<vtkTable> SPDAnalysisModel::GetDataMatrix()
{
	return this->DataTable;
}

bool SPDAnalysisModel::ReadCellTraceFile(std::string fileName)
{
	this->DataTable = ftk::LoadTable(fileName);
	if( this->DataTable != NULL)
	{
		ParseTraceFile(this->DataTable);
		return true;
	}
	else
	{
		return false;
	}
}
	//int line = LineNum(fileName);

	//std::ifstream file(fileName, std::ifstream::in);
	//if( !file.is_open() || line == -1)
	//{
	//	return false;
	//}

	//std::string feature;

	//bool bfirst = true;
	//int rowIndex = 0; 
	//std::vector<std::string> rowValue;
	//while( !file.eof())
	//{
	//	getline(file, feature);
	//	if( feature.length() > 5)
	//	{
	//		split(feature, '\t', &rowValue);
	//		if( bfirst)
	//		{
	//			this->FeatureNames = rowValue;
	//			this->FeatureNames.erase(this->FeatureNames.begin());     // first trace index
	//			(this->FeatureNames).pop_back();
	//			(this->FeatureNames).pop_back();
	//			(this->FeatureNames).pop_back();                            // distance to device && cell names
	//			 bfirst = false;
	//			(this->DataMatrix).set_size(line - 2, this->FeatureNames.size());
	//		}
	//		else
	//		{
	//			std::vector<std::string>::iterator iter = rowValue.begin();
	//			for( int colIndex = 0, col = 0; col < rowValue.size() - 1; col++, iter++)
	//			{
	//				if ( col == 0)
	//				{
	//					this->CellTraceIndex.push_back( atoi((*iter).c_str()));    // trace index
	//				}
	//				else if( col == rowValue.size() - 3)		// strings to be eliminated
	//				{
	//					//break;
	//				}
	//				else if( col == rowValue.size() - 2)	
	//				{
	//					this->DistanceToDevice.push_back( atof((*iter).c_str()));  // save distance to device in the array
	//				}
	//				else if( col > 0 && col < rowValue.size() - 3)
	//				{
	//					(this->DataMatrix)(rowIndex, colIndex++) = atof( (*iter).c_str());   // save the feature data in the matrix for analysis
	//				}
	//			}
	//			rowIndex++;
	//		}
	//		rowValue.clear();
	//	}
	//	feature.clear();
	//}
	////std::ofstream ofs("readdata.txt", std::ofstream::out);
	////ofs<<this->DataMatrix<<endl;


void SPDAnalysisModel::ParseTraceFile(vtkSmartPointer<vtkTable> table)
{
	for( long int i = 0; i < table->GetNumberOfRows(); i++)
	{
		long int var = table->GetValue( i, 0).ToLong();
		this->indMapFromIndToVertex.push_back( var);
	}

	this->DataMatrix.set_size( table->GetNumberOfRows(), table->GetNumberOfColumns() - 2);

	for( int i = 0, rowIndex = 0; i < table->GetNumberOfRows(); i++, rowIndex++)
	{
		int colIndex = 0;
		for( int j = 0; j < table->GetNumberOfColumns(); j++)
		{
			if( j == 0)
			{
				this->CellTraceIndex.push_back( table->GetValue(i, j).ToInt());
			}
			else if( j == table->GetNumberOfColumns() - 1)
			{
				this->DistanceToDevice.push_back( table->GetValue(i, j).ToDouble());
			}
			else
			{
				(this->DataMatrix)(rowIndex, colIndex++) = table->GetValue(i, j).ToDouble();
			}
		}
	}
}

int SPDAnalysisModel::LineNum(const char* fileName)
{
	std::ifstream file(fileName, std::ifstream::in);

	if( !file.is_open())
	{
		return -1;
	}

	std::string temp;
	int line = 0;
	while( getline(file, temp) && !file.eof())
	{
		line++;
	}
	file.close();
	return line;
}

unsigned int SPDAnalysisModel::GetSampleNum()
{
	return this->DataMatrix.rows();
}

unsigned int SPDAnalysisModel::GetFeatureNum()
{
	return this->DataMatrix.cols();
}

void SPDAnalysisModel::split(std::string& s, char delim, std::vector< std::string >* ret)
{
	size_t last = 0;
	size_t index = s.find_first_of(delim,last);
	while( index!=std::string::npos)
	{
		ret->push_back(s.substr(last,index-last));
		last = index+1;
		index=s.find_first_of(delim,last);
	}
	if( index-last>0)
	{
		ret->push_back(s.substr(last,index-last));
	}
} 

void SPDAnalysisModel::NormalizeData()
{
	vnl_vector<double> std_vec;
	vnl_vector<double> mean_vec;

	GetMatrixRowMeanStd(this->DataMatrix, mean_vec, std_vec);
	
	for( unsigned int i = 0; i < this->DataMatrix.columns(); ++i)
	{
		vnl_vector<double> temp_col = this->DataMatrix.get_column(i);
		if( std_vec(i) > 0)
		{	
			for( unsigned int j =0; j < temp_col.size(); ++j)
			{
				temp_col[j] = (temp_col[j] - mean_vec(i))/std_vec(i);
			}
			//temp_col = temp_col.normalize();
		}
		this->DataMatrix.set_column(i,temp_col);
	}
}

int SPDAnalysisModel::ClusterAgglomerate(double cor)
{
	vnl_vector<unsigned int> clusterIndex;
	clusterIndex.set_size( this->DataMatrix.cols());

	vnl_matrix<double> moduleMean = this->DataMatrix;
	moduleMean.normalize_columns();

	std::ofstream ofs("result.txt", std::ofstream::out);
	ofs<<"The clustering process:"<<endl;

	for( unsigned int i = 0; i < this->DataMatrix.cols(); i++)
	{
		clusterIndex[i] = i;
	}

	int old_cluster_num = this->DataMatrix.cols();
	int new_cluster_num = this->DataMatrix.cols();

	do
	{
		old_cluster_num = new_cluster_num;
		new_cluster_num = ClusterAggFeatures( clusterIndex, moduleMean, cor, 2);
		//ofs<<clusterIndex<<endl<<endl;
	}
	while( old_cluster_num != new_cluster_num);

	this->ClusterIndex = clusterIndex;
	this->ModuleMean = moduleMean;
	this->ModuleSize = GetModuleSize( clusterIndex);

	for( unsigned int i = 0; i <= this->ClusterIndex.max_value(); i++)
	{
		ofs<<"Cluster:"<<i + 1<<endl;
		for( unsigned int j = 0; j < this->ClusterIndex.size(); j++)
		{
			if( ClusterIndex[j] == i)
			{
				ofs<<j + 1<<endl;
			}
		}
	}

	ofs.close();

	std::cout<<"The cluster size after cluster: "<<this->ClusterIndex.max_value() + 1<<endl;
	std::cout<<clusterIndex<<endl;

	return new_cluster_num;
}


int SPDAnalysisModel::ClusterAggFeatures(vnl_vector<unsigned int>& index, vnl_matrix<double>& mean, double cor, int fold = 2)
{
	vnl_vector<int> moduleSize = GetModuleSize(index);
	vnl_vector<int> isActiveModule;
	vnl_vector<double> moduleCenter;
	vnl_vector<double> moduleCor;   
	isActiveModule.set_size( moduleSize.size());
	unsigned int i = 0;
	for(  i = 0; i < moduleSize.size(); i++)
	{
		isActiveModule[i] = 1;
	}
	
	vnl_vector<double> zeroCol;
	zeroCol.set_size( mean.rows());
	for( i = 0; i < mean.rows(); i++)
	{
		zeroCol[i] = 0;
	}

	while( isActiveModule.sum() > 1)
	{
		vnl_vector<int> tmp;
		tmp.set_size( moduleSize.size());

		vnl_matrix<double> newModule;
		vnl_vector<double> newModuleMeans;
		vnl_vector<double> newModuleStd;
		vnl_vector<double> newModuleCor;

		for( i = 0; i < moduleSize.size(); i++)
		{
			if( moduleSize[i] == 0 || isActiveModule[i] == 0)
			{
				tmp[i] = moduleSize.max_value() + 1;
			}
			else
			{
				tmp[i] = moduleSize[i];
			}
		}

		unsigned int moduleId = tmp.arg_min();
		moduleCenter = mean.get_column(moduleId);
		moduleCor = moduleCenter * mean;
		moduleCor[moduleId] = 0;

		unsigned int moduleToDeleteId = moduleCor.arg_max();
		int newModuleSize = moduleSize[moduleId] + moduleSize[moduleToDeleteId];

		newModule.set_size(this->DataMatrix.rows(), newModuleSize);

		if( moduleToDeleteId != moduleId && isActiveModule[moduleToDeleteId] == 1)
		{
			GetCombinedMatrix(index, moduleId, moduleToDeleteId, newModule);   
			newModule.normalize_columns();

			vnl_matrix<double> mat = newModule.transpose();
			GetMatrixRowMeanStd(mat, newModuleMeans, newModuleStd);

			newModuleMeans.normalize();
			newModuleCor = newModuleMeans * newModule;

			if( abs(newModuleCor.mean()) > cor)
			{
				isActiveModule[moduleToDeleteId] = 0;
				moduleSize[moduleId] += moduleSize[moduleToDeleteId];
				moduleSize[moduleToDeleteId] = 0;
				mean.set_column(moduleId, newModuleMeans);
				mean.set_column(moduleToDeleteId, zeroCol);

				for( unsigned int j = 0; j < index.size(); j++)
				{
					if( index[j] == moduleToDeleteId)
					{
						index[j] = moduleId;
					}
				}
			}
		}
		isActiveModule[moduleId] = 0;
	}

	StandardizeIndex(index);

	int max = index.max_value() + 1;
	EraseZeroCol( mean);

	return index.max_value() + 1;   // how many clusters have been generated
}

void SPDAnalysisModel::GetCombinedMatrix(vnl_vector<unsigned int>& index, unsigned int moduleId, unsigned int moduleDeleteId, 
										 vnl_matrix<double>& mat)
{
	unsigned int i = 0;
	unsigned int ind = 0;

	for( i = 0; i < index.size(); i++)
	{
		if( index[i] == moduleId || index[i] == moduleDeleteId)
		{
			mat.set_column( ind++, this->DataMatrix.get_column(i));
		}
	}
}

vnl_vector<int> SPDAnalysisModel::GetModuleSize(vnl_vector<unsigned int>& index)   // index starts from 0
{
	vnl_vector<int> moduleSize;
	moduleSize.set_size(index.max_value() + 1);
	for( unsigned int i = 0; i < index.size(); i++)
	{
		moduleSize[index[i]]++;
	}
	return moduleSize;
}

int SPDAnalysisModel::GetSingleModuleSize(vnl_vector<unsigned int>& index, unsigned int ind)
{
	int size = 0;
	for( unsigned int i = 0; i < index.size(); i++)
	{
		if( index[i] == ind)
		{
			size++;
		}
	}
	return size;
}

void SPDAnalysisModel::GetMatrixRowMeanStd(vnl_matrix<double>& mat, vnl_vector<double>& mean, vnl_vector<double>& std)
{
	mbl_stats_nd stats;

	for( unsigned int i = 0; i < mat.rows() ; ++i)
	{
		vnl_vector<double> temp_row = mat.get_row(i);
		stats.obs(temp_row);	
	}

	std = stats.sd();
	mean = stats.mean();
}

void SPDAnalysisModel::StandardizeIndex(vnl_vector<unsigned int>& index)
{
	vnl_vector<unsigned int> newIndex = index;
	vnl_vector<unsigned int> sortY( index.size());
	vnl_vector<unsigned int> sortI( index.size());
	
	for( unsigned int i = 0; i < index.size(); i++)
	{
		sortY[i] = newIndex.min_value();
		sortI[i] = newIndex.arg_min();

		newIndex[ sortI[i]] = -1;
	}

	vnl_vector<unsigned int> indicator = index;
	indicator[0] = 0;
	for( unsigned int i = 1; i < index.size(); i++)
	{
		indicator[i] = (sortY[i] != sortY[i-1]);
 	}

	for( unsigned int i = 0; i < index.size(); i++)
	{
		newIndex[i] = 0;
	}
	for( unsigned int i = 0; i < index.size(); i++)
	{
		for( unsigned int j = 0; j <= i; j++)
		{
			newIndex[sortI[i]] += indicator[j];
		}
	}

	index = newIndex;
}

void SPDAnalysisModel::EraseZeroCol(vnl_matrix<double>& mat)
{
	int colNum = 0;
	for( unsigned int i = 0; i < mat.cols(); i++)
	{
		if( mat.get_column(i).is_zero() == false)
		{
			colNum++;
		}
	}

	vnl_matrix<double> newMat( mat.rows(), colNum);

	unsigned int newMatInd = 0;

	for( unsigned int j = 0; j < mat.cols(); j++)
	{
		if( mat.get_column(j).is_zero() == false)
		{
			newMat.set_column( newMatInd++, mat.get_column(j));
		}
	}

	mat = newMat;
}

void SPDAnalysisModel::ClusterMerge(double cor, double mer)
{
	vnl_matrix<double> coherence = this->ModuleMean.transpose() * this->ModuleMean;
	vnl_matrix<double> oneMat(coherence.rows(), coherence.cols());
	for( unsigned int i = 0; i < coherence.rows(); i++)
	{
		for( unsigned int j = 0; j < coherence.rows(); j++)
		{
			if( i == j)
			{
				oneMat(i,j) = 1;
			}
			else
			{
				oneMat(i,j) = 0;
			}
		}
	}
	coherence -= oneMat;

	//std::ofstream ofs("coherence.txt", std::ofstream::out);
	//ofs<< coherence<<endl;

	while( coherence.rows() > 0 && coherence.max_value() >= cor)
	{
		unsigned int maxId = coherence.arg_max();
		unsigned int j = maxId / coherence.rows();
		unsigned int i = maxId - coherence.rows() * j;

		vnl_matrix<double> dataTmp( this->DataMatrix.rows(), this->ModuleSize[i] + this->ModuleSize[j]);
		GetCombinedMatrix( this->ClusterIndex, i, j, dataTmp);
			
		vnl_vector<double> dataMean;
		vnl_vector<double> dataSd;
		vnl_matrix<double> mat = dataTmp.transpose();
		GetMatrixRowMeanStd( mat, dataMean, dataSd);
		dataMean.normalize();

		double moduleCor = (dataMean * dataTmp).mean();

		if( moduleCor >= cor || coherence(i, j) >= mer)
		{
			unsigned int min = i > j ? j : i;
			unsigned int max = i > j ? i : j;
			this->ModuleMean.set_column( min, dataMean);

			DeleteMatrixColumn( this->ModuleMean, max);
			DeleteMatrixColumn( coherence, max);

			vnl_matrix<double> coherenceTmp = coherence.transpose();
			DeleteMatrixColumn( coherenceTmp, max);
			coherence = coherenceTmp.transpose();

			SubstitudeVectorElement( this->ClusterIndex, max, min);

			unsigned int clusternum = this->ClusterIndex.max_value();

			for( unsigned int ind = max + 1; ind < clusternum + 1; ind++)
			{
				SubstitudeVectorElement( this->ClusterIndex, ind, ind - 1);
			}
		}
		else
		{
			break;
		}
	}
	
	//ofs.close();

	std::cout<<"The cluster size after merge: "<<this->ClusterIndex.max_value() + 1<<endl;
	std::cout<<this->ClusterIndex<<endl<<endl;

	std::ofstream ofs("result.txt", std::ofstream::out | ios::app);
	ofs<<"The cluster size after merge: "<<this->ClusterIndex.max_value() + 1<<endl;
	for ( unsigned int i = 0; i < this->ClusterIndex.size(); i++)
	{
		ofs<<this->ClusterIndex[i] + 1<<endl;
	}
	ofs<<this->ClusterIndex<<endl<<endl;
	ofs.close();
}

void SPDAnalysisModel::SubstitudeVectorElement( vnl_vector<unsigned int>& vector, unsigned int ori, unsigned int newValue)
{
	for( unsigned int i = 0; i < vector.size(); i++)
	{
		if( vector[i] == ori)
		{
			vector[i] = newValue;
		}
	}
}

void SPDAnalysisModel::DeleteMatrixColumn( vnl_matrix<double>& mat, unsigned int col)
{
	vnl_matrix<double> matTmp( mat.rows(), mat.cols()-1);
	unsigned int index = 0;
	for( unsigned int i = 0; i < mat.cols() ; i++)
	{
		if( i != col)
		{
			matTmp.set_column( index, mat.get_column(i));
		}
	}
	mat = matTmp;
}

void SPDAnalysisModel::GenerateMST()
{
	std::ofstream ofs("MST.txt", std::ofstream::out);
	//std::ofstream ofsmat("CLusterMat.txt", std::ofstream::out);

	unsigned int num_nodes = this->DataMatrix.rows();

	for( unsigned int i = 0; i <= this->ClusterIndex.max_value(); i++)
	{
		std::cout<<"Building MST for module "<<i<<endl;


		vnl_matrix<double> clusterMat( this->DataMatrix.rows(), GetSingleModuleSize( this->ClusterIndex, i));
		GetCombinedMatrix( this->ClusterIndex, i, i, clusterMat);

		//ofsmat<<"Module:"<< i + 1<<endl;
		//ofsmat<< clusterMat<<endl<<endl;

		Graph graph(num_nodes);

		for( unsigned int k = 0; k < num_nodes; k++)
		{
			for( unsigned int j = k + 1; j < num_nodes; j++)
			{
				double dist = CityBlockDist( clusterMat, k, j);   // need to be modified 
				//ofs<<dist<<endl;
				boost::add_edge(k, j, dist, graph);
			}
		}

		std::vector< boost::graph_traits< Graph>::vertex_descriptor> vertex(this->DataMatrix.rows());

		try
		{
			boost::prim_minimum_spanning_tree(graph, &vertex[0]);

		}
		catch(...)
		{
			std::cout<< "MST construction failure!"<<endl;
			exit(111);
		}

		this->ModuleMST.push_back(vertex);

		// construct vtk table and show it
		vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();

		for(int i = 0; i < this->headers.size(); i++)
		{		
			column = vtkSmartPointer<vtkDoubleArray>::New();
			column->SetName( (this->headers)[i].c_str());
			table->AddColumn(column);
		}
		
		double mapweight = 0;   // construct the vtk table and compute the overall weight
		for( int i = 0; i < vertex.size(); i++)
		{
			if( i != vertex[i])
			{
				double dist = CityBlockDist( clusterMat, i, vertex[i]);
				vtkSmartPointer<vtkVariantArray> DataRow = vtkSmartPointer<vtkVariantArray>::New();
				DataRow->InsertNextValue(i);
				DataRow->InsertNextValue( vertex[i]);
				DataRow->InsertNextValue(dist);
				table->InsertNextRow(DataRow);
				mapweight += dist;
			}
		}
		this->MSTTable.push_back(table);
		this->MSTWeight.push_back( mapweight);
	}

	std::vector<std::vector< boost::graph_traits<Graph>::vertex_descriptor> >::iterator iter = this->ModuleMST.begin();
	int moduleIndex = 1;
	while( iter != this->ModuleMST.end())
	{
		ofs<<"Module:"<<moduleIndex<<endl;
		ofs<<"Distance:"<< MSTWeight[moduleIndex - 1]<<endl;
		moduleIndex++;

		std::vector< boost::graph_traits<Graph>::vertex_descriptor> vertex = *iter;
		std::multimap<int,int> tree;

		for( unsigned int i = 0; i < vertex.size(); i++)
		{
			if( i != vertex[i])
			{
				tree.insert(std::pair<int,int>(i, vertex[i]));
				tree.insert(std::pair<int,int>(vertex[i], i));
			}
			//ofs<<i<<"\t"<< vertex[i]<<endl;
		}

		std::multimap<int,int>::iterator iterTree = tree.begin();
		while( iterTree!= tree.end())
		{
			ofs<<(*iterTree).second + 1<<"\t"<<((*iterTree).first) + 1<<endl;
			iterTree++;
		}

		iter++;	
	}

	ofs.close();
}

double SPDAnalysisModel::CityBlockDist( vnl_matrix<double>& mat, unsigned int ind1, unsigned int ind2)
{
	vnl_vector<double> module1 = mat.get_row(ind1);
	vnl_vector<double> module2 = mat.get_row(ind2);
	vnl_vector<double> subm = module1 - module2;
	double dis = 0;
	for( unsigned int i = 0; i < subm.size(); i++)
	{
		dis += abs(subm[i]);
	}
	return dis;
}

vtkSmartPointer<vtkTable> SPDAnalysisModel::GetMSTTable( int MSTIndex)
{
	if( MSTIndex >= 0 && MSTIndex < this->MSTTable.size())
	{
		vtkSmartPointer<vtkTable> table = this->MSTTable[MSTIndex];
		for( unsigned int i = 0; i < table->GetNumberOfRows(); i++)
		{
			long ind1 = this->MSTTable[MSTIndex]->GetValue(i,0).ToLong();
			long ind2 = this->MSTTable[MSTIndex]->GetValue(i,1).ToLong();
			table->SetValue(i, 0, this->indMapFromIndToVertex[ind1]);
			table->SetValue(i, 1, this->indMapFromIndToVertex[ind2]);
		}
		return table;
	}
	else
	{
		return NULL;
	}
}

void SPDAnalysisModel::GetTableHeaders(std::vector<std::string> &headers)
{
	headers = this->headers;
}

void SPDAnalysisModel::GetMatrixDistance( vnl_matrix<double>& data, vnl_vector<double>&distance, DISTANCE_TYPE type)
{
	unsigned int num = data.rows() * ( data.rows() - 1) / 2;
	distance.set_size( num);
	unsigned int ind = 0;
	
	switch( type)
	{
	case CITY_BLOCK:
		for( unsigned int i = 0; i < data.rows(); i++)
		{
			for( unsigned int j = i + 1; j < data.rows(); j++)
			{
				distance[ind] = 0;
				vnl_vector<double> sub = data.get_row(j) - data.get_row(i);
				for( unsigned int k = 0; k < sub.size(); k++)
				{
					distance[ind] += abs( sub[k]);
				}
				ind++;
			}
		}
		break;
	case SQUARE_DISTANCE:
		break;
	default:
		break;
	}
}

void SPDAnalysisModel::GetMSTMatrixDistance( vnl_vector<double>& distance, std::vector< boost::graph_traits<Graph>::vertex_descriptor>& vertex, vnl_vector<double>& MSTdistance)
{
	int size = vertex.size();
	MSTdistance.set_size( size - 1);   // vertex num is "size", MST edges num is "size - 1"
	int j = 0;
	for( int i = 0; i < size; i++)
	{
		if( i != vertex[i])
		{
			int min = vertex[i] > i ? i : vertex[i];
			int max = vertex[i] > i ? vertex[i] : i;
			int index = min * size - ( 1 + min) * min / 2 + ( max - min - 1);
			MSTdistance[j++] = distance[ index];
		}
	}
}

vnl_vector<unsigned int> SPDAnalysisModel::Hist(vnl_vector<double>&distance, int num_bin, vnl_vector<double>& interval)
{
	double inter = ( distance.max_value() - distance.min_value()) / num_bin;
	interval.set_size( num_bin + 1);
	for( unsigned int i = 0; i <= num_bin; i++)
	{
		interval[i] = distance.min_value() + i * inter;
	}

	vnl_vector<unsigned int> amount( num_bin);
	for( unsigned int i = 0; i < num_bin; i++)
	{
		amount[i] = 0;
	}

	for( unsigned int i = 0; i < distance.size(); i++)
	{
		int ind = floor((distance[i] - distance.min_value()) / inter);
		if( ind == num_bin)
		{
			ind =  num_bin - 1;
		}
		amount[ind] += 1;
	}
	return amount;
}

vnl_vector<unsigned int> SPDAnalysisModel::Hist(vnl_vector<double>&distance, vnl_vector<double>& interval)
{
	double inter = interval[1] - interval[0];   // interval value
	vnl_vector<unsigned int> distribution( interval.size() - 1);
	for( unsigned int i = 0; i < interval.size() - 1; i++)
	{
		distribution[i] = 0;
	}

	for( unsigned int i = 0; i < distance.size(); i++)
	{
		if( distance[i] >= interval.min_value() && distance[i] <= interval.max_value())
		{
			unsigned int index = floor( (distance[i] - interval.min_value()) / inter);
			if( index >= interval.size() - 1)
			{
				index = interval.size() - 2;
			}
			distribution[index]++;
		}
	}
	return distribution;
}

double Dist(int *first, int *second)
{
	return abs(*first - *second);
}

double SPDAnalysisModel::EarthMoverDistance(vnl_vector<unsigned int>& first, vnl_vector<unsigned int>& second)
{
	using namespace dt_simplex;
	int *src = new int[first.size()];
	double *src_num = new double[first.size()];
	unsigned int first_sum = first.sum();
	for( int i = 0; i < first.size(); i++)
	{
		src[i] = i;
		src_num[i] = (double)first[i] / first_sum;    // convert to the distribution percentage
	}

	int *snk = new int[second.size()];
	unsigned int second_sum =  second.sum();
	double *snk_num = new double[second.size()];
	for( int i = 0; i < second.size(); i++)
	{
		snk[i] = i;
		snk_num[i] = (double)second[i] / second_sum;    // convert to the distribution percentage
	}

	TsSignature<int> srcSig(first.size(), src, src_num);
	TsSignature<int> snkSig(second.size(), snk, snk_num);

	double result = transportSimplex(&srcSig, &snkSig, Dist);
	
	delete src;
	delete src_num;
	delete snk;
	delete snk_num;

	return result;
}

void SPDAnalysisModel::RunEMDAnalysis()
{
	vnl_vector<int> moduleSize = GetModuleSize( this->ClusterIndex);
	vnl_vector<double> distance;
	vnl_vector<double> mstdis;

	/** Generating the distance vector for module data and its MST data */
	std::vector<vnl_vector<double> > moduleDistance;
	std::vector<vnl_vector<double> > MSTDistance;

	std::vector<std::vector< boost::graph_traits<Graph>::vertex_descriptor> >::iterator iter = this->ModuleMST.begin();

	for( int i = 0; i <= this->ClusterIndex.max_value() && iter != this->ModuleMST.end(); i++, iter++)
	{
		vnl_matrix<double> clusterData( this->DataMatrix.rows(), moduleSize[i]);
		GetCombinedMatrix( this->ClusterIndex, i, i, clusterData);    /// Get the module data for cluster i
		
		GetMatrixDistance( clusterData, distance, CITY_BLOCK);
		moduleDistance.push_back( distance);
		
		std::vector< boost::graph_traits<Graph>::vertex_descriptor> vertex = *iter;
		GetMSTMatrixDistance( distance, vertex, mstdis);
		MSTDistance.push_back( mstdis);
	}

	std::ofstream histOfs("emd.txt", std::ofstream::out);
	this->EMDMatrix.set_size( this->ClusterIndex.max_value() + 1, this->ClusterIndex.max_value() + 1);

	for( int i = 0; i < moduleDistance.size(); i++)
	{
		vnl_vector<double> hist_interval;
		vnl_vector<unsigned int> moduleHist = Hist( moduleDistance[i], NUM_BIN, hist_interval);
		//histOfs << moduleHist << endl;
		for( int j = 0; j < MSTDistance.size(); j++)
		{	
			std::cout<< "matching module " << i << " with MST " << j<<endl;
			vnl_vector<unsigned int> mstHist = Hist( MSTDistance[j], hist_interval);   
			//histOfs << mstHist << endl;
			this->EMDMatrix( i, j) = EarthMoverDistance( moduleHist, mstHist);
		}
	}
	histOfs << EMDMatrix << endl;
	histOfs.close();
	std::cout<< "EMD matrix saved in file emd.txt"<<endl;
}