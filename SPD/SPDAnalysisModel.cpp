#include "SPDAnalysisModel.h"
#include <fstream>
#include <string>
#include <vtkVariant.h>
#include <math.h>
#include <mbl/mbl_stats_nd.h>
#include <boost/graph/prim_minimum_spanning_tree.hpp>

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
	graphWindow = new GraphWindow();
	headers.push_back("node1");
	headers.push_back("node2");
	headers.push_back("weight");
}

SPDAnalysisModel::~SPDAnalysisModel()
{
	if( graphWindow)
	{
		delete graphWindow;
	}
}

bool SPDAnalysisModel::ReadCellTraceFile(const char* fileName)
{
	int line = LineNum(fileName);

	std::ifstream file(fileName, std::ifstream::in);
	if( !file.is_open() || line == -1)
	{
		return false;
	}

	std::string feature;

	bool bfirst = true;
	int rowIndex = 0; 
	std::vector<std::string> rowValue;
	while( !file.eof())
	{
		getline(file, feature);
		if( feature.length() > 5)
		{
			split(feature, '\t', &rowValue);
			if( bfirst)
			{
				this->FeatureNames = rowValue;
				this->FeatureNames.erase(this->FeatureNames.begin());     // first trace index
				(this->FeatureNames).pop_back();
				(this->FeatureNames).pop_back();
				(this->FeatureNames).pop_back();                            // distance to device && cell names
				 bfirst = false;
				(this->DataMatrix).set_size(line - 2, this->FeatureNames.size());
			}
			else
			{
				std::vector<std::string>::iterator iter = rowValue.begin();
				for( int colIndex = 0, col = 0; col < rowValue.size() - 1; col++, iter++)
				{
					if ( col == 0)
					{
						this->CellTraceIndex.push_back( atoi((*iter).c_str()));    // trace index
					}
					else if( col == rowValue.size() - 3)		// strings to be eliminated
					{
						//break;
					}
					else if( col == rowValue.size() - 2)	
					{
						this->DistanceToDevice.push_back( atof((*iter).c_str()));  // save distance to device in the array
					}
					else if( col > 0 && col < rowValue.size() - 3)
					{
						(this->DataMatrix)(rowIndex, colIndex++) = atof( (*iter).c_str());   // save the feature data in the matrix for analysis
					}
				}
				rowIndex++;
			}
			rowValue.clear();
		}
		feature.clear();
	}
	//std::ofstream ofs("readdata.txt", std::ofstream::out);
	//ofs<<this->DataMatrix<<endl;
	return true;
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

			GetMatrixRowMeanStd(newModule.transpose(), newModuleMeans, newModuleStd);

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
		GetMatrixRowMeanStd( dataTmp.transpose(), dataMean, dataSd);
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
		//this->graphWindow->SetGraphTable(table, this->headers[0], this->headers[1], this->headers[2]);
		//this->graphWindow->ShowGraphWindow();
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

void SPDAnalysisModel::ShowMST()
{

}
