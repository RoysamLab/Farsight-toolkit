#include "SPDAnalysisModel.h"
#include <fstream>
#include <string>
#include <vtkVariant.h>
#include <math.h>
#include <mbl/mbl_stats_nd.h>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include "ftkUtils.h"
#include "transportSimplex.h"
#include <iomanip>
#ifdef _OPENMP
#include "omp.h"
#endif
#define NUM_THREAD 4
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
	filename = "";
	cc1 = new clusclus();
	cc2 = new clusclus();
	this->bProgression = false;
}

SPDAnalysisModel::~SPDAnalysisModel()
{
	if( cc1 != NULL)
	{
		delete cc1;
		cc1 = NULL;
	}
	if( cc2 != NULL)
	{
		delete cc2;
		cc2 = NULL;
	}
}

vtkSmartPointer<vtkTable> SPDAnalysisModel::GetDataTable()
{
	if( this->DataTable->GetNumberOfRows() > 0)
	{
		return this->DataTable;
	}
	else
	{
		return NULL;
	}
}

QString SPDAnalysisModel::GetFileName()
{
	return this->filename;
}

bool SPDAnalysisModel::ReadCellTraceFile(std::string fileName, bool btest)
{
	std::ifstream file(fileName.c_str(), std::ifstream::in);
	int line = LineNum(fileName.c_str());
	if( !file.is_open() || line == -1)
	{
		return false;
	}

	if( btest)
	{
		int rowIndex = 0;
		std::string feature;
		std::vector<std::string> rowValue;
		bool bfirst = true;
		while( !file.eof())
		{
			getline(file, feature);
			if( feature.length() > 3)
			{
				split(feature, '\t', &rowValue);
				std::vector<std::string>::iterator iter = rowValue.begin();
				if ( bfirst)
				{
					this->DataMatrix.set_size( line, rowValue.size() - 1);
					bfirst = false;
				}
				for( int colIndex = 0; colIndex < rowValue.size() - 1; iter++, colIndex++)
				{
					this->DataMatrix( rowIndex, colIndex) = atof( (*iter).c_str());   
				}
				rowIndex++;
			}
			rowValue.clear();
			feature.clear();
		}

		// build the index and its mapping for the test data
		for( unsigned int k = 0; k <= this->DataMatrix.cols(); k++)
		{
			vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
			this->DataTable->AddColumn(column);
		}

		for( unsigned int i = 0; i < this->DataMatrix.rows(); i++)
		{
			vtkSmartPointer<vtkVariantArray> row = vtkSmartPointer<vtkVariantArray>::New();

			for( unsigned int j = 0; j <= this->DataMatrix.cols(); j++)
			{ 
				if( j == 0)
				{
					row->InsertNextValue( vtkVariant( i));
				}
				else
				{
					row->InsertNextValue( vtkVariant(this->DataMatrix( i, j - 1)));
				}
			}
			this->DataTable->InsertNextRow(row);
		}
	}
	else
	{
		this->DataTable = ftk::LoadTable(fileName);
		if( this->DataTable != NULL)
		{
			ParseTraceFile(this->DataTable);
		}
		else
		{
			return false;
		}
	}	
	return true;
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
	std::cout<< table->GetNumberOfRows()<< "\t"<< table->GetNumberOfColumns()<<endl;
	for( long int i = 0; i < table->GetNumberOfRows(); i++)
	{
		long int var = table->GetValue( i, 0).ToLong();
		this->indMapFromIndToVertex.push_back( var);
	}

	this->DataTable = table;

	this->DataMatrix.set_size( this->DataTable->GetNumberOfRows(), this->DataTable->GetNumberOfColumns() - 2);
	for( int i = 0, rowIndex = 0; i < this->DataTable->GetNumberOfRows(); i++, rowIndex++)
	{
		int colIndex = 0;
		for( int j = 0; j < this->DataTable->GetNumberOfColumns(); j++)
		{
			if( j == 0 )
			{
				this->CellTraceIndex.push_back( this->DataTable->GetValue(i, j).ToInt());
			}
			else if( j == this->DataTable->GetNumberOfColumns() - 1)
			{
				this->DistanceToDevice.push_back( this->DataTable->GetValue(i, j).ToDouble());
			}
			else
			{
				double var = this->DataTable->GetValue(i, j).ToDouble();
				if( !boost::math::isnan(var))
				{
					(this->DataMatrix)(rowIndex, colIndex++) = this->DataTable->GetValue(i, j).ToDouble();
				}
				else
				{
					(this->DataMatrix)(rowIndex, colIndex++) = 0;
				}
			}
		}
	}

	this->CellClusterIndex.set_size(this->DataMatrix.rows());
	for( int i = 0; i < CellCluster.size(); i++)
	{
		this->CellCluster[i].clear();
	}
	this->CellCluster.clear();

	for( int i = 0; i < this->DataMatrix.rows(); i++)
	{
		this->CellClusterIndex[i] = i;
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
		if( abs(std_vec(i)) > 0)
		{	
			for( unsigned int j =0; j < temp_col.size(); ++j)
			{
				temp_col[j] = (temp_col[j] - mean_vec(i))/std_vec(i);
			}
			//temp_col = temp_col.normalize();
		}
		else
		{
			for( unsigned int j =0; j < temp_col.size(); ++j)
			{
				temp_col[j] = temp_col[j] - mean_vec(i);
			}
		}
		this->DataMatrix.set_column(i,temp_col);
	}
}

void SPDAnalysisModel::ClusterCells( double cor)
{
	this->filename = "C++_" + QString::number(this->DataMatrix.cols())+ "_" + QString::number(this->DataMatrix.rows()) + 
					"_" + QString::number( cor, 'g', 4) + "_";
	QString filenameCluster = this->filename + "Cellclustering.txt";
	std::ofstream ofs(filenameCluster.toStdString().c_str(), std::ofstream::out);

	vnl_matrix<double> tmpMat = this->DataMatrix;

	NormalizeData();   // eliminate the differences of different features

	this->DataMatrix =  this->DataMatrix.transpose();

	NormalizeData();   // for calculating covariance

	vnl_matrix<double> moduleMean = this->DataMatrix;
	moduleMean.normalize_columns();

	for( int i = 0; i < CellCluster.size(); i++)
	{
		CellCluster[i].clear();
	}
	CellCluster.clear();
	
	vnl_vector<unsigned int> clusterIndex;
	clusterIndex.set_size( this->DataMatrix.cols());
	vnl_vector<int> TreeIndex(this->DataMatrix.cols());
	for( unsigned int i = 0; i < this->DataMatrix.cols(); i++)
	{
		clusterIndex[i] = i;
		TreeIndex[i] = i;
	}

	ofs<<"cell clustering:"<<endl;

	TreeData.clear();

	int old_cluster_num = this->DataMatrix.cols();
	int new_cluster_num = this->DataMatrix.cols();
	do
	{
		old_cluster_num = new_cluster_num;

		new_cluster_num = ClusterAggFeatures( clusterIndex, moduleMean, cor, TreeIndex, 2);
		std::cout<< "new cluster num:"<< new_cluster_num<<" vs "<<"old cluster num:"<<old_cluster_num<<endl;
	}
	while( old_cluster_num != new_cluster_num);

	this->CellClusterIndex = clusterIndex;
	std::cout<< "Cell Cluster Size:"<<new_cluster_num<<endl;

	/// rebuild the datamatrix for MST 
	std::vector<int> clus;
	for( int i = 0; i < new_cluster_num; i++)
	{
		CellCluster.push_back(clus);
	}
	for( int i = 0; i < this->DataMatrix.cols(); i++)
	{
		int index = this->CellClusterIndex[i];
		CellCluster[index].push_back(i);
	}

	this->MatrixAfterCellCluster.set_size( CellCluster.size(), tmpMat.cols());
	for( int i = 0; i < CellCluster.size(); i++)
	{
		ofs<< "Cluster "<<i<<endl;
		vnl_vector<double> tmpVec( tmpMat.cols());
		tmpVec.fill(0);
		for( int j = 0; j < CellCluster[i].size(); j++)
		{
			ofs<< CellCluster[i][j]<<endl;
			tmpVec += tmpMat.get_row(CellCluster[i][j]);
		}
		tmpVec = tmpVec / CellCluster[i].size();
		MatrixAfterCellCluster.set_row( i, tmpVec);
	}	
	this->DataMatrix = MatrixAfterCellCluster;  // for normalization
	NormalizeData();
	MatrixAfterCellCluster = this->DataMatrix;
	
	this->DataMatrix = tmpMat;   // restore datamatrix
	ofs.close();
}

void SPDAnalysisModel::GetClusterMapping( std::vector<int> &indToclusInd)
{
	indToclusInd.clear();
	if( this->CellClusterIndex.max_value() + 1 < this->DataMatrix.rows())    // cell cluster has been on
	{
		for( int i = 0; i < this->CellClusterIndex.size(); i++)
		{
			indToclusInd.push_back( this->CellClusterIndex[i]);
		}
	}
}

int SPDAnalysisModel::ClusterAgglomerate(double cor, double mer)
{
	NormalizeData();
	vnl_vector<unsigned int> clusterIndex;
	clusterIndex.set_size( this->DataMatrix.cols());
	vnl_matrix<double> moduleMean = this->DataMatrix;
	moduleMean.normalize_columns();

	this->filename = "C++_" + QString::number(this->DataMatrix.cols())+ "_" + QString::number(this->DataMatrix.rows()) + 
					"_" + QString::number( cor, 'g', 4) + "_" + QString::number( mer, 'g', 4)+ "_";
	QString filenameCluster = this->filename + "clustering.txt";
	std::ofstream ofs(filenameCluster.toStdString().c_str(), std::ofstream::out);
	ofs<<"first clustering:"<<endl;

	vnl_vector<int> TreeIndex(this->DataMatrix.cols());
	for( unsigned int i = 0; i < this->DataMatrix.cols(); i++)
	{
		clusterIndex[i] = i;
		TreeIndex[i] = i;
	}
	TreeData.clear();

	int old_cluster_num = this->DataMatrix.cols();
	int new_cluster_num = this->DataMatrix.cols();

	do
	{
		old_cluster_num = new_cluster_num;
		new_cluster_num = ClusterAggFeatures( clusterIndex, moduleMean, cor, TreeIndex, 2);
		//ofs<<clusterIndex<<endl<<endl;
	}
	while( old_cluster_num != new_cluster_num);

	this->ClusterIndex = clusterIndex;
	this->ModuleMean = moduleMean;
	this->TreeIndex = TreeIndex;

	for( unsigned int i = 0; i <= this->ClusterIndex.max_value(); i++)
	{
		ofs<<"Cluster:"<< i + 1<<endl;
		for( unsigned int j = 0; j < this->ClusterIndex.size(); j++)
		{
			if( ClusterIndex[j] == i)
			{
				ofs<<j + 1<<endl;
			}
		}
	}
	ofs<<endl;
	ofs.close();

	std::cout<<"The cluster size after cluster: "<<this->ClusterIndex.max_value() + 1<<endl;
	std::cout<<clusterIndex<<endl;

	return new_cluster_num;
}


int SPDAnalysisModel::ClusterAggFeatures(vnl_vector<unsigned int>& index, vnl_matrix<double>& mean, double cor, vnl_vector<int>& TreeIndex, int fold = 2)
{
	vnl_vector<int> moduleSize = GetModuleSize(index);
	vnl_vector<int> isActiveModule;
	vnl_vector<double> moduleCenter;
	int newIndex = TreeIndex.max_value() + 1;

	isActiveModule.set_size( moduleSize.size());
	int i = 0;
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

	//ofstream ofs("moduleCor.txt");
	//ofs<< mean<<endl;
	std::vector<unsigned int> delColsVec;
	while( isActiveModule.sum() > 1)
	{
		vnl_vector<int> tmp;
		tmp.set_size( moduleSize.size());

		vnl_matrix<double> newModule;
		vnl_vector<double> newModuleMeans;
		vnl_vector<double> newModuleStd;
		vnl_vector<double> newModuleCor;
		vnl_vector<double> moduleCor; 

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
		moduleCor.set_size( mean.cols());
		moduleCor = moduleCenter * mean;
		moduleCor[moduleId] = 0;

		unsigned int moduleToDeleteId = moduleCor.arg_max();
		int newModuleSize = moduleSize[moduleId] + moduleSize[moduleToDeleteId];

		newModule.set_size(this->DataMatrix.rows(), newModuleSize);

		if( moduleToDeleteId != moduleId && isActiveModule[moduleToDeleteId] == 1)
		{
			GetCombinedMatrix( this->DataMatrix, index, moduleId, moduleToDeleteId, newModule);   
			newModule.normalize_columns();

			vnl_matrix<double> mat = newModule.transpose();
			GetMatrixRowMeanStd(mat, newModuleMeans, newModuleStd);

			newModuleMeans = newModuleMeans.normalize();
			newModuleCor = newModuleMeans * newModule;
			double newCor = abs( newModuleCor.mean());                          //??? there should be no abs???

			if( newCor > cor)
			{
				isActiveModule[moduleToDeleteId] = 0;
				moduleSize[moduleId] += moduleSize[moduleToDeleteId];
				moduleSize[moduleToDeleteId] = 0;

				mean.set_column(moduleId, newModuleMeans);
				mean.set_column(moduleToDeleteId, zeroCol);
				
				delColsVec.push_back(moduleToDeleteId);
				
				for( unsigned int j = 0; j < index.size(); j++)
				{
					if( index[j] == moduleToDeleteId)
					{
						index[j] = moduleId;
					}
				}

				Tree tr( TreeIndex[moduleId], TreeIndex[moduleToDeleteId], newCor, newIndex);
				TreeIndex[moduleId] = newIndex;
				TreeIndex[moduleToDeleteId] = -1;
				TreeData.push_back( tr);
				newIndex++;
			}
		}
		
		//ofs<< isActiveModule.sum()<<endl;
		//ofs<< moduleCenter<<endl;
		//ofs<< moduleCor<<endl<<endl;
		isActiveModule[moduleId] = 0;
		//moduleCor.clear();
	}
	//ofs.close();

	StandardizeIndex(index);
	EraseCols( mean, delColsVec);

    vnl_vector<int> newTreeIndex(index.max_value() + 1);
	int j = 0;
	for( int i = 0; i < TreeIndex.size(); i++)
	{
		if( TreeIndex[i] != -1 && j <= index.max_value())
		{
			newTreeIndex[j++] = TreeIndex[i];
		}
	}
	TreeIndex = newTreeIndex;

	return index.max_value() + 1;   // how many clusters have been generated
}

void SPDAnalysisModel::GetCombinedMatrix( vnl_matrix<double> &datamat, vnl_vector<unsigned int>& index, unsigned int moduleId, unsigned int moduleDeleteId, 
										 vnl_matrix<double>& mat)
{
	unsigned int i = 0;
	unsigned int ind = 0;

	for( i = 0; i < index.size(); i++)
	{
		if( index[i] == moduleId || index[i] == moduleDeleteId)
		{
			mat.set_column( ind++, datamat.get_column(i));
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

	if( mat.cols()> 0 && mean.size() <= 0)
	{
		std::cout<< "meanstd ERROR"<<endl;
	}
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

void SPDAnalysisModel::EraseCols(vnl_matrix<double>& mat, std::vector<unsigned int> vec)
{
	vnl_matrix<double> newMat( mat.rows(), mat.cols() - vec.size());
	unsigned int index = 0;
	for( unsigned int i = 0; i < mat.cols(); i++)
	{
		if( IsExist(vec, i) == false)
		{
			newMat.set_column( index++, mat.get_column(i));
		}
	}
	mat = newMat;
}

void SPDAnalysisModel::ClusterMerge(double cor, double mer)
{
	vnl_matrix<double> coherence = this->ModuleMean.transpose() * this->ModuleMean;
	vnl_vector<int> moduleSize = GetModuleSize( this->ClusterIndex);
	vnl_matrix<double> oneMat(coherence.rows(), coherence.cols());

	int newIndex = this->TreeIndex.max_value() + 1;

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

	//std::ofstream oft("coherence.txt", std::ofstream::out);
	//oft<< setprecision(4)<<coherence<<endl;
	//oft.close();

	while( coherence.rows() > 0 && coherence.max_value() >= cor)
	{
		unsigned int maxId = coherence.arg_max();
		unsigned int j = maxId / coherence.rows();
		unsigned int i = maxId - coherence.rows() * j;

		vnl_matrix<double> dataTmp( this->DataMatrix.rows(), moduleSize[i] + moduleSize[j]);
		GetCombinedMatrix( this->DataMatrix, this->ClusterIndex, i, j, dataTmp);
		dataTmp.normalize_columns();

		vnl_vector<double> dataMean;
		vnl_vector<double> dataSd;
		vnl_matrix<double> mat = dataTmp.transpose();    // bugs lying here
		GetMatrixRowMeanStd( mat, dataMean, dataSd);
		dataMean = dataMean.normalize();

		double moduleCor = abs( (dataMean * dataTmp).mean());

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

			moduleSize[ min] += moduleSize[max];
			moduleSize[ max] = 0;

			Tree tr( TreeIndex[min], TreeIndex[max], moduleCor, newIndex);
			TreeIndex[min] = newIndex;
			TreeData.push_back( tr);
			newIndex++;

			SubstitudeVectorElement( this->ClusterIndex, max, min);
			unsigned int clusternum = this->ClusterIndex.max_value();

			vnl_vector<int> newTreeIndex( TreeIndex.size() - 1);
			for( int k = 0; k < max; k++)
			{
				newTreeIndex[k] = TreeIndex[k];
			}
			for( unsigned int ind = max + 1; ind < clusternum + 1; ind++)
			{
				SubstitudeVectorElement( this->ClusterIndex, ind, ind - 1);
				moduleSize[ind - 1] = moduleSize[ind];
				newTreeIndex[ind - 1] = TreeIndex[ind];
			}
			moduleSize[ clusternum] = 0;
			TreeIndex = newTreeIndex;
		}
		else
		{
			break;
		}
	}

	std::cout<<"The cluster size after merge: "<< this->ClusterIndex.max_value() + 1<<endl;
	std::cout<<this->ClusterIndex<<endl<<endl;

	QString filenameCluster = this->filename + "clustering.txt";
	std::ofstream ofs( filenameCluster.toStdString().c_str(), std::ofstream::out | ios::app);
	ofs<< "clustering result after merge:"<<endl;
	for( unsigned int i = 0; i <= this->ClusterIndex.max_value(); i++)
	{
		ofs<<"Cluster: "<<i + 1<<endl;
		for( unsigned int j = 0; j < this->ClusterIndex.size(); j++)
		{
			if( ClusterIndex[j] == i)
			{
				ofs<< j + 1 <<endl;
			}
		}
	}

	ofs<<endl<<"TreeData:"<<endl;
	for( int i = 0; i < TreeData.size(); i++)
	{
		Tree tr = TreeData[i];
		ofs<< tr.first<< "\t"<< tr.second<< "\t"<< tr.cor<< "\t"<< tr.parent<<endl;
	}
	ofs<< TreeIndex<<endl;
	ofs.close();
}

void SPDAnalysisModel::HierachicalClustering()    
{
	assert(TreeIndex.size() == this->ClusterIndex.max_value() + 1);
	PublicTreeData = TreeData;

	vnl_vector<int> moduleSize = GetModuleSize(this->ClusterIndex);
	vnl_matrix<double> matrix = DataMatrix;
	matrix.normalize_columns();
	vnl_matrix<double> CorMat(TreeIndex.size(), TreeIndex.size());
	CorMat.fill(0);
	for( int i = 0; i < TreeIndex.size(); i++)
	{
		for( int j = i + 1; j < TreeIndex.size(); j++)
		{
			vnl_matrix<double> newModule;
			vnl_vector<double> newModuleMeans;
			vnl_vector<double> newModuleStd;
			vnl_vector<double> newModuleCor;

			int newModuleSize = moduleSize[i] + moduleSize[j];
			newModule.set_size(this->DataMatrix.rows(), newModuleSize);
			GetCombinedMatrix(this->DataMatrix, this->ClusterIndex, i, j, newModule);   
			newModule.normalize_columns();
			vnl_matrix<double> mat = newModule.transpose();
			GetMatrixRowMeanStd(mat, newModuleMeans, newModuleStd);
			newModuleMeans = newModuleMeans.normalize();
			newModuleCor = newModuleMeans * newModule;
			double newCor = abs( newModuleCor.mean());
			CorMat(i, j) = newCor;
			CorMat(j, i) = newCor;
		}
	}

	vnl_vector<double> zeroCol;
	zeroCol.set_size( CorMat.rows());
	for( int i = 0; i < CorMat.rows(); i++)
	{
		zeroCol[i] = 0;
	}

	std::vector<unsigned int> ids;
	TreeIndexNew = TreeIndex;

	int newIndex = TreeIndexNew.max_value() + 1;
	vnl_vector<unsigned int> cIndex = this->ClusterIndex;

	for( int count = 0; count < TreeIndexNew.size() - 1; count++)
	{
		int maxId = CorMat.arg_max();
		double maxCor = CorMat.max_value();
		int j = maxId / CorMat.rows();
		int i = maxId - CorMat.rows() * j;
		int min = i > j ? j : i;
		int max = i > j ? i : j;

		Tree tr( TreeIndexNew[min], TreeIndexNew[max], maxCor, newIndex);
		PublicTreeData.push_back( tr);
		TreeIndexNew[min] = newIndex;
		TreeIndexNew[max] = -1;
		newIndex += 1;

		SubstitudeVectorElement(cIndex, max, min);       
		moduleSize[min] += moduleSize[max];
		moduleSize[max] = 0;

		CorMat.set_column(min, zeroCol);
		CorMat.set_column(max, zeroCol);
		CorMat.set_row(min, zeroCol);
		CorMat.set_row(max, zeroCol);

		for( int k = 0; k < CorMat.cols(); k++)
		{
			if( k != min && TreeIndexNew[k] != -1)
			{
				vnl_matrix<double> newModule(DataMatrix.rows(), moduleSize(min) + moduleSize(k));
				vnl_vector<double> newModuleMeans;
				vnl_vector<double> newModuleStd;
				vnl_vector<double> newModuleCor;

				GetCombinedMatrix( this->DataMatrix, cIndex, min, k, newModule);   
				newModule.normalize_columns();

				vnl_matrix<double> mat = newModule.transpose();
				GetMatrixRowMeanStd(mat, newModuleMeans, newModuleStd);

				newModuleMeans = newModuleMeans.normalize();
				newModuleCor = newModuleMeans * newModule;
				double newCor = abs( newModuleCor.mean());

				CorMat( min, k) = newCor;
				CorMat( k, min) = newCor;
			}	
		}
	}

	QString str = "dendrogram_feature.txt";
	ofstream ofs(str.toStdString().c_str());
	ofs<< "Tree Index:"<<endl;
	ofs<< TreeIndexNew<<endl;
	for( int i = 0; i < PublicTreeData.size(); i++)
	{
		Tree tr = PublicTreeData[i];
		ofs<< tr.first<< "\t"<< tr.second<< "\t"<< tr.cor<< "\t"<< tr.parent<<endl;
	}
	ofs.close();
}

/// Hierachical Clustering for cols of data
void SPDAnalysisModel::HierachicalClustering(vtkSmartPointer<vtkTable> table, bool bcol)    
{
	ParseTraceFile(table);
	TreeData.clear(); 
	//vnl_matrix<double> tmpData;
	if( bcol ==  false)
	{
		DataMatrix = DataMatrix.transpose();
	}
	DataMatrix.normalize_columns();

	vnl_matrix<double> CorMat(DataMatrix.cols(), DataMatrix.cols());
	CorMat.fill(0);
	for( int i = 0; i < DataMatrix.cols(); i++)
	{
		vnl_vector<double> tmpColi = DataMatrix.get_column(i);
		for( int j = i + 1; j < DataMatrix.cols(); j++)
		{
			vnl_vector<double> tmpColj = DataMatrix.get_column(j);
			double cor = VnlVecMultiply(tmpColi, tmpColj);
			CorMat(i, j) = cor;
			CorMat(j, i) = cor;
		}
	}

	vnl_vector<int> index( DataMatrix.cols());    // node index
	vnl_vector<unsigned int> cIndex(DataMatrix.cols());    // cluster index
	vnl_vector<int> size(DataMatrix.cols());      // cluster size
	int newIndex = DataMatrix.cols();
	for( int i = 0; i < DataMatrix.cols(); i++)
	{
		index[i] = i;
		cIndex[i] = i;
		size[i] = 1;
	}

	vnl_vector<double> zeroCol;
	zeroCol.set_size( CorMat.rows());
	for( int i = 0; i < CorMat.rows(); i++)
	{
		zeroCol[i] = 0;
	}
	
	std::vector<unsigned int> ids;
	for( int count = 0; count < DataMatrix.cols() - 1; count++)
	{
		int maxId = CorMat.arg_max();
		double maxCor = CorMat.max_value();
		int j = maxId / CorMat.rows();
		int i = maxId - CorMat.rows() * j;
		int min = i > j ? j : i;
		int max = i > j ? i : j;

		Tree tr( index[min], index[max], maxCor, newIndex);
		TreeData.push_back( tr);

		index[min] = newIndex;
		index[max] = -1;
		SubstitudeVectorElement(cIndex, max, min);        // !!! get all index which is max to min

		CorMat.set_column(min, zeroCol);
		CorMat.set_column(max, zeroCol);
		CorMat.set_row(min, zeroCol);
		CorMat.set_row(max, zeroCol);

		for( int k = 0; k < CorMat.cols(); k++)
		{
			if( k != min && index[k] != -1)
			{
				ids.clear();
				ids.push_back(min);
				ids.push_back(max);
				ids.push_back(k);

				vnl_matrix<double> newModule(DataMatrix.rows(), size(min) + size(max) + size(k));
				vnl_vector<double> newModuleMeans;
				vnl_vector<double> newModuleStd;
				vnl_vector<double> newModuleCor;
				vnl_vector<double> moduleCor; 

				GetCombinedMatrix( this->DataMatrix, cIndex, ids, newModule);   
				//newModule.normalize_columns();

				vnl_matrix<double> mat = newModule.transpose();
				GetMatrixRowMeanStd(mat, newModuleMeans, newModuleStd);

				newModuleMeans = newModuleMeans.normalize();
				newModuleCor = newModuleMeans * newModule;
				double newCor = abs( newModuleCor.mean());

				CorMat( min, k) = newCor;
				CorMat( k, min) = newCor;
			}	
		}
		newIndex += 1;
		size(min) += size(max);
		size(max) = 0;
	}

	QString str = "dendrogram";
	if(bcol)
	{
		str += "_feature.txt";
	}
	else
	{
		str += "_sample.txt";
	}
	ofstream ofs(str.toStdString().c_str());
	for( int i = 0; i < TreeData.size(); i++)
	{
		Tree tr = TreeData[i];
		ofs<< tr.first<< "\t"<< tr.second<< "\t"<< tr.cor<< "\t"<< tr.parent<<endl;
	}
	ofs.close();
}

double SPDAnalysisModel::VnlVecMultiply(vnl_vector<double> const &vec1, vnl_vector<double> const &vec2)
{
	if( vec1.size() != vec2.size())
	{
		std::cout<< "vector size doesn't match for multiplication"<<endl;
		return 0;
	}
	double sum = 0;
	for( int i = 0; i < vec1.size(); i++)
	{
		sum += vec1[i] * vec2[i];
	}
	return sum;
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

void SPDAnalysisModel::SetProgressionTag(bool bProg)
{
	this->bProgression = bProg;
}

bool SPDAnalysisModel::GetProgressionTag()
{
	return this->bProgression;
}

void SPDAnalysisModel::GenerateMST()
{
	this->ModuleMST.clear();
	this->MSTWeight.clear();
	this->ModuleGraph.clear();

	if( this->MatrixAfterCellCluster.rows() <= 0)
	{
		this->MatrixAfterCellCluster = this->DataMatrix;
		std::cout<< "MatrixAfterCellCluster hasn't been generated, given datamatrix instead!"<<endl;
	}

	unsigned int num_nodes = this->MatrixAfterCellCluster.rows();

	std::vector< vnl_matrix<double> > clusterMatList;
	for( int i = 0; i <= this->ClusterIndex.max_value(); i++)
	{
		vnl_matrix<double> clusterMat( this->MatrixAfterCellCluster.rows(), GetSingleModuleSize( this->ClusterIndex, i));
		GetCombinedMatrix( this->MatrixAfterCellCluster, this->ClusterIndex, i, i, clusterMat);
		clusterMatList.push_back(clusterMat);
	}

	this->ModuleMST.resize( this->ClusterIndex.max_value() + 1);

	#ifdef _OPENMP
		omp_lock_t lock;
		omp_init_lock(&lock);
		omp_lock_t lock2;
		omp_init_lock(&lock2);
	#endif

	#pragma omp parallel for num_threads(NUM_THREAD)
	for( int i = 0; i <= this->ClusterIndex.max_value(); i++)
	{
		std::cout<<"Building MST for module " << i<<endl;

		//ofsmat<<"Module:"<< i + 1<<endl;
		//ofsmat<< clusterMat<<endl<<endl;
		#ifdef _OPENMP
			omp_set_lock(&lock);
		#endif
		vnl_matrix<double> clusterMat = clusterMatList[i];
		#ifdef _OPENMP
			omp_unset_lock(&lock);
		#endif

		Graph graph(num_nodes);

		for( unsigned int k = 0; k < num_nodes; k++)
		{
			for( unsigned int j = k + 1; j < num_nodes; j++)
			{
				double dist = CityBlockDist( clusterMat, k, j);   // need to be modified 
				//double dist = EuclideanBlockDist( clusterMat, k, j);
				boost::add_edge(k, j, dist, graph);
			}
		}

		std::vector< boost::graph_traits< Graph>::vertex_descriptor> vertex(this->MatrixAfterCellCluster.rows());

		try
		{
			boost::prim_minimum_spanning_tree(graph, &vertex[0]);

		}
		catch(...)
		{
			std::cout<< "MST construction failure!"<<endl;
			exit(111);
		}

		#ifdef _OPENMP
			omp_set_lock(&lock2);
		#endif
		this->ModuleMST[i] = vertex;
		#ifdef _OPENMP
			omp_unset_lock(&lock2);
		#endif
	}

	#ifdef _OPENMP
		omp_destroy_lock(&lock);
		omp_destroy_lock(&lock2);
	#endif

	for( int n = 0; n <= this->ClusterIndex.max_value(); n++)
	{
		vnl_matrix<double> clusterMat = clusterMatList[n];

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
		for( int i = 0; i < this->ModuleMST[n].size(); i++)
		{
			if( i != this->ModuleMST[n][i])
			{
				double dist = CityBlockDist( clusterMat, i, this->ModuleMST[n][i]);
				vtkSmartPointer<vtkVariantArray> DataRow = vtkSmartPointer<vtkVariantArray>::New();
				DataRow->InsertNextValue(i);
				DataRow->InsertNextValue( this->ModuleMST[n][i]);
				DataRow->InsertNextValue(dist);
				table->InsertNextRow(DataRow);
				mapweight += dist;
			}
		}
		this->MSTTable.push_back(table);
		this->MSTWeight.push_back( mapweight);
	}

	//std::vector<std::vector< boost::graph_traits<Graph>::vertex_descriptor> >::iterator iter = this->ModuleMST.begin();
	//int moduleIndex = 1;

	//QString filenameMST = this->filename + "mst.txt";
	//std::ofstream ofs(filenameMST.toStdString().c_str(), std::ofstream::out);
	//while( iter != this->ModuleMST.end())
	//{
	//	ofs<<"Module:"<<moduleIndex<<endl;
	//	ofs<<"Distance:"<< MSTWeight[moduleIndex - 1]<<endl;
	//	moduleIndex++;

	//	std::vector< boost::graph_traits<Graph>::vertex_descriptor> vertex = *iter;
	//	std::multimap<int,int> tree;

	//	for( unsigned int i = 0; i < vertex.size(); i++)
	//	{
	//		if( i != vertex[i])
	//		{
	//			tree.insert(std::pair<int,int>(i, vertex[i]));
	//			tree.insert(std::pair<int,int>(vertex[i], i));
	//		}
	//		//ofs<<i<<"\t"<< vertex[i]<<endl;
	//	}

	//	std::multimap<int,int>::iterator iterTree = tree.begin();
	//	int oldValue = (*iterTree).first;
	//	while( iterTree!= tree.end())
	//	{
	//		std::map<int, int> sortTree;
	//		while( iterTree!= tree.end() && (*iterTree).first == oldValue)
	//		{
	//			sortTree.insert( std::pair<int, int>((*iterTree).second, (*iterTree).first));
	//			iterTree++;
	//		}

	//		if( iterTree != tree.end())
	//		{
	//			oldValue = (*iterTree).first;
	//		}
	//		
	//		std::map<int,int>::iterator iterSort = sortTree.begin();
	//		while( iterSort != sortTree.end())
	//		{
	//			ofs<<(*iterSort).first + 1<<"\t"<<((*iterSort).second) + 1<<endl;
	//			iterSort++;
	//		}
	//	}
	//	iter++;	
	//}
	//ofs.close();
}

void SPDAnalysisModel::GenerateDistanceMST()
{
	if( this->MatrixAfterCellCluster.rows() <= 0)
	{
		this->MatrixAfterCellCluster = this->DataMatrix;
	}

	unsigned int num_nodes = this->MatrixAfterCellCluster.rows();

	std::cout<<"Building MST for distance to device "<<endl;

	vnl_matrix<double> clusterMat( this->MatrixAfterCellCluster.rows(), 1);
	for( int i = 0; i < this->MatrixAfterCellCluster.rows(); i++)
	{
		clusterMat(i,0) = DistanceToDevice[i];
	}

	Graph graph(num_nodes);

	for( unsigned int k = 0; k < num_nodes; k++)
	{
		for( unsigned int j = k + 1; j < num_nodes; j++)
		{
			double dist = CityBlockDist( clusterMat, k, j);   // need to be modified 
				//double dist = EuclideanBlockDist( clusterMat, k, j);
			boost::add_edge(k, j, dist, graph);
		}
	}

	std::vector< boost::graph_traits< Graph>::vertex_descriptor> vertex(this->MatrixAfterCellCluster.rows());

	try
	{
		boost::prim_minimum_spanning_tree(graph, &vertex[0]);
	}
	catch(...)
	{
		std::cout<< "MST construction failure!"<<endl;
		exit(111);
	}

	DistanceMST = vertex;
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

double SPDAnalysisModel::EuclideanBlockDist( vnl_matrix<double>& mat, unsigned int ind1, unsigned int ind2)
{
	vnl_vector<double> module1 = mat.get_row(ind1);
	vnl_vector<double> module2 = mat.get_row(ind2);
	vnl_vector<double> subm = module1 - module2;
	double dis = subm.magnitude();
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
			double dis = this->MSTTable[MSTIndex]->GetValue(i,2).ToDouble();
			if( this->indMapFromIndToVertex.size() > 0)
			{
				table->SetValue(i, 0, this->indMapFromIndToVertex[ind1]);
				table->SetValue(i, 1, this->indMapFromIndToVertex[ind2]);
				table->SetValue(i, 2, dis);
			}
			else
			{
				table->SetValue(i, 0, ind1);
				table->SetValue(i, 1, ind2);
				table->SetValue(i, 2, dis);
			}
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

void SPDAnalysisModel::Hist(vnl_vector<double>&distance, int num_bin, vnl_vector<double>& interval, vnl_vector<unsigned int>& histDis)
{
	double min = distance.min_value();
	double inter = ( distance.max_value() - min) / num_bin;
	interval.set_size( num_bin + 1);
	for( unsigned int i = 0; i <= num_bin; i++)
	{
		interval[i] = distance.min_value() + i * inter;
	}

	histDis.set_size(num_bin);
	for( unsigned int i = 0; i < num_bin; i++)
	{
		histDis[i] = 0;
	}

	for( unsigned int i = 0; i < distance.size(); i++)
	{
		int ind = 0;
		if( inter != 0)
		{
			ind = floor((distance[i] - min) / inter);
		}
		else
		{
			ind = floor(distance[i] - min);
		}

		if( ind == num_bin)
		{
			ind =  num_bin - 1;
		}
		histDis[ind] += 1;
	}
}

void SPDAnalysisModel::Hist(vnl_vector<double>&distance, vnl_vector<double>& interval, vnl_vector<unsigned int>& histDis)
{
	double inter = interval[1] - interval[0];   // interval value
	histDis.set_size( interval.size() - 1);
	for( unsigned int i = 0; i < interval.size() - 1; i++)
	{
		histDis[i] = 0;
	}

	for( unsigned int i = 0; i < distance.size(); i++)
	{
		if( distance[i] >= interval.min_value() && distance[i] <= interval.max_value())
		{
			unsigned int index = 0;
			if( inter != 0)
			{
				index = floor( (distance[i] - interval.min_value()) / inter);
			}
			else
			{
				index = floor( distance[i] - interval.min_value());
			}

			if( index >= interval.size() - 1)
			{
				index = interval.size() - 2;
			}
			histDis[index]++;
		}
	}
}

double Dist(int *first, int *second)
{
	if( *first > *second)
	{
		return *first - *second;
	}
	else
	{
		return 0;
	}
	//return abs(*first - *second);
}

double SPDAnalysisModel::EarthMoverDistance(vnl_vector<unsigned int>& first, vnl_vector<unsigned int>& second,
											vnl_matrix<double> &flowMatrix)
{
	using namespace dt_simplex;
	int *src = new int[first.size()];
	double *src_num = new double[first.size()];
	unsigned int first_sum = first.sum();
	for( int i = 0; i < first.size(); i++)
	{
		src[i] = i;
		if( first_sum != 0)
		{
			src_num[i] = (double)first[i] / first_sum;    // convert to the distribution percentage
		}
		else
		{
			src_num[i] = (double)first[i];
		}
	}

	int *snk = new int[second.size()];
	unsigned int second_sum =  second.sum();
	double *snk_num = new double[second.size()];
	for( int i = 0; i < second.size(); i++)
	{
		snk[i] = i;
		if( second_sum != 0)
		{
			snk_num[i] = (double)second[i] / second_sum;    // convert to the distribution percentage
		}
		else
		{
			snk_num[i] = (double)second[i];
		}
	}

	TsSignature<int> srcSig(first.size(), src, src_num);
	TsSignature<int> snkSig(second.size(), snk, snk_num);

	TsFlow flow[NUM_BIN * NUM_BIN];
	int flowVars = 0;
	double result = transportSimplex(&srcSig, &snkSig, Dist, flow, &flowVars);
	
	for( int i = 0; i < NUM_BIN; i++)
	{
		for( int j = 0; j < NUM_BIN; j++)
		{
			flowMatrix( i, j) = 0;
		}
	}
	for( int i = 0; i < flowVars; i++)
	{
		flowMatrix( flow[i].to, flow[i].from) = flow[i].amount;
	}
	
	delete src;
	delete src_num;
	delete snk;
	delete snk_num;

	return result;
}

void SPDAnalysisModel::RunEMDAnalysis()
{
	vnl_vector<int> moduleSize = GetModuleSize( this->ClusterIndex);

	/** Generating the distance vector for module data and its MST data */
	std::vector<vnl_vector<double> > moduleDistance;
	moduleDistance.resize( this->ModuleMST.size());

	#pragma omp parallel for num_threads(NUM_THREAD)
	for( int i = 0; i <= this->ClusterIndex.max_value(); i++)
	{
		vnl_vector<double> distance;
		std::cout<< "build module distance "<< i<<endl;
		vnl_matrix<double> clusterData( this->MatrixAfterCellCluster.rows(), moduleSize[i]);
		GetCombinedMatrix( this->MatrixAfterCellCluster, this->ClusterIndex, i, i, clusterData);    /// Get the module data for cluster i
		GetMatrixDistance( clusterData, distance, CITY_BLOCK);
		moduleDistance[i] = distance;
	}

	//QString filenameEMD = this->filename + "emd.txt";
	//std::ofstream histOfs(filenameEMD.toStdString().c_str(), std::ofstream::out);
	//histOfs.precision(4);
	this->EMDMatrix.set_size( this->ClusterIndex.max_value() + 1, this->ClusterIndex.max_value() + 1);
	this->DistanceEMDVector.set_size( this->ClusterIndex.max_value() + 1);
	this->EMDMatrix.fill(0);
	this->DistanceEMDVector.fill(0);

	if( bProgression)
	{
		vnl_vector<double> distance;
		vnl_matrix<double> clusterMat( this->MatrixAfterCellCluster.rows(), 1);
		for( int i = 0; i < this->MatrixAfterCellCluster.rows(); i++)
		{
			clusterMat(i,0) = DistanceToDevice[i];
		}
		GetMatrixDistance( clusterMat, distance, CITY_BLOCK);

		vnl_vector<double> hist_interval;
		vnl_vector<unsigned int> moduleHist;
		Hist( distance, NUM_BIN, hist_interval, moduleHist);

		for( int j = 0; j < this->ModuleMST.size(); j++)
		{
			vnl_vector<double> mstDistance; 
			GetMSTMatrixDistance( distance, this->ModuleMST[j], mstDistance);  // get mst j distance from module i's distance matrix
			vnl_vector<unsigned int> mstHist;
			Hist( mstDistance, hist_interval, mstHist); 
			vnl_matrix<double> flowMatrix( NUM_BIN, NUM_BIN);		
			DistanceEMDVector[j] = EarthMoverDistance( moduleHist, mstHist, flowMatrix);
		}
		double max = DistanceEMDVector.max_value();
		DistanceEMDVector = DistanceEMDVector / max;
	}

	#pragma omp parallel for num_threads(NUM_THREAD)
	for( int i = 0; i < moduleDistance.size(); i++)
	{
		vnl_vector<double> hist_interval;
		vnl_vector<unsigned int> moduleHist;
		Hist( moduleDistance[i], NUM_BIN, hist_interval, moduleHist);

		//histOfs << "Module "<< i + 1 <<endl;
		
		// test the distance of the modules
		//for( int s = 0; s < moduleDistance[i].size(); s = s + 10)
		//{
		//	for( int t = 0; t < 10; t++)
		//	{
		//		if( s + t < moduleDistance[i].size())
		//		{
		//			histOfs << moduleDistance[i][s + t]<<"\t";
		//		}
		//	}
		//	histOfs<<endl;		
		//}
		//histOfs <<endl;
		//histOfs << moduleHistPer<< endl<<endl;
		for( int j = 0; j < moduleDistance.size(); j++)
		{	
			std::cout<< "matching module " << i<< " with MST " << j<<endl;
			//histOfs << "MST "<< j + 1 <<endl;

			vnl_vector<double> mstDistance; 
			GetMSTMatrixDistance( moduleDistance[i], this->ModuleMST[j], mstDistance);  // get mst j distance from module i's distance matrix
			vnl_vector<unsigned int> mstHist;
			Hist( mstDistance, hist_interval, mstHist); 

			//test the distance of the MSTs
			//for( int s = 0; s < mstDistance.size(); s = s + 10)
			//{
			//	for( int t = 0; t < 10; t++)
			//	{
			//		if( s + t < mstDistance.size())
			//		{
			//			histOfs << mstDistance[s + t]<<"\t";
			//		}
			//	}
			//	histOfs<<endl;		
			//}
			//histOfs << setiosflags(ios::fixed) << mstHistPer<< endl<<endl;

			vnl_matrix<double> flowMatrix( NUM_BIN, NUM_BIN);
			
			this->EMDMatrix( i, j) = EarthMoverDistance( moduleHist, mstHist, flowMatrix);

			/*histOfs << "flows:"<<endl;
			vnl_matrix<double> tranMatrix = flowMatrix.transpose();
			histOfs << setiosflags(ios::fixed) << tranMatrix <<endl <<endl;
			
			histOfs << "sums:" <<endl;
			vnl_vector<double> sumvector( tranMatrix.cols());
			for( int i = 0; i < tranMatrix.cols(); i++)
			{
				sumvector[i] = tranMatrix.get_column(i).sum();
			}

			histOfs << setiosflags(ios::fixed) << sumvector <<endl<<endl;*/
		}
	}

	//histOfs << setiosflags(ios::fixed) << this->EMDMatrix <<endl <<endl;

	for( unsigned int i = 0; i < this->EMDMatrix.rows(); i++)
	{
		double max = this->EMDMatrix.get_row(i).max_value();
		for( unsigned int j = 0; j < this->EMDMatrix.cols(); j++)
		{
			if( max != 0)
			{
				this->EMDMatrix( i, j) = this->EMDMatrix( i, j) / max;
			}
		}
	}

	//histOfs << setiosflags(ios::fixed) << this->EMDMatrix << endl;
	//histOfs.close();
	std::cout<< "EMD matrix has been built successfully"<<endl;
}

void SPDAnalysisModel::GetEMDMatrixDivByMax(vnl_matrix<double> &emdMatrix)
{
	emdMatrix =  this->EMDMatrix;
}

void SPDAnalysisModel::GetClusClusData(clusclus& c1, clusclus& c2, double threshold, double magFactor)
{
	this->heatmapMatrix.set_size( this->EMDMatrix.rows(), this->EMDMatrix.cols());
	this->heatmapMatrix.fill(0);
	std::vector< unsigned int> simModIndex;
	std::vector< unsigned int> disModIndex;

	if( bProgression)
	{
		for(unsigned int i = 0; i < DistanceEMDVector.size(); i++)
		{
			if( DistanceEMDVector[i] >= threshold)
			{
				disModIndex.push_back(i);
			}
		}
	}

	QString filenameSM = this->filename + "similarity_matrix.txt";
	std::ofstream ofSimatrix(filenameSM .toStdString().c_str(), std::ofstream::out);
	ofSimatrix.precision(4);

	for( unsigned int i = 0; i < this->EMDMatrix.cols(); i++)
	{
		ofSimatrix <<"MST "<<i<<endl;
		for( unsigned int j = 0; j < this->EMDMatrix.rows(); j++)
		{
			if( this->EMDMatrix( j, i) >= threshold)     // find the modules that support common mst
			{
				simModIndex.push_back(j);
				ofSimatrix << j <<"\t";
			}
		}
		ofSimatrix <<endl;

		if( simModIndex.size() > 0)
		{
			for( unsigned int j = 0; j < simModIndex.size(); j++)
			{
				for( unsigned int k = j + 1; k < simModIndex.size(); k++)
				{
					this->heatmapMatrix( simModIndex[j], simModIndex[k]) = heatmapMatrix( simModIndex[j], simModIndex[k]) + 1;
					this->heatmapMatrix( simModIndex[k], simModIndex[j]) = heatmapMatrix( simModIndex[k], simModIndex[j]) + 1;
				}
				this->heatmapMatrix(simModIndex[j], simModIndex[j]) = this->heatmapMatrix(simModIndex[j], simModIndex[j]) + 1;
			}
		}
		simModIndex.clear();
	}

	if( bProgression && disModIndex.size() > 0 && magFactor > 0)
	{
		for( unsigned int j = 0; j < disModIndex.size(); j++)
		{
			for( unsigned int k = j + 1; k < disModIndex.size(); k++)
			{
				this->heatmapMatrix( disModIndex[j], disModIndex[k]) = heatmapMatrix( disModIndex[j], disModIndex[k]) + magFactor;
				this->heatmapMatrix( disModIndex[k], disModIndex[j]) = heatmapMatrix( disModIndex[k], disModIndex[j]) + magFactor;
			}
			this->heatmapMatrix(disModIndex[j], disModIndex[j]) = this->heatmapMatrix(disModIndex[j], disModIndex[j]) + magFactor;
		}
	}

	ofSimatrix << "EMD matrix( devided by max value of each row):"<<endl;
	ofSimatrix << setiosflags(ios::fixed) << this->EMDMatrix <<endl<<endl;
	ofSimatrix << "Heatmap Matrix:"<<endl;
	ofSimatrix << this->heatmapMatrix <<endl<<endl;

	this->cc1->Initialize( heatmapMatrix.data_array(), this->EMDMatrix.rows(), this->EMDMatrix.cols());
	this->cc1->RunClusClus();
	this->cc1->Transpose();

	ofSimatrix << "Optimal leaf order:"<<endl;
	for(int i = 0; i < this->EMDMatrix.rows(); i++)
	{
		ofSimatrix << cc1->optimalleaforder[i]<<endl;
	}

	this->cc2->Initialize( cc1->transposefeatures,cc1->num_features, cc1->num_samples);
	this->cc2->RunClusClus();

	c1 = *this->cc1;
	c2 = *this->cc2;

	ofSimatrix.close();
}

void SPDAnalysisModel::GetCombinedMatrix( vnl_matrix<double> &datamat, vnl_vector< unsigned int> & index, std::vector< unsigned int> moduleID, vnl_matrix<double>& mat)
{
	unsigned int i = 0;
	unsigned int ind = 0;

	if( mat.rows() == datamat.rows())
	{
		for( i = 0; i < index.size(); i++)
		{
			if( IsExist(moduleID, index[i]))
			{
				mat.set_column( ind++, datamat.get_column(i));
			}
		}
	}
	else
	{
		std::cout<< "Row number does not match!"<<endl;
	}
}

bool SPDAnalysisModel::IsExist(std::vector<unsigned int> vec, unsigned int value)
{
	for( int i = 0; i < vec.size(); i++)
	{
		if( value == vec[i])
		{
			return true;
		}
	}
	return false;
}

long int SPDAnalysisModel::GetSelectedFeatures( vnl_vector< unsigned int> & index, std::vector<unsigned int> moduleID, std::set<long int>& featureSelectedIDs)
{
	long int count = 0;
	featureSelectedIDs.clear();
	for( long int i = 0; i < index.size(); i++)
	{
		if( IsExist(moduleID, index[i]))
		{
			count++;
			featureSelectedIDs.insert(i);   // feature id is one bigger than the original feature id 
		}
	}
	return count;
}

vtkSmartPointer<vtkTable> SPDAnalysisModel::GenerateProgressionTree( std::string& selectedModules)
{
	QString filenameMST = this->filename + "progression_mst.txt";
	std::ofstream ofs(filenameMST.toStdString().c_str(), std::ofstream::out);
	ofs.precision(4);

	std::vector<std::string> moduleIDStr;
	std::vector<unsigned int> moduleID;
	split( selectedModules, ',', &moduleIDStr);
	int coln = 0;
	for( int i = 0; i < moduleIDStr.size(); i++)
	{
		unsigned int index = atoi( moduleIDStr[i].c_str());
		moduleID.push_back( index);
	}
	coln = GetSelectedFeatures( this->ClusterIndex, moduleID, this->selectedFeatureIDs);
	
	vnl_matrix<double> clusterMat( this->MatrixAfterCellCluster.rows(), coln);
	GetCombinedMatrix( this->MatrixAfterCellCluster, this->ClusterIndex, moduleID, clusterMat);
	ofs<< setiosflags(ios::fixed)<< clusterMat<<endl<<endl;

	std::cout<<"Build progression tree"<<endl;
	int num_nodes = this->MatrixAfterCellCluster.rows();
	Graph graph( num_nodes);

	for( unsigned int k = 0; k < num_nodes; k++)
	{
		for( unsigned int j = k + 1; j < num_nodes; j++)
		{
			double dist = CityBlockDist( clusterMat, k, j);   // need to be modified 
			//double dist = EuclideanBlockDist( clusterMat, k, j);
			boost::add_edge(k, j, dist, graph);
		}
	}

	std::vector< boost::graph_traits< Graph>::vertex_descriptor> vertex(this->MatrixAfterCellCluster.rows());

	try
	{
		boost::prim_minimum_spanning_tree(graph, &vertex[0]);

	}
	catch(...)
	{
		std::cout<< "MST construction failure!"<<endl;
		exit(111);
	}

	// construct vtk table and show it

	std::cout<<"Construct vtk table"<<endl;
	vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
	vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();

	for(int i = 0; i < this->headers.size(); i++)
	{		
		column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName( (this->headers)[i].c_str());
		table->AddColumn(column);
	}
	
	for( int i = 0; i < vertex.size(); i++)
	{
		if( i != vertex[i])
		{
			double dist = CityBlockDist( clusterMat, i, vertex[i]);
			vtkSmartPointer<vtkVariantArray> DataRow = vtkSmartPointer<vtkVariantArray>::New();
			if( this->CellCluster.size() > 0)    // cluster index continuous, no mapping needed
			{
				DataRow->InsertNextValue( i);
				DataRow->InsertNextValue( vertex[i]);
			}
			else if( this->indMapFromIndToVertex.size() > 0)
			{
				DataRow->InsertNextValue( this->indMapFromIndToVertex[i]);
				DataRow->InsertNextValue( this->indMapFromIndToVertex[vertex[i]]);
			}
			else
			{
				DataRow->InsertNextValue( i);
				DataRow->InsertNextValue( vertex[i]);
			}
			DataRow->InsertNextValue(dist);
			table->InsertNextRow(DataRow);
		}
	}

	// write into files
	std::cout<<"Write into files"<<endl;
	std::multimap<int,int> tree;

	for( unsigned int i = 0; i < vertex.size(); i++)
	{
		if( i != vertex[i])
		{
			tree.insert(std::pair<int,int>(i, vertex[i]));
			tree.insert(std::pair<int,int>(vertex[i], i));
		}
	}

	std::multimap<int,int>::iterator iterTree = tree.begin();
	int oldValue = (*iterTree).first;
	while( iterTree!= tree.end())
	{
		std::map<int, int> sortTree;
		while( iterTree!= tree.end() && (*iterTree).first == oldValue)
		{
			sortTree.insert( std::pair<int, int>((*iterTree).second, (*iterTree).first));
			iterTree++;
		}

		if( iterTree != tree.end())
		{
			oldValue = (*iterTree).first;
		}
			
		std::map<int,int>::iterator iterSort = sortTree.begin();
		while( iterSort != sortTree.end())
		{
			ofs<<(*iterSort).first + 1<<"\t"<<((*iterSort).second) + 1<<endl;
			iterSort++;
		}
	}
	ofs.close();

	return table;
}

void SPDAnalysisModel::GetSelectedFeatures(std::set<long int>& selectedFeatures)
{
	selectedFeatures = this->selectedFeatureIDs;
}

void SPDAnalysisModel::SaveSelectedFeatureNames(QString filename, std::set<long int>& selectedFeatures)
{
	std::set<long int>::iterator iter;
	QString file = this->filename + filename;
	std::ofstream ofs( file.toStdString().c_str(), std::ofstream::out);

	for( iter = selectedFeatures.begin(); iter != selectedFeatures.end(); iter++)
	{
		long int index = *iter;
		char* name = this->DataTable->GetColumn(index + 1)->GetName();
		std::string featureName(name);
		ofs<< index<<"\t"<<featureName<<endl;
	}
	ofs.close();
}

double SPDAnalysisModel::GetEMDSelectedPercentage(double thres)
{
	unsigned int count = 0;
	for( unsigned int i = 0; i < this->EMDMatrix.cols(); i++)
	{
		for( unsigned int j = 0; j < this->EMDMatrix.rows(); j++)
		{
			if( this->EMDMatrix( j, i) >= thres)     // find the modules that support common mst
			{
				count++;
			}
		}
	}
	unsigned int allnum = this->EMDMatrix.cols() * this->EMDMatrix.rows();
	return (double)count / allnum;
}

double SPDAnalysisModel::GetEMDSelectedThreshold( double per)
{
	return 0;
}

void SPDAnalysisModel::GetMatrixData(vnl_matrix<double> &mat)
{
	mat = DataMatrix.normalize_columns();
}

void SPDAnalysisModel::GetDistanceOrder(std::vector<long int> &order)
{
	order.clear();
	vnl_vector<double> distance( DistanceToDevice.size());
	for( long int i = 0; i < DistanceToDevice.size(); i++)
	{
		distance[i] = DistanceToDevice[i];
	}
	double max = distance.max_value();
	for( long int i = 0; i < DistanceToDevice.size(); i++)
	{
		long int minind = distance.arg_min();
		order.push_back(minind);
		distance[minind] = max;
	}
}