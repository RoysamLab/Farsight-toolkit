#include "SPDAnalysisModel.h"
#include <fstream>
#include <string>
#include <vtkVariant.h>
#include <stdlib.h>
#include <math.h>
#include <mbl/mbl_stats_nd.h>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include "ftkUtils.h"
#include "transportSimplex.h"
#include <iomanip>
#include "ClusClus/Biclustering.h"

#ifdef _OPENMP
#include "omp.h"
#endif
#define NUM_THREAD 4
#define NUM_BIN 30
#define DISTANCE_PRECISION 10
#define LMEASURETABLE 0
#define MYINTMAX 10240000

SPDAnalysisModel::SPDAnalysisModel()
{
	DataTable = vtkSmartPointer<vtkTable>::New();
	headers.push_back("node1");
	headers.push_back("node2");
	headers.push_back("weight");
	filename = "";
	this->bProgression = false;
	disCor = 0;
	ContrastDataTable = NULL;
	m_kNeighbor =0;
}

SPDAnalysisModel::~SPDAnalysisModel()
{
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

bool SPDAnalysisModel::ReadCellTraceFile(std::string fileName, bool bContrast)
{
	std::ifstream file(fileName.c_str(), std::ifstream::in);
	int line = LineNum(fileName.c_str());
	if( !file.is_open() || line == -1)
	{
		return false;
	}

	vtkSmartPointer<vtkTable> table = ftk::LoadTable(fileName);
	if( table != NULL)
	{
		if( bContrast)
		{
			this->ContrastDataTable = table;
			ParseTraceFile(this->ContrastDataTable, bContrast);
		}
		else
		{
			this->DataTable = table;
			ParseTraceFile(this->DataTable, bContrast);
		}
	}
	else
	{
		return false;
	}	
	return true;
}

bool SPDAnalysisModel::ReadRawData(std::string fileName)
{
	vtkSmartPointer<vtkTable> table = ftk::LoadTable(fileName);
	if( table != NULL)
	{
		DataTable = table;
		std::cout<< table->GetNumberOfRows()<<"\t"<<table->GetNumberOfColumns()<<std::endl;
		DataMatrix.set_size(table->GetNumberOfRows(), table->GetNumberOfColumns());
		for( int i = 0; i < table->GetNumberOfRows(); i++)
		{
			for( int j = 0; j < table->GetNumberOfColumns(); j++)
			{
				double var = table->GetValue(i, j).ToDouble();
				if( !boost::math::isnan(var))
				{
					DataMatrix( i, j) = table->GetValue(i, j).ToDouble();
				}
				else
				{
					DataMatrix( i, j) = 0;
				}
			}
		}

		maxVertexId = DataMatrix.rows() - 1;
		indMapFromIndToVertex.resize(DataMatrix.rows());
		for( int i = 0; i < DataMatrix.rows(); i++)
		{
			indMapFromIndToVertex[i] = i;
		}

		for( int i = 0; i < CellCluster.size(); i++)
		{
			this->CellCluster[i].clear();
		}
		this->CellCluster.clear();

		this->CellClusterIndex.set_size(this->DataMatrix.rows());
		combinedCellClusterIndex.clear();
		for( int i = 0; i < this->DataMatrix.rows(); i++)
		{
			this->CellClusterIndex[i] = i;
		}
		combinedCellClusterIndex = CellClusterIndex;

		indMapFromVertexToClus.clear();
		combinedIndexMapping.clear();
		for( int i = 0; i < CellClusterIndex.size(); i++)
		{
			indMapFromVertexToClus.insert( std::pair<int, int>(indMapFromIndToVertex[i], CellClusterIndex[i]));
		}
		combinedIndexMapping = indMapFromVertexToClus;

		UNMatrixAfterCellCluster = DataMatrix;
		MatrixAfterCellCluster = UNMatrixAfterCellCluster;
		NormalizeData(MatrixAfterCellCluster);
		return true;
	}
	else
	{
		return false;
	}	
}
		//int rowIndex = 0;
		//std::string feature;
		//std::vector<std::string> rowValue;
		//bool bfirst = true;
		//while( !file.eof())
		//{
		//	getline(file, feature);
		//	if( feature.length() > 3)
		//	{
		//		split(feature, '\t', &rowValue);
		//		std::vector<std::string>::iterator iter = rowValue.begin();
		//		if ( bfirst)
		//		{
		//			this->DataMatrix.set_size( line, rowValue.size() - 1);
		//			bfirst = false;
		//		}
		//		for( int colIndex = 0; colIndex < rowValue.size() - 1; iter++, colIndex++)
		//		{
		//			this->DataMatrix( rowIndex, colIndex) = atof( (*iter).c_str());   
		//		}
                //		rowIndex++;
		//	}
		//	rowValue.clear();
		//	feature.clear();
		//}

		//// build the index and its mapping for the test data
		//for( unsigned int k = 0; k <= this->DataMatrix.cols(); k++)
		//{
		//	vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		//	this->DataTable->AddColumn(column);
		//}

		//for( unsigned int i = 0; i < this->DataMatrix.rows(); i++)
		//{
		//	vtkSmartPointer<vtkVariantArray> row = vtkSmartPointer<vtkVariantArray>::New();

		//	for( unsigned int j = 0; j <= this->DataMatrix.cols(); j++)
		//	{ 
		//		if( j == 0)
		//		{
		//			row->InsertNextValue( vtkVariant( i));
		//		}
		//		else
		//		{
		//			row->InsertNextValue( vtkVariant(this->DataMatrix( i, j - 1)));
		//		}
		//	}
		//	this->DataTable->InsertNextRow(row);
		//}

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

bool SPDAnalysisModel::MergeMatrix( vnl_matrix<double> &firstMat, vnl_matrix<double> &secondMat, vnl_matrix<double> &combinedMat)
{
	if( firstMat.cols() == secondMat.cols())
	{
		combinedMat.set_size( firstMat.rows() + secondMat.rows(), firstMat.cols());
		unsigned int i = 0;
		unsigned int rownum = firstMat.rows();
		for( i = 0; i < rownum; i++)
		{
			combinedMat.set_row( i, firstMat.get_row(i));
		}
		for( i = rownum; i < firstMat.rows() + secondMat.rows(); i++)
		{
			combinedMat.set_row( i, secondMat.get_row(i - rownum));
		}
		return true;
	}
	else
	{
		combinedMat = firstMat;   // if two matrix's columns are not equal, then it equals first matrix
		return false;
	}
}

void SPDAnalysisModel::ParseTraceFile(vtkSmartPointer<vtkTable> table, bool bContrast)
{
	table->RemoveColumnByName("Soma_X");
	table->RemoveColumnByName("Soma_Y");
	table->RemoveColumnByName("Soma_Z");
	table->RemoveColumnByName("Trace_File");

	table->RemoveColumnByName("Soma_X_Pos");
	table->RemoveColumnByName("Soma_Y_Pos");
	table->RemoveColumnByName("Soma_Z_Pos");
	table->RemoveColumnByName("Distance_to_Device");
	//for( long int i = 0; i < table->GetNumberOfRows(); i++)
	//{
	//	long int var = table->GetValue( i, 0).ToLong();
	//	this->indMapFromIndToVertex.push_back( var);
	//}
	
	if(bContrast)
	{
		std::cout<< "Constrast table:"<<table->GetNumberOfRows()<<"\t"<< table->GetNumberOfColumns()<<endl;
		this->ContrastDataTable = table;
		ConvertTableToMatrix( ContrastDataTable, this->ContrastDataMatrix, indMapFromIndToContrastVertex, contrastDistance);
		MergeMatrix(DataMatrix, ContrastDataMatrix, UNMatrixAfterCellCluster);
		this->constrastCellClusterIndex.set_size(this->ContrastDataMatrix.rows());
		for( int i = 0; i < this->ContrastDataMatrix.rows(); i++)
		{
			constrastCellClusterIndex[i] = i;
		}

		unsigned int maxclusInd = CellClusterIndex.max_value();
		int clusterIndexSize = CellClusterIndex.size();
		combinedCellClusterIndex.set_size(DataMatrix.rows() + ContrastDataMatrix.rows());
		for( int i = 0; i < clusterIndexSize; i++)
		{
			combinedCellClusterIndex[i] = CellClusterIndex[i];
		}
		for( int i = 0; i < constrastCellClusterIndex.size(); i++)
		{
			combinedIndexMapping.insert( std::pair< int, int>(indMapFromIndToContrastVertex[i] + maxVertexId + 1, constrastCellClusterIndex[i] + maxclusInd + 1));
			combinedCellClusterIndex[ i + clusterIndexSize] = constrastCellClusterIndex[i] + maxclusInd + 1;
		}
	}
	else
	{
		std::cout<< "Data table:"<<table->GetNumberOfRows()<<"\t"<< table->GetNumberOfColumns()<<endl;
		this->DataTable = table;
		this->indMapFromIndToVertex.clear();

#if LMEASURETABLE
		ConvertTableToMatrix( this->DataTable, this->DataMatrix, this->indMapFromIndToVertex, DistanceToDevice);
		UNDistanceToDevice = DistanceToDevice;
		
		double disMean = DistanceToDevice.mean();
		disCor = sqrt((DistanceToDevice - disMean).squared_magnitude() / DistanceToDevice.size());
		std::cout<< "Distance cor: "<<disCor<<endl;
		if( disCor > 1e-6)
		{
			DistanceToDevice = (DistanceToDevice - disMean) / disCor;
		}
		else
		{
			DistanceToDevice = DistanceToDevice - disMean;
		}

#else 
		//std::cout<< "Layer Data"<<std::endl;
		ConvertTableToMatrixForValidation(this->DataTable, this->DataMatrix, this->indMapFromIndToVertex, clusNo);
#endif
		UNMatrixAfterCellCluster = this->DataMatrix;
		maxVertexId = 0;
		for( int i = this->indMapFromIndToVertex.size() - 1; i >= 0; i--)
		{
			if( this->indMapFromIndToVertex[i] > maxVertexId)
			{
				maxVertexId = this->indMapFromIndToVertex[i];
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
		combinedCellClusterIndex = CellClusterIndex;

		indMapFromVertexToClus.clear();
		combinedIndexMapping.clear();
		for( int i = 0; i < CellClusterIndex.size(); i++)
		{
			indMapFromVertexToClus.insert( std::pair<int, int>(indMapFromIndToVertex[i], CellClusterIndex[i]));
		}
		combinedIndexMapping = indMapFromVertexToClus;
	}
	MatrixAfterCellCluster = UNMatrixAfterCellCluster;
	NormalizeData(MatrixAfterCellCluster);

	//vtkSmartPointer<vtkTable> normalizedTable = vtkSmartPointer<vtkTable>::New();
	//ConvertMatrixToTable(normalizedTable, MatrixAfterCellCluster, DistanceToDevice);
	//ftk::SaveTable("NormalizedTable.txt", normalizedTable);
}

void SPDAnalysisModel::ConvertTableToMatrix(vtkSmartPointer<vtkTable> table, vnl_matrix<double> &mat, std::vector<int> &index, vnl_vector<double> &distance)
{
	distance.set_size( table->GetNumberOfRows());
	for( int i = 0; i < distance.size(); i++)
	{
		distance[i] = 0;
	}

	vtkVariantArray *distanceArray = vtkVariantArray::SafeDownCast(table->GetColumnByName("Distance to Device"));
	if( distanceArray == NULL)
	{
		distanceArray = vtkVariantArray::SafeDownCast(table->GetColumnByName("Distance_to_Device"));
	}

	vtkDoubleArray *distanceDoubleArray = vtkDoubleArray::SafeDownCast(table->GetColumnByName("Distance to Device"));
	if( distanceDoubleArray == NULL)
	{
		distanceDoubleArray = vtkDoubleArray::SafeDownCast(table->GetColumnByName("Distance_to_Device"));
	}

	if( distanceArray)
	{
		for( int i = 0; i < distanceArray->GetNumberOfValues(); i++)
		{
			distance[i] = distanceArray->GetValue(i).ToDouble();
		}
	}
	else if( distanceDoubleArray)
	{
		for( int i = 0; i < distanceDoubleArray->GetNumberOfTuples(); i++)
		{
			distance[i] = distanceDoubleArray->GetValue(i);
		}
	}

	mat.set_size( table->GetNumberOfRows(), table->GetNumberOfColumns() - 2);
	mat.fill(0);
	for( int i = 0; i < table->GetNumberOfRows(); i++)
	{
		int colIndex = 0;
		for( int j = 0; j < table->GetNumberOfColumns() - 1; j++)
		{
			if( j == 0 )
			{
				index.push_back( table->GetValue(i, j).ToInt());
			}
			else
			{
				double var = table->GetValue(i, j).ToDouble();
				if( !boost::math::isnan(var))
				{
					mat( i, colIndex++) = table->GetValue(i, j).ToDouble();
				}
				else
				{
					mat( i, colIndex++) = 0;
				}
			}
		}
	}
}

void SPDAnalysisModel::ConvertTableToMatrixForValidation(vtkSmartPointer<vtkTable> table, vnl_matrix<double> &mat, 
        std::vector<int> &index, vnl_vector<int> &clusNo)
{
	std::cout<< "Table input: "<<table->GetNumberOfRows()<<"\t"<<table->GetNumberOfColumns()<<std::endl;

	clusNo.set_size( table->GetNumberOfRows());
	for( int i = 0; i < clusNo.size(); i++)
	{
		clusNo[i] = 0;
	}

	vtkVariantArray *tmpArray = vtkVariantArray::SafeDownCast(table->GetColumnByName("Label"));
	vtkDoubleArray *tmpDoubleArray = vtkDoubleArray::SafeDownCast(table->GetColumnByName("Label"));

	if( tmpArray)
	{
		for( int i = 0; i < tmpArray->GetNumberOfValues(); i++)
		{
			clusNo[i] = tmpArray->GetValue(i).ToDouble();
		}
	}
	else if( tmpDoubleArray)
	{
		for( int i = 0; i < tmpDoubleArray->GetNumberOfTuples(); i++)
		{
			clusNo[i] = tmpDoubleArray->GetValue(i);
		}
	}

	table->RemoveColumnByName("Label");

	std::cout<< "Table input: "<<table->GetNumberOfRows()<<"\t"<<table->GetNumberOfColumns()<<std::endl;

	mat.set_size( table->GetNumberOfRows(), table->GetNumberOfColumns() - 1);
	
	for( int i = 0; i < table->GetNumberOfRows(); i++)
	{
		int colIndex = 0;
		for( int j = 0; j < table->GetNumberOfColumns(); j++)
		{
			if( j == 0 )
			{
				index.push_back( table->GetValue(i, j).ToInt());
			}
			else
			{
				double var = table->GetValue(i, j).ToDouble();
				if( !boost::math::isnan(var))
				{
					mat( i, colIndex++) = table->GetValue(i, j).ToDouble();
				}
				else
				{
					mat( i, colIndex++) = 0;
				}
			}
		}
	}
}

void SPDAnalysisModel::ConvertMatrixToTable(vtkSmartPointer<vtkTable> table, vnl_matrix<double> &mat, vnl_vector<double> &distance)
{
	std::cout<< DataTable->GetNumberOfColumns()<< "\t"<< mat.cols()<<std::endl;

	for(int i = 0; i < DataTable->GetNumberOfColumns(); i++)
	{		
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName( DataTable->GetColumn(i)->GetName());
		table->AddColumn(column);
	}
	
	bool bindex = true;
	if(mat.cols() == DataTable->GetNumberOfColumns())
	{
		bindex = false;
	}

	for( int i = 0; i < mat.rows(); i++)
	{
		vtkSmartPointer<vtkVariantArray> DataRow = vtkSmartPointer<vtkVariantArray>::New();
		if( bindex)
		{
			DataRow->InsertNextValue(i);
		}
		for( int j = 0; j < mat.cols(); j++)
		{
			DataRow->InsertNextValue( mat(i,j));
		}
		if( distance.data_block() != NULL)
		{
			DataRow->InsertNextValue(distance[i]);
		}
		table->InsertNextRow(DataRow);
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

unsigned int SPDAnalysisModel::GetContrastDataSampleNum()
{
	return this->ContrastDataMatrix.rows();
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

void SPDAnalysisModel::NormalizeData( vnl_matrix<double> &mat)
{
	vnl_vector<double> std_vec;
	vnl_vector<double> mean_vec;

	GetMatrixRowMeanStd(mat, mean_vec, std_vec);
	
	for( unsigned int i = 0; i < mat.columns(); ++i)
	{
		vnl_vector<double> temp_col = mat.get_column(i);
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
		mat.set_column(i,temp_col);
	}
}

void SPDAnalysisModel::ClusterSamples(double cor)
{
	this->filename = "C++_" + QString::number(this->DataMatrix.cols())+ "_" + QString::number(this->DataMatrix.rows()) + 
					"_" + QString::number( cor, 'g', 4) + "_";
	QString filenameCluster = this->filename + "Cellclustering.txt";
	std::ofstream ofs(filenameCluster.toStdString().c_str(), std::ofstream::app);

	vnl_matrix<double> tmp = this->DataMatrix;
	NormalizeData(tmp);   // eliminate the differences of different features
	vnl_matrix<double> tmpMat =  tmp.transpose();
	ofs << tmpMat<<endl;
	int new_cluster_num = ClusterSamples(cor, tmpMat, CellClusterIndex);
	GetClusterIndexFromVnlVector( CellCluster, CellClusterIndex);

	vnl_vector<double> distance;
	
	GetAverageModule( this->DataMatrix, DistanceToDevice, CellCluster, DataMatrixAfterCellCluster, distance);
	DistanceToDevice = distance;

	indMapFromVertexToClus.clear();
	combinedIndexMapping.clear();
	for( int i = 0; i < CellClusterIndex.size(); i++)
	{
		indMapFromVertexToClus.insert( std::pair<int, int>( indMapFromIndToVertex[i], CellClusterIndex[i]));
		combinedIndexMapping.insert( std::pair< int, int>(indMapFromIndToVertex[i], CellClusterIndex[i]));
	}

	if( this->ContrastDataMatrix.empty() == false)
	{
		tmp = this->ContrastDataMatrix;
		NormalizeData(tmp);   // eliminate the differences of different features
		vnl_matrix<double> tmpMat =  tmp.transpose();
		ClusterSamples(cor, tmpMat, constrastCellClusterIndex);
		GetClusterIndexFromVnlVector( contrastCellCluster, constrastCellClusterIndex);
		GetAverageModule( this->ContrastDataMatrix, contrastDistance, contrastCellCluster, ContrastMatrixAfterCellCluster, distance);

		// this->MatrixAfterCellCluster
		vnl_matrix<double> combinedMatrixAfterCellCluster;
		MergeMatrix(DataMatrixAfterCellCluster, ContrastMatrixAfterCellCluster, combinedMatrixAfterCellCluster);
		MatrixAfterCellCluster = combinedMatrixAfterCellCluster; 
		
		// combined vertex - clusindex mapping
		unsigned int maxclusInd = CellClusterIndex.max_value();
		int clusterIndexSize = CellClusterIndex.size();
		combinedCellClusterIndex.set_size(DataMatrix.rows() + ContrastDataMatrix.rows());
		for( int i = 0; i < clusterIndexSize; i++)
		{
			combinedCellClusterIndex[i] = CellClusterIndex[i];
		}
		for( int i = 0; i < constrastCellClusterIndex.size(); i++)
		{
			combinedIndexMapping.insert( std::pair< int, int>(indMapFromIndToContrastVertex[i] + maxVertexId + 1, constrastCellClusterIndex[i] + maxclusInd + 1));
			combinedCellClusterIndex[ i + clusterIndexSize] = constrastCellClusterIndex[i] + maxclusInd + 1;
		}
	}
	else
	{
		MatrixAfterCellCluster = DataMatrixAfterCellCluster; 
		combinedCellClusterIndex = CellClusterIndex;
	}

	this->UNMatrixAfterCellCluster = MatrixAfterCellCluster;
	NormalizeData(MatrixAfterCellCluster);
	ofs.close();
}

void SPDAnalysisModel::GetCombinedDataTable(vtkSmartPointer<vtkTable> table)
{
	if( ContrastDataTable)
	{
		if( ContrastDataTable->GetNumberOfRows() > 0)
		{
			MergeTables(DataTable, ContrastDataTable, table);
		}
	}
	else
	{
		CopyTable( DataTable, table);
	}
}

void SPDAnalysisModel::CopyTable( vtkSmartPointer<vtkTable> oriTable, vtkSmartPointer<vtkTable> targetTable)
{
	for(int i = 0; i < oriTable->GetNumberOfColumns(); i++)
	{		
		vtkSmartPointer<vtkVariantArray> column = vtkSmartPointer<vtkVariantArray>::New();
		column->SetName( oriTable->GetColumn(i)->GetName());
		targetTable->AddColumn(column);
	}
	
	for( vtkIdType i = 0; i < oriTable->GetNumberOfRows(); i++)
	{
		vtkSmartPointer<vtkVariantArray> row = oriTable->GetRow(i);
		targetTable->InsertNextRow(row);
	}
}

bool SPDAnalysisModel::MergeTables(vtkSmartPointer<vtkTable> firstTable, vtkSmartPointer<vtkTable> secondTable, vtkSmartPointer<vtkTable> table)
{
	if( firstTable->GetNumberOfColumns() != secondTable->GetNumberOfColumns())
	{
		CopyTable( table, firstTable);
		return false;
	}

	for(int i = 0; i < firstTable->GetNumberOfColumns(); i++)
	{		
		vtkSmartPointer<vtkVariantArray> column = vtkSmartPointer<vtkVariantArray>::New();
		column->SetName( firstTable->GetColumn(i)->GetName());
		table->AddColumn(column);
	}
	
	vtkIdType firstRowNum = firstTable->GetNumberOfRows();
	int maxId = 0;
	for( vtkIdType i = 0; i < firstTable->GetNumberOfRows(); i++)
	{
		table->InsertNextRow(firstTable->GetRow(i));
		int id = firstTable->GetValue( i, 0).ToInt();
		if( id > maxId)
		{
			maxId = id;
		}
	}
	for( vtkIdType i = 0; i < secondTable->GetNumberOfRows(); i++)
	{
		table->InsertNextRow(secondTable->GetRow(i));
		int id = secondTable->GetValue( i, 0).ToInt();
		table->SetValue( firstRowNum + i, 0, vtkVariant(id + maxId + 1));
	}
	return true;
}

void SPDAnalysisModel::GetClusterIndexFromVnlVector(std::vector< std::vector<int> > &clusterIndex, vnl_vector<unsigned int> &index)
{
	for( int i = 0; i < clusterIndex.size(); i++)
	{
		clusterIndex[i].clear();
	}
	clusterIndex.clear();

	/// rebuild the datamatrix for MST 
	std::vector<int> clus;
	for( int i = 0; i <= index.max_value(); i++)
	{
		clusterIndex.push_back(clus);
	}
	for( int i = 0; i < index.size(); i++)
	{
		int k = index[i];
		clusterIndex[k].push_back(i);
	}
}

void SPDAnalysisModel::GetAverageModule( vnl_matrix<double> &mat, vnl_vector<double> &distance, std::vector< std::vector<int> > &index, 
										vnl_matrix<double> &averageMat, vnl_vector<double> &averageDistance)
{
	averageMat.set_size( index.size(), mat.cols());
	averageDistance.set_size( index.size());
	for( int i = 0; i < index.size(); i++)
	{
		vnl_vector<double> tmpVec( mat.cols());
		tmpVec.fill(0);
		double tmpDis;
		tmpDis = 0;
		for( int j = 0; j < index[i].size(); j++)
		{
			tmpVec += mat.get_row(index[i][j]);
			tmpDis += distance[ index[i][j]];
		}
		tmpVec = tmpVec / CellCluster[i].size();
		tmpDis = tmpDis / CellCluster[i].size();
		averageMat.set_row( i, tmpVec);
		averageDistance[i] = tmpDis;
	}	
}

int SPDAnalysisModel::ClusterSamples( double cor, vnl_matrix<double> &mat, vnl_vector<unsigned int> &index)									
{
	vnl_matrix<double> mainMatrix = mat;
	NormalizeData(mainMatrix);
	vnl_matrix<double> moduleMean = mainMatrix;

	TreeIndex.set_size(mainMatrix.cols());
	vnl_vector<unsigned int> clusterIndex( moduleMean.cols());

	for( unsigned int i = 0; i < moduleMean.cols(); i++)
	{
		clusterIndex[i] = i;
	}

	int old_cluster_num = moduleMean.cols();
	int new_cluster_num = moduleMean.cols();
	do
	{
		old_cluster_num = new_cluster_num;

		new_cluster_num = ClusterAggFeatures( mainMatrix, clusterIndex, moduleMean, cor);
		std::cout<< "new cluster num:"<< new_cluster_num<<" vs "<<"old cluster num:"<<old_cluster_num<<endl;
	}
	while( old_cluster_num > new_cluster_num);

	index = clusterIndex; 
	mat = moduleMean;
	std::cout<< "Cell Cluster Size:"<<new_cluster_num<<endl;
	return new_cluster_num;
}

void SPDAnalysisModel::GetCellClusterSize( std::vector<int> &clusterSize)
{
	if( CellCluster.size() > 0)
	{
		for( int i = 0; i < CellCluster.size(); i++)
		{
			int num = CellCluster[i].size();
			clusterSize.push_back( num);
		}
	}
	else
	{
		for( int i = 0; i < CellClusterIndex.size(); i++)
		{
			int num = 1;
			clusterSize.push_back( num);
		}
	}
}

vtkSmartPointer<vtkTable> SPDAnalysisModel::GetDataTableAfterCellCluster()
{
	vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
	ConvertMatrixToTable( table, MatrixAfterCellCluster, DistanceToDevice);
	return table;
}

int SPDAnalysisModel::ClusterAgglomerate(double cor, double mer)
{
	vnl_vector<unsigned int> clusterIndex;
	clusterIndex.set_size( MatrixAfterCellCluster.cols());
	vnl_matrix<double> moduleMean = MatrixAfterCellCluster;
	moduleMean.normalize_columns();

	this->filename = "C++_" + QString::number(MatrixAfterCellCluster.cols())+ "_" + QString::number(MatrixAfterCellCluster.rows()) + 
					"_" + QString::number( cor, 'g', 4) + "_" + QString::number( mer, 'g', 4)+ "_";
	QString filenameCluster = this->filename + "clustering.txt";
	std::ofstream ofs(filenameCluster.toStdString().c_str(), std::ofstream::out);

	//NormalizeData( this->DataMatrix);

	TreeIndex.set_size(this->MatrixAfterCellCluster.cols());
	for( unsigned int i = 0; i < this->MatrixAfterCellCluster.cols(); i++)
	{
		clusterIndex[i] = i;
		TreeIndex[i] = i;
	}
	TreeData.clear();

	int old_cluster_num = MatrixAfterCellCluster.cols();
	int new_cluster_num = MatrixAfterCellCluster.cols();
	std::cout<< "Start Clustering:" <<std::endl;
	ofs<<"Left cluster number:"<<endl;
	do
	{
		old_cluster_num = new_cluster_num;
		//ofs<<old_cluster_num<<endl;
		new_cluster_num = ClusterAggFeatures( MatrixAfterCellCluster, clusterIndex, moduleMean, cor);
		std::cout<< new_cluster_num<<"\t";
	}
	while( old_cluster_num > new_cluster_num);
	//ofs<<endl<<endl;
	std::cout<<std::endl;

	this->ClusterIndex = clusterIndex;
	this->ModuleMean = moduleMean;

	ClusterMerge(cor, mer);

	ofs << this->ClusterIndex <<std::endl;

	for( unsigned int i = 0; i <= this->ClusterIndex.max_value(); i++)
	{
		ofs<<"Cluster:"<< i<<endl;
		for( unsigned int j = 0; j < this->ClusterIndex.size(); j++)
		{
			if( ClusterIndex[j] == i)
			{
				char* name = this->DataTable->GetColumn(j + 1)->GetName();
				std::string featureName(name);
				ofs<<j<<"\t"<<featureName<<endl;
			}
		}
	}

	ofs.close();
	std::cout<<"The cluster size after cluster: "<<this->ClusterIndex.max_value() + 1<<endl;
	//std::cout<<clusterIndex<<endl;

	return new_cluster_num;
}


int SPDAnalysisModel::ClusterAggFeatures(vnl_matrix<double>& mainmatrix, vnl_vector<unsigned int>& index, vnl_matrix<double>& mean, double cor)
{
	vnl_vector<int> moduleSize = GetModuleSize(index);
	vnl_vector<int> isActiveModule;
	vnl_vector<double> moduleCenter;
	vnl_vector<unsigned int> indexmap(moduleSize.size());
	for( unsigned int i = 0; i < indexmap.size(); i++)
	{
		indexmap[i] = i;
	}
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
	vnl_vector<int> tmp(moduleSize.size());
	unsigned int activeModuleNum = isActiveModule.size();
	while( activeModuleNum > 1)
	{
		//std::cout<< activeModuleNum<<"\t";
		if( activeModuleNum % 1000 == 0)
		{
			std::cout<< activeModuleNum << "\t/"<<isActiveModule.size()<<"\n";
		}
		vnl_matrix<double> newModule;
		vnl_vector<double> newModuleMeans;
		vnl_vector<double> newModuleStd;
		vnl_vector<double> newModuleCor;
		vnl_vector<double> moduleCor; 

		unsigned int moduleId = moduleSize.arg_min();
		moduleCenter = mean.get_column(moduleId);
		moduleCor.set_size( mean.cols());
		moduleCor = moduleCenter * mean;
		moduleCor[moduleId] = 0;

		unsigned int moduleToDeleteId = moduleCor.arg_max();
		int newModuleSize = moduleSize[moduleId] + moduleSize[moduleToDeleteId];

		newModule.set_size(mainmatrix.rows(), newModuleSize);

		if( moduleToDeleteId != moduleId && isActiveModule[moduleToDeleteId] == 1)
		{
			GetCombinedMatrix( mainmatrix, index, moduleId, moduleToDeleteId, newModule, &indexmap);   
			newModule.normalize_columns();

			vnl_matrix<double> mat = newModule.transpose();
			GetMatrixRowMeanStd(mat, newModuleMeans, newModuleStd);

			newModuleMeans = newModuleMeans.normalize();
			newModuleCor = newModuleMeans * newModule;
			double newCor = newModuleCor.mean();                       

			if( newCor > cor)
			{
				isActiveModule[moduleToDeleteId] = 0;
				activeModuleNum--;
				moduleSize[moduleId] += moduleSize[moduleToDeleteId];
				moduleSize[moduleToDeleteId] = MYINTMAX;

				mean.set_column(moduleId, newModuleMeans);
				mean.set_column(moduleToDeleteId, zeroCol);
				
				delColsVec.push_back(moduleToDeleteId);
				
				indexmap[moduleToDeleteId] = moduleId;

				Tree tr( TreeIndex[moduleId], TreeIndex[moduleToDeleteId], newCor, newIndex);
				TreeIndex[moduleId] = newIndex;
				TreeIndex[moduleToDeleteId] = -1;
				TreeData.push_back( tr);
				newIndex++;
			}
		}

		isActiveModule[moduleId] = 0;
		activeModuleNum--;
		moduleSize[moduleId] = MYINTMAX;
	}

	StandardizeIndex(index, indexmap);
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
										 vnl_matrix<double>& mat, vnl_vector<unsigned int>* indexMap)
{
	unsigned int ind = 0;
	if(indexMap == NULL)
	{
		for(unsigned int i = 0; i < index.size(); i++)
		{
			if( index[i] == moduleId || index[i]== moduleDeleteId)
			{

				mat.set_column( ind++, datamat.get_column(i));
			}
		}
	}
	else
	{
		for(unsigned int i = 0; i < index.size(); i++)
		{
			unsigned int k = index[i];
			if( (*indexMap)[k] == moduleId || (*indexMap)[k] == moduleDeleteId)
			{

				mat.set_column( ind++, datamat.get_column(i));
			}
		}
	}
}

vnl_vector<int> SPDAnalysisModel::GetModuleSize(vnl_vector<unsigned int>& index)   // index starts from 0
{
	vnl_vector<int> moduleSize;
	moduleSize.set_size(index.max_value() + 1);
	moduleSize.fill(0);
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

void SPDAnalysisModel::StandardizeIndex(vnl_vector<unsigned int>& index, vnl_vector<unsigned int> &indexmap)
{
	for(unsigned int i = 0; i < index.size(); i++)
	{
		unsigned int k = index[i];
		index[i] = indexmap[k];
	}

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

		vnl_matrix<double> dataTmp( MatrixAfterCellCluster.rows(), moduleSize[i] + moduleSize[j]);
		GetCombinedMatrix( MatrixAfterCellCluster, this->ClusterIndex, i, j, dataTmp);
		dataTmp.normalize_columns();

		vnl_vector<double> dataMean;
		vnl_vector<double> dataSd;
		vnl_matrix<double> mat = dataTmp.transpose();    // bugs lying here
		GetMatrixRowMeanStd( mat, dataMean, dataSd);
		dataMean = dataMean.normalize();

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
	//std::cout<<this->ClusterIndex<<endl<<endl;

	//QString filenameCluster = this->filename + "clustering.txt";
	//std::ofstream ofs( filenameCluster.toStdString().c_str(), std::ofstream::out | ios::app);
	//ofs<< "clustering result after merge:"<<endl;
	//for( unsigned int i = 0; i <= this->ClusterIndex.max_value(); i++)
	//{
	//	ofs<<"Cluster: "<<i<<endl;
	//	for( unsigned int j = 0; j < this->ClusterIndex.size(); j++)
	//	{
	//		if( ClusterIndex[j] == i)
	//		{
	//			char* name = this->DataTable->GetColumn(j + 1)->GetName();
	//			std::string featureName(name);
	//			ofs<< j<<"\t"<<featureName<<endl;
	//		}
	//	}
	//}

	//ofs<<endl<<"TreeData:"<<endl;
	//for( int i = 0; i < TreeData.size(); i++)
	//{
	//	Tree tr = TreeData[i];
	//	ofs<< tr.first<< "\t"<< tr.second<< "\t"<< tr.cor<< "\t"<< tr.parent<<endl;
	//}
	//ofs<< TreeIndex<<endl;
	//ofs.close();
}

void SPDAnalysisModel::GetFeatureIdbyModId(std::vector<unsigned int> &modID, std::vector<unsigned int> &featureID)
{
	for( int i = 0; i < ClusterIndex.size(); i++)
	{
		if( IsExist( modID, ClusterIndex[i]))
		{
			featureID.push_back(i);
		}
	}
}

double SPDAnalysisModel::GetANOVA(std::vector< std::vector< long int> > &index, std::vector< unsigned int> &selFeatureId)
{
	if( MatrixAfterCellCluster.rows() == index.size() || index.size() == 1)
	{
		std::cout<< "No clusters. "<<std::endl;
		return 0;
	}
	vnl_matrix<double> subMatFeaSel;
	vnl_vector<double> overallMean;
	vnl_matrix<double> clusMean;

	GetSubFeatureMatrix( MatrixAfterCellCluster, selFeatureId, subMatFeaSel, overallMean);
	clusMean.set_size(index.size(), selFeatureId.size());

	double MSB = 0;
	double MSW = 0;
	double K = 0;
	double n = MatrixAfterCellCluster.rows();
	double m = index.size();

	for( int i = 0; i < index.size(); i++)
	{
		std::vector<int> cindex(index[i].size());
		for( int j = 0; j < index[i].size(); j++)
		{
			int id = (int)index[i][j];
			int cId = combinedIndexMapping.find(id)->second;
			cindex[j]= cId;
		}
		vnl_vector<double> imean;
		vnl_matrix<double> imat;
		GetSubSampleMatrix(subMatFeaSel, cindex, imat, imean);
		for( int j = 0; j < cindex.size(); j++)
		{
			vnl_vector<double> tmp = subMatFeaSel.get_row(cindex[j]) - imean;
			MSW += 1.0 / ( n - m) * tmp.squared_magnitude();
		}
		vnl_vector<double> tmp = imean - overallMean;
		MSB += ((double)index[i].size() /( m - 1)) * tmp.squared_magnitude();
		K = - 1.0 / (m - 1) * (double)index[i].size() * index[i].size() / n;
	}
	K += n / ( m - 1);

	double anova = (MSB - MSW) / (MSB + (K - 1) * MSW);
	return anova;
}

void SPDAnalysisModel::GetSubFeatureMatrix(vnl_matrix<double> &mat, std::vector< unsigned int> &featureId, vnl_matrix<double> &subMat, vnl_vector<double> &mean)
{
	subMat.set_size( mat.rows(), featureId.size());
	mean.set_size(featureId.size());
	mean.fill(0);

	for(unsigned int i = 0; i < featureId.size(); i++)
	{
		unsigned int ind = featureId[i];
		vnl_vector<double> colVal = mat.get_column(ind);
		subMat.set_column( (unsigned int)i, colVal);
	}

	for( unsigned int i = 0; i < subMat.rows(); i++)
	{
		vnl_vector<double> rowVal = subMat.get_row(i);
		mean = (mean * i + rowVal) / (i + 1);
	}
}

void SPDAnalysisModel::GetSubSampleMatrix(vnl_matrix<double> &mat, std::vector< int> &sampleId, vnl_matrix<double> &subMat, vnl_vector<double> &mean)
{
	subMat.set_size(sampleId.size(), mat.cols());
	mean.set_size(mat.cols());
	mean.fill(0);

	for( unsigned int i = 0; i < sampleId.size(); i++)
	{
		vnl_vector<double> rowVal = mat.get_row(sampleId[i]);
		subMat.set_row(i, rowVal);
		mean = (mean * i + rowVal) / (i + 1);
	}
}

void SPDAnalysisModel::ConvertClusIndexToSampleIndex(std::vector< std::vector< long int> > &clusIndex, std::vector< std::vector< long int> > &sampleIndex)
{
	sampleIndex.resize(clusIndex.size());
	for( int i = 0; i < clusIndex.size(); i++)
	{
		for( int j = 0; j < clusIndex[i].size(); j++)
		{
			long int id = clusIndex[i][j];
			sampleIndex[i].push_back((long int)indMapFromIndToVertex[id]);
		}
	}
}

void SPDAnalysisModel::GetSingleLinkageClusterAverage(std::vector< std::vector< long int> > &index, vnl_matrix<double> &clusAverageMat)  // after single linkage clustering
{
	clusAverageMat.set_size( index.size(), MatrixAfterCellCluster.cols());
	for( int i = 0; i < index.size(); i++)
	{
		vnl_vector<double> tmp( MatrixAfterCellCluster.cols());
		tmp.fill(0);
		for( int j = 0; j < index[i].size(); j++)
		{
			long int id = index[i][j];
			int clusId = combinedIndexMapping.find(id)->second;
			vnl_vector<double> row = MatrixAfterCellCluster.get_row(clusId);
			tmp = tmp + row;
		}
		tmp = tmp / index[i].size();
		clusAverageMat.set_row(i, tmp);
	}
}

void SPDAnalysisModel::GetDataMatrix( vnl_matrix<double> &mat)
{
	mat = MatrixAfterCellCluster;
}

void SPDAnalysisModel::GetArrangedMatrixByConnectedComponent(std::vector< std::vector< long int> > &index, vnl_matrix<double> &mat)
{
	mat.set_size(MatrixAfterCellCluster.rows(), MatrixAfterCellCluster.cols());
	unsigned int count = 0;
	for( int i = 0; i < index.size(); i++)
	{
		for( int j = 0; j < index[i].size(); j++)
		{
			mat.set_row( count++,MatrixAfterCellCluster.get_row(index[i][j]));
		}
	}
}

void SPDAnalysisModel::GetClusterMapping( std::map< int, int> &index)
{
	index = combinedIndexMapping;
}

void SPDAnalysisModel::GetSingleLinkageClusterMapping(std::vector< std::vector< long int> > &index, std::vector<int> &newIndex)
{
	//std::ofstream ofs("clusIndexMapping.txt");
	//ofs<< "Cell cluster index before single linkage:"<<endl;
	//ofs<< CellClusterIndex<<endl;

	std::vector< unsigned int> tmpInd;
	tmpInd.resize( combinedCellClusterIndex.max_value() + 1);
	for( unsigned int i = 0; i < index.size(); i++)
	{
		for( unsigned int j = 0; j < index[i].size(); j++)
		{
			unsigned int id = index[i][j];
			int clusId = combinedIndexMapping.find(id)->second;
			tmpInd[clusId] = i;
		}
	}
	
	newIndex.resize( combinedCellClusterIndex.size());
	for( unsigned int i = 0; i < combinedCellClusterIndex.size(); i++)
	{
		unsigned int ind = combinedCellClusterIndex[i];
		newIndex[i] = tmpInd[ind];
	}

	//ofs<< "Cell cluster index after single linkage:"<<endl;
	//for( int i = 0; i < newIndex.size(); i++)
	//{
	//	ofs<< newIndex[i]<<"\t";
	//}
	//ofs.close();
}

void SPDAnalysisModel::HierachicalClustering()    
{
	assert(TreeIndex.size() == this->ClusterIndex.max_value() + 1);
	PublicTreeData = TreeData;

	vnl_vector<int> moduleSize = GetModuleSize(this->ClusterIndex);
	vnl_matrix<double> matrix = MatrixAfterCellCluster;
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
			newModule.set_size(this->MatrixAfterCellCluster.rows(), newModuleSize);
			GetCombinedMatrix(this->MatrixAfterCellCluster, this->ClusterIndex, i, j, newModule);   
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
				vnl_matrix<double> newModule(MatrixAfterCellCluster.rows(), moduleSize(min) + moduleSize(k));
				vnl_vector<double> newModuleMeans;
				vnl_vector<double> newModuleStd;
				vnl_vector<double> newModuleCor;

				GetCombinedMatrix( this->MatrixAfterCellCluster, cIndex, min, k, newModule);   
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

bool SPDAnalysisModel::SetProgressionType(bool bProg)
{
	if( disCor != 0)
	{
		this->bProgression = bProg;
		return true;
	}
	else
	{
		return false;
	}
}

bool SPDAnalysisModel::GetProgressionType()
{
	return this->bProgression;
}

void SPDAnalysisModel::GenerateMST()
{
	this->ModuleMST.clear();
	this->MSTWeight.clear();
	this->ModuleGraph.clear();

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

	#pragma omp parallel for
	for( int i = 0; i <= this->ClusterIndex.max_value(); i++)
	{
		std::cout<<"Building MST for module " << i<<endl;

		//ofsmat<<"Module:"<< i + 1<<endl;
		//ofsmat<< clusterMat<<endl<<endl;
		vnl_matrix<double> clusterMat = clusterMatList[i];

		Graph graph(num_nodes);

		for( unsigned int k = 0; k < num_nodes; k++)
		{
			for( unsigned int j = k + 1; j < num_nodes; j++)
			{
				double dist = EuclideanBlockDist( clusterMat, k, j);   // need to be modified 
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

		this->ModuleMST[i] = vertex;
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
				double dist = EuclideanBlockDist( clusterMat, i, this->ModuleMST[n][i]);
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

bool SPDAnalysisModel::IsConnected(std::multimap<int, int> &neighborGraph, std::vector<long int> &index1, std::vector<long int> &index2)
{
    for( int i = 0; i < index1.size(); i++)
    {
        int ind1 = index1[i];
        for( int j = 0; j < index2.size(); j++)
        {
            int ind2 = index2[j];
            unsigned int count = neighborGraph.count(ind1);
            if( count > 0)
            {
                std::multimap<int, int>::iterator iter = neighborGraph.find(ind1);
                for( unsigned int k = 0; k < count; k++)
                {
                    if(iter->second == ind2)
                    {
                        return true;
                    }
                    else
                    {
                        iter++;
                    }
                }
            }
        }
    }
    return false;
}

// clusterNum: for connected component
vtkSmartPointer<vtkTable> SPDAnalysisModel::GenerateMST( vnl_matrix<double> &mat, std::vector< unsigned int> &selFeatures, std::vector<int> &clusterNum)
{
	std::vector<int> nstartRow;
	nstartRow.resize(clusterNum.size());
	for( int i = 0; i < clusterNum.size(); i++)
	{
		nstartRow[i] = 0;
		for( int j = 0; j < i; j++)
		{
			nstartRow[i] += clusterNum[j];
		}
	}

	vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
	for(int i = 0; i < this->headers.size(); i++)
	{		
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName( (this->headers)[i].c_str());
		table->AddColumn(column);
	}

	for( int i = 0; i < clusterNum.size(); i++)
	{
		vnl_matrix<double> clusterMat;
		GetCombinedMatrix( mat, nstartRow[i], clusterNum[i], selFeatures, clusterMat);

		int num_nodes = clusterMat.rows();
		Graph graph(num_nodes);

		for( unsigned int k = 0; k < num_nodes; k++)
		{
			for( unsigned int j = k + 1; j < num_nodes; j++)
			{
				double dist = EuclideanBlockDist( clusterMat, k, j);   // need to be modified 
				//double dist = EuclideanBlockDist( clusterMat, k, j);
				boost::add_edge(k, j, dist, graph);
			}
		}

		std::vector< boost::graph_traits< Graph>::vertex_descriptor> vertex( num_nodes);

		try
		{
			boost::prim_minimum_spanning_tree(graph, &vertex[0]);
		}
		catch(...)
		{
			std::cout<< "MST construction failure!"<<endl;
			exit(111);
		}

		for( int k = 0; k < vertex.size(); k++)
		{
			if( k != vertex[k])
			{
				double dist = EuclideanBlockDist( clusterMat, k, vertex[k]);
				vtkSmartPointer<vtkVariantArray> DataRow = vtkSmartPointer<vtkVariantArray>::New();
				DataRow->InsertNextValue(k + nstartRow[i]);
				DataRow->InsertNextValue( vertex[k] + nstartRow[i]);
				DataRow->InsertNextValue(dist);
				table->InsertNextRow(DataRow);
			}
		}
	}

	/// build MST for the modules based on their center distance
	if(clusterNum.size() > 1)
	{
		std::cout<< "Build MST for modules"<<std::endl;
		Graph globe_graph(clusterNum.size());
		vnl_matrix< double> dist(clusterNum.size(), clusterNum.size());
		vnl_matrix< int> rowi(clusterNum.size(), clusterNum.size());
		vnl_matrix< int> rowj(clusterNum.size(), clusterNum.size());
		for( int i = 0; i < clusterNum.size(); i++)
		{
			vnl_matrix<double> clusterMati;
			GetCombinedMatrix( mat, nstartRow[i], clusterNum[i], selFeatures, clusterMati);
			for( int j = i + 1; j < clusterNum.size(); j++)
			{
				vnl_matrix<double> clusterMatj;
				GetCombinedMatrix( mat, nstartRow[j], clusterNum[j], selFeatures, clusterMatj);
				dist(i,j) = ComputeModuleDistanceAndConnection(clusterMati, clusterMatj, rowi(i,j), rowj(i,j));
				dist(j,i) = dist(i,j);
				rowi(j,i) = rowi(i,j);
				rowj(j,i) = rowj(i,j);
				boost::add_edge(i, j, dist(i,j), globe_graph);
			}
		}

		std::vector< boost::graph_traits< Graph>::vertex_descriptor> global_vertex( clusterNum.size());
		try
		{
			boost::prim_minimum_spanning_tree(globe_graph, &global_vertex[0]);
		}
		catch(...)
		{
			std::cout<< "MST construction failure!"<<endl;
			exit(111);
		}
		for( int i = 0; i < global_vertex.size(); i++)
		{
			if( i != global_vertex[i])
			{
				int j = global_vertex[i];
				int min = i < j ? i : j;
				int max = i > j ? i : j;

				int modulei = rowi(min, max);
				int modulej = rowj(min, max);
				
				vtkSmartPointer<vtkVariantArray> DataRow = vtkSmartPointer<vtkVariantArray>::New();
				DataRow->InsertNextValue( modulei + nstartRow[min]);
				DataRow->InsertNextValue( modulej + nstartRow[max]);
				DataRow->InsertNextValue(dist(i,j));
				table->InsertNextRow(DataRow);
			}
		}
	}
	
	//ftk::SaveTable("MST.txt", table);
	return table;
}

void SPDAnalysisModel::WriteGraphToGDF( std::vector< unsigned int> &selFeatures)
{
	std::ofstream ofs("test.gdf");
	ofs<<"nodedef> name\n";
	int num_nodes = MatrixAfterCellCluster.rows();
	for( int i = 0; i < num_nodes; i++)
	{
		ofs<<i<<std::endl;
	}
	ofs<<"edgedef> node1,node2,weight DOUBLE\n";

	vnl_matrix<double> clusterAllMat;
    GetCombinedMatrix( MatrixAfterCellCluster, 0, num_nodes, selFeatures, clusterAllMat);
	std::vector<unsigned int> nearIndex;
	FindNearestKSample(clusterAllMat, nearIndex, m_kNeighbor);
	vnl_matrix<unsigned char> kNearMat(num_nodes, num_nodes);
	kNearMat.fill(0);
	for( int ind = 0; ind < num_nodes; ind++)
	{
		for( int k = 0; k < m_kNeighbor; k++)
		{
			int neighborIndex = ind * m_kNeighbor + k;
			neighborIndex = nearIndex[neighborIndex];
			if( kNearMat(ind, neighborIndex) == 0)
			{
				kNearMat(ind, neighborIndex) = 1;
				kNearMat(neighborIndex, ind) = 1;
				double dist = EuclideanBlockDist( clusterAllMat, ind, neighborIndex);
				ofs<<ind<<","<<neighborIndex<<","<<dist<<"\n";
			}
		}
	}
	ofs.close();
}

vtkSmartPointer<vtkTable> SPDAnalysisModel::GenerateSubGraph( vnl_matrix<double> &mat, std::vector< std::vector< long int> > &clusIndex, std::vector< unsigned int> &selFeatures, std::vector<int> &clusterNum)
{
	std::vector<int> nstartRow;
	nstartRow.resize(clusterNum.size());
	for( int i = 0; i < clusterNum.size(); i++)
	{
		nstartRow[i] = 0;
		for( int j = 0; j < i; j++)
		{
			nstartRow[i] += clusterNum[j];
		}
	}

	vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
	for(int i = 0; i < this->headers.size(); i++)
	{		
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName( (this->headers)[i].c_str());
		table->AddColumn(column);
	}

	int num_nodes = MatrixAfterCellCluster.rows();
	vnl_matrix<unsigned char> kNearMat(num_nodes, num_nodes);
	kNearMat.fill(0);
	vnl_matrix<double> clusterAllMat;
    GetCombinedMatrix( MatrixAfterCellCluster, 0, num_nodes, selFeatures, clusterAllMat);
	std::vector<unsigned int> nearIndex;
	FindNearestKSample(clusterAllMat, nearIndex, m_kNeighbor);
	std::multimap<int, int> nearNeighborGraph;
	for( int ind = 0; ind < num_nodes; ind++)
	{
		for( int k = 0; k < m_kNeighbor; k++)
		{
			int neighborIndex = ind * m_kNeighbor + k;
			neighborIndex = nearIndex[neighborIndex];
			if( kNearMat(ind, neighborIndex) == 0)
			{
				kNearMat(ind, neighborIndex) = 1;
				kNearMat(neighborIndex, ind) = 1;
				nearNeighborGraph.insert( std::pair<int, int>(ind, neighborIndex));
				nearNeighborGraph.insert( std::pair<int, int>(neighborIndex, ind));
			}
		}
	}

	for( int i = 0; i < clusterNum.size(); i++)
	{ 
        vnl_matrix<double> clusterMat;
        GetCombinedMatrix( mat, nstartRow[i], clusterNum[i], selFeatures, clusterMat);
        for( int k = nstartRow[i]; k < nstartRow[i] + clusterNum[i]; k++)
        {
            std::vector<long int> tmpVeck = clusIndex[k];
            for( int j = k + 1; j < nstartRow[i] + clusterNum[i]; j++)
			{
                std::vector<long int> tmpVecj = clusIndex[j];

                if(IsConnected(nearNeighborGraph, tmpVeck, tmpVecj))
                {
                    double dist = EuclideanBlockDist( clusterMat, k - nstartRow[i], j - nstartRow[i]);
                    vtkSmartPointer<vtkVariantArray> DataRow = vtkSmartPointer<vtkVariantArray>::New();
                    DataRow->InsertNextValue(k);
                    DataRow->InsertNextValue(j);
                    DataRow->InsertNextValue(dist);
                    table->InsertNextRow(DataRow);
                }
			}
        }
	}

	//ftk::SaveTable("NearestGraph1.txt", table);

	/// build MST for the modules based on their center distance
	Graph globe_graph(clusterNum.size());
	vnl_matrix< double> dist(clusterNum.size(), clusterNum.size());
	vnl_matrix< int> rowi(clusterNum.size(), clusterNum.size());
	vnl_matrix< int> rowj(clusterNum.size(), clusterNum.size());
	for( int i = 0; i < clusterNum.size(); i++)
	{
		vnl_matrix<double> clusterMati;
		GetCombinedMatrix( mat, nstartRow[i], clusterNum[i], selFeatures, clusterMati);
		for( int j = i + 1; j < clusterNum.size(); j++)
		{
			vnl_matrix<double> clusterMatj;
			GetCombinedMatrix( mat, nstartRow[j], clusterNum[j], selFeatures, clusterMatj);
			dist(i,j) = ComputeModuleDistanceAndConnection(clusterMati, clusterMatj, rowi(i,j), rowj(i,j));
			dist(j,i) = dist(i,j);
			rowi(j,i) = rowi(i,j);
			rowj(j,i) = rowj(i,j);
			boost::add_edge(i, j, dist(i,j), globe_graph);
		}
	}

	std::vector< boost::graph_traits< Graph>::vertex_descriptor> global_vertex( clusterNum.size());
	try
	{
		boost::prim_minimum_spanning_tree(globe_graph, &global_vertex[0]);
	}
	catch(...)
	{
		std::cout<< "MST construction failure!"<<endl;
		exit(111);
	}
	for( int i = 0; i < global_vertex.size(); i++)
	{
		if( i != global_vertex[i])
		{
			int j = global_vertex[i];
			int min = i < j ? i : j;
			int max = i > j ? i : j;

			int modulei = rowi(min, max);
			int modulej = rowj(min, max);
			
			vtkSmartPointer<vtkVariantArray> DataRow = vtkSmartPointer<vtkVariantArray>::New();
			DataRow->InsertNextValue( modulei + nstartRow[min]);
			DataRow->InsertNextValue( modulej + nstartRow[max]);
			DataRow->InsertNextValue(dist(i,j));
			table->InsertNextRow(DataRow);
		}
	}
	
    //ftk::SaveTable("NearestGraph.txt", table);
	return table;
}

// rowi and rowj are the nearest between two module, distance is the distance between two centers
double SPDAnalysisModel::ComputeModuleDistanceAndConnection(vnl_matrix<double> &mati, vnl_matrix<double> &matj, int &rowi, int &rowj)
{
	vnl_vector<double> averageVeci;
	GetAverageVec(mati, averageVeci);
	vnl_vector<double> averageVecj;
	GetAverageVec(matj, averageVecj);
	double distance = EuclideanBlockDist(averageVeci, averageVecj);
	double min_dist = 1e9;
	for( unsigned int i = 0; i < mati.rows(); i++)
	{
		vnl_vector<double> veci = mati.get_row(i);
		for( unsigned int j = 0; j < matj.rows(); j++)
		{
			vnl_vector<double> vecj = matj.get_row(j);
			double dist = EuclideanBlockDist(veci, vecj);
			if( dist < min_dist)
			{
				rowi = i;
				rowj = j;
				min_dist = dist;
			}
		}
	}
	return distance;
}

void SPDAnalysisModel::GetAverageVec(vnl_matrix<double> &mat, vnl_vector<double> &vec)
{
	vec.set_size( mat.cols());
	vec.fill(0);
	for( unsigned int i = 0; i < mat.rows(); i++)
	{
		vec += mat.get_row(i);
	}
	vec = vec / mat.rows();
}

double SPDAnalysisModel::EuclideanBlockDist( vnl_vector<double>& vec1, vnl_vector<double>& vec2)
{
	vnl_vector<double> subm = vec1 - vec2;
	double dis = subm.magnitude();
	return dis;
}

void SPDAnalysisModel::GetCombinedMatrix( vnl_matrix<double> &datamat, int nstart, int nrow, std::vector< unsigned int> selFeatureIDs, vnl_matrix<double>& mat)
{
	int count = 0;
	mat.set_size( nrow, selFeatureIDs.size());
	for( unsigned int j = 0; j < nrow; j++)
	{
		count = 0;
		for( unsigned int i = 0; i < datamat.cols(); i++)
		{
			if( IsExist(selFeatureIDs, i))
			{
				mat(j, count) = datamat( nstart + j, i);
				count++;
			}
		}
	}
}

void SPDAnalysisModel::GetCombinedMatrixByModuleId( vnl_matrix<double> &datamat, std::vector< std::vector< unsigned int> > &featureClusterIndex, std::vector< unsigned int> &selModuleIDs, vnl_matrix<double>& mat)
{
	mat.set_size(datamat.rows(), selModuleIDs.size());
	unsigned int count = 0;
	for( unsigned int i = 0; i < selModuleIDs.size(); i++)
	{
		std::vector<unsigned int> vec = featureClusterIndex[ selModuleIDs[i]];
		vnl_vector<double> modMean(datamat.rows());
		modMean.fill(0);
		for( unsigned int j = 0; j < vec.size(); j++)
		{
			vnl_vector<double> colVec = datamat.get_column(vec[j]);
			modMean += colVec;
		}
		modMean = modMean / vec.size();
		mat.set_column( i, modMean);
	}
}

void SPDAnalysisModel::GenerateDistanceMST()
{
	if(disCor <= 10)
	{
		return;
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
			double dist = EuclideanBlockDist( clusterMat, k, j);   // need to be modified 
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

void SPDAnalysisModel::CityBlockDist( vnl_matrix<double>& mat, vnl_matrix<double>& matDis)
{
	unsigned int nrows = mat.rows();
	matDis.set_size( nrows, nrows);
	matDis.fill(1e6);

	for( int ind1 = 0; ind1 < mat.rows(); ind1++)
	{
//#pragma omp parallel for
		for( int ind2 = ind1 + 1; ind2 < mat.rows(); ind2++)
		{
			vnl_vector<double> module1 = mat.get_row(ind1);
			vnl_vector<double> module2 = mat.get_row(ind2);
			vnl_vector<double> subm = module1 - module2;
			double dis = 0;
			for( unsigned int i = 0; i < subm.size(); i++)
			{
				dis += abs(subm[i]);
			}
			matDis(ind1, ind2) = dis;
			matDis(ind2, ind1) = dis;
		}
	}
}

void SPDAnalysisModel::EuclideanBlockDist( vnl_matrix<double>& mat, vnl_matrix<double>& matDis)
{
	unsigned int nrows = mat.rows();
	matDis.set_size( nrows, nrows);
	matDis.fill(1e6);

	for( int ind1 = 0; ind1 < mat.rows(); ind1++)
	{
//#pragma omp parallel for
		for( int ind2 = ind1 + 1; ind2 < mat.rows(); ind2++)
		{
			vnl_vector<double> module1 = mat.get_row(ind1);
			vnl_vector<double> module2 = mat.get_row(ind2);
			vnl_vector<double> subm = module1 - module2;
			double dis = subm.magnitude();
			matDis(ind1, ind2) = dis;
			matDis(ind2, ind1) = dis;
		}
	}
}

void SPDAnalysisModel::EuclideanBlockDist( vnl_vector<double>& vec, vnl_matrix<double>& matDis)
{
    unsigned int size = vec.size();
    matDis.set_size( size, size);
    matDis.fill(1e6);
    for( int ind1 = 0; ind1 < size; ind1++)
    {
//#pragma omp parallel for
            for( int ind2 = ind1 + 1; ind2 <size; ind2++)
            {
                    double m1 = vec[ind1];
                    double m2 = vec[ind2];
                    double dis = abs( m2 - m1);
                    matDis(ind1, ind2) = dis;
                    matDis(ind2, ind1) = dis;
            }
    }
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
	histDis.fill(0);

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
	histDis.fill(0);

	if( inter < 1e-9)
	{
		return;
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

			if( index >= interval.size() - 1)
			{
				index = interval.size() - 2;
			}
			histDis[index]++;
		}
	}
}

void SPDAnalysisModel::Hist(vnl_vector<double>&distance, double interVal, double minVal, vnl_vector<unsigned int>& histDis, int num_bin = NUM_BIN)
{
	histDis.set_size( num_bin);
	histDis.fill(0);

	if( interVal < 1e-9)
	{
		return;
	}

	for( unsigned int i = 0; i < distance.size(); i++)
	{
		int index = floor( (distance[i] - minVal) / interVal);
		if(index >= 0)
		{
			if( index > num_bin - 1)
			{
				index = num_bin - 1;
			}
			histDis[index]++;
		}
	}
}

void SPDAnalysisModel::Hist(vnl_matrix<double>&distance, double interVal, double minVal, vnl_vector<unsigned int>& histDis, int num_bin = NUM_BIN)
{
	histDis.set_size( num_bin);
	histDis.fill(0);

	if( interVal < 1e-9)
	{
		return;
	}

	for( unsigned int i = 0; i < distance.rows(); i++)
	{
		for( unsigned int j = i + 1; j < distance.cols(); j++)
		{
			int index = floor( (distance(i,j) - minVal) / interVal);
			if( index >= 0)
			{
				if( index > num_bin - 1)
				{
					index = num_bin - 1;
				}
				histDis[index]++;
			}
		}
	}
}

double Dist(int *first, int *second)
{
	return abs(*first - *second);
}

double OnedirectDist(int *first, int *second)
{
	if( *first > *second)
	{
		return *first - *second;
	}
	else
	{
		return 0;
	}
}

double SPDAnalysisModel::EarthMoverDistance(vnl_vector<unsigned int>& first, vnl_vector<unsigned int>& second,
											vnl_matrix<double> &flowMatrix, int num_bin = NUM_BIN, bool bOneDirection)
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

	TsFlow *flow = new TsFlow[num_bin * num_bin];
	int flowVars = 0;

	double result = 0;
	if(bOneDirection)
	{
		result = transportSimplex(&srcSig, &snkSig, OnedirectDist, flow, &flowVars);
	}
	else
	{
		result = transportSimplex(&srcSig, &snkSig, Dist, flow, &flowVars);
	}
	
	//double cost = 0;
	//flowMatrix.fill(0);
	//if( bOneDirection)
	//{
	//	for( int i = 0; i < flowVars; i++)
	//	{
	//		if( flow[i].from > flow[i].to)
	//		{
	//			flowMatrix( flow[i].to, flow[i].from) = flow[i].amount;
	//			cost += flow[i].amount *( flow[i].from - flow[i].to);
	//		}
	//	}
	//}
	//else
	//{
	//	for( int i = 0; i < flowVars; i++)
	//	{
	//		flowMatrix( flow[i].to, flow[i].from) = flow[i].amount;
	//		cost += flow[i].amount * abs( flow[i].from - flow[i].to);
	//	}

	//}

	delete src;
	delete src_num;
	delete snk;
	delete snk_num;
	delete flow;

	return result;
}

void SPDAnalysisModel::RunEMDAnalysis(int numBin)
{
	std::ofstream ofs("EMD.txt");
	ofs.precision(4);
	vnl_vector<int> moduleSize = GetModuleSize( this->ClusterIndex);

	/** Generating the distance vector for module data and its MST data */
	std::vector<vnl_vector<double> > moduleDistance;
	moduleDistance.resize( this->ModuleMST.size());

	#pragma omp parallel for
	for( int i = 0; i <= this->ClusterIndex.max_value(); i++)
	{
		vnl_vector<double> distance;
		std::cout<< "build module distance "<< i<<endl;
		vnl_matrix<double> clusterData( this->MatrixAfterCellCluster.rows(), moduleSize[i]);
		GetCombinedMatrix( this->MatrixAfterCellCluster, this->ClusterIndex, i, i, clusterData);    /// Get the module data for cluster i
		GetMatrixDistance( clusterData, distance, CITY_BLOCK);
		moduleDistance[i] = distance;
	}

	this->EMDMatrix.set_size( this->ClusterIndex.max_value() + 1, this->ClusterIndex.max_value() + 1);
	this->EMDMatrix.fill(0);
	this->DistanceEMDVector.set_size( this->ClusterIndex.max_value() + 1);
	this->DistanceEMDVector.fill(0);

	for( int i = 0; i < moduleDistance.size(); i++)
	{
		if( i % 100 == 0)
		{
			std::cout<< "matching module " << i<<endl;
		}
		vnl_vector<double> hist_interval;
		vnl_vector<unsigned int> moduleHist;
		Hist( moduleDistance[i], numBin, hist_interval, moduleHist);
		ofs<< moduleHist<< endl;

		bool bMatch = true;
		unsigned int sumHist = moduleHist.sum();
		for(unsigned int k = 0; k < moduleHist.size(); k++)
		{
			if(moduleHist[k] >= sumHist * 0.9)
			{
				bMatch = false;
				std::cout<< "Remove module "<<i<<std::endl;
				break;
			}
		}

		if( bMatch)
		{
			#pragma omp parallel for
			for( int j = 0; j < moduleDistance.size(); j++)
			{	
				vnl_vector<double> mstDistance; 
				GetMSTMatrixDistance( moduleDistance[i], this->ModuleMST[j], mstDistance);  // get mst j distance from module i's distance matrix
				vnl_vector<unsigned int> mstHist;
				Hist( mstDistance, hist_interval, mstHist); 

				ofs<< mstHist<< endl;

				vnl_matrix<double> flowMatrix( numBin, numBin);
				double earth = EarthMoverDistance( moduleHist, mstHist, flowMatrix, numBin);
				this->EMDMatrix( i, j) = earth > 0 ? earth : 0;
			}

			///// matching with distance to device
			//if( disCor != 0)
			//{
			//	std::cout<< "matching module " << i<< " with distance MST "<<endl;
			//	vnl_vector<double> mstDistance; 
			//	GetMSTMatrixDistance( moduleDistance[i], this->DistanceMST, mstDistance);  // get mst j distance from module i's distance matrix
			//	vnl_vector<unsigned int> mstHist;    
			//	Hist( mstDistance, hist_interval, mstHist); 
			//	vnl_matrix<double> flowMatrix( numBin, numBin);
			//	this->DistanceEMDVector[ i] = EarthMoverDistance( moduleHist, mstHist, flowMatrix, numBin);
			//}
		}
	}

	moduleForSelection.clear();
	for( int i = 0; i < EMDMatrix.rows(); i++)
	{
		if( EMDMatrix(i,0) > -1e-9)
		{
			moduleForSelection.push_back(i);
		}
	}

	for(int i = 0; i < this->EMDMatrix.rows(); i++)
	{
		this->EMDMatrix(i,i) = 1;
		for( int j = 0; j < this->EMDMatrix.rows(); j++)
		{
			if( j != i)
			{
				double max = this->EMDMatrix(i,j) > this->EMDMatrix(j,i) ? this->EMDMatrix(i,j) : this->EMDMatrix(j,i);
				this->EMDMatrix(i,j) = max;
				this->EMDMatrix(j,i) = max;
			}
		}
	}
	this->EMDMatrix = this->EMDMatrix / this->EMDMatrix.max_value();

	ofs<< this->EMDMatrix<< std::endl;
	ofs.close();
	//for( unsigned int i = 0; i < this->EMDMatrix.rows(); i++)
	//{
	//	double max = this->EMDMatrix.get_row(i).max_value();
	//	for( unsigned int j = 0; j < this->EMDMatrix.cols(); j++)
	//	{
	//		if( max != 0)
	//		{
	//			this->EMDMatrix( i, j) = this->EMDMatrix( i, j) / max;
	//		}
	//	}
	//	if( max != 0)
	//	{
	//		this->DistanceEMDVector[ i] = this->DistanceEMDVector[ i] / max;
	//	}
	//}
	std::cout<< "EMD matrix has been built successfully"<<endl;
}

void SPDAnalysisModel::GetEMDMatrixDivByMax(vnl_matrix<double> &emdMatrix)
{
	emdMatrix =  this->EMDMatrix;
}

//void SPDAnalysisModel::GetClusClusData(clusclus *c1, clusclus *c2, double threshold, std::vector< unsigned int> *disModIndex)
//{
//	QString filenameSM = this->filename + "similarity_matrix.txt";
//	std::ofstream ofSimatrix(filenameSM .toStdString().c_str(), std::ofstream::out);
//	ofSimatrix.precision(4);
//
//	this->heatmapMatrix.set_size( this->EMDMatrix.rows(), this->EMDMatrix.cols());
//	this->heatmapMatrix.fill(0);
//	std::vector< unsigned int> simModIndex;
//
//	if( bProgression)
//	{
//		ofSimatrix<< "Progression over distance"<<endl;
//		for(unsigned int i = 0; i < DistanceEMDVector.size(); i++)
//		{
//			if( DistanceEMDVector[i] <= threshold)    // change to smaller than
//			{
//				disModIndex->push_back(i);
//				ofSimatrix<< i<< "\t"<< DistanceEMDVector[i]<<endl;
//			}
//		}
//	}
//	else
//	{
//		ofSimatrix<< "Overal Progression"<<endl;
//		for( unsigned int i = 0; i < this->EMDMatrix.cols(); i++)
//		{
//			ofSimatrix <<"MST "<<i<<endl;
//			for( unsigned int j = 0; j < this->EMDMatrix.rows(); j++)
//			{
//				if( this->EMDMatrix( j, i) <= threshold)     // find the modules that support common mst    change to smaller than
//				{
//					simModIndex.push_back(j);
//					ofSimatrix << j <<"\t";
//				}
//			}
//			ofSimatrix <<endl;
//
//			if( simModIndex.size() > 0)
//			{
//				for( unsigned int j = 0; j < simModIndex.size(); j++)
//				{
//					for( unsigned int k = j + 1; k < simModIndex.size(); k++)
//					{
//						this->heatmapMatrix( simModIndex[j], simModIndex[k]) = heatmapMatrix( simModIndex[j], simModIndex[k]) + 1;
//						this->heatmapMatrix( simModIndex[k], simModIndex[j]) = heatmapMatrix( simModIndex[k], simModIndex[j]) + 1;
//					}
//					this->heatmapMatrix(simModIndex[j], simModIndex[j]) = this->heatmapMatrix(simModIndex[j], simModIndex[j]) + 1;
//				}
//			}
//			simModIndex.clear();
//		}
//
//		c1->Initialize( heatmapMatrix.data_array(), this->EMDMatrix.rows(), this->EMDMatrix.cols());
//		c1->RunClusClus();
//		c1->Transpose();
//
//		c2->Initialize( c1->transposefeatures,c1->num_features, c1->num_samples);
//		c2->RunClusClus();
//	}
//	ofSimatrix.close();	
//}

void SPDAnalysisModel::GetClusClusDataMST(clusclus *c1, double threshold, std::vector< unsigned int> *disModIndex)
{
	QString filenameSM = this->filename + "similarity_matrix.txt";
	std::ofstream ofSimatrix(filenameSM .toStdString().c_str(), std::ofstream::out);
	ofSimatrix.precision(4);

	this->heatmapMatrix.set_size( this->EMDMatrix.rows(), this->EMDMatrix.cols());
	this->heatmapMatrix.fill(0);
	std::vector< unsigned int> simModIndex;

	if( bProgression)
	{
		//ofSimatrix<< "Progression over distance"<<endl;
		for(unsigned int i = 0; i < DistanceEMDVector.size(); i++)
		{
			if( DistanceEMDVector[i] <= threshold)    // change to smaller than
			{
				disModIndex->push_back(i);
				//ofSimatrix<< i<< "\t"<< DistanceEMDVector[i]<<endl;
			}
		}
	}
	else
	{
		//ofSimatrix<< "Overal Progression"<<endl;
		//vnl_vector<int> moduleSize = GetModuleSize( this->ClusterIndex);
		for( unsigned int i = 0; i < moduleForSelection.size(); i++)
		{
			int moduleI = moduleForSelection[i];
			//ofSimatrix <<"MST "<<moduleI<<endl;
			for( unsigned int j = 0; j < moduleForSelection.size(); j++)
			{
				if( this->EMDMatrix( moduleI, moduleForSelection[j]) >= threshold)     // find the modules that support common mst    change to smaller than
				{
					simModIndex.push_back(moduleForSelection[j]);
					//ofSimatrix << moduleForSelection[j] <<"\t";
				}
			}
			//ofSimatrix <<endl;

			if( simModIndex.size() > 1)   // not only include itself
			{
				for( unsigned int j = 0; j < simModIndex.size(); j++)
				{
					for( unsigned int k = j + 1; k < simModIndex.size(); k++)
					{
						this->heatmapMatrix( simModIndex[j], simModIndex[k]) = heatmapMatrix( simModIndex[j], simModIndex[k]) + 1;//moduleSize[moduleI];
						this->heatmapMatrix( simModIndex[k], simModIndex[j]) = heatmapMatrix( simModIndex[k], simModIndex[j]) + 1;//moduleSize[moduleI];
					}
					this->heatmapMatrix(simModIndex[j], simModIndex[j]) = this->heatmapMatrix(simModIndex[j], simModIndex[j]) + 1;//moduleSize[moduleI];
				}
			}
			simModIndex.clear();
		}

		c1->Initialize( this->heatmapMatrix.data_array(), this->EMDMatrix.rows(), this->EMDMatrix.cols());
		c1->RunClusClus();
	}

	ofSimatrix << this->heatmapMatrix<< std::endl;
	ofSimatrix.close();	
}

void SPDAnalysisModel::GetBiClusData(clusclus *c1, vnl_vector<double> *diagVec)
{
	if( c1 == NULL)
	{
		return;
	}
	c1->Initialize( this->EMDMatrix.data_array(), this->EMDMatrix.rows(), this->EMDMatrix.cols());

	double max = this->EMDMatrix.max_value();
	double min = this->EMDMatrix.min_value();

	std::set< unsigned int> mIds;
	std::vector< std::vector< unsigned int> > order;
	std::vector< std::vector< unsigned int> > tmpOrder;

	for(double selThreshold = max; selThreshold >= min; selThreshold -= 0.01)
	{
		std::vector<std::vector<unsigned int> > modID;

		GetSelectedFeaturesModulesForBlockVisualization(selThreshold, modID); // the opposite way, output include the previous!!

		for( size_t i = 0; i < modID.size(); i++)
		{
			bool tag = false;
			std::vector<unsigned int> tmp;
			for( size_t j = 0; j < order.size(); j++)
			{
				if( IsExist(modID[i], order[j][0]))
				{
					tag = true;
					
					for( size_t k = 0; k< order[j].size(); k++) 
					{
						tmp.push_back(order[j][k]);
					}
				}
			}

			if( modID[i].size() > tmp.size())
			{
				for( size_t k = 0; k < modID[i].size(); k++)
				{
					if( !IsExist(tmp, modID[i][k]))
					{
						tmp.push_back( modID[i][k]);
					}
				}
			}
		
			tmpOrder.push_back(tmp);
		}
		order = tmpOrder;
		tmpOrder.clear();
	}

	std::vector< unsigned int> lastVec;

	if( order[0].size() < this->EMDMatrix.rows())
	{
		std::cout<< "Last searching for the missed component."<<std::endl;
		for( unsigned int i = 0; i < this->EMDMatrix.rows(); i++)
		{
			
			bool bexist = false;
			for( size_t j = 0; j < order.size(); j++)
			{
				if( IsExist(order[j], i))
				{
					bexist = true;
					break;
				}
			}
			if( false == bexist)
			{
				lastVec.push_back(i);
			}
		}
	}

	int count = 0;
	for( size_t i = 0; i < order.size(); i++)
	{
		for( size_t j = 0; j < order[i].size(); j++)
		{
			if( count < this->EMDMatrix.rows())
			{
				c1->optimalleaforder[count++] = order[i][j];
			}
			else
			{
				std::cout<< "Exceed matrix dim!!"<<std::endl;
			}
		}
	}

	if( count <  this->EMDMatrix.rows())
	{
		for( size_t k = 0; k < lastVec.size(); k++)
		{
			if( count < this->EMDMatrix.rows())
			{
				c1->optimalleaforder[count++] = lastVec[k];
			}
			else
			{
				std::cout<< "Exceed matrix dim!!"<<std::endl;
			}
		}
	}

	if(diagVec)
	{
		diagVec->set_size(this->EMDMatrix.rows());
		for(unsigned int i = 0; i < this->EMDMatrix.rows(); i++)
		{
			(*diagVec)[i] = 1;
		}
	}
}

/// Find largest fully connected component in EMDMatrix.
void SPDAnalysisModel::GetSelectedFeaturesModulesTest(double selThreshold, std::vector<unsigned int> &selModules, std::vector<unsigned int> &size)
{
	selModules.clear();
	size.clear();
	unsigned int maxSize = 0;
	std::vector< std::vector<unsigned int> > tmpSelModules;
	std::set< unsigned int> processedModules;
	std::vector< unsigned int> tmpModules;

	for( unsigned int i = 0; i < this->EMDMatrix.rows(); i++)
	{
		if(processedModules.find(i) == processedModules.end())
		{
			processedModules.insert(i);
			tmpModules.clear();
			tmpModules.push_back(i);
			for( unsigned int j = i + 1; j < this->EMDMatrix.cols(); j++)
			{
				if( this->EMDMatrix(i,j) >= selThreshold)
				{
					tmpModules.push_back(j);
					processedModules.insert(j);
				}
			}
			tmpSelModules.push_back(tmpModules);
			size.push_back(tmpModules.size());
		}
	}
	
	//std::ofstream ofs("SelectedFeaturesModulesTest.txt");
	//std::cout<< "Module size > 2:"<<std::endl;
	for( unsigned int i = 0; i < tmpSelModules.size(); i++)
	{
		std::vector<unsigned int> module = tmpSelModules[i];
		if( module.size() > 2)
		{
			std::cout<< module.size()<< "\t";
			for( unsigned int j = 0; j < module.size(); j++)
			{
				selModules.push_back(module[j]);			
			}
		}
	}
	std::cout<< std::endl;
	//ofs.close();
}

void SPDAnalysisModel::GetSelectedFeaturesModulesForBlockVisualization(double selThreshold, std::vector< std::vector<unsigned int> > &tmpSelModules)
{
	tmpSelModules.clear();
	std::set< unsigned int> processedModules;
	std::vector< unsigned int> tmpModules;

	for( unsigned int i = 0; i < this->EMDMatrix.rows(); i++)
	{
		if(processedModules.find(i) == processedModules.end())
		{
			processedModules.insert(i);
			tmpModules.clear();
			tmpModules.push_back(i);
			for( unsigned int j = i + 1; j < this->EMDMatrix.cols(); j++)
			{
				if( this->EMDMatrix(i,j) >= selThreshold)
				{
					tmpModules.push_back(j);
					processedModules.insert(j);
				}
			}
			if( tmpModules.size() >= 2)
			{
				tmpSelModules.push_back(tmpModules);
			}
		}
	}
}

double SPDAnalysisModel::GetConnectionAccuracy( vtkSmartPointer<vtkTable> treeTable, vnl_matrix<double> &disMat, vnl_vector<double> &accuracyVec,
											    vnl_vector<double> &aggDegree, double &aggDegreeValue, int neighborScope = 1, int clusterScope = 0)
{
	if( clusNo.max_value() <= 0)
	{
		std::cout<< "Validation information not available"<<std::endl;
		return -1;
	}

	vtkAbstractArray *arrayID1 = treeTable->GetColumnByName((this->headers)[0].c_str());
	vtkAbstractArray *arrayID2 = treeTable->GetColumnByName((this->headers)[1].c_str());
	vnl_matrix<int> connectionCount( clusNo.max_value(),2); //first row count correct connection, second row count wrong connection
	connectionCount.fill(0);

	std::vector< std::vector<int> > connectionNode( 2 * clusNo.max_value());
	std::vector< std::vector<int> > ncNode( clusNo.max_value());

	Graph graph(treeTable->GetNumberOfRows() + 1);

	for( vtkIdType i = 0; i < treeTable->GetNumberOfRows(); i++) 
	{
		int ver1 = arrayID1->GetVariantValue(i).ToInt();
		int ver2 = arrayID2->GetVariantValue(i).ToInt();
		int tag1 = clusNo[ver1];
		int tag2 = clusNo[ver2];
		
		if( abs(tag1 - tag2) <= clusterScope)
		{
			boost::add_edge( ver1, ver2, graph);
		}

		if( abs(tag1 - tag2) == 0)
		{
			ncNode[tag1-1].push_back(ver1);
			ncNode[tag2-1].push_back(ver2);
		}
		else
		{
			if( abs(tag1 - tag2) <= neighborScope) // right connection
			{
				connectionCount(tag1-1,0) += 1;
				connectionCount(tag2-1,0) += 1;	
				connectionNode[ 2 * (tag1-1)].push_back(ver1);
				connectionNode[ 2 * (tag2-1)].push_back(ver2);
			}
			else if( abs(tag1 - tag2) > neighborScope) // wrong connection
			{
				connectionCount(tag1-1,1) += 1;
				connectionCount(tag2-1,1) += 1;
				connectionNode[ 2 * (tag1-1) + 1].push_back(ver1);
				connectionNode[ 2 * (tag2-1) + 1].push_back(ver2);
			}
		}
	}

	for( size_t i = 0; i < ncNode.size(); i++) 
	{
		std::vector<int> nodes = ncNode[i];
		for( size_t j = 0; j < nodes.size(); j++)
		{
			int ucn = nodes[j];
			double minDis = 1e9;
			bool brightConnection = true;
			for( size_t k = 0; k < connectionNode[2*i].size(); k++)
			{
				int cn = connectionNode[2*i][k];
				if( cn != ucn && disMat(ucn, cn) < minDis)   // probably problem here~~
				{
					minDis = disMat(ucn, cn);
					brightConnection = true;
				}
			}
			for(size_t k = 0; k < connectionNode[2*i + 1].size(); k++)
			{
				int cn = connectionNode[2*i + 1][k];
				if( cn != ucn && disMat(ucn, cn) < minDis)
				{
					minDis = disMat(ucn, cn);
					brightConnection = false;
					break;
				}
			}
			if( brightConnection)
			{
				connectionCount(i,0) += 1;
			}
			else
			{
				connectionCount(i,1) += 1;
			}
		}
	}

	std::cout<< connectionCount<<std::endl;

	accuracyVec.set_size(clusNo.max_value());
	double totalConnection = 0;
	double accurateConnection = 0;
	for( vtkIdType i = 0; i < accuracyVec.size(); i++)
	{
		accuracyVec[i] = (double)(connectionCount(i,0) - connectionCount(i,1))/(connectionCount(i,0) + connectionCount(i,1));
		accurateConnection += connectionCount(i,0);
		totalConnection += (connectionCount(i,0) + connectionCount(i,1));
	}
	double rtnVal = accurateConnection / totalConnection;

	aggDegree.set_size(clusNo.max_value());
	aggDegree.fill(0);
	
	if( boost::num_vertices(graph) > 0)
	{
		std::vector< int> component;
		component.resize(boost::num_vertices(graph));
		int num = boost::connected_components(graph, &component[0]);
		std::cout<< "Connected clusters: "<< num<<std::endl;
		std::vector< std::vector< int> > connectedNodes(num);
		vnl_vector< int> state(num);
		state.fill(-1);
		for( size_t i = 0; i < component.size(); i++)
		{
			int index = component[i];
			connectedNodes[index].push_back(i);
			if( state[index] == -1)
			{
				state[index] = clusNo[i] - 1;   /// problem!!!!
			}
		}

		vnl_vector<double> size(clusNo.max_value());
		size.fill(0);
		vnl_vector<double> maxAgg(clusNo.max_value());
		maxAgg.fill(1);
		for( size_t i = 0; i < connectedNodes.size(); i++)
		{
			int clus = state[ i];
			size[clus] += connectedNodes[i].size();
			if( connectedNodes[i].size() > maxAgg[clus])
			{
				maxAgg[clus] = connectedNodes[i].size();	
			}
		}
		std::cout<< "Size: "<<size<<std::endl;
		double totalSize = 0;
		double aggSize = 0;
		for( unsigned int i = 0; i < aggDegree.size(); i++)
		{
			aggDegree[i] = maxAgg[i] / size[i];
			aggSize += maxAgg[i];
			totalSize += size[i];
		}
		aggDegreeValue = aggSize / totalSize;
	}
	else
	{
		aggDegreeValue = 1;
	}
	
	return rtnVal;
}

void SPDAnalysisModel::GetValidationVec(vnl_vector<int> &validationVec)
{
	validationVec = clusNo;
}

void SPDAnalysisModel::GetClusClusDataKNNG(clusclus *c1,vnl_vector<double> *diagVec, std::vector< unsigned int> *disModIndex)
{
	//QString filenameSM = this->filename + "similarity_matrix.txt";
	//std::ofstream ofSimatrix(filenameSM .toStdString().c_str(), std::ofstream::out);
	//ofSimatrix.precision(4);
	if( c1 == NULL)
	{
		return;
	}
	c1->Initialize( this->EMDMatrix.data_array(), this->EMDMatrix.rows(), this->EMDMatrix.cols());
	c1->RunClusClus();

	unsigned int nrows = this->EMDMatrix.rows();
	vnl_matrix<double> EMDMat(nrows, nrows);
	vnl_vector<int> order(nrows);

	for( unsigned int i = 0; i < nrows; i++)
	{
		order[i] = c1->optimalleaforder[i];
	}

	vnl_vector<double> vec;
	bool iterate = true;
	int count = 0;
	while(iterate)
	{
		for(int i = 0; i < nrows; i++)
		{
			for(int j = 0; j < nrows; j++)
			{
				EMDMat(i,j) = EMDMatrix(order[i],order[j]);
			}
		}
		iterate = DiagnalIteration(EMDMat, order, vec);

		count++;
		if(count >= 10)
		{
			break;
		}
	}
	std::cout<<"Converge after "<<count<<std::endl;

	// arrange according to the first row
	//for(int i = 0; i < nrows; i++)
	//{
	//	for(int j = 0; j < nrows; j++)
	//	{
	//		EMDMat(i,j) = EMDMatrix(order[i],order[j]);
	//	}
	//}
	//for( unsigned int i = 0; i < order.size(); i++)
	//{
	//	for( unsigned int j = i + 1; j < order.size(); j++)
	//	{
	//		if( EMDMat(0,i) < EMDMat(0,j))
	//		{
	//			std::swap(vec[i],vec[j]);
	//			std::swap(order[i],order[j]);
	//			std::swap(EMDMat(0,i), EMDMat(0,j));
	//		}
	//	}
	//}
	//for(int i = 0; i < nrows; i++)
	//{
	//	for(int j = 0; j < nrows; j++)
	//	{
	//		EMDMat(i,j) = EMDMatrix(order[i],order[j]);
	//	}
	//}

	// save the order
	for(unsigned int i = 0; i < nrows; i++)
	{
		c1->optimalleaforder[i] = order[i];
	}

	if(diagVec)
	{
		diagVec->set_size(nrows);
		for(unsigned int i = 0; i < nrows; i++)
		{
			(*diagVec)[i] = vec(i);
		}
	}
	//ofSimatrix << EMDMat<< std::endl;
	//ofSimatrix.close();	
}

bool SPDAnalysisModel::DiagnalIteration(vnl_matrix<double> &mat, vnl_vector<int> &order, vnl_vector<double> &vec)
{
	unsigned int nrows = mat.rows();
	vec.set_size(nrows);
	for( unsigned int i = 0; i < nrows; i++)
	{
		if( i == 0 )
		{
			mat(0,0) = mat(1,0);
		}
		else if( i == nrows - 1)
		{
			mat(nrows - 1, nrows - 1) = mat(nrows - 1, nrows - 2);
		}
		else
		{
			mat(i,i) = mat(i,i+1) > mat(i-1,i) ? mat(i,i+1) : mat(i-1,i);
		}
		vec[i] = mat(i,i);
	}

	bool bstate = false;
	for( unsigned int i = 0; i < vec.size(); i++)
	{
		for( unsigned int j = i + 1; j < vec.size(); j++)
		{
			if( vec[i] < vec[j])
			{
				std::swap(vec[i],vec[j]);
				std::swap(order[i],order[j]);
				bstate = true;
			}
		}
	}
	return bstate;
}

void SPDAnalysisModel::GetClusClusPSCWithoutIterData(clusclus* c1, double threshold)
{
	QString filenameSM = this->filename + "similarity_matrix.txt";
	std::ofstream ofSimatrix(filenameSM .toStdString().c_str(), std::ofstream::out);
	ofSimatrix.precision(4);

	this->heatmapMatrix.set_size( this->EMDMatrix.rows(), this->EMDMatrix.cols());
	this->heatmapMatrix.fill(0);

	for( unsigned int i = 0; i < this->EMDMatrix.rows(); i++)
	{
		std::vector< unsigned int> simModIndex;
		for( unsigned int j = 0; j < this->EMDMatrix.cols(); j++)
		{
			if( this->EMDMatrix( i, j) >= threshold)
			{
				simModIndex.push_back(j);
			}
		}
		//ofSimatrix <<endl;

		if( simModIndex.size() > 1)   // not only include itself
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
	}

	c1->Initialize( EMDMatrix.data_array(), this->EMDMatrix.rows(), this->EMDMatrix.cols());
	c1->RunClusClus();

	ofSimatrix << this->EMDMatrix<< std::endl;
	ofSimatrix.close();	
}

void SPDAnalysisModel::GetClusClusPSCData(clusclus* c1)
{
	if(c1)
	{
		c1->Initialize( this->EMDMatrix.data_array(), this->EMDMatrix.rows(), this->EMDMatrix.cols());
		c1->RunClusClus();
	}
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
			double dist = EuclideanBlockDist( clusterMat, k, j);   // need to be modified 
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
	vtkSmartPointer<vtkDoubleArray> column;

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
			double dist = EuclideanBlockDist( clusterMat, i, vertex[i]);
			vtkSmartPointer<vtkVariantArray> DataRow = vtkSmartPointer<vtkVariantArray>::New();
 
				DataRow->InsertNextValue( i);
				DataRow->InsertNextValue( vertex[i]);
			//}
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

void SPDAnalysisModel::SaveSelectedFeatureNames(QString filename, std::vector<int>& selectedFeatures)
{
	std::ofstream ofs( filename.toStdString().c_str(), std::ofstream::out);

	for( int i = 0; i < selectedFeatures.size(); i++)
	{
		char* name = this->DataTable->GetColumn(selectedFeatures[i] + 1)->GetName();
		std::string featureName(name);
		ofs<<featureName<<endl;
	}
	ofs.close();
}

void SPDAnalysisModel::SaveSelectedFeatureNames(QString filename, std::vector<unsigned int>& selectedFeatures)
{
	QString file = this->filename + filename;
	std::ofstream ofs( file.toStdString().c_str(), std::ofstream::out);

	for( int i = 0; i < selectedFeatures.size(); i++)
	{
		char* name = this->DataTable->GetColumn(selectedFeatures[i] + 1)->GetName();
		std::string featureName(name);
		ofs<< selectedFeatures[i]<<"\t"<<featureName<<endl;
	}
	ofs.close();
}

void SPDAnalysisModel::SaveNormalizedTableAfterFeatureSelection(std::string filename, std::vector<int>& selectedFeatures)
{
	vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
	//for(int i = 0; i < selectedFeatures.size(); i++)
	//{		
	//	vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
	//	column->SetName( DataTable->GetColumn(selectedFeatures[i] + 1)->GetName());
	//	table->AddColumn(column);
	//}

	std::string fileNameOri = filename + "_ori.txt";
	std::string fileNameSel = filename + ".txt";

	vtkSmartPointer<vtkTable> normalTable = GetDataTableAfterCellCluster();
	ftk::SaveTable(fileNameOri,normalTable);

	table->AddColumn(normalTable->GetColumn(0));
	for(int i = 0; i < selectedFeatures.size(); i++)
	{		
		table->AddColumn(normalTable->GetColumn(selectedFeatures[i] + 1));
	}
	ftk::SaveTable(fileNameSel,table);
}

double SPDAnalysisModel::GetEMDSelectedPercentage(double thres)
{
	unsigned int count = 0;
	double per = 0;
	if( false == bProgression)
	{
		for( unsigned int i = 0; i < this->EMDMatrix.rows(); i++)
		{
			for( unsigned int j = 0; j < this->EMDMatrix.cols(); j++)
			{
				if( this->EMDMatrix( i, j) >= thres && this->EMDMatrix( i, j) > 0)     // find the modules that support common mst
				{
					count++;
				}
			}
		}
		unsigned int allnum = this->EMDMatrix.cols() * this->EMDMatrix.rows();
		per = (double)count / allnum;
	}
	else
	{
		for( unsigned int i = 0; i < this->DistanceEMDVector.size(); i++)
		{
			if( this->DistanceEMDVector[i] <= thres)     // find the modules that support common mst
			{
				count++;
			}
		}
		per = (double) count / this->DistanceEMDVector.size();
	}
	return per;
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

void SPDAnalysisModel::GetClusterMatrixData(vnl_matrix<double> &mat)
{
	mat = UNMatrixAfterCellCluster;
}

void SPDAnalysisModel::GetClusterOrder(std::vector< std::vector< long int> > &clusIndex, std::vector<long int> &treeOrder, std::vector< int> &clusterOrder)
{
	clusterOrder.clear();
	for( int i = 0; i < treeOrder.size(); i++)
	{
		std::vector< long int> order = clusIndex[ treeOrder[i]];
		for( int j = 0; j < order.size(); j++)
		{
			clusterOrder.push_back(order[j]);
		}
	}
}

void SPDAnalysisModel::GetValidationOrder(std::vector< int> &clusterOrder)
{
	if(clusNo.max_value() >= 1)
	{
		clusterOrder.clear();
		std::vector< std::vector<int> > index;
		index.resize(clusNo.max_value());
		for(unsigned int i = 0; i < clusNo.size(); i++)
		{
			int clus = clusNo[i] - 1;
			index[clus].push_back(i);
		}
		for(size_t i = 0; i < index.size(); i++)
		{
			for( size_t j = 0; j < index[i].size(); j++)
			{
				clusterOrder.push_back(index[i][j]);
			}
		}
	}
}

/// for spdtestwindow
void SPDAnalysisModel::ModuleCoherenceMatchAnalysis()
{
	int size = ClusterIndex.max_value() + 1;

	vnl_matrix<double> mat( this->MatrixAfterCellCluster.rows(), size);
	std::vector< std::vector< unsigned int> > featureClusterIndex;
	featureClusterIndex.resize( size);
	for( int i = 0; i < ClusterIndex.size(); i++)
	{
		unsigned int index = ClusterIndex[i];
		featureClusterIndex[index].push_back(i);
	}
	for( int i = 0; i < featureClusterIndex.size(); i++)
	{
		vnl_vector< double> vec( this->MatrixAfterCellCluster.rows());
		vec.fill(0);
		for( int j = 0; j < featureClusterIndex[i].size(); j++)
		{
			vec += this->MatrixAfterCellCluster.get_column( featureClusterIndex[i][j]);
		}
		vec = vec / featureClusterIndex[i].size();
		mat.set_column(i, vec);
	}

	QString filenameSM = this->filename + "module_coherence_match.txt";
	std::ofstream ofSimatrix(filenameSM .toStdString().c_str(), std::ofstream::out);
	ofSimatrix.precision(4);
	
	ModuleCompareCorMatrix.set_size(size, size);
	ModuleCompareCorMatrix.fill(0);
	DistanceCorVector.set_size(size);
	DistanceCorVector.fill(0);

	ofSimatrix<< mat<<endl<<endl;

	for( int i = 0; i < featureClusterIndex.size(); i++)
	{
		if( disCor > 0)
		{
			vnl_matrix< double> module(this->MatrixAfterCellCluster.rows(), featureClusterIndex[i].size());
			GetCombinedMatrix( this->MatrixAfterCellCluster, ClusterIndex, i, i, module);
			module.normalize_columns();
			vnl_vector< double> disVec = DistanceToDevice.normalize();
			vnl_vector< double> discorvec = disVec * module;
			DistanceCorVector[i] = discorvec.mean();
			ofSimatrix << "module "<< i<< "\t"<< DistanceCorVector[i]<<endl;
		}
		
		for( int j = 0; j < featureClusterIndex.size(); j++)
		{
			if( i != j)
			{
				ofSimatrix << "module "<<i<<"\t"<<j<<endl;
				vnl_vector< double> veci = mat.get_column(i);
				vnl_vector< double> vecj = mat.get_column(j);
				int ijsize = featureClusterIndex[i].size() + featureClusterIndex[j].size();

				//vnl_vector< double> mean = (veci * featureClusterIndex[i].size() + vecj * featureClusterIndex[j].size()) / ijsize;
				//mean = mean.normalize();
				//vnl_matrix< double> newmodule(this->MatrixAfterCellCluster.rows(), ijsize);
				//GetCombinedMatrix( this->MatrixAfterCellCluster, ClusterIndex, i, j, newmodule);
				//newmodule.normalize_columns();
				//vnl_vector< double> corvec = mean * newmodule;
				//CorMatrix(i,j) = corvec.mean();
				//CorMatrix(j,i) = CorMatrix(i,j);

				///// calculate negative cormatrix, inverse vector j
				//vecj = -vecj;
				//mean = (veci * featureClusterIndex[i].size() + vecj * featureClusterIndex[j].size()) / ijsize;
				//mean = mean.normalize();
				//GetCombinedInversedMatrix( this->MatrixAfterCellCluster, ClusterIndex, i, j, newmodule);  // module j is inversed
				//newmodule.normalize_columns();
				//corvec = mean * newmodule;
				//ofSimatrix << corvec<<endl;
				//ModuleCompareCorMatrix(i,j) = corvec.mean();
				//ModuleCompareCorMatrix(j,i) = ModuleCompareCorMatrix(i,j);
				veci = veci.normalize();
				ijsize = featureClusterIndex[j].size();
				vnl_matrix< double> newmodule(this->MatrixAfterCellCluster.rows(), ijsize);
				GetCombinedMatrix( this->MatrixAfterCellCluster, ClusterIndex, j, j, newmodule);  // module j is inversed
				newmodule.normalize_columns();
				vnl_vector< double> corvec = veci * newmodule;
				ofSimatrix << corvec<<endl;
				ModuleCompareCorMatrix(i,j) = corvec.mean();
			}
		}
	}
	ofSimatrix.close();
}

unsigned int SPDAnalysisModel::GetKNeighborNum()
{
	return m_kNeighbor;
}

void SPDAnalysisModel::GetDiagnalMinMax(vnl_matrix<double> &mat, double &min, double &max)
{
	min = mat(0,1);
	max = mat(0,1);
	for(unsigned int i = 0; i < mat.rows(); i++)
	{
		for( unsigned int j = i + 1; j < mat.cols(); j++)
		{
			if(mat(i,j) > max)
			{
				max = mat(i,j);
			}
			if(mat(i,j) < min)
			{
				min = mat(i,j);
			}
		}
	}
}

void SPDAnalysisModel::ModuleCorrelationMatrixMatch2(unsigned int kNeighbor, int nbins)
{
	m_kNeighbor = kNeighbor;
	std::cout<< m_kNeighbor<<std::endl;
	std::ofstream ofs("ModuleCorrelationMatrixMatch2.txt");
	int size = ClusterIndex.max_value() + 1;
	std::vector< std::vector< unsigned int> > featureClusterIndex;
	featureClusterIndex.resize( size);
	for( int i = 0; i < ClusterIndex.size(); i++)
	{
		unsigned int index = ClusterIndex[i];
		featureClusterIndex[index].push_back(i);
	}
	this->EMDMatrix.set_size( featureClusterIndex.size(),  featureClusterIndex.size());
	this->EMDMatrix.fill(-1);

	std::set< int> excludeModule;
	std::vector< vnl_vector<int> > nearIndex(featureClusterIndex.size());
	vnl_vector< unsigned char> state(featureClusterIndex.size());
	state.fill(1);

	#pragma omp parallel for
	for( int i = 0; i < featureClusterIndex.size(); i++)
	{
		vnl_matrix< double> modulei(this->MatrixAfterCellCluster.rows(), featureClusterIndex[i].size());
		GetCombinedMatrix( this->MatrixAfterCellCluster, ClusterIndex, i, i, modulei);

		vnl_vector< double> meanVec( modulei.rows());
		meanVec.fill(0);
		for( unsigned int j = 0; j < modulei.cols(); j++)
		{
			meanVec += modulei.get_column(j);
		}
		meanVec = meanVec / modulei.cols();
		double mmin = meanVec.min_value();
		double mmax = meanVec.max_value();
		
		if( mmax - mmin > 1e-9)
		{
			bool bstate = true;
			vnl_vector<unsigned int> histMean;
			unsigned int binTest = meanVec.size();
			if( binTest * 0.1 > 10)
			{
				binTest = 10;
			}
			else if( binTest * 0.1 < 3)
			{
				binTest = 3;
			}
			double minterval = (mmax - mmin) / binTest;
			Hist(meanVec, minterval, mmin, histMean, binTest);
			unsigned int tmpMax = histMean.max_value();
			unsigned int maxHisti = (unsigned int)(histMean.sum() * 0.9);   // if the module's full graph edge distribution highly agglomerate at one bin, then discard this module.
			if(tmpMax >= maxHisti)
			{
				std::cout<< "1 Remove module "<<i<<std::endl;
				state[i] = 0;
			}
			else
			{
				unsigned int zeroSum = 0;
				for( unsigned int s = 0; s < histMean.size(); s++)
				{
					if(histMean[s] < 1)
					{
						zeroSum++;
					}
				}
				if( zeroSum >= 2.0 / 3.0 * histMean.size())
				{
					std::cout<< "Do not use module "<<i<<" for iteration.\t"<<std::endl;
					excludeModule.insert(i);
					state[i] = 0;
				}
			}

			if(state[i])
			{
				double * vecArray = meanVec.data_block();
				vnl_vector<int> idArray(meanVec.size());
				for( int c = 0; c < meanVec.size(); c++)
				{
					idArray[c] = c;
				}
				quickSort(vecArray, 0, meanVec.size() - 1, idArray.data_block());
				nearIndex[i] = idArray;
			}
		}
	}

	///neighbor match process
	
	for( int i = 0; i < featureClusterIndex.size(); i++)
	{
		if( state[i])
		{
			vnl_vector< double> row;
			row.set_size(featureClusterIndex.size());
			row.fill(0);

			vnl_matrix< double> modulei(this->MatrixAfterCellCluster.rows(), featureClusterIndex[i].size());
			GetCombinedMatrix( this->MatrixAfterCellCluster, ClusterIndex, i, i, modulei);
			vnl_matrix< double> modDisti;
			EuclideanBlockDist(modulei, modDisti);

			double min = 0;
			double max = 0;
			GetDiagnalMinMax(modDisti, min, max);
			double interval = ( max - min) / nbins;

			vnl_vector<double> nearWeights;
			GetWeightsFromSorting( modDisti, nearIndex[i], nearWeights, kNeighbor);

			vnl_vector<unsigned int> histMod, histNear;
			Hist(modDisti, interval, min, histMod, nbins); 
			Hist(nearWeights, interval, min, histNear, nbins);
			vnl_matrix<double> flowMatrix( nbins, nbins);
			double disScale = EarthMoverDistance( histMod, histNear, flowMatrix, nbins, true);

			if( i % 100 == 0)
			{
				std::cout<< "Matching module "<<i<<std::endl; 
			}

			#pragma omp parallel for
			for( int j = 0; j < featureClusterIndex.size(); j++)
			{
				if(state[j])
				{
					if( j == i)
					{
						row[j] = 1;
					}
					else 
					{
						vnl_vector<double> matchWeights;
						vnl_vector<unsigned int> histj;
						
						GetWeightsFromSorting( modDisti, nearIndex[j], matchWeights, kNeighbor);
						Hist( matchWeights, interval, min, histj, nbins);
						vnl_matrix<double> flowMatrix( nbins, nbins);
						double earth = EarthMoverDistance( histMod, histj, flowMatrix, nbins, true);
						row[j] = earth > 0 ? earth : 0;
						if( abs(disScale) > 1e-6)
						{
							row[j] = row[j] / disScale;
						}
					}
				}
			}
			this->EMDMatrix.set_row(i, row);
		}
	}

	moduleForSelection.clear();
	for( int i = 0; i < EMDMatrix.rows(); i++)
	{
		if( EMDMatrix(i,0) == -1)
		{
			this->EMDMatrix.set_row(i,(double)0);
			this->EMDMatrix.set_column(i,(double)0);
		}
		else
		{
			moduleForSelection.push_back(i);
		}
	}
	ofs<< this->EMDMatrix<< std::endl;

	for(int i = 0; i < this->EMDMatrix.rows(); i++)
	{
		this->EMDMatrix(i,i) = 1;
		for( int j = 0; j < this->EMDMatrix.rows(); j++)
		{
			if( j != i)
			{
				double max = this->EMDMatrix(i,j) > this->EMDMatrix(j,i) ? this->EMDMatrix(i,j) : this->EMDMatrix(j,i);
				this->EMDMatrix(i,j) = max;
				this->EMDMatrix(j,i) = max;
			}
		}
	}

	ofs<< this->EMDMatrix<< std::endl;
	ofs.close();
	std::cout<< "EMD matrix has been built successfully"<<endl;
}

void SPDAnalysisModel::GetWeightsFromSorting(vnl_matrix< double> &modDist, vnl_vector<int>& index, vnl_vector<double>& weights, unsigned int kNeighbor)
{
	weights.set_size(index.size() - 1);
	for( unsigned int i = 0; i < index.size() - 1; ++i)
	{
		weights[i] = modDist(index[i], index[i+1]);
	}
}

void SPDAnalysisModel::ModuleCorrelationMatrixMatch(unsigned int kNeighbor, int nbins)
{
	m_kNeighbor = kNeighbor;
	std::cout<< m_kNeighbor<<std::endl;
	std::ofstream ofs("ModuleCorrelationMatrixMatch.txt");
	int size = ClusterIndex.max_value() + 1;
	std::vector< std::vector< unsigned int> > featureClusterIndex;
	featureClusterIndex.resize( size);
	for( int i = 0; i < ClusterIndex.size(); i++)
	{
		unsigned int index = ClusterIndex[i];
		featureClusterIndex[index].push_back(i);
	}
	this->EMDMatrix.set_size( featureClusterIndex.size(),  featureClusterIndex.size());
	this->EMDMatrix.fill(-1);

	std::vector< std::vector< unsigned int> > kNearIndex(featureClusterIndex.size());
	vnl_vector< unsigned char> state(featureClusterIndex.size());
	std::set< int> excludeModule;
	state.fill(1);

	omp_lock_t my_lock;
	omp_init_lock(&my_lock);

	#pragma omp parallel for
	for( int i = 0; i < featureClusterIndex.size(); i++)
	{
		vnl_matrix< double> modulei(this->MatrixAfterCellCluster.rows(), featureClusterIndex[i].size());
		GetCombinedMatrix( this->MatrixAfterCellCluster, ClusterIndex, i, i, modulei);

		vnl_vector< double> meanVec( modulei.rows());
		meanVec.fill(0);
		for( unsigned int j = 0; j < modulei.cols(); j++)
		{
			meanVec += modulei.get_column(j);
		}
		meanVec = meanVec / modulei.cols();
		double mmin = meanVec.min_value();
		double mmax = meanVec.max_value();
		
		if( mmax - mmin > 1e-9)
		{
			bool bstate = true;
			vnl_vector<unsigned int> histMean;
			unsigned int binTest = meanVec.size();
			if( binTest * 0.1 > 10)
			{
				binTest = 10;
			}
			else if( binTest * 0.1 < 3)
			{
				binTest = 3;
			}
			double minterval = (mmax - mmin) / binTest;
			Hist(meanVec, minterval, mmin, histMean, binTest);
			unsigned int tmpMax = histMean.max_value();
			unsigned int maxHisti = (unsigned int)(histMean.sum() * 0.9);   // if the module's full graph edge distribution highly agglomerate at one bin, then discard this module.
			if(tmpMax >= maxHisti)
			{
				std::cout<< "1 Remove module "<<i<<std::endl;
				state[i] = 0;
			}
			else
			{
				unsigned int zeroSum = 0;
				for( unsigned int s = 0; s < histMean.size(); s++)
				{
					if(histMean[s] < 1)
					{
						zeroSum++;
					}
				}
				if( zeroSum >= 2.0 / 3.0 * histMean.size())
				{
					std::cout<< "Do not use module "<<i<<" for iteration.\t"<<std::endl;
					omp_set_lock(&my_lock);
						excludeModule.insert(i);
					omp_unset_lock(&my_lock);
				}
			}

			if(state[i])
			{
				vnl_matrix< double> modDisti;
				EuclideanBlockDist(modulei, modDisti);
				double min = 0;
				double max = 0;
				GetDiagnalMinMax(modDisti, min, max);
				double interval = ( max - min) / nbins;

				std::vector<unsigned int> nearIndex;
				vnl_vector<double> nearWeights;
				FindNearestKSample(modDisti, nearIndex, kNeighbor);
				GetKWeights( modDisti, nearIndex, nearWeights, kNeighbor);
				kNearIndex[i] = nearIndex;
			}
		}
	}

	std::cout<< "Modules excluded: "<< excludeModule.size()<<std::endl;

	/// k nearest neighbor graph match process
	
	for( int i = 0; i < featureClusterIndex.size(); i++)
	{
		if( state[i] && ((excludeModule.size() > 0 && excludeModule.find(i) != excludeModule.end()) || excludeModule.size() == 0))
		{
			vnl_vector< double> row;
			row.set_size(featureClusterIndex.size());
			row.fill(0);

			if( i % 100 == 0)
			{
				std::cout<< "Matching module "<<i<<std::endl; 
			}

			vnl_matrix< double> modulei(this->MatrixAfterCellCluster.rows(), featureClusterIndex[i].size());
			GetCombinedMatrix( this->MatrixAfterCellCluster, ClusterIndex, i, i, modulei);
			vnl_matrix< double> modDisti;
			EuclideanBlockDist(modulei, modDisti);
			double min = 0;
			double max = 0;
			GetDiagnalMinMax(modDisti, min, max);
			double interval = ( max - min) / nbins;
			vnl_vector<unsigned int> histModi;
			vnl_vector<double> knnWeights;
			vnl_vector<unsigned int> histi;
			Hist(modDisti, interval, min, histModi, nbins); 
			GetKWeights( modDisti, kNearIndex[i], knnWeights, kNeighbor);
			Hist( knnWeights, interval, min, histi, nbins);
			vnl_matrix<double> flowMatrix( nbins, nbins);
			double discale = EarthMoverDistance( histModi, histi, flowMatrix, nbins, true);
			
			#pragma omp parallel for
			for( int j = 0; j < featureClusterIndex.size(); j++)
			{
				if( state[j])
				{
					if(i == j)
					{
						row[j] = 1;
					}
					else
					{
						vnl_vector<double> matchWeights;
						vnl_vector<unsigned int> histj;
						
						GetKWeights( modDisti, kNearIndex[j], matchWeights, kNeighbor);
						Hist( matchWeights, interval, min, histj, nbins);
						vnl_matrix<double> flowMatrix( nbins, nbins);
						double earth = EarthMoverDistance( histModi, histj, flowMatrix, nbins, true);
						row[j] = earth > 0 ? earth : 0;
						row[j] = row[j] / discale;
					}
				}
			}
			this->EMDMatrix.set_row(i, row);
		}
	}

	moduleForSelection.clear();
	for( int i = 0; i < EMDMatrix.rows(); i++)
	{
		if( EMDMatrix(i,0) == -1)
		{
			this->EMDMatrix.set_row(i,(double)0);
			this->EMDMatrix.set_column(i,(double)0);
		}
		else
		{
			moduleForSelection.push_back(i);
		}
	}
	ofs<< this->EMDMatrix<< std::endl;

	for(int i = 0; i < this->EMDMatrix.rows(); i++)
	{
		this->EMDMatrix(i,i) = 1;
		for( int j = 0; j < this->EMDMatrix.rows(); j++)
		{
			if( j != i)
			{
				double max = this->EMDMatrix(i,j) > this->EMDMatrix(j,i) ? this->EMDMatrix(i,j) : this->EMDMatrix(j,i);
				this->EMDMatrix(i,j) = max;
				this->EMDMatrix(j,i) = max;
			}
		}
	}

	//this->EMDMatrix = this->EMDMatrix / this->EMDMatrix.max_value();

	//ofs<< kNeighbor<<"\t"<<max1<<std::endl;

	ofs<< this->EMDMatrix<< std::endl;
	ofs.close();
	std::cout<< "EMD matrix has been built successfully"<<endl;
}

bool SPDAnalysisModel::SearchSubsetsOfFeatures(std::vector< unsigned int> &selModules)   // feature id 
{
	if( disCor > 0)
	{
		if(selModules.size() == 0)  // search all modules
		{
			for( unsigned int i = 0; i <= ClusterIndex.max_value(); i++)
			{
				selModules.push_back(i);
			}
		}
		std::ofstream ofs("SubsetsOfFeatures_hist_tmp.txt");
        vnl_matrix< double> distanceMat;
		vnl_vector< double> disVec = DistanceToDevice;

		EuclideanBlockDist(disVec, distanceMat);
		double min = distanceMat.min_value();
		std::vector< unsigned int> distancekNearIndex;
		FindNearestKSample(distanceMat, distancekNearIndex, m_kNeighbor);

		for( int k = 0; k < distanceMat.rows(); k++)
		{
			distanceMat(k,k) = 0;
		}
		double max = distanceMat.max_value();
		std::vector< unsigned int> distancekFarIndex;
		FindFarthestKSample(distanceMat, distancekFarIndex, m_kNeighbor);
        
        double interval = ( max - min) / NUM_BIN;
        
        vnl_vector<unsigned int> histDist;
        vnl_vector<unsigned int> histDist2;
        vnl_vector<double> distWeights;
        vnl_vector<double> distWeights2;
       
        GetKWeights( distanceMat, distancekNearIndex, distWeights, m_kNeighbor);
        GetKWeights( distanceMat, distancekFarIndex, distWeights2, m_kNeighbor);
        Hist( distWeights, interval, min, histDist);
        Hist( distWeights2, interval, min, histDist2);
        vnl_matrix<double> flowMatrix( NUM_BIN, NUM_BIN);
        double disScale = EarthMoverDistance( histDist2, histDist, flowMatrix);
		ofs << disScale<<std::endl;
		int size = ClusterIndex.max_value() + 1;
		std::vector< std::vector< unsigned int> > featureClusterIndex;
		featureClusterIndex.resize( size);
		for( int i = 0; i < ClusterIndex.size(); i++)
		{
			unsigned int index = ClusterIndex[i];
			featureClusterIndex[index].push_back(i);
		}

        int searchSize = (int)pow( 2.0, (double)selModules.size());
        vnl_vector<double> fitVec(searchSize - 1);
		std::cout<< searchSize - 1 <<std::endl;

		clock_t start_time = clock();
#pragma omp parallel for
        for( int i = 0; i < searchSize - 1; i++)
        {
			if( i % 10000 == 0)
			{
				std::cout<<i<<std::endl;
			}
            std::vector<unsigned int> modIndex;
            unsigned int count = 0;
            unsigned int tmpi = i;
            while(count <= i / 2)
            {
                unsigned char bit = tmpi % 2;
                if( bit == 1)
                {
                       modIndex.push_back( selModules[count]);
                }
                tmpi = tmpi >> 1;
                count++;
            }

            vnl_matrix<double> comMat;
            GetCombinedMatrixByModuleId( this->MatrixAfterCellCluster, featureClusterIndex, modIndex, comMat);
            vnl_matrix< double> modDisti;
            EuclideanBlockDist(comMat, modDisti);
            std::vector< unsigned int> nearIndex;
            FindNearestKSample(modDisti, nearIndex, m_kNeighbor);
            vnl_vector<double> matchWeights;
            GetKWeights( distanceMat, nearIndex, matchWeights, m_kNeighbor);
            vnl_vector<unsigned int> histi;
            Hist( matchWeights, interval, min, histi);
            vnl_matrix<double> flowMatrix( NUM_BIN, NUM_BIN);
            fitVec[i] = EarthMoverDistance( histi, histDist, flowMatrix);
            //fitVec[i] = fitVec[i] / disScale;
        }

		clock_t duration_time = clock() - start_time;
		double hours = duration_time / ((double) CLOCKS_PER_SEC * 3600);
		ofs <<"Duration hours: "<<hours<<std::endl;

		double max1 = fitVec.max_value();
		double min1 = fitVec.min_value();
		double inter = (max1 - min1) / 90;
		ofs << min1<<"\t"<<max1<<"\t"<<inter<<std::endl;
		vnl_vector<unsigned int> histk;
		Hist(fitVec, inter, min1, histk, 90);
		ofs << histk<<std::endl;
	
		for( unsigned int ks = 0; ks < histk[0]; ks++)
		{
			unsigned int min_arg = fitVec.arg_min();
			ofs << min_arg<<"\t"<<fitVec[min_arg]<<std::endl;
			int ct = 0;
			unsigned tmp = min_arg + 1;
			while(ct <= (min_arg + 1) / 2)
			{
				unsigned int bit = tmp % 2;
				if( bit == 1)
				{
					ofs<< selModules[ct]<<"\t";
				}
				tmp = tmp >> 1;
				ct++;
			}
			ofs<<std::endl;
			fitVec[min_arg] = 1e9;
		}
        return true;
	}
	else
    {
        return false;
    }
}

void SPDAnalysisModel::FindNearestKSample(vnl_matrix< double> &modDist, std::vector< unsigned int>& index, unsigned int kNeighbor)
{
	index.resize( modDist.rows() * kNeighbor);
	for( int i = 0; i < modDist.rows(); i++)
	{
		vnl_vector<double> row = modDist.get_row(i);
		double max = row.max_value();
		
		for( int j = 0; j < kNeighbor; j++)
		{
			unsigned int argMin = row.arg_min();
			index[i * kNeighbor + j] = argMin;
			row[argMin] = max;
		}
	}
}

void SPDAnalysisModel::FindFarthestKSample(vnl_matrix< double> &modDist, std::vector< unsigned int>& index, unsigned int kNeighbor)
{
	index.resize( modDist.rows() * kNeighbor);
	for( int i = 0; i < modDist.rows(); i++)
	{
		vnl_vector<double> row = modDist.get_row(i);
		
		for( int j = 0; j < kNeighbor; j++)
		{
			unsigned int argMax = row.arg_max();
			index[i * kNeighbor + j] = argMax;
			row[argMax] = 0;
		}
	}
}

void SPDAnalysisModel::GetKWeights(vnl_matrix< double> &modDist, std::vector< unsigned int>& index, 
										vnl_vector<double>& weights, unsigned int kNeighbor)
{
	weights.set_size( index.size());
	weights.fill(0);

	vnl_matrix< unsigned char> tagConnected(modDist.rows(), modDist.cols());
	tagConnected.fill(0);

	//#pragma omp parallel for
	for( int i = 0; i < index.size(); i++)
	{
		unsigned int nrow = i / kNeighbor;
		unsigned int indexi = index[i];
		if( tagConnected(nrow, indexi) == 0) 
		{
			weights[i] = modDist(nrow, indexi);
			tagConnected(nrow, indexi) = 1;
			tagConnected(indexi, nrow) = 1;
		}
		else
		{
			weights[i] = -1;
		}
	}
}

void SPDAnalysisModel::GetKWeightsComplement(vnl_matrix< double> &modDist, std::vector< unsigned int>& index, 
										vnl_vector<double>& weightsComplement, unsigned int kNeighbor)
{
	vnl_matrix< unsigned char> tagConnected(modDist.rows(), modDist.cols());
	tagConnected.fill(0);

	//#pragma omp parallel for
	for( int i = 0; i < index.size(); i++)
	{
		unsigned int nrow = i / kNeighbor;
		unsigned int indexi = index[i];
		tagConnected(nrow, indexi) = 1;
		tagConnected(indexi, nrow) = 1;
	}
	
	std::vector<double> weights;
	for( unsigned int i = 0; i < modDist.rows(); i++)
	{
		for( unsigned int j = i + 1; j < modDist.cols(); j++)
		{
			if( tagConnected(i,j) == 0)
			{
				weights.push_back(modDist(i,j));
			}
		}
	}
	weightsComplement.set_size( weights.size());
	for( int i = 0; i < weights.size(); i++)
	{
		weightsComplement[i] = weights[i];
	}
}

int SPDAnalysisModel::GetConnectedComponent(std::vector< unsigned int> &selFeatureID, std::vector<int> &component)
{
	vnl_matrix<double> selMat;
	GetCombinedMatrix( this->MatrixAfterCellCluster, 0, this->MatrixAfterCellCluster.rows(), selFeatureID, selMat);
	vnl_matrix< double> distMat;
	EuclideanBlockDist(selMat, distMat);
	std::vector<unsigned int> nearIndex;
	FindNearestKSample(distMat, nearIndex, m_kNeighbor);

	component.clear();
	int connectedNum = GetConnectedComponent(nearIndex, component, m_kNeighbor);
	return connectedNum;
}

void SPDAnalysisModel::WriteKNNGConnectionMatrix(const char *filename, std::vector< unsigned int> selFeatureID)
{
	vnl_matrix<double> selMat;
	GetCombinedMatrix( this->MatrixAfterCellCluster, 0, this->MatrixAfterCellCluster.rows(), selFeatureID, selMat);
	/*vnl_matrix< double> distMat;
	EuclideanBlockDist(selMat, distMat);
	std::vector<unsigned int> nearIndex;
	FindNearestKSample(distMat, nearIndex, m_kNeighbor);
	vnl_matrix<int> connectedMatrix(this->MatrixAfterCellCluster.rows(),this->MatrixAfterCellCluster.rows());
	connectedMatrix.fill(0);
	for( int i = 0; i < nearIndex.size() / m_kNeighbor; i++)
	{
		for(int j = 0; j < m_kNeighbor; j++)
		{
			int index = nearIndex[i * m_kNeighbor + j];
			connectedMatrix(i,index) = 1;
			connectedMatrix(index,i) = 1;
		}
	}*/

	std::ofstream ofs(filename);
	unsigned int lastCol = selMat.cols() - 1;
	for( unsigned int i = 0; i < selMat.rows(); i++)
	{
		for( unsigned int j = 0; j < selMat.cols() - 1; j++)
		{
			ofs<< selMat(i,j)<<",";
		}
		ofs<< selMat(i,lastCol)<<std::endl;
	}	
	ofs.close();
}

void SPDAnalysisModel::BuildMSTForConnectedComponent(std::vector< unsigned int> &selFeatureID, std::vector<int> &component, int connectedNum)
{
	vnl_matrix<double> selMat;
	GetCombinedMatrix( this->MatrixAfterCellCluster, 0, this->MatrixAfterCellCluster.rows(), selFeatureID, selMat);
	vnl_matrix< double> distMat;
	EuclideanBlockDist(selMat, distMat);
	std::vector<unsigned int> nearIndex;
	FindNearestKSample(distMat, nearIndex, m_kNeighbor);
	std::vector< std::vector<unsigned int> >neighborIndex;
	neighborIndex.resize( nearIndex.size() / m_kNeighbor);
	for( int i = 0; i < nearIndex.size() / m_kNeighbor; i++)
	{
		std::vector<unsigned int> index;
		index.resize(m_kNeighbor);
		for(int j = 0; j < m_kNeighbor; j++)
		{
			index[j] = nearIndex[ i * m_kNeighbor + j];

		}
		neighborIndex[i] = index;
	}

	std::vector< std::vector<int> > graphVertex;
	std::vector< std::map<int, int> > graphMap;
	vnl_vector< int> size(connectedNum);
	graphVertex.resize(connectedNum);
	graphMap.resize(connectedNum);
	size.fill(0);

	for( int i = 0; i < component.size(); i++)
	{
		int n = component[i];
		graphVertex[n].push_back(i);
		graphMap[n][i] = size[n];
		size[n]++;
	}

	std::vector< std::multimap< double, std::pair<int, int> > > edgeMapVec;
	edgeMapVec.resize(connectedNum);

	//std::ofstream ofs("BuildMSTForConnectedComponent.txt");
	double minEdge = 1e6;
	for( int i = 0; i < graphVertex.size(); i++)
	{
		Graph graph(graphVertex[i].size());
		vnl_matrix< unsigned char> mat;
		mat.set_size(this->MatrixAfterCellCluster.rows(), this->MatrixAfterCellCluster.rows());
		mat.fill(0);
		for( int j = 0; j < graphVertex[i].size(); j++)
		{
			int node = graphVertex[i][j];
			for( int k = 0; k < neighborIndex[node].size(); k++)
			{
				int node2 = neighborIndex[node][k];
				if( mat(node, node2) == 0)
				{
					//ofs << graphMap[i][node]<<"\t"<< graphMap[i][node2]<< "\t"<<distMat(node, node2)<<std::endl;
					boost::add_edge(graphMap[i][node], graphMap[i][node2], distMat(node, node2), graph);
					mat(node, node2) = 1;
					mat(node2, node) = 1;
				}
			}
		}
		
		std::vector< boost::graph_traits< Graph>::vertex_descriptor> vertex(graphVertex[i].size());

		try
		{
			boost::prim_minimum_spanning_tree(graph, &vertex[0]);

		}
		catch(...)
		{
			std::cout<< "MST construction failure!"<<endl;
			exit(111);
		}
		
		std::multimap< double, std::pair<int, int> > edgeMap;
		for( int k = 0; k < vertex.size(); k++)
		{
			int nodeInd = vertex[k];
			if( k != nodeInd)
			{
				std::pair< int, int> pair = std::pair< int, int>( graphVertex[i][k], graphVertex[i][nodeInd]);
				double dis = distMat(graphVertex[i][k], graphVertex[i][nodeInd]);
				if( dis < minEdge)
				{
					minEdge = dis;
				}
				edgeMap.insert( std::pair< double, std::pair<int, int> >(dis, pair));
				//ofs << dis<<"\t"<< graphVertex[i][k]<< "\t"<<graphVertex[i][nodeInd]<<std::endl;
			}
		}
		edgeMapVec[i] = edgeMap;
	}

	//ofs <<std::endl;

	int clusNo = MatrixAfterCellCluster.rows();
	std::vector< int> parentIndex;	// save tree root index for each node
	parentIndex.resize(clusNo);
	for( int i = 0; i < clusNo; i++)
	{
		parentIndex[i] = i;
	}
	
	mstTreeList.clear();
	vnl_vector<int> parent(edgeMapVec.size()); // save the parent tree index for each connected graph
	for( int i = 0; i < edgeMapVec.size(); i++)
	{
		std::multimap< double, std::pair<int, int> > map = edgeMapVec[i];
		std::multimap< double, std::pair<int, int> >::iterator iter;
		for( iter = map.begin(); iter != map.end(); iter++)
		{
			std::pair< int, int> pair = iter->second;
			int fir = pair.first;
			int sec = pair.second;
			Tree tree(parentIndex[fir], parentIndex[sec], iter->first - minEdge, clusNo);
			mstTreeList.push_back( tree);
			int firstIndex = parentIndex[fir];
			int secIndex = parentIndex[sec];
			//ofs << iter->first<<"\t"<< parentIndex[fir]<< "\t"<<parentIndex[sec]<<"\t"<<clusNo<<std::endl;
			for( int j = 0; j < parentIndex.size(); j++)
			{
				if( parentIndex[j] == firstIndex || parentIndex[j] == secIndex)
				{
					parentIndex[j] = clusNo;
				}
			}
			clusNo++;
		}
		parent[i] = clusNo - 1;
	}
	//ofs << parentIndex <<std::endl;
	//ofs <<std::endl;

	vnl_matrix< double> subTreeDistance;
	GetComponentMinDistance(distMat, component, connectedNum, subTreeDistance);
	//ofs<< subTreeDistance<<std::endl;
	while(mstTreeList.size() < MatrixAfterCellCluster.rows() - 1)  
	{
		unsigned int location = subTreeDistance.arg_min();
		int nrow = location / subTreeDistance.rows();
		int ncol = location % subTreeDistance.rows();
		if( parent[nrow] != parent[ncol])
		{
			Tree tree( parent[nrow], parent[ncol], subTreeDistance(nrow, ncol) - minEdge, clusNo);
			//ofs << subTreeDistance(nrow, ncol)<<"\t"<< parent[nrow]<< "\t"<< parent[ncol]<<"\t"<<clusNo<<std::endl;
			mstTreeList.push_back(tree);
			subTreeDistance(ncol, nrow)= 1e9;
			subTreeDistance(nrow, ncol)= 1e9;

			for( int i = 0; i < parent.size(); i++)
			{
				if( parent[i] == parent[nrow] || parent[i] == parent[ncol])
				{
					parent[i] = clusNo;
				}
			}
			clusNo++;
		}
		else
		{
			subTreeDistance(ncol, nrow)= 1e9;
			subTreeDistance(nrow, ncol)= 1e9;
		}
	}
	//ofs.close();
}

void SPDAnalysisModel::GetComponentMinDistance( vnl_matrix<double> &distMat, std::vector<int> &component, int connectedNum, vnl_matrix<double> &dis)
{
	dis.set_size(connectedNum, connectedNum);
	dis.fill(1e9);
	
	std::vector< std::vector<int> > graphVertex;
	graphVertex.resize(connectedNum);

	for( int i = 0; i < component.size(); i++)
	{
		int n = component[i];
		graphVertex[n].push_back(i);
	}

	for( int i = 0; i < graphVertex.size(); i++)
	{
		for( int j = i + 1; j < graphVertex.size(); j++)
		{
			double min = FindMinBetweenConnectComponent(distMat, graphVertex[i], graphVertex[j]);
			dis( i, j) = min;
			dis( j, i) = min;
		}
	}
}

void SPDAnalysisModel::GetComponentMinDistance(std::vector< unsigned int> selFeatureID, std::vector<int> &component, int connectedNum, vnl_matrix<double> &dis)
{
	vnl_matrix<double> selMat;
	GetCombinedMatrix( this->MatrixAfterCellCluster, 0, this->MatrixAfterCellCluster.rows(), selFeatureID, selMat);
	vnl_matrix< double> distMat;
	EuclideanBlockDist(selMat, distMat);

	dis.set_size(connectedNum, connectedNum);
	dis.fill(1e9);
	
	std::vector< std::vector<int> > graphVertex;
	graphVertex.resize(connectedNum);

	for( int i = 0; i < component.size(); i++)
	{
		int n = component[i];
		graphVertex[n].push_back(i);
	}

	for( int i = 0; i < graphVertex.size(); i++)
	{
		for( int j = i + 1; j < graphVertex.size(); j++)
		{
			double min = FindMinBetweenConnectComponent(distMat, graphVertex[i], graphVertex[j]);
			dis( i, j) = min;
			dis( j, i) = min;
		}
	}
}

double SPDAnalysisModel::FindMinBetweenConnectComponent(vnl_matrix<double> &dis, std::vector<int> ver1, std::vector<int> ver2)
{
	double min = 1e9;
	for( int i = 0; i < ver1.size(); i++)
	{
		int m = ver1[i];
		for( int j = 0; j < ver2.size(); j++)
		{
			int n = ver2[j];
			if( dis(m, n) < min)
			{
				min = dis(m,n);
			}
		}
	}
	return min;
}

int SPDAnalysisModel::GetConnectedComponent(std::vector< unsigned int>& index, std::vector< int>& component, unsigned int kNeighbor)
{
	Graph graph;
	for( int i = 0; i < index.size(); i++)
	{
		int nrow = i / kNeighbor;
		boost::add_edge( nrow, index[i], graph);
	}
	component.resize(boost::num_vertices(graph));
	int num = boost::connected_components(graph, &component[0]);
	return num;
}

void SPDAnalysisModel::GetCombinedInversedMatrix(vnl_matrix<double> &datamat, vnl_vector<unsigned int>& index, unsigned int moduleId, unsigned int moduleInversedId, vnl_matrix<double>& mat)
{
	unsigned int i = 0;
	unsigned int ind = 0;

	for( i = 0; i < index.size(); i++)
	{
		vnl_vector<double> vec = datamat.get_column(i);
		if( index[i] == moduleId)
		{
			mat.set_column( ind++, vec);
		}
		else if( index[i] == moduleInversedId)
		{
			mat.set_column( ind++, -vec);
		}
	}
}

void SPDAnalysisModel::GetClusClusDataForCorMatrix( clusclus* c1, clusclus* c2, double threshold, std::vector< unsigned int> *disModIndex)
{
	QString filenameSM = this->filename + "similarity_cor_matrix.txt";
	std::ofstream ofSimatrix(filenameSM .toStdString().c_str(), std::ofstream::out);
	ofSimatrix.precision(4);
	//ofSimatrix <<"CorMatrix:"<<endl;
	//ofSimatrix << this->CorMatrix<<endl<<endl;
	//ofSimatrix <<"Module Compare CorMatrix:"<<endl;
	//ofSimatrix << this->ModuleCompareCorMatrix<<endl<<endl;
	this->heatmapMatrix.set_size( this->ModuleCompareCorMatrix.rows(), this->ModuleCompareCorMatrix.cols());
	this->heatmapMatrix.fill(0);
	std::vector< unsigned int> simModIndex;

	if( bProgression)
	{
		ofSimatrix<< "Progression over distance"<<endl;
		for(unsigned int i = 0; i < DistanceCorVector.size(); i++)
		{
			if( abs(DistanceCorVector[i]) >= threshold)
			{
				disModIndex->push_back(i);
				ofSimatrix<< i<< "\t"<< DistanceCorVector[i]<<endl;
			}
		}
	}
	else
	{
		ofSimatrix<< "Overal Progression"<<endl;
		for( unsigned int i = 0; i < this->ModuleCompareCorMatrix.cols(); i++)
		{
			ofSimatrix << "matching with module " <<i<<endl;
			for( unsigned int j = 0; j < this->ModuleCompareCorMatrix.rows(); j++)
			{
				if( abs(this->ModuleCompareCorMatrix( i, j)) >= threshold)     // find the modules that support common mst
				{
					simModIndex.push_back(j);
					ofSimatrix << j <<"\t"<< ModuleCompareCorMatrix( i, j)<<endl;;
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

		c1->Initialize( heatmapMatrix.data_array(), this->ModuleCompareCorMatrix.rows(), this->ModuleCompareCorMatrix.cols());
		c1->RunClusClus();
		c1->Transpose();

		c2->Initialize( c1->transposefeatures,c1->num_features, c1->num_samples);
		c2->RunClusClus();
	}
	ofSimatrix.close();	
}

double SPDAnalysisModel::GetCorMatSelectedPercentage(double thres)
{
	unsigned int count = 0;
	double per = 0;
	if( false == bProgression)
	{
		for( unsigned int i = 0; i < this->ModuleCompareCorMatrix.cols(); i++)
		{
			for( unsigned int j = 0; j < this->ModuleCompareCorMatrix.rows(); j++)
			{
				if( abs(this->ModuleCompareCorMatrix( i, j)) >= thres)     // find the modules that support common mst
				{
					count++;
				}
			}
		}
		unsigned int allnum = this->ModuleCompareCorMatrix.cols() * this->ModuleCompareCorMatrix.rows();
		per = (double)count / allnum;
	}
	else
	{
		for( unsigned int i = 0; i < this->DistanceCorVector.size(); i++)
		{
			if( this->DistanceCorVector[i] >= thres)     // find the modules that support common mst
			{
				count++;
			}
		}
		per = (double) count / this->DistanceCorVector.size();
	}
	return per;
}

void SPDAnalysisModel::SetMaxVertexID(int verId)
{
	maxVertexId = verId;
}

void SPDAnalysisModel::GetPercentage(std::vector< std::vector< long int> > &clusIndex, std::vector< double> &colorVec)
{
	colorVec.clear();
	for( int i = 0; i < clusIndex.size(); i++)
	{
		int count = 0;
		for( int j = 0; j < clusIndex[i].size(); j++)
		{
			if( clusIndex[i][j] <= maxVertexId)
			{
				count++;
			}
		}
		colorVec.push_back( (double)count / clusIndex[i].size());
	}
}

void SPDAnalysisModel::GetCloseToDevicePercentage( std::vector< std::vector< long int> > &clusIndex, std::vector< double> &disPer, double disThreshold)
{
	disPer.clear();
	//std::ofstream ofs("Distance.txt");
	if( disCor > 1)
	{
		for( int i = 0; i < clusIndex.size(); i++)
		{
			//ofs<< "Cluster "<<i<<"\t"<<clusIndex[i].size()<<std::endl;
			int count = 0;
			int distanceCount = 0;
			for( int j = 0; j < clusIndex[i].size(); j++)
			{
				long int index = clusIndex[i][j];
				if( index <= maxVertexId)
				{
					count++;
					std::map< int, int>::iterator iter = indMapFromVertexToClus.find( index);  
					//ofs<< UNDistanceToDevice[iter->second]<<"\t";
					if( iter != indMapFromVertexToClus.end() && UNDistanceToDevice[iter->second] <= disThreshold)
					{
						distanceCount++;
					}
				}
			}
			//ofs<<std::endl<<std::endl;
			if( count > 0)
			{
				disPer.push_back( (double)distanceCount / count);
			}
		}
	}
	else
	{
		for( int i = 0; i < clusIndex.size(); i++)
		{
			disPer.push_back( 0);
		}
	}
}

void SPDAnalysisModel::GetClusterFeatureValue(std::vector< std::vector< long int> > &clusIndex, int nfeature, vnl_vector<double> &featureValue, std::string &featureName)
{
	vnl_vector<double> selCol = UNMatrixAfterCellCluster.get_column( nfeature);
	featureValue.set_size(clusIndex.size());
	for(int i = 0; i < clusIndex.size(); i++)
	{
		double averFeature = 0;
		for( int j = 0; j < clusIndex[i].size(); j++)
		{
			averFeature += selCol[clusIndex[i][j] ];
		}
		averFeature /= clusIndex[i].size();
		featureValue[i] = averFeature;
	}

	char* name = this->DataTable->GetColumn(nfeature + 1)->GetName();
	featureName = "Colored by " + std::string(name);
}

vtkSmartPointer<vtkTable> SPDAnalysisModel::GetAverModuleTable(std::vector< std::vector< long int> > &clusIndex, std::vector<long int> &TreeOrder, std::vector< int> &selFeatureOrder, std::vector< int> &unselFeatureOrder)
{
	vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();

	vtkSmartPointer<vtkVariantArray> column;
	column = vtkSmartPointer<vtkVariantArray>::New();
	column->SetName( "Id");
	table->AddColumn(column);
	column = vtkSmartPointer<vtkVariantArray>::New();
	column->SetName( "Progression Order");
	table->AddColumn(column);
	column = vtkSmartPointer<vtkVariantArray>::New();
	column->SetName( "Progression Tag");
	table->AddColumn(column);
	
	unsigned int colNo = DataTable->GetNumberOfColumns();
	for(int i = 0; i < 10; i++)
	{	
		column = vtkSmartPointer<vtkVariantArray>::New();
		column->SetName( DataTable->GetColumn(colNo - 1 - i)->GetName());
		table->AddColumn(column);
	}

	//for(int i = 0; i < selFeatureOrder.size(); i++)
	//{		
	//	column = vtkSmartPointer<vtkVariantArray>::New();
	//	column->SetName( DataTable->GetColumn(selFeatureOrder[i] + 1)->GetName());
	//	table->AddColumn(column);
	//}
	//for(int i = 0; i < unselFeatureOrder.size(); i++)
	//{		
	//	column = vtkSmartPointer<vtkVariantArray>::New();
	//	column->SetName( DataTable->GetColumn(unselFeatureOrder[i] + 1)->GetName());
	//	table->AddColumn(column);
	//}

	int count = 0;
	for( int i = 0; i < TreeOrder.size(); i++)
	{
		long int n = TreeOrder[i];
		for( int j = 0; j < clusIndex[n].size(); j++)
		{
			int id = clusIndex[n][j];
			vtkSmartPointer<vtkVariantArray> DataRow = vtkSmartPointer<vtkVariantArray>::New();
			DataRow->InsertNextValue( indMapFromIndToVertex[id]);
			DataRow->InsertNextValue( count);
			DataRow->InsertNextValue( i);

			for( int k = 0; k < 10; k++)
			{
				DataRow->InsertNextValue( UNMatrixAfterCellCluster(id, colNo - 1 - k));
			}

			//for( int k = 0; k < selFeatureOrder.size(); k++)
			//{
			//	DataRow->InsertNextValue( UNMatrixAfterCellCluster(id, selFeatureOrder[k]));
			//}
			//for( int k = 0; k < unselFeatureOrder.size(); k++)
			//{
			//	DataRow->InsertNextValue( UNMatrixAfterCellCluster(id, unselFeatureOrder[k]));
			//}

			table->InsertNextRow(DataRow);
			count++;
		}
	}

	//ftk::SaveTable("AverTable.txt",table);
	return table;
}

vtkSmartPointer<vtkTable> SPDAnalysisModel::GetTableForHist(std::vector< int> &selFeatureOrder, std::vector< int> &unselFeatureOrder)
{
	vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();

	vtkSmartPointer<vtkVariantArray> column = vtkSmartPointer<vtkVariantArray>::New();
	column->SetName( "Distance To Device");
	table->AddColumn(column);
	
	for(int i = 0; i < selFeatureOrder.size(); i++)
	{		
		column = vtkSmartPointer<vtkVariantArray>::New();
		column->SetName( DataTable->GetColumn(selFeatureOrder[i] + 1)->GetName());
		table->AddColumn(column);
	}
	for(int i = 0; i < unselFeatureOrder.size(); i++)
	{		
		column = vtkSmartPointer<vtkVariantArray>::New();
		column->SetName( DataTable->GetColumn(unselFeatureOrder[i] + 1)->GetName());
		table->AddColumn(column);
	}

	for( int i =0; i < indMapFromIndToVertex.size(); i++)
	{
		if(indMapFromIndToVertex[i] < maxVertexId)
		{
			vtkSmartPointer<vtkVariantArray> DataRow = vtkSmartPointer<vtkVariantArray>::New();

			DataRow->InsertNextValue( UNDistanceToDevice[i]);

			for( int k = 0; k < selFeatureOrder.size(); k++)
			{
				DataRow->InsertNextValue( UNMatrixAfterCellCluster(i, selFeatureOrder[k]));
			}
			for( int k = 0; k < unselFeatureOrder.size(); k++)
			{
				DataRow->InsertNextValue( UNMatrixAfterCellCluster(i, unselFeatureOrder[k]));
			}
			table->InsertNextRow(DataRow);
		}
		else
		{
			break;  // all the left vertex id larger than maxVertexId
		}
	}

	return table;
}

/// For milti-level demo


/// read the feature distribution table
bool SPDAnalysisModel::RunSPDforFeatureDistributionTable(std::string fileName)
{
	std::ifstream file(fileName.c_str(), std::ifstream::in);
	int line = LineNum(fileName.c_str());
	if( !file.is_open() || line == -1)
	{
		return false;
	}

	vtkSmartPointer<vtkTable> table = ftk::LoadTable(fileName);
	if( table != NULL)
	{
		nBinNum = table->GetNumberOfColumns();
		nFeatureSize = table->GetNumberOfRows() / nSampleSize;
		for( int k = 0; k < nFeatureSize; k++)
		{
			vnl_matrix<unsigned int> mat( nSampleSize, table->GetNumberOfColumns());
			vnl_matrix<double> dismat( nSampleSize, nSampleSize);
			for( int i = k * nSampleSize; i < ( k + 1) * nSampleSize; i++)
			{
				for( int j = 0; j < table->GetNumberOfColumns(); j++)
				{
					unsigned int var = table->GetValue(i, j).ToUnsignedInt();
					if( !boost::math::isnan(var))
					{
						mat( i, j) = var;
					}
					else
					{
						mat( i, j) = 0;
					}
				}
			}
			ComputeDistributionDistance(mat, dismat);
			if( k == 0)
			{
				GenerateMST(dismat, true);  // first time clean the 
			}
			else
			{
				GenerateMST(dismat, false);
			}
		}

		this->EMDMatrix.set_size( nFeatureSize, nFeatureSize);
		this->EMDMatrix.fill(0);
		for( int k = 0; k < nFeatureSize; k++)
		{
			vnl_matrix<unsigned int> mat( nSampleSize, table->GetNumberOfColumns());
			vnl_vector<double> moduleDistance( nSampleSize, nSampleSize);
			for( int i = k * nSampleSize; i < ( k + 1) * nSampleSize; i++)
			{
				for( int j = 0; j < table->GetNumberOfColumns(); j++)
				{
					unsigned int var = table->GetValue(i, j).ToUnsignedInt();
					if( !boost::math::isnan(var))
					{
						mat( i, j) = var;
					}
					else
					{
						mat( i, j) = 0;
					}
				}
			}
			ComputeDistributionDistance(mat, moduleDistance);
			RunEMDAnalysis(moduleDistance, k);
		}
	}
	else
	{
		return false;
	}	
	return true;
}

void SPDAnalysisModel::ComputeDistributionDistance(vnl_matrix<unsigned int> &mat, vnl_matrix<double> &dismat)
{
	dismat.set_size(mat.rows(), mat.rows());
	dismat.fill(0);
	for( int i = 0; i < mat.rows(); i++)
	{
		vnl_vector<unsigned int> dis1 = mat.get_row(i);
		for( int j = i + 1; j < mat.rows(); j++)
		{
			vnl_matrix<double> flowMatrix( nBinNum, nBinNum);
			vnl_vector<unsigned int> dis2 = mat.get_row(j);
			dismat( i, j) = EarthMoverDistance( dis1, dis2, flowMatrix);
		}
	}
}

void SPDAnalysisModel::ComputeDistributionDistance(vnl_matrix<unsigned int> &mat, vnl_vector<double> &moduleDistance)
{
	int rows = mat.rows();
	moduleDistance.set_size( rows * (rows - 1) / 2);
	moduleDistance.fill(0);
	int ind = 0;
	for( unsigned int i = 0; i < mat.rows(); i++)
	{
		vnl_vector<unsigned int> dis1 = mat.get_row(i);
		for( unsigned int j = i + 1; j < mat.rows(); j++)
		{
			vnl_matrix<double> flowMatrix( nBinNum, nBinNum);
			vnl_vector<unsigned int> dis2 = mat.get_row(j);
			moduleDistance[ind++] = EarthMoverDistance( dis1, dis2, flowMatrix);
		}
	}
}

/// mat: distance matrix, nSize: Feature Size, how many msts need to be generated
bool SPDAnalysisModel::GenerateMST( vnl_matrix<double> &mat, bool bfirst)
{
	if( mat.rows() != nSampleSize)
	{
		return false;
	}

	if( bfirst)
	{
		this->ModuleMST.clear();
	}

	unsigned int num_nodes = mat.rows();

	Graph graph(num_nodes);

	for( unsigned int k = 0; k < num_nodes; k++)
	{
		for( unsigned int j = k + 1; j < num_nodes; j++)
		{
			boost::add_edge(k, j, mat(k,j), graph);
		}
	}

	std::vector< boost::graph_traits< Graph>::vertex_descriptor> vertex(num_nodes);

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
	
	for( unsigned int i = 0; i < num_nodes; i++)
	{
		if( i != vertex[i])
		{
			double dist = mat(i,vertex[i]);
			vtkSmartPointer<vtkVariantArray> DataRow = vtkSmartPointer<vtkVariantArray>::New();
			DataRow->InsertNextValue(i);
			DataRow->InsertNextValue(vertex[i]);
			DataRow->InsertNextValue(dist);
			table->InsertNextRow(DataRow);
		}
	}
	this->MSTTable.push_back(table);
	return true;
}

void SPDAnalysisModel::RunEMDAnalysis( vnl_vector<double> &moduleDistance, int ind)
{
	if( this->ModuleMST.size() < nFeatureSize)
	{
		std::cout<< "MSTs haven't been built!"<<std::endl;
		return;
	}

	vnl_vector<double> hist_interval;
	vnl_vector<unsigned int> moduleHist;
	Hist( moduleDistance, NUM_BIN, hist_interval, moduleHist);

	for( int j = 0; j < nFeatureSize; j++)
	{	
		std::cout<< "matching module " << ind<< " with MST " << j<<endl;

		vnl_vector<double> mstDistance; 
		GetMSTMatrixDistance( moduleDistance, this->ModuleMST[j], mstDistance);  // get mst j distance from module i's distance matrix
		vnl_vector<unsigned int> mstHist;
		Hist( mstDistance, hist_interval, mstHist); 

		vnl_matrix<double> flowMatrix( NUM_BIN, NUM_BIN);
		this->EMDMatrix( ind, j) = EarthMoverDistance( moduleHist, mstHist, flowMatrix);
	}


	if( ind == nFeatureSize - 1)
	{
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
			if( max != 0)
			{
				this->DistanceEMDVector[ i] = this->DistanceEMDVector[ i] / max;
			}
		}
		std::cout<< "EMD matrix has been built successfully"<<endl;
	}
}

// compute PS distance between two variables
double SPDAnalysisModel::CaculatePS(unsigned int kNeighbor, unsigned int nbins, vnl_vector<double> &vec1, vnl_vector<double> &vec2, bool debug)
{
	if( abs(vec1.max_value() - vec1.min_value()) < 1e-6 ||  abs(vec2.max_value() - vec2.min_value()) < 1e-6)
	{
		return 0;
	}

	vnl_matrix< double> contextDis;
	EuclideanBlockDist(vec1, contextDis);

	std::vector<unsigned int> nearIndex;
	vnl_vector<double> nearWeights;
	FindNearestKSample(contextDis, nearIndex, kNeighbor);
	GetKWeights( contextDis, nearIndex, nearWeights, kNeighbor);
	double min = contextDis.min_value();

	for( int k = 0; k < contextDis.rows(); k++)
	{
		contextDis(k,k) = 0;
	}

	double max = contextDis.max_value();
	double interval = ( max - min) / nbins;
	double range = 0;
	double rangeNoise = 0;

	vnl_matrix< double> dis2;
	std::vector<unsigned int> nearIndex2;
	vnl_vector<double> matchWeights2;
	vnl_vector<unsigned int> hist2;

	EuclideanBlockDist(vec2, dis2);
	FindNearestKSample(dis2, nearIndex2, kNeighbor);
	GetKWeights( contextDis, nearIndex2, matchWeights2, kNeighbor);
	Hist( matchWeights2, interval, min, hist2, nbins);
	
	double movedEarth = 0;
	double ps = 0;

	if( interval > 1e-6)
	{
		vnl_vector<unsigned int> histNear;
		vnl_vector<unsigned int> histFar;
		Hist(nearWeights, interval, min, histNear, nbins);

		vnl_vector<double> noiseVec(vec1.size());
		for( unsigned int i = 0; i < noiseVec.size(); i++)
		{
			noiseVec[i] = gaussrand(0,1);
		}
		vnl_matrix< double> noiseDis;
		EuclideanBlockDist(noiseVec, noiseDis);
		std::vector<unsigned int> noisenearIndex;
		vnl_vector<double> noiseWeights;
		vnl_vector<unsigned int> histNoise;
		FindNearestKSample(noiseDis, noisenearIndex, kNeighbor);
		GetKWeights(contextDis, noisenearIndex, noiseWeights, kNeighbor);
		Hist(noiseWeights, interval, min, histNoise, nbins);
		std::vector<unsigned int> farIndex;
		vnl_vector<double> farWeights;
		FindFarthestKSample(contextDis, farIndex, kNeighbor);
		GetKWeights(contextDis, farIndex, farWeights, kNeighbor);
		Hist(farWeights, interval, min, histFar, nbins);

		vnl_matrix<double> flowMatrix( nbins, nbins);
		movedEarth = EarthMoverDistance( histNoise, hist2, flowMatrix, nbins);
		double sum = 0;
		for( unsigned int n = 0; n < flowMatrix.rows(); n++)
		{
			for( unsigned int m = n + 1; m < flowMatrix.cols(); m++)
			{
				if( flowMatrix(n,m) > 1e-6)
				{
					sum += flowMatrix(n,m);
				}
			}
		}
		//std::cout<< "Movied earth from upper bins to lower bins:"<<sum<<std::endl;
		if( sum >= 0.5)  // between k-NNG and noise-NNG
		{
			rangeNoise = EarthMoverDistance( histNoise, histNear, flowMatrix, nbins);
			ps = movedEarth / rangeNoise;
		}
		else  // between noise-NNG and k-FNG
		{
			rangeNoise = EarthMoverDistance( histFar, histNoise, flowMatrix, nbins);
			ps = movedEarth / rangeNoise;
		}
		if(debug)
		{
			std::ofstream ofs("debug.txt",std::ios_base::app);
			ofs<< "kNNG1:"<<std::endl;
			ofs<< histNear<<std::endl;
			ofs<< "kFNG1:"<<std::endl;
			ofs<< histFar<<std::endl;
			ofs<< "kNNGn:"<<std::endl;
			ofs<<histNoise<<std::endl;
			ofs<< "kNNG2:"<<std::endl;
			ofs<<hist2<<std::endl;
			if( sum >= 0.5)
			{
				ofs<< "PS(Between kNNG1 and kNNGn):"<<std::endl;
			}
			else
			{
				ofs<< "PS(Between kFNG1 and kNNGn):"<<std::endl;
				
			}
			ofs<<ps<<std::endl<<std::endl;
			ofs.close();
		}
	}
	else
	{
		return 0;
	}

	return ps;
}

void SPDAnalysisModel::ModuleCorrelationPSC(unsigned int kNeighbor, int nbins)
{
	m_kNeighbor = kNeighbor;
	std::cout<< "K equals: "<< m_kNeighbor<<std::endl;

	vnl_vector<int> moduleSize = GetModuleSize( this->ClusterIndex);

	this->EMDMatrix.set_size( moduleSize.size(),  moduleSize.size());
	this->EMDMatrix.fill(0);

	std::vector< std::vector< unsigned int> > kNearIndex;
	std::vector< vnl_matrix< double> > disMatrixPtList;
	vnl_vector<double> disScale;
	disMatrixPtList.resize(moduleSize.size());
	kNearIndex.resize(moduleSize.size());
	disScale.set_size(moduleSize.size());
	disScale.fill(0);

	#pragma omp parallel for
	for( int i = 0; i < moduleSize.size(); i++)
	{
		std::cout<< i<<std::endl;
		vnl_matrix< double> modulei(this->MatrixAfterCellCluster.rows(), moduleSize[i]);
		GetCombinedMatrix( this->MatrixAfterCellCluster, ClusterIndex, i, i, modulei);
		vnl_matrix< double> modDisti;
		EuclideanBlockDist(modulei, modDisti);

		std::vector<unsigned int> nearIndex;
		vnl_vector<double> nearWeights;

		FindNearestKSample(modDisti, nearIndex, kNeighbor);
		GetKWeights( modDisti, nearIndex, nearWeights, kNeighbor);
		double min = modDisti.min_value();

		for( int k = 0; k < modDisti.rows(); k++)
		{
			modDisti(k,k) = 0;
		}
		
		double max = modDisti.max_value();
		double interval = ( max - min) / nbins;

		if( interval > 1e-9)
		{
			vnl_vector<unsigned int> histMod, histNear, histComplement;
			Hist(modDisti, interval, min, histMod, nbins);
			unsigned int maxHist = (unsigned int)(histMod.sum() * 0.9);   // if the module's full graph edge distribution highly agglomerate at one bin, then discard this module.
			bool bmatch = true;
			for( int id = 0; id < histMod.size(); id++) 
			{
				if(histMod[id] > maxHist)
				{
					std::cout<<"Remove "<< i<<std::endl;
					bmatch = false;
					break;
				}
			}
			if(bmatch)
			{
				Hist(nearWeights, interval, min, histNear, nbins);
				histComplement = histMod - histNear;
				vnl_matrix<double> flowMatrix( nbins, nbins);
				disScale[i] = EarthMoverDistance( histComplement, histNear, flowMatrix, nbins);
			}
		}

		disMatrixPtList[i] = modDisti;
		kNearIndex[i] = nearIndex;
	}

	//std::ofstream ofs1("match.txt");
	/// k nearest neighbor graph match process
	for( int i = 0; i < moduleSize.size(); i++)
	{
		//std::cout<< "Matching module "<< i<<std::endl;
		double min = disMatrixPtList[i].min_value();
		double max = disMatrixPtList[i].max_value();
		double interval = ( max - min) / nbins;

		vnl_vector< double> row;
		row.set_size(moduleSize.size());
		row.fill(0);
		if( abs(disScale[i]) > 1e-9)
		{
			vnl_vector<unsigned int> histModi;
			Hist(disMatrixPtList[i], interval, min, histModi, nbins);

			//#pragma omp parallel for
			for( int j = 0; j < moduleSize.size(); j++)
			{
				//std::cout<< j<< "\t";
				if( abs(disScale[j]) > 1e-9)
				{
					vnl_vector<double> matchWeights;
					vnl_vector<unsigned int> histj;
					vnl_vector<unsigned int> histjComplement;
					GetKWeights( disMatrixPtList[i], kNearIndex[j], matchWeights, kNeighbor);
					Hist( matchWeights, interval, min, histj, nbins);
					histjComplement = histModi - histj;

					//ofs1 << i<<"\t"<<j<< std::endl;
					//ofs1 << histjComplement<<std::endl;
					//ofs1 << histj<<std::endl;
					vnl_matrix<double> flowMatrix( nbins, nbins);
					double movedEarth = EarthMoverDistance( histjComplement, histj, flowMatrix, nbins);
					row[j] = movedEarth > 0 ? movedEarth : 0;
				}
			}
			std::cout<< "\n";
		}
	}
	//ofs1.close();

	for(int i = 0; i < this->EMDMatrix.rows(); i++)
	{
		this->EMDMatrix(i,i) = 0;
		for( int j = 0; j < this->EMDMatrix.rows(); j++)
		{
			if( j != i)
			{
				double max = this->EMDMatrix(i,j) > this->EMDMatrix(j,i) ? this->EMDMatrix(i,j) : this->EMDMatrix(j,i);
				this->EMDMatrix(i,j) = max;
				this->EMDMatrix(j,i) = max;
			}
		}
	}
	this->EMDMatrix = this->EMDMatrix / this->EMDMatrix.max_value();

	std::cout<< "EMD matrix has been built successfully."<<endl;

	std::ofstream ofs("PSCEMD.txt");
	ofs<< EMDMatrix<<std::endl<<std::endl;
	ofs.close();
}

void SPDAnalysisModel::EMDMatrixIteration()
{
	std::cout<< "Begin PSC iterations."<<endl;
	vnl_matrix< double> newMatrix;
	bool changed = true;
	int count = 0;
	while(changed)
	{
		std::cout<< "Iteration: "<< count<<std::endl;

		changed = PSCIterationRadius2(EMDMatrix, newMatrix);
		count++;
		EMDMatrix = newMatrix;
		if( count >= 1000)
		{
			std::cout<< "Not able to converge"<<std::endl;
			break;
		}
	}
	std::cout<< "Converged after "<<count<<" iterations."<<std::endl;
	std::ofstream ofs("PSCEMDIter.txt");
	ofs<< EMDMatrix<<std::endl<<std::endl;
	ofs.close();
}

bool SPDAnalysisModel::PSCIterationRadius2(vnl_matrix< double> &input, vnl_matrix< double> &output)
{
	bool bchanged = false;
	output.set_size(input.rows(), input.cols());
	output.fill(0);
	for( unsigned int i = 0; i < output.rows(); i++)
	{
		output(i,i) = input(i,i);
	}

	#ifdef _OPENMP
		omp_lock_t lock;
		omp_init_lock(&lock);
	#endif

	#pragma omp parallel for
	for( int i = 0; i < output.rows(); i++)
	{
		for( int j = i + 1; j < output.cols(); j++)
		{
			double val = input(i,j);
			for( int k = 0; k < output.rows(); k++)
			{
				if( k != i && k != j)
				{
					double min = input(i,k) > input(k,j) ? input(k,j) : input(i,k);
					if(min > val)
					{
						bchanged = true;
						omp_set_lock(&lock);
						val = min;
						omp_unset_lock(&lock);
					}
				}
			}
			output(i,j) = val;
			output(j,i) = val;
		}
	}
	#ifdef _OPENMP
		omp_destroy_lock(&lock);
	#endif
	return bchanged;
}

double SPDAnalysisModel::CaculatePSComplement(unsigned int kNeighbor, unsigned int nbins, vnl_vector<double> &vec1, vnl_vector<double> &vec2, bool debug)
{
	if( abs(vec1.max_value() - vec1.min_value()) < 1e-6 ||  abs(vec2.max_value() - vec2.min_value()) < 1e-6)
	{
		return 0;
	}

	vnl_matrix< double> contextDis;
	EuclideanBlockDist(vec1, contextDis);

	std::vector<unsigned int> nearIndex;
	vnl_vector<double> nearWeights;
	
	FindNearestKSample(contextDis, nearIndex, kNeighbor);
	GetKWeights( contextDis, nearIndex, nearWeights, kNeighbor);
	
	double min = contextDis.min_value();

	for( int k = 0; k < contextDis.rows(); k++)
	{
		contextDis(k,k) = 0;
	}
	
	double max = contextDis.max_value();
	double interval = ( max - min) / nbins;

	vnl_vector<unsigned int> histContext;
	vnl_vector<unsigned int> histNear;

	Hist(contextDis, interval, min, histContext, nbins);
	Hist(nearWeights, interval, min, histNear, nbins);

	vnl_matrix< double> dis2;
	std::vector<unsigned int> nearIndex2;
	vnl_vector<double> matchWeights2;
	
	vnl_vector<unsigned int> hist2;
	EuclideanBlockDist(vec2, dis2);
	FindNearestKSample(dis2, nearIndex2, kNeighbor);
	GetKWeights( contextDis, nearIndex2, matchWeights2, kNeighbor);
	Hist( matchWeights2, interval, min, hist2, nbins);

	bool bstate = false;
	unsigned int maxHist = (unsigned int)(histContext.sum() * 0.9);   // if the module's full graph edge distribution highly agglomerate at one bin, then discard this module.
	for( int id = 0; id < histContext.size(); id++) 
	{
		if(histContext[id] > maxHist)
		{	
			bstate = true;
			break;
		}
	}

	if(bstate)
	{
		return 0;
	}
	
	double ps = 0;
	vnl_matrix<double> flowMatrix( nbins, nbins);
	flowMatrix.fill(0);
	double range = EarthMoverDistance( histContext, histNear, flowMatrix, nbins);
	double movedEarth = EarthMoverDistance( histContext, hist2, flowMatrix, nbins);
	double sum = 0;
	
	for( unsigned int n = 0; n < flowMatrix.rows(); n++)
	{
		for( unsigned int m = n; m < flowMatrix.cols(); m++)
		{
			if( flowMatrix(n,m) > 1e-6)
			{
				sum += flowMatrix(n,m);
			}
		}
	}
	//std::cout<< "Moved earth from upper bins to lower bins:"<<sum<<std::endl;
	if( sum >= 0.5)  
	{
		ps = movedEarth / range;
	}
	else  
	{
		ps = -movedEarth / range;
	}
	
	if(debug)
	{
		std::ofstream ofs("debug.txt",std::ios_base::app);
		ofs<< "kNNG1:"<<std::endl;
		ofs<< histNear<<std::endl;
		ofs<< "full graph:"<<std::endl;
		ofs<< histContext<<std::endl;
		ofs<< "kNNG2:"<<std::endl;
		ofs<< hist2<<std::endl;
		ofs<<ps<<std::endl<<std::endl;
		ofs.close();
	}

	return ps;
}

double SPDAnalysisModel::CaculatePSComplementUsingShortestPath(unsigned int kNeighbor, unsigned int nbins, vnl_vector<double> &vec1, vnl_vector<double> &vec2, double ratio, bool debug)
{
	if( abs(vec1.max_value() - vec1.min_value()) < 1e-6 ||  abs(vec2.max_value() - vec2.min_value()) < 1e-6)
	{
		return 0;
	}
	
	vnl_matrix<double> shortestPathDisRatio;
	bool bstate = GetKDistanceMetricRatio(vec1, vec2, shortestPathDisRatio, vec1.size() * 0.1);
	
	if( !bstate)
	{
		std::cout<<"Failed to build shortest path"<<std::endl;
		return 0;
	}
	//ofs<< shortestPathDisRatio<<std::endl<<std::endl;

	vnl_matrix< double> contextDis;
	EuclideanBlockDist(vec1, contextDis);

	std::vector<unsigned int> nearIndex;
	vnl_vector<double> nearWeights;
	vnl_vector<double> nearWeightsComplement;
	
	FindNearestKSample(contextDis, nearIndex, kNeighbor);
	GetKWeights( contextDis, nearIndex, nearWeights, kNeighbor);
	GetKWeightsComplement( contextDis, nearIndex, nearWeightsComplement, kNeighbor);
	
	double min = contextDis.min_value();

	for( int k = 0; k < contextDis.rows(); k++)
	{
		contextDis(k,k) = 0;
	}
	
	double max = contextDis.max_value();
	double interval = ( max - min) / nbins;

	vnl_matrix< double> dis2;
	std::vector<unsigned int> nearIndex2;
	vnl_vector<double> matchWeights2;
	vnl_vector<double> matchWeights2Complement;
	vnl_vector<unsigned int> hist2;
	vnl_vector<unsigned int> hist2Complement;
	vnl_vector<unsigned int> histNear;
	vnl_vector<unsigned int> histComplement;
	
	EuclideanBlockDist(vec2, dis2);
	//ofs<<dis2<<std::endl<<std::endl;
	for( unsigned int i = 0; i < dis2.rows(); i++)
	{
		for( unsigned int j = i + 1; j < dis2.cols(); j++)
		{
			if( shortestPathDisRatio( i, j) > ratio)
			{
				dis2(i,j) = 1e6;
				dis2(j,i) = 1e6;
			}
		}
	}
	//ofs<<dis2<<std::endl<<std::endl;

	FindNearestKSample(dis2, nearIndex2, kNeighbor);
	GetKWeights( contextDis, nearIndex2, matchWeights2, kNeighbor);
	GetKWeightsComplement( contextDis, nearIndex2, matchWeights2Complement, kNeighbor);
	Hist( matchWeights2, interval, min, hist2, nbins);
	Hist( matchWeights2Complement,  interval, min, hist2Complement, nbins);
	Hist(nearWeights, interval, min, histNear, nbins);
	Hist(nearWeightsComplement, interval, min, histComplement, nbins);
	
	double ps = 0;
	vnl_matrix<double> flowMatrix( nbins, nbins);
	flowMatrix.fill(0);
	double range = EarthMoverDistance( histComplement, histNear, flowMatrix, nbins);
	double movedEarth = EarthMoverDistance( hist2Complement, hist2, flowMatrix, nbins);
	double sum = 0;
	
	for( unsigned int n = 0; n < flowMatrix.rows(); n++)
	{
		for( unsigned int m = n; m < flowMatrix.cols(); m++)
		{
			if( flowMatrix(n,m) > 1e-6)
			{
				sum += flowMatrix(n,m);
			}
		}
	}
	std::cout<< "Moved earth from upper bins to lower bins:"<<sum<<std::endl;
	if( sum >= 0.5)  
	{
		ps = movedEarth / range;
	}
	else  
	{
		ps = -movedEarth / range;
	}
	
	if(debug)
	{
		std::ofstream ofs("debug.txt");
		ofs<< "kNNG1:"<<std::endl;
		ofs<< histNear<<std::endl;
		ofs<< "kNNG1Complement:"<<std::endl;
		ofs<< histComplement<<std::endl;
		ofs<< "kNNG2:"<<std::endl;
		ofs<< hist2<<std::endl;
		ofs<< "kNNG2Complement:"<<std::endl;
		ofs<< hist2Complement<<std::endl;
		ofs<<ps<<std::endl<<std::endl;
		ofs.close();
	}

	return ps;
}

bool SPDAnalysisModel::GetKDistanceMetricRatio(vnl_vector<double> &vec1, vnl_vector<double> &vec2, vnl_matrix<double> &shortestPath, unsigned int kNeighbor)
{	
	//std::ofstream ofs("dis.txt");
	int num_nodes = vec1.size();
	vnl_matrix<double> combMat(num_nodes,2);
	vnl_matrix<double> disMat;
	combMat.set_column(0, vec1);
	combMat.set_column(1, vec2);
	EuclideanBlockDist(combMat, disMat);
	//ofs<< disMat<<std::endl;
	
	std::vector<unsigned int> nearIndex;
	FindNearestKSample(disMat, nearIndex, kNeighbor);
	
	double **distance = new double *[num_nodes];
	for( unsigned i = 0; i < num_nodes; i++)
	{
		distance[i] = new double[num_nodes];
	}
	
	Graph graph(num_nodes);

	for( unsigned int k = 0; k < nearIndex.size(); k++)
	{
		unsigned int ind1 = k / kNeighbor;
		unsigned int ind2 = nearIndex[k];
		boost::add_edge(ind1, ind2, disMat(ind1,ind2), graph);
	}

	bool breturn = false;
	try
	{
		breturn = boost::johnson_all_pairs_shortest_paths(graph, distance);
	}
	catch(...)
	{
		std::cout<< "GetKDistanceMetric failed!"<<endl;
		exit(111);
	}
	
	if( breturn)
	{
		shortestPath.set_size(num_nodes, num_nodes);
		shortestPath.fill(1e6);
		for( unsigned int i = 0; i < num_nodes; i++)
		{
			for( unsigned int j = i + 1; j < num_nodes; j++)
			{
				//ofs<< distance[i][j]<<"\t";
				if( disMat(i,j) > 1e-6)
				{
					shortestPath(i,j) = distance[i][j] / disMat(i,j);
					shortestPath(j,i) = distance[j][i] / disMat(j,i);
				}
				else
				{
					shortestPath(i,j) = 0;
					shortestPath(j,i) = 0;
				}
			}
			//ofs<< std::endl;
		}
	}
	//ofs.close();

	for( unsigned i = 0; i < num_nodes; i++)
	{
		delete distance[i];
	}
	delete distance;
	return breturn;
}

double SPDAnalysisModel::CaculatePSAveragebin(unsigned int kNeighbor, unsigned int nbins, vnl_vector<double> &vec1, vnl_vector<double> &vec2, bool debug)
{
	vnl_matrix< double> contextDis;
	vnl_vector<double> binInterval;
	EuclideanBlockDist(vec1, contextDis);
	AvarageBinHistogram(contextDis, nbins, binInterval);

	std::vector<unsigned int> nearIndex;
	vnl_vector<double> nearWeights;
	FindNearestKSample(contextDis, nearIndex, kNeighbor);
	GetKWeights( contextDis, nearIndex, nearWeights, kNeighbor);
	double min = contextDis.min_value();

	for( int k = 0; k < contextDis.rows(); k++)
	{
		contextDis(k,k) = 0;
	}

	double max = contextDis.max_value();
	double interval = ( max - min) / nbins;
	double range = 0;
	double rangeNoise = 0;

	vnl_matrix< double> dis2;
	std::vector<unsigned int> nearIndex2;
	vnl_vector<double> matchWeights2;
	vnl_vector<unsigned int> hist2;

	EuclideanBlockDist(vec2, dis2);
	FindNearestKSample(dis2, nearIndex2, kNeighbor);
	GetKWeights( contextDis, nearIndex2, matchWeights2, kNeighbor);
	AverageHist( matchWeights2, binInterval, hist2);
	
	double movedEarth = 0;
	double ps = 0;

	if( interval > 1e-6)
	{
		vnl_vector<unsigned int> histNear;
		vnl_vector<unsigned int> histFar;
		AverageHist( nearWeights, binInterval, histNear);

		vnl_vector<double> noiseVec(vec1.size());
		for( unsigned int i = 0; i < noiseVec.size(); i++)
		{
			noiseVec[i] = gaussrand(0,1);
		}
		vnl_matrix< double> noiseDis;
		EuclideanBlockDist(noiseVec, noiseDis);
		std::vector<unsigned int> noisenearIndex;
		vnl_vector<double> noiseWeights;
		vnl_vector<unsigned int> histNoise;
		FindNearestKSample(noiseDis, noisenearIndex, kNeighbor);
		GetKWeights(contextDis, noisenearIndex, noiseWeights, kNeighbor);
		AverageHist( noiseWeights, binInterval, histNoise);
		std::vector<unsigned int> farIndex;
		vnl_vector<double> farWeights;
		FindFarthestKSample(contextDis, farIndex, kNeighbor);
		GetKWeights(contextDis, farIndex, farWeights, kNeighbor);
		AverageHist( farWeights, binInterval, histFar);

		vnl_matrix<double> flowMatrix( nbins, nbins);
		movedEarth = EarthMoverDistance( histNoise, hist2, flowMatrix, nbins);
		double sum = 0;
		for( unsigned int n = 0; n < flowMatrix.rows(); n++)
		{
			for( unsigned int m = n + 1; m < flowMatrix.cols(); m++)
			{
				if( flowMatrix(n,m) > 1e-6)
				{
					sum += flowMatrix(n,m);
				}
			}
		}
		std::cout<< "Moved earth from upper bins to lower bins:"<<sum<<std::endl;
		if( sum >= 0.5)  // between k-NNG and noise-NNG
		{
			rangeNoise = EarthMoverDistance( histNoise, histNear, flowMatrix, nbins);
			ps = movedEarth / rangeNoise;
		}
		else  // between noise-NNG and k-FNG
		{
			rangeNoise = EarthMoverDistance( histFar, histNoise, flowMatrix, nbins);
			ps = movedEarth / rangeNoise;
		}
		if(debug)
		{
			std::ofstream ofs("debug.txt", std::ios_base::app);
			ofs<< "kNNG1:"<<std::endl;
			ofs<< histNear<<std::endl;
			ofs<< "kFNG1:"<<std::endl;
			ofs<< histFar<<std::endl;
			ofs<< "kNNGn:"<<std::endl;
			ofs<<histNoise<<std::endl;
			ofs<< "kNNG2:"<<std::endl;
			ofs<<hist2<<std::endl;
			if( sum >= 0.5)
			{
				ofs<< "PS(Between kNNG1 and kNNGn):"<<std::endl;
			}
			else
			{
				ofs<< "PS(Between kFNG1 and kNNGn):"<<std::endl;
				
			}
			ofs<<ps<<std::endl<<std::endl;
			ofs.close();
		}
	}
	else
	{
		return 0;
	}

	return ps;
}

double SPDAnalysisModel::gaussrand(double exp, double std)
{
    static double V1, V2, S;
    static int phase = 0;
    double X;
     
    if ( phase == 0 ) {
        do {
            double U1 = (double)rand() / RAND_MAX;
            double U2 = (double)rand() / RAND_MAX;
             
            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
        } while(S >= 1 || S == 0);
         
        X = V1 * sqrt(-2 * log(S) / S);
    } else
        X = V2 * sqrt(-2 * log(S) / S);
         
    phase = 1 - phase;
 	X = X * std + exp;
    return X;
}

void SPDAnalysisModel::AvarageBinHistogram(vnl_matrix<double> &disMetric, int nbins, vnl_vector<double> &binInterval)
{
	unsigned int rown = disMetric.rows();
	unsigned int coln = disMetric.cols();
	int edgeNum = rown * (coln - 1);
	double *distance = new double[edgeNum];
	int k = 0;
	for( unsigned int i = 0; i < rown; i++)
	{
		for( unsigned int j = 0; j < coln; j++)
		{
			if( i != j)
			{
				distance[k++] = disMetric(i,j);
			}
		}
	}
	quickSort( distance, 0, edgeNum - 1);

	int num = ceil(double(edgeNum) / nbins);
	int binnum = ceil(double(edgeNum) / num);
	binInterval.set_size( binnum + 1);
	int count = 0;
	for( int i = 0; i < edgeNum; i++)
	{
		if( i % num == 0)
		{
			binInterval[count++] = distance[i];
		}
	}
	if( count == binnum)
	{
		binInterval[count] = distance[edgeNum - 1];
	}
	delete distance;
}

void SPDAnalysisModel::AverageHist(vnl_vector<double>&distance, vnl_vector<double>& binInterval, vnl_vector<unsigned int>& histDis)
{
	histDis.set_size(binInterval.size() - 1);
	histDis.fill(0);
	for( unsigned int i = 0; i < distance.size(); i++)
	{
		for( unsigned int j = 0; j < binInterval.size() - 1; j++)
		{
			if( distance[i] >= binInterval[j] && distance[i] < binInterval[j + 1])
			{
				histDis[j]++;
			}
		}
	}
}

void SPDAnalysisModel::quickSort(double *arr, int left, int right, int *id) 
{
      int i = left;
	  int j = right;
      double tmp = 0;
      double pivot = arr[(left + right) / 2];
 
      /* partition */
      while (i <= j) {
            while (arr[i] < pivot)
                  i++;
            while (arr[j] > pivot)
                  j--;
            if (i <= j) {
                  tmp = arr[i];
                  arr[i] = arr[j];
                  arr[j] = tmp;
				
				  if( id)
				  {
					  int idtmp = id[i];
					  id[i] = id[j];
					  id[j] = idtmp;
				  }

                  i++;
                  j--;
            }
      }
 
      /* recursion */
      if (left < j)
            quickSort(arr, left, j, id);
      if (i < right)
            quickSort(arr, i, right, id);
}

double SPDAnalysisModel::SimulateVec2(unsigned int kNeighbor, unsigned int nbins, vnl_vector<double> &vec1, vnl_vector<double> &vec2, bool debug)
{
	if( abs(vec1.max_value() - vec1.min_value()) < 1e-6)
	{
		return 0;
	}

	vnl_matrix< double> contextDis;
	EuclideanBlockDist(vec1, contextDis);

	std::vector<unsigned int> nearIndex;
	vnl_vector<double> nearWeights;
	FindNearestKSample(contextDis, nearIndex, kNeighbor);
	GetKWeights( contextDis, nearIndex, nearWeights, kNeighbor);
	double min = contextDis.min_value();

	for( int k = 0; k < contextDis.rows(); k++)
	{
		contextDis(k,k) = 0;
	}

	double max = contextDis.max_value();
	double interval = ( max - min) / nbins;
	double range = 0;
	double rangeNoise = 0;
	
	double movedEarth = 0;
	double ps = 0;

	if( interval > 1e-6)
	{
		vnl_vector<unsigned int> histNear;
		vnl_vector<unsigned int> histFar;
		Hist(nearWeights, interval, min, histNear, nbins);

		vnl_vector<double> noiseVec(vec1.size());
		for( unsigned int i = 0; i < noiseVec.size(); i++)
		{
			noiseVec[i] = gaussrand(0,1);
		}
		vnl_matrix< double> noiseDis;
		EuclideanBlockDist(noiseVec, noiseDis);
		std::vector<unsigned int> noisenearIndex;
		vnl_vector<double> noiseWeights;
		vnl_vector<unsigned int> histNoise;
		FindNearestKSample(noiseDis, noisenearIndex, kNeighbor);
		GetKWeights(contextDis, noisenearIndex, noiseWeights, kNeighbor);
		Hist(noiseWeights, interval, min, histNoise, nbins);
		std::vector<unsigned int> farIndex;
		vnl_vector<double> farWeights;
		FindFarthestKSample(contextDis, farIndex, kNeighbor);
		GetKWeights(contextDis, farIndex, farWeights, kNeighbor);
		Hist(farWeights, interval, min, histFar, nbins);

		vnl_vector<double> newY(vec1.size());
		newY.fill(-1);
		double yNum = 0.5;
		for( unsigned int i = 0; i < vec1.size(); i++)
		{
			bool bfind = false;
			double yN = yNum;
			for( unsigned int k = 0; k < kNeighbor;k++)
			{
				unsigned int ind = i * kNeighbor + k;
				unsigned int fN = farIndex[ind];
				if( newY[i] > 0)
				{
					bfind = true;
					yN = newY[i];
                    break;
				}
				else if( newY[fN] > 0)
				{
                    bfind = true;
                    yN = newY[fN];
                    break;
				}
			}
            if( bfind)
            {
				newY[i] = yN;
                for( unsigned int k = 0; k < kNeighbor;k++)
                {
                        unsigned int ind = i * kNeighbor + k;
                        unsigned int fN = farIndex[ind];
                        newY[fN]=yN;
                }
            }
            else
            {
				newY[i] = yN;
                for( unsigned int k = 0; k < kNeighbor;k++)
                {
                        unsigned int ind = i * kNeighbor + k;
                        unsigned int fN = farIndex[ind];
                        newY[fN]=yN;
                }
                yNum += 0.5;
            }
		}

	vnl_matrix< double> dis2;
	std::vector<unsigned int> nearIndex2;
	vnl_vector<double> matchWeights2;
	vnl_vector<unsigned int> hist2;

	EuclideanBlockDist(newY, dis2);
	FindNearestKSample(dis2, nearIndex2, kNeighbor);
	GetKWeights( contextDis, nearIndex2, matchWeights2, kNeighbor);
	Hist( matchWeights2, interval, min, hist2, nbins);

		vnl_matrix<double> flowMatrix( nbins, nbins);
		movedEarth = EarthMoverDistance( histNoise, hist2, flowMatrix, nbins);
		double sum = 0;
		for( unsigned int n = 0; n < flowMatrix.rows(); n++)
		{
			for( unsigned int m = n + 1; m < flowMatrix.cols(); m++)
			{
				if( flowMatrix(n,m) > 1e-6)
				{
					sum += flowMatrix(n,m);
				}
			}
		}
		std::cout<< "Movied earth from upper bins to lower bins:"<<sum<<std::endl;
		if( sum >= 0.5)  // between k-NNG and noise-NNG
		{
			rangeNoise = EarthMoverDistance( histNoise, histNear, flowMatrix, nbins);
			ps = movedEarth / rangeNoise;
		}
		else  // between noise-NNG and k-FNG
		{
			rangeNoise = EarthMoverDistance( histFar, histNoise, flowMatrix, nbins);
			ps = movedEarth / rangeNoise;
		}

                if(debug)
                {
                    std::ofstream ofs("newY.txt");
                    ofs<< "x:"<<std::endl;
                    ofs<< vec1<<std::endl;
                    ofs<< "y:"<<std::endl;
                    ofs<< newY<<std::endl;
                }
	}
	else
	{
		return 0;
	}

	return ps;
}
