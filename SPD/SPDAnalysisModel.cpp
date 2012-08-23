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
#define DISTANCE_PRECISION 10
#define LMEASURETABLE 1

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

	//table->RemoveColumnByName("centroid_x");
	//table->RemoveColumnByName("centroid_y");
	//table->RemoveColumnByName("centroid_z");

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
		disCor = (DistanceToDevice - disMean).two_norm();
		std::cout<< "Distance cor: "<<disCor<<endl;
#else 
		ConvertTableToMatrixForLayerData(this->DataTable, this->DataMatrix, this->indMapFromIndToVertex, clusNo);
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

	table->RemoveColumnByName("Distance to Device");
	table->RemoveColumnByName("Distance_to_Device");

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

void SPDAnalysisModel::ConvertTableToMatrixForLayerData(vtkSmartPointer<vtkTable> table, vnl_matrix<double> &mat, 
        std::vector<int> &index, vnl_vector<int> &clusNo)
{
	clusNo.set_size( table->GetNumberOfRows());
	for( int i = 0; i < clusNo.size(); i++)
	{
		clusNo[i] = 0;
	}

	vtkVariantArray *distanceArray = vtkVariantArray::SafeDownCast(table->GetColumnByName("prediction_active"));
	if( distanceArray == NULL)
	{
		distanceArray = vtkVariantArray::SafeDownCast(table->GetColumnByName("prediction_active"));
	}

	vtkDoubleArray *distanceDoubleArray = vtkDoubleArray::SafeDownCast(table->GetColumnByName("prediction_active"));
	if( distanceDoubleArray == NULL)
	{
		distanceDoubleArray = vtkDoubleArray::SafeDownCast(table->GetColumnByName("prediction_active"));
	}

	if( distanceArray)
	{
		for( int i = 0; i < distanceArray->GetNumberOfValues(); i++)
		{
			clusNo[i] = distanceArray->GetValue(i).ToDouble();
		}
	}
	else if( distanceDoubleArray)
	{
		for( int i = 0; i < distanceDoubleArray->GetNumberOfTuples(); i++)
		{
			clusNo[i] = distanceDoubleArray->GetValue(i);
		}
	}

	table->RemoveColumnByName("prediction_active");

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
	for(int i = 0; i < DataTable->GetNumberOfColumns(); i++)
	{		
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName( DataTable->GetColumn(i)->GetName());
		table->AddColumn(column);
	}
		
	for( int i = 0; i < mat.rows(); i++)
	{
		vtkSmartPointer<vtkVariantArray> DataRow = vtkSmartPointer<vtkVariantArray>::New();
		DataRow->InsertNextValue(i);
		for( int j = 0; j < mat.cols(); j++)
		{
			DataRow->InsertNextValue( mat(i,j));
		}
		//DataRow->InsertNextValue(distance[i]);
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
	ofs<<"first clustering:"<<endl;

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

	ofs<<"Left cluster number:"<<endl;
	do
	{
		old_cluster_num = new_cluster_num;
		ofs<<old_cluster_num<<endl;
		new_cluster_num = ClusterAggFeatures( MatrixAfterCellCluster, clusterIndex, moduleMean, cor);
		//ofs<<clusterIndex<<endl<<endl;
	}
	while( old_cluster_num != new_cluster_num);
	ofs<<endl<<endl;

	this->ClusterIndex = clusterIndex;
	this->ModuleMean = moduleMean;

	for( unsigned int i = 0; i <= this->ClusterIndex.max_value(); i++)
	{
		ofs<<"Cluster:"<< i<<endl;
		for( unsigned int j = 0; j < this->ClusterIndex.size(); j++)
		{
			if( ClusterIndex[j] == i)
			{
				char* name = this->DataTable->GetColumn(j + 1)->GetName();
				std::string featureName(name);
				ofs<<featureName<<endl;
			}
		}
	}
	ofs<<endl;
	ofs.close();

	std::cout<<"The cluster size after cluster: "<<this->ClusterIndex.max_value() + 1<<endl;
	std::cout<<clusterIndex<<endl;

	return new_cluster_num;
}


int SPDAnalysisModel::ClusterAggFeatures(vnl_matrix<double>& mainmatrix, vnl_vector<unsigned int>& index, vnl_matrix<double>& mean, double cor)
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
	vnl_vector<int> tmp(moduleSize.size());
	while( isActiveModule.sum() > 1)
	{
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

		newModule.set_size(mainmatrix.rows(), newModuleSize);

		if( moduleToDeleteId != moduleId && isActiveModule[moduleToDeleteId] == 1)
		{
			GetCombinedMatrix( mainmatrix, index, moduleId, moduleToDeleteId, newModule);   
			newModule.normalize_columns();

			vnl_matrix<double> mat = newModule.transpose();
			GetMatrixRowMeanStd(mat, newModuleMeans, newModuleStd);

			newModuleMeans = newModuleMeans.normalize();
			newModuleCor = newModuleMeans * newModule;
			double newCor = newModuleCor.mean();                       

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

		isActiveModule[moduleId] = 0;
	}

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
	std::cout<<this->ClusterIndex<<endl<<endl;

	QString filenameCluster = this->filename + "clustering.txt";
	std::ofstream ofs( filenameCluster.toStdString().c_str(), std::ofstream::out | ios::app);
	ofs<< "clustering result after merge:"<<endl;
	for( unsigned int i = 0; i <= this->ClusterIndex.max_value(); i++)
	{
		ofs<<"Cluster: "<<i<<endl;
		for( unsigned int j = 0; j < this->ClusterIndex.size(); j++)
		{
			if( ClusterIndex[j] == i)
			{
				char* name = this->DataTable->GetColumn(j + 1)->GetName();
				std::string featureName(name);
				ofs<< j<<"\t"<<featureName<<endl;
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

vtkSmartPointer<vtkTable> SPDAnalysisModel::GenerateMST( vnl_matrix<double> &mat, std::vector< unsigned int> &selFeatures)
{
	NormalizeData(mat);
	vnl_matrix<double> clusterMat( mat.rows(), selFeatures.size());
	GetCombinedMatrix( mat, selFeatures, clusterMat);

	int num_nodes = clusterMat.rows();
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

	vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();

	for(int i = 0; i < this->headers.size(); i++)
	{		
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName( (this->headers)[i].c_str());
		table->AddColumn(column);
	}

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
		}
	}
	return table;
}

void SPDAnalysisModel::GetCombinedMatrix( vnl_matrix<double> &datamat, std::vector< unsigned int> selFeatureIDs, vnl_matrix<double>& mat)
{
	int count = 0;
	for( unsigned int i = 0; i < datamat.cols(); i++)
	{
		if( IsExist(selFeatureIDs, i))
		{
			vnl_vector<double> vec = datamat.get_column(i);
			mat.set_column(count, vec);
			count++;
		}
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

	#ifdef _OPENMP
		omp_lock_t lock;
		omp_init_lock(&lock);
	#endif

	#pragma omp parallel for
	for( int i = 0; i <= this->ClusterIndex.max_value(); i++)
	{
		vnl_vector<double> distance;
		std::cout<< "build module distance "<< i<<endl;
		vnl_matrix<double> clusterData( this->MatrixAfterCellCluster.rows(), moduleSize[i]);
		GetCombinedMatrix( this->MatrixAfterCellCluster, this->ClusterIndex, i, i, clusterData);    /// Get the module data for cluster i
		GetMatrixDistance( clusterData, distance, CITY_BLOCK);
		#ifdef _OPENMP
			omp_set_lock(&lock);
		#endif
			moduleDistance[i] = distance;
		#ifdef _OPENMP
			omp_unset_lock(&lock);
		#endif
	}
	#ifdef _OPENMP
		omp_destroy_lock(&lock);
	#endif

	this->EMDMatrix.set_size( this->ClusterIndex.max_value() + 1, this->ClusterIndex.max_value() + 1);
	this->EMDMatrix.fill(0);
	this->DistanceEMDVector.set_size( this->ClusterIndex.max_value() + 1);
	this->DistanceEMDVector.fill(0);

	//#pragma omp parallel for
	for( int i = 0; i < moduleDistance.size(); i++)
	{
		vnl_vector<double> hist_interval;
		vnl_vector<unsigned int> moduleHist;
		Hist( moduleDistance[i], NUM_BIN, hist_interval, moduleHist);

		for( int j = 0; j < moduleDistance.size(); j++)
		{	
			std::cout<< "matching module " << i<< " with MST " << j<<endl;

			vnl_vector<double> mstDistance; 
			GetMSTMatrixDistance( moduleDistance[i], this->ModuleMST[j], mstDistance);  // get mst j distance from module i's distance matrix
			vnl_vector<unsigned int> mstHist;
			Hist( mstDistance, hist_interval, mstHist); 

			vnl_matrix<double> flowMatrix( NUM_BIN, NUM_BIN);
			this->EMDMatrix( i, j) = EarthMoverDistance( moduleHist, mstHist, flowMatrix);
		}

		/// matching with distance to device
		if( disCor != 0)
		{
			std::cout<< "matching module " << i<< " with distance MST "<<endl;
			vnl_vector<double> mstDistance; 
			GetMSTMatrixDistance( moduleDistance[i], this->DistanceMST, mstDistance);  // get mst j distance from module i's distance matrix
			vnl_vector<unsigned int> mstHist;
			Hist( mstDistance, hist_interval, mstHist); 
			vnl_matrix<double> flowMatrix( NUM_BIN, NUM_BIN);
			this->DistanceEMDVector[ i] = EarthMoverDistance( moduleHist, mstHist, flowMatrix);
		}
	}

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

void SPDAnalysisModel::GetEMDMatrixDivByMax(vnl_matrix<double> &emdMatrix)
{
	emdMatrix =  this->EMDMatrix;
}

void SPDAnalysisModel::GetClusClusData(clusclus *c1, clusclus *c2, double threshold, std::vector< unsigned int> *disModIndex)
{
	QString filenameSM = this->filename + "similarity_matrix.txt";
	std::ofstream ofSimatrix(filenameSM .toStdString().c_str(), std::ofstream::out);
	ofSimatrix.precision(4);

	this->heatmapMatrix.set_size( this->EMDMatrix.rows(), this->EMDMatrix.cols());
	this->heatmapMatrix.fill(0);
	std::vector< unsigned int> simModIndex;

	if( bProgression)
	{
		ofSimatrix<< "Progression over distance"<<endl;
		for(unsigned int i = 0; i < DistanceEMDVector.size(); i++)
		{
			if( DistanceEMDVector[i] >= threshold)
			{
				disModIndex->push_back(i);
				ofSimatrix<< i<< "\t"<< DistanceEMDVector[i]<<endl;
			}
		}
	}
	else
	{
		ofSimatrix<< "Overal Progression"<<endl;
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

		c1->Initialize( heatmapMatrix.data_array(), this->EMDMatrix.rows(), this->EMDMatrix.cols());
		c1->RunClusClus();
		c1->Transpose();

		c2->Initialize( c1->transposefeatures,c1->num_features, c1->num_samples);
		c2->RunClusClus();
	}
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
			double dist = CityBlockDist( clusterMat, i, vertex[i]);
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

double SPDAnalysisModel::GetEMDSelectedPercentage(double thres)
{
	unsigned int count = 0;
	double per = 0;
	if( false == bProgression)
	{
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
		per = (double)count / allnum;
	}
	else
	{
		for( unsigned int i = 0; i < this->DistanceEMDVector.size(); i++)
		{
			if( this->DistanceEMDVector[i] >= thres)     // find the modules that support common mst
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
	std::ofstream ofs("Distance.txt");
	if( disCor > 1)
	{
		for( int i = 0; i < clusIndex.size(); i++)
		{
			ofs<< "Cluster "<<i<<"\t"<<clusIndex[i].size()<<std::endl;
			int count = 0;
			int distanceCount = 0;
			for( int j = 0; j < clusIndex[i].size(); j++)
			{
				long int index = clusIndex[i][j];
				if( index <= maxVertexId)
				{
					count++;
					std::map< int, int>::iterator iter = indMapFromVertexToClus.find( index);  
					ofs<< UNDistanceToDevice[iter->second]<<"\t";
					if( iter != indMapFromVertexToClus.end() && UNDistanceToDevice[iter->second] <= disThreshold)
					{
						distanceCount++;
					}
				}
			}
			ofs<<std::endl<<std::endl;
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

vtkSmartPointer<vtkTable> SPDAnalysisModel::GetAverModuleTable(std::vector< std::vector< long int> > &clusIndex, std::vector<long int> &TreeOrder, std::vector< double> &percentageOfSamples,
				    std::vector< double> &percentageOfNearDeviceSamples, std::vector< int> &selFeatureOrder, std::vector< int> &unselFeatureOrder)
{
	vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
	//vnl_matrix<double> averFeatureModule( clusIndex.size(), UNMatrixAfterCellCluster.cols());
	//for( int i = 0; i < clusIndex.size(); i++)
	//{
	//	vnl_vector<double> tmp( UNMatrixAfterCellCluster.cols());
	//	tmp.fill(0);
	//	for( int j = 0; j < clusIndex[i].size(); j++)
	//	{
	//		vnl_vector<double> row = UNMatrixAfterCellCluster.get_row(clusIndex[i][j]);
	//		tmp += row;
	//	}
	//	tmp /= clusIndex[i].size();
	//	averFeatureModule.set_row(i, tmp);
	//}

	//std::ofstream ofs("averFeatureModule.txt");
	//ofs<< averFeatureModule<<endl;
	//vnl_vector<long int> order( TreeOrder.size());
	//for( int i = 0; i < TreeOrder.size(); i++)
	//{
	//	order[ TreeOrder[i]] = i;
	//}
	//ofs<< order<<endl;
	//ofs.close();

	vtkSmartPointer<vtkVariantArray> column;
	column = vtkSmartPointer<vtkVariantArray>::New();
	column->SetName( "Id");
	table->AddColumn(column);
	column = vtkSmartPointer<vtkVariantArray>::New();
	column->SetName( "Progression Order");
	table->AddColumn(column);
	column = vtkSmartPointer<vtkVariantArray>::New();
	column->SetName( "Device Sample Percentage");
	table->AddColumn(column);
	column = vtkSmartPointer<vtkVariantArray>::New();
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
			DataRow->InsertNextValue( percentageOfSamples[i]);
			DataRow->InsertNextValue( UNDistanceToDevice[id]);

			for( int k = 0; k < selFeatureOrder.size(); k++)
			{
				DataRow->InsertNextValue( UNMatrixAfterCellCluster(id, selFeatureOrder[k]));
			}
			for( int k = 0; k < unselFeatureOrder.size(); k++)
			{
				DataRow->InsertNextValue( UNMatrixAfterCellCluster(id, unselFeatureOrder[k]));
			}

			table->InsertNextRow(DataRow);
			count++;
		}
	}
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