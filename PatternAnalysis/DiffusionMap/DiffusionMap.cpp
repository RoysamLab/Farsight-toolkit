#include "DiffusionMap.h"

DiffusionMap::DiffusionMap()
{	
}

DiffusionMap::~DiffusionMap()
{
}


void DiffusionMap::Initialize(vtkSmartPointer< vtkTable > tbl, bool val)
{
	featureTable = tbl;
	multi_cell_struct = val;
	
	vnl_matrix<double> temp;
	if(multi_cell_struct)
	{
		temp.set_size((int)featureTable->GetNumberOfRows(), (int)featureTable->GetNumberOfColumns()-4);
		for(int row=0; row<(int)featureTable->GetNumberOfRows(); ++row)
		{
			for(int col=4; col<(int)featureTable->GetNumberOfColumns(); ++col)
			{
				temp(row, col-4) = featureTable->GetValue(row,col).ToDouble();
			}
			idRowMap[featureTable->GetValue(row,0).ToInt()] = row;
			rowIdMap[row] = featureTable->GetValue(row,0).ToInt();
			std::vector<double> c;
			c.push_back(featureTable->GetValue(row,1).ToDouble());
			c.push_back(featureTable->GetValue(row,2).ToDouble());
			c.push_back(featureTable->GetValue(row,3).ToDouble());
			centroidMap[featureTable->GetValue(row,0).ToInt()] = c;	
		}
	}
	else
	{
		temp.set_size((int)featureTable->GetNumberOfRows(), (int)featureTable->GetNumberOfColumns()-1);
		for(int row=0; row<(int)featureTable->GetNumberOfRows(); ++row)
		{
			for(int col=1; col<(int)featureTable->GetNumberOfColumns(); ++col)
			{
				temp(row, col-1) = featureTable->GetValue(row,col).ToDouble();
			}
			idRowMap[featureTable->GetValue(row,0).ToInt()] = row;
			rowIdMap[row] = featureTable->GetValue(row,0).ToInt();
		}
	}

	mbl_stats_nd stats;
	
	for(int i = 0; i<temp.rows() ; ++i)
	{
		vnl_vector<double> temp_row = temp.get_row(i);
		stats.obs(temp_row);	
	}

	vnl_vector<double> std_vec = stats.sd();
	vnl_vector<double> mean_vec = stats.mean();
	
	data_matrix.set_size((int)featureTable->GetNumberOfRows(), (int)featureTable->GetNumberOfColumns()-1);

	for(int i = 0; i<temp.columns() ; ++i)
	{
		vnl_vector<double> temp_col = temp.get_column(i);
		if(std_vec(i) > 0)
		{	
			for(int j =0; j<temp_col.size() ; ++j)
				temp_col[j] = (temp_col[j] - mean_vec(i))/std_vec(i) ;
		}
		data_matrix.set_column(i,temp_col);
	}
}

void DiffusionMap::ComputeDiffusionMap(void)
{
	std::cout << "Computing Random Walk matrix...";
	similarity_matrix.set_size(data_matrix.rows(), data_matrix.rows());
	vnl_vector<double> weight_vector(data_matrix.rows(), 0);

	if(multi_cell_struct)
	{
		std::cout << "I'm here\n";
		similarity_matrix.fill(0);
		kNearestObjects<3>* KNObj = new kNearestObjects<3>(centroidMap);
		//KNObj->setFeatureTable(table);
		std::vector<std::vector< std::pair<unsigned int, double> > > kNeighborIDs;
		kNeighborIDs = KNObj->k_nearest_neighbors_All(7, 0, 0);

		//#pragma omp parallel for
		for(int i=0; i<kNeighborIDs.size(); ++i)
		{
			int row = idRowMap[kNeighborIDs[i][0].first];
			#pragma omp parallel for
			for(int j=1; j<kNeighborIDs[i].size(); ++j)
			{
				int col = idRowMap[kNeighborIDs[i][j].first];
				if(similarity_matrix(row,col) == 0)
				{
					vnl_vector<double> dt1 = data_matrix.get_row(row);
					vnl_vector<double> dt2 = data_matrix.get_row(col);
					similarity_matrix(row,col) = exp( -((dt1-dt2).squared_magnitude()) / (2*pow((double)2/*sigma*/, 2)) );
					weight_vector(row) += similarity_matrix(row,col);
					similarity_matrix(col,row) = similarity_matrix(row,col);
					weight_vector(col) += similarity_matrix(col,row);
				}
			}
		}
	}
	else
	{
		#pragma omp parallel for
		for(int row=0; row<(int)data_matrix.rows(); ++row)
		{
			for(int col=0; col<(int)data_matrix.rows(); ++col)
			{
				if(row == col)
					similarity_matrix(row,col) = 0;
				else if(row > col)
				{
					similarity_matrix(row,col) = similarity_matrix(col,row);
					weight_vector(row) += similarity_matrix(row,col);
				}
				else
				{
					vnl_vector<double> dt1 = data_matrix.get_row(row);
					vnl_vector<double> dt2 = data_matrix.get_row(col);
					similarity_matrix(row,col) = exp( -((dt1-dt2).squared_magnitude()) / (2*pow((double)2/*sigma*/, 2)) );
					weight_vector(row) += similarity_matrix(row,col);
				}
			}
		}
	}
	
	ofstream outFile; 
	outFile.open("mmmmmmmmmmmmmmmm.txt", ios::out | ios::trunc );
	random_walk_matrix.set_size(data_matrix.rows(), data_matrix.rows());
	for(int row=0; row<(int)similarity_matrix.rows(); ++row)
	{
		for(int col=0; col<(int)similarity_matrix.rows(); ++col)
		{
			random_walk_matrix(row,col) = similarity_matrix(row,col)/weight_vector(row);
			outFile << random_walk_matrix(row,col) << "\t";
		}
		outFile << "\n";
	}
	outFile.close();

	std::cout << " Done\nComputing Diffusion Map...";
	
	vnl_real_eigensystem Eyegun( random_walk_matrix );
	
	vnl_matrix<vcl_complex<double> > EVals = Eyegun.D;
	vnl_vector<double> temp_vec(similarity_matrix.rows());
	for(int row=0; row<EVals.rows(); ++row)
	{
		temp_vec(row) = vnl_real(EVals)(row,row);
	}

	vnl_matrix<double> temp_mat = Eyegun.Vreal;

	EigVals.set_size(data_matrix.rows());
	EigVecs.set_size(data_matrix.rows(), data_matrix.rows());

	for(int i=0; i<(int)temp_vec.size(); ++i)
	{
		EigVals(i) = temp_vec.max_value();
		int pos = temp_vec.arg_max();
		temp_vec(pos) = temp_vec.arg_min() - 1;
		EigVecs.set_column(i, temp_mat.get_column(pos));
	}

	std::cout << " Done\n";
}


vtkSmartPointer< vtkTable > DiffusionMap::GetKDiffusionNeighbors(std::vector<unsigned int> IDs, unsigned int k)
{
	std::map< unsigned int, std::vector<double> > eigenMap;
	for(int row=0; row<(int)EigVecs.rows(); ++row)
	{
		unsigned int id = (int)rowIdMap[row];		
		std::vector<double> c;
		for ( int i=0; i<10; ++i)
		{		
			c.push_back(EigVecs(row,i));	
		}		
		eigenMap[id] = c;
	}

	kNearestObjects<10>* KNObj = new kNearestObjects<10>(eigenMap);
	//KNObj->setFeatureTable(table);
	std::vector<std::vector< std::pair<unsigned int, double> > > kNeighborIDs;
	if(IDs.at(0) == 0)
		kNeighborIDs = KNObj->k_nearest_neighbors_All(k, 0, 0);
	else
		kNeighborIDs = KNObj->k_nearest_neighbors_IDs(IDs, k, 0);

	for(int i=0; i<kNeighborIDs.size(); ++i)
	{
		for(int j=0; j<kNeighborIDs[i].size(); ++j)
			std::cout << kNeighborIDs[i][j].first << "_";
		std::cout << "\n";
	}

	return KNObj->vectorsToGraphTable(kNeighborIDs);
}