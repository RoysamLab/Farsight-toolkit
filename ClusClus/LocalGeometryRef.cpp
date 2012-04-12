#include"LocalGeometryRef.h"

LocalGeometryRef::LocalGeometryRef()
{
}

LocalGeometryRef::~LocalGeometryRef()
{
}

void LocalGeometryRef::Initialize(vnl_matrix<double> data)
{
	this->num_row = data.rows();
	this->num_col = data.columns();
	this->data_matrix.set_size(this->num_row, this->num_col);
	for(int row = 0; row < this->num_row; row++)
	{
		for(int col = 0; col < this->num_col; col++)
		{
			this->data_matrix(row, col) = data(row, col);
		}
	}
}

void LocalGeometryRef::Initialize(vtkSmartPointer< vtkTable > table)
{
	this->num_row = table->GetNumberOfRows();
	this->num_col = table->GetNumberOfColumns();
	this->data_matrix.set_size(this->num_row, this->num_col);

	for(int row = 0; row < this->num_row; row++)
	{
		for(int col = 1; col < this->num_col; col++)
		{
			double var = table->GetValue(row, col).ToDouble();
			if( !boost::math::isnan(var))
			{
				this->data_matrix(row, col - 1) = var;
			}
			else
			{
				this->data_matrix(row, col - 1) = 0;
			}	
		}
	}
}

void LocalGeometryRef::Normalize()
{
	this->nor_data_matrix.set_size(this->num_row, this->num_col);
	this->nor_data_matrix.fill(0);
	for(int col = 0; col < this->num_col; col++)
	{
		double mean = 0.0;
		for(int row = 0; row < this->num_row; row++)
		{
			mean += this->data_matrix(row,col);
		}
		mean /= this->num_row;

		double sum = 0.0;
		for(int row = 0; row < this->num_row; row++)
		{
			sum += ( this->data_matrix(row,col) - mean ) * ( this->data_matrix(row,col) - mean );
		}
		double std = sqrt(sum / this->num_row);

		for(int row = 0; row < this->num_row; row++)
		{
			this->nor_data_matrix(row, col) = (this->data_matrix(row,col) - mean)/std;
		}
	}
}
void LocalGeometryRef::ComputeSimilarityMatrix(char dis, int k)
{
	this->Normalize();
	//initialize similarity matrix
	this->similarity_matrix.set_size(this->num_row, this->num_row);
	this->similarity_matrix.fill(0);

	//searching for the k nearest neighbours
	//std::vector<std::vector< std::pair<unsigned int, double> > > kNeighborIDs;
	//kNeighborIDs = this->KnnSearch(k);

	//starting compute similarity matrix
	for(int i = 0; i < this->num_row; i++)
	{
		for(int j = i + 1; j < this->num_row; j++)
		{
			vnl_vector<double> dt1 = this->nor_data_matrix.get_row(i);
			vnl_vector<double> dt2 = this->nor_data_matrix.get_row(j);
			double distance = ComputeSimilarity(dt1,dt2,dis);
			similarity_matrix(i,j) = distance;
			similarity_matrix(j,i) = distance;
		}
	}
}

double LocalGeometryRef::ComputeSimilarity(vnl_vector<double> dt1, vnl_vector<double> dt2, char dis)
{
	switch(dis)
	{
	case 'G':
		return exp( -((dt1-dt2).squared_magnitude()) / (2*pow((double)2/*sigma*/, 2)) );
		break;
	case 'E':
		break;
	default:
		break;
	}
}

std::vector<std::vector< std::pair<unsigned int, double> > > LocalGeometryRef::KnnSearch(int k)
{
	std::map< unsigned int, std::vector<double> > centroidMap;
	for(int row = 0; row < this->num_row; row++)
	{
		std::vector<double> c;
		for(int col = 0; col < this->num_col; col++)
		{
			c.push_back(this->nor_data_matrix(row, col));
		}
		centroidMap[row] = c;	
	}

	kNearestObjects<10>* KNObj = new kNearestObjects<10>(centroidMap);
	std::vector<std::vector< std::pair<unsigned int, double> > > kNeighborIDs;
	kNeighborIDs = KNObj->k_nearest_neighbors_All(7, 0, 0);
	return kNeighborIDs; 
}

void LocalGeometryRef::ComputeProbabilityMatrix()
{
	this->probability_matrix.set_size(this->num_row, this->num_row);
	vnl_matrix< double > D;
	D.set_size(this->num_row, this->num_row);
	D.fill(0);
	for(int i = 0; i < this->num_row; i++)
	{
		D(i, i) = this->similarity_matrix.get_row(i).sum(); 
	}

	for(int row=0; row <this->num_row; row++)
	{
		for(int col=0; col<this->num_row; col++)
		{
			this->probability_matrix(row,col) = this->similarity_matrix(row,col)/D(row, row);
		}
	}

	vnl_symmetric_eigensystem<double> sysD(D);
	this->D12 = sysD.square_root();
	this->Dn12 = sysD.inverse_square_root();	

	this->probability_matrix = D12 * this->probability_matrix;
	this->probability_matrix = this->probability_matrix * Dn12;

	//output files for debugging
	FILE *fp1 = fopen("D.txt","w");
	for(int i=0; i<this->num_row; i++)
	{
		for(int j=0; j<this->num_row; j++)
		{
			fprintf(fp1,"%f\t",D(i, j));
		}
		fprintf(fp1,"\n");
	}
	fclose(fp1);

	FILE *fp2 = fopen("D12.txt","w");
	for(int i=0; i<this->num_row; i++)
	{
		for(int j=0; j<this->num_row; j++)
		{
			fprintf(fp2,"%f\t",D12(i, j));
		}
		fprintf(fp2,"\n");
	}
	fclose(fp2);

	FILE *fp3 = fopen("similarity.txt","w");
	for(int i=0; i<this->num_row; i++)
	{
		for(int j=0; j<this->num_row; j++)
		{
			fprintf(fp3,"%f\t",this->similarity_matrix(i, j));
		}
		fprintf(fp3,"\n");
	}
	fclose(fp3);

	FILE *fp4 = fopen("probability.txt","w");
	for(int i=0; i<this->num_row; i++)
	{
		for(int j=0; j<this->num_row; j++)
		{
			fprintf(fp4,"%f\t",this->probability_matrix(i, j));
		}
		fprintf(fp4,"\n");
	}
	fclose(fp4);

	FILE *fp5 = fopen("Dn12.txt","w");
	for(int i=0; i<this->num_row; i++)
	{
		for(int j=0; j<this->num_row; j++)
		{
			fprintf(fp5,"%f\t",Dn12(i, j));
		}
		fprintf(fp5,"\n");
	}
	fclose(fp5);
}


void LocalGeometryRef::SVD(int num_eigen)
{
	vnl_symmetric_eigensystem<double> eigen(this->probability_matrix);
	vnl_matrix<double> eig_vec = eigen.V;
	vnl_diag_matrix<double> eig_val = eigen.D;
	
	vnl_vector<double> temp_vec(this->num_row);
	for(int row = 0; row<eig_val.rows(); row++)
	{
		temp_vec(row) = eig_val(row,row);
	}

	this->eigVals.set_size(num_eigen);
	this->eigVecs.set_size(this->num_row, num_eigen);
	for(int i=0; i<num_eigen; i++)
	{
		this->eigVals(i) = temp_vec.max_value();
		int pos = temp_vec.arg_max();
		temp_vec(pos) = temp_vec.min_value() - 1;
		this->eigVecs.set_column(i, eig_vec.get_column(pos) * this->eigVals(i) );
	}

	//output files for debugging
	FILE *fp1 = fopen("eigenvecters.txt","w");
	for(int i=0; i<this->num_row; i++)
	{
		for(int j=0; j<num_eigen; j++)
		{
			fprintf(fp1,"%f\t",this->eigVecs(i, j));
		}
		fprintf(fp1,"\n");
	}
	fclose(fp1);

	FILE *fp2 = fopen("eigenvalues.txt","w");
	for(int i=0; i<num_eigen; i++)
	{
		fprintf(fp2,"%f\t",this->eigVals(i));
		std::cout<<this->eigVals(i)<<std::endl;
		fprintf(fp2,"\n");
	}
	fclose(fp2);
}

vnl_matrix<double> LocalGeometryRef::GetSimilarityMatrix()
{
	return this->similarity_matrix;
}

vnl_matrix<double> LocalGeometryRef::GetProbabilityMatrix()
{
	return this->probability_matrix;
}

vnl_vector<double> LocalGeometryRef::GetEigenValues()
{
	return this->eigVals;
}

vnl_matrix<double> LocalGeometryRef::GetEigenVectors()
{
	return this->eigVecs;
}