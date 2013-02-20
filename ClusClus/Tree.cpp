#include "Tree.h"

BiTree::BiTree()
{
	this->reorderFlag = false;
}

BiTree::~BiTree()
{
}

void BiTree::setDataToTree(std::vector<std::vector<double > > & originalmatrix, std::vector<std::vector<double > > & matrix)
{
	if(originalmatrix.empty())
	{
		std::cout<<"No data to construct tree !......."<<std::endl;
		return;
	}
	this->matrix = matrix;
	this->originalmatrix = originalmatrix;
	this->num_rows = matrix.size();
	this->num_cols = matrix[0].size();
	this->num_rows_original = originalmatrix.size();
	this->num_cols_original = originalmatrix[0].size();
}

void BiTree::setFirstIterationClusterNumber(int firstK)
{ 
	this->firstK = firstK;
}

void BiTree::constructTree()
{
	this->InitialLevelId();
	int num_cluster = this->firstK;
	std::vector<std::vector<double > > datatobekmeans = this->matrix;
	int L = 1;

	while(num_cluster > 0)
	{
		std::vector<std::set<int > > clustersToPoints;
		std::vector<std::vector<double > > normalized;
		this->normalize(datatobekmeans, normalized);

		//cluster with normalized data 
		clustersToPoints = this->KmeansClustering(normalized,num_cluster);

		//cluster with unnormalized data
		//clustersToPoints = this->KmeansClustering(datatobekmeans,num_cluster);

		std::vector<std::vector<double > > meanmatrix;
		for(int i = 0; i < num_cluster; i++)
		{
			std::vector<std::vector<double > > submatrix;
			std::vector<int > ctp;
			std::set<int >::iterator it;
			for(it = clustersToPoints[i].begin(); it != clustersToPoints[i].end(); it++)
			{
				ctp.push_back(*it);
				submatrix.push_back(datatobekmeans[*it]);
			}
		
			if(submatrix.empty())
			{
				continue;
			}
			//function clustering reorder maybe
			if((this->reorderFlag == true) && (submatrix.size() > 2))
			{
				this->reOrder(submatrix, ctp);
			}

			std::vector<double > mean;
			mean = this->computerMeanSubmatrix(submatrix);
			meanmatrix.push_back(mean);
			this->levels.level_id.resize(L + 1);
			this->levels.level_id[L].push_back(ctp);
		}
		this->levels.numm_levels++;
		datatobekmeans = meanmatrix;

		//function to append meanmatrix to original matrix
		this->append(datatobekmeans);
		if(num_cluster == 1)
			num_cluster = 0;
		else
		{
			float last_num_cluster = (float)num_cluster;
			num_cluster = ceil(last_num_cluster/(float)3) ;
		}
		L++;
	}
}

void BiTree::constructSpecTree()
{
	this->InitialLevelId();
	int num_cluster = this->firstK;
	std::vector<std::vector<double > > datatobekmeans = this->matrix;
	int L = 1;

	while(num_cluster > 0)
	{
		std::vector<std::set<int > > clustersToPoints;
		std::vector<std::vector<double > > normalized;
		this->normalize(datatobekmeans, normalized);

		//cluster with normalized data 
		clustersToPoints = this->KmeansClustering(normalized,num_cluster);

		//cluster with unnormalized data
		//clustersToPoints = this->KmeansClustering(datatobekmeans,num_cluster);

		std::vector<std::vector<double > > meanmatrix;
		for(int i = 0; i < num_cluster; i++)
		{
			std::vector<std::vector<double > > submatrix;
			std::vector<int > ctp;
			std::set<int >::iterator it;
			for(it = clustersToPoints[i].begin(); it != clustersToPoints[i].end(); it++)
			{
				ctp.push_back(*it);
				submatrix.push_back(datatobekmeans[*it]);
			}
		
			if(submatrix.empty())
			{
				continue;
			}
			//function clustering reorder maybe
			if((this->reorderFlag == true) && (submatrix.size() > 2))
			{
				this->reOrder(submatrix, ctp);
			}

			std::vector<double > mean;
			mean = this->computerMeanSubmatrix(submatrix);
			meanmatrix.push_back(mean);
			this->levels.level_id.resize(L + 1);
			this->levels.level_id[L].push_back(ctp);
		}
		this->levels.numm_levels++;
		datatobekmeans = meanmatrix;

		//function to append meanmatrix to original matrix
		this->append(datatobekmeans);
		if(num_cluster == 1)
			num_cluster = 0;
		else
		{
			float last_num_cluster = (float)num_cluster;
			num_cluster = ceil(last_num_cluster/(float)3) ;
		}
		L++;
	}
}

void BiTree::InitialLevelId()
{
	this->levels.numm_levels = 0;
	this->levels.level_id.resize(1);
	for(int i = 0; i < this->num_rows; i ++)
	{
		std::vector<int > temp;
		temp.push_back(i);
		this->levels.level_id[0].push_back(temp);
	}
}

void BiTree::append(std::vector<std::vector<double > > & datatobekmeans)
{
	int p = datatobekmeans.size();
	for(int i = 0; i < p; i++)
	{
		std::vector<double > temp;
		temp.resize(this->num_cols_original);
		for(int j = 0; j < this->num_cols_original; j++)
		{
			temp[j] = datatobekmeans[i][j];
		}
		this->originalmatrix.push_back(temp);
	}
}
void BiTree::WriteFile(const char *filename1)
{
	FILE *fp1 = fopen(filename1,"w");
	for(int i=0; i<this->levels.numm_levels + 1; i++)
	{
		for(int j=0; j<this->levels.level_id[i].size(); j++)
		{
			for(int k=0; k<this->levels.level_id[i][j].size(); k++)
			{
				fprintf(fp1,"%d",levels.level_id[i][j][k]);
				fprintf(fp1,"%c",'*');
			}				
			fprintf(fp1,"\t");
		}
		fprintf(fp1,"\n");
	}
	fclose(fp1);
}

std::vector<std::set<int > > BiTree::KmeansClustering(std::vector<std::vector<double > > & datatobekmeans, int num_cluster)
{
	std::vector<std::set<int > > clustersToPoints;
	Kmeans* kmeans = new Kmeans();
	kmeans->setDataToKmeans(datatobekmeans);
	kmeans->setClusterNumber(num_cluster);
	kmeans->Clustering();
	clustersToPoints = kmeans->getClustersToPoints();
	delete kmeans;
	return clustersToPoints;
}

std::vector<double > BiTree::computerMeanSubmatrix(std::vector<std::vector<double > > & submatrix)
{
	std::vector<double > mean;
	mean.resize(submatrix[0].size());
	for(int j = 0; j < submatrix[0].size(); j++)
	{
		for(int i = 0; i < submatrix.size(); i++)
			mean[j]+= submatrix[i][j];
		mean[j] /= submatrix.size();
	}
	return mean;		
}

std::vector<std::vector<double > > BiTree::getPartitionTree()
{
	return this->originalmatrix;
}

Level_id BiTree::getLevels()
{
	return this->levels;
}

void BiTree::reOrder(std::vector<std::vector<double > > submatrix, std::vector<int > & ctp)
{
	double** data;
	int nrow = submatrix.size();
	int ncol = submatrix[0].size();
	data = new double* [nrow];
	for(int i = 0; i < nrow; i++)
	{
		data[i] = new double [ncol];
		for(int j = 0; j < ncol; j++)
			data[i][j] = submatrix[i][j];
	}

	clusclus *cc = new clusclus(data, nrow, ncol);
	cc->RunClusClus();
	std::cout<<"Finish clusclus now! "<<std::endl;

	std::vector<int > temp;
	for(int i = 0; i < nrow; i++)
		temp.push_back(ctp[cc->optimalleaforder[i]]);

	ctp = temp;

	for(int i = 0; i < nrow; i++)
	{
		delete data[i];
	}
	delete data;
	delete cc;
}

void BiTree::setReorderFlag(bool value)
{
	this->reorderFlag = true;
}

void BiTree::normalize(std::vector<std::vector<double > > datatobenormalized, std::vector<std::vector<double > > & normalized)
{
	int numr = datatobenormalized.size();
	int numc = datatobenormalized[0].size();
	normalized.resize (numr);

	for(int j = 0; j < numc; j++)
	{
		double mean = 0.0;
		for(int i = 0; i < numr; i++)
		{
			mean += datatobenormalized[i][j];
		}
		mean /= numr;

		double sum = 0.0;
		for(int i = 0; i < numr; i++)
		{
			sum += ( datatobenormalized[i][j] - mean ) * ( datatobenormalized[i][j] - mean );
		}
		double std = sqrt(sum / numr);

		for(int i = 0; i < numr; i++)
		{
			normalized[i].push_back ((datatobenormalized[i][j] - mean) / std);
		}
	}
}