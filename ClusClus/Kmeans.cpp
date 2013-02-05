#include "Kmeans.h"

//Constructor
Kmeans::Kmeans()
{	
	this->num_clusters = 1;
	this->num_iterations = 0;
	this->distance_mode = 1;
	
}

//Deconstructor
Kmeans::~Kmeans()
{
}

//Set Data to the Kmeans object
void Kmeans::setDataToKmeans(std::vector<std::vector<double > > & points)
{
	if(points.empty())
	{
		std::cout<<"No data to clustering !......."<<std::endl;
		return;
	}
	this->points = points;
	this->num_rows = points.size();
	this->num_cols = points[0].size();
}

//set Data to the Kmeans object in vnl form
void Kmeans::setDataToKmeans(vnl_matrix<double> points)
{
	if(points.empty())
	{
		std::cout<<"No data to clustering !......."<<std::endl;
		return;
	}
	
	this->num_rows = points.rows();
	this->num_cols = points.columns();
	this->points.resize(points.rows());
	for(int row = 0; row < this->num_rows; row++)
	{
		this->points[row].resize(this->num_cols);
		for (int col = 0; col < this->num_cols; col++)
		{
			this->points[row][col] = points(row,col);
		}
	}
}

//Set distance mode(Euclidean or Cosine)
void Kmeans::setDistasnce(char distancestring)
{
	switch(distancestring)
	{
	case 'E':
		this->distance_mode = 1;
		break;
	case 'C':
		this->distance_mode = 2;
		break;
	case 'G':
		this->distance_mode = 3;
		break;
	default:
		break;
	}
}

//Set the number of clusters would like 
void Kmeans::setClusterNumber(int clusternumber)
{
	this->num_clusters = clusternumber;
}

//Clustering iteratiom
void Kmeans::Clustering()
{
	//whether it si improved or not in this iteration
	bool move = true ;
	if(this->num_clusters < 1)
	{
		std::cout<<" cluster number is less than 1, no need to clustering"<<std::endl;
		return;
	}

	int num_iteration = 0;

	this->initialClusters();
	this->initialCentroids();

	int counter = 0;
	while(move && counter < 50)
	{
		move = false;

		this->updateCentroids();

		//repositioning each point to new cluster
		this->updateClusters(move);

		//update each cluster with new point ids
		if(move)
		{
			this->updatePointIdInClusters();
		}

		counter++;
	}
}

//Allocate each point to a cluster by % operation
void Kmeans::initialClusters()
{
	int pointId;
	int clusterId;
	this->clustersToPoints.resize(this->num_clusters);
	for(pointId = 0; pointId < this->num_rows; pointId++)
	{
		clusterId = pointId % this->num_clusters;
		this->pointsToClusters.push_back(clusterId);
		this->clustersToPoints[clusterId].insert(pointId);
	}
}

//Resize the centroids 
void Kmeans::initialCentroids()
{
	this->centroids.resize(this->num_clusters);
	for(int clusterId = 0; clusterId < this->num_clusters; clusterId++)
		this->centroids[clusterId].resize(this->num_cols);

}

void Kmeans::updateCentroids()
{
	for(int clusterId = 0; clusterId < this->num_clusters; clusterId++)
	{
		this->updateCentroid(clustersToPoints[clusterId], centroids[clusterId]);
	}
}

//Reallocate each point to a new cluster
void Kmeans::updateClusters(bool & move)
{
	for(int pointId = 0; pointId < this->num_rows; pointId++)
	{
		int point_to_cluster = this->pointsToClusters[pointId];
		double min_distance = computeDistance(points[pointId], centroids[pointsToClusters[pointId]]);
		for(int k = 0; k < this->num_clusters; k++)
		{	
			double temp_distance = computeDistance(points[pointId], centroids[k]);
			if (temp_distance < min_distance)
			{
				min_distance = temp_distance;
				point_to_cluster = k;
				move = true;
			}
		}
		this->pointsToClusters[pointId] = point_to_cluster;
	}
}
void Kmeans::updatePointIdInClusters()
{
	for(int clusterId = 0; clusterId < this->num_clusters; clusterId++)
	{
		this->clustersToPoints[clusterId].clear();
	}
	for(int pointId = 0; pointId < this->num_rows; pointId++)
	{
		int point_to_cluster = this->pointsToClusters[pointId];
		this->clustersToPoints[point_to_cluster].insert(pointId);
	}
}

//Recompute each centroid
void Kmeans::updateCentroid(std::set<int > & clusterToPoints, std::vector<double > & point)
{
	std::set<int >::iterator it;
	int num_point_in_cluster = clusterToPoints.size(); 

	for(int j = 0; j < this->num_cols; j++)
	{
		point[j] = 0;
		for(it = clusterToPoints.begin(); it != clusterToPoints.end(); it++)
		{
			point[j] += points[*it][j]; 
		}
		point[j] /= num_point_in_cluster;		
	}
}


double Kmeans::computeDistance(std::vector<double > & point1,std::vector<double > & point2)
{
	if(distance_mode == 1)
	{
		double sum = 0.0;
		for(int i = 0; i < this->num_cols; i++)
		{
			sum += (point1[i] - point2[i]) * (point1[i] - point2[i]);
		}

		return sqrt(sum);
	}

	if(distance_mode == 3)
	{
		double sum = 0.0;
		for(int i = 0; i < this->num_cols; i++)
		{
			sum += (point1[i] - point2[i]) * (point1[i] - point2[i]);
		}
		long double exponent = sum;
		long double result = exp(- exponent / (2 * 20));
		return result;
	}
}


void Kmeans::WriteFile(const char *filename1, const char *filename2)
{
	FILE *fp1 = fopen(filename1,"w");
	for(int i=0; i<this->num_clusters; i++)
	{
		std::set<int >::iterator it;
		for(it = this->clustersToPoints[i].begin(); it!= this->clustersToPoints[i].end(); it++)
			fprintf(fp1,"%d\t",*it);
		fprintf(fp1,"\n");
	}
	fclose(fp1);

	FILE *fp2 = fopen(filename2,"w");
	for(int i=0; i<this->num_rows; i++)
	{
		fprintf(fp2,"%d\n",pointsToClusters[i]);
	}
	fclose(fp2);
}

std::vector<int > Kmeans::getPointsToClusters()
{
	return this->pointsToClusters;
}

std::vector<std::set<int > > Kmeans::getClustersToPoints()
{
	return this->clustersToPoints;
}

std::vector<std::vector<double > > Kmeans::getClusterCenters()
{
	return this->centroids;
}
std::vector<std::vector<int > > Kmeans::getClustersToPointsvector()
{
	this->clustersToPointsvector.clear();
	for(int i = 0; i < this->num_clusters; i++)
	{
		std::vector<int > tempvector;
		std::set<int >::iterator it = this->clustersToPoints[i].begin();
		for(int j = 0; j < this->clustersToPoints[i].size(); j++)
		{
			tempvector.push_back(*it++);
		}
		this->clustersToPointsvector.push_back(tempvector);
	}
	return this->clustersToPointsvector;
}
