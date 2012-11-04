/*=========================================================================
Copyright 2009 Rensselaer Polytechnic Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 
=========================================================================*/
#if defined(_MSC_VER)
#pragma warning(disable : 4996)
#endif
#include "fuzzy_clustering.h"
#include <math.h>
#include <stdlib.h> /* exit */
#include <cstdio>


// CTOR1
FuzzyClustering::FuzzyClustering(int numC)
{
	num_Clusters = numC;
	num_Iterations = 100;
	precision = 0.001;
	num_Features = 0;
	num_Samples = 0;
	weighting_exp = 2.0;
	trainingSetSize = 0;
	Features = NULL;
	Centers = NULL;
	new_Centers = NULL;
	Memberships = NULL;
	sum_Memberships = NULL;
	Distances = NULL;
	sum_Distances = NULL;
	Singularities = NULL;
	clustering_Output = NULL;
	cluster_Distributions = NULL;
	trainingSet = NULL;
}
// CTOR2
FuzzyClustering::FuzzyClustering(int numC, int numIter, int prec)
{
	num_Clusters = numC;
	num_Iterations = numIter;
	precision = prec;
	num_Features = 0;
	num_Samples = 0;
	weighting_exp = 2.0;
	trainingSetSize = 0;
	Features = NULL;
	Centers = NULL;
	new_Centers = NULL;
	Memberships = NULL;
	sum_Memberships = NULL;
	Distances = NULL;
	sum_Distances = NULL;
	Singularities = NULL;
	clustering_Output = NULL;
	cluster_Distributions = NULL;
	trainingSet = NULL;
}
// DTOR
FuzzyClustering::~FuzzyClustering()
{
	for(int i=0;i<num_Samples;i++)
	{
		delete Features[i];
		delete Memberships[i];
		delete Distances[i];
		delete Singularities[i];		
	}
	for(int i=0; i<num_Clusters; i++)
	{
		delete Centers[i];
		delete new_Centers[i];
	}
	for(int i=0; i<trainingSetSize; i++)
		delete trainingSet[i];

	delete Features;
	delete Memberships;
	delete Centers;
	delete new_Centers;
	delete Distances;
	delete sum_Distances;
	delete sum_Memberships;
	delete Singularities;
	delete clustering_Output;
	delete cluster_Distributions;
	delete trainingSet;
	delete labels;
}
//by Yousef (04-02-2009):
//this is copied from the read_problem function available in svm-train.c
//I did that because we will be using the same format used by libsvm
void FuzzyClustering::ReadProblem(const char *filename)
{
	int n, i, j;
	FILE *fp = fopen(filename,"r");
	
	if(fp == NULL)
	{
		fprintf(stderr,"can't open input file %s\n",filename);
		exit(1);
	}

	n=0;
	while(1)
	{
		int c = fgetc(fp);
		switch(c)
		{
			case '\n':				
				++num_Samples;				
				if(num_Features == 0)
					num_Features = n;
			case ':':
				++n;
				break;
			case EOF:
				goto out;
			default:
				;
		}
	}
out:
	rewind(fp);
	
	Features = new double*[num_Samples];
	Distances = new double*[num_Samples];
	Memberships = new double*[num_Samples];
	sum_Memberships = new double[num_Clusters];
	sum_Distances = new double[num_Samples];
	Singularities = new int*[num_Samples];
	//added by Yousef on 10/13/2009
	//To hold the input labels (first column of the file)
	labels = new int[num_Samples];
			
	
	
	for(i=0;i<num_Samples;i++)
	{
		j=0;
		Features[i] = new double[num_Features];
		Distances[i] = new double[num_Clusters];
		Memberships[i] = new double[num_Clusters];
		Singularities[i] = new int[num_Clusters];

		double label;		
		if( fscanf(fp,"%lf",&label) == EOF )
      {
      std::cerr << "fscanf encountered end of file!" << std::endl;
      }
		labels[i] = (int) label;
		
		int dd;
		while(1)
		{
			int c;
			do {
				c = getc(fp);
				if(c=='\n') goto out2;
			} while(isspace(c));
			ungetc(c,fp);
			if (fscanf(fp,"%d:%lf",&dd,&Features[i][j]) < 2)
			{
				fprintf(stderr,"Wrong input format at line %d\n", i+1);
				exit(1);
			}
			++j;
		}			
	out2:
		//do nothing
		dd=0;
	}	

	fclose(fp);

	//Normalize the features
	NormalizeFeatures();
}

//This function initializes the centers of the clusters
//For now, I will just use the first num_Clusters points in the list
void FuzzyClustering::InitializeCenters()
{		
	//get the first num_Clusters points
	Centers = new double*[num_Clusters];
	new_Centers = new double*[num_Clusters];
	for(int i=0; i<num_Clusters; i++)
	{
		Centers[i] = new double[num_Features];
		new_Centers[i] = new double[num_Features];
		for(int j=0; j<num_Features; j++)
		{
			new_Centers[i][j] = Centers[i][j] = Features[i][j];
		}
	}
}
//This function compute the Euclidean distances between the points and the current centers of the clusters
void FuzzyClustering::ComputeDistances()
{
	for(int i=0; i<num_Samples; i++)
	{
		sum_Distances[i] = 0;
		for(int j=0; j<num_Clusters; j++)
		{
			Distances[i][j] = 0;
			for(int k=0; k<num_Features; k++)
				Distances[i][j]+= ((Features[i][k]-Centers[j][k])*(Features[i][k]-Centers[j][k]));
			Distances[i][j] = sqrt(Distances[i][j]);
			if(Distances[i][j]!=0)
			{
				Distances[i][j] = pow(Distances[i][j],(-2/(weighting_exp-1)));
				Singularities[i][j] = 0;
				sum_Distances[i]+=Distances[i][j];
			}
			else
				Singularities[i][j] = 1; //this means that the point is a center point
		}
	}
}
//This function computes the mebership degrees of each point to each cluster
void FuzzyClustering::ComputeMembershipValues()
{
	for(int j=0; j<num_Clusters; j++)
		sum_Memberships[j] = 0;
	for(int i=0; i<num_Samples; i++)
	{		
		double sum_memb = 0.0;
		for(int j=0; j<num_Clusters; j++)
		{
			if(Singularities[i][j] == 0)
				Memberships[i][j] = Distances[i][j]/sum_Distances[i];
			else
				Memberships[i][j] = 1;
			Memberships[i][j] = pow(Memberships[i][j],weighting_exp);
			sum_Memberships[j] += Memberships[i][j];
			sum_memb+=Memberships[i][j];
		}

		//try this: Normalize membership values		
		for(int j=0; j<num_Clusters; j++)
			Memberships[i][j]/=sum_memb;
	}	
}
//This function updates the cluster centers
void FuzzyClustering::FindNewCenters()
{
	double S;
	for(int j=0; j<num_Clusters; j++)
	{
		for(int k=0; k<num_Features; k++)
		{
			S = 0;
			for(int i=0; i<num_Samples; i++)
			{
				S+= Memberships[i][j]*Features[i][k];
			}
			new_Centers[j][k] = S/sum_Memberships[j];
		}
	}
}
//This function computes the error in between the new and the old cluster centers
//This is done by computing the maximum error between the corresponding old and new centers
/*double FuzzyClustering::ComputeCurrentError()
{
	double Err, S;
	
	for(int i=0; i<num_Clusters; i++)
	{
		S = 0;
		for(int j=0; j<num_Features; j++)
		{
			S += ((new_Centers[i][j]-Centers[i][j])*(new_Centers[i][j]-Centers[i][j]));		
			//also copy new to old
			Centers[i][j] = new_Centers[i][j];
		}
		S = sqrt(S);
		if(i==0)
			Err = S;
		else		
			Err = (Err>S)? Err : S;		
	}
	return Err;
}*/

double FuzzyClustering::ComputeCurrentError()
{
  double Err = 0.0;

  for(int i=0; i<num_Clusters; i++)
    {
    double S = 0.0;
    for(int j=0; j<num_Features; j++)
      {
      double feature_distance = new_Centers[i][j]-Centers[i][j];
      S += (feature_distance*feature_distance);		
      //also copy new to old
      Centers[i][j] = new_Centers[i][j];
      }
    S = sqrt(S);
    if (S > Err)
      Err = S;		
    }
  return Err;
}


//This function runs one step of the fuzzy c-means clustering
void FuzzyClustering::RunOneStep()
{
	//compute distances between the points and the current cluster centers
	ComputeDistances();
	
	//calculate new membership values
	ComputeMembershipValues();

	//find new centers
	FindNewCenters();
}

//This function finds the class (cluster) label for each point based on the maximum fuzzy membership value
void FuzzyClustering::AssignClusters()
{
	int ind = 0;
	double memb = 0.0;
	clustering_Output = new int[num_Samples];
	cluster_Distributions = new int[num_Clusters];
	for(int i=0; i<num_Clusters; i++)
		cluster_Distributions[i] = 0;

	for(int i=0; i<num_Samples; i++)
	{
		memb =  Memberships[i][0];
		ind = 0;
		for(int j=1; j<num_Clusters; j++)
		{
			if(Memberships[i][j] > memb)
			{
				memb = Memberships[i][j];
				ind = j;
			}
		}
		//assign label
		clustering_Output[i] = ind + 1;
		//increase the number of samples in that cluster
		cluster_Distributions[ind]++;
	}
}

//This is the main clustering function 
void FuzzyClustering::Run()
{
	double Err;
	//initialize centers
	InitializeCenters();
	//here is the main clustering loop
	for(int i=0; i<num_Iterations; i++)
	{		
		//run one clustering step
		RunOneStep();
		
		//compute errors (max difference between old and new centers)
		Err = ComputeCurrentError();		

		//check the termination condition
		if(Err <= precision)
		{
			std::cout<< "Number of Iterations is: "<<i+1<<"\nCurrent Error is "<< Err <<std::endl;
			break;
		}
	}

	//assign clusters
	AssignClusters();
}
//This function extracts a training set using the top C percent
double** FuzzyClustering::ExtractTrainingSet(double C)
{ 
	//Allocate memory and initialize the training set
	trainingSetSize = 0;
	for(int i=0; i<num_Clusters; i++)
	{
		int cc = (int)(cluster_Distributions[i]*C/100);		
		if(cluster_Distributions[i]>0 && cc==0)
			cc=1;
		trainingSetSize+=cc;
	}
	trainingSet = new double*[trainingSetSize];
	for(int i=0; i<trainingSetSize; i++)
		trainingSet[i] = new double[num_Features+1];

	int training_index = -1;
	int err = 0;
	for(int i=0; i<num_Clusters; i++)
	{
		//get the membership values of the samples assigned to the ith cluster
		double* samples = new double[cluster_Distributions[i]];
		int* ids = new int[cluster_Distributions[i]];
		int ind = -1;
		for(int j=0; j<num_Samples; j++)
		{			
			int cl = clustering_Output[j] - 1;
			if(cl==i)
			{
				ind++;
				samples[ind] = Memberships[j][i];
				ids[ind] = j;
			}
		}
		if(ind<0)
			continue;

		//sort based on the membership values
		SortMembershipValues(samples, ids, cluster_Distributions[i]);
		
		//get the top C percent
		int numTrnSmpls = (int) (cluster_Distributions[i]*C/100);
		if(cluster_Distributions[i]>0 && numTrnSmpls<1)
			numTrnSmpls = 1;
		for(int k1=0; k1<numTrnSmpls; k1++)
		{
			training_index++;			
			trainingSet[training_index][0] = i+1;
			for(int k2=0; k2<num_Features; k2++)
				trainingSet[training_index][k2+1] = Features[ids[k1]][k2];
			if(labels[ids[k1]] != i+1)
				err++;
		}

		//deallocate memory
		delete samples;
		delete ids;
	}

	std::cout<<"The training set has "<<training_index<<" elements with "<<err<<" errors"<<std::endl;

	return trainingSet; 
}

//this function will sort the membership values and will change the ids respectively
void FuzzyClustering::SortMembershipValues(double *vals, int* ids, int l)
{
	for(int i=0; i<l-1; i++)
	{
		double v = vals[i];
		int id = ids[i];
		for(int j=i+1; j<l; j++)
		{
			if(vals[j]>v)
			{
				vals[i] = vals[j];
				vals[j] = v;
				v = vals[i];
				ids[i] = ids[j];
				ids[j] = id;
				id = ids[i];
			}
		}
	}
}

//The next function will write the clustering output to a file
//each point in the output is represented by a row
//the first element in the row is the cluster assignemnt, while the others
//represent the cluster membership degrees
void FuzzyClustering::WriteClusteringOutputToFile(const char *filename)
{
	FILE *fp = fopen(filename,"w");
	for(int i=0; i<num_Samples; i++)
	{
		//write the cluster number for which the ith point is assigned
		fprintf(fp,"%d ",clustering_Output[i]);
		//write the membership values
		for(int j=0; j<num_Clusters; j++)
			fprintf(fp,"%f ",Memberships[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);
}
//The next function will write the training set to a file using the libsvm format
void FuzzyClustering::WriteTrainingSetToFile(const char *filename)
{
	FILE *fp = fopen(filename,"w");

	for(int i=0; i<trainingSetSize; i++)
	{
		//write the label
		fprintf(fp,"%d ",(int) trainingSet[i][0]);
		//write the feature vector
		for(int j=1; j<=num_Features; j++)		
			fprintf(fp,"%d:%f ", j,trainingSet[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);
}

//Added by Yousef on 10/13/2009
//This function is used for validation purposes
//It assumes that the first column of the file contains the labels (classes)
//It relables the output in order to minimize the error and it retures the percentage of accurately classified objects
int FuzzyClustering::Validate()
{
	//we need to find the best arrangement of the labels (classes) such that the number of errors will be minimized
	int factorial = 1;
	for(int i=num_Clusters; i>=1; i--)
		factorial *= i;
	int **combs = new int*[factorial];
	for(int i=0; i<factorial; i++)
	{
		combs[i] = new int[num_Clusters];
		for(int j=0; j<num_Clusters; j++)
			combs[i][j] = 0;
	}

	int *s;
	int lenth = num_Clusters;		
	s= new int[lenth];
	for(int i=0; i<lenth; i++)
		s[i] = i+1;	

	int ll = 0;


	rotrec(0,lenth,lenth, s, combs, &ll); 

	int* tmpLbl = new int[num_Samples];
	for(int i=0; i<num_Samples; i++)
		tmpLbl[i] = clustering_Output[i];

	int num_err = num_Samples;
	int best_comb = 0;
	for(int i=0; i<factorial; i++)
	{
		rearrangeLabels(tmpLbl, combs[i]);		
		int err = 0;
		for(int j=0; j<num_Samples; j++)
		{
			if(tmpLbl[j] != labels[j])
			{
				err++;
			}
		}
		if(err<num_err)
		{
			num_err = err;
			best_comb = i;
		}
	}
	rearrangeLabels(clustering_Output, combs[best_comb]);	
	double* tmpMembVals = new double[num_Clusters];
	for(int i=0; i<num_Clusters; i++)
	{
		cluster_Distributions[i]=0;		
	}
	for(int i=0; i<num_Samples; i++)
	{		
		cluster_Distributions[clustering_Output[i]-1]++;	
		for(int j=0; j<num_Clusters; j++)
			tmpMembVals[j] = Memberships[i][combs[best_comb][j]-1];
		for(int j=0; j<num_Clusters; j++)
			Memberships[i][j] = tmpMembVals[j];
	}
	float err_perc = ((float) num_err / (float) num_Samples)*100;
	std::cout<<"Clustering accuracy = "<<100-err_perc<<"%"<<std::endl;

	for(int i=0; i<factorial; i++)
	{				
		delete combs[i];
	}
	delete [] combs;
	delete [] s;
	delete [] tmpLbl;	
	delete [] tmpMembVals;
	
	return num_err;
}

//Added by Yousef on 10/13/2009
//Relabel the classes
void FuzzyClustering::rearrangeLabels(int* tmpLables, int* arrangement)
{
	for(int i=0; i<num_Samples; i++)
	{
		int cl = clustering_Output[i];	
		for(int j=0; j<num_Clusters; j++)
		{
			if(cl == arrangement[j])
			{
				tmpLables[i] = j+1;
				continue;
			}
		}
	}
}
//Added by Yousef on 10/14/2009
void FuzzyClustering::NormalizeFeatures()
{
	//Normalize by subtracting the mean and then dividing by the standard deviation	
	for(int i=0; i<num_Features; i++)
	{
		//1. Compute the mean (average) of each feature
		double avg = 0.0;
		for(int j=0; j<num_Samples; j++)
		{
			avg +=Features[j][i];
		}	
		avg/=num_Samples;

		//1. Compute the standard deviation of each feature
		double stdd = 0.0;
		for(int j=0; j<num_Samples; j++)
		{
			stdd+= sqrt((Features[j][i]-avg)*(Features[j][i]-avg));
		}
		stdd/=num_Samples;

		//3. Normalize
		for(int j=0; j<num_Samples; j++)
		{
			Features[j][i] = (Features[j][i]-avg)/stdd;
		}
	}
}

void FuzzyClustering::rotate(int *c,int start, int len)
{     
	int x;     
	int temp; 	
	temp=c[start]; 
	x=len-start-1;     
	while(x>0)     
	{
		c[start]=c[start+1];
		x--; 
		start++;     
	}     
	c[start]=temp;
}
void FuzzyClustering::rotrec(int start,int lenth,int len, int *c, int** combs, int* ll)
{     
	int i;     
	i=0;     
	if(lenth==0) 
	{
		//printf("%d\t",ll[0]);		
		for(int k=0; k<len; k++)
			//printf("%d",c[k]);
			combs[ll[0]][k] = c[k];
		ll[0]++;
		//printf("\n");
	}
	else     
	{        
		while(i < lenth)         
		{           
			rotrec(start+1,lenth-1,len, c, combs, ll);           
			rotate(c,start,len);           
			i++;         
		}     
	}
}
