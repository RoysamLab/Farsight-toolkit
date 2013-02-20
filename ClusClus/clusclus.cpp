#include"clusclus.h"

clusclus::clusclus(int linkmode)
{
	this->num_samples = 0;
	this->num_features = 0;	
	this->linkmode = linkmode;
	this->num_gaps = 5;
	this->gap = NULL;
	this->mergers = NULL;
	this->treedata = NULL;
	this->features = NULL;
	this->members = NULL;
	this->progress = NULL;
	this->optimalleaforder = NULL;
	this->sample_distances = NULL;
	this->transposefeatures = NULL;
	this->cluster_distances = NULL;
	this->num_cluster_samples = NULL;	
}

clusclus::clusclus(double** feature,int numsamples, int numfeatures, int linkmode)
{
	this->linkmode = linkmode;
	this->num_gaps = 5;
	this->num_features = numfeatures;
	this->num_samples = numsamples;
	this->gap = new double*[num_samples-1];
	this->progress = new int*[num_samples];
	this->members = new int*[num_samples];
	this->mergers = new double*[num_samples-1];
	this->treedata = new double*[num_samples-1];
	this->features = new double*[num_samples];
	this->num_cluster_samples = new int[num_samples];
	this->sample_distances = new double[num_samples*(num_samples+1)/2];
	this->cluster_distances = new double[num_samples*(num_samples+1)/2];
	this->transposefeatures = new double*[num_features];
	this->optimalleaforder = new int[num_samples];

	for(int i=0; i<num_samples; i++)
	{
		this->progress[i] = new int[num_samples];
		this->members[i] = new int[2];
		this->features[i] = new double[num_features + 2];
		for(int j = 0 ; j<num_features; j++)
			this->features[i][j] = feature[i][j];
	}
	for(int i=0; i<num_samples-1; i++)
	{
		this->mergers[i] = new double[5];
		this->gap[i] = new double[8];
		this->treedata[i] = new double[4];
	}
	for(int i=0; i<num_features-1; i++)
		this->transposefeatures[i] = new double[num_samples + 2];
}

clusclus::~clusclus()
{
	if(this->features)
	{
		for(int i=0; i<num_samples; i++)
		{
			delete this->progress[i];
			delete this->members[i];
			delete this->features[i];
		}
		for(int i=0; i<num_samples-1; i++)
		{
			delete this->gap[i];
			delete this->mergers[i];
			delete this->treedata[i];
		}
		for(int i=0; i<num_features-1; i++)
			delete this->transposefeatures[i];

		delete this->gap;
		delete this->progress;
		delete this->members;
		delete this->mergers;
		delete this->treedata;
		delete this->features;
		delete this->num_cluster_samples;
		delete this->sample_distances;
		delete this->cluster_distances;
		delete this->transposefeatures;
		delete this->optimalleaforder;
	}
	else if( this->treedata)   // the clusclus is intialized by tree data and used for getting the optimal order
	{
		for( int i = 0; i < num_samples-1; i++)
		{
			delete this->treedata[i];
		}
		delete this->treedata;
		this->treedata = NULL;
		delete this->optimalleaforder;
		this->optimalleaforder = NULL;
	}
}

void clusclus::Initialize(double** feature,int numsamples, int numfeatures)
{
	if(this->features)
	{
		for(int i=0; i<num_samples; i++)
		{
			delete this->progress[i];
			delete this->members[i];
			delete this->features[i];
		}
		for(int i=0; i<num_samples-1; i++)
		{
			delete this->gap[i];
			delete this->mergers[i];
			delete this->treedata[i];
		}
		for(int i=0; i<num_features-1; i++)
			delete this->transposefeatures[i];

		delete this->gap;
		delete this->progress;
		delete this->members;
		delete this->mergers;
		delete this->treedata;
		delete this->features;
		delete this->num_cluster_samples;
		delete this->sample_distances;
		delete this->cluster_distances;
		delete this->transposefeatures;
		delete this->optimalleaforder;
	}

	this->num_features = numfeatures;
	this->num_samples = numsamples;

	this->gap = new double*[num_samples-1];
	this->progress = new int*[num_samples];
	this->members = new int*[num_samples];
	this->mergers = new double*[num_samples-1];
	this->treedata = new double*[num_samples-1];
	this->features = new double*[num_samples];
	this->num_cluster_samples = new int[num_samples];
	this->sample_distances = new double[num_samples*(num_samples+1)/2];
	this->cluster_distances = new double[num_samples*(num_samples+1)/2];
	this->transposefeatures = new double*[num_features];
	this->optimalleaforder = new int[num_samples];

	for(int i=0; i<num_samples; i++)
	{
		this->progress[i] = new int[num_samples];
		this->members[i] = new int[2];
		this->features[i] = new double[num_features + 2];
		for(int j = 0 ; j<num_features; j++)
			this->features[i][j] = feature[i][j];
	}

	for(int i=0; i<num_samples-1; i++)
	{
		this->mergers[i] = new double[5];
		this->gap[i] = new double[8];
		this->treedata[i] = new double[4];
	}
	for(int i=0; i<num_features-1; i++)
		this->transposefeatures[i] = new double[num_samples + 2];
}

void clusclus::Initialize( double** treedata, int numsamples)
{
	this->num_samples = numsamples;
	if( this->treedata)
	{
		for( int i = 0; i < numsamples - 1; i++)
		{
			delete this->treedata[i];
		}
		delete this->treedata;
		this->treedata = NULL;
	}

	this->treedata = new double*[numsamples - 1];
	this->optimalleaforder = new int[num_samples];

	for(int i = 0; i < numsamples - 1; i++)
	{
		this->treedata[i] = new double[4];
		for( int j = 0; j < 4; j++)
		{
			this->treedata[i][j] = treedata[i][j];
		}
	}

	//FILE *fp = fopen("ClusInitialize.txt","w");
	//for(int i=0; i<num_samples - 1; i++)
	//{
	//	for(int j=0; j<4; j++)
	//		fprintf(fp,"%f\t",this->treedata[i][j]);
	//	fprintf(fp,"\n");
	//}
	//fclose(fp);
}

void clusclus::Clustering()
{
	int num_currcluster = num_samples;
	int pivot1,pivot2;

	ComputeSampleDistances();
	for (int i = 0; i < num_samples; i++)
	{
		features[i][num_features]=i;
		features[i][num_features+1]=i;
		num_cluster_samples[i]=1;
		for(int j = 0; j< i+1; j++)
		{
			cluster_distances[i*(1+i)/2+j]=sample_distances[i*(1+i)/2+j];
		}
	}
	for(int i=0; i<num_samples-1; i++)
	{
		if( num_currcluster % 100 == 0)
		{
			std::cout<<"left: "<<num_currcluster<< std::endl;
		}
		mergers[i][0] = i;
		mergers[i][4]=MergeClusters(num_currcluster, &pivot1, &pivot2);
		mergers[i][1]=pivot1;
		mergers[i][2]=pivot2;
		mergers[i][3]=UpdateClusterDistances(num_currcluster,pivot1,pivot2);

		num_currcluster--;
	}
}

void clusclus::ComputeSampleDistances()
{
	#pragma omp parallel for
	for(int i=0; i<num_samples; i++)
	{
		for(int j=0; j<i+1; j++)
		{
			sample_distances[i*(i+1)/2+j]= ComputeDistance(features[i],features[j]);
		}
	}
}

double clusclus::MergeClusters(int num_currcluster, int *pivot1, int *pivot2)
{     
	int		 pivot;
	double   minn = 1.0e30;

	for (int k = 0; k < num_currcluster; k++)
	{
		for (int kk = 0; kk < k; kk++)
		{ 
			if (cluster_distances[k*(k+1)/2+kk] < minn)
			{
				minn = cluster_distances[k*(k+1)/2+kk];
				*pivot2 = k;
				*pivot1 = kk;
			}
		}
	}

	for (int i = 0; i < num_samples; i++) 
	{
		pivot = (int) features[i][num_features+1];
		if (pivot == *pivot2)
		{
			features[i][num_features+1] = *pivot1;
		}
		else if (pivot > *pivot2)
		{
			features[i][num_features+1] = pivot-1;
		}
	}

	num_cluster_samples[*pivot1] += num_cluster_samples[*pivot2];
	for(int i=*pivot2; i<num_currcluster - 1; i++)
	{
		num_cluster_samples[i] = num_cluster_samples[i+1];
	}
	return minn;
}

/*void clusclus::NormalizeFeatures()
{
	for(int i=0; i<this->num_features; i++)
	{
		double avg = 0.0;
		for(int j=0; j<this->num_samples; j++)
		{
			avg += this->features[j][i];
		}
		avg /= this->num_samples;

		double sum = 0.0;
		for(int j=0; j<this->num_samples; j++)
		{
			sum += (this->features[j][i]-avg)*(this->features[j][i]-avg);
		}

		double stdd = sqrt(sum / this->num_samples);
		for(int j=0; j<this->num_samples; j++)
		{
			this->features[j][i] = (this->features[j][i]-avg)/stdd;
		}
	}
}*/

double clusclus::ComputeDistance(double* vector1, double* vector2)
{
	double   sum, distance;

	sum = 0.0;
	for (int i = 0; i < num_features; i++) 
	{
		sum += (vector1[i]-vector2[i])*(vector1[i]-vector2[i]);
	}
	distance = sqrt(sum);
	return distance;
}

double clusclus::UpdateClusterDistances(int num_currcluster, int pivot1, int pivot2)
{
	int counter =0;
	int num_sample_p1 = num_cluster_samples[pivot1];
	double dispersion;
	int * clusp1_ID;
	clusp1_ID = new int[num_sample_p1];

	for(int i=pivot2; i<num_currcluster-1; i++)
	{
		for(int j=0; j<pivot2; j++)
		{
			cluster_distances[i*(i+1)/2+j]=cluster_distances[(i+2)*(i+1)/2+j];
		}
		for(int j=pivot2; j<i+1; j++)
		{
			cluster_distances[i*(i+1)/2+j]=cluster_distances[(i+2)*(i+1)/2+j+1];
		}
	}
	for(int i=0; i<num_samples; i++)
	{
		if(features[i][num_features+1] == pivot1)
		{
			clusp1_ID[counter++]= i;
		}
	}
	for(int k=0; k<num_currcluster-1; k++)
	{
		int counter=0;		
		int num_sample_k = num_cluster_samples[k];
		double tempdistance=0.0,finaldistance=0.0,temp;
		double minn = 1.0e30, maxx = -1.0e30;
		int* cluspk_ID;
		cluspk_ID = new int[num_sample_k];
		
		for(int j=0; j<num_samples; j++)
		{
			if(features[j][num_features+1] == k)
			{
				cluspk_ID[counter++]= j;
			}
		}
		for(int m=0; m<num_sample_p1; m++)
		{
			int clusp1_m = clusp1_ID[m];
			for(int n=0; n<num_sample_k; n++)
			{
				int cluspk_n = cluspk_ID[n];
				if(clusp1_m>cluspk_n)
				{
					temp = sample_distances[clusp1_m*(clusp1_m+1)/2+cluspk_n];
					tempdistance += temp;
				}
				else 
				{
					temp = sample_distances[cluspk_n*(cluspk_n+1)/2+clusp1_m];
					tempdistance += temp;
				}
				if(temp < minn)
					minn = temp;
				if(temp > maxx)
					maxx = temp;
			}
		}
		if(this->linkmode == 2)
		{
			finaldistance = sqrt(tempdistance*tempdistance/num_sample_p1/num_sample_k);
		}
		else if(this->linkmode == 1)
			finaldistance = minn;
		else if(this->linkmode == 3)
			finaldistance = maxx;

		if(k < pivot1) cluster_distances[pivot1*(pivot1+1)/2+k] = finaldistance;
		else if(k > pivot1) cluster_distances[k*(k+1)/2+pivot1] = finaldistance;
		else cluster_distances[k*(k+1)/2+k] = tempdistance*tempdistance/num_sample_p1/num_sample_k;

		delete cluspk_ID;
	}
	delete clusp1_ID;

	dispersion = 0.0;
	for(int k=0; k<num_currcluster-1; k++)
	{
		if(num_cluster_samples[k]!= 0)
		{
			dispersion += cluster_distances[k*(k+1)/2+k]*num_cluster_samples[k];
		}
	}
	return sqrt(dispersion/num_samples);
}

void clusclus::MergersToProgress()
{
	int pivot1, pivot2;
	for(int i=0; i<num_samples; i++)
	{
		progress[i][0] = i;
	}
	for(int k=0; k<num_samples-1; k++)
	{
		pivot1 = mergers[k][1];
		pivot2 = mergers[k][2];
		for(int i=0; i<num_samples; i++)
		{
			progress[i][k+1] = progress[i][k];
		}
		for(int i=0; i<num_samples; i++)
		{
			if(progress[i][k+1] == pivot2) 
			{
				progress[i][k+1] = pivot1; 
			}
			if(progress[i][k+1] > pivot2)
			{
				progress[i][k+1] -= 1;
			}
		}
	}
}

void clusclus::GetMembers(int num_cluster)
{
	int counter = 0;
	for(int k=0; k<num_cluster; k++)
	{
		for(int i=0; i<num_samples; i++) 
		{
			if (progress[i][num_samples-num_cluster]==k)
			{
				members[counter][0]=k;
				members[counter++][1]=i;
			}
		}
	}
}
void clusclus::WriteClusteringOutputToFile(const char *filename1, const char *filename2,const char *filename3, 
										   const char *filename4, const char* filename5, const char *filename6, const char *filename7)
{
	//FILE *fp1 = fopen(filename1,"w");
	//for(int i=0; i<num_samples-1; i++)
	//{
	//	for(int j=0; j<5; j++)
	//		fprintf(fp1,"%14.7e\t",mergers[i][j]);
	//	fprintf(fp1,"\n");
	//}
	//fclose(fp1);

	//FILE *fp2 = fopen(filename2,"w");
	//for(int i=0; i<num_samples; i++)
	//{
	//	for(int j=0; j<num_features+2; j++)
	//		fprintf(fp2,"%f\t ",features[i][j]);
	//	fprintf(fp2,"\n");
	//}
	//fclose(fp2);

	//FILE *fp3 = fopen(filename3,"w");
	//for(int i=0; i<num_samples; i++)
	//{
	//	for(int j=0; j<num_samples; j++)
	//		fprintf(fp3,"%d\t ",progress[i][j]);
	//	fprintf(fp3,"\n");
	//}
	//fclose(fp3);

	//FILE *fp4 = fopen(filename4,"w");
	//for(int i=0; i<num_samples; i++)
	//{
	//	for(int j=0; j<2; j++)
	//		fprintf(fp4,"%d\t",members[i][j]);
	//	fprintf(fp4,"\n");
	//}
	//fclose(fp4);

	//FILE *fp5 = fopen(filename5,"w");
	//for(int i=0; i<num_gaps; i++)
	//{
	//	for(int j=0; j<8; j++)
	//		fprintf(fp5,"%14.7e\t",gap[i][j]);
	//	fprintf(fp5,"\n");
	//}
	//fclose(fp5);

	FILE *fp6 = fopen(filename6,"w");
	for(int i=0; i<num_samples -1; i++)
	{
		for(int j=0; j<4; j++)
			fprintf(fp6,"%f\t",treedata[i][j]);
		fprintf(fp6,"\n");
	}
	fclose(fp6);

	FILE *fp7 = fopen(filename7,"w");
	for(int i=0; i<num_samples; i++)
	{
		fprintf(fp7,"%d\t",optimalleaforder[i]);
		fprintf(fp7,"\n");
	}
	fclose(fp7);

/*	FILE *fp6 = fopen(filename6,"w");
	for(int i=0; i<num_samples; i++)
	{
		for(int j=0; j<2; j++)
			fprintf(fp6,"%d\t",metrix[i][j]);
		fprintf(fp6,"\n");
	}
	fclose(fp6);*/
}
void clusclus::WriteClusteringOutputToFile(const char *filename)
{
	FILE *fp = fopen(filename,"w");
	for(int i=0; i<num_gaps; i++)
	{
		for(int j=0; j<8; j++)
			fprintf(fp,"%14.7e\t",gap[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);
}

int  clusclus::ComputeGapStatistics()
{
	double ratio, maxm ;
	int cluster_number ;

	ratio = 2.0*mergers[num_samples-2][3]/mergers[num_samples-2][4]/3.0;
	for(int k=0; k<num_samples-1; k++)
	{
		gap[k][0] = k+1;
		gap[k][1] = mergers[num_samples-2-k][3];		
	}
		gap[0][7] = ratio*mergers[num_samples-2][4];
	for(int k=0; k<num_samples-2; k++)
	{
		gap[k+1][7] = ratio*mergers[num_samples-2-k][4];
		gap[k][2] = gap[k+1][1] - gap[k][1];
		gap[k][5] = gap[k][7]-gap[k+1][7];
	}

	maxm = -1000;
	for(int k=0; k < num_gaps; k++)
	{
		if (gap[k][5]>maxm)
		{
			maxm = gap[k][5];
			cluster_number = k+1;
		}
	}
	return cluster_number;
}

void clusclus::ReadFile(const char *filename)
{
	FILE *fp = fopen(filename,"r");
	if(fp == NULL)
	{
		fprintf(stderr,"can't open input file %s\n",filename);
		exit(1);
	}

	int n=0;
	while(1)
	{
		int c = fgetc(fp);
		switch(c)
		{
			case '\n':				
				(num_samples)++;	
				n++;
				if(num_features == 0)num_features = n;
				break;	
			case '\t':
				n++;
				break;
			case EOF:
				goto out;
			default:
				;
		}
	}
out:
	rewind(fp);

	features = new double*[num_samples];
	for(int i=0; i<num_samples; i++)
	{
		double temp;
		features[i] = new double[num_features+2];

		for(int j=0; j<num_features; j++)
		{
			fscanf(fp, "%lf", &temp);
			features[i][j] = temp;
		}
	}

	this->mergers = new double*[num_samples-1];
	this->progress = new int*[num_samples];
	this->num_cluster_samples = new int[num_samples];
	this->members = new int*[num_samples];
	this->sample_distances = new double[num_samples*(num_samples+1)/2];
	this->cluster_distances = new double[num_samples*(num_samples+1)/2];
	this->gap = new double*[num_samples-1];
	this->transposefeatures = new double*[num_features];
	this->treedata = new double*[num_samples-1];
	this->optimalleaforder = new int[num_samples];

	for(int i=0; i<num_samples; i++)
	{
		this->progress[i] = new int[num_samples];
		this->members[i] = new int[2];
	}
	for(int i=0; i<num_samples-1; i++)
	{
		this->mergers[i] = new double[5];
		this->gap[i] = new double[8];
		this->treedata[i] = new double[4];
	}
	for(int i=0; i<num_features-1; i++)
		this->transposefeatures[i] = new double[num_samples + 2];

	return ;
}

void clusclus::Transpose()
{

	transposefeatures = new double*[num_features];
	for( int i= 0; i<num_features; i++)
	{
		transposefeatures[i] = new double[num_samples+2];
		for( int j= 0; j<num_samples; j++)
		{
			transposefeatures[i][j] = features[j][i];
		}
	}
	return ;
}

void clusclus::PrepareTreeData()
{
	int      i, ii, j, jj, jjlast, indj, k, kk, kklast, indk;
	int 	 found, found1, found2;
	double*  distances = new double[num_samples - 1];

	for (i = 0; i < num_samples - 1; i++) 
		distances[i] = mergers[i][4];               //get distances

	for (i = 0; i < num_samples-1; i++) 
	{
		treedata[i][0] = -1;					    //placeholder for cluster number, all real cluster numbers are non-neg
		treedata[i][1] = -1;						//placeholder for cluster number, all real cluster numbers are non-neg
	}
	//main loop: work through columns of progress
	for (i = 0; i < (num_samples-1); i++) 
	{
		j = (int)mergers[i][1];                       //low cluster number
		k = (int)mergers[i][2];                       //high cluster number

		//Find parent for cluster number j
		for (ii = 0; ii < num_samples; ii++)
			if ((int)progress[ii][i] == j)
				indj = ii;

		//find row in matlab_chris that includes indj
		jj = -1;
		for (ii = 0; ii < num_samples-1; ii++) 
		{
			if (treedata[ii][0] == indj) jj = ii;
			if (treedata[ii][1] == indj) jj = ii;
		}
		//find ultimate parent --------------------
		if (jj == -1) jjlast = indj;                  //if not in tree, add new cluster number
		else 
		{ 								                 //if currently in tree, walk through tree
			found = -1;									//found1 = -1; found2 = -1;
			while (found == -1) 
			{
				jj = num_samples + jj;
				jjlast = jj;                              //save jj in jjlast incase jj is not currently in tree

				//look in tree, see if jj is in either column
				found1 = FindRowIndex(0, jj); //find jj in matlab_chris col 0
				found2 = FindRowIndex(1, jj); //find jj in matlab_chris col 1
				if (found1 != -1) jj = found1;                        //if found, save index
				if (found2 != -1) jj = found2;                        //if found, save index
				if ((found1 == -1) && (found2 == -1)) found = -2;     //if not found in either, exit loop
			}
		}

		//Find parent for cluster number k --------
		for (ii = 0; ii < num_samples; ii++) 
			if ((int)progress[ii][i] == k) indk = ii;
		//find row in matlab_chris that includes indk
		kk = -1;
		for (ii = 0; ii < num_samples-1; ii++) 
		{
			if (treedata[ii][0] == indk) kk = ii;
			if (treedata[ii][1] == indk) kk = ii;
		}
		//find ultimate parent ----------------------
		if (kk == -1) kklast = indk; 				                //if not in tree, add new cluster number
		else 
		{ 									                         //if currently in tree, walk through tree
			found = -1; //found1 = -1; found2 = -1;
			while (found == -1) 
			{
				kk = num_samples + kk;
				kklast = kk;                            //save kk in kklast incase kk is not currently in tree
				//look in tree, see if kk is in either column
				found1 = FindRowIndex(0, kk);//find kk in matlab_chris col 0
				found2 = FindRowIndex(1, kk);//find kk in matlab_chris col 1
				if (found1 != -1) kk = found1;                       //if found, save index
				if (found2 != -1) kk = found2;                       //if found, save index
				if ((found1 == -1) && (found2 == -1)) found = -2;    //if not found in either, exit loop
			}
		}

		//add cluster numbers to matlab_chris ------
		treedata[i][0] = jjlast < kklast ? jjlast : kklast;
		treedata[i][1] = jjlast > kklast ? jjlast : kklast;
	}
	for (ii = 0; ii < num_samples-1; ii++) 
	{
		treedata[ii][2] = distances[ii];                    //add distance column
		treedata[ii][3] = num_samples+ii;                      //add cluster number
	}
	
	delete distances;
	return;
}

int  clusclus::FindRowIndex(int col, int c) 
{
	int     i, ii;

	ii = -1;
	 for (i = 0; i < num_samples-1; i++) 
	 {
		 if (treedata[i][col] == c) ii = i;
  }
  return ii;
}

void clusclus::GetOptimalLeafOrderD() 
{
	int  i, k;
	int* Tnums;
	int* pickedup;

	Tnums = new int[num_samples - 1];
	pickedup = new int[num_samples];

	//Initialize variables
	k = 2*(num_samples-1);                                //maximum cluster number
	for (i = 0; i < num_samples-1; i++) Tnums[i] = num_samples+i;
	for (i = 0; i < num_samples; i++) pickedup[i] = -1;
	for (i = 0; i < num_samples; i++) optimalleaforder[i] = -1;

	//Get kids
	GetKids(k, Tnums, pickedup, optimalleaforder);

	delete Tnums;
	delete pickedup;

	return;
}

void clusclus::GetKids(int k, int* Tnums, int* pickedup, int* kids) 
{
	int   i = 0, j = 0, k1, k2;

	if (k < num_samples)                        //if node number is a leaf
	{     
		while (j == 0) 
		{                        //find first open position in pickedup
			if (kids[i] == -1) 
			{
				j = 1;
				kids[i] = k;
				pickedup[i] = k;   //may not need this
			}
			i++;
		}
	}
	else                     //if node number is not a leaf
	{                  
		for (j = 0; j < num_samples-1; j++) 
		{
			if (Tnums[j] == k) i = j;
		}                     //find index of k in Tnums
		k1 = (int)treedata[i][0];
		k2 = (int)treedata[i][1];
		GetKids(k1, Tnums, pickedup, kids); //pick up kids from node k1
		GetKids(k2, Tnums, kids, kids);     //pick up kids from node k2
	}
	return;
}

void clusclus::RunClusClus()
{
	this->Clustering();
	this->MergersToProgress();
	this->PrepareTreeData();

	//double max1 = 0;
	//int maxInd = 0;
	//double max2 = 0;
	//double min = 1e10;
	//for (int i = 0; i < num_samples-1; i++) 
	//{
	//	if( treedata[i][2] > max1) 
	//	{
	//		max2 = max1;
	//		max1 = treedata[i][2];
	//		maxInd = i;
	//	}

	//	if( treedata[i][2] < min)
	//	{
	//		min = treedata[i][2];
	//	}
	//}

	//treedata[maxInd][2] = max2+10;
	//if(max2 + 10 > min)
	//{
	//	for(int i = 0; i < num_samples-1; i++)
	//	{
	//		treedata[i][2] = (treedata[i][2] - min) / ( max2 + 10 - min);
	//	}
	//}
	this->GetOptimalLeafOrderD();
}

void clusclus::GetTreeStructure(std::vector< ClusterTree> &treeVec)
{
	treeVec.clear();
	for(int i = 0; i < num_samples - 1; i++)
	{
		ClusterTree tree((int)treedata[i][0], (int)treedata[i][1], treedata[i][2], (int)treedata[i][3]);
		treeVec.push_back(tree);
	}
}