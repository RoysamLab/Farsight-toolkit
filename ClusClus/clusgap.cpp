#include"clusgap.h"

clusgap::clusgap(clusclus* object, int numtrials, int numgaps)
{
	this->features = object->features;
	this->num_samples = object->num_samples;
	this->num_features = object->num_features;
	this->num_trials = numtrials;
	this->num_gaps = numgaps;
	this->cc = object; 

	this->syntheticfeatures = new double*[num_samples];	
	this->dispersionmatrix = new double*[num_samples-1];

	for(int i = 0; i<num_samples; i++)
		this->syntheticfeatures[i] = new double[num_features+2];	
	for(int i = 0; i<num_samples-1; i++)
		this->dispersionmatrix[i] = new double[num_trials];
}

clusgap::clusgap(clusclus* object)
{
	this->features = object->features;
	this->num_samples = object->num_samples;
	this->num_features = object->num_features;
	this->num_trials = 20;
	this->num_gaps = 20;
	this->cc = object; 

	this->syntheticfeatures = new double*[num_samples];	
	this->dispersionmatrix = new double*[num_samples-1];

	for(int i = 0; i<num_samples; i++)
		this->syntheticfeatures[i] = new double[num_features+2];	
	for(int i = 0; i<num_samples-1; i++)
		this->dispersionmatrix[i] = new double[num_trials];
}
clusgap::~clusgap()
{
	for(int i = 0; i<num_samples; i++)
		delete this->syntheticfeatures[i];	
	for(int i = 0; i<num_samples-1; i++)
		delete this->dispersionmatrix[i];

	delete this->syntheticfeatures;
	delete this->dispersionmatrix;
}

int clusgap::ComputeGap()
{

	int seed = 5;
	int num_cluster;
	double avg, std, offset, maxm;

	for(int i = 0; i< num_trials; i++)
	{
		GenerateSyntheticData(seed);
		clusclus *ccsynthetic = new clusclus(syntheticfeatures, num_samples, num_features);
		ccsynthetic->RunClusClus();
		for(int k = 0; k<num_samples-1; k++)
		{
			dispersionmatrix[k][i] = ccsynthetic->mergers[num_samples-2-k][3];
		}			
	}
	for(int i = 0; i < num_gaps; i++)
	{
		Stand_Devv(&avg, &std, dispersionmatrix[i]);
		cc->gap[i][3] = avg;
	}
	offset = cc->gap[0][3] - cc->gap[0][1];
	for (int i = 0; i < num_gaps; i++) cc->gap[i][3] = cc->gap[i][3] - offset;   
    for (int i = 0; i < num_gaps; i++) cc->gap[i][4] = cc->gap[i][3] - cc->gap[i][1];   
    for (int i = 0; i < num_gaps-1; i++) cc->gap[i][5] = cc->gap[i+1][4] - cc->gap[i][4]; 
    for (int i = 0; i < num_gaps-2; i++) cc->gap[i][6] = cc->gap[i+1][5] - cc->gap[i][5]; 
	
	maxm = 0.0;
	for(int k = 0; k<num_gaps - 2; k++)
	{
		if(fabs(cc->gap[k][6])>maxm)
		{
			if(cc->gap[k][6]<0)
			{
				maxm = fabs(cc->gap[k][6]);
				num_cluster = (int)cc->gap[k][0]+1;
			}
		}
	}

	//cout<<"......................"<< num_cluster <<endl;
	return num_cluster;
}

void clusgap::GenerateSyntheticData(int seed)
{
	int      i, j;
	double   minn, maxx;

	srand(seed);
	for (j = 0; j < num_features; j++)
	{
		minn = 999;
		maxx = -999;
		for (i = 0; i < num_samples; i++)
		{
			if (features[i][j] < minn) minn = features[i][j];
			if (features[i][j] > maxx) maxx = features[i][j];
		}
		for (i = 0; i < num_samples; i++) syntheticfeatures[i][j] = minn + (rand()%32000)/32000.0*(maxx-minn);
	}
	for (i = 0; i < num_samples; i++) 
	{
		syntheticfeatures[i][num_features] = features[i][num_features];
	}
}

void clusgap::Stand_Devv(double *avg, double *std, double* my_vector) 
{
	double   tempd, sum = 0.0;
	for (int i = 0; i < num_trials; i++) 
		sum += my_vector[i];

	*avg = sum/num_trials;

	sum = 0.0;
	for (int i = 0; i < num_trials; i++) 
	{
		tempd = my_vector[i]-*avg;
		sum += tempd*tempd;
	}
	*std = sqrt(sum/(num_trials-1));
}