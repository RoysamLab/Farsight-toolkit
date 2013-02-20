#include "ActiveValidation.h"

//Constructor
ActiveValidation::ActiveValidation()
{	
	this->data = NULL;
	this->lable = NULL;
	this->table = NULL;
	this->N = 0;
	this->numfeat = 0;
	this->numbin = 0;
	this->phat = 0;
	this->delta = 0;
	this->varphat = 0;
	this->numiteration = 0;
	this->stratas.clear();
	this->numsampled = 0;
}

//Deconstructor
ActiveValidation::~ActiveValidation()
{
}

//Initialize an instance
void ActiveValidation::Initializing(vnl_matrix<double> data, vnl_vector<int> lable, 
									int K, double delta, unsigned long seed)
{
	//check input data
	if(data.empty() || data.rows() < 1 || data.rows() != lable.size())
	{ 
		std::cout<<"invalid input data!!!"<<std::endl;
		return;
	}

	//clear existing variables
	if(!this->data.empty())
	{
		this->data = NULL;
		this->lable = NULL;
		this->table = NULL;
		this->N = 0;
		this->numfeat = 0;
		this->numbin = 0;
		this->phat = 0;
		this->varphat = 0;
		this->delta = 0;
		this->numiteration = 0;
		this->stratas.clear();
		this->numsampled = 0;
	}

	//reload and setup variables
	this->seed = seed;
	this->N = data.rows();
	this->numbin = K;
	this->delta = delta;
	this->numfeat = data.columns();
	this->data.set_size(this->N, this->numfeat);
	for(int row = 0; row < data.rows(); row++)
	{
		for(int col = 0; data.columns(); col++)
		{
			this->data(row, col) = data(row, col);
		}
	}
	this->lable.set_size(lable.size());
	for(int i = 0; i < lable.size(); i++)
	{
		this->lable(i) = lable(i);
	}
}

void ActiveValidation::Initializing(std::vector<std::vector<double > > datavector, 
								  std::vector<int > lablevector, int K, double delta, unsigned long seed)
{
	//check input data
	if(datavector.empty() || (datavector.size() < 1) || datavector.size() != lablevector.size())
	{ 
		std::cout<<"invalid input data!!!"<<std::endl;
		return;
	}

	//clear existing variables
	if(!this->data.empty())
	{
		this->data = NULL;
		this->lable = NULL;
		this->table = NULL;
		this->N = 0;
		this->numfeat = 0;
		this->numbin = 0;
		this->phat = 0;
		this->varphat = 0;
		this->delta = 0;
		this->numiteration = 0;
		this->stratas.clear();
		this->numsampled = 0;
	}

	//reload and setup variables
	this->seed = seed;
	this->N = datavector.size();
	this->numbin = K;
	this->delta = delta;
	this->numfeat = datavector[0].size();
	this->data.set_size(this->N, this->numfeat);
	for(int row = 0; row < this->N; row++)
	{
		for(int col = 0; col < this->numfeat; col++)
		{
			this->data(row, col) = datavector[row][col];
		}
	}
	this->lable.set_size(lablevector.size());
	for(int i = 0; i < lablevector.size(); i++)
	{
		this->lable(i) = lablevector[i];
	}
}
void ActiveValidation::Initializing(vtkSmartPointer<vtkTable > table, int K, double delta)
{
	//check input data
	if(table->GetNumberOfRows() < 1 )
	{ 
		std::cout<<"invalid input data!!!"<<std::endl;
		return;
	}

	//clear existing variables
	if(!this->data.empty())
	{
		this->data = NULL;
		this->lable = NULL;
		this->table = NULL;
		this->N = 0;
		this->numfeat = 0;
		this->numbin = 0;
		this->phat = 0;
		this->varphat = 0;
		this->delta = 0;
		this->numiteration = 0;
		this->stratas.clear();
		this->numsampled = 0;
	}

	//reload and setup variables
	this->current = 0;
	this->convergepub = false;
	this->t = 0;
	this->seed = 0;
	this->N = table->GetNumberOfRows();
	this->numbin = K;
	this->delta = delta;
	this->numfeat = table->GetNumberOfColumns() - 4;
	this->data.set_size(this->N, this->numfeat);
	for(int row = 0; row < this->N; row++)
	{
		for(int col = 0; col < this->numfeat; col++)
		{
			this->data(row, col) = table->GetValue(row, col + 1).ToDouble();
		}
	}
	this->table = table;
	this->lable.set_size(this->N);
}
//Get homogenous groups
void ActiveValidation::Stratifing()
{
	//Stratifing data set into numbin clusters using kmeans
	std::vector<std::vector<int > > clustersToPointsvector;
	Kmeans* kmeans = new Kmeans();
	kmeans->setDataToKmeans(this->data);
	kmeans->setClusterNumber(this->numbin);
	kmeans->Clustering();
	clustersToPointsvector = kmeans->getClustersToPointsvector();
	delete kmeans;

	//set up structure element for each strata
	for(int i = 0; i < this->numbin; i++)
	{
		if(clustersToPointsvector[i].size() == 0)
		{}
		else if( (i < numbin - 1) && (clustersToPointsvector[i + 1].size() == 0) )
		{
			std::vector<int >  clustersToPointsvector1;
			std::vector<int >  clustersToPointsvector2;
			int num1 = clustersToPointsvector[i].size() / 2;
			int num2 = clustersToPointsvector[i].size()  - num1;
			int k = 0;
			for(; k < num1; k++)
				clustersToPointsvector1.push_back(clustersToPointsvector[i][k]);
			for(; k < clustersToPointsvector[i].size(); k++)
				clustersToPointsvector2.push_back(clustersToPointsvector[i][k]);

			strata temp1;
			temp1.Nk = num1;
			temp1.idsk = clustersToPointsvector1;	
			vnl_matrix<double> tempdata1;
			tempdata1.set_size(temp1.Nk, this->numfeat);
			vnl_vector<int> templable1; 
			templable1.set_size(temp1.Nk);
			std::vector<int > tempunselectedidsk1;
			for(int j = 0; j < temp1.Nk; j++)
			{
				vnl_vector<double> temprow1 = this->data.get_row(clustersToPointsvector1[j]);
				tempdata1.set_row(j,temprow1);
				templable1[j] = this->lable[clustersToPointsvector1[j]];
				tempunselectedidsk1.push_back(j);
			}
			temp1.strdata = tempdata1;
			temp1.strlable = templable1;
			temp1.unselectedidsk = tempunselectedidsk1;
			temp1.numleft = temp1.Nk;
			temp1.hk = 0;
			temp1.nk = 0;
			this->stratas.push_back(temp1);

			strata temp2;
			temp2.Nk = num2;
			temp2.idsk = clustersToPointsvector2;	
			vnl_matrix<double> tempdata2;
			tempdata2.set_size(temp2.Nk, this->numfeat);
			vnl_vector<int> templable2; 
			templable2.set_size(temp2.Nk);
			std::vector<int > tempunselectedidsk2;
			for(int j = 0; j < temp2.Nk; j++)
			{
				vnl_vector<double> temprow2 = this->data.get_row(clustersToPointsvector1[j]);
				tempdata2.set_row(j,temprow2);
				templable2[j] = this->lable[clustersToPointsvector1[j]];
				tempunselectedidsk2.push_back(j);
			}
			temp2.strdata = tempdata2;
			temp2.strlable = templable2;
			temp2.unselectedidsk = tempunselectedidsk2;
			temp2.numleft = temp2.Nk;
			temp2.hk = 0;
			temp2.nk = 0;
			this->stratas.push_back(temp2);
		}
		else
		{
			strata temp;
			temp.Nk = clustersToPointsvector[i].size();
			temp.idsk = clustersToPointsvector[i];	

			vnl_matrix<double> tempdata;
			tempdata.set_size(temp.Nk, this->numfeat);
			vnl_vector<int> templable; 
			templable.set_size(temp.Nk);
			std::vector<int > tempunselectedidsk;
			for(int j = 0; j < temp.Nk; j++)
			{
				vnl_vector<double> temprow = this->data.get_row(clustersToPointsvector[i][j]);
				tempdata.set_row(j,temprow);
				templable[j] = this->lable[clustersToPointsvector[i][j]];
				tempunselectedidsk.push_back(j);
			}
			temp.strdata = tempdata;
			temp.strlable = templable;
			temp.unselectedidsk = tempunselectedidsk;
			temp.numleft = temp.Nk;
			temp.hk = 0;
			temp.nk = 0;
			this->stratas.push_back(temp);
		}
	}
}

//Get groups acording to claasification result
void ActiveValidation::StratifingwithoutLabel()
{
	this->stratas.resize(this->numbin);
	for(int i = 0; i< this->table->GetNumberOfRows(); i++)
	{	
		vtkVariant cn = table->GetValueByName(i, "prediction_active" );
		//std::cout<<table->GetColumnByName("prediction_active")->GetSize()<<std::endl;
		int strataid = cn.ToInt();
		this->stratas[strataid].idsk.push_back(i);
	}

	for(int i = 0; i < this->numbin; i++)
	{
		this->stratas[i].Nk = this->stratas[i].idsk.size();

		vnl_matrix<double> tempdata;
		tempdata.set_size(this->stratas[i].Nk, this->numfeat);
		std::vector<int > tempunselectedidsk;
		for(int j = 0; j < this->stratas[i].Nk; j++)
		{
			for(int  k = 0 ; k < this->numfeat; k++)
			{
				tempdata(j, k ) = this->table->GetValue(this->stratas[i].idsk[j], k + 1).ToDouble();
			}

			tempunselectedidsk.push_back(j);
		}
		this->stratas[i].strdata = tempdata;
		this->stratas[i].unselectedidsk = tempunselectedidsk;
		this->stratas[i].numleft = this->stratas[i].Nk;
		this->stratas[i].hk = 0;
		this->stratas[i].nk = 0;
		this->stratas[i].strlable.set_size(this->stratas[i].Nk);
	}
}

//Random sampling m variables from n entries
void ActiveValidation::RandomSampling(int m, int n, int kth_strata)
{

	std::vector<int > tempselecting;

	//psudo-random number generator using boost generator 
	//boost::mt19937 gen;                                    
    //boost::uniform_int<> dist(1,n);
    //boost::variate_generator<boost::mt19937&, boost::uniform_int<>> multi(gen,dist);
	//for(int i = 0; i < m; i++) 
	//{
	//	int tempint = multi() - 1;
	//	for(int k = 0; k < tempselecting.size(); k++)
	//	{
	//		if(tempint == tempselecting[k])
	//		{
	//			tempint = multi() - 1;
	//		}
	//	}
	//	tempselecting.push_back(tempint);
	//}

	//real random number generator using time as seed
	srand( (unsigned)time( NULL ) + this->seed );
	for(int i = 0; i < m; i++) 
	{
		int tempint = rand()%n;
		for(int k = 0; k < tempselecting.size(); k++)
		{
			if(tempint == tempselecting[k])
			{
				tempint = rand()%n;
			}
		}
		tempselecting.push_back(tempint);
	}
	this->stratas[kth_strata].selectingidsk = tempselecting;
	this->stratas[kth_strata].nk += m;
}

std::vector<int > ActiveValidation::RandomSamplingwithoutLabel(int m, int n, int kth_strata)
{
	std::vector<int > tempselecting;
	std::vector<int > queries;
	//real random number generator using time as seed
	srand( (unsigned)time( NULL ) + this->seed );
	for(int i = 0; i < m; i++) 
	{
		int tempint = rand()%n;
		for(int k = 0; k < tempselecting.size(); k++)
		{
			if(tempint == tempselecting[k])
			{
				tempint = rand()%n;
			}
		}
		tempselecting.push_back(tempint);
	}
	this->stratas[kth_strata].selectingidsk = tempselecting;
	this->stratas[kth_strata].nk += m;

	for(int i = 0; i < m; i++)
	{
		int id = tempselecting[i];
		queries.push_back(this->stratas[kth_strata].idsk[this->stratas[kth_strata].unselectedidsk [id] ]);
	}

	return queries;
}


//Split strata with variance larger than threshold
void ActiveValidation::StrataSplit(int kth_strata)
{
	//Split kth strata into 2 stratas using kmeans
	std::vector<std::vector<int > > clustersToPointsvector;
	Kmeans* kmeans = new Kmeans();
	kmeans->setDataToKmeans(this->stratas[kth_strata].strdata);
	kmeans->setClusterNumber(2);
	kmeans->Clustering();
	clustersToPointsvector = kmeans->getClustersToPointsvector();
	delete kmeans;

	//update structure element for splited strata
	int tempNk = clustersToPointsvector[0].size();
	vnl_matrix<double> tempmatrix(tempNk, this->numfeat);
	std::vector<int > tempidsk;
	std::vector<int > tempselectedidsk;
	vnl_vector<int> templable;
	templable.set_size(tempNk);

	//update idsk for cluster 1
	for(int i = 0; i < tempNk; i++)
	{
		int tempid = clustersToPointsvector[0][i];
		vnl_vector<double> temprow = this->stratas[kth_strata].strdata.get_row(tempid);
		tempmatrix.set_row(i,temprow);
		tempidsk.push_back(this->stratas[kth_strata].idsk[i]);
		templable[i] = this->stratas[kth_strata].strlable[i];
	}
	//update selected idsk for cluster 1
	for(int i = 0; i < this->stratas[kth_strata].selectedidsk.size(); i++)
	{
		int tempid = this->stratas[kth_strata].selectedidsk[i];
		for(int j = 0; j < clustersToPointsvector[0].size(); j++)
		{
			if(tempid == clustersToPointsvector[0][j])
			{
				tempselectedidsk.push_back(j);
				break;
			}
		}
	}

	this->stratas[kth_strata].Nk = tempNk;
	this->stratas[kth_strata].idsk = tempidsk;
	this->stratas[kth_strata].strdata = tempmatrix;	
	this->stratas[kth_strata].strlable = templable;
	this->stratas[kth_strata].selectedidsk = tempselectedidsk;
	this->stratas[kth_strata].nk = tempselectedidsk.size();

	//update structure element for cluster 2
	strata tempstrata;
	int tempNk2 = clustersToPointsvector[1].size();
	vnl_matrix<double> tempmatrix2(tempNk2, this->numfeat);
	std::vector<int > tempidsk2;
	std::vector<int > tempselectedidsk2;
	vnl_vector<int> templable2;
	templable2.set_size(tempNk2);

	//update idsk for cluster 2
	for(int i = 0; i < tempNk2; i++)
	{
		int tempid = clustersToPointsvector[1][i];
		vnl_vector<double> temprow = this->stratas[kth_strata].strdata.get_row(tempid);
		tempmatrix2.set_row(i,temprow);
		tempidsk2.push_back(this->stratas[kth_strata].idsk[i]);
		templable2[i] = this->stratas[kth_strata].strlable[i];	
	}
	//update selected idsk for cluster 2
	for(int i = 0; i < this->stratas[kth_strata].selectedidsk.size(); i++)
	{
		int tempid = this->stratas[kth_strata].selectedidsk[i];
		for(int j = 0; j < clustersToPointsvector[1].size(); j++)
		{
			if(tempid == clustersToPointsvector[1][j])
			{
				tempselectedidsk2.push_back(j);
				break;
			}
		}
	}

	tempstrata.Nk = tempNk2;
	tempstrata.idsk = tempidsk2;
	tempstrata.strdata = tempmatrix2;	
	tempstrata.selectedidsk = tempselectedidsk2;
	tempstrata.nk = tempselectedidsk2.size();
	tempstrata.strlable = templable2;

	this->stratas.push_back(tempstrata);
}

//Draw elements from entire dataset
void ActiveValidation::Sampling()
{
	bool converge = false;
	int t = 0;                                    //avoid spurious

	//Split data into numbin groups
	this->Stratifing();

	while(t < 2)
	{
		//Split check for each bin
		//for(int i = 0; i < this->numbin; i++)
		//{	
		//	bool check = this->SplitCheck(i);
		//	if(check)
		//	{
		//		this->StrataSplit(i);
		//	}
		//}

		//sampling from each bin and calculate parameter for each strata
		this->CalculatenumToBeDrawn();
		for(int i = 0; i < this->numbin; i++)
		{
			this->RandomSampling(this->stratas[i].numtobeDrawn, this->stratas[i].numleft, i);		
			this->Calculatehk(i);
			this->Calculatephatak(i);
			this->CalculateVarphatk(i);
			this->CalculateSigmak(i);
			this->UpdateBins(i);
		}

		//calculate result and check converge
		this->Calculatephata();
		this->CalculateVarphat();
		converge = this->ConvergeCheck();

		//avoid spurious
		if(converge)
		{
			t++;
		}
		else
		{
			t = 0;
		}
		this->numiteration++;
	}
	//for debuging 
	//std::cout<<"number of iteration = "<<this->numiteration<<std::endl;
	//std::cout<<"number of samples = "<<this->numsampled<<std::endl;
	//std::cout<<"p_hat = "<<this->phat<<std::endl;
	//std::cout<<"delta = "<<sqrt((double)this->varphat)<<std::endl;
	//this->FileWrite();
}

//Check strata whether it need to to be split
bool ActiveValidation::SplitCheck(int kth_strata)
{
	bool split = false;
	if(this->numiteration == 0)
	{
		return split;
	}

	double nk = this->stratas[kth_strata].nk;
	double Nk = this->stratas[kth_strata].Nk;
	double fpcvarphatk = this->stratas[kth_strata].fpcvarphatk;

	double maxfpc = 0.25 * (1 - nk/Nk) / (nk - 1);
	if( (0.5*maxfpc < fpcvarphatk) && (Nk > 0.1* this->N ))
	{
		split = true;
	}
	return split;
}

//Calculate sigmak for each strata
void ActiveValidation::CalculateSigmak(int kth_strata)
{
	double temp = this->stratas[kth_strata].fpcvarphatk;
	this->stratas[kth_strata].sigmak = sqrt((double)temp);
}

//Calculate phatk for each strata
void ActiveValidation::Calculatephatak(int kth_strata)
{
	double pk = 0.5;
	double m;
	double tempphatk;

	if(this->stratas[kth_strata].nk == 0)
	{
		m = 2;
	}
	else
	{
		m = 1 / sqrt((double)this->stratas[kth_strata].nk);
	}

	tempphatk = (this->stratas[kth_strata].hk + m*pk) / (this->stratas[kth_strata].nk + m);
	this->stratas[kth_strata].phatk = tempphatk;
}

//Calculate variance of phatk
void ActiveValidation::CalculateVarphatk(int kth_strata)
{
	if(this->stratas[kth_strata].nk <= 1)
	{
		std::cout<<"nk is less or equal than 1 !!!!!"<<std::endl;
	}
	else
	{
		double tempphatk = this->stratas[kth_strata].phatk;
		double tempvarphatk;
		double tempfpcvarphatk;
		tempvarphatk = tempphatk * (1 - tempphatk) / (this->stratas[kth_strata].nk - 1);
		tempfpcvarphatk = tempvarphatk * (1 - (double)this->stratas[kth_strata].nk / (double)this->stratas[kth_strata].Nk);
		this->stratas[kth_strata].varphatk = tempvarphatk;
		this->stratas[kth_strata].fpcvarphatk = tempfpcvarphatk;
	}
}

//Calculate number of elements to be drawn
void ActiveValidation::CalculatenumToBeDrawn()
{
	double init = 0.02;
	double tempnumtobedrawn;
	if(this->numiteration == 0)
	{
		for(int i = 0; i < this->numbin; i++)
		{
			tempnumtobedrawn = init * this->stratas[i].Nk;
			int temint = floor( (double)tempnumtobedrawn);
			if( (tempnumtobedrawn - temint) >= 0.5)
				temint += 1;
			this->stratas[i].numtobeDrawn = tempnumtobedrawn;
			this->numsampled += tempnumtobedrawn;
		}
	}
	else
	{
		//double tempnum = 0.1 * this->N * 1.96 *sqrt((double)this->varphat);
		//double tempnummax = (std::max<double >)(40, tempnum);
		//double coe = tempnummax / tempnum;

		for(int i = 0; i < this->numbin; i++)
		{
			tempnumtobedrawn = 0.15 * (this->stratas[i].Nk) * 1.96 * (this->stratas[i].sigmak);
			int temint = floor( (double)tempnumtobedrawn);
			if( (tempnumtobedrawn - temint) >= 0.5)
				temint += 1;
				
			/*tempnumtobedrawn = tempnumtobedrawn * coe;*/
			//using floor
			//this->stratas[i].numtobeDrawn = tempnumtobedrawn;
			//this->numsampled += tempnumtobedrawn;
			//using around
			this->stratas[i].numtobeDrawn = temint;
			this->numsampled += temint;
		}
	}
}

//Check whether convergence has been derived
bool ActiveValidation::ConvergeCheck()
{
	bool converge = false;
	
	double co = 1 / sqrt( ( double)this->N);
	double c1 = sqrt( (double)(this->phat / (1 - this->phat)));
	double c2 = sqrt( (double)((1 - this->phat) / this->phat));
	double norapproxcri = co * abs(c1 - c2);
	double var = 1.96 * sqrt((double)this->varphat);
	if((var < this->delta) && norapproxcri < 0.3 )
	{
		converge = true;
	}
	return converge;
}

//Calculate phat
void ActiveValidation::Calculatephata()
{
	double tempphat = 0;
	for(int i = 0; i < this->numbin; i++)
	{
		tempphat += this->stratas[i].phatk * (double)this->stratas[i].Nk / (double)this->N;
	}
	this->phat = tempphat;
	//std::cout<<"phat = "<<this->phat<<std::endl;
}

//Calculate variance of phat
void ActiveValidation::CalculateVarphat()
{
	double tempvarphat = 0;
	for(int i = 0; i < this->numbin; i++)
	{
		double coe = (double)this->stratas[i].Nk / (double)this->N;
		//tempvarphat += this->stratas[i].varphatk  * coe * coe ;
		tempvarphat += this->stratas[i].fpcvarphatk  * coe * coe ;
	}
	this->varphat = tempvarphat;
	//std::cout<<"varphat = "<<this->varphat<<std::endl;
}

//Calculate hk for each strata
void ActiveValidation::Calculatehk(int kth_strata)
{
	int counter = 0;
	for(int i = 0; i < this->stratas[kth_strata].selectingidsk.size(); i++)
	{
		int tempid = this->stratas[kth_strata].selectingidsk[i];
		if(this->stratas[kth_strata].strlable[this->stratas[kth_strata].unselectedidsk [tempid]] == 1)
		{
			counter++;
		}
	}
	this->stratas[kth_strata].hk += counter;
}

//Calculate hk for each strata without labels
void ActiveValidation::CalculatehkwithoutLabel(int kth_strata, std::vector<int > learned)
{
	int counter = 0;
	for(int i = 0; i < learned.size(); i++)
	{
		int tempid = learned[i];
		if(tempid == kth_strata )
		{
			counter++;
		}
	}
	this->stratas[kth_strata].hk += counter;
	std::cout<<"hk = "<<this->stratas[kth_strata].hk<<std::endl;
}

//Update bins after random sampling from each strata
void ActiveValidation::UpdateBins(int kth_strata)
{
	std::vector<int > tempunselectedidsk = this->stratas[kth_strata].unselectedidsk;
	for(int i = 0; i < this->stratas[kth_strata].numtobeDrawn; i++)
	{
		int templocalid = this->stratas[kth_strata].selectingidsk[i];
		this->stratas[kth_strata].selectedidsk.push_back(tempunselectedidsk[templocalid]);
	}

	std::vector<int > tempunselectedidskcurrent;
	for(int p = 0; p < tempunselectedidsk.size(); p++)
	{
		bool exist = false;
		for(int q = 0; q < this->stratas[kth_strata].selectingidsk.size(); q++)
		{
			if(this->stratas[kth_strata].selectingidsk[q] == p)
			{
				exist = true;
				break;
			}
		}
		if(exist == false)
		{
			tempunselectedidskcurrent.push_back(tempunselectedidsk[p]);
		}
	}
	
	this->stratas[kth_strata].unselectedidsk = tempunselectedidskcurrent;
	this->stratas[kth_strata].numleft = tempunselectedidskcurrent.size();
}
void ActiveValidation::FileWrite()
{
	const char* filename = "strabel";
	FILE *fp = fopen(filename,"w");
	for(int i = 0; i < this->numbin; i++)
	{
		fprintf(fp,"%s", "stratalabel");
		fprintf(fp,"\t");
		fprintf(fp,"%d", i);
		fprintf(fp,"\n");
		for(int j = 0; j < this->stratas[i].Nk; j++)
		{
			fprintf(fp,"%d",j);
			fprintf(fp,"\t");
			fprintf(fp,"%d",this->stratas[i].strlable(j));
			fprintf(fp,"\t");
			fprintf(fp,"%d",this->lable (this->stratas[i].idsk[j]));
			fprintf(fp,"\n");
		}
	}
	fclose(fp);
}
//Calculate the true accuracy of entire dataset
void ActiveValidation::TrueAccuracy()
{
	this->trueaccuracy = 0;
	for(int i = 0; i < this->N; i++)
	{
		if(this->lable(i) == 1)
		{
			this->trueaccuracy +=1;
		}
	}

	this->trueaccuracy = this->trueaccuracy / (double)this->N ;
}
//Query and sample from outside
void ActiveValidation::QuearyandSample()
{
	//while(this->t < 2)
	//{
	//	std::vector<int > queries;
	//	std::vector<int > learned;
	//	a_v->CalculatenumToBeDrawn();
	//	for(int i = 0; i < numbin; i++)
	//	{
	//		queries = a_v->RandomSamplingwithoutLabel
	//			(a_v->stratas[i].numtobeDrawn, a_v->stratas[i].numleft, i);

	//		learned = Labellearning(queries, i, numbin);
	//		a_v->CalculatehkwithoutLabel(i, learned);
	//		a_v->Calculatephatak(i);
	//		a_v->CalculateVarphatk(i);
	//		a_v->CalculateSigmak(i);
	//		a_v->UpdateBins(i);
	//	}

	//	//calculate result and check converge
	//	a_v->Calculatephata();
	//	a_v->CalculateVarphat();
	//	converge = a_v->ConvergeCheck();

	//	//avoid spurious
	//	if(converge)
	//	{
	//		t++;
	//	}
	//	else
	//	{
	//		t = 0;
	//	}
	//	a_v->numiteration++;
	//}
}