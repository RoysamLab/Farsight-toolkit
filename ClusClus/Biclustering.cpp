#include "Biclustering.h"

Bicluster::Bicluster()
{
	this->num_iteration = 9;
	this->K1 = 15;
	this->K2 = 25;
}

Bicluster::~Bicluster()
{
}

void Bicluster::setDataToBicluster(std::vector<std::vector<double > > & data)
{
	if(data.empty())
	{
		std::cout<<"No data to biclustering !......."<<std::endl;
		return;
	}
	this->data = data;
	this->num_rows = data.size();
	this->num_cols = data[0].size();
}

//void Bicluster::bispectralclustering()
//{
//	int iteration = 3;
//	this->Data1 = data;
//	std::vector<std::vector<double > > dataT;
//	dataT = this->transpose(data);
//
//	for(int i = 0; i < iteration; i++)
//	{
//		//compute local geometry reference
//		LocalGeometryRef* lgr1 = new LocalGeometryRef();
//		lgr1->Initialize(this->Data1);
//		lgr1->ComputeSimilarityMatrix();
//		lgr1->ComputeProbabilityMatrix();
//		lgr1->SVD(20);
//		this->Data1 = lgr1->GetEigenVectors(true);
//		delete lgr1;
//
//		//compute tree structure
//		BiTree* tree1 = new BiTree();
//		tree1->setDataToTree(Data1,Data1);
//		tree1->setFirstIterationClusterNumber(10);
//		if(i == iteration - 1)
//			tree1->setReorderFlag(true);
//		tree1->constructTree();
//		this->levels1 = tree1->getLevels();
//		delete tree1;
//		
//		//get order after clustering
//		this->reOrder(this->levels1, 1);
//		const int num1 = this->qorder1.size();
//		for(int jj = 0; jj < num1; ++jj)
//		{
//			int temp = this->qorder1.front();
//			this->order1.push_back(temp);
//			this->qorder1.pop();
//		}
//
//		//reorganize data for the other direction
//		this->Data2 = reorganize(data, order1);
//
//		//compute local geometry reference for the other direction
//		LocalGeometryRef* lgr2 = new LocalGeometryRef();
//		lgr2->Initialize(this->Data2);
//		lgr2->ComputeSimilarityMatrix();
//		lgr2->ComputeProbabilityMatrix();
//		lgr2->SVD(10);
//		this->Data2 = lgr2->GetEigenVectors(true);
//		delete lgr2;
//
//		//compute tree structure
//		BiTree* tree2 = new BiTree();
//		tree2->setDataToTree(Data2,Data2);
//		tree2->setFirstIterationClusterNumber(10);
//		if(i == iteration - 1)
//			tree2->setReorderFlag(true);
//		tree2->constructTree();
//		this->levels2 = tree2->getLevels();
//		delete tree2;
//
//		//get order after clustering		
//		this->reOrder(this->levels2, 2);
//		const int num2 = this->qorder2.size();
//		for(int jj = 0; jj < num2; ++jj)
//		{
//			int temp = this->qorder2.front();
//			this->order2.push_back(temp);
//			this->qorder2.pop();
//		}
//
//		//reorganize data for the other direction
//		this->Data1 = reorganize(dataT, order2);
//	}
//}

void Bicluster::bispectralclustering()
{
	this->Data1 = data;
	std::vector<std::vector<double > > dataT;
	dataT = this->transpose(data);
	for(int i = 0; i < this->num_iteration; i++)
	{
		Data1 = this->transpose(Data1);
		BiTree* tree1 = new BiTree();
		tree1->setDataToTree(dataT,Data1);
		tree1->setFirstIterationClusterNumber(this->K1);
		if(i == this->num_iteration - 1)
			tree1->setReorderFlag(true);
		tree1->constructTree();
		this->Data2 = tree1->getPartitionTree();
		this->levels2 = tree1->getLevels();
		delete tree1;

		Data2 = this->transpose(Data2);
		BiTree* tree2 = new BiTree();
		tree2->setDataToTree(data,Data2);
		tree2->setFirstIterationClusterNumber(this->K2);
		if(i == this->num_iteration - 1)
			tree2->setReorderFlag(true);
		tree2->constructTree();
		this->Data1 = tree2->getPartitionTree();
		this->levels1 = tree2->getLevels();
		delete tree2;
		std::cout<<"****..."<<i<<"th...iteration"<<std::endl;
	}
	this->reOrder(this->levels1, 1);
	this->reOrder(this->levels2, 2);

	int num1 = this->qorder1.size();
	for(int i = 0; i<num1; i++)
	{
		int temp = this->qorder1.front();
		this->order1.push_back(temp);
		this->qorder1.pop();
	}
	int num2 = this->qorder2.size();
	for(int i = 0; i < num2; i++)
	{
		int temp = this->qorder2.front();
		this->order2.push_back(temp);
		this->qorder2.pop();
	}
}

void Bicluster::biclustering()
{
	this->Data1 = data;
	std::vector<std::vector<double > > dataT;
	dataT = this->transpose(data);
	for(int i = 0; i < this->num_iteration; i++)
	{
		Data1 = this->transpose(Data1);
		BiTree* tree1 = new BiTree();
		tree1->setDataToTree(dataT,Data1);
		tree1->setFirstIterationClusterNumber(this->K1);
		if(i == this->num_iteration - 1)
			tree1->setReorderFlag(true);
		tree1->constructTree();
		this->Data2 = tree1->getPartitionTree();
		this->levels2 = tree1->getLevels();
		delete tree1;

		Data2 = this->transpose(Data2);
		BiTree* tree2 = new BiTree();
		tree2->setDataToTree(data,Data2);
		tree2->setFirstIterationClusterNumber(this->K2);
		if(i == this->num_iteration - 1)
			tree2->setReorderFlag(true);
		tree2->constructTree();
		this->Data1 = tree2->getPartitionTree();
		this->levels1 = tree2->getLevels();
		delete tree2;
		std::cout<<"****..."<<i<<"th...iteration"<<std::endl;
	}
	this->reOrder(this->levels1, 1);
	this->reOrder(this->levels2, 2);

	int num1 = this->qorder1.size();
	for(int i = 0; i<num1; i++)
	{
		int temp = this->qorder1.front();
		this->order1.push_back(temp);
		this->qorder1.pop();
	}
	int num2 = this->qorder2.size();
	for(int i = 0; i < num2; i++)
	{
		int temp = this->qorder2.front();
		this->order2.push_back(temp);
		this->qorder2.pop();
	}
}

std::vector<std::vector<double > > Bicluster::transpose(std::vector<std::vector<double > > & M)
{
	std::vector<std::vector<double > > transposeddata;
	int m = M.size();
	int n = M[0].size();
	transposeddata.resize(n);
	for(int i = 0; i < n; i++)
	{
		transposeddata[i].resize(m);
		for(int j = 0; j < m; j++)
			transposeddata[i][j] = M[j][i];
	}
	return transposeddata;
}

void Bicluster::reOrder(Level_id levels, int flag)
{
	if(flag == 1)
	{
		int level = levels.numm_levels;
		for(int i = 0; i < levels.level_id[level][0].size(); i++)
			this->qorder1.push(levels.level_id[level][0][i]);

		for(int counter = 1; counter < level; counter++)
		{
			int num = this->qorder1.size();
			for(int k = 0; k < num; k++)
			{
				int temp = this->qorder1.front();
				this->qorder1.pop();
				for(int i = 0; i < levels.level_id[level - counter][temp].size(); i++)
					this->qorder1.push(levels.level_id[level - counter][temp][i]);
			}
		}
	}
	else if(flag == 2)
	{
		int level = levels.numm_levels;
		for(int i = 0; i < levels.level_id[level][0].size(); i++)
			this->qorder2.push(levels.level_id[level][0][i]);

		for(int counter = 1; counter < level; counter++)
		{
			int num = this->qorder2.size();
			for(int k = 0; k < num; k++)
			{
				int temp = this->qorder2.front();
				this->qorder2.pop();
				for(int i = 0; i < levels.level_id[level - counter][temp].size(); i++)
					this->qorder2.push(levels.level_id[level - counter][temp][i]);
			}
		}
	}
}

std::vector<std::vector<double > > Bicluster::reorganize(std::vector<std::vector<double > > datatoreorganize, std::vector<int > order)
{
	//reorder data
	std::vector<std::vector<double > > redata;
	for(int row = 0; row < datatoreorganize.size(); row++)
	{
		redata.push_back(datatoreorganize[order[row]]);
	}

	//transpose
	std::vector<std::vector<double > > trandata;
	for(int col = 0; col < datatoreorganize[0].size(); col++)
	{
		std::vector<double > temp;
		for(int row = 0; row < datatoreorganize.size(); row++)
		{
			temp.push_back(redata[row][col]);
		}
		trandata.push_back(temp);
	}
	return trandata;
}

void Bicluster::WriteFile(const char *filename1, const char *filename2)
{
	FILE *fp1 = fopen(filename1,"w");
	for(int i = 0; i < this->order1.size(); i++)
		fprintf(fp1,"%d\t",this->order1[i]);			
	fclose(fp1);

	FILE *fp2 = fopen(filename2,"w");
	for(int i = 0; i < this->order2.size(); i++)
		fprintf(fp2,"%d\t",this->order2[i]);			
	fclose(fp2);
}

void Bicluster::Setparameter(int k1, int k2)
{
	this->K1 = k1;
	this->K2 = k2;
}