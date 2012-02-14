#include "Biclustering.h"

Bicluster::Bicluster()
{
	this->num_iteration = 4;
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
		tree1->setFirstIterationClusterNumber(15);
		if(i == this->num_iteration - 1)
			tree1->setReorderFlag(true);
		tree1->constructTree();
		this->Data2 = tree1->getPartitionTree();
		this->levels2 = tree1->getLevels();
		delete tree1;

		Data2 = this->transpose(Data2);
		BiTree* tree2 = new BiTree();
		tree2->setDataToTree(data,Data2);
		tree2->setFirstIterationClusterNumber(25);
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