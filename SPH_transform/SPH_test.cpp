#include "Spherical_Harmonic_transform.h"

#include <fstream>
#include <stdio.h>

std::vector<std::vector<double> > ReadFile(const char *filename);

int  main(int argc, char * argv [])
{
	std::vector<std::vector<double> > data;
	std::cout<<argv [1]<<std::endl;
	const char* filename = "949.txt";
	data = ReadFile("207.txt");

	int L = 20;
	SPH_Trasform* sph_trans =    new SPH_Trasform();
	sph_trans->Initialize(data, L);
	sph_trans->Transform();

	sph_trans->WriteFile();
	return 0;
}

std::vector<std::vector<double> > ReadFile(const char *filename)
{
	FILE *fp = fopen(filename,"r");
	int num_samples = 0;
	int num_features = 0;
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

	std::vector<std::vector<double> > tempdata;
	for(int i=0; i<num_samples; i++)
	{
		std::vector<double > tempvector;
		double temp;

		for(int j=0; j<num_features; j++)
		{
			fscanf(fp, "%lf", &temp);
			tempvector.push_back(temp);
		}
		tempdata.push_back(tempvector);
	}

	std::cout<<"Data file loaded !!!"<<std::endl;
	return tempdata;	
}

