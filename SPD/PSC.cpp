#include "SPDAnalysisModel.h"
#include <vnl/vnl_vector.h>
#include <fstream>

int main(int argc, char* argv[])
{
	unsigned int k = 3;
	unsigned int nbin = 20;
	bool bnoise = true;
	if( argc < 2)
	{
		std::cout<< "PSC: PSC <Table of 2 columns> <K for KNN(default 3)> <Number of bins(default 20)> < noise range 1 or fNNG range 0(default 1)>"<<std::endl;
		return 0;
	}
	else if( argc == 3)
	{
		k = atoi(argv[2]);
	}
	else if( argc == 4)
	{
		k = atoi(argv[2]);
		nbin = atoi(argv[3]);
	}
	else if( argc == 5)
	{
		k = atoi(argv[2]);
		nbin = atoi(argv[3]);
		if( atoi(argv[4]) < 1e-6)
		{
			bnoise = false;
		}
	}

	ifstream file(argv[1]);
	unsigned int n = 0;
	char str[100];
	while(!file.eof())
	{
		file.getline(str, sizeof(str));
		if( std::string(str).size() > 0)
		{
			n++;
		}
	}
	std::cout<< "Read in "<< n<<" lines."<<std::endl;
	if( n <= 1)
	{
		return -1;
	}

    double x,y;
	char delim;
	vnl_vector<double> vecx(n - 1);
	vnl_vector<double> vecy(n - 1);

	file.clear();
	file.seekg(0, std::ios_base::beg);
	file.getline(str, sizeof(str));   // remove the header

	unsigned int i = 0;
    while(!file.eof())
    {
        file >> x;
        file >> delim;
        file >> y;

        vecx[i] = x;
        vecy[i] = y;
		i++;
    }


	double ps = SPDAnalysisModel::CaculatePS(bnoise, k, nbin,vecx, vecy);
	double ps2 = SPDAnalysisModel::CaculatePS(bnoise, k, nbin, vecy, vecx);
	std::ofstream ofs("PSC-result.csv");
	ofs<< ps<<","<<ps2<<std::endl;
	ofs.close();
	std::cout<< "PS: "<<ps<<"\t"<<ps2<<std::endl;
}