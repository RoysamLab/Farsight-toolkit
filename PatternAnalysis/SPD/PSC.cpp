#include "SPDAnalysisModel.h"
#include <vnl/vnl_vector.h>
#include <vul/vul_arg.h>
#include <fstream>

int main(int argc, char* argv[])
{
	vul_arg< vcl_string > arg_file(0, "A file containing x and y with headers(x,y), seperated by','");
	vul_arg< int > arg_bin("-bin", "bin size for the histogram of the full graph distance metric", 20);
	vul_arg< int > arg_k("-neighbor", "k for k-NNG", 4);
	vul_arg< bool > arg_debug("-debug", "Print out the histgram information", false);
	vul_arg_parse(argc, argv);

	ifstream file(arg_file().c_str());
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

	double ps = 0;
	double ps2 = 0;

	ps = SPDAnalysisModel::CaculatePSC(arg_k(), arg_bin(),vecx, vecy, arg_debug());
	ps2 = SPDAnalysisModel::CaculatePSC(arg_k(), arg_bin(),vecy, vecx, arg_debug());

	std::ofstream ofs("PSC-result.csv");
	ofs<< ps<<","<<ps2<<std::endl;
	ofs.close();
	std::cout<< "PS: "<<ps<<"\t"<<ps2<<std::endl;
}