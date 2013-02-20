#include <time.h>
#include "Haralick.h"
#include "itkImageFileReader.h"
# include "omp.h"

using namespace std;



int main(int argc, char * argv [])
{
	
	
	//argv[1] = "C:\\SACHIN\\Farsight_bin\\exe\\Release\\sachin61.tif";
	if( argc < 2 )
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0];
		std::cerr << "Haralick_test.exe <Input File Name> " << std::endl;
		return EXIT_FAILURE;
	}
	std::string imageFileName = argv[1];
    Haralick* myHaralick =    new Haralick(imageFileName);
	double* features=(double *) malloc(13 * sizeof(double)); 
	clock_t start_time = clock();

	//int i=0;
	//#pragma omp parallel for
	for (int i = 0; i < 1; i++)
	{
		features=myHaralick->GetHaralik();	
		for(int k=0;k<14;++k)
		{
			cout<<features[k]<<"   ";

		}
	}        



	cout << "Total time taken: " << ( (clock() - start_time)/(float) CLOCKS_PER_SEC) << endl; 
	return 0;
}
