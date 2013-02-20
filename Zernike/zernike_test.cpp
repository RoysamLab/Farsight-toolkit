
#include "zernike.h"
#include <time.h>
#include "itkImageFileReader.h"
using namespace std;



int main(int argc, char * argv [])
{

	
	//argv[1] = "C:\\SACHIN\\Farsight_bin\\exe\\Release\\test1.tif";
	if( argc < 3 )
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0];
		std::cerr << "zernike_test.exe <Input File Name> <order of zernike moments>  " << std::endl;
		return EXIT_FAILURE;
	}
	std::string imageFileName = argv[1];
	int orderofmoments = atoi(argv[2]);
	zernike* myzernike =    new zernike( imageFileName,orderofmoments);
	//now lets say we have the input ITK image pointer inputImage
	//zernike* myzernike =    new zernike(ImageType::Pointer inputImage,int orderofmoments);
		
	
	clock_t start_time = clock();

    myzernike-> GetZernike();

	cout << "Total time taken: " << ( (clock() - start_time)/(float) CLOCKS_PER_SEC) << endl; 
	return 0;


}

