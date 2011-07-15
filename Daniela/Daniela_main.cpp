
#if defined(_MSC_VER)
#pragma warning (disable : 4786)
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif
#
#include <iostream>
#include <math.h>
#include "Daniela.h"
using namespace std;

int main(int argc, char * argv[])
{
	//argv[1] = "C:\\Lab\\SuperBuild\\Bin\\exe\\Release\\s05-68478.tif";
	
	if( argc < 3 )
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0];
		std::cerr << "daniela.exe <Input File Name> <Number of classes> " << std::endl;
		return EXIT_FAILURE;
	}

	//parse command line arguments
	std::string imageFileName = std::string(argv[1]);
	
	int numberOfInitialClasses = atoi( argv[2] );
	
	//Construct Daniela
	Daniela* myDaniela = new Daniela(imageFileName, numberOfInitialClasses);
	
	//Calculate labeled bone image
	myDaniela->Get_Daniela();

	return 0;
}
