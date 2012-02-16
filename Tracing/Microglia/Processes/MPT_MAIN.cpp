#include "MicrogliaProcessTracer.h"

int main(int argc, char* argv[])
{	
	if(argc < 6)
	{
		std::cout << argv[0] << " input-image seedfile.txt soma-image distance-threshold output.swc [spacingXY spacingZ]" << std::endl;
		return 1;
	}

	clock_t start_time = clock();
	std::string InputFilename = std::string(argv[1]);
	std::string SWCFilename = InputFilename;
	SWCFilename.erase(SWCFilename.length()-4,SWCFilename.length());
	SWCFilename.append("_ANT.swc");

	MicrogliaProcessTracer * MPT = new MicrogliaProcessTracer();
  
  //check if image spacing was provided on command-line
  if( argc > 6)
    {
    MPT->SetXSpacing( atof(argv[6]) );
    MPT->SetZSpacing( atof(argv[7]) );
    }

	std::cout << "Loading the input image" << std::endl;
	MPT->LoadInputImage(std::string(argv[1]));

	std::cout << "Reading the start points" << std::endl;
	MPT->LoadSeedPoints(std::string(argv[2]));

	std::cout << "Reading the soma image" << std::endl;
	MPT->LoadSomaImage(std::string(argv[3]));

	double threshold = atof(argv[4]);
	std::cout << "Distance Threshold: " << threshold << std::endl;
	MPT->SetMaxDistance( threshold );

	std::cout << "Beginning Tracing" << std::endl;
	MPT->RunTracing();

	MPT->SetSeparateFilePerCell(false);
	MPT->WriteToSWC( std::string(argv[5]) );

	std::cout << "Total time to segmentation is : " << (clock() - start_time)/(float) CLOCKS_PER_SEC << std::endl;

	delete MPT;
	return 0;
}
