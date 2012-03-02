#include "SomaExtraction.h"

int main(int argc, char* argv[])
{	
	clock_t SomaExtraction_start_time = clock();

	if(argc < 7)
	{
		std::cout<<"Usage: SomaExtraction <InputImageFileName> <SomaCentroids.txt> <SomaImage.tif>\
			<Time Threshold (typically 100)> <Curvature Scaling (typically 0.5)> <RMS Error (typically 0.02)>\n";
		return 0;
	}
	
	//if(argc < 8)
	//{
	//	std::cout<<"Usage: SomaExtraction <InputImageFileName> <SomaCentroids.txt> <ParamFile> <SomaImage.mhd> <open morph filter size (typically) 8> <minimum object size (typically 1500)> <number of bins>\n";
	//	return 0;
	//}

	SomaExtractor *Somas = new SomaExtractor();
	
	std::cout << "Entering SetInputImage" << std::endl;
	Somas->SetInputImage(argv[1]);
	
	Somas->SegmentSoma(argv[2], argv[3], atoi(argv[4]), atof(argv[5]), atof(argv[6]));

	//GenerateSeedPoints(argv[3], atoi(argv[7]));
	//std::cout << "Entering LoadSegParams" << std::endl;
	//Somas->LoadSegParams(atoi(argv[5]), atoi(argv[6]));
	
	//std::cout << "Entering binarizeImage" << std::endl;
	//Somas->binarizeImage(argv[3], atoi(argv[7]));

	//std::cout << "Entering relabelBinaryImage" << std::endl;
	//Somas->relabelBinaryImage();
	//
	//std::cout << "Entering GetSomaCentroids" << std::endl;
	//Somas->GetSomaCentroids();

	//std::cout << "Entering writeSomaCentroids" << std::endl;
	//Somas->writeSomaCentroids(argv[2]);

	//std::cout << "Entering writeSomaImage" << std::endl;
	//Somas->writeSomaImage(argv[4]);

	std::cout << "Total time for SomaExtraction is: " << (clock() - SomaExtraction_start_time) / (float) CLOCKS_PER_SEC << std::endl;

	return 0;
}