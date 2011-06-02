#include "SomaExtraction.h"

int main(int argc, char* argv[])
{	

	//if(argc < 4)
	//{
	//	std::cout<<"Usage: SomaExtraction <InputImageFileName> <OutputImageFileName> <ParametersFileName> <SegParams>\n";
	//	return 0;
	//}


    /*argv[1] = "C:\\ROYSAMLAB\\FARSIGHT\\BINARY\\exe\\Debug\\R2080_6wk_Crop_sigma_0.030000_CV.tif";
	argv[2] = "C:\\ROYSAMLAB\\FARSIGHT\\BINARY\\exe\\Debug\\11111111.tif";
	argv[3] = "C:\\ROYSAMLAB\\FARSIGHT\\BINARY\\exe\\Debug\\ParameterFile.ini";
	argv[4] = "C:\\ROYSAMLAB\\FARSIGHT\\BINARY\\exe\\Debug\\SegParams.ini";*/
	SomaExtractor * Somas = new SomaExtractor();
	Somas->SetInputImage(argv[1]);
	Somas->LoadSegParams(atoi(argv[4]), atoi(argv[5]));
	Somas->binarizeImage(argv[3]);
	Somas->relabelBinaryImage();
	Somas->GetSomaCentroids();
	Somas->writeSomaCentroids(argv[2]);
	Somas->writeSomaImage(argv[2]);

	return 0;

}