#include "MultipleNeuronTracer.h"

int main(int argc, char* argv[])
{	

    /*argv[1] = "C:\\ROYSAMLAB\\FARSIGHT\\BINARY\\exe\\Debug\\R2080_6wk_Crop_sigma_0.030000_CV.tif";
	argv[2] = "C:\\ROYSAMLAB\\FARSIGHT\\BINARY\\exe\\Debug\\11111111.tif";
	argv[3] = "C:\\ROYSAMLAB\\FARSIGHT\\BINARY\\exe\\Debug\\ParameterFile.ini";
	argv[4] = "C:\\ROYSAMLAB\\FARSIGHT\\BINARY\\exe\\Debug\\SegParams.ini";*/
	
	clock_t start_time = clock();
	std::string InputFilename = std::string(argv[1]);
	std::string SWCFilename = InputFilename;
	SWCFilename.erase(SWCFilename.length()-4,SWCFilename.length());
	SWCFilename.append("_ANT.swc");

	MultipleNeuronTracer * MNT = new MultipleNeuronTracer();
	MNT->LoadCurvImage(std::string(argv[1]), 1);
	MNT->ReadStartPoints(std::string(argv[2]), 1);
	MNT->SetCostThreshold(atof(argv[3]));
	MNT->LoadSomaImage(std::string(argv[4]));
	MNT->RunTracing();
	MNT->WriteSWCFile(SWCFilename, 1);

	std::cout << "Total time to segmentation is : " << (clock() - start_time)/(float) CLOCKS_PER_SEC << std::endl;
	

	return 0;

}