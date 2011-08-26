#include "MultipleNeuronTracer.h"
#include "time.h"

int main(int argc, char* argv[])
{	

    /*argv[1] = "C:\\ROYSAMLAB\\FARSIGHT\\BINARY\\exe\\Debug\\R2080_6wk_Crop_sigma_0.030000_CV.tif";
	argv[2] = "C:\\ROYSAMLAB\\FARSIGHT\\BINARY\\exe\\Debug\\11111111.tif";
	argv[3] = "C:\\ROYSAMLAB\\FARSIGHT\\BINARY\\exe\\Debug\\ParameterFile.ini";
	argv[4] = "C:\\ROYSAMLAB\\FARSIGHT\\BINARY\\exe\\Debug\\SegParams.ini";*/
	
	if (argc < 4)
	{
		std::cout << "MultipleNeuronTracer.exe <InputFileName> <SomaCentroids.txt> <CostThreshold (Usually 1000)> <SomaImageFileName>" << std::endl;
		return -1;
	}
	clock_t start_time = clock();
	std::string InputFilename = std::string(argv[1]);
	std::string SWCFilename = InputFilename;
	SWCFilename.erase(SWCFilename.length()-4,SWCFilename.length());
	SWCFilename.append("_ANT.swc");

	MultipleNeuronTracer * MNT = new MultipleNeuronTracer();
	clock_t LoadCurvImage_start_time = clock();
	MNT->LoadCurvImage(std::string(argv[1]), 1);
	std::cout << "LoadCurvImage took: " << (clock() - LoadCurvImage_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	clock_t ReadStartPoints_start_time = clock();
	MNT->ReadStartPoints(std::string(argv[2]), 1);
	std::cout << "ReadStartPoints took: " << (clock() - ReadStartPoints_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	clock_t SetCostThreshold_start_time = clock();
	MNT->SetCostThreshold(atof(argv[3]));
	std::cout << "SetCostThreshold took: " << (clock() - SetCostThreshold_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	//clock_t LoadSomaImage_start_time = clock();
	//MNT->LoadSomaImage(std::string(argv[4]));
	//std::cout << "LoadSomaImage took: " << (clock() - LoadSomaImage_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	clock_t RunTracing_start_time = clock();
	MNT->RunTracing();
	std::cout << "RunTracing took: " << (clock() - RunTracing_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	MNT->WriteSWCFile(SWCFilename, 1);

	std::cout << "Total time to tracing is : " << (clock() - start_time)/(float) CLOCKS_PER_SEC << std::endl;
	

	return 0;

}