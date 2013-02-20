#include "MultipleNeuronTracer.h"
#include <ctime>

int main(int argc, char* argv[])
{	

    /*argv[1] = "C:\\Lab\\data\\forMultineuron\\montage_8bitkt11156_w311GFPdsu.raw_Tiles\\Chan3_3_2100_0.tif"; //"C:\\ROYSAMLAB\\FARSIGHT\\BINARY\\exe\\Debug\\R2080_6wk_Crop_sigma_0.030000_CV.tif";
	argv[2] = "C:\\Lab\\data\\forMultineuron\\montage_8bitkt11156_w311GFPdsu.raw_Tiles\\MG_centroids_3_3_2100_0.txt"; //"C:\\ROYSAMLAB\\FARSIGHT\\BINARY\\exe\\Debug\\11111111.tif";
	argv[3] = "1000"; //"C:\\ROYSAMLAB\\FARSIGHT\\BINARY\\exe\\Debug\\ParameterFile.ini";
	argv[4] = "C:\\Lab\\data\\forMultineuron\\montage_8bitkt11156_w311GFPdsu.raw_Tiles\\op_3_3_2100_0.tif"; //"C:\\ROYSAMLAB\\FARSIGHT\\BINARY\\exe\\Debug\\SegParams.ini";
	*/

	if(argc < 4 || argc > 6)
	{
		std::cout << "MultipleNeuronTracer.exe <InputFileName> <SomaCentroids.txt> <SomaImageFile> [options file] [Optimize Coverage (0 or 1)]" << std::endl;
		return -1;
	}

	std::cout<<"argc="<<argc<<std::endl;

	int opt_coverage = atoi(argv[5]);

	MultipleNeuronTracer * MNT = new MultipleNeuronTracer();

	//Setting _indxDice to zero, when MNT is executed as standalone
	itk::Index<3> zero_index;
	zero_index[0]=0;
	zero_index[1]=0;
	zero_index[2]=0;
	MNT->setDiceIndex(zero_index);

	clock_t LoadCurvImage_start_time = clock();
	MNT->LoadCurvImage(std::string(argv[1]), 1);//
	std::cout << "LoadCurvImage took: " << (clock() - LoadCurvImage_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	clock_t ReadStartPoints_start_time = clock();
	MNT->ReadStartPoints(std::string(argv[2]), 1);//
	std::cout << "ReadStartPoints took: " << (clock() - ReadStartPoints_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	clock_t LoadSomaImage_start_time = clock();
	MNT->LoadSomaImage(std::string(argv[3]));
	std::cout << "LoadSomaImage took: " << (clock() - LoadSomaImage_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	 
	std::string coverageFileName = std::string(argv[1]);
	coverageFileName.erase(coverageFileName.length()-4, coverageFileName.length());
	coverageFileName.append("_coverage.txt");
	
	MNT->Set_isCoverageOptimized(false);

	if(opt_coverage == 1)
		MNT->OptimizeCoverage(coverageFileName, true);

	clock_t LoadParameters_start_time = clock();
	MNT->LoadParameters(argv[4], argc);
	std::cout << "LoadParameters took: " << (clock() - LoadParameters_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	clock_t start_time = clock();
	std::string InputFilename = std::string(argv[1]);
	std::string SWCFilename = InputFilename;
	SWCFilename.erase(SWCFilename.length()-4,SWCFilename.length());
	SWCFilename.append("_ANT.swc");

	//MultipleNeuronTracer * MNT = new MultipleNeuronTracer();
	
	clock_t SetCostThreshold_start_time = clock();
	MNT->SetCostThreshold(MNT->cost_threshold);//instead of the previous argv[3]
	std::cout << "SetCostThreshold took: " << (clock() - SetCostThreshold_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	
	clock_t RunTracing_start_time = clock();
	int tracing_type = MNT->tracing_type;
	if(tracing_type == 1){
		MNT->RunTracing();
	}else{
		MNT->RunGVFTracing(false);
	}
	
	//
	std::cout << "RunTracing took: " << (clock() - RunTracing_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	MNT->WriteSWCFile(SWCFilename, 1);

	//MNT->WriteMultipleSWCFiles(SWCFilename, 1);

	std::cout << "Total time to tracing is : " << (clock() - start_time)/(float) CLOCKS_PER_SEC << std::endl;
	

	return 0;

}
