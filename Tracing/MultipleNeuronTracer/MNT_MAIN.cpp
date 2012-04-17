#include "MultipleNeuronTracer.h"
#include "time.h"

int main(int argc, char* argv[])
{	

    /*argv[1] = "C:\\Lab\\data\\forMultineuron\\montage_8bitkt11156_w311GFPdsu.raw_Tiles\\Chan3_3_2100_0.tif"; //"C:\\ROYSAMLAB\\FARSIGHT\\BINARY\\exe\\Debug\\R2080_6wk_Crop_sigma_0.030000_CV.tif";
	argv[2] = "C:\\Lab\\data\\forMultineuron\\montage_8bitkt11156_w311GFPdsu.raw_Tiles\\MG_centroids_3_3_2100_0.txt"; //"C:\\ROYSAMLAB\\FARSIGHT\\BINARY\\exe\\Debug\\11111111.tif";
	argv[3] = "1000"; //"C:\\ROYSAMLAB\\FARSIGHT\\BINARY\\exe\\Debug\\ParameterFile.ini";
	argv[4] = "C:\\Lab\\data\\forMultineuron\\montage_8bitkt11156_w311GFPdsu.raw_Tiles\\op_3_3_2100_0.tif"; //"C:\\ROYSAMLAB\\FARSIGHT\\BINARY\\exe\\Debug\\SegParams.ini";
	*/

	if(argc < 4 || argc > 5)
	{
		std::cout << "MultipleNeuronTracer.exe <InputFileName> <SomaCentroids.txt> <SomaImageFile> [options file]" << std::endl;
		return -1;
	}
	std::cout<<"argc="<<argc<<std::endl;
	MultipleNeuronTracer * MNT = new MultipleNeuronTracer();

	// code for loading parameters from txt file
	std::map<std::string, std::string> opts;  
	if(argc == 5)
		MNT->optionsCreate(argv[4], opts);
	
	std::map<std::string,std::string>::iterator mi;

	mi = opts.find("-intensity_threshold"); 
	if(mi!=opts.end())
	{ std::istringstream ss((*mi).second); ss>>MNT->intensity_threshold; 
	}
	else
	{ MNT->intensity_threshold = 0.005; printf("Chose intensity_threshold = 0.005 as default\n");}

	mi = opts.find("-contrast_threshold");
	if(mi!=opts.end())
	{ std::istringstream ss((*mi).second); ss>>MNT->contrast_threshold; }
	else
	{	  MNT->contrast_threshold = 0.0003; printf("Chose contrast_threshold = 0.0003 as default\n"); }

	mi = opts.find("-cost_threshold"); 
	if(mi!=opts.end())
	{ std::istringstream ss((*mi).second); ss>>MNT->cost_threshold; }
	else
	{ MNT->cost_threshold = 700; printf("Chose cost_threshold = 700 as default\n");}

	mi = opts.find("-debris_threshold"); 
	if(mi!=opts.end())
	{ std::istringstream ss((*mi).second); ss>>MNT->debris_threshold; }
	else
	{ MNT->debris_threshold = 0.8; printf("Chose debris_threshold = 0.8 as default\n"); }

	mi = opts.find("-offshoot"); 
	if(mi!=opts.end())
	{ std::istringstream ss((*mi).second); ss>>MNT->offshoot; }
	else
	{ MNT->offshoot = 10; printf("Chose offshoot = 10 as default\n"); }

	mi = opts.find("-device"); 
	if(mi!=opts.end())
	{ std::istringstream ss((*mi).second); ss>>MNT->device; }
	else
	{ MNT->device = 1; printf("Chose device = 0 as default\n"); }

	std::cout<<"intensity_threshold="<<MNT->intensity_threshold<<std::endl;
	std::cout<<"contrast_threshold="<<MNT->contrast_threshold<<std::endl;
	std::cout<<"cost_threshold="<<MNT->cost_threshold<<std::endl;
	std::cout<<"debris_threshold="<<MNT->debris_threshold<<std::endl;
	std::cout<<"offshoot="<<MNT->offshoot<<std::endl;
	std::cout<<"device="<<MNT->device<<std::endl;


	// end of code for loading parameters from txt file


	clock_t start_time = clock();
	std::string InputFilename = std::string(argv[1]);
	std::string SWCFilename = InputFilename;
	SWCFilename.erase(SWCFilename.length()-4,SWCFilename.length());
	SWCFilename.append("_ANT.swc");

	//MultipleNeuronTracer * MNT = new MultipleNeuronTracer();
	clock_t LoadCurvImage_start_time = clock();
	MNT->LoadCurvImage(std::string(argv[1]), 1);
	std::cout << "LoadCurvImage took: " << (clock() - LoadCurvImage_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	clock_t ReadStartPoints_start_time = clock();
	MNT->ReadStartPoints(std::string(argv[2]), 1);
	std::cout << "ReadStartPoints took: " << (clock() - ReadStartPoints_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	clock_t SetCostThreshold_start_time = clock();
	MNT->SetCostThreshold(MNT->cost_threshold);//instead of the previous argv[3]
	std::cout << "SetCostThreshold took: " << (clock() - SetCostThreshold_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	clock_t LoadSomaImage_start_time = clock();
	MNT->LoadSomaImage(std::string(argv[3]));
	std::cout << "LoadSomaImage took: " << (clock() - LoadSomaImage_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	clock_t RunTracing_start_time = clock();
	MNT->RunTracing();
	std::cout << "RunTracing took: " << (clock() - RunTracing_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	MNT->WriteSWCFile(SWCFilename, 1);

	//MNT->WriteMultipleSWCFiles(SWCFilename, 1);

	std::cout << "Total time to tracing is : " << (clock() - start_time)/(float) CLOCKS_PER_SEC << std::endl;
	

	return 0;

}