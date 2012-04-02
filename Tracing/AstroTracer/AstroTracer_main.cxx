#include "AstroTracer.h"
#include "time.h"

int main(int argc, char* argv[])
//int main(void)
{	
	
    argv[1] = "C:\\Prathamesh\\Astrocytes\\Cropped_Experiment\\CroppedAstroTRITC.tif"; //"C:\\ROYSAMLAB\\FARSIGHT\\BINARY\\exe\\Debug\\R2080_6wk_Crop_sigma_0.030000_CV.tif";
	argv[2] = "C:\\Prathamesh\\Astrocytes\\Cropped_Experiment\\roots_cropped.txt"; //"C:\\ROYSAMLAB\\FARSIGHT\\BINARY\\exe\\Debug\\11111111.tif";
	argv[4] = "C:\\Prathamesh\\Astrocytes\\Cropped_Experiment\\CroppedAstroLABEL.tif"; //"C:\\ROYSAMLAB\\FARSIGHT\\BINARY\\exe\\Debug\\SegParams.ini";
	
	argv[3] = "400"; //"C:\\ROYSAMLAB\\FARSIGHT\\BINARY\\exe\\Debug\\ParameterFile.ini";

	argc = 5;

	bool justComputeRootFeatures = true; //false;
	bool startWithCandidateRoots = false;

	std::string tracePointsFileName = "C:\\Prathamesh\\Astrocytes\\Cropped_Experiment\\trace_points.txt";
	std::string featureVectorFileName = "C:\\Prathamesh\\Astrocytes\\Cropped_Experiment\\feature_vector.txt";
	std::string IDImageFileName = "C:\\Prathamesh\\Astrocytes\\Cropped_Experiment\\IDImage.tif";
	std::string rootPointsFileName = "C:\\Prathamesh\\Astrocytes\\Cropped_Experiment\\feature_vector_with_classes_backup.txt";
	std::string nucleiFeaturesFileName = "C:\\Prathamesh\\Astrocytes\\Cropped_Experiment\\nucleus_features_intrinsic_1.txt";
	std::string nucleiFeaturesAppendedFileName = "C:\\Prathamesh\\Astrocytes\\Cropped_Experiment\\nucleus_features_appended_1.txt";
	
	if (argc < 5)
	{
		std::cout << "AstroTracer.exe <InputFileName> <SomaCentroids.txt> <CostThreshold (Usually 1000)> <SomaImageFile>" << std::endl;
		return -1;
	}

	clock_t start_time = clock();
	std::string InputFilename = std::string(argv[1]);
	std::string SWCFilename = InputFilename;
	SWCFilename.erase(SWCFilename.length()-4,SWCFilename.length());
	SWCFilename.append("_ANT.swc");

	AstroTracer * AT = new AstroTracer();
	clock_t LoadCurvImage_start_time = clock();
	AT->LoadCurvImage(std::string(argv[1]), 0);////
	std::cout << "LoadCurvImage took: " << (clock() - LoadCurvImage_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	clock_t ReadStartPoints_start_time = clock();
	AT->ReadStartPoints(std::string(argv[2]), 1);
	std::cout << "ReadStartPoints took: " << (clock() - ReadStartPoints_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	clock_t SetCostThreshold_start_time = clock();
	AT->SetCostThreshold(atof(argv[3]));
	std::cout << "SetCostThreshold took: " << (clock() - SetCostThreshold_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	clock_t LoadSomaImage_start_time = clock();
	AT->LoadSomaImage(std::string(argv[4]));
	std::cout << "LoadSomaImage took: " << (clock() - LoadSomaImage_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	if(justComputeRootFeatures){
		
	//	clock_t RunTracing_start_time = clock();

		AT->SetScaleRange(2, 5); //(4, 5)

		//MNT->RunTracing();
		AT->CallFeatureMainExternal();

		//std::cout << "RunTracing took: " << (clock() - RunTracing_start_time)/(float) CLOCKS_PER_SEC << std::endl;

		AT->ComputeAstroFeatures(featureVectorFileName, IDImageFileName, 0, std::string("LOG"));
		//MNT->WriteSWCFile(SWCFilename, 1);

	}
	
	if(startWithCandidateRoots){

		AT->ReadRootPointsExternal(rootPointsFileName);
		AT->ReadNucleiFeaturesExternal(nucleiFeaturesFileName);
		AT->ComputeFeaturesFromCandidateRoots();
		AT->WriteNucleiFeatures(nucleiFeaturesAppendedFileName);
	}
		

	//MNT->WriteMultipleSWCFiles(SWCFilename, 1);

	//std::cout << "Total time to tracing is : " << (clock() - start_time)/(float) CLOCKS_PER_SEC << std::endl;

	return 0;

}