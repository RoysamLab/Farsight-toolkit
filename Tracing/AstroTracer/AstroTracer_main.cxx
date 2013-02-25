
#include "AstroTracer.h"
#include "time.h"

int main(int argc, char* argv[]){	
    
	if(argc < 2 || argc > 8){
		std::cout << "AstroTracer.exe <InputFileName> <SomaImageFile> <DoPreprocessing?> <Step_no: 0: Optimize coverage 1:Compute root features, 2:Compute nuclei features";
		std::cout << " 3:Generate centroids for tracing and run tracing 4. Only run tracing (like MNT)> <OptionsFileName> <RootPointsFileName>";
		std::cout << " <NucleiFeaturesFileName>" << std::endl;
		return -1;
	}

	int do_preprocessing = atoi(argv[3]);
	
	int step_no = atoi(argv[4]);
	
	if(step_no < 0 || step_no > 4){
		std::cout << "Incorrect parameter: step_no. " << std::endl;
		return -1;
	}
	
	bool getCentroids = false;
	bool doTracing = false;

	std::string InputFilename = std::string(argv[1]);

	std::string coverageFileName = InputFilename;
	coverageFileName.erase(coverageFileName.length()-4, coverageFileName.length());
	coverageFileName.append("_coverage.txt");

	std::string nucleiFeaturesAppendedFileName = InputFilename;
	nucleiFeaturesAppendedFileName.erase(nucleiFeaturesAppendedFileName.length()-4, nucleiFeaturesAppendedFileName.length());
	nucleiFeaturesAppendedFileName.append("_nuc_features_all.txt");

	std::string centroidsForTracingFileName = InputFilename;
	centroidsForTracingFileName.erase(centroidsForTracingFileName.length()-4, centroidsForTracingFileName.length());
	centroidsForTracingFileName.append("_centroids.txt");

	std::string finalIDImageFileName = InputFilename;
	finalIDImageFileName.erase(finalIDImageFileName.length()-4, finalIDImageFileName.length());
	finalIDImageFileName.append("_finalTracingCentroids.tif");

	std::string SWCFilename = InputFilename;
	SWCFilename.erase(SWCFilename.length()-4, SWCFilename.length());
	SWCFilename.append("_AstroTraces.swc");
	
	AstroTracer * AT = new AstroTracer();

	AT->SetInputDataPath(InputFilename.erase(InputFilename.length()-4, InputFilename.length()));

	clock_t LoadCurvImage_start_time = clock();
	AT->LoadCurvImageFromPath(std::string(argv[1]), 0);
	std::cout << "LoadCurvImage took: " << (clock() - LoadCurvImage_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	clock_t LoadSomaImage_start_time = clock();
	AT->LoadSomaImage(std::string(argv[2]));
	std::cout << "LoadSomaImage took: " << (clock() - LoadSomaImage_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	if(do_preprocessing){

		clock_t preprocessing_start_time = clock();
		AT->DoPreprocessing();
		std::cout << "Anisotropic diffusion and GVF took: " << (clock() - preprocessing_start_time)/(float) CLOCKS_PER_SEC << std::endl;	
	}
	
	std::cout << "No preprocessing. Loading from disk. " << std::endl;
	AT->LoadPreprocessingResults();
	

	// step 0 is for testing whatever you want
	if(step_no == 0){
		AT->OptimizeCoverage(coverageFileName, true);	
		std::cout << std::endl << "Done with step 0. " << std::endl;
	}

	if(step_no == 1){

		//bool inPipeline = false;
		AT->LoadParameters(argv[5]);	
		AT->SetScaleRange(4, 4); //(2, 5); //(2, 2)
		AT->CallFeatureMainExternal();
		AT->ComputeRootPointFeatures();
		std::cout << "Done with computing root-based features. " << std::endl;
	}
	if(step_no == 2){

		AT->ReadRootPointsExternal(std::string(argv[6]));
		AT->ReadNucleiFeaturesExternal(std::string(argv[7]));
		AT->ComputeFeaturesFromCandidateRoots();
		AT->WriteNucleiFeatures(nucleiFeaturesAppendedFileName);
		std::cout << "Done with computing nuclei-based features. " << std::endl;
	}
	if(step_no == 3){

		//For this part, make sure:
		// 1. The soma file should have only astrocyte somas	2. The nucleiFeaturesFile and the rootPointsFile should have the respective class labels
		
		AT->LoadParameters(argv[5]);	
		AT->SetCostThreshold(AT->cost_threshold);

		AT->ReadRootPointsExternal(std::string(argv[6]));
		AT->ReadFinalNucleiTable(std::string(argv[7])); 
		
		AT->GetCentroidsForTracing(centroidsForTracingFileName, finalIDImageFileName);
		
		std::cout << "GetCentroidsForTracing finished. " << std::endl; 

		if(true){
			
			AT->ReadStartPointsInternal();

			clock_t RunTracing_start_time = clock();
			AT->RunTracing();
			std::cout << "RunTracing took: " << (clock() - RunTracing_start_time)/(float) CLOCKS_PER_SEC << std::endl;

			AT->WriteSWCFile(SWCFilename, 1);
	
		}
	}
	if(step_no == 4){

		//For this part, make sure:
		// 1. Root points file name is a simple list of centroids

		AT->LoadParameters(argv[5]);	
		AT->SetCostThreshold(AT->cost_threshold);
		
		AT->ReadStartPointsFromPath(std::string(argv[6]), 0);
		
		clock_t RunTracing_start_time = clock();
		AT->RunTracing();
		std::cout << "RunTracing took: " << (clock() - RunTracing_start_time)/(float) CLOCKS_PER_SEC << std::endl;

		AT->WriteSWCFile(SWCFilename, 1);
	}
		
	return 0;
}