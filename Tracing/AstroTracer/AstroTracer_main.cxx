#include "AstroTracer.h"
#include "time.h"
#include <conio.h>

int main(int argc, char* argv[]){	
    
	/*argv[1] = "C:\\Prathamesh\\Astrocytes\\ControlExp\\Testing\\_CROPPED_montage_kt01341_w212TRITCdsu-2_pre.tif";
	argv[6] = "C:\\Prathamesh\\Astrocytes\\TestPipelineExp1\\idealRoots.txt";
	
	argv[2] = "C:\\Prathamesh\\Astrocytes\\ControlExp\\Testing\\nuc_labels.tif";*/


	if(argc < 2 || argc > 7)
	{
		std::cout << "AstroTracer.exe <InputFileName> <SomaImageFile> <Step_no: 0: Optimize coverage 1:Compute root features, 2:Compute nuclei features";
		std::cout << " 3:Generate centroids for tracing and run tracing 4. Only run tracing (like MNT)> <OptionsFileName> <RootPointsFileName> <NucleiFeaturesFileName>" << std::endl;
		return -1;
	}
	
	int step_no = atoi(argv[3]);
	
	if(step_no < 0 || step_no > 4){
		std::cout << "Incorrect parameter: step_no. " << std::endl;
		return -1;
	}
	
	bool getCentroids = false;
	bool doTracing = false;

	std::string InputFilename = std::string(argv[1]);
	
	std::string featureVectorFileName = InputFilename;
	featureVectorFileName.erase(featureVectorFileName.length()-4, featureVectorFileName.length());
	featureVectorFileName.append("_feature_vector_roots.txt");
	
	std::string IDImageFileName = InputFilename;
	IDImageFileName.erase(IDImageFileName.length()-4, IDImageFileName.length());
	IDImageFileName.append("_RootsImage.tif");

	std::string nucleiFeaturesAppendedFileName = InputFilename;
	nucleiFeaturesAppendedFileName.erase(nucleiFeaturesAppendedFileName.length()-4, nucleiFeaturesAppendedFileName.length());
	nucleiFeaturesAppendedFileName.append("_nuc_features_all.txt");

	std::string coverageFileName = InputFilename;
	coverageFileName.erase(coverageFileName.length()-4, coverageFileName.length());
	coverageFileName.append("_coverage.txt");

	std::string centroidsForTracingFileName = InputFilename;
	centroidsForTracingFileName.erase(centroidsForTracingFileName.length()-4, centroidsForTracingFileName.length());
	centroidsForTracingFileName.append("_centroids.txt");

	std::string finalIDImageFileName = InputFilename;
	finalIDImageFileName.erase(finalIDImageFileName.length()-4, finalIDImageFileName.length());
	finalIDImageFileName.append("_finalTracingCentroids.tif");

	std::string SWCFilename = InputFilename;
	SWCFilename.erase(SWCFilename.length()-4, SWCFilename.length());
	SWCFilename.append("_traces.swc");
	
	
	//std::string featureVectorFileName = "C:\\Prathamesh\\Astrocytes\\ControlExp\\Training\\feature_vector_roots.txt";
	//std::string IDImageFileName = "C:\\Prathamesh\\Astrocytes\\ControlExp\\Testing\\RootsImage.tif";
	
	//std::string rootPointsFileName = "C:\\Prathamesh\\Astrocytes\\ControlExp\\Testing\\feature_vector_roots_classified.txt"; //argv[4]
	//std::string nucleiFeaturesFileName = "C:\\Prathamesh\\Astrocytes\\ControlExp\\Testing\\nuc_intrinsic_asso_table.txt"; //argv[5]
	//std::string nucleiFeaturesAppendedFileName = "C:\\Prathamesh\\Astrocytes\\ControlExp\\Testing\\nuc_features_all.txt";
	
	//std::string finalNucleiTableFileName, finalIDImageFileName;
	//std::string centroidsForTracingFileName = "C:\\Users\\msavelon\\Desktop\Astro\\IdealRootExperiment\\centroids.txt";
	
	//const char* optionsFileName = "C:\\Prathamesh\\Astrocytes\\TestPipelineExp1\\data\\options_mnt";
	

	AstroTracer * AT = new AstroTracer();

	
	clock_t LoadCurvImage_start_time = clock();
	AT->LoadCurvImage(std::string(argv[1]), 0);////
	std::cout << "LoadCurvImage took: " << (clock() - LoadCurvImage_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	clock_t LoadSomaImage_start_time = clock();
	AT->LoadSomaImage(std::string(argv[2]));
	std::cout << "LoadSomaImage took: " << (clock() - LoadSomaImage_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	

	// step 0 is for testing whatever you want
	if(step_no == 0){

		//ObjectnessMeasures obj_measures;
		//AT->ComputeObjectnessImage(obj_measures);
		AT->OptimizeCoverage(coverageFileName);	

		std::cout << std::endl << "Done with step 0. " << std::endl;
	}

	if(step_no == 1){

		AT->LoadParameters(argv[4]);	
		AT->SetScaleRange(4, 4); //(2, 5); //(2, 2)
		AT->CallFeatureMainExternal();
		AT->ComputeAstroFeatures(featureVectorFileName, IDImageFileName, 0, std::string("LOG"));
		std::cout << "Done with computing root-based features. " << std::endl;
	}
	if(step_no == 2){

		AT->ReadRootPointsExternal(std::string(argv[5]));
		AT->ReadNucleiFeaturesExternal(std::string(argv[6]));
		AT->ComputeFeaturesFromCandidateRoots();
		AT->WriteNucleiFeatures(nucleiFeaturesAppendedFileName);
		std::cout << "Done with computing nuclei-based features. " << std::endl;
	}
	if(step_no == 3){

		//For this part, make sure:
		// 1. The soma file should have only astrocyte somas	2. The nucleiFeaturesFile and the rootPointsFile should have the respective class labels
		
		AT->LoadParameters(argv[4]);	
		AT->SetCostThreshold(AT->cost_threshold);

		AT->ReadRootPointsExternal(std::string(argv[5]));
		AT->ReadFinalNucleiTable(std::string(argv[6])); 
		
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

		AT->LoadParameters(argv[4]);	
		AT->SetCostThreshold(AT->cost_threshold);
		
		AT->ReadStartPoints(std::string(argv[5]), 0);
		
		clock_t RunTracing_start_time = clock();
		AT->RunTracing();
		std::cout << "RunTracing took: " << (clock() - RunTracing_start_time)/(float) CLOCKS_PER_SEC << std::endl;

		AT->WriteSWCFile(SWCFilename, 1);
	}
		
	return 0;
}