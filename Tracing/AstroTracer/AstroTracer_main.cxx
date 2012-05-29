#include "AstroTracer.h"
#include "time.h"
#include <conio.h>

int main(int argc, char* argv[]){	
    
	/*argv[1] = "C:\\Prathamesh\\Astrocytes\\ControlExp\\Testing\\_CROPPED_montage_kt01341_w212TRITCdsu-2_pre.tif";
	argv[6] = "C:\\Prathamesh\\Astrocytes\\TestPipelineExp1\\idealRoots.txt";
	
	argv[2] = "C:\\Prathamesh\\Astrocytes\\ControlExp\\Testing\\nuc_labels.tif";*/


	if(argc < 2 || argc > 6)
	{
		std::cout << "AstroTracer.exe <InputFileName> <SomaImageFile> <Step_no: 1:Compute root features, 2:Compute nuclei features> <RootPointsFileName> <NucleiFeaturesFileName>" << std::endl;
		return -1;
	}
	
	int step_no = atoi(argv[3]);
	
	if(step_no < 1 || step_no > 2){
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
	
	
	//std::string featureVectorFileName = "C:\\Prathamesh\\Astrocytes\\ControlExp\\Training\\feature_vector_roots.txt";
	//std::string IDImageFileName = "C:\\Prathamesh\\Astrocytes\\ControlExp\\Testing\\RootsImage.tif";
	
	//std::string rootPointsFileName = "C:\\Prathamesh\\Astrocytes\\ControlExp\\Testing\\feature_vector_roots_classified.txt"; //argv[4]
	//std::string nucleiFeaturesFileName = "C:\\Prathamesh\\Astrocytes\\ControlExp\\Testing\\nuc_intrinsic_asso_table.txt"; //argv[5]
	//std::string nucleiFeaturesAppendedFileName = "C:\\Prathamesh\\Astrocytes\\ControlExp\\Testing\\nuc_features_all.txt";
	
	//std::string finalNucleiTableFileName, finalIDImageFileName;
	//std::string centroidsForTracingFileName = "C:\\Users\\msavelon\\Desktop\Astro\\IdealRootExperiment\\centroids.txt";
	
	//const char* optionsFileName = "C:\\Prathamesh\\Astrocytes\\TestPipelineExp1\\data\\options_mnt";
	

	AstroTracer * AT = new AstroTracer();

	/*clock_t LoadParameters_start_time = clock();
	AT->LoadParameters(optionsFileName);
	std::cout << "LoadParameters took: " << (clock() - LoadParameters_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	clock_t start_time = clock();
	std::string InputFilename = std::string(argv[1]);
	std::string SWCFilename = InputFilename;
	SWCFilename.erase(SWCFilename.length()-4,SWCFilename.length());
	SWCFilename.append("_ANT.swc");*/

	clock_t LoadCurvImage_start_time = clock();
	AT->LoadCurvImage(std::string(argv[1]), 0);////
	std::cout << "LoadCurvImage took: " << (clock() - LoadCurvImage_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	clock_t LoadSomaImage_start_time = clock();
	AT->LoadSomaImage(std::string(argv[2]));
	std::cout << "LoadSomaImage took: " << (clock() - LoadSomaImage_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	
	if(step_no == 1){
		
		AT->SetScaleRange(4, 4); //(2, 5); //(2, 2)
		AT->CallFeatureMainExternal();
		AT->ComputeAstroFeatures(featureVectorFileName, IDImageFileName, 0, std::string("LOG"));
		std::cout << "Done with computing root-based features. " << std::endl;
	}
	if(step_no == 2){

		AT->ReadRootPointsExternal(std::string(argv[4]));
		AT->ReadNucleiFeaturesExternal(std::string(argv[5]));
		AT->ComputeFeaturesFromCandidateRoots();
		AT->WriteNucleiFeatures(nucleiFeaturesAppendedFileName);
		std::cout << "Done with computing nuclei-based features. " << std::endl;
	}
		

	// This part is for tracing
	/*clock_t ReadStartPoints_start_time = clock();
	AT->ReadStartPoints(std::string(argv[2]), 0);///
	std::cout << "ReadStartPoints took: " << (clock() - ReadStartPoints_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	clock_t SetCostThreshold_start_time = clock();
	AT->SetCostThreshold(AT->cost_threshold);
	std::cout << "SetCostThreshold took: " << (clock() - SetCostThreshold_start_time)/(float) CLOCKS_PER_SEC << std::endl;
	*/

	/*if (getCentroids) {
		AT->ReadRootPointsExternal(rootPointsFileName);
		AT->ReadFinalNucleiTable(finalNucleiTableFileName); 
		AT->GetCentroidsForTracing(centroidsForTracingFileName, finalIDImageFileName);
		
		std::cout << "GetCentroidsForTracing finished!" << std::endl; 
	}

	if (doTracing) {

		clock_t RunTracing_start_time = clock();
		AT->RunTracing();
		std::cout << "RunTracing took: " << (clock() - RunTracing_start_time)/(float) CLOCKS_PER_SEC << std::endl;

		AT->WriteSWCFile(SWCFilename, 1);

	}*/


	return 0;
}