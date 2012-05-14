#include "AstroTracer.h"
#include "time.h"
#include <conio.h>

int main(int argc, char* argv[])
//int main(void)
{	
	
    argv[1] = "C:\\Prathamesh\\Astrocytes\\TestPipelineExp1\\data\\kt11780_w212TRITCdsu_pre.tif"; 
	//argv[1] = "C:\\Prathamesh\\Astrocytes\\Cropped_Experiment\\CroppedAstroTRITC.tif"; 

	argv[2] = "C:\\Prathamesh\\Astrocytes\\TestPipelineExp1\\roots_cropped.txt"; 
	
	//argv[4] = "C:\\Prathamesh\\Astrocytes\\TestPipelineExp1\\data\\label_nuc_1.tif"; 			
	argv[4] = "C:\\Prathamesh\\Astrocytes\\TestPipelineExp1\\data\\class_1_labels_somas_0.tif"; 

	//argv[4] = "C:\\Prathamesh\\Astrocytes\\Cropped_Experiment\\FourClassExp\\class_1_soma_somas_0.tif";

	argv[3] = "400";

	argc = 5;

	bool justComputeRootFeatures = false; //false;
	bool startWithCandidateRoots = false;//true;
	bool doTracing = true;

	// Output of part 1
	std::string featureVectorFileName = "C:\\Prathamesh\\Astrocytes\\TestPipelineExp1\\feature_vector.txt";
	std::string IDImageFileName = "C:\\Prathamesh\\Astrocytes\\TestPipelineExp1\\rootsImage.tif";


	std::string rootPointsFileName = "C:\\Prathamesh\\Astrocytes\\TestPipelineExp1\\Part1_result.txt";
	//std::string rootPointsFileName = "C:\\Prathamesh\\Astrocytes\\Cropped_Experiment\\feature_vector_with_classes_backup.txt";

	std::string nucleiFeaturesFileName = "C:\\Prathamesh\\Astrocytes\\TestPipelineExp1\\data\\table_nuc_1.txt";
	//std::string nucleiFeaturesFileName = "C:\\Prathamesh\\Astrocytes\\Cropped_Experiment\\nucleus_features_intrinsic_4.txt";

	std::string nucleiFeaturesAppendedFileName = "C:\\Prathamesh\\Astrocytes\\TestPipelineExp1\\data\\table_nuc_1_appended.txt";

	std::string finalNucleiTableFileName = "C:\\Prathamesh\\Astrocytes\\TestPipelineExp1\\data\\table_final.txt";
	std::string finalIDImageFileName = "C:\\Prathamesh\\Astrocytes\\TestPipelineExp1\\finalRootsImage.tif";
	std::string centroidsForTracingFileName="C:\\Prathamesh\\Astrocytes\\TestPipelineExp1\\data\\centroids.txt";
	
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
	
	// USE THIS FOR READING THE REAL ROOT POINTSI IN PART3
	/*clock_t ReadStartPoints_start_time = clock();
	AT->ReadStartPoints(std::string(argv[2]), 1);
	std::cout << "ReadStartPoints took: " << (clock() - ReadStartPoints_start_time)/(float) CLOCKS_PER_SEC << std::endl;*/

	clock_t SetCostThreshold_start_time = clock();
	AT->SetCostThreshold(atof(argv[3]));
	std::cout << "SetCostThreshold took: " << (clock() - SetCostThreshold_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	clock_t LoadSomaImage_start_time = clock();
	AT->LoadSomaImage(std::string(argv[4]));
	std::cout << "LoadSomaImage took: " << (clock() - LoadSomaImage_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	if(justComputeRootFeatures){
		
	//	clock_t RunTracing_start_time = clock();

		AT->SetScaleRange(2, 2); //(2, 5)

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
		

	if (doTracing) {
		
		AT->ReadRootPointsExternal(rootPointsFileName);
		AT->ReadFinalNucleiTable(finalNucleiTableFileName);

		AT->GetCentroidsForTracing(centroidsForTracingFileName, finalIDImageFileName);
		
		std::cout << "GetCentroidsForTracing finished!" << std::endl; 

		//clock_t ReadStartPoints_start_time = clock();
		//AT->ReadStartPoints(std::string("C:\\Users\\msavelon\\Desktop\\Astro\\TrainingWithBill\\centroids.txt"), 1);
		//std::cout << "ReadStartPoints took: " << (clock() - ReadStartPoints_start_time)/(float) CLOCKS_PER_SEC << std::endl;


		//clock_t RunTracing_start_time = clock();
		//AT->RunTracing();
		//std::cout << "RunTracing took: " << (clock() - RunTracing_start_time)/(float) CLOCKS_PER_SEC << std::endl;

		//AT->WriteSWCFile(SWCFilename, 1);

	}

	//MNT->WriteMultipleSWCFiles(SWCFilename, 1);

	//std::cout << "Total time to tracing is : " << (clock() - start_time)/(float) CLOCKS_PER_SEC << std::endl;
	
	getch();

	return 0;

}