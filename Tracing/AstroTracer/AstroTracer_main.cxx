#include "AstroTracer.h"
#include "time.h"
#include <conio.h>

int main(int argc, char* argv[])
//int main(void)
{	
	
    //argv[1] = "C:\\Users\\msavelon\\Desktop\\Astro\\TrainingWithBill\\CroppedAstroTRITC.tif"; //"C:\\ROYSAMLAB\\FARSIGHT\\BINARY\\exe\\Debug\\R2080_6wk_Crop_sigma_0.030000_CV.tif";
	argv[1] = "C:\\Users\\msavelon\\Desktop\\Astro\\IdealRootExperiment\\kt11780_w212TRITCdsu_pre.tif";
	//argv[2] = "C:\\Users\\msavelon\\Desktop\\Astro\\TrainingWithBill\\roots_cropped.txt"; //"C:\\ROYSAMLAB\\FARSIGHT\\BINARY\\exe\\Debug\\11111111.tif";
	argv[2] = "C:\\Users\\msavelon\\Desktop\\Astro\\IdealRootExperiment\\centroids.txt";
	//argv[4] = "C:\\Users\\msavelon\\Desktop\\Astro\\TrainingWithBill\\labels_2_label_nuc_edited_somas_0.tif"; 
	argv[4] = "C:\\Users\\msavelon\\Desktop\\Astro\\IdealRootExperiment\\Black.tif"; 

	

		//CroppedAstroLABEL.tif //"C:\\ROYSAMLAB\\FARSIGHT\\BINARY\\exe\\Debug\\SegParams.ini";
	
	//argv[3] = "400"; //"C:\\ROYSAMLAB\\FARSIGHT\\BINARY\\exe\\Debug\\ParameterFile.ini";

	//argc = 5;

	bool justComputeRootFeatures = false; //false;
	bool startWithCandidateRoots = false;//true;
	bool getCentroids=false;
	bool doTracing=true;

	std::string tracePointsFileName = "C:\\Users\\msavelon\\Desktop\\Astro\\TrainingWithBill\\trace_points.txt";
	std::string featureVectorFileName = "C:\\Users\\msavelon\\Desktop\\Astro\\TrainingWithBill\\feature_vector.txt";
	std::string IDImageFileName = "C:\\Users\\msavelon\\Desktop\\Astro\\TrainingWithBill\\IDImage.tif";
	std::string rootPointsFileName = "C:\\Users\\msavelon\\Desktop\\Astro\\TrainingWithBill\\feature_vector_with_classes_backup.txt";
	std::string nucleiFeaturesFileName = "C:\\Users\\msavelon\\Desktop\\Astro\\TrainingWithBill\\nucleus_features_intrinsic_4.txt";
	std::string nucleiFeaturesAppendedFileName = "C:\\Users\\msavelon\\Desktop\\Astro\\TrainingWithBill\\nucleus_features_appended_1.txt";
	std::string finalNucleiTableFileName,finalIDImageFileName;
	std::string centroidsForTracingFileName="C:\\Users\\msavelon\\Desktop\Astro\\IdealRootExperiment\\centroids.txt";
	const char* optionsFileName="C:\\Users\\msavelon\\Desktop\\Astro\\IdealRootExperiment\\options_mnt";
	
	//if (argc < 5)
	//{
	//	std::cout << "AstroTracer.exe <InputFileName> <SomaCentroids.txt> <CostThreshold (Usually 1000)> <SomaImageFile>" << std::endl;
	//	return -1;
	//}

	AstroTracer * AT = new AstroTracer();

	clock_t LoadParameters_start_time = clock();
	AT->LoadParameters(optionsFileName);
	std::cout << "LoadParameters took: " << (clock() - LoadParameters_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	clock_t start_time = clock();
	std::string InputFilename = std::string(argv[1]);
	std::string SWCFilename = InputFilename;
	SWCFilename.erase(SWCFilename.length()-4,SWCFilename.length());
	SWCFilename.append("_ANT.swc");

	//AstroTracer * AT = new AstroTracer();
	clock_t LoadCurvImage_start_time = clock();
	AT->LoadCurvImage(std::string(argv[1]), 0);////
	std::cout << "LoadCurvImage took: " << (clock() - LoadCurvImage_start_time)/(float) CLOCKS_PER_SEC << std::endl;
	
	// USE THIS FOR READING THE REAL ROOT POINTSI IN PART3
	clock_t ReadStartPoints_start_time = clock();
	AT->ReadStartPoints(std::string(argv[2]), 0);///
	std::cout << "ReadStartPoints took: " << (clock() - ReadStartPoints_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	clock_t SetCostThreshold_start_time = clock();
	AT->SetCostThreshold(AT->cost_threshold);
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
		
	if (getCentroids) {
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

	}

	//MNT->WriteMultipleSWCFiles(SWCFilename, 1);

	//std::cout << "Total time to tracing is : " << (clock() - start_time)/(float) CLOCKS_PER_SEC << std::endl;
	
	getch();

	return 0;

}