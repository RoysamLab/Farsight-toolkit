
#include <iostream>

#include "ftkVesselTracer.h"

#include "Common.h"



/* To do:
 . Test preprocessing with itk-only filters. This might be faster.
 . Use vesselness-enhanced diffusion filter for preprocessing http://www.insight-journal.org/browse/publication/163
 . Add vesselness to the likelihood or to the optimization http://www.insight-journal.org/browse/publication/175
 . Switch back to the old priority_queue structure to compare performance (DONE - makes no difference)
 . Compare the code performance with Amit's 
 . Implement pruning on the MSF and loop completion
 . Get a binary mask for the final segmentation
 . Retrace for filling tracing gaps
 . Add coverage idea to the cost function
 . Add smart editing capabilities (game theoretic/tensor voting)
 . Write SWC file output (DONE)
 . Compute vessel network features (new class)
 . Classify networks
 */

int main(int argc, char* argv[]){

	if(argc < 2 || argc > 3){
		std::cout << "ftkVesselTracer.exe <InputFileName> <preProcessData?>" << std::endl;
		return -1;
	}


	// IO for now..
	//std::string input_data_path = "C:\\Lab\\data\\Amit\\25_6_1099_AJ_left_b_Ch1_cropped.tif";
	//std::string input_data_path = "C:\\Prathamesh\\Vessels\\CroppedExp\\25_6_1099_AJ_left_b_Ch1_cropped_8bit.mhd";
	//std::string input_data_path = "C:\\Lab\\data\\Naren_data\\9.5_1\\ContourImage__TiffStack_pre_1.tif";

	//std::string input_data_path = "C:\\Prathamesh\\Vessels\\EriksenExp\\_CROP_Lectin_IM2011-12-21_15-18_pre.tif";
	//std::string input_data_path = "C:\\Prathamesh\\Vessels\\BigExp\\25_0_1101_AJ_left_Ch1_8bit.tif";
	//std::string input_data_path = "C:\\Lab\\data\\forAstrocytes\\ASTRO_Cropped2_montage_8bitkt11410_w212TRITCdsu.tif";

	std::string input_data_path = std::string(argv[1]);

	
	std::string preprocess_str = argv[2];

	/*if(strcmp(preprocess_str.c_str(), "0") != 0 || strcmp(preprocess_str.c_str(), "1") != 0){
		std::cout << "Incorrect option for preprocessData, should be 0 or 1. " << std::endl;
		return -1;
	}*/

	bool preprocess = false;
	if(strcmp(preprocess_str.c_str(), "1") == 0)
		preprocess = true; 

	bool useVesselness = true;

	bool startWithMST = false; 
	ftkVesselTracer *Tracer = new ftkVesselTracer(input_data_path, preprocess, startWithMST, useVesselness);


	// Please do not use this, it is not cross platform.
	//_getch();

	return 0;
}
