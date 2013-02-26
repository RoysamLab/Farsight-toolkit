
#include <iostream>

#include "ftkVesselTracer.h"

#include "Common.h"



/* To do:
 . Test preprocessing with itk-only filters. This might be faster.
 . Use vesselness-enhanced diffusion filter for preprocessing http://www.insight-journal.org/browse/publication/163
 . Add vesselness to the likelihood or to the optimization http://www.insight-journal.org/browse/publication/175 (DONE - not much improvement)
 . Switch back to the old priority_queue structure to compare performance (DONE - makes no difference)
 . Compare the code performance with Amit's (DONE)
 . Implement pruning on the MSF and loop completion
 . Get a binary mask for the final segmentation (DONE)
 . Retrace for filling tracing gaps
 . Add coverage idea to the cost function  
 . Add smart editing capabilities (game theoretic/tensor voting)
 . Write SWC file output (DONE)
 . Compute vessel network features (Read papers from Audrey. To be computed after final editing. Get branch-based features from TraceEditor.)
 . Classify networks
 . Set inside region for vessel tracing and improve vessel mask by using spheres
 . Write label segmentation mask image as output
 */

int main(int argc, char* argv[]){

	if(argc < 2 || argc > 7){
		std::cout << "ftkVesselTracer.exe <Mode_no: 1. Tracing mode 2. Vessel features for nuclei> ";
		std::cout << " <(1)InputFileName, (2)NeucleiFeatureTable> <(1)preProcessData?, (2) SkeletonImage> ";
		std::cout << " <(1)useVesselness?, (2)NodeFeaturesFile> <(2)NucleiLabelImage> <(2)VesselMaskImage>" << std::endl;
		return -1;
	}

	int step_no = atoi(argv[1]);
	
	if(step_no < 1 || step_no > 2){
		std::cout << "Incorrect parameter: step_no. Returning. " << std::endl;
		return -1;
	}

	if(step_no == 1){

		std::string input_data_path = std::string(argv[2]);
		
		std::string preprocess_str = argv[3];

		std::string vesselness_str = argv[4];

		/*if(strcmp(preprocess_str.c_str(), "0") != 0 || strcmp(preprocess_str.c_str(), "1") != 0){
			std::cout << "Incorrect option for preprocessData, should be 0 or 1. " << std::endl;
			return -1;
		}*/

		bool preprocess = false;
		if(strcmp(preprocess_str.c_str(), "1") == 0)
			preprocess = true; 

		int useVesselness = 1;
		if(strcmp(vesselness_str.c_str(), "0") == 0)
			useVesselness = 0;
		else if(strcmp(vesselness_str.c_str(), "1") == 0)
			useVesselness = 1;
		else if(strcmp(vesselness_str.c_str(), "2") == 0)
			useVesselness = 2;
		
		std::cout << "Use vesselness: " << useVesselness << std::endl;

		bool startWithMST = false; 
		ftkVesselTracer *Tracer = new ftkVesselTracer(input_data_path, preprocess, startWithMST, useVesselness);	
	}

	if(step_no == 2){

		std::string nuc_table_path = std::string(argv[2]);
		std::string skeleton_img_path = std::string(argv[3]);
		std::string node_prop_path = std::string(argv[4]);
		std::string nuc_label_img_path = std::string(argv[5]);
		std::string vessel_mask_path = std::string(argv[6]);

		VesselBasedNucleiFeatures *VesselNucfeatures = new VesselBasedNucleiFeatures(nuc_table_path, skeleton_img_path, node_prop_path, nuc_label_img_path, vessel_mask_path);
	}

	return 0;
}
