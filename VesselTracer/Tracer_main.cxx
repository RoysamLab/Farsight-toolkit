
#include <iostream>

#include "ftkVesselTracer.h"

#include "Common.h"



/* To do:
 . Ask user if preprocessing is required for a given data
 . Make parameters available in function calls from  main function
 . Side-by-side visualization of preprocessing result
 */

int main(int argc, char* argv[]){


	// IO for now..
	std::string input_data_path = "C:\\Lab\\data\\Amit\\25_6_1099_AJ_left_b_Ch1_cropped.tif";
	//std::string input_data_path = "C:\\Lab\\data\\forAstrocytes\\ASTRO_Cropped2_montage_8bitkt11410_w212TRITCdsu.tif";

	bool preprocess = false; //false;
	bool startWithMST = false; //true; //false;

	ftkVesselTracer *Tracer = new ftkVesselTracer(input_data_path, preprocess, startWithMST);

	// Please do not use this, it is not cross platform.
	//_getch();

	return 0;
}
