
#include <iostream>
#include <conio.h>

#include "ftkVesselTracer.h"

#include "Common.h"



/* To do:
 . Ask user if preprocessing is required for a given data
 . Make parameters available in function calls from  main function
 . Side-by-side visualization of preprocessing result
 */

void main(int argc, char* argv[]){

	
	// IO for now..
	std::string input_data_path = "C:\\Lab\\data\\Amit\\25_6_1099_AJ_left_b_Ch1_cropped.tif";
	
	bool preprocess = false;
	bool startWithMST = true; //false;
	
	ftkVesselTracer *Tracer = new ftkVesselTracer(input_data_path, preprocess, startWithMST);
	
	_getch();
}