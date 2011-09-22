
#include <iostream>
#include <conio.h>

#include "ftkVesselTracer.h"

#include "Common.h"



/* To do:
 1. Add visualization - MIP, 3D
 2. Ask user if preprocessing is required for a given data
 3. Make parameters available in function calls
 4. Compute Oribin in the preprocessing
 5. Side-by-side visualization of preprocessing result
 */

void main(int argc, char* argv[]){

	
	// IO for now..
	std::string input_data_path = "C:\\Lab\\data\\Amit\\25_6_1099_AJ_left_b_Ch1.tif";
	
	bool preprocess = false;
	
	ftkVesselTracer *Tracer = new ftkVesselTracer(input_data_path, preprocess);
	
	_getch();
}