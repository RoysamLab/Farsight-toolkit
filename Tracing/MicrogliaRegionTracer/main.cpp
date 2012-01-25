#include "MicrogliaRegionTracer.h"

int main(int argc, char* argv[])
{
	MicrogliaRegionTracer *MRT = new MicrogliaRegionTracer();
	
	std::cout << "Loading image" << std::endl;
	MRT->LoadImage("D:/farsight_images/MicrogliaRegionTracer/input.tif");
	
	std::cout << "Loading seed points" << std::endl;
	MRT->LoadSeedPoints("D:/farsight_images/MicrogliaRegionTracer/DAPI/8bitkt06041_w410DAPIdsu_seedPoints.txt");
}