#include "MicrogliaRegionTracer.h"

int main(int argc, char* argv[])
{
	MicrogliaRegionTracer *MRT = new MicrogliaRegionTracer();
	
	/*std::cout << "Loading image" << std::endl;
	MRT->LoadImage("E:/farsight_images/MicrogliaRegionTracer/input.tif");*/
	
	std::cout << "Loading seed points" << std::endl;
	MRT->LoadSeedPoints("E:/farsight_images/MicrogliaRegionTracer/DAPI/seedpoints.txt");

	std::cout << "Entering WriteInitialMicrogliaImages" << std::endl;
	MRT->WriteInitialMicrogliaImages();
}