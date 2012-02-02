#include "MicrogliaRegionTracer.h"

int main(int argc, char* argv[])
{
	MicrogliaRegionTracer *MRT = new MicrogliaRegionTracer("E:/Farsight_Images/MicrogliaRegionTracer/GFP/joint_transforms.xml", "E:/Farsight_Images/MicrogliaRegionTracer/GFP/", "8bitkt06045_w311GFPdsu.TIF");
	
	/*std::cout << "Loading image" << std::endl;
	MRT->LoadImage("E:/farsight_images/MicrogliaRegionTracer/input.tif");*/
	
	std::cout << "Entering LoadSeedPoints" << std::endl;
	//MRT->LoadSeedPoints("E:/farsight_images/MicrogliaRegionTracer/GFP/singleseedpoint.txt");
	MRT->LoadSeedPoints("E:/farsight_images/MicrogliaRegionTracer/GFP/seedpoints.txt");

	//std::cout << "Entering LoadSeedImages" << std::endl;
	//MRT->WriteSeedImages();

	//std::cout << "Entering Trace" << std::endl;
	MRT->Trace();
}