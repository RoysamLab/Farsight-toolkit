#include "MicrogliaRegionTracer.h"

#include "time.h"

int main(int argc, char* argv[])
{
	MicrogliaRegionTracer *MRT = new MicrogliaRegionTracer("E:/Farsight_Images/MicrogliaRegionTracer/GFP/joint_transforms.xml", "E:/Farsight_Images/MicrogliaRegionTracer/GFP/", "8bitkt06045_w311GFPdsu.TIF");
	
	clock_t start_time = clock();

	/*std::cout << "Loading image" << std::endl;
	MRT->LoadImage("E:/farsight_images/MicrogliaRegionTracer/input.tif");*/
	
	std::cout << "Entering LoadCellPoints" << std::endl;
	//MRT->LoadCellPoints("E:/farsight_images/MicrogliaRegionTracer/GFP/singleseedpoint.txt", "E:/Farsight_Images/MicrogliaRegionTracer/GFP/montage_8bitkt06041_w311GFPdsu_soma.mhd" );
	//MRT->LoadCellPoints("E:/farsight_images/MicrogliaRegionTracer/GFP/bottomleftseedpoint.txt", "E:/Farsight_Images/MicrogliaRegionTracer/GFP/montage_8bitkt06041_w311GFPdsu_soma.mhd");
	//MRT->LoadCellPoints("E:/farsight_images/MicrogliaRegionTracer/GFP/16seedpoints.txt", "E:/Farsight_Images/MicrogliaRegionTracer/GFP/montage_8bitkt06041_w311GFPdsu_soma.mhd");
	MRT->LoadCellPoints("E:/farsight_images/MicrogliaRegionTracer/GFP/seedpoints.txt", "E:/Farsight_Images/MicrogliaRegionTracer/GFP/montage_8bitkt06041_w311GFPdsu_soma.mhd");
	//MRT->LoadCellPoints("E:/farsight_images/MicrogliaRegionTracer/GFP/850_444_36_seedPoint.txt", "E:/Farsight_Images/MicrogliaRegionTracer/GFP/montage_8bitkt06041_w311GFPdsu_soma.mhd");
	
	std::cout << "Entering Trace" << std::endl;
	MRT->Trace();
	//MRT->Trace2();

	std::cout << "Total time for MicrogliaRegionTracing is: " << (clock() - start_time) / CLOCKS_PER_SEC << " seconds" << std::endl;

	//getchar();
}
