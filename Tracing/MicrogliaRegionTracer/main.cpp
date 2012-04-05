#include "MicrogliaRegionTracer.h"

#include "time.h"

#include "itkMultiThreader.h"

int main(int argc, char* argv[])
{
	MicrogliaRegionTracer *MRT;
	
	if (argc == 1)	//Only for development purposes
		MRT = new MicrogliaRegionTracer("E:/Farsight_Images/MicrogliaRegionTracer/GFP/joint_transforms.xml", "E:/Farsight_Images/MicrogliaRegionTracer/GFP/", "8bitkt06045_w311GFPdsu.TIF", "E:/Farsight_Images/MicrogliaRegionTracer/DAPI/montage_8bitkt06045_w410DAPIdsu_soma.mhd");
	else if( argc < 5 )
	{
		std::cerr << "Usage: "
			<< "<joint_transforms.xml> "
			<< "<image series path> "
			<< "<filename from image series> "
			<< "<seedpoints file> "
			<< "[mask image]"
			<< std::endl;
		return 1;
	}
	std::string mask_image;
	if( argc > 5 )
	{
		mask_image = argv[5];
		MRT = new MicrogliaRegionTracer(argv[1], argv[2], argv[3], mask_image.c_str() );
	}

	

	clock_t start_time = clock();

	/*std::cout << "Loading image" << std::endl;
	MRT->LoadImage("E:/farsight_images/MicrogliaRegionTracer/input.tif");*/

	std::cout << "Entering LoadCellPoints" << std::endl;
	MRT->LoadCellPoints("E:/farsight_images/MicrogliaRegionTracer/GFP/singleseedpoint.txt");
	//MRT->LoadCellPoints("E:/farsight_images/MicrogliaRegionTracer/GFP/bottomleftseedpoint.txt");
	//MRT->LoadCellPoints("E:/farsight_images/MicrogliaRegionTracer/GFP/16seedpoints.txt");
	//MRT->LoadCellPoints("E:/farsight_images/MicrogliaRegionTracer/GFP/seedpoints.txt");
	//MRT->LoadCellPoints("E:/farsight_images/MicrogliaRegionTracer/GFP/1262_1723_169_seedPoint.txt");
	//MRT->LoadCellPoints("E:/farsight_images/MicrogliaRegionTracer/DAPI/SomaCentroids.txt");

	//MRT->LoadCellPoints( argv[4] );

	std::cout << "Entering Trace" << std::endl;
	try
	{
		MRT->Trace();
	}
	catch ( const std::exception & e )
	{
		std::cerr << "Error: " << e.what() << std::endl;
		delete MRT;
		return 1;
	}

	std::cout << "Total time for MicrogliaRegionTracing is: " << (clock() - start_time) / CLOCKS_PER_SEC << " seconds" << std::endl;

	delete MRT;
	return 0;
}
