#include "MicrogliaRegionTracer.h"
#include "time.h"
#include "itkMultiThreader.h"
#include <cstdlib>

int main(int argc, char* argv[])
{
	MicrogliaRegionTracer MRT;
	
	if (argc == 1)	//Only for development purposes
	{		
		MRT.SetJointTransformsFile("E:/Farsight_Images/MicrogliaRegionTracer/GFP/joint_transforms.xml");
		MRT.SetImageSeriesPath("E:/Farsight_Images/MicrogliaRegionTracer/GFP/");
		MRT.SetAnchorImage("8bitkt06045_w311GFPdsu.TIF");
		MRT.SetSomaImage("E:/Farsight_Images/MicrogliaRegionTracer/DAPI/montage_8bitkt06045_w410DAPIdsu_soma.mhd");
				
		try
		{
			//MRT.LoadSeedPoints("E:/farsight_images/MicrogliaRegionTracer/GFP/singleseedpoint.txt");
			//MRT.LoadSeedPoints("E:/farsight_images/MicrogliaRegionTracer/GFP/bottomleftseedpoint.txt");
			//MRT.LoadSeedPoints("E:/farsight_images/MicrogliaRegionTracer/GFP/16seedpoints.txt");
			MRT.LoadSeedPoints("E:/farsight_images/MicrogliaRegionTracer/DAPI/SomaCentroids.txt");
			//MRT.LoadSeedPoints("E:/farsight_images/MicrogliaRegionTracer/GFP/763_690_19_seedPoint.txt");
			//MRT.LoadSeedPoints("E:/farsight_images/MicrogliaRegionTracer/GFP/850_444_36_seedpoint.txt");
			//MRT.LoadSeedPoints("E:/farsight_images/MicrogliaRegionTracer/GFP/726_78_238_seedpoint.txt");
			//MRT.LoadSeedPoints("E:/farsight_images/MicrogliaRegionTracer/GFP/140_1559_164_seedpoint.txt");
			//MRT.LoadSeedPoints("E:/Farsight_Images/MicrogliaRegionTracer/DAPI/20seedpoint.txt");
		}
		catch (std::exception & e)
		{
			std::cerr << "Error detected in LoadSeedPoints: " << std::endl;
			std::cerr << e.what() << std::endl;
			return -1;
		}
		MRT.SetAspectRatio(3.0);
	}
	else if( argc < 6 )
	{
		std::cerr << "Usage: "
			<< "<joint_transforms.xml> "
			<< "<image series path> "
			<< "<filename from image series> "
			<< "<seedpoints file> "
			<< "<aspect_ratio> "
			<< "[mask image]"
			<< std::endl;
		return 1;
	}
    else //( argc >= 6 )
	{
		std::string mask_image;
        mask_image = argv[5];
		MRT.SetJointTransformsFile(argv[1]);
		MRT.SetImageSeriesPath(argv[2]);
		MRT.SetAnchorImage(argv[3]);
		MRT.LoadSeedPoints(argv[4]);
		MRT.SetAspectRatio(atof(argv[5]));
		MRT.SetSomaImage(argv[6]);
	}

	

	clock_t start_time = clock();
	std::cout << "Entering Trace" << std::endl;
	try
	{
		MRT.Trace();
	}
	catch ( const std::exception & e )
	{
		std::cerr << "Error: " << e.what() << std::endl;
		return 1;
	}

	std::cout << "Total time for MicrogliaRegionTracing is: " << (clock() - start_time) / CLOCKS_PER_SEC << " seconds" << std::endl;

	return 0;
}
