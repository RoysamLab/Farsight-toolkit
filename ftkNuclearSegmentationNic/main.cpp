// ############################################################################################################################################################################

//#includes goes here

#include "ftkNuclearSegmentationNic.h"

// FTK INCLUDES
#include <ftkObject.h>
#include <ftkImage/ftkImage.h>
#include <ftkCommon/ftkUtils.h>


// STANDARD INCLUES
#include <vector>
#include <iostream>
#include <map>

// ITK INCLUDES
#include <itkImageIOBase.h>


// ############################################################################################################################################################################

int main( int argc, char * argv[] ){

	// rm -rf exe/ftkNuclearSegmentationNic; cmake .; make -j16; ./exe/ftkNuclearSegmentationNic
	
	std::string path = std::string(argv[0]);
	int lastslash;
	path = path.substr(0, path.rfind('\\'));
	std::cout << path << std::endl;
	
	
	
	
// 	std::cout << std::endl << "	step 1";
// 	std::cout << std::endl << "	step 1";
	ftk::Image::Pointer myImg;
	//std::string inputName = "testImages\Histo_Input.xml";
// 	std::cout << std::endl << "	step 1";
// 	myImg = ftk::LoadXMLImage( "/media/sf_11_SharedWithUbuntu/Results/Histo_Input_big2.xml" );
// 	myImg = ftk::LoadXMLImage( "/media/sf_11_SharedWithUbuntu/Results/Histo_Input_big_invert2.xml" );
// 	myImg = ftk::LoadXMLImage( "/media/sf_11_SharedWithUbuntu/Results/Histo_Input_big.xml" );
// 	myImg = ftk::LoadXMLImage( "/media/sf_11_SharedWithUbuntu/Results/Histo_Input_big_invert.xml" );
	myImg = ftk::LoadXMLImage( "/media/sf_11_SharedWithUbuntu/Results/Histo_Input_toy.xml" );
// 	myImg = ftk::LoadXMLImage( "/media/sf_11_SharedWithUbuntu/Results/Histo_Input_noise.xml" );
	
//	ftk::Image::DataType dataType;
// 	itk::ImageIOBase::UCHAR
// 	std::cout << std::endl << "	step 2";
// 	std::cout << std::endl << "	step 2";
	std::cout << myImg->GetImageInfo()->dataType;
	
// 	switch(myImg->GetImageInfo()->dataType)
// 	{
// 		case itk::ImageIOBase::USHORT:
// 		{
// 			typedef ftk::nucSecNic::ftkNuclearSegmentationNic<unsigned short> my_nucSecNic;
// 			std::cout << "1111111";
// 			
// 		break;
// 		}
// 		case itk::ImageIOBase::UCHAR:
// 		{
// 			typedef ftk::nucSecNic::ftkNuclearSegmentationNic<unsigned char> my_nucSecNic;
// 			std::cout << "2222222";
// 		break;
// 		}
// 	}
	
// 	if( myImg->GetImageInfo()->dataType == itk::ImageIOBase::USHORT)
// 	{
// 		typedef ftk::nucSecNic::ftkNuclearSegmentationNic<unsigned short> my_nucSecNic;
// 	}
// 	else if(myImg->GetImageInfo()->dataType == itk::ImageIOBase::UCHAR)
// 	{
// 		typedef ftk::nucSecNic::ftkNuclearSegmentationNic<unsigned char> my_nucSecNic;
// 	}	
	
	typedef ftk::nucSecNic::ftkNuclearSegmentationNic<unsigned char, unsigned short> my_nucSecNic;
// 	typedef ftk::nucSecNic::ftkNuclearSegmentationNic<unsigned char, float> my_nucSecNic;
// 	std::cout << std::endl << "	step 3";
// 	std::cout << std::endl << "	step 3";
	my_nucSecNic nucSecNic;
	bool out = nucSecNic.setInput( myImg, "test_data", 0, 0 );
	
// 	std::cout << std::endl << "	step 4";
// 	std::cout << std::endl << "	step 4";
	nucSecNic.ResetRealeaseAll();
// 	std::cout << std::endl << "	step 5";
	nucSecNic.ReleaseMemory();
// 	std::cout << std::endl << "	step 6";
	nucSecNic.setBinarizeMixPoissonParameters_1(128);
	
// 	std::cout << std::endl << "	step 4";
// 	std::cout << std::endl << "	step 4";

	nucSecNic.runBinarizeMixPoisson_1();
	
// 	std::cout << std::endl << "	step 55";
// 	std::cout << std::endl << "	step 54";
	
	nucSecNic.runSeedDetection_1();
	
	nucSecNic.setMaxClusteringParameters_1();
	nucSecNic.runMaxClustering_1();
	

	int * nic;
	nic = new int[2];
	delete [] nic;
	
	
	std::cout << std::endl;
	return 0;
};

// ############################################################################################################################################################################