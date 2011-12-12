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

	ftk::Image::Pointer myImg;
	//std::string inputName = "testImages\Histo_Input.xml";
	myImg = ftk::LoadXMLImage( "/home/nicolasreyv/farsight/src/farsight-src/ftkNuclearSegmentationNic/testImages/Histo_Input.xml" );
	
//	ftk::Image::DataType dataType;
// 	itk::ImageIOBase::UCHAR
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
	
	typedef ftk::nucSecNic::ftkNuclearSegmentationNic<unsigned char> my_nucSecNic;
	my_nucSecNic nucSecNic;
	
	nucSecNic.ResetRealeaseAll();
	nucSecNic.ReleaseMemory();
	nucSecNic.runBinarizeMixPoisson_1();
	
	std::cout << myImg->GetImageInfo()->numColumns;
	
	bool out = nucSecNic.setInput( myImg, "test_data", 0, 0 );
	nucSecNic.runBinarizeMixPoisson_1();

	int * nic;
	nic = new int[2];
	delete [] nic;
	
	
	std::cout << std::endl;
	return 0;
};

// ############################################################################################################################################################################