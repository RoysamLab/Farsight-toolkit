/* 
* Copyright 2009 Rensselaer Polytechnic Institute
* This program is free software; you can redistribute it and/or modify 
* it under the terms of the GNU General Public License as published by 
* the Free Software Foundation; either version 2 of the License, or 
* (at your option) any later version.
* 
* This program is distributed in the hope that it will be useful, but 
* WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
* or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
* for more details.
* 
* You should have received a copy of the GNU General Public License along 
* with this program; if not, write to the Free Software Foundation, Inc., 
* 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*/

#include "yousef_core/yousef_seg.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include <time.h>



int main(int argc, char* argv[])
{
	//This executable is for 16-bit grayscale images only, and assuming a gaussian mixture model instead of a poisson mixture:

	if(argc <4)
	{
	//	std::cout<<"Usage: segment_nuclei <InputImageFileName> <BinaryImageFileName> <OutputImageFileName> <ParametersFileName>\n";
		std::cout<<"Usage: segment_nuclei <InputImageFileName>  <OutputImageFileName> <ParametersFileName>\n";
		return 0;
	}
	typedef itk::Image< unsigned short, 3 > InputImageType16;
	typedef itk::ImageFileReader< InputImageType16 > ReaderType16;
	typedef itk::ImageRegionConstIterator< InputImageType16 > ConstIteratorType16;
	

	// Read Input Image:////////////////////////////////////////////////////////////////////////////////////
	std::cout<<"reading input image...\n";
	ReaderType16::Pointer reader = ReaderType16::New();
	reader->SetFileName (argv[1]);
	try{
		reader->Update();	
	}
	catch( itk::ExceptionObject &err)
	{
		  std::cerr<<err<<std::endl;
		  return EXIT_FAILURE;
	 }	
	std::cout<<"done reading input image..\n";
	InputImageType16::Pointer imagedata = reader->GetOutput();
	//typedef itk::ImageFileWriter< InputImageType16 > WriterType1;
	//WriterType1::Pointer writer1 = WriterType1::New();
	//writer1->SetFileName("C:/Lab/AminFiles/Debug/SegmentationNavin/bin_im.tif");
	//writer1->SetInput( imagedata );
	//writer1->Update();
	//scanf("%*d");

	//////// Read Initial Binary Image for GMM Parameter Estimation://///////////////////////////////////////////
	//////std::cout<<"reading binary image...\n";
	//////ReaderType16::Pointer reader2 = ReaderType16::New();

	//////reader2->SetFileName (argv[2]);
	//////try{
	//////	reader2->Update();	
	//////}
	//////catch( itk::ExceptionObject &err)
	//////{
	//////	  std::cerr<<err<<std::endl;
	//////	  return EXIT_FAILURE;
	////// }	
	//////std::cout<<"done reading binary image..\n";

	//////InputImageType16::Pointer binary_imagedata = reader2->GetOutput();

  	//Copy The Images into  Arrays://////////////////////////////////////////////////////////////////////
	size_t size1=imagedata->GetLargestPossibleRegion().GetSize()[0];
	size_t size2=imagedata->GetLargestPossibleRegion().GetSize()[1];
	size_t size3=imagedata->GetLargestPossibleRegion().GetSize()[2];
	std::cout<<"size_x:\t"<<size1<<std::endl;
	std::cout<<"size_y:\t"<<size2<<std::endl;
	std::cout<<"size_z:\t"<<size3<<std::endl;

	unsigned short *array_imagedata;
	array_imagedata = (unsigned short *) malloc(size1*size2*size3*sizeof(array_imagedata[0]));
	std::cout<<"trying to allocate memory to input image Array..\n";
	std::cout<<size1*size2*size3<<std::endl;
	if( array_imagedata == NULL )
	{
		std::cout<<"Memory allocation for image failed\n";
		return 0;
	}
	//std::cout << array_imagedata[3] << std::endl;
	//std::cout<<"memset is not working..\n";
	std::cout << sizeof(array_imagedata[0]) << std::endl;
	memset(array_imagedata/*destination*/,0/*value*/,size1*size2*size3*sizeof(array_imagedata[0])/*num bytes to move*/);
	std::cout<<"done copying input image..\n" << std::flush;


	//unsigned short *array_binary_imagedata;
	//array_binary_imagedata = (unsigned short *) malloc(size1*size2*size3*sizeof(array_binary_imagedata[0]));
	//std::cout<<"trying to allocate memory to binary image Array..\n";


	//if( array_binary_imagedata == NULL )
	//{
	//	std::cout<<"Memory allocation for image failed\n";
	//	return 0;
	//}
	//memset(array_binary_imagedata/*destination*/,0/*value*/,size1*size2*size3*sizeof(array_binary_imagedata[0])/*num bytes to move*/);
	//std::cout<<"done copying binary image..\n";
	
	ConstIteratorType16 imagedata_iter( imagedata, imagedata->GetRequestedRegion() );
	//ConstIteratorType16 binary_imagedata_iter( binary_imagedata, binary_imagedata->GetRequestedRegion() );
	size_t index=0;
	//for ( imagedata_iter.GoToBegin(); !imagedata_iter.IsAtEnd(); ++imagedata_iter, ++index,++binary_imagedata_iter )
	for ( imagedata_iter.GoToBegin(); !imagedata_iter.IsAtEnd(); ++imagedata_iter, ++index)
	{
		array_imagedata[index]=(imagedata_iter.Get());
		//array_binary_imagedata[index] = (binary_imagedata_iter.Get());
	}
	reader = 0;
	imagedata=0;
	//binary_imagedata=0;

	// Create an instance of Nuclear Segmentation and Start Segmentation: //////////////////////////////////

	yousef_nucleus_seg *NucleusSeg = new yousef_nucleus_seg();
	if(argc == 4)
		NucleusSeg->readParametersFromFile("");
	else
		NucleusSeg->readParametersFromFile(argv[4]);

	// Set Input Image
	NucleusSeg->setDataImage(array_imagedata,size1,size2,size3,argv[1]);
	//NucleusSeg->setInitialBinaryImage(array_binary_imagedata);
	
	if(!NucleusSeg->EstimateGMMParameters())
	{
		NucleusSeg->RunKmeansClustering();
		//std::cout<<"Something is up with the initial binary image\n";
		//return 0;
	}

	unsigned short *output_img;
	//segmentation steps
	//1-Binarization
	NucleusSeg->runBinarization16();
	//2-Seeds Detection
	std::cout << "about to start seed detection........" << std::endl;
	try{
		NucleusSeg->runSeedDetection16();
	}
	catch( std::bad_alloc & excp ){
		std::cout<<"You have requested more memory than "
			 <<"what is currently available in this "
			 <<"system, please try again with a smaller "
			 <<"input image\n";
		return EXIT_FAILURE;
	}
	catch( itk::ExceptionObject & excp ){
		std::cout<<"Error: " << excp <<std::endl;
		return EXIT_FAILURE;
	}
	catch( std::exception & excp ){
		std::cout<<"Error: " << excp.what() << std::endl;
		return EXIT_FAILURE;
	}

	////std::cout << "zackdebug: clustering" << std::endl;
	////3-Initial Segmentation (CLustering)
	NucleusSeg->runClustering16();
	if(NucleusSeg->isSegmentationFinEnabled())
	{
		//4-Optional Segmentation Refinement (Alpha-expansion) 
		//std::cout << "zackdebug: alpha expansion" << std::endl;
		NucleusSeg->runAlphaExpansion16();		
		output_img=NucleusSeg->getSegImage();
	}
	else
		output_img=NucleusSeg->getClustImage();


	typedef itk::Image< unsigned short, 3 > SegmentedImageType;
	SegmentedImageType::Pointer image;
	image = SegmentedImageType::New();
	SegmentedImageType::PointType origin;
	origin[0] = 0; 
	origin[1] = 0;    
	origin[2] = 0;    
	image->SetOrigin( origin );

	SegmentedImageType::IndexType start;
	start[0] = 0;  // first index on X
	start[1] = 0;  // first index on Y    
	start[2] = 0;  // first index on Z    
	SegmentedImageType::SizeType  size;
	size[0] = size1;  // size along X
	size[1] = size2;  // size along Y
	size[2] = size3;  // size along Z

	SegmentedImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );

	image->SetRegions( region );
	image->Allocate();
	image->FillBuffer(0);
	image->Update();

	//copy the output image into the ITK image
	typedef itk::ImageRegionIteratorWithIndex< SegmentedImageType > IteratorType;
	IteratorType iterator1(image,image->GetRequestedRegion());		
	for(size_t i=0; i<(size1*size2*size3); i++)
	{				
		unsigned short val = (unsigned short)output_img[i];
		iterator1.Set(val);			
		++iterator1;				
	}
	std::cout<< "About to write result..\n";
	typedef itk::ImageFileWriter< SegmentedImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(argv[3]);
	writer->SetInput( image );
	writer->Update();
	image = 0;
	writer = 0;

//	NucleusSeg->outputSeeds();
//	NucleusSeg-> saveIntoIDLFormat(argv[1]);

	delete NucleusSeg;
	//delete output_img; //This is deleted when I delete NucleusSeg
	delete array_imagedata;
	//}

	//std::cout<<"Time elapsed is: "<<(((double)clock() - startTimer) / CLOCKS_PER_SEC)<<" seconds"<<std::endl;

	return 0;
}
