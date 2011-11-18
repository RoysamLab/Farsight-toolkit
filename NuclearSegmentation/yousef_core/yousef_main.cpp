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

#include "yousef_seg.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include <time.h>
#include <limits.h>

int main(int argc, char* argv[])
{
	if(argc <4)
	{
		std::cout<<"Usage: segment_nuclei <InputImageFileName> <OutputImageFileName> <ParametersFileName>\n";
		return 0;
	}
	clock_t startTimer = clock();

	std::cout<<"Starting segmentation\n";

	//For 8-bit grayscale images only 
	typedef itk::Image< unsigned char, 3 > OutputImageType;
	typedef itk::ImageFileReader< OutputImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName (argv[1]);
	reader->Update();
	//read the input image into an ITK image
	OutputImageType::Pointer img = reader->GetOutput();	
	//get the image dimensions
	int size1=img->GetLargestPossibleRegion().GetSize()[0];
	int size2=img->GetLargestPossibleRegion().GetSize()[1];
	int size3=img->GetLargestPossibleRegion().GetSize()[2];
	//size3 +=1;
	unsigned char *in_Image;
	in_Image = (unsigned char *) malloc (size1*size2*size3);
	
	if( in_Image == NULL )
	{
		std::cout<<"Memory allocation for image failed\n";
		return 0;
	}
	memset(in_Image/*destination*/,0/*value*/,size1*size2*size3*sizeof(unsigned char)/*num bytes to move*/);

	typedef itk::ImageRegionConstIterator< OutputImageType > ConstIteratorType;
	ConstIteratorType pix_buf( img, img->GetRequestedRegion() );

	int ind=0;
	for ( pix_buf.GoToBegin(); !pix_buf.IsAtEnd(); ++pix_buf, ++ind )
	in_Image[ind]=(pix_buf.Get());

	reader = 0;
	img=0;

	yousef_nucleus_seg *NucleusSeg = new yousef_nucleus_seg();
	NucleusSeg->readParametersFromFile(argv[3]);
	NucleusSeg->setDataImage(in_Image,size1,size2,size3,argv[1]);
	
	
	int *output_img;
	//segmentation steps
	//1-Binarization
	NucleusSeg->runBinarization();
	//2-Seeds Detection
	try {
		std::cout<<"Starting seed detection\n";
		NucleusSeg->runSeedDetection();
	}
	catch( bad_alloc & excp ){
		std::cout<<"You have requested more memory than "
			 <<"what is currently available in this "
			 <<"system, please try again with a smaller "
			 <<"input image\n";
	}
	catch( itk::ExceptionObject & excp ){
		std::cout<<"Error: " << excp <<std::endl;
	}
	catch( std::excption & excp ){
		std::cout<<"Error: " << excp.what() << std::endl;
	}
	//3-Initial Segmentation (CLustering)
	NucleusSeg->runClustering();
	if(NucleusSeg->isSegmentationFinEnabled())
	{
		//4-Optional Segmentation Refinement (Alpha-expansion) 
		NucleusSeg->runAlphaExpansion();		
		output_img=NucleusSeg->getSegImage();
	}
	else
		output_img=NucleusSeg->getClustImage();

	
	//Hard coded parameters, add cin for each if required..
	/*int params[6];
	params[0]=0;
	params[1]=30;
	params[2]=5;
	params[3]=8;
	params[4]=5;
	params[5]=2;

	int choice;
	
	cout<<"Use default parameters?\n1->\t\tNo, I will define the parameters\nAny other key->\tYes, Use default parameters\n";
	cin>>choice;
	if(choice==1)
	{
		cout<<"Enter the parameters for shift sigma scaleMin scaleMax regionXY and regionZ\n";
		cin>>params[0]>>params[1]>>params[2]>>params[3]>>params[4]>>params[5];
	}
	
	cout<<"1->Run Binarization\n2->Run Seed Detection\n3->Run Clustering\n4->Run Aplha Expansion\nAny other key->Exit\n";
	cin>>choice;
	
	if( choice==1 || choice==2 || choice==3 || choice==4 )
	{
		yousef_nucleus_seg *NucleusSeg = new yousef_nucleus_seg();
		NucleusSeg->setParams(params);
		NucleusSeg->setDataImage(in_Image,size1,size2,size3);
		int *output_img;
		NucleusSeg->runBinarization();
		choice--;
		if ( choice )
		{
			NucleusSeg->runSeedDetection();
			choice--;
			if ( choice )
			{
				NucleusSeg->runClustering();
				choice--;
				if ( choice )
				{
					NucleusSeg->runAlphaExpansion3D();		
					output_img=NucleusSeg->getSegImage();			
				}
				else
					output_img=NucleusSeg->getClustImage();
			}
			else
				output_img=NucleusSeg->getSeedImage();
		}
		else
			output_img=NucleusSeg->getBinImage();
		*/

		//size3 -= 1;
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
		int max_pix, min_pix; //Kedar Testing ****************
		min_pix=USHRT_MAX; max_pix=0;
		for(int i=0; i<(size1*size2*size3); i++)
		{	
			unsigned short val = (unsigned short)output_img[i];
			if( min_pix > val ) min_pix = val;
			if( max_pix < val ) max_pix = val;
			iterator1.Set(val);			
			++iterator1;	
		}

		printf("Max pixel value: %i\nMin pixel value: %i\n",max_pix,min_pix);

		typedef itk::ImageFileWriter< SegmentedImageType > WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetFileName(argv[2]);
		writer->SetInput( image );
		writer->Update();
		image = 0;
		writer = 0;

		NucleusSeg->outputSeeds();
        //NucleusSeg-> saveIntoIDLFormat(argv[1]);
        
		delete NucleusSeg;
		//delete output_img; //This is deleted when I delete NucleusSeg
		delete in_Image;
	//}
	
	std::cout<<"Time elapsed is: "<<(((double)clock() - startTimer) / CLOCKS_PER_SEC)<<" seconds"<<std::endl;
      
	return 1;
}
