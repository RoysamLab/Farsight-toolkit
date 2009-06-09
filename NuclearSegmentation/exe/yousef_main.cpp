/*=========================================================================
Copyright 2009 Rensselaer Polytechnic Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 
=========================================================================*/

#include "yousef_core/yousef_seg.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include <time.h>

using namespace std;

int main(int argc, char* argv[])
{
	if(argc <4)
	{
		std::cout<<"Usage: segment_nuclei <InputImageFileName> <OutputImageFileName> <ParametersFileName>\n";
		return 0;
	}
	clock_t startTimer = clock();
	
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
	NucleusSeg->runSeedDetection();
	//3-Initial Segmentation (CLustering)
	NucleusSeg->runClustering();
	if(NucleusSeg->isSegmentationFinEnabled())
	{
		//4-Optional Segmentation Refinement (Alpha-expansion) 
		NucleusSeg->runAlphaExpansion3D();		
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
		for(int i=0; i<(size1*size2*size3); i++)
		{	
			unsigned short val = (unsigned short)output_img[i];
			iterator1.Set(val);			
			++iterator1;	
		}
		
		typedef itk::ImageFileWriter< SegmentedImageType > WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetFileName(argv[2]);
		writer->SetInput( image );
		writer->Update();
		image = 0;
		writer = 0;

		NucleusSeg->outputSeeds();
        NucleusSeg-> saveIntoIDLFormat(argv[1]);
        
		delete NucleusSeg;
		//delete output_img; //This is deleted when I delete NucleusSeg
		delete in_Image;
	//}
	
	std::cout<<"Time elapsed is: "<<(((double)clock() - startTimer) / CLOCKS_PER_SEC)<<" seconds"<<std::endl;
      
	return 1;
}
