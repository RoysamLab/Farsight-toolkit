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
#include "itkBinaryBallStructuringElement.h"
#include "itkGrayscaleMorphologicalOpeningImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include <iostream>
#include <time.h>


int main(int argc, char* argv[])
{	

	if(argc <6)
	{
		std::cout<<"Usage: segment_nuclei <InputImageFileName> <OutputImageFileName> <ParametersFileName> <Opening Kernel Size> <Min Object Size>\n";
		return 0;
	}
	clock_t startTimer = clock();
	
	std::cout<<"reading input image...";
	//For 8-bit grayscale images only 
	typedef itk::Image< unsigned char, 3 > OutputImageType;
	typedef itk::ImageFileReader< OutputImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName (argv[1]);
	reader->Update();	
	std::cout<<"done"<<std::endl;
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
	if(argc == 3)
		NucleusSeg->readParametersFromFile("");
	else
		NucleusSeg->readParametersFromFile(argv[3]);
	NucleusSeg->setDataImage(in_Image,size1,size2,size3,argv[1]);
	
	
	unsigned short *output_img;
	//segmentation steps
	//1-Binarization
	NucleusSeg->runBinarization();
	output_img=NucleusSeg->getBinImage();
	
 
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
		

		

	//Binary Morphological closing....
	//1-Create the structuring element
	typedef itk::BinaryBallStructuringElement< unsigned short, 3 > KernelType;
	KernelType ball;
	KernelType::SizeType ballSize;
	//ballSize.Fill( atoi(argv[4]) ); //for now, set the radius to 4
	ballSize[0] = ( atoi(argv[4]) );
	ballSize[1] = ( atoi(argv[4]) );
	ballSize[2] = ( 2 );
	ball.SetRadius( ballSize );
	ball.CreateStructuringElement();
	
	typedef itk::GrayscaleMorphologicalOpeningImageFilter< SegmentedImageType, SegmentedImageType, KernelType > FilterType;
	FilterType::Pointer filter = FilterType::New();
	filter->SetInput( image );
	filter->SetKernel( ball );
	//filter->SetForegroundValue( 255 ); 

    try
    {
		filter->Update();
    } 
    catch ( itk::ExceptionObject & excp )
    {
		std::cerr << excp << std::endl;
		return EXIT_FAILURE;
    }
		
	
	typedef itk::RelabelComponentImageFilter< SegmentedImageType, SegmentedImageType > RelabelFilterType;
	RelabelFilterType::Pointer relabel = RelabelFilterType::New();
	relabel->SetInput( filter->GetOutput() );
	relabel->SetMinimumObjectSize( atoi(argv[5]) );
	relabel->InPlaceOn();

	try
    {
		relabel->Update();
    }
    catch( itk::ExceptionObject & excep )
    {
		std::cerr << "Relabel: exception caught !" << std::endl;
		std::cerr << excep << std::endl;
		return 1;
    }
	
		typedef itk::ImageFileWriter< SegmentedImageType > WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetFileName(argv[2]);
		writer->SetInput( filter->GetOutput() );
		writer->Update();
		image = 0;
		writer = 0;

		     
		delete NucleusSeg;
		//delete output_img; //This is deleted when I delete NucleusSeg
		delete in_Image;
	//}
	
	std::cout<<"Time elapsed is: "<<(((double)clock() - startTimer) / CLOCKS_PER_SEC)<<" seconds"<<std::endl;
      
	return 0;
}
