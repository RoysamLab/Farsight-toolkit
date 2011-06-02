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

#include "SomaExtraction.h"


SomaExtractor::SomaExtractor()
{

}

void SomaExtractor::SetInputImage(char * fileName)
{
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName (fileName);
	reader->Update();	
	//read the input image into an ITK image
	OutputImageType::Pointer img = reader->GetOutput();	
	this->SetInputImage(img);
	
}

void SomaExtractor::SetInputImage(OutputImageType::Pointer img)
{
	if(!img) return;

	inputImage = img;
	size1 = inputImage->GetLargestPossibleRegion().GetSize()[0];
	size2 = inputImage->GetLargestPossibleRegion().GetSize()[1];
	size3 = inputImage->GetLargestPossibleRegion().GetSize()[2];
	//size3 +=1;
	in_Image = (unsigned char *) malloc (size1*size2*size3);
	
	if( in_Image == NULL )
	{
		std::cout<<"Memory allocation for image failed\n";
		return ;
	}

	memset(in_Image/*destination*/,0/*value*/,size1*size2*size3*sizeof(unsigned char)/*num bytes to move*/);

	ConstIteratorType pix_buf( inputImage, inputImage->GetRequestedRegion() );

	int ind=0;
	for ( pix_buf.GoToBegin(); !pix_buf.IsAtEnd(); ++pix_buf, ++ind )
	in_Image[ind]=(pix_buf.Get());


}

bool SomaExtractor::LoadSegParams(int kernel, int minObj)
{
	KernelSize = kernel;
	minObjSize = minObj;
	//TiXmlDocument doc;
	//if ( !doc.LoadFile( filename ) )
	//	return false;

	//TiXmlElement* rootElement = doc.FirstChildElement();
	//const char* docname = rootElement->Value();
	//if ( strcmp( docname, "SomaExtractionDefinition" ) != 0 )
	//	return false;

	//std::string name = rootElement->Attribute("name");//default

	//TiXmlElement * parentElement = rootElement->FirstChildElement();
	//const char * parent = parentElement->Value();
	//	if ( strcmp(parent,"SegParameters") == 0 )
	//	{
	//		TiXmlElement * parameterElement = parentElement->FirstChildElement();
	//		const char * parameter = parameterElement ->Value();
	//		this->KernelSize = atoi(parameterElement->Attribute("KERNEL_SIZE"));
	//		this->minObjSize = atoi(parameterElement->Attribute("MIN_OBJ_SIZE"));
	//		std::cout<<"minObjSize = "<<minObjSize<<std::endl;

	//	}
	//	else 
	//	{
	//		printf("Wrong Format for Parameter File\n");
	//	}
	//return true;
	return true;
}

int SomaExtractor::binarizeImage(char* paramFile)
{
	NucleusSeg = new yousef_nucleus_seg(); //
	if(!paramFile)
		NucleusSeg->readParametersFromFile("");
	else
	NucleusSeg->readParametersFromFile(paramFile);
	NucleusSeg->setDataImage(in_Image, size1, size2, size3, "");
	
	
	unsigned short *output_img;
	//segmentation steps
	//1-Binarization
	
	NucleusSeg->runBinarization();
	output_img = NucleusSeg->getBinImage();	
	
	binImage = SegmentedImageType::New();
	SegmentedImageType::PointType origin;
   	origin[0] = 0; 
    origin[1] = 0;    
	origin[2] = 0;    
    binImage->SetOrigin( origin );

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

	binImage->SetRegions( region );
   	binImage->Allocate();
   	binImage->FillBuffer(0);
	binImage->Update();

	//copy the output image into the ITK image
	IteratorType iterator1(binImage, binImage->GetRequestedRegion());		
	for(int i=0; i<(size1*size2*size3); i++)
	{				
		unsigned short val = (unsigned short)output_img[i];
		iterator1.Set(val);			
		++iterator1;				
	}

	KernelType ball;
	KernelType::SizeType ballSize;
	//ballSize.Fill( atoi(argv[4]) ); //for now, set the radius to 4
	ballSize[0] = ( KernelSize );
	ballSize[1] = ( KernelSize );
	ballSize[2] = ( 2 );
	ball.SetRadius( ballSize );
	ball.CreateStructuringElement();
	
	morphOpenFilter = morphOpenFilterType::New();
	morphOpenFilter->SetInput( binImage );
	morphOpenFilter->SetKernel( ball );
	//filter->SetForegroundValue( 255 ); 

    try
    {
		morphOpenFilter->Update();
    } 
    catch ( itk::ExceptionObject & excp )
    {
		std::cerr << excp << std::endl;
		return EXIT_FAILURE;
    }
	binImage = morphOpenFilter->GetOutput();
		
}

int SomaExtractor::relabelBinaryImage(void)
{
	RelabelFilterType::Pointer relabel = RelabelFilterType::New();
	relabel->SetInput( binImage );
	relabel->SetMinimumObjectSize( minObjSize );  //Still have to read the segParams file before this
	//std::cout<<"minObjSize = "<<minObjSize<<std::endl;
	relabel->InPlaceOn();

	try
    {
		relabel->Update();
    }
    catch( itk::ExceptionObject & excep )
    {
		std::cerr << "Relabel: soma extraction exception caught !" << std::endl;
		std::cerr << excep << std::endl;
		return EXIT_FAILURE;
    }

	BinaryThresholdImageType::Pointer  BinaryThreshold = BinaryThresholdImageType::New();
	BinaryThreshold->SetInput( binImage );
	BinaryThreshold->SetLowerThreshold( 1 );
	//BinaryThresholdImage->SetUpperThreshold( 255 );	//save for potential later use
	BinaryThreshold->SetInsideValue(255);
	BinaryThreshold->SetOutsideValue(0);
	BinaryThreshold->Update();
	outputImage = BinaryThreshold->GetOutput();
	
}

SomaExtractor::centroidVectorType SomaExtractor::GetSomaCentroids()
{
	ConverterType::Pointer converter = ConverterType::New();
	converter->SetInputForegroundValue(255);
	converter->SetInput(outputImage);
	LabelMapType::Pointer Somas = converter->GetOutput();
	converter->Update();
	unsigned int numSomas = Somas->GetNumberOfLabelObjects();

	for(unsigned int label=1; label<= numSomas; ++label)
	{
		const LabelObjectType * labelObject = Somas->GetLabelObject(label);
		if(labelObject->GetPhysicalSize() < minObjSize)
		{	
			//skip small blobs: they probably aren't real somas
			std::cout<< "rejected cell "<< label <<" of size " << labelObject->GetPhysicalSize() << std::endl;
			continue;
		}
		const LabelObjectType::CentroidType centroid = labelObject->GetCentroid();
		OutputImageType::IndexType pixelIndex;
		outputImage->TransformPhysicalPointToIndex( centroid, pixelIndex );
		Centroids.push_back(pixelIndex);
	}
	return Centroids;
}

void SomaExtractor::writeSomaCentroids(char* writeFileName)
{
	std::string FileName = writeFileName;
	FileName.erase(FileName.length() -4, FileName.length());
	FileName.append(".txt");
	std::ofstream outfile(FileName.c_str());
	outfile.precision(1);
	for(int i=0; i<(int)Centroids.size(); ++i)
	{
		OutputImageType::IndexType Index = Centroids.at(i);
		outfile << std::fixed << Index[0] << " " << Index[1] << " " << Index[2] << std::endl;
	}
	outfile.close();
}

void SomaExtractor::writeSomaImage(char* writeFileName)
{
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(writeFileName);
	writer->SetInput(outputImage);
	writer->Update();
}



//
//int main(int argc, char* argv[])
//{	
//
//	if(argc <6)
//	{
//		std::cout<<"Usage: segment_nuclei <InputImageFileName> <OutputImageFileName> <ParametersFileName> <Opening Kernel Size> <Min Object Size>\n";
//		return 0;
//	}
//	clock_t startTimer = clock();
//	
//	std::cout<<"reading input image...";
//	//For 8-bit grayscale images only 
//	/*typedef itk::Image< unsigned char, 3 > OutputImageType;
//	typedef itk::ImageFileReader< OutputImageType > ReaderType;*/
//	ReaderType::Pointer reader = ReaderType::New();
//	reader->SetFileName (argv[1]);
//	reader->Update();	
//	std::cout<<"done"<<std::endl;
//	//read the input image into an ITK image
//	OutputImageType::Pointer img = reader->GetOutput();	
//	//get the image dimensions
//	int size1=img->GetLargestPossibleRegion().GetSize()[0];
//	int size2=img->GetLargestPossibleRegion().GetSize()[1];
//	int size3=img->GetLargestPossibleRegion().GetSize()[2];
//	//size3 +=1;
//	unsigned char *in_Image;
//	in_Image = (unsigned char *) malloc (size1*size2*size3);
//	
//	if( in_Image == NULL )
//	{
//		std::cout<<"Memory allocation for image failed\n";
//		return 0;
//	}
//	memset(in_Image/*destination*/,0/*value*/,size1*size2*size3*sizeof(unsigned char)/*num bytes to move*/);
//
//	//typedef itk::ImageRegionConstIterator< OutputImageType > ConstIteratorType;
//	ConstIteratorType pix_buf( img, img->GetRequestedRegion() );
//
//	int ind=0;
//	for ( pix_buf.GoToBegin(); !pix_buf.IsAtEnd(); ++pix_buf, ++ind )
//	in_Image[ind]=(pix_buf.Get());
//
//	//reader = 0;
//	//img=0;
//
//	yousef_nucleus_seg *NucleusSeg = new yousef_nucleus_seg();
//	if(argc == 3)
//		NucleusSeg->readParametersFromFile("");
//	else
//		NucleusSeg->readParametersFromFile(argv[3]);
//	NucleusSeg->setDataImage(in_Image,size1,size2,size3,argv[1]);
//	
//	
//	unsigned short *output_img;
//	//segmentation steps
//	//1-Binarization
//	NucleusSeg->runBinarization();
//	output_img=NucleusSeg->getBinImage();
//	
// 
//		//typedef itk::Image< unsigned short, 3 > SegmentedImageType;
//		SegmentedImageType::Pointer image;
//		image = SegmentedImageType::New();
//		SegmentedImageType::PointType origin;
//   		origin[0] = 0; 
//    	origin[1] = 0;    
//		origin[2] = 0;    
//    	image->SetOrigin( origin );
//
//    	SegmentedImageType::IndexType start;
//    	start[0] = 0;  // first index on X
//    	start[1] = 0;  // first index on Y    
//		start[2] = 0;  // first index on Z    
//    	SegmentedImageType::SizeType  size;
//    	size[0] = size1;  // size along X
//    	size[1] = size2;  // size along Y
//		size[2] = size3;  // size along Z
//
//    	SegmentedImageType::RegionType region;
//    	region.SetSize( size );
//    	region.SetIndex( start );
//
//    	image->SetRegions( region );
//    	image->Allocate();
//    	image->FillBuffer(0);
//		image->Update();
//
//		//copy the output image into the ITK image
//		//typedef itk::ImageRegionIteratorWithIndex< SegmentedImageType > IteratorType;
//		IteratorType iterator1(image,image->GetRequestedRegion());		
//		for(int i=0; i<(size1*size2*size3); i++)
//		{				
//			unsigned short val = (unsigned short)output_img[i];
//			iterator1.Set(val);			
//			++iterator1;				
//		}
//
//
//		//////////////////////////////////////////////////////////////////////////////////////////
//		
//
//		
//
//	//Binary Morphological closing....
//	//1-Create the structuring element
//	//typedef itk::BinaryBallStructuringElement< unsigned short, 3 > KernelType;
//	KernelType ball;
//	KernelType::SizeType ballSize;
//	//ballSize.Fill( atoi(argv[4]) ); //for now, set the radius to 4
//	ballSize[0] = ( atoi(argv[4]) );
//	ballSize[1] = ( atoi(argv[4]) );
//	ballSize[2] = ( 2 );
//	ball.SetRadius( ballSize );
//	ball.CreateStructuringElement();
//	
//	//typedef itk::GrayscaleMorphologicalOpeningImageFilter< SegmentedImageType, SegmentedImageType, KernelType > FilterType;
//	FilterType::Pointer filter = FilterType::New();
//	filter->SetInput( image );
//	filter->SetKernel( ball );
//	//filter->SetForegroundValue( 255 ); 
//
//    try
//    {
//		filter->Update();
//    } 
//    catch ( itk::ExceptionObject & excp )
//    {
//		std::cerr << excp << std::endl;
//		return EXIT_FAILURE;
//    }
//		
//	
//	//typedef itk::RelabelComponentImageFilter< SegmentedImageType, SegmentedImageType > RelabelFilterType;
//	RelabelFilterType::Pointer relabel = RelabelFilterType::New();
//	relabel->SetInput( filter->GetOutput() );
//	relabel->SetMinimumObjectSize( atoi(argv[5]) );
//	relabel->InPlaceOn();
//
//	try
//    {
//		relabel->Update();
//    }
//    catch( itk::ExceptionObject & excep )
//    {
//		std::cerr << "Relabel: exception caught !" << std::endl;
//		std::cerr << excep << std::endl;
//		return 1;
//    }
//	//typedef itk::BinaryThresholdImageFilter<SegmentedImageType, OutputImageType> BinaryThresholdImageType;
//	BinaryThresholdImageType::Pointer  BinaryThresholdImage = BinaryThresholdImageType::New();
//	BinaryThresholdImage->SetInput( filter->GetOutput() );
//	BinaryThresholdImage->SetLowerThreshold( 1 );
//	//BinaryThresholdImage->SetUpperThreshold( 255 );	//save for potential later use
//	BinaryThresholdImage->SetInsideValue(255);
//	BinaryThresholdImage->SetOutsideValue(0);
//
//	//typedef unsigned long LabelType;
//	//typedef itk::ShapeLabelObject< LabelType, 3 > LabelObjectType;
//	//typedef itk::LabelMap< LabelObjectType > LabelMapType;
//	//typedef itk::BinaryImageToShapeLabelMapFilter< OutputImageType, LabelMapType >	ConverterType;
//	ConverterType::Pointer converter = ConverterType::New();
//	converter->SetInputForegroundValue(255);
//	converter->SetInput(BinaryThresholdImage->GetOutput());
//	LabelMapType::Pointer Somas = converter->GetOutput();
//	converter->Update();
//	unsigned int numSomas = Somas->GetNumberOfLabelObjects();
//	std::string FileName = argv[2];
//	FileName.erase(FileName.length() -4, FileName.length());
//	FileName.append(".txt");
//	std::ofstream outfile(FileName.c_str());
//	outfile.precision(1);
//
//	for(unsigned int label=1; label<= numSomas; ++label)
//	{
//		const LabelObjectType * labelObject = Somas->GetLabelObject(label);
//		if(labelObject->GetPhysicalSize() < 100)
//		{	//skip small blobs: they probably aren't real somas
//			std::cout<< "rejected cell "<< label <<" of size " << labelObject->GetPhysicalSize() << std::endl;
//			continue;
//		}
//		const LabelObjectType::CentroidType centroid = labelObject->GetCentroid();
//		OutputImageType::IndexType pixelIndex;
//		BinaryThresholdImage->GetOutput()->TransformPhysicalPointToIndex( centroid, pixelIndex );
//		outfile << std::fixed << pixelIndex[0] << " " << pixelIndex[1] << " " << pixelIndex[2] << endl;
//	}
//	outfile.close();
////save soma image
//	//typedef itk::ImageFileWriter< OutputImageType > WriterType;
//	WriterType::Pointer writer = WriterType::New();
//	writer->SetFileName(argv[2]);
//	writer->SetInput( BinaryThresholdImage->GetOutput() );
//	writer->Update();
//	image = 0;
//	writer = 0;
//
//	     
//	delete NucleusSeg;
//	//delete output_img; //This is deleted when I delete NucleusSeg
//	delete in_Image;
//	std::cout<<"Time elapsed is: "<<(((double)clock() - startTimer) / CLOCKS_PER_SEC)<<" seconds"<<std::endl;
//      
//	return 0;
//}
