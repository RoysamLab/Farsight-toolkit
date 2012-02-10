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
	SetInputImage(img);
	
}

void SomaExtractor::SetInputImage(OutputImageType::Pointer img)
{
	if(!img) return;

	inputImage = img;
	size1 = inputImage->GetLargestPossibleRegion().GetSize()[0];
	size2 = inputImage->GetLargestPossibleRegion().GetSize()[1];
	size3 = inputImage->GetLargestPossibleRegion().GetSize()[2];
	//size3 +=1;

	std::cout << "Allocating memory of: " << size1*size2*size3 / (1024.0 * 1024.0) << " MB"  << std::endl;
	in_Image = (unsigned char *) malloc (size1*size2*size3*sizeof(unsigned char));
	
	if( in_Image == NULL )
	{
		std::cout<<"Memory allocation for image failed\n";
		return ;
	}
	std::cout << "Done allocating memory" << std::endl;

	std::cout << "Copying 0s into new array" << std::endl;
	memset(in_Image/*destination*/,0/*value*/,size1*size2*size3*sizeof(unsigned char)/*num bytes to move*/);
	std::cout << "Done copying 0s into image" << std::endl;

	ConstIteratorType pix_buf( inputImage, inputImage->GetRequestedRegion() );
	
	std::cout << "Copying image into in_Image" << std::endl;
	size_t ind=0;
	for ( pix_buf.GoToBegin(); !pix_buf.IsAtEnd(); ++pix_buf, ++ind )
		in_Image[ind]=(pix_buf.Get());

	std::cout << "Copied " << ind << " pixels into in_Image" << std::endl;


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

int SomaExtractor::binarizeImage(char* paramFile, unsigned short num_bins)
{
	std::cout << "Constructing youself_nucleus_seg" << std::endl;
	NucleusSeg = new yousef_nucleus_seg(); //
	
	std::cout << "Reading parameters file" << std::endl;
	if(!paramFile)
		NucleusSeg->readParametersFromFile("");
	else
	NucleusSeg->readParametersFromFile(paramFile);
	
	std::cout << "Entering setDataImage" << std::endl;
	NucleusSeg->setDataImage(in_Image, size1, size2, size3, "");
	
	
	unsigned short *output_img;
	//segmentation steps
	//1-Binarization
	
	std::cout << "Entering runBinarization" << std::endl;
	NucleusSeg->runBinarization(num_bins);
	
	std::cout << "Entering getBinImage" << std::endl;
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
	for(size_t i=0; i<(size1*size2*size3); i++)
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

	std::cout << "Running GrayscaleMorphologicalOpeningImageFilter" << std::endl;
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
	std::cout << "Done with binarizeImage" << std::endl;	
}

int SomaExtractor::relabelBinaryImage(void)
{
	RelabelFilterType::Pointer relabel = RelabelFilterType::New();
	relabel->SetInput( binImage );
	relabel->SetMinimumObjectSize( minObjSize );  //Still have to read the segParams file before this
	//std::cout<<"minObjSize = "<<minObjSize<<std::endl;
	
	//relabel->InPlaceOn();

	//Calculate labels
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

	somaImage = relabel->GetOutput();
	
	std::cout << "Originally there were " << relabel->GetOriginalNumberOfObjects() << " objects" << std::endl;
	std::cout << "After relabel there are now " << relabel->GetNumberOfObjects() << " objects" << std::endl;

	/*typedef itk::LabelGeometryImageFilter<SegmentedImageType> LabelGeometryType;
	LabelGeometryType::Pointer labelGeometryFilter = LabelGeometryType::New();
	labelGeometryFilter->SetInput(binImage);
	labelGeometryFilter->SetCalculateOrientedBoundingBox(false);
	labelGeometryFilter->SetCalculateOrientedIntensityRegions(false);
 	labelGeometryFilter->SetCalculateOrientedLabelRegions(false);
	labelGeometryFilter->SetCalculatePixelIndices(true);
	std::cout << "Running labelGeometryFilter->Update" << std::endl;
	try
    {
		labelGeometryFilter->Update();
    }
    catch( itk::ExceptionObject & excep )
    {
		std::cerr << "labelGeometryFilter: soma extraction exception caught !" << std::endl;
		std::cerr << excep << std::endl;
		return EXIT_FAILURE;
    }
	
	//Reject small cells

	//Get the vector containing all the labels of the image
	std::cout << "Getting labelVector" << std::endl;
	std::vector<LabelGeometryType::LabelPixelType> labelVector = labelGeometryFilter->GetLabels();
	std::vector<LabelGeometryType::LabelPixelType>::iterator labelVectorIterator;

	//For each label
	std::cout << "Iterating through all the labels" << std::endl;
	for (labelVectorIterator = labelVector.begin(); labelVectorIterator != labelVector.end(); labelVectorIterator++)
	{
		LabelGeometryType::LabelPixelType label = *labelVectorIterator;
		//Get the volume of that label
		unsigned long volume = labelGeometryFilter->GetVolume(label);
		
		//If the volume is smaller than minObjSize
		if (volume < minObjSize)
		{
			std::cout << "Rejected cell " << label << " of size " << volume << std::endl;
			//Get the index of all the pixels in that label with small value
			LabelGeometryType::LabelIndicesType indices = labelGeometryFilter->GetPixelIndices(label);
			LabelGeometryType::LabelIndicesType::iterator indexIterator;

			//And Reject it by setting all the pixel values to 0
			for (indexIterator = indices.begin(); indexIterator != indices.end(); indexIterator++)
			{
				LabelGeometryType::LabelIndexType index = *indexIterator;
				binImage->SetPixel(index, 0);
			}
		}
	}*/


	std::cout << "Creating Binary Threshold Image" << std::endl;
	BinaryThresholdImageType::Pointer  BinaryThreshold = BinaryThresholdImageType::New();
	BinaryThreshold->SetInput( binImage );
	BinaryThreshold->SetLowerThreshold( 1 );
	//BinaryThresholdImage->SetUpperThreshold( 255 );	//save for potential later use
	BinaryThreshold->SetInsideValue(255);
	BinaryThreshold->SetOutsideValue(0);
	BinaryThreshold->Update();
	outputImage = BinaryThreshold->GetOutput();

	return 1;	
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
		/*if(labelObject->GetPhysicalSize() < minObjSize)
		{	
			//skip small blobs: they probably aren't real somas
			std::cout<< "rejected cell "<< label <<" of size " << labelObject->GetPhysicalSize() << std::endl;
			continue;
		}*/
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
	try
	{
		writer->Update();
	}
	catch ( itk::ExceptionObject &err)
	{
		std::cout << "Error in writer: " << err << std::endl; 
	}
	
	
	somaImageWriter::Pointer soma_image_writer = somaImageWriter::New();
	soma_image_writer->SetFileName("somaImage.mhd");
	soma_image_writer->SetInput(somaImage);
	try
	{
		soma_image_writer->Update();
	}
	catch ( itk::ExceptionObject &err)
	{
		std::cout << "Error in bin_image_writer: " << err << std::endl; 
	}
}