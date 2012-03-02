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
#include "ftkUtils.h"
#include "itkShapeDetectionLevelSetImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkFastMarchingImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

SomaExtractor::SomaExtractor()
{
}

void SomaExtractor::SetInputImage(char * fileName)
{
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName (fileName);
	reader->Update();	
	//read the input image into an ITK image
	OutputImageType::Pointer inputImage = reader->GetOutput();	

	typedef itk::CastImageFilter<OutputImageType, SegmentedImageType> CasterType;
    CasterType::Pointer caster = CasterType::New();
    caster->SetInput(inputImage);
	caster->Update();
	binImage = caster->GetOutput();

	SM = binImage->GetLargestPossibleRegion().GetSize()[0];
    SN = binImage->GetLargestPossibleRegion().GetSize()[1];
    SZ = binImage->GetLargestPossibleRegion().GetSize()[2];
}

bool SomaExtractor::LoadSegParams(int kernel, int minObj)
{
	KernelSize = kernel;
	minObjSize = minObj;

	return true;
}

int SomaExtractor::GenerateSeedPoints(char* paramFile, unsigned short num_bins)
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
	std::cout << "Entering runBinarization" << std::endl;
	NucleusSeg->runBinarization(num_bins);

	NucleusSeg->runSeedDetection();
	NucleusSeg->outputSeeds();
	seeds = NucleusSeg->getSeedsList();

	std::cout << "Writing Seeds into image file"<<endl;
	SegmentedImageType::PointType origin;
   	origin[0] = 0; 
    origin[1] = 0;    
	origin[2] = 0;    

    SegmentedImageType::IndexType start;
    start[0] = 0;  // first index on X
    start[1] = 0;  // first index on Y    
	start[2] = 0;  // first index on Z    
    SegmentedImageType::SizeType  size;
    size[0] = size1;  // size along X
    size[1] = size2;  // size along Y
	size[2] = size3;  // size along Z
	//
	SegmentedImageType::RegionType region;
   	region.SetSize( size );
   	region.SetIndex( start );

	OutputImageType::Pointer seedImage = OutputImageType::New();

    seedImage->SetOrigin( origin );
	seedImage->SetRegions( region );
   	seedImage->Allocate();
   	seedImage->FillBuffer(0);
	seedImage->Update();

	//set the seedpoints into the image	
	for(int i=0; i < seeds.size(); i++)
	{	
		Seed seedIndex = seeds[i];
		OutputImageType::IndexType pixelIndex;
		pixelIndex[0] = seedIndex.x();
		pixelIndex[1] = seedIndex.y();
		pixelIndex[2] = seedIndex.z();
			
		seedImage->SetPixel(pixelIndex, 255);
	}
	WriterType::Pointer image_writer = WriterType::New();
	image_writer->SetFileName("SeedImage.mhd");
	image_writer->SetInput(seedImage);
	try
	{
		image_writer->Update();
	}
	catch ( itk::ExceptionObject &err)
	{
		std::cout << "Error in bin_image_writer: " << err << std::endl; 
	} 
	return EXIT_SUCCESS;
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

void SomaExtractor::SegmentSoma(char* centroidsFileName, const char * somafileName, int timethreshold, double curvatureScaling, double rmsThres)
{
	std::cout<< "Load Centroids Table"<<endl;
	vtkSmartPointer<vtkTable> centroidsTable = ftk::LoadXYZTable( std::string(centroidsFileName));
	
	PointList3D centroidsList;

	for( vtkIdType i = 0; i < centroidsTable->GetNumberOfRows(); i++)
	{
		float x = centroidsTable->GetValue( i, 0).ToDouble();
		float y = centroidsTable->GetValue( i, 1).ToDouble();
		float z = centroidsTable->GetValue( i, 2).ToDouble();
		centroidsList.AddPt(x, y, z);
	}

	ImFastMarching_Soma(centroidsList, timethreshold, curvatureScaling, rmsThres, somafileName);
}

void SomaExtractor::ImFastMarching_Soma(PointList3D seg_seeds, int timeThreshold, double curvatureScaling, double rmsError, const char *somaFileName)
{
	typedef itk::NearestNeighborInterpolateImageFunction< SegmentedImageType, float>  InterpolatorType;
	InterpolatorType::Pointer I_Interpolator = InterpolatorType::New();
	I_Interpolator->SetInputImage(binImage);

	//move the picked points along z axis
    for( int i = 0; i < seg_seeds.NP; i++ )
	{
		seg_seeds.Pt[i].check_out_of_range_3D(SM,SN,SZ);
		SegmentedImageType::IndexType index;
		SegmentedImageType::IndexType index1;
		index[0] = index1[0] = seg_seeds.Pt[i].x;
		index[1] = index1[1] = seg_seeds.Pt[i].y;
		index[2] = seg_seeds.Pt[i].z;
		for( int j = 0; j < SZ; j++ )
		{
			index1[2] = j;
			if( I_Interpolator->EvaluateAtIndex(index1) > I_Interpolator->EvaluateAtIndex(index) )
			{
				seg_seeds.Pt[i].z = j;
				index[2] = j;
			}
		}
	}

	std::cout << "RescaleIntensity"<<endl;
	typedef itk::RescaleIntensityImageFilter< SegmentedImageType, ProbImageType> RescaleFilterType;
    RescaleFilterType::Pointer rescale = RescaleFilterType::New();
	rescale->SetInput( binImage); 
    rescale->SetOutputMinimum( 0 );
    rescale->SetOutputMaximum( 1 );
    rescale->Update();

	std::cout<< "Generating Distance Map..." <<endl;

	clock_t SomaExtraction_start_time = clock();
	typedef  itk::FastMarchingImageFilter< ProbImageType, ProbImageType >    FastMarchingFilterType;
	FastMarchingFilterType::Pointer  fastMarching = FastMarchingFilterType::New();

	typedef FastMarchingFilterType::NodeContainer           NodeContainer;
    typedef FastMarchingFilterType::NodeType                NodeType;
    NodeContainer::Pointer seeds = NodeContainer::New();
    seeds->Initialize();

	for( int i = 0; i< seg_seeds.NP; i++ )
	{
		ProbImageType::IndexType  seedPosition;
		seedPosition[0] = seg_seeds.Pt[i].x;
		seedPosition[1] = seg_seeds.Pt[i].y;
		seedPosition[2] = seg_seeds.Pt[i].z;
		NodeType node;
		const double seedValue = -3;
		node.SetValue( seedValue );
		node.SetIndex( seedPosition );
		seeds->InsertElement( i, node );
	}

	fastMarching->SetTrialPoints(  seeds);
    fastMarching->SetOutputSize( binImage->GetBufferedRegion().GetSize() );
	const double stoppingTime = timeThreshold * 1.1;
    fastMarching->SetStoppingValue(  stoppingTime);
	fastMarching->SetSpeedConstant( 1.0 );
	fastMarching->Update();
	std::cout<< "Total time for Distance Map is: " << (clock() - SomaExtraction_start_time) / (float) CLOCKS_PER_SEC << std::endl;

	/// Shape Detection
	SomaExtraction_start_time = clock();
	std::cout<< "Shape Detection..." <<std::endl;
	typedef itk::ShapeDetectionLevelSetImageFilter< ProbImageType, ProbImageType> ShapeDetectionFilterType;
	ShapeDetectionFilterType::Pointer shapeDetection = ShapeDetectionFilterType::New();
	
	shapeDetection->SetPropagationScaling(1.0);
	shapeDetection->SetCurvatureScaling(curvatureScaling);
	shapeDetection->SetMaximumRMSError( rmsError);
	shapeDetection->SetNumberOfIterations( 800);

	shapeDetection->SetInput( fastMarching->GetOutput());
	shapeDetection->SetFeatureImage( rescale->GetOutput());
	shapeDetection->Update();
	std::cout << "No. elpased iterations: " << shapeDetection->GetElapsedIterations() << std::endl;
	std::cout << "RMS change: " << shapeDetection->GetRMSChange() << std::endl;
	std::cout<< "Total time for Shape Detection is: " << (clock() - SomaExtraction_start_time) / (float) CLOCKS_PER_SEC << std::endl;

	std::cout<< "Thresholding..."<<endl;
	typedef itk::BinaryThresholdImageFilter< ProbImageType, SegmentedImageType> ShapeThresholdingFilterType;
    ShapeThresholdingFilterType::Pointer thresholder = ShapeThresholdingFilterType::New();
	thresholder->SetLowerThreshold( -10000);
	thresholder->SetUpperThreshold(0.0);
	thresholder->SetOutsideValue( 0);
	thresholder->SetInsideValue(255);
	thresholder->SetInput( shapeDetection->GetOutput());
	thresholder->Update();

	typedef itk::CastImageFilter<SegmentedImageType, OutputImageType> CasterType;
	CasterType::Pointer caster = CasterType::New(); 
	caster->SetInput(thresholder->GetOutput());

	WriterType::Pointer writer1 = WriterType::New();
	writer1->SetInput( caster->GetOutput());
	writer1->SetFileName(somaFileName);
	writer1->Update();
}