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
#include "yousef_core/yousef_seg.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include <itkShapeLabelObject.h>
#include "itkCurvatureAnisotropicDiffusionImageFilter.h"

#include <itkSigmoidImageFilter.h>
#include <math.h>
#include "itkImageSliceIteratorWithIndex.h"
#include "itkDivideImageFilter.h"
#include "itkAddImageFilter.h"
#include <itkNormalizeImageFilter.h>
#include "itkDiscreteGaussianImageFilter.h"
#include <itkGetAverageSliceImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkBinaryMorphologicalOpeningImageFilter.h>
#include <sstream>
#include <itkMultiplyImageFilter.h>
#include "itkMedianImageFunction.h"
#include <itkImageDuplicator.h>
#include <itkGradientAnisotropicDiffusionImageFilter.h>

#ifdef _OPENMP
#include "omp.h"
#endif

#define N 11
#define BOUNDARY_EXPAND 5
std::string SomaInfo[N]={"ID", "centroid_x", "centroid_y", "centroid_z", "soma_volume", "eccentricity", "elongation", "orientation", 
"majorAxisLength", "minorAxisLength", "soma_surface_area"};


SomaExtractor::SomaExtractor()
{
}

SomaExtractor::~SomaExtractor()
{
}

SomaExtractor::ProbImageType::Pointer SomaExtractor::SetInputImage(const char * fileName)
{
	ushortImageReader::Pointer reader = ushortImageReader::New();
	reader->SetFileName (fileName);	

	typedef itk::CastImageFilter<UShortImageType, ProbImageType> CasterType;
	CasterType::Pointer caster = CasterType::New();
	caster->SetInput(reader->GetOutput());
	caster->Update();
	ProbImageType::Pointer inputImage = caster->GetOutput();
	width = inputImage->GetLargestPossibleRegion().GetSize()[0];
	height = inputImage->GetLargestPossibleRegion().GetSize()[1];
	depth = inputImage->GetLargestPossibleRegion().GetSize()[2];
	//std::cout<<width<<"\t"<<height<<"\t"<<depth<<std::endl;
	return inputImage;
}

SomaExtractor::ProbImageType::Pointer SomaExtractor::SetInputImage8bit(const char * fileName)
{
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName (fileName);	

	typedef itk::CastImageFilter<OutputImageType, ProbImageType> CasterType;
	CasterType::Pointer caster = CasterType::New();
	caster->SetInput(reader->GetOutput());
	caster->Update();
	ProbImageType::Pointer inputImage = caster->GetOutput();
	width = inputImage->GetLargestPossibleRegion().GetSize()[0];
	height = inputImage->GetLargestPossibleRegion().GetSize()[1];
	depth = inputImage->GetLargestPossibleRegion().GetSize()[2];
	std::cout<<width<<"\t"<<height<<"\t"<<depth<<std::endl;
	return inputImage;
}

SomaExtractor::ProbImageType::Pointer SomaExtractor::SetInputImageByPortion(const char * fileName)
{
	std::cout<< sizeof(UShortPixel)<<std::endl;
	std::cout<<"Read region:"<<startX<<"\t"<<startY<<"\t"<<startZ<<"\t"<<sizeX<<"\t"<<sizeY<<"\t"<<sizeZ<<std::endl;

	UShortImageType::IndexType index;
	index[0] = startX;
	index[1] = startY;
	index[2] = startZ;

	UShortImageType::SizeType size;
	size[0] = sizeX;
	size[1] = sizeY;
	size[2] = sizeZ;

	UShortImageType::RegionType region;
	region.SetIndex( index);
	region.SetSize( size);

	ushortImageReader::Pointer reader = ushortImageReader::New();

	reader->SetFileName(fileName);
	reader->GetOutput()->SetRequestedRegion(region);

	try
	{
		reader->Update();
	}
	catch ( itk::ExceptionObject &err)
	{
		std::cout << "Error in bin_image_writer: " << err << std::endl; 
	}

	typedef itk::CastImageFilter<UShortImageType, ProbImageType> CasterType;
	CasterType::Pointer caster = CasterType::New();
	caster->SetInput(reader->GetOutput());
	caster->Update();

	ProbImageType::Pointer inputImage = caster->GetOutput();
	width = inputImage->GetLargestPossibleRegion().GetSize()[0];
	height = inputImage->GetLargestPossibleRegion().GetSize()[1];
	depth = inputImage->GetLargestPossibleRegion().GetSize()[2];
	std::cout<<"Region:"<<width<<"\t"<<height<<"\t"<<depth<<std::endl;
	return inputImage;
}

SomaExtractor::OutputImageType::Pointer SomaExtractor::Read8BitImage(const char * fileName)
{
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(fileName);	
	reader->Update();
	OutputImageType::Pointer imagePtr = reader->GetOutput();
	int SM = imagePtr->GetLargestPossibleRegion().GetSize()[0];
	int SN = imagePtr->GetLargestPossibleRegion().GetSize()[1];
	int SZ = imagePtr->GetLargestPossibleRegion().GetSize()[2];
	std::cout<<SM<<"\t"<<SN<<"\t"<<SZ<<std::endl;
	return imagePtr;
}

SomaExtractor::OutputImageType::Pointer SomaExtractor::Read16BitImage(const char * fileName)
{
	ushortImageReader::Pointer reader = ushortImageReader::New();
	reader->SetFileName(fileName);	

	RescaleUshortFilterType::Pointer rescaleFilter = RescaleUshortFilterType::New();
	rescaleFilter->SetInput(reader->GetOutput());
	rescaleFilter->SetOutputMinimum(0);
	rescaleFilter->SetOutputMaximum(255);
	try
	{
		rescaleFilter->Update();
	}
	catch ( itk::ExceptionObject &err)
	{
		std::cout << "Error in Read16BitImage " << err << std::endl; 
	}
	OutputImageType::Pointer imagePtr = rescaleFilter->GetOutput();
	return imagePtr;
}

SomaExtractor::SegmentedImageType::Pointer SomaExtractor::SetInitalContourImage(const char * fileName)
{
	somaImageReader::Pointer reader = somaImageReader::New();
	reader->SetFileName(fileName);	
	try
	{
		reader->Update();
	}
	catch ( itk::ExceptionObject &err)
	{
		std::cout << "Error in bin_image_writer: " << err << std::endl; 
	}
	SegmentedImageType::Pointer imagePtr = reader->GetOutput();
	//OutputImageType::Pointer image = Read8BitImage(fileName);
	//typedef itk::CastImageFilter<OutputImageType, SegmentedImageType> CasterType;
	//   CasterType::Pointer caster = CasterType::New();
	//caster->SetInput(image);
	//caster->Update();
	return imagePtr;
}

SomaExtractor::SegmentedImageType::Pointer SomaExtractor::SetInitalContourImage16bit(const char * fileName)
{
	ushortImageReader::Pointer reader = ushortImageReader::New();
	reader->SetFileName(fileName);	
	try
	{
		reader->Update();
	}
	catch ( itk::ExceptionObject &err)
	{
		std::cout << "Error in bin_image_writer: " << err << std::endl; 
	}
	UShortImageType::Pointer imagePtr = reader->GetOutput();

	typedef itk::CastImageFilter<UShortImageType, SegmentedImageType> CasterType;
	CasterType::Pointer caster = CasterType::New();
	caster->SetInput(imagePtr);
	caster->Update();
	SegmentedImageType::Pointer imagePtr2 = caster->GetOutput();
	return imagePtr2;
}

SomaExtractor::SegmentedImageType::Pointer SomaExtractor::SetInitalContourImageByPortion(const char * fileName)
{
	std::cout<<"Read region:"<<startX<<"\t"<<startY<<"\t"<<startZ<<"\t"<<sizeX<<"\t"<<sizeY<<"\t"<<sizeZ<<std::endl;
	SegmentedImageType::IndexType index;
	//index[0] = startX;
	//index[1] = startY;
	//index[2] = startZ;
	index.Fill(5);
	std::cout<<index<<std::endl;

	SegmentedImageType::SizeType size;
	//size[0] = sizeX;
	//size[1] = sizeY;
	//size[2] = sizeZ;
	size.Fill(5);
	std::cout<<size<<std::endl;

	SegmentedImageType::RegionType region;
	region.SetIndex( index);
	region.SetSize( size);

	ushortImageReader::Pointer reader = ushortImageReader::New();
	reader->SetFileName(fileName);
	//reader->GetOutput()->SetRequestedRegion(region);

	/*try
	{
	reader->Update();
	reader->GetOutput()->Print(std::cout);
	}
	catch(itk::ExceptionObject &err)
	{
	std::cerr << "ExceptionObject caught!" <<std::endl;
	std::cerr << err << std::endl;
	}*/
	reader->Print(std::cout);

	typedef itk::CastImageFilter<UShortImageType, SegmentedImageType> CasterType;
	CasterType::Pointer caster = CasterType::New();
	caster->SetInput(reader->GetOutput());

	typedef itk::ExtractImageFilter< SegmentedImageType, SegmentedImageType> ROIFilterType;
	ROIFilterType::Pointer ROIfilter = ROIFilterType::New();
	ROIfilter->SetExtractionRegion(region);
	ROIfilter->SetInput( caster->GetOutput());
	ROIfilter->SetDirectionCollapseToIdentity();
	try
	{
		ROIfilter->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "ExceptionObject caught!" <<std::endl;
		std::cerr << err << std::endl;
	}

	//caster->Update();
	SegmentedImageType::Pointer imagePtr = ROIfilter->GetOutput();
	return imagePtr;
}

void SomaExtractor::ReadSeedpoints(const char * fileName, std::vector< itk::Index<3> > &seedVec, bool bNucleusTable)
{
	std::cout << "ReadSeedpoints" << std::endl;
	seedVec.clear();
	//int endX = startX + width;
	//int endY = startY + height;
	//int endZ = startZ + depth;
	//std::cout<< "Bounding Area:"<<startX<<"\t"<<startY<<"\t"<<endX<<"\t"<<endY<<std::endl;
	if( false == bNucleusTable)
	{
		vtkSmartPointer<vtkTable> table = ftk::LoadXYZTable( std::string(fileName));
		for( vtkIdType i = 0; i < table->GetNumberOfRows(); i++)
		{
			itk::Index<3> index;
			index[0] = table->GetValue(i,0).ToUnsignedInt();
			index[1] = table->GetValue(i,1).ToUnsignedInt();
			index[2] = table->GetValue(i,2).ToUnsignedInt();
			//if( index[0] >= startX && index[0] <= endX && index[1] >= startY && index[1] <= endY && index[2] >= startZ && index[2] <= endZ)
			//{
			//	itk::Index<3> newIndex;
			//	newIndex[0] = index[0] - startX;
			//	newIndex[1] = index[1] - startY;
			//	newIndex[2] = index[2] - startZ;
			seedVec.push_back(index);
			//}
		}
	}
	else
	{
		vtkSmartPointer<vtkTable> centroidsTable = ftk::LoadTable( std::string(fileName));
		vtkSmartPointer<vtkDoubleArray> centroidX = vtkDoubleArray::SafeDownCast(centroidsTable->GetColumnByName("centroid_x"));
		vtkSmartPointer<vtkDoubleArray> centroidY = vtkDoubleArray::SafeDownCast(centroidsTable->GetColumnByName("centroid_y"));
		vtkSmartPointer<vtkDoubleArray> centroidZ = vtkDoubleArray::SafeDownCast(centroidsTable->GetColumnByName("centroid_z"));

		std::cout<< "Seeds Size:"<<centroidX->GetSize()<<std::endl;
		if( centroidX && centroidY && centroidZ)
		{
			for( vtkIdType i = 0; i < centroidX->GetSize(); i++)
			{
				itk::Index<3> index;
				index[0] = (unsigned int)centroidX->GetValue( i);
				index[1] = (unsigned int)centroidY->GetValue( i);
				index[2] = (unsigned int)centroidZ->GetValue( i);

				//if( index[0] >= startX && index[0] <= endX && index[1] >= startY && index[1] <= endY && index[2] >= startZ && index[2] <= endZ)
				//{
				//itk::Index<3> newIndex;
				//newIndex[0] = index[0] - startX;
				//newIndex[1] = index[1] - startY;
				//newIndex[2] = index[2] - startZ;
				if( index[0] == 0 || index[1] == 0 || index[2] == 0)
				{
					continue;
				}
				std::cout<< index[0]<<"\t"<<index[1]<<"\t"<<index[2]<<std::endl;
				seedVec.push_back(index);
				//}
			}
		}
	}
}

template <class T>
bool SomaExtractor::SetParamValue(std::map<std::string,std::string> &opts, std::string str, T &value, T defVal)
{
	std::map<std::string,std::string>::iterator mi;
	mi = opts.find(str); 
	if( mi != opts.end())
	{ 
		std::istringstream ss((*mi).second); 
		ss>>value;
		return true;
	}
	else
	{ 
		value = defVal; 
		return false;
	}
}

void SomaExtractor::LoadOptions(const char* paramFileName)
{
	std::map<std::string,std::string> opts;
	ifstream fin(paramFileName); 
	std::string name; 
	fin>>name;
	while(fin.good()) 
	{
		char cont[100];	 
		fin.getline(cont, 99);
		opts[name] = std::string(cont);
		fin>>name;
	}
	fin.close();

	SetParamValue<double>(opts, "-sigmoid_alfa", alfa, -2);
	SetParamValue<double>(opts, "-sigmoid_beta", beta, 70);
	SetParamValue<double>(opts, "-open_radius", open_radius, 2);
	SetParamValue<double>(opts, "-fill_radius", fill_radius, 2);
	SetParamValue<int>(opts, "-time_threhold", timethreshold, 5);
	SetParamValue<double>(opts, "-seed_value", seedValue, -3);
	SetParamValue<double>(opts, "-curvature_scaling", curvatureScaling, 0.5);
	SetParamValue<double>(opts, "-advection_scaling", advectScaling, 0.5);
	SetParamValue<double>(opts, "-rmsThreshold", rmsThres, 0.02);
	SetParamValue<int>(opts, "-max_iterations", maxIterations, 800);
	SetParamValue<int>(opts, "-hole_size", holeSize, 10);
	SetParamValue<int>(opts, "-min_object_size", minObjSize, 1000);

	SetParamValue<double>(opts, "-sigma", sigma, 3);
	SetParamValue<unsigned int>(opts, "-numberOfIterations", numberOfIterations, 5);
	SetParamValue<double>(opts, "-noiseLevel", noiseLevel, 2000);
	SetParamValue<double>(opts, "-outlierExpand",outlierExpandValue, 5);
	SetParamValue<unsigned int>(opts, "-smoothIteration",smoothIteration, 5);
	SetParamValue<double>(opts, "-conductance",conductance, 9);
	SetParamValue<int>(opts, "-start_x", startX, 0);
	SetParamValue<int>(opts, "-start_y", startY, 0);
	SetParamValue<int>(opts, "-start_z", startZ, 0);
	SetParamValue<int>(opts, "-size_x", sizeX, 0);
	SetParamValue<int>(opts, "-size_y", sizeY, 0);
	SetParamValue<int>(opts, "-size_z", sizeZ, 0);

	SetParamValue<int>(opts, "-num_bins", num_bins, 128);
	SetParamValue<int>(opts, "-shift", shift, 0);
	SetParamValue<double>(opts, "-scaleMin", scaleMin, 30);
	SetParamValue<double>(opts, "-scaleMax", scaleMax, 35);
	SetParamValue<double>(opts, "-regionXY", regionXY, 30);
	SetParamValue<double>(opts, "-regionZ", regionZ, 20);
	SetParamValue<int>(opts, "-useDistMap", useDistMap, 1);
	SetParamValue<int>(opts, "-sampling_ratio_XY_to_Z", sampling_ratio_XY_to_Z, 2);
	SetParamValue<int>(opts, "-radius", radius, 10);
	SetParamValue<int>(opts, "-rerun", brerun, 0);
}

void SomaExtractor::writeImage(const char* writeFileName, OutputImageType::Pointer image)
{
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(writeFileName);
	writer->SetInput(image);
	try
	{
		writer->Update();
	}
	catch ( itk::ExceptionObject &err)
	{
		std::cout << "Error in bin_image_writer: " << err << std::endl; 
	}
}

void SomaExtractor::writeImage(const char* writeFileName, SegmentedImageType::Pointer image)
{
	//typedef itk::CastImageFilter<SegmentedImageType, OutputImageType> CasterType;
	//   CasterType::Pointer caster = CasterType::New();
	//   caster->SetInput(image);

	somaImageWriter::Pointer soma_image_writer = somaImageWriter::New();
	soma_image_writer->SetFileName(writeFileName);
	soma_image_writer->SetInput(image);

	try
	{
		soma_image_writer->Update();
	}
	catch ( itk::ExceptionObject &err)
	{
		std::cout << "Error in bin_image_writer: " << err << std::endl; 
	}
}

void SomaExtractor::writeImage(const char* writeFileName, ProbImageType::Pointer image, bool bscale)
{
	typedef itk::ImageFileWriter< ProbImageType > ProbImageWriter;
	ProbImageWriter::Pointer writer = ProbImageWriter::New();
	writer->SetFileName( writeFileName);
	writer->SetInput( image);
	writer->Update();
}

void SomaExtractor::writeImage(const char* writeFileName, GradientImageType::Pointer image)
{
	VectorImageWriterType::Pointer writer = VectorImageWriterType::New();
	writer->SetFileName( writeFileName);
	writer->SetInput( image);
	try
	{
		writer->Update();
	}
	catch ( itk::ExceptionObject &err)
	{
		std::cout << "Error in bin_image_writer: " << err << std::endl; 
	}
}

void SomaExtractor::writeCentroids(const char* writeFileName, std::vector< itk::Index<3> > &seedVec)
{
	std::ofstream ofs(writeFileName);
	ofs<<"centroid_x"<<"\t"<<"centroid_y"<<"\t"<<"centroid_z"<<std::endl;
	for( int i = 0; i< seedVec.size(); i++ )
	{
		ofs<< seedVec[i][0]<< "\t"<<seedVec[i][1]<<"\t"<<seedVec[i][2]<<std::endl;
	}
	ofs.close();
}

/// Get gradient edge map
SomaExtractor::ProbImageType::Pointer SomaExtractor::GetEdgePotentialMap(ProbImageType::Pointer inputImage, double sigma)
{
	std::cout<< "Sharpen the images..."<<std::endl;
	LaplacianSharpeningImageFilterType::Pointer laplacianSharpeningImageFilter = LaplacianSharpeningImageFilterType::New();
	laplacianSharpeningImageFilter->SetInput( inputImage);

	//std::cout<< "Smoothing the image with edges preserved."<<std::endl;
	//typedef itk::CurvatureAnisotropicDiffusionImageFilter<ProbImageType, ProbImageType > SmoothFilter;
	//SmoothFilter::Pointer smooth = SmoothFilter::New();
	//smooth->SetInput(inputImage);
	//smooth->SetNumberOfIterations(3);
	//smooth->SetTimeStep(0.0625);
	//smooth->SetConductanceParameter(3);
	//smooth->Update();

	std::cout<< "Gradient Image:"<<sigma<<std::endl;
	GradientFilterType::Pointer gradientMagnitude = GradientFilterType::New();
	gradientMagnitude->SetSigma(sigma);
	gradientMagnitude->SetInput(laplacianSharpeningImageFilter->GetOutput());
	gradientMagnitude->Update();
	ProbImageType::Pointer image = gradientMagnitude->GetOutput();

	//SigmoidImageFilterType::Pointer sigmoidFilter = SigmoidImageFilterType::New();
	//sigmoidFilter->SetInput(image);
	//sigmoidFilter->SetOutputMinimum(0);
	//sigmoidFilter->SetOutputMaximum(255);
	//sigmoidFilter->SetAlpha(alfa);
	//sigmoidFilter->SetBeta(beta);
	//sigmoidFilter->Update();
	//ProbImageType::Pointer floatImage = sigmoidFilter->GetOutput();

	//this->writeImage("GradientImage.tif", image);

	return image;
}

/// Huang Auto Thresholding
SomaExtractor::ProbImageType::Pointer SomaExtractor::EnhanceContrast( ProbImageType::Pointer inputImage, int sliceNum, double alfa, double beta, double &threshold)
{
	ProbImageType::IndexType start;
	start[0] = 0;
	start[1] = 0;
	start[2] = sliceNum;
	//std::cout<< "Slide: "<<start[2]<<std::endl;
	ProbImageType::SizeType size;
	size[0] = inputImage->GetLargestPossibleRegion().GetSize()[0];
	size[1] = inputImage->GetLargestPossibleRegion().GetSize()[1];
	size[2] = 0;
	ProbImageType::RegionType desiredRegion;
	desiredRegion.SetSize(size);
	desiredRegion.SetIndex(start);

	ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
	extractFilter->SetInput(inputImage);
	extractFilter->SetExtractionRegion(desiredRegion);
#if ITK_VERSION_MAJOR >= 4
	extractFilter->SetDirectionCollapseToIdentity(); // This is required.
#endif
	extractFilter->Update();
	ProbImageType::Pointer slice = extractFilter->GetOutput();
	//std::cout<< "Slice Done... "<< std::endl;

	HuangThresholdFilter::Pointer huangThresholdFilter = HuangThresholdFilter::New();
	huangThresholdFilter->SetInput(slice);
	huangThresholdFilter->SetNumberOfHistogramBins( 256);
	huangThresholdFilter->Update();
	threshold = huangThresholdFilter->GetThreshold();

	//std::cout<< "HuangThreshold: "<< threshold<<std::endl;

	BinaryProbThresholdingFilterType::Pointer thresholder = BinaryProbThresholdingFilterType::New();
	thresholder->SetInput(inputImage);
	thresholder->SetUpperThreshold(threshold);                                                           
	thresholder->SetOutsideValue( 255);
	thresholder->SetInsideValue( 0);
	thresholder->Update();
	ProbImageType::Pointer probImage = thresholder->GetOutput();
	//this->writeImage("ThresholdImage.tif", probImage);
	//SigmoidImageFilterType::Pointer sigmoidFilter = SigmoidImageFilterType::New();
	//sigmoidFilter->SetInput(inputImage);
	//sigmoidFilter->SetOutputMinimum(0);
	//sigmoidFilter->SetOutputMaximum(255);
	//sigmoidFilter->SetAlpha(alfa);
	//sigmoidFilter->SetBeta(beta);
	//sigmoidFilter->Update();
	//ProbImageType::Pointer floatImage = sigmoidFilter->GetOutput();

	return probImage;
}

/// Adaptive histogram
SomaExtractor::ProbImageType::Pointer SomaExtractor::EnhanceContrast( ProbImageType::Pointer inputImage, double alfa, double beta, double radius)
{
	AdaptiveHistogramEqualizationImageFilterType::Pointer adaptiveHistogramEqualizationImageFilter = AdaptiveHistogramEqualizationImageFilterType::New();
	adaptiveHistogramEqualizationImageFilter->SetInput(inputImage);
	adaptiveHistogramEqualizationImageFilter->SetRadius(radius);
	adaptiveHistogramEqualizationImageFilter->SetAlpha(alfa);
	adaptiveHistogramEqualizationImageFilter->SetBeta(beta);
	adaptiveHistogramEqualizationImageFilter->Update();
	ProbImageType::Pointer imagePt = adaptiveHistogramEqualizationImageFilter->GetOutput();
	return imagePt;
}

void SomaExtractor::GetSeedpointsInRegion(std::vector< itk::Index<3> > &seedVec, std::vector< itk::Index<3> > &seedInRegion, int startX, int startY, int width, int height)
{
	seedInRegion.clear();
	int endX = startX + width;
	int endY = startY + height;
	for( int i = 0; i < seedVec.size(); i++ )
	{
		itk::Index<3> ind = seedVec[i];
		if( ind[0] >= startX && ind[1] >= startY && ind[0] <= endX && ind[1] <= endY)
		{
			ind[0] -= startX;
			ind[1] -= startY;
			seedInRegion.push_back( ind);
		}
	}
}

// Shape Detection Active Contour Without GVF:
//	double alfa;
//	double beta;
//	int timethreshold; 
//	double curvatureScaling; 
//	double rmsThres;
//	int holeSize;
//	int minObjSize;
SomaExtractor::SegmentedImageType::Pointer SomaExtractor::SegmentSoma( std::vector< itk::Index<3> > &somaCentroids, ProbImageType::Pointer binImagePtr)
{
	//typedef itk::NearestNeighborInterpolateImageFunction< ProbImageType, float>  InterpolatorType;
	//InterpolatorType::Pointer I_Interpolator = InterpolatorType::New();
	//I_Interpolator->SetInputImage(input);

	int SM = binImagePtr->GetLargestPossibleRegion().GetSize()[0];
	int SN = binImagePtr->GetLargestPossibleRegion().GetSize()[1];
	int SZ = binImagePtr->GetLargestPossibleRegion().GetSize()[2];
	std::cout<<SM<<"\t"<<SN<<"\t"<<SZ<<std::endl;

	//move the seed points along z axis
	//   for( int i = 0; i < somaCentroids.size(); i++ )
	//{
	//	ProbImageType::IndexType index;
	//	ProbImageType::IndexType index1;
	//	index[0] = index1[0] = somaCentroids[i][0];
	//	index[1] = index1[1] = somaCentroids[i][1];
	//	index[2] = somaCentroids[i][2];
	//	for( int j = 0; j < SZ; j++ )
	//	{
	//		index1[2] = j;
	//		if( I_Interpolator->EvaluateAtIndex(index1) > I_Interpolator->EvaluateAtIndex(index) )
	//		{
	//			somaCentroids[i][2] = j;
	//			index[2] = j;
	//		}
	//	}
	//}

	//   SigmoidImageFilterType::Pointer sigmoidFilter = SigmoidImageFilterType::New();
	//sigmoidFilter->SetInput(input);
	//sigmoidFilter->SetOutputMinimum(0);
	//sigmoidFilter->SetOutputMaximum(1);
	//sigmoidFilter->SetAlpha(alfa);
	//sigmoidFilter->SetBeta(beta);
	//sigmoidFilter->Update();
	//writeImage("Sigmoid.tif",sigmoidFilter->GetOutput());

	BinaryProbThresholdingFilterType::Pointer binaryFilterPointer = BinaryProbThresholdingFilterType::New();
	binaryFilterPointer->SetInput(binImagePtr);
	binaryFilterPointer->SetLowerThreshold(1);
	binaryFilterPointer->SetInsideValue(1);
	binaryFilterPointer->SetOutsideValue(0);

	FastMarchingFilterType::Pointer fastMarching = FastMarchingFilterType::New();
	NodeContainer::Pointer seeds = NodeContainer::New();
	seeds->Initialize();

	std::cout<< "Seed Size"<< somaCentroids.size()<<std::endl;
	for( int i = 0; i< somaCentroids.size(); i++ )
	{
		ProbImageType::IndexType  seedPosition;
		seedPosition[0] = somaCentroids[i][0];
		seedPosition[1] = somaCentroids[i][1];
		seedPosition[2] = somaCentroids[i][2];
		NodeType node;
		node.SetValue( seedValue );
		node.SetIndex( seedPosition );
		seeds->InsertElement( i, node );
	}

	fastMarching->SetTrialPoints(  seeds);
	fastMarching->SetOutputSize( binImagePtr->GetBufferedRegion().GetSize() );
	fastMarching->SetStoppingValue(  timethreshold);
	fastMarching->SetSpeedConstant( 1.0 );

	std::cout<< "Shape Detection " <<curvatureScaling<<"\t"<<rmsThres<<std::endl;

	ShapeDetectionFilterType::Pointer shapeDetection = ShapeDetectionFilterType::New();

	shapeDetection->SetPropagationScaling(1.0);
	shapeDetection->SetCurvatureScaling(curvatureScaling);
	shapeDetection->SetMaximumRMSError( rmsThres);
	shapeDetection->SetNumberOfIterations( maxIterations);

	shapeDetection->SetInput( fastMarching->GetOutput());
	shapeDetection->SetFeatureImage( binaryFilterPointer->GetOutput());

	//try
	//{
	shapeDetection->Update();
	std::cout << "No. elpased iterations: " << shapeDetection->GetElapsedIterations() << std::endl;
	std::cout << "RMS change: " << shapeDetection->GetRMSChange() << std::endl;

	if( brerun == 1 && shapeDetection->GetElapsedIterations() >= maxIterations)
	{
		curvatureScaling += 0.05;
	}
	while( brerun == 1 && shapeDetection->GetElapsedIterations() >= maxIterations && curvatureScaling < 1)
	{
		std::cout<< "rerun shape detection"<<"\t"<<curvatureScaling<<std::endl;
		shapeDetection->SetCurvatureScaling(curvatureScaling);
		shapeDetection->Update();
		std::cout << "No. elpased iterations: " << shapeDetection->GetElapsedIterations() << std::endl;
		std::cout << "RMS change: " << shapeDetection->GetRMSChange() << std::endl;
		curvatureScaling += 0.05;
	}

	if(curvatureScaling < 1)
	{
		std::cout<<"converge!"<<std::endl;
	}
	else
	{
		std::cout<<"fail to converge!!!!!!!!!!!!!!!!!!!"<<std::endl;
	}
	//}
	//catch( itk::ExceptionObject &err)
	//{
	//	std::cout << "Error in shape detection filter: " << err << std::endl; 
	//	return NULL;
	//}

	std::cout<< "Thresholding..."<<endl;

	BinaryThresholdingFilterType::Pointer thresholder = BinaryThresholdingFilterType::New();
	thresholder->SetLowerThreshold( -10000);
	thresholder->SetUpperThreshold(0.0);                                                           
	thresholder->SetOutsideValue( 0);
	thresholder->SetInsideValue(255);
	thresholder->SetInput( shapeDetection->GetOutput());

	///// morphological closing
	//KernelType ball;
	//KernelType::SizeType ballSize;
	//ballSize[0] = holeSize;
	//ballSize[1] = holeSize;
	//ballSize[2] = 1;
	//ball.SetRadius(ballSize);
	//ball.CreateStructuringElement();
	//CloseFilterType::Pointer closeFilter = CloseFilterType::New();
	//closeFilter->SetInput( thresholder->GetOutput());
	//closeFilter->SetKernel( ball );
	//closeFilter->SetForegroundValue( 255);

	/// Label image
	LabelFilterType::Pointer label = LabelFilterType::New();
	label->SetInput(thresholder->GetOutput());

	RelabelFilterType::Pointer relabel = RelabelFilterType::New();
	relabel->SetInput( label->GetOutput());
	relabel->SetMinimumObjectSize( minObjSize);  

	relabel->Update();
	SegmentedImageType::Pointer somas = relabel->GetOutput();

	LabelGeometryImageFilterType::Pointer labelGeometryImageFilter = LabelGeometryImageFilterType::New();
	labelGeometryImageFilter->SetInput( somas);
	labelGeometryImageFilter->CalculatePixelIndicesOff();
	labelGeometryImageFilter->CalculateOrientedBoundingBoxOff();
	labelGeometryImageFilter->CalculateOrientedLabelRegionsOff();
	labelGeometryImageFilter->CalculateOrientedIntensityRegionsOff();
	labelGeometryImageFilter->Update();
	LabelGeometryImageFilterType::LabelsType allLabels = labelGeometryImageFilter->GetLabels();
	LabelGeometryImageFilterType::LabelsType::iterator allLabelsIt =  allLabels.begin();

	somaCentroids.clear();
	for( allLabelsIt++; allLabelsIt != allLabels.end(); allLabelsIt++ )
	{
		LabelGeometryImageFilterType::LabelPixelType labelValue = *allLabelsIt;
		vtkSmartPointer<vtkVariantArray> row = vtkSmartPointer<vtkVariantArray>::New();
		LabelGeometryImageFilterType::LabelPointType point = labelGeometryImageFilter->GetCentroid(labelValue);
		itk::Index<3> index;
		index[0] = point[0];
		index[1] = point[1];
		index[2] = point[2];
		somaCentroids.push_back(index);
	}
	return somas;
}

// Shape Detection Active Contour Without GVF:
//	double alfa;
//	double beta;
//	int timethreshold; 
//	double curvatureScaling; 
//	double rmsThres;
//	int holeSize;
//	int minObjSize;

SomaExtractor::SegmentedImageType::Pointer SomaExtractor::SegmentHeart( const char *imageName, const char *fileName, ProbImageType::Pointer inputImage, vnl_vector<int> &seperator, vnl_vector<double> &curvature)
{
	int SM = inputImage->GetLargestPossibleRegion().GetSize()[0];
	int SN = inputImage->GetLargestPossibleRegion().GetSize()[1];
	int SZ = inputImage->GetLargestPossibleRegion().GetSize()[2];
	std::cout<<SM<<"\t"<<SN<<"\t"<<SZ<<std::endl;

	int count = seperator.size();

	std::string inputName = std::string(imageName);
	inputName.erase(inputName.length()-4,inputName.length());
	inputName.append("_smooth.tif");
	typedef itk::CurvatureAnisotropicDiffusionImageFilter<ProbImageType, ProbImageType> SmoothingFilterType;
	SmoothingFilterType::Pointer smoothing = SmoothingFilterType::New();
	smoothing->SetInput(inputImage);
	smoothing->SetTimeStep( 0.0625);
	smoothing->SetNumberOfIterations( smoothIteration);
	smoothing->SetConductanceParameter( conductance);
	
	typedef itk::CastImageFilter<ProbImageType, OutputImageType> CasterType;
	CasterType::Pointer caster = CasterType::New();
	caster->SetInput(smoothing->GetOutput());
	caster->Update();
	OutputImageType::Pointer smoothImage = caster->GetOutput();
	writeImage(inputName.c_str(), smoothImage);

	BinaryProbThresholdingFilterType::Pointer thresholdFilter = BinaryProbThresholdingFilterType::New();
	thresholdFilter->SetLowerThreshold( 0);
	thresholdFilter->SetUpperThreshold(beta);                                                           
	thresholdFilter->SetOutsideValue( 0);
	thresholdFilter->SetInsideValue(1);
	thresholdFilter->SetInput( smoothing->GetOutput());
	//thresholdFilter->Update();
	//writeImage("binary.nrrd", thresholdFilter->GetOutput());

	typedef itk::BinaryBallStructuringElement< ProbImageType::PixelType, Dim> KType;
	typedef itk::BinaryMorphologicalOpeningImageFilter< ProbImageType, ProbImageType, KType > OpenMorphFilterType;
	OpenMorphFilterType::Pointer openFilter = OpenMorphFilterType::New();
	KType ball;
	KType::SizeType ballSize;
	ballSize.Fill( open_radius);
	ball.SetRadius(ballSize);
	ball.CreateStructuringElement();
	openFilter->SetInput( thresholdFilter->GetOutput());
	openFilter->SetKernel( ball );
	openFilter->SetForegroundValue( 1);
	typedef itk::VotingBinaryHoleFillingImageFilter< ProbImageType, ProbImageType> HoleFillFilterType;
	HoleFillFilterType::Pointer holeFillFilter = HoleFillFilterType::New();
	holeFillFilter->SetInput(openFilter->GetOutput());
	ProbImageType::SizeType indexRadius;
	indexRadius[0] = fill_radius; 
	indexRadius[1] = fill_radius;
	indexRadius[2] = fill_radius; 
	holeFillFilter->SetRadius( indexRadius);
	holeFillFilter->SetBackgroundValue( 0);
	holeFillFilter->SetForegroundValue( 1);
	holeFillFilter->Update();
	//writeImage("openandfill.nrrd", holeFillFilter->GetOutput());
	//typedef itk::SigmoidImageFilter<ProbImageType, ProbImageType> SigmoidFilterType;
	//SigmoidFilterType::Pointer sigmoid = SigmoidFilterType::New();
	//sigmoid->SetAlpha( alfa);
	//sigmoid->SetBeta(  beta);
	//sigmoid->SetOutputMinimum( 0.0);
	//sigmoid->SetOutputMaximum( 1.0);
	//sigmoid->SetInput( smoothing->GetOutput());
	//sigmoid->Update();
	ProbImageType::Pointer sigmoidImage = holeFillFilter->GetOutput();
	//writeImage("sigmoid.nrrd", sigmoidImage);

	SegmentedImageType::Pointer segImage = SegmentedImageType::New();
	segImage->SetLargestPossibleRegion( inputImage->GetLargestPossibleRegion());
	segImage->SetBufferedRegion( inputImage->GetLargestPossibleRegion());
	segImage->SetRequestedRegion( inputImage->GetLargestPossibleRegion());
	segImage->Allocate();
	segImage->FillBuffer(0);
	ProbImageType::PixelType * inputBuffer = sigmoidImage->GetBufferPointer();
	SegmentedImageType::PixelType * segBuffer = segImage->GetBufferPointer();

	for(int n = 0; n < count; n++)
	{
		/// seperate the images into N small stacks
		ProbImageType::Pointer sepStack = ProbImageType::New();
		ProbImageType::SizeType size;
		size[0] = SM;
		size[1] = SN;
		int sliceIndex = 0;
		if( n == 0)
		{
			size[2] = seperator[0] + 2;
		}
		else
		{
			sliceIndex = seperator[n-1];
			size[2] = seperator[n] - seperator[n-1] + 2;
		}
		ProbImageType::IndexType start;
		start.Fill(0);
		ProbImageType::RegionType region;
		region.SetIndex( start);
		region.SetSize( size);

		sepStack->SetLargestPossibleRegion( region);
		sepStack->SetBufferedRegion( region);
		sepStack->SetRequestedRegion( region);
		sepStack->Allocate();
		sepStack->FillBuffer(0);
		ProbImageType::PixelType * sepBuffer = sepStack->GetBufferPointer();

		for( int kk = 0; kk < size[2]; ++kk)
		{
			for( int jj = 0; jj< size[1]; ++jj)
			{
				for( int ii = 0; ii < size[0]; ++ii)
				{
					itk::Index<1> offset;
					itk::Index<1> offset2;
					offset[0] = (sliceIndex + kk) * size[0] * size[1] + jj * size[0] + ii;
					offset2[0] = kk * size[0] * size[1] + jj * size[0] + ii;
					sepBuffer[offset2[0]] = inputBuffer[offset[0]];
				}
			}
		}

		/// read the centroids
		std::stringstream ss;
		ss << (n+1);
		//std::string oriImageName = ss.str();
		//oriImageName.append("sm.nrrd");
		//std::string segImageName = ss.str();
		//segImageName.append("seg.nrrd");
		//writeImage(oriImageName.c_str(), sepStack);

		std::vector< itk::Index<3> > somaCentroids;
		std::string seedFileName = std::string(fileName);
		seedFileName += ss.str();
		seedFileName.append(".txt");
		std::cout<< "Seed File "<<seedFileName<<std::endl;
		ReadSeedpoints(seedFileName.c_str(), somaCentroids, true);
		FastMarchingFilterType::Pointer fastMarching = FastMarchingFilterType::New();
		NodeContainer::Pointer seeds = NodeContainer::New();
		seeds->Initialize();

		std::cout<< "Seed Size"<< somaCentroids.size()<<std::endl;
		for( int i = 0; i< somaCentroids.size(); i++ )
		{
			ProbImageType::IndexType seedPosition;
			seedPosition[0] = somaCentroids[i][0];
			seedPosition[1] = somaCentroids[i][1];
			seedPosition[2] = somaCentroids[i][2];
			NodeType node;
			node.SetValue( seedValue );
			node.SetIndex( seedPosition );
			seeds->InsertElement( i, node );
		}

		fastMarching->SetTrialPoints(  seeds);
		fastMarching->SetOutputSize( sepStack->GetBufferedRegion().GetSize() );
		fastMarching->SetStoppingValue(  timethreshold);
		fastMarching->SetSpeedConstant( 1.0 );

		std::cout<< n<<": Shape Detection " <<curvature[n]<<"\t"<<rmsThres<<std::endl;

		ShapeDetectionFilterType::Pointer shapeDetection = ShapeDetectionFilterType::New();

		shapeDetection->SetPropagationScaling(1.0);
		shapeDetection->SetCurvatureScaling(curvature[n]);
		shapeDetection->SetMaximumRMSError( rmsThres);
		shapeDetection->SetNumberOfIterations( maxIterations);

		shapeDetection->SetInput( fastMarching->GetOutput());
		shapeDetection->SetFeatureImage( sepStack);

		//try
		//{
		shapeDetection->Update();
		std::cout << "No. elpased iterations: " << shapeDetection->GetElapsedIterations() << std::endl;
		std::cout << "RMS change: " << shapeDetection->GetRMSChange() << std::endl;

		if( brerun == 1 && shapeDetection->GetElapsedIterations() >= maxIterations)
		{
			curvatureScaling += 0.05;
		}
		while( brerun == 1 && shapeDetection->GetElapsedIterations() >= maxIterations && curvatureScaling < 1)
		{
			std::cout<< "rerun shape detection"<<"\t"<<curvatureScaling<<std::endl;
			shapeDetection->SetCurvatureScaling(curvatureScaling);
			shapeDetection->Update();
			std::cout << "No. elpased iterations: " << shapeDetection->GetElapsedIterations() << std::endl;
			std::cout << "RMS change: " << shapeDetection->GetRMSChange() << std::endl;
			curvatureScaling += 0.05;
		}

		if(curvatureScaling < 1)
		{
			std::cout<<"converge!"<<std::endl;
		}
		else
		{
			std::cout<<"fail to converge!!!!!!!!!!!!!!!!!!!"<<std::endl;
		}
		//}
		//catch( itk::ExceptionObject &err)
		//{
		//	std::cout << "Error in shape detection filter: " << err << std::endl; 
		//	return NULL;
		//}

		std::cout<< "Thresholding..."<<endl;

		BinaryProbThresholdingFilterType::Pointer thresholder = BinaryProbThresholdingFilterType::New();
		thresholder->SetLowerThreshold( -10000);
		thresholder->SetUpperThreshold(0.0);                                                           
		thresholder->SetOutsideValue( 0);
		thresholder->SetInsideValue(255);
		thresholder->SetInput( shapeDetection->GetOutput());
		thresholder->Update();
		ProbImageType::Pointer segContour = thresholder->GetOutput();
		ProbImageType::PixelType * segContourBuffer = segContour->GetBufferPointer();

		for( int kk = 1; kk < size[2] - 1; ++kk)
		{
			for( int jj=0; jj < size[1]; ++jj)
			{
				for( int ii = 0; ii < size[0]; ++ii)
				{
					itk::Index<1> offset;
					itk::Index<1> offset2;
					offset[0] = kk * size[0] * size[1] + jj * size[0] + ii;
					offset2[0] = (sliceIndex + kk) * size[0] * size[1] + jj * size[0] + ii;
					segBuffer[offset2[0]] = (SegmentedImageType::PixelType)segContourBuffer[offset[0]];
				}
			}
		}
		//writeImage(segImageName.c_str(), segContour);
	}

	/// Label image

	LabelFilterType::Pointer label = LabelFilterType::New();
	label->SetInput(segImage);
	label->Update();

	//RelabelFilterType::Pointer relabel = RelabelFilterType::New();
	//relabel->SetInput( label->GetOutput());
	//relabel->SetMinimumObjectSize( minObjSize);  

	//relabel->Update();
	SegmentedImageType::Pointer somas = label->GetOutput();

	return somas;
}

/// Generate Expanded Inital Contours within the boundary of labeled object by Daniel Distance Map
SomaExtractor::ProbImageType::Pointer SomaExtractor::GetInitalContourByDistanceMap(SegmentedImageType::Pointer labelImage, double outlierExpand)
{
	int SX = labelImage->GetLargestPossibleRegion().GetSize()[0];
	int SY = labelImage->GetLargestPossibleRegion().GetSize()[1];
	int SZ = labelImage->GetLargestPossibleRegion().GetSize()[2];
	std::cout<< SX<<"\t"<<SY<<"\t"<<SZ<<std::endl;

	SegThresholdingFilterType::Pointer thresholder = SegThresholdingFilterType::New();
	thresholder->SetInput(labelImage);
	thresholder->SetLowerThreshold(1); 
	thresholder->SetUpperThreshold(100000); 
	thresholder->SetOutsideValue( 0);
	thresholder->SetInsideValue(255);
	std::cout<< "SegThresholdingFilterType"<<std::endl;

	/////generate bounding box
	//LabelFilterType::Pointer label = LabelFilterType::New();
	//label->SetInput(thresholder->GetOutput());

	//LabelGeometryImageFilterType::Pointer labelGeometryImageFilter = LabelGeometryImageFilterType::New();
	//labelGeometryImageFilter->CalculateOrientedBoundingBoxOff();
	//labelGeometryImageFilter->CalculateOrientedIntensityRegionsOff();
	//labelGeometryImageFilter->CalculateOrientedLabelRegionsOff();
	//labelGeometryImageFilter->CalculatePixelIndicesOff();
	//labelGeometryImageFilter->SetInput(label->GetOutput());
	//labelGeometryImageFilter->Update();

	//LabelGeometryImageFilterType::LabelsType allLabels = labelGeometryImageFilter->GetLabels();
	//LabelGeometryImageFilterType::LabelsType::iterator allLabelsIt =  allLabels.begin();
	SegmentedImageType::IndexType start;
	SegmentedImageType::IndexType end;
	start[0] = 0;
	start[1] = 0;
	start[2] = 0;
	end[0] = SX;
	end[1] = SY;
	end[2] = SZ;
	//end[0] = end[1] = end[2] = 0;

	//for( allLabelsIt++; allLabelsIt != allLabels.end(); allLabelsIt++ )
	//   {
	//	LabelGeometryImageFilterType::LabelPixelType labelValue = *allLabelsIt;
	//	LabelGeometryImageFilterType::BoundingBoxType boundingBox = labelGeometryImageFilter->GetBoundingBox(labelValue);
	//	start[0] = boundingBox[0] < start[0] ? boundingBox[0] : start[0];
	//	start[1] = boundingBox[2] < start[1] ? boundingBox[2] : start[1];
	//	start[2] = boundingBox[4] < start[2] ? boundingBox[4] : start[2];

	//	end[0] = boundingBox[1] > end[0] ? boundingBox[1] : end[0];
	//	end[1] = boundingBox[2] > end[1] ? boundingBox[3] : end[1];
	//	end[2] = boundingBox[5] > end[2] ? boundingBox[5] : end[2];
	//}

	//start[0] -= (outlierExpand + BOUNDARY_EXPAND);
	//start[1] -= (outlierExpand + BOUNDARY_EXPAND);
	//start[2] -= (outlierExpand + BOUNDARY_EXPAND);
	//end[0] += (outlierExpand + BOUNDARY_EXPAND);
	//end[1] += (outlierExpand + BOUNDARY_EXPAND);
	//end[2] += (outlierExpand + BOUNDARY_EXPAND);
	//CheckBoundary(start, end, SX, SY, SZ);
	//std::cout<< "bounding box"<<std::endl;

	/// Generate Distance Map within the bounding box
	SegmentedImageType::SizeType size;
	size[0] = end[0] - start[0];
	size[1] = end[1] - start[1];
	size[2] = end[2] - start[2];

	SegmentedImageType::RegionType desiredRegion;
	desiredRegion.SetSize(size);
	desiredRegion.SetIndex(start);

	RegionOfInterestFilter::Pointer regionFilter = RegionOfInterestFilter::New();
	regionFilter->SetInput(thresholder->GetOutput());
	regionFilter->SetRegionOfInterest(desiredRegion);

	typedef itk::CastImageFilter<SegmentedImageType, ProbImageType> CasterType;
	CasterType::Pointer caster = CasterType::New();
	caster->SetInput(regionFilter->GetOutput());

	//DanielssonDistanceMapFilterType::Pointer distanceMapFilter = DanielssonDistanceMapFilterType::New();
	MaurerDistanceMapFilterType::Pointer distanceMapFilter = MaurerDistanceMapFilterType::New();
	distanceMapFilter->SetInput(caster->GetOutput());
	SubtractImageFilterType::Pointer substractImageFilter = SubtractImageFilterType::New();   // expand the contour
	substractImageFilter->SetInput1( distanceMapFilter->GetOutput());
	substractImageFilter->SetConstant2( outlierExpand);

	try
	{
		substractImageFilter->Update();
	}
	catch(itk::ExceptionObject exp)
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << exp << std::endl;
	}

	/// Write the distance map back to the original size image filled with the maximum value of the distance map
	ProbImageType::Pointer initialContourRegionImage = substractImageFilter->GetOutput();
	ProbImageType::Pointer initialContourImage = ProbImageType::New();
	initialContourImage->SetRegions(labelImage->GetLargestPossibleRegion());
	initialContourImage->Allocate();

	ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New();
	imageCalculatorFilter->SetImage(initialContourRegionImage);
	imageCalculatorFilter->Compute();
	ProbImageType::PixelType maxPixel = imageCalculatorFilter->GetMaximum();
	initialContourImage->FillBuffer(maxPixel);

	ProbConstIteratorType inputIt( initialContourRegionImage, initialContourRegionImage->GetLargestPossibleRegion());
	ProbIteratorType outputIt( initialContourImage, desiredRegion);
	inputIt.GoToBegin();
	outputIt.GoToBegin();
	while( !inputIt.IsAtEnd())
	{
		outputIt.Set( inputIt.Get());
		++inputIt;
		++outputIt;
	}
	return initialContourImage;
}


void SomaExtractor::CheckBoundary(SegmentedImageType::IndexType &start, SegmentedImageType::IndexType &end, int SX, int SY, int SZ)
{
	start[0] = start[0] > 0 ? start[0] : 0;
	start[1] = start[1] > 0 ? start[1] : 0;
	start[2] = start[2] > 0 ? start[2] : 0;
	end[0] = end[0] < SX ? end[0] : SX;
	end[1] = end[1] < SY ? end[1] : SY;
	end[2] = end[2] < SZ ? end[2] : SZ;
}

// GVF Active Contour with advection field and curvature constraint:
// edge image:
// double sigma
// GVF:
// unsigned int numberOfIterations; 
// double noiseLevel;
// double outlierExpandValue to expand the initial contour
// activecontour:
// double curvatureScaling; 
// double advectScaling
// double rmsThres;
// int minObjSize;
SomaExtractor::SegmentedImageType::Pointer SomaExtractor::SegmentSomaUsingGradient( ProbImageType::Pointer input, SegmentedImageType::Pointer initialContour, std::vector< itk::Index<3> > &somaCentroids) 																	  
{
	int SM = input->GetLargestPossibleRegion().GetSize()[0];
	int SN = input->GetLargestPossibleRegion().GetSize()[1];
	int SZ = input->GetLargestPossibleRegion().GetSize()[2];
	std::cout<<SM<<"\t"<<SN<<"\t"<<SZ<<std::endl;

	FastMarchingFilterType::Pointer fastMarching = FastMarchingFilterType::New();
	NodeContainer::Pointer seeds = NodeContainer::New();
	seeds->Initialize();
	std::cout<< "Seed Size"<< somaCentroids.size()<<std::endl;
	for( int i = 0; i< somaCentroids.size(); i++ )
	{
		ProbImageType::IndexType  seedPosition;
		seedPosition[0] = somaCentroids[i][0];
		seedPosition[1] = somaCentroids[i][1];
		seedPosition[2] = somaCentroids[i][2];
		NodeType node;
		node.SetValue( seedValue);
		node.SetIndex( seedPosition);
		seeds->InsertElement( i, node);
	}

	fastMarching->SetTrialPoints(  seeds);
	fastMarching->SetOutputSize( input->GetBufferedRegion().GetSize() );
	fastMarching->SetStoppingValue(  timethreshold);
	fastMarching->SetSpeedConstant( 1.0 );

	typedef itk::CurvatureAnisotropicDiffusionImageFilter<ProbImageType, ProbImageType> SmoothingFilterType;
	typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<ProbImageType, ProbImageType>  GradientFilterType;
	typedef itk::SigmoidImageFilter<ProbImageType, ProbImageType> SigmoidFilterType;

	SmoothingFilterType::Pointer smoothing = SmoothingFilterType::New();
	smoothing->SetInput(input);
	smoothing->SetTimeStep( 0.0625);
	smoothing->SetNumberOfIterations( smoothIteration);
	smoothing->SetConductanceParameter( conductance);
	smoothing->Update();
	writeImage("smooth.nrrd", smoothing->GetOutput());

	GradientFilterType::Pointer  gradientMagnitude = GradientFilterType::New();
	gradientMagnitude->SetInput( smoothing->GetOutput());
	gradientMagnitude->SetSigma( sigma);
	gradientMagnitude->Update();
	writeImage("gradient.nrrd", gradientMagnitude->GetOutput());

	SigmoidFilterType::Pointer sigmoid = SigmoidFilterType::New();
	sigmoid->SetAlpha( alfa);
	sigmoid->SetBeta(  beta);
	sigmoid->SetOutputMinimum( 0.0);
	sigmoid->SetOutputMaximum( 1.0);
	sigmoid->SetInput( gradientMagnitude->GetOutput());
	sigmoid->Update();
	writeImage("sigmoid.nrrd", sigmoid->GetOutput());
	return NULL;

	std::cout<<"Active Contour: "<<curvatureScaling<<"\t"<<advectScaling<<"\t"<<rmsThres<<std::endl;
	GeodesicActiveContourFilterType::Pointer GVF_snake = GeodesicActiveContourFilterType::New();
	GVF_snake->SetInput(fastMarching->GetOutput());
	GVF_snake->SetFeatureImage(sigmoid->GetOutput());
	GVF_snake->SetPropagationScaling( 1.0);
	GVF_snake->SetCurvatureScaling( curvatureScaling);
	GVF_snake->SetAdvectionScaling( advectScaling);
	GVF_snake->SetMaximumRMSError( rmsThres);
	GVF_snake->SetNumberOfIterations( maxIterations);

	try
	{
		GVF_snake->Update();
	}
	catch(...)
	{
		cout << "An exception occurred" << std::endl;
		exit(1);
	}

	std::cout << "No. elpased iterations: " << GVF_snake->GetElapsedIterations() << std::endl;
	std::cout << "RMS Error: " << GVF_snake->GetRMSChange() << std::endl;
	std::cout<< "Thresholding..."<<endl;

	BinaryThresholdingFilterType::Pointer thresholder = BinaryThresholdingFilterType::New();
	thresholder->SetLowerThreshold( -10000);
	thresholder->SetUpperThreshold(0.0);                                                           
	thresholder->SetOutsideValue( 0);
	thresholder->SetInsideValue(255);
	thresholder->SetInput( GVF_snake->GetOutput());

	/// Label image, recaculate centroids

	LabelFilterType::Pointer label = LabelFilterType::New();
	label->SetInput(thresholder->GetOutput());

	RelabelFilterType::Pointer relabel = RelabelFilterType::New();
	relabel->SetInput( label->GetOutput());
	relabel->SetMinimumObjectSize( minObjSize);  

	relabel->Update();
	SegmentedImageType::Pointer somas = relabel->GetOutput();

	LabelGeometryImageFilterType::Pointer labelGeometryImageFilter = LabelGeometryImageFilterType::New();
	labelGeometryImageFilter->SetInput( somas);
	labelGeometryImageFilter->CalculatePixelIndicesOff();
	labelGeometryImageFilter->CalculateOrientedBoundingBoxOff();
	labelGeometryImageFilter->CalculateOrientedLabelRegionsOff();
	labelGeometryImageFilter->CalculateOrientedIntensityRegionsOff();
	labelGeometryImageFilter->Update();
	LabelGeometryImageFilterType::LabelsType allLabels = labelGeometryImageFilter->GetLabels();
	LabelGeometryImageFilterType::LabelsType::iterator allLabelsIt =  allLabels.begin();

	somaCentroids.clear();
	for( allLabelsIt++; allLabelsIt != allLabels.end(); allLabelsIt++ )
	{
		LabelGeometryImageFilterType::LabelPixelType labelValue = *allLabelsIt;
		vtkSmartPointer<vtkVariantArray> row = vtkSmartPointer<vtkVariantArray>::New();
		LabelGeometryImageFilterType::LabelPointType point = labelGeometryImageFilter->GetCentroid(labelValue);
		itk::Index<3> index;
		index[0] = point[0];
		index[1] = point[1];
		index[2] = point[2];
		somaCentroids.push_back(index);
	}
	return somas;
}

//// write the x component of gvf
//GradientImageType::Pointer gradientImage = gvfFilter->GetOutput();
//GradientImageType::RegionType region = gradientImage->GetLargestPossibleRegion();
//GradientImageType::SizeType imageSize = region.GetSize();
//
//ProbImageType::Pointer ximage = ProbImageType::New();
//ximage->SetRegions(region);
//ximage->Allocate();
//ximage->FillBuffer(0);

//for ( unsigned int z = 0; z < imageSize[2]; z++)
//{
//	for(unsigned int j = 0; j < imageSize[1]; j++)
//	{
//		for(unsigned int i = 0; i < imageSize[0]; i++)
//		{
//			GradientImageType::IndexType index;
//			index[0] = i;
//			index[1] = j;
//			index[2] = z;

//			GradientImageType::PixelType pixel = gradientImage->GetPixel(index);
//			ProbImageType::IndexType index2;
//			index2[0] = i;
//			index2[1] = j;
//			index2[2] = z;
//			ximage->SetPixel(index2, pixel[0]);
//		 }
//	}
//}
//writeImage("GVFX.nrrd", ximage);

/// For soma surface feature
void SomaExtractor::SomaBoundaryScan(SegmentedImageType::Pointer labelImage, std::map< TLPixel, int> &LtoIMap, std::vector< int> &boundaryPixSize)
{
	if(!labelImage) return;

	typedef itk::ConstantBoundaryCondition< SegmentedImageType> boundaryConditionType;
	typedef itk::ConstNeighborhoodIterator< SegmentedImageType, boundaryConditionType > NeighborhoodIteratorType;

	// The offsets for the neighboring pixels for 4-connectivity
	bool cyto_image = false;
	unsigned int dim;
	if ( cyto_image == true)
		dim = 4 * Dim;
	else
		dim = 2 * Dim;

	std::vector< NeighborhoodIteratorType::OffsetType > offsets(dim);
	for ( unsigned int i=0; i < dim; ++i)
	{
		offsets[i].Fill(0);
	}
	unsigned int p = 0, o = 0;
	while ( p < Dim)
	{
		if ( cyto_image == true)
			offsets[o++][p] = -2;
		offsets[o++][p] = -1;
		offsets[o++][p] = 1;
		if ( cyto_image == true)
			offsets[o++][p] = 2;
		p++;
	}

	NeighborhoodIteratorType::RadiusType radius;
	if ( cyto_image == true)
		radius.Fill(2);
	else
		radius.Fill(1);	

	boundaryPixSize.clear();
	int numLabels = (int)LtoIMap.size();
	boundaryPixSize.resize( numLabels );
	for( int i = 0; i < boundaryPixSize.size(); i++)
	{
		boundaryPixSize[i] = 0;
	}

	NeighborhoodIteratorType it( radius, labelImage, labelImage->GetRequestedRegion() );
	for ( it.GoToBegin(); !it.IsAtEnd(); ++it ) 
	{
		TLPixel v = it.GetCenterPixel();  // in the mask

		if ( v <= 0 ) continue;

		bool allSame = true;
		for (unsigned int i=0; i<dim; ++i)
		{
			TLPixel p = it.GetPixel( offsets[i] );

			if ( v != p )
			{
				allSame = false;
				boundaryPixSize[LtoIMap[v]] += 1;
				break;
			}
		}
	}
}

vtkSmartPointer<vtkTable> SomaExtractor::ComputeSomaFeatures(SegmentedImageType::Pointer inputImage)
{
	LabelGeometryImageFilterType::Pointer labelGeometryImageFilter = LabelGeometryImageFilterType::New();
	labelGeometryImageFilter->SetInput( inputImage);

	// These generate optional outputs.
	//labelGeometryImageFilter->CalculatePixelIndicesOn();
	//labelGeometryImageFilter->CalculateOrientedBoundingBoxOn();
	//labelGeometryImageFilter->CalculateOrientedLabelRegionsOn();

	labelGeometryImageFilter->Update();

	/// mapping from label to continuous index
	std::map< TLPixel, int> LtoIMap;
	LabelGeometryImageFilterType::LabelsType allLabels = labelGeometryImageFilter->GetLabels();
	for (int i = 0; i < allLabels.size(); ++i)
	{
		TLPixel l = allLabels.at(i);
		LtoIMap[l] = i;
	}

	std::vector< int> boundaryPixSize; //boundary pixels size for each label
	SomaBoundaryScan(inputImage, LtoIMap, boundaryPixSize);

	/// Calculating geometry features and build a vtkTable
	LabelGeometryImageFilterType::LabelsType::iterator allLabelsIt =  allLabels.begin();
	std::cout << "Calculating geometry labels, number of labels: " << labelGeometryImageFilter->GetNumberOfLabels() << std::endl;

	vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();	

	for( int i = 0; i < N; i++)
	{
		vtkSmartPointer<vtkVariantArray> column = vtkSmartPointer<vtkVariantArray>::New();
		column->SetName(SomaInfo[i].c_str());
		table->AddColumn(column);
	}

	// remove the first labeled object which is the background
	int count = 0;
	for( allLabelsIt++; allLabelsIt != allLabels.end(); allLabelsIt++ )
	{
		LabelGeometryImageFilterType::LabelPixelType labelValue = *allLabelsIt;
		vtkSmartPointer<vtkVariantArray> row = vtkSmartPointer<vtkVariantArray>::New();
		LabelGeometryImageFilterType::LabelPointType point = labelGeometryImageFilter->GetCentroid(labelValue);
		unsigned int volume = labelGeometryImageFilter->GetVolume(labelValue);
		row->InsertNextValue(vtkVariant(count++));

		row->InsertNextValue(vtkVariant(point[0]));
		row->InsertNextValue(vtkVariant(point[1]));
		row->InsertNextValue(vtkVariant(point[2]));

		row->InsertNextValue(vtkVariant(volume));
		row->InsertNextValue(vtkVariant(labelGeometryImageFilter->GetEccentricity(labelValue)));
		row->InsertNextValue(vtkVariant(labelGeometryImageFilter->GetElongation(labelValue)));
		row->InsertNextValue(vtkVariant(labelGeometryImageFilter->GetOrientation(labelValue)));
		row->InsertNextValue(vtkVariant(labelGeometryImageFilter->GetMajorAxisLength(labelValue)));
		row->InsertNextValue(vtkVariant(labelGeometryImageFilter->GetMinorAxisLength(labelValue)));
		row->InsertNextValue(vtkVariant(boundaryPixSize[LtoIMap[labelValue] ]));
		table->InsertNextRow(row);
	}

	return table;
}

vtkSmartPointer<vtkTable> SomaExtractor::ComputeHeartFeatures(SegmentedImageType::Pointer inputImage)
{
	LabelGeometryImageFilterType::Pointer labelGeometryImageFilter = LabelGeometryImageFilterType::New();
	labelGeometryImageFilter->SetInput( inputImage);
	labelGeometryImageFilter->Update();
	LabelGeometryImageFilterType::LabelsType allLabels = labelGeometryImageFilter->GetLabels();
	/// Calculating geometry features and build a vtkTable
	LabelGeometryImageFilterType::LabelsType::iterator allLabelsIt =  allLabels.begin();
	std::cout << "Calculating geometry labels, number of labels: " << labelGeometryImageFilter->GetNumberOfLabels() << std::endl;

	vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();	
	vtkSmartPointer<vtkVariantArray> column1 = vtkSmartPointer<vtkVariantArray>::New();
	column1->SetName("Id");
	table->AddColumn(column1);
	vtkSmartPointer<vtkVariantArray> column2 = vtkSmartPointer<vtkVariantArray>::New();
	column2->SetName("Volume");
	table->AddColumn(column2);

	// remove the first labeled object which is the background
	int count = 0;
	unsigned int volume = 0;
	for( allLabelsIt++; allLabelsIt != allLabels.end(); allLabelsIt++ )
	{
		LabelGeometryImageFilterType::LabelPixelType labelValue = *allLabelsIt;
		volume += labelGeometryImageFilter->GetVolume(labelValue);
	}
	vtkSmartPointer<vtkVariantArray> row = vtkSmartPointer<vtkVariantArray>::New();
	row->InsertNextValue(vtkVariant(0));
	row->InsertNextValue(vtkVariant(volume));
	table->InsertNextRow(row);
	return table;
}

/// Debris Association
void SomaExtractor::GetDebrisCentroids( OutputImageType::Pointer inputImage, std::vector< itk::Index<3> > &debrisSeeds)
{
	debrisSeeds.clear();
	int SM = inputImage->GetLargestPossibleRegion().GetSize()[0];
	int SN = inputImage->GetLargestPossibleRegion().GetSize()[1];
	int SZ = inputImage->GetLargestPossibleRegion().GetSize()[2];
	std::cout<<SM<<"\t"<<SN<<"\t"<<SZ<<std::endl;
	for(int x=0; x<SM; x++) // for all Columns
	{
		for(int y=0; y<SN; y++) // for all Rows
		{
			for(int z = 0; z < SZ; z++)
			{
				OutputImageType::IndexType ind;
				ind[0] = x;
				ind[1] = y;
				ind[2] = z;
				if(inputImage->GetPixel(ind)> 0)
				{
					itk::Index<3> index;
					index[0] = x;
					index[1] = y;
					index[2] = z;
					debrisSeeds.push_back(index);
				}
			}
		}
	}
}

void SomaExtractor::AssociateDebris(OutputImageType::Pointer inputImage, std::vector< itk::Index<3> > &somaCentroids, std::vector< itk::Index<3> > &debrisSeeds)
{
	std::vector< int> debrisIntensity;
	debrisIntensity.resize( somaCentroids.size());
	for( int i = 0; i < debrisIntensity.size(); i++)
	{
		debrisIntensity[i] = 0;
	}

#pragma omp parallel for
	for( int i = 0; i < debrisSeeds.size(); i++)
	{
		int x = debrisSeeds[i][0];
		int y = debrisSeeds[i][1];
		int z = debrisSeeds[i][2];
		double minDistance = 10e5;
		int tag = 0;
		for( int j = 0; j < somaCentroids.size(); j++)
		{
			double distance = sqrt( pow( (double)somaCentroids[j][0] - x, 2) + pow( (double)somaCentroids[j][1] - y, 2) + pow( (double)somaCentroids[j][2] - z, 2));
			if( distance < 20)
			{
				tag = j;
				break;
			}
			if( distance < minDistance )
			{
				minDistance = distance;
				tag = j;
			}
		}

#pragma omp atomic
		debrisIntensity[tag] += (int)inputImage->GetPixel( debrisSeeds[i]);
	}

	vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();	
	vtkSmartPointer<vtkDoubleArray> column1 = vtkSmartPointer<vtkDoubleArray>::New();
	column1->SetName("centroid_x");
	table->AddColumn(column1);

	vtkSmartPointer<vtkDoubleArray> column2 = vtkSmartPointer<vtkDoubleArray>::New();
	column2->SetName("centroid_y");
	table->AddColumn(column2);

	vtkSmartPointer<vtkDoubleArray> column3 = vtkSmartPointer<vtkDoubleArray>::New();
	column3->SetName("centroid_z");
	table->AddColumn(column3);

	vtkSmartPointer<vtkDoubleArray> column4 = vtkSmartPointer<vtkDoubleArray>::New();
	column4->SetName("debris_intensity");
	table->AddColumn(column4);

	for( int i = 0; i < debrisIntensity.size(); i++)
	{
		vtkSmartPointer<vtkVariantArray> row = vtkSmartPointer<vtkVariantArray>::New();
		row->InsertNextValue( vtkVariant( somaCentroids[i][0]));
		row->InsertNextValue( vtkVariant( somaCentroids[i][1]));
		row->InsertNextValue( vtkVariant( somaCentroids[i][2]));
		row->InsertNextValue( vtkVariant( debrisIntensity[i]));
		table->InsertNextRow(row);
	}
	ftk::SaveTable( "DebrisTable.txt", table);
}

SomaExtractor::ProbImageType::Pointer SomaExtractor::GenerateSeedPoints(OutputImageType::Pointer inputImgPt, std::vector< itk::Index<3> > &somaCentroids)
{
	int size1 = inputImgPt->GetLargestPossibleRegion().GetSize()[0];
	int size2 = inputImgPt->GetLargestPossibleRegion().GetSize()[1];
	int size3 = inputImgPt->GetLargestPossibleRegion().GetSize()[2];

	unsigned char *in_Image;
	in_Image = (unsigned char *) malloc (size1*size2*size3);
	memset(in_Image, 0, size1*size2*size3*sizeof(unsigned char));

	typedef itk::ImageRegionConstIterator< OutputImageType > ConstIteratorType;
	ConstIteratorType pix_buf( inputImgPt, inputImgPt->GetRequestedRegion() );

	int ind=0;
	for ( pix_buf.GoToBegin(); !pix_buf.IsAtEnd(); ++pix_buf, ++ind )
	{
		in_Image[ind]=(pix_buf.Get());
	}

	yousef_nucleus_seg * NucleusSeg = new yousef_nucleus_seg();
	NucleusSeg->setDataImage(in_Image, size1, size2, size3, "");
	NucleusSeg->setParamsForSeedDetection(shift, scaleMin, scaleMax, regionXY, regionZ, useDistMap, sampling_ratio_XY_to_Z, minObjSize);
	std::cout<<std::endl << "Run Binarization"<<std::endl;
	NucleusSeg->runBinarization(num_bins);
	std::cout<<std::endl << "Run SeedDetection"<<std::endl;
	NucleusSeg->runSeedDetection();

	std::vector<Seed> seeds = NucleusSeg->getSeeds();
	somaCentroids.clear();
	somaCentroids.resize(seeds.size());
	for( int i = 0; i < seeds.size(); i++)
	{
		itk::Index<3> index;
		index[0] = seeds[i].x();
		index[1] = seeds[i].y();
		index[2] = seeds[i].z();
		somaCentroids[i] = index;
	}

	/// save binarized image as the speed image
	unsigned short *binImage = NucleusSeg->getBinImage();
	ProbImageType::Pointer binImagePtr = ProbImageType::New();
	ProbImageType::PointType origin;
	origin[0] = 0; 
	origin[1] = 0;    
	origin[2] = 0;    
	binImagePtr->SetOrigin( origin );

	ProbImageType::IndexType start;
	start[0] = 0;  
	start[1] = 0;     
	start[2] = 0;     
	ProbImageType::SizeType  size;
	size[0] = size1;  
	size[1] = size2; 
	size[2] = size3;  
	ProbImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );

	binImagePtr->SetRegions( region );
	binImagePtr->Allocate();
	binImagePtr->FillBuffer(0);
	binImagePtr->Update();

	//copy the bin image into the ITK image
	typedef itk::ImageRegionIteratorWithIndex< ProbImageType > IteratorType;
	IteratorType iterator1( binImagePtr, binImagePtr->GetRequestedRegion());		
	for( int i = 0; i < size1 * size2 * size3; i++)
	{				
		unsigned short val = (unsigned short)binImage[i];
		if( val > 0)
		{
			iterator1.Set(1);		
		}
		++iterator1;
	}

	//writeImage("binarize_image.nrrd",binImagePtr);
	delete NucleusSeg;
	return binImagePtr;
}

SomaExtractor::ProbImageType::Pointer SomaExtractor::GenerateSeedPoints( unsigned char* inputBuffer, int size1, int size2, int size3, std::vector< itk::Index<3> > &somaCentroids)
{
	unsigned char *in_Image = inputBuffer;

	/// save binarized image as the speed image
	yousef_nucleus_seg * NucleusSeg = new yousef_nucleus_seg();
	NucleusSeg->setDataImage(in_Image, size1, size2, size3, "");
	NucleusSeg->setParamsForSeedDetection(shift, scaleMin, scaleMax, regionXY, regionZ, useDistMap, sampling_ratio_XY_to_Z, 0);
	std::cout<<std::endl << "Run Binarization"<<std::endl;
	NucleusSeg->runBinarization(num_bins);
	std::cout<<std::endl << "Run SeedDetection"<<std::endl;
	NucleusSeg->runSeedDetection();

	std::vector<Seed> seeds = NucleusSeg->getSeeds();
	somaCentroids.clear();
	somaCentroids.resize(seeds.size());
	for( int i = 0; i < seeds.size(); i++)
	{
		itk::Index<3> index;
		index[0] = seeds[i].x();
		index[1] = seeds[i].y();
		index[2] = seeds[i].z();
		somaCentroids[i] = index;
	}

	/// save binarized image as the speed image
	unsigned short *binImage = NucleusSeg->getBinImage();
	ProbImageType::Pointer binImagePtr = ProbImageType::New();
	ProbImageType::PointType origin;
	origin[0] = 0; 
	origin[1] = 0;    
	origin[2] = 0;    
	binImagePtr->SetOrigin( origin );

	ProbImageType::IndexType start;
	start[0] = 0;  
	start[1] = 0;     
	start[2] = 0;     
	ProbImageType::SizeType  size;
	size[0] = size1;  
	size[1] = size2; 
	size[2] = size3;  
	ProbImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );

	binImagePtr->SetRegions( region );
	binImagePtr->Allocate();
	binImagePtr->FillBuffer(0);
	binImagePtr->Update();

	//copy the bin image into the ITK image
	typedef itk::ImageRegionIteratorWithIndex< ProbImageType > IteratorType;
	IteratorType iterator1( binImagePtr, binImagePtr->GetRequestedRegion());		
	for( int i = 0; i < size1 * size2 * size3; i++)
	{				
		unsigned short val = (unsigned short)binImage[i];
		if( val > 0)
		{
			iterator1.Set(1);		
		}
		++iterator1;				
	}

	//writeImage("binarize_image.nrrd",binImagePtr);
	delete NucleusSeg;
	return binImagePtr;
}

SomaExtractor::ProbImageType2D::Pointer SomaExtractor::SetInputImage2D(const char * fileName)
{
	ushortImage2DReader::Pointer reader = ushortImage2DReader::New();
	reader->SetFileName (fileName);	

	typedef itk::CastImageFilter<UShortImageType2D, ProbImageType2D> CasterType;
	CasterType::Pointer caster = CasterType::New();
	caster->SetInput(reader->GetOutput());
	try
	{
		caster->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cout << "Error in read image filter: " << err << std::endl;
	}
	ProbImageType2D::Pointer image = caster->GetOutput();
	return image;
}

SomaExtractor::ProbImageType2D::Pointer SomaExtractor::GetAverage(const char * channelName, int n, double sigma)
{
	typedef itk::AddImageFilter< ProbImageType2D, ProbImageType2D, ProbImageType2D > AddImageFilterType;
	AddImageFilterType::Pointer addFilter = AddImageFilterType::New();

	std::string str(channelName);
	str += std::string("1");
	str.append(".tif");
	ProbImageType2D::Pointer sumProbImage = SetInputImage2D(str.c_str());

	for(int i = 2; i <= n; i++)
	{	
		std::string str(channelName);
		std::stringstream ss;
		ss<< i;
		str += ss.str();
		str.append(".tif");
		std::cout<< str<<std::endl;
		ProbImageType2D::Pointer postprobImage = SetInputImage2D(str.c_str());
		addFilter->SetInput1(sumProbImage);
		addFilter->SetInput2(postprobImage);
		try
		{
			addFilter->Update();
		}
		catch( itk::ExceptionObject &err)
		{
			std::cout << "Error in SumFilter: " << err << std::endl; 
		}
		sumProbImage = addFilter->GetOutput();
	}

	typedef itk::DivideImageFilter <ProbImageType2D, ProbImageType2D, ProbImageType2D > DivideImageFilterType;
	DivideImageFilterType::Pointer divideImageFilter = DivideImageFilterType::New();
	divideImageFilter->SetInput1(sumProbImage);
	divideImageFilter->SetConstant2(n);
	try
	{
		divideImageFilter->Update();
	}
	catch( itk::ExceptionObject &err)
	{
		std::cout << "Error in SumFilter: " << err << std::endl; 
	}

	typedef itk::DiscreteGaussianImageFilter< ProbImageType2D, ProbImageType2D >  GaussinafilterType;
	GaussinafilterType::Pointer gaussianFilter = GaussinafilterType::New();
	gaussianFilter->SetInput( divideImageFilter->GetOutput());
	gaussianFilter->SetVariance(sigma);
	gaussianFilter->Update();
	ProbImageType2D::Pointer probImage = gaussianFilter->GetOutput();

	typedef itk::ImageFileWriter< ProbImageType2D > ufloatImage2DWriter;
	ufloatImage2DWriter::Pointer writer = ufloatImage2DWriter::New();
	writer->SetInput(probImage);
	writer->SetFileName("Average.nrrd");
	writer->Update();
	return probImage;
}

SomaExtractor::ProbImageType2D::Pointer SomaExtractor::GetBackgroundImage(ProbImageType::Pointer image, double sigma)
{
	std::cout<< "Get Background Image."<<std::endl;
	typedef itk::GetAverageSliceImageFilter< ProbImageType, ProbImageType2D> AverageSliceImageFilterType;
	AverageSliceImageFilterType::Pointer averageFilter = AverageSliceImageFilterType::New();
	averageFilter->SetInput(image);
	averageFilter->SetAccumulateDimension(2);
	averageFilter->SetAverage( true);

	std::cout<< "DiscreteGaussianImageFilter."<<std::endl;
	typedef itk::DiscreteGaussianImageFilter< ProbImageType2D, ProbImageType2D >  GaussinafilterType;
	GaussinafilterType::Pointer gaussianFilter = GaussinafilterType::New();
	gaussianFilter->SetInput( averageFilter->GetOutput());
	gaussianFilter->SetVariance(sigma);
	try
	{
		gaussianFilter->Update();
	}
	catch( itk::ExceptionObject &err)
	{
		std::cout << "Error in DiscreteGaussianImageFilter: " << err << std::endl; 
	}

	ProbImageType2D::Pointer averageImage = gaussianFilter->GetOutput();
	return averageImage;
}

SomaExtractor::ProbImageType2D::Pointer SomaExtractor::ExtractSlice(ProbImageType::Pointer image, int sliceId)
{
	ProbImageType::IndexType start;
	start[0] = 0;
	start[1] = 0;
	start[2] = sliceId;
	ProbImageType::SizeType size;
	size[0] = image->GetLargestPossibleRegion().GetSize()[0];
	size[1] = image->GetLargestPossibleRegion().GetSize()[1];
	size[2] = 0;
	ProbImageType::RegionType desiredRegion;
	desiredRegion.SetSize(size);
	desiredRegion.SetIndex(start);

	typedef itk::ExtractImageFilter< ProbImageType, ProbImageType2D> ExtractFilterType2D;
	ExtractFilterType2D::Pointer extractFilter = ExtractFilterType2D::New();
	extractFilter->SetInput(image);
	extractFilter->SetExtractionRegion(desiredRegion);
#if ITK_VERSION_MAJOR >= 4
	extractFilter->SetDirectionCollapseToIdentity(); // This is required.
#endif
	extractFilter->Update();
	ProbImageType2D::Pointer slice = extractFilter->GetOutput();
	return slice;
}

SomaExtractor::ProbImageType2D::Pointer SomaExtractor::GetBackgroundImageByFirstSlice(ProbImageType::Pointer image, double sigma)
{
	ProbImageType2D::Pointer probImage1 = ExtractSlice(image, 0);
	ProbImageType2D::Pointer probImage2 = ExtractSlice(image, 1);

	typedef itk::AddImageFilter< ProbImageType2D, ProbImageType2D, ProbImageType2D > AddImageFilterType;
	AddImageFilterType::Pointer addFilter = AddImageFilterType::New();
	addFilter->SetInput1(probImage1);
	addFilter->SetInput2(probImage2);

	typedef itk::DivideImageFilter <ProbImageType2D, ProbImageType2D, ProbImageType2D > DivideImageFilterType;
	DivideImageFilterType::Pointer divideImageFilter = DivideImageFilterType::New();
	divideImageFilter->SetInput1(addFilter->GetOutput());
	divideImageFilter->SetConstant2(2);

	typedef itk::DiscreteGaussianImageFilter< ProbImageType2D, ProbImageType2D >  GaussinafilterType;
	GaussinafilterType::Pointer gaussianFilter = GaussinafilterType::New();
	gaussianFilter->SetInput( divideImageFilter->GetOutput());
	gaussianFilter->SetVariance(sigma);
	try
	{
		gaussianFilter->Update();
	}
	catch( itk::ExceptionObject &err)
	{
		std::cout << "Error in DiscreteGaussianImageFilter: " << err << std::endl; 
	}

	ProbImageType2D::Pointer imagePt = gaussianFilter->GetOutput();

	//WriteFloat2DImage("BackgroundImage1.nrrd", imagePt);
	return imagePt;
}

void SomaExtractor::WriteFloat2DImage(const char* writeFileName, ProbImageType2D::Pointer image)
{
	probImage2DWriter::Pointer writer = probImage2DWriter::New();
	writer->SetFileName(writeFileName);
	writer->SetInput(image);
	writer->Update();
}

SomaExtractor::UShortImageType::Pointer SomaExtractor::DevideAndScaleToOriginalMean(ProbImageType::Pointer image, ProbImageType2D::Pointer backgroundImage, int border)
{
	//typedef itk::StatisticsImageFilter<ProbImageType> StatisticsImageFilterType;
	//StatisticsImageFilterType::Pointer statisticsImageFilter = StatisticsImageFilterType::New();
	//statisticsImageFilter->SetInput(image);
	//statisticsImageFilter->Update();
	//double mean1 = statisticsImageFilter->GetMean();
	////double std = statisticsImageFilter->GetSigma();
	//
	//typedef itk::StatisticsImageFilter<ProbImageType2D> StatisticsImageFilterType2;
	//StatisticsImageFilterType2::Pointer statisticsImageFilter2 = StatisticsImageFilterType2::New();
	//statisticsImageFilter2->SetInput(backgroundImage);
	//statisticsImageFilter2->Update();
	//double mean2 = statisticsImageFilter2->GetMean();
	//std::cout<< "Mean: "<<mean1<<"\t"<<mean2<<std::endl;

	////typedef itk::NormalizeImageFilter< ProbImageType, ProbImageType > NormalizeImageFilter;
	////NormalizeImageFilter::Pointer normalizeFilter = NormalizeImageFilter::New();
	////normalizeFilter->SetInput(image);
	////normalizeFilter->Update();
	////image = normalizeFilter->GetOutput();

	////typedef itk::CastImageFilter<ProbImageType2D, UShortImageType2D> CasterType;
	////CasterType::Pointer castFilter = CasterType::New();
	////castFilter->SetInput(backgroundImage);
	////ushortImage2DWriter::Pointer writer = ushortImage2DWriter::New();
	////writer->SetInput(castFilter->GetOutput());
	////writer->SetFileName("BackgroundImage.tif");
	////writer->Update();

	////typedef itk::NormalizeImageFilter< ProbImageType2D, ProbImageType2D > NormalizeImage2DFilter;
	////NormalizeImage2DFilter::Pointer normalize2DFilter = NormalizeImage2DFilter::New();
	////normalize2DFilter->SetInput(backgroundImage);
	////normalize2DFilter->Update();
	////backgroundImage = normalize2DFilter->GetOutput();
	//typedef itk::SubtractImageFilter <ProbImageType2D, ProbImageType2D > SubtractImageFilterType;
	//SubtractImageFilterType::Pointer subtractFilter = SubtractImageFilterType::New ();
	//subtractFilter->SetInput1(backgroundImage);
	//subtractFilter->SetConstant2(mean2-mean1);
	//subtractFilter->Update();
	//backgroundImage = subtractFilter->GetOutput();

	std::cout<< "Devide Image."<<std::endl;
	int width = image->GetLargestPossibleRegion().GetSize()[0];
	int height = image->GetLargestPossibleRegion().GetSize()[1];
	int depth = image->GetLargestPossibleRegion().GetSize()[2];

	int width2 = backgroundImage->GetLargestPossibleRegion().GetSize()[0];
	int height2 = backgroundImage->GetLargestPossibleRegion().GetSize()[1];
	if( width != width2 || height != height2)
	{
		std::cout<< "Input image and background image do not match!"<<std::endl;
		return NULL;
	}
	else
	{
		std::cout<< "Input image and background image match!"<<std::endl;
	}

	//typedef itk::DiscreteGaussianImageFilter< ProbImageType2D, ProbImageType2D >  GaussinafilterType;
	//GaussinafilterType::Pointer gaussianFilter = GaussinafilterType::New();
	//gaussianFilter->SetInput(backgroundImage);
	//gaussianFilter->SetVariance(5);
	//try
	//{
	//	gaussianFilter->Update();
	//}
	//catch( itk::ExceptionObject &err)
	//{
	//	std::cout << "Error in DiscreteGaussianImageFilter: " << err << std::endl; 
	//}

	//backgroundImage = gaussianFilter->GetOutput();
	//ShrinkPixel(backgroundImage, border);
	//WriteFloat2DImage("BackgroundImage.nrrd", backgroundImage);

	StatisticsImageFilterType::Pointer statisticsImageFilter1 = StatisticsImageFilterType::New();
	statisticsImageFilter1->SetInput(image);
	statisticsImageFilter1->Update();
	double mean1 = statisticsImageFilter1->GetMean();
	std::cout<< "Pre Mean: "<< mean1<<std::endl;

	//std::ofstream ofpre("PreScaleMax.txt",ios::app); 
	//ofpre<<maxPixelPreDevide<<"\t";
	//ofpre.close();

	typedef itk::ImageRegionConstIterator< ProbImageType2D > Prob2DConstIteratorType;

	for(int i = 0; i < depth; i++)
	{
		Prob2DConstIteratorType backgroundIt( backgroundImage, backgroundImage->GetLargestPossibleRegion());
		backgroundIt.GoToBegin();
		ProbImageType::IndexType desiredStart;
		desiredStart[0] = 0;
		desiredStart[1] = 0;
		desiredStart[2] = i;

		ProbImageType::SizeType desiredSize;
		desiredSize[0] = width;
		desiredSize[1] = height;
		desiredSize[2] = 1;

		ProbImageType::RegionType desiredRegion(desiredStart, desiredSize);
		ProbIteratorType imageIt( image, desiredRegion);
		imageIt.GoToBegin();

		while( !backgroundIt.IsAtEnd() && !imageIt.IsAtEnd())
		{
			double val = imageIt.Get();
			if( backgroundIt.Get() > 0)
			{
				val = (double)imageIt.Get() / (double)backgroundIt.Get();
			}
			else
			{
				val = 1;
			}
			imageIt.Set( val);
			++imageIt;
			++backgroundIt;
		}
	}

	StatisticsImageFilterType::Pointer statisticsImageFilter2 = StatisticsImageFilterType::New();
	statisticsImageFilter2->SetInput(image);
	statisticsImageFilter2->Update();
	double mean2 = statisticsImageFilter2->GetMean();
	double max2 = statisticsImageFilter2->GetMaximum();

	std::cout<< "Max: "<<max2<<"\t"<<"Mean: "<< mean2<<std::endl;
	std::cout<< "Multiply by "<< mean1 / mean2<<std::endl;

	typedef itk::MultiplyImageFilter<ProbImageType, ProbImageType, ProbImageType> MultiplyImageFilterType;
	MultiplyImageFilterType::Pointer multiplyImageFilter = MultiplyImageFilterType::New();
	multiplyImageFilter->SetInput(image);
	multiplyImageFilter->SetConstant(mean1 / mean2);

	typedef itk::CastImageFilter<ProbImageType, UShortImageType> CasterType2;
	CasterType2::Pointer castFilter2 = CasterType2::New();
	castFilter2->SetInput(multiplyImageFilter->GetOutput());
	castFilter2->Update();
	UShortImageType::Pointer rtnImage = castFilter2->GetOutput();
	return rtnImage;
}

SomaExtractor::UShortImageType::Pointer SomaExtractor::DevideAndScale(ProbImageType::Pointer oriImage, ProbImageType2D::Pointer backgroundImage, double median, double ratio_threshold)
{
	bool bAutoThreshold = false;
	if(ratio_threshold > 1e-6)
	{
		StatisticsImageFilterType::Pointer statisticsImageFilter = StatisticsImageFilterType::New();
		statisticsImageFilter->SetInput(oriImage);
		statisticsImageFilter->Update();
		double mean = statisticsImageFilter->GetMean();
		double std = statisticsImageFilter->GetSigma();
		double ratio = mean / std;
		if( ratio < ratio_threshold)
		{
			bAutoThreshold = true;
		}
	}

	typedef itk::MedianImageFilter< OutputImageType, OutputImageType >  MedianFilterType;
	MedianFilterType::Pointer medianFilter = MedianFilterType::New();
	OutputImageType::Pointer BinarizeImagePtr;
	if(bAutoThreshold)
	{
		HuangThresholdFilter::Pointer huangThresholdFilter = HuangThresholdFilter::New();
		huangThresholdFilter->SetInput(oriImage);
		huangThresholdFilter->SetNumberOfHistogramBins( 256);
		huangThresholdFilter->SetOutsideValue(1);
		huangThresholdFilter->SetInsideValue(0);
		medianFilter->SetInput(huangThresholdFilter->GetOutput());
		medianFilter->SetRadius(2);
		medianFilter->Update();
		BinarizeImagePtr = medianFilter->GetOutput();
	}

	ProbImageType::Pointer image = oriImage;

	int width = image->GetLargestPossibleRegion().GetSize()[0];
	int height = image->GetLargestPossibleRegion().GetSize()[1];
	int depth = image->GetLargestPossibleRegion().GetSize()[2];

	int width2 = backgroundImage->GetLargestPossibleRegion().GetSize()[0];
	int height2 = backgroundImage->GetLargestPossibleRegion().GetSize()[1];
	if( width != width2 || height != height2)
	{
		std::cout<< "Input image and background image do not match!"<<std::endl;
		return NULL;
	}
	else
	{
		//std::cout<< "Input image and background image match!"<<std::endl;
	}

	for(int i = 0; i < depth; i++)
	{
		Prob2DConstIteratorType backgroundIt( backgroundImage, backgroundImage->GetLargestPossibleRegion());
		backgroundIt.GoToBegin();
		ProbImageType::IndexType desiredStart;
		desiredStart[0] = 0;
		desiredStart[1] = 0;
		desiredStart[2] = i;

		ProbImageType::SizeType desiredSize;
		desiredSize[0] = width;
		desiredSize[1] = height;
		desiredSize[2] = 1;

		ProbImageType::RegionType desiredRegion(desiredStart, desiredSize);
		ProbIteratorType imageIt( image, desiredRegion);
		imageIt.GoToBegin();

		while( !backgroundIt.IsAtEnd() && !imageIt.IsAtEnd())
		{
			double val = imageIt.Get();
			if( backgroundIt.Get() > 0)
			{
				val = (double)imageIt.Get() / (double)backgroundIt.Get();
			}
			else
			{
				val = 1;
			}
			imageIt.Set( val);
			++imageIt;
			++backgroundIt;
		}
	}

	ProbImageType2D::Pointer imageSlice = ExtractSlice(image, 0);
	typedef itk::MedianImageFunction< ProbImageType2D > MedianImageFunctionType;
	MedianImageFunctionType::Pointer medianImageFunction = MedianImageFunctionType::New();
	medianImageFunction->SetInputImage( imageSlice);
	int r = width > height ? height : width;
	unsigned int radius = r / 2;
	itk::Index<2> index;
	index[0] = width / 2;
	index[1] = height / 2;
	//std::cout<< "Estimated Radius: "<<radius<<std::endl;
	medianImageFunction->SetNeighborhoodRadius(radius);
	double median2 = medianImageFunction->EvaluateAtIndex(index);
	//std::cout<< "Multiply by "<< median / median2<<std::endl;

	typedef itk::MultiplyImageFilter<ProbImageType, ProbImageType, ProbImageType> MultiplyImageFilterType;
	MultiplyImageFilterType::Pointer multiplyImageFilter = MultiplyImageFilterType::New();
	multiplyImageFilter->SetInput(image);
	multiplyImageFilter->SetConstant(median / median2);
	multiplyImageFilter->Update();
	image = multiplyImageFilter->GetOutput();

	if( bAutoThreshold)
	{
		//std::cout<< "AutoThreshold..."<<std::endl;
		ProbIteratorType oriImageIt( image, image->GetLargestPossibleRegion());
		UcharIteratorType imageIt( BinarizeImagePtr, BinarizeImagePtr->GetLargestPossibleRegion());
		while( !imageIt.IsAtEnd())
		{
			if( imageIt.Get() == 0)   // background
			{
				oriImageIt.Set(0);
			}
			++oriImageIt;
			++imageIt;
		}
	}

	/// remove border by 10 pixel
	ProbImageType::Pointer extractImage = RemoveImageBorderByPixel(image, 10);

	typedef itk::CastImageFilter<ProbImageType, UShortImageType> CasterType2;
	CasterType2::Pointer castFilter2 = CasterType2::New();
	castFilter2->SetInput(extractImage);
	castFilter2->Update();
	UShortImageType::Pointer rtnImage = castFilter2->GetOutput();
	return rtnImage;
}

SomaExtractor::ProbImageType::Pointer SomaExtractor::RemoveImageBorderByPixel(ProbImageType::Pointer image, int border)
{
	int width = image->GetLargestPossibleRegion().GetSize()[0];
	int height = image->GetLargestPossibleRegion().GetSize()[1];
	int depth = image->GetLargestPossibleRegion().GetSize()[2];
	ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
	ProbImageType::IndexType start;
	start[0] = border;
	start[1] = border;
	start[2] = 0;
	ProbImageType::SizeType size;
	size[0] = width - 2 * border;
	size[1] = height - 2 * border;
	size[2] = depth;
	ProbImageType::RegionType desiredRegion(start, size);
	extractFilter->SetExtractionRegion(desiredRegion);
	extractFilter->SetInput(image);
#if ITK_VERSION_MAJOR >= 4
	extractFilter->SetDirectionCollapseToIdentity(); 
#endif
	extractFilter->Update();

	ProbImageType::Pointer output = extractFilter->GetOutput();
	int awidth = output->GetLargestPossibleRegion().GetSize()[0];
	int aheight = output->GetLargestPossibleRegion().GetSize()[1];
	int adepth = output->GetLargestPossibleRegion().GetSize()[2];
	//std::cout<<width<<"\t"<<awidth<<std::endl<<height<<"\t"<<aheight<<std::endl<<depth<<"\t"<<adepth<<std::endl;
	return output;
}

SomaExtractor::ProbImageType2D::Pointer SomaExtractor::adjustMeanStd(ProbImageType2D::Pointer image, double globalMean, double globalStd)
{
	typedef itk::StatisticsImageFilter<ProbImageType2D> StatisticsImageFilter2DType;
	StatisticsImageFilter2DType::Pointer statisticsImageFilter = StatisticsImageFilter2DType::New();
	statisticsImageFilter->SetInput(image);
	statisticsImageFilter->Update();
	double mean = statisticsImageFilter->GetMean();
	double std = statisticsImageFilter->GetSigma();
	std::cout<< "mean: "<< mean<<" sigma "<< std<<std::endl;

	SubtractImageFilterType2D::Pointer subtractFilter = SubtractImageFilterType2D::New();
	subtractFilter->SetInput(image);
	subtractFilter->SetConstant2(mean);

	DivideImageFilterType2D::Pointer divideImageFilter = DivideImageFilterType2D::New();
	divideImageFilter->SetInput1(subtractFilter->GetOutput());
	divideImageFilter->SetConstant2(std / globalStd);

	AddImageFilterType2D::Pointer addImageFilter = AddImageFilterType2D::New();
	addImageFilter->SetInput(divideImageFilter->GetOutput());
	addImageFilter->SetConstant2(globalMean);
	addImageFilter->Update(); 

	ProbImageType2D::Pointer rtnImage = addImageFilter->GetOutput();
	return rtnImage;
}

void SomaExtractor::CaculateMeanStd(std::string fileName, ProbImageType::Pointer image)
{
	StatisticsImageFilterType::Pointer statisticsImageFilter2 = StatisticsImageFilterType::New();
	statisticsImageFilter2->SetInput(image);
	statisticsImageFilter2->Update();
	double mean = statisticsImageFilter2->GetMean();
	double std = statisticsImageFilter2->GetSigma();
	std::ofstream ofs(fileName.c_str(), fstream::app);
	ofs<< mean <<"\t"<<std<<"\t"<< mean/std<<std::endl;
	ofs.close();
}

void SomaExtractor::ShrinkPixel(ProbImageType2D::Pointer image, int border)
{
	int width = image->GetLargestPossibleRegion().GetSize()[0];
	int height = image->GetLargestPossibleRegion().GetSize()[1];
	for( int i = 0; i < width; i++)
	{
		for( int j = 0; j < border; j++)
		{
			ProbImageType2D::IndexType pixelIndex1;
			pixelIndex1[0] = i; // x position
			pixelIndex1[1] = j; // y position

			ProbImageType2D::IndexType pixelIndex2;
			pixelIndex2[0] = i; // x position
			pixelIndex2[1] = border; // y position

			image->SetPixel(pixelIndex1, image->GetPixel(pixelIndex2));
		}
	}

	for( int i = width - border; i < width; i++)
	{
		for( int j = 0; j < height; j++)
		{
			ProbImageType2D::IndexType pixelIndex1;
			pixelIndex1[0] = i; // x position
			pixelIndex1[1] = j; // y position

			ProbImageType2D::IndexType pixelIndex2;
			pixelIndex2[0] = width - border - 1; // x position
			pixelIndex2[1] = j; // y position

			image->SetPixel(pixelIndex1, image->GetPixel(pixelIndex2));
		}
	}

	for( int i = 0; i < width; i++)
	{
		for( int j = height - border; j < height; j++)
		{
			ProbImageType2D::IndexType pixelIndex1;
			pixelIndex1[0] = i; // x position
			pixelIndex1[1] = j; // y position

			ProbImageType2D::IndexType pixelIndex2;
			pixelIndex2[0] = i; // x position
			pixelIndex2[1] = height - border - 1; // y position
			image->SetPixel(pixelIndex1, image->GetPixel(pixelIndex2));
		}
	}

	for( int i = 0; i < border; i++)
	{
		for( int j = 0; j < height; j++)
		{
			ProbImageType2D::IndexType pixelIndex1;
			pixelIndex1[0] = i; // x position
			pixelIndex1[1] = j; // y position

			ProbImageType2D::IndexType pixelIndex2;
			pixelIndex2[0] = border ; // x position
			pixelIndex2[1] = j; // y position
			image->SetPixel(pixelIndex1, image->GetPixel(pixelIndex2));
		}
	}
}

void SomaExtractor::NormalizeUsingBackgroundImage(ProbImageType2D::Pointer image, ProbImageType2D::Pointer backgroundimage, double sigma)
{
	StatisticsImageFilterType2D::Pointer statisticsImageFilter = StatisticsImageFilterType2D::New();
	statisticsImageFilter->SetInput(image);
	statisticsImageFilter->Update();
	double meanImage = statisticsImageFilter->GetMean();
	double stdImage = statisticsImageFilter->GetSigma();
	std::cout<< "mean: "<< meanImage<<" sigma "<< stdImage<<std::endl;

	ShrinkPixel(backgroundimage, (int)sigma);

	statisticsImageFilter->SetInput(backgroundimage);
	statisticsImageFilter->Update();
	double meanBackground = statisticsImageFilter->GetMean();
	double stdBackground = statisticsImageFilter->GetSigma();
	std::cout<< "mean: "<< meanBackground<<" sigma "<< stdBackground<<std::endl;

	ProbImageType2D::Pointer newBackground = adjustMeanStd(backgroundimage, meanImage, stdImage);

	ThresholdImageFilterType::Pointer thresholdImageFilter = ThresholdImageFilterType::New();
	thresholdImageFilter->SetInput(newBackground);
	thresholdImageFilter->SetLower(0);
	thresholdImageFilter->SetOutsideValue(0);

	typedef itk::DiscreteGaussianImageFilter< ProbImageType2D, ProbImageType2D >  GaussinafilterType;
	GaussinafilterType::Pointer gaussianFilter = GaussinafilterType::New();
	gaussianFilter->SetInput( thresholdImageFilter->GetOutput());
	gaussianFilter->SetVariance(sigma);
	try
	{
		gaussianFilter->Update();
	}
	catch( itk::ExceptionObject &err)
	{
		std::cout << "Error in DiscreteGaussianImageFilter: " << err << std::endl; 
	}

	ProbImageType2D::Pointer imagePt = gaussianFilter->GetOutput();
	ShrinkPixel(imagePt, (int)sigma);

	WriteFloat2DImage("BackgroundImage.nrrd", imagePt);
}

void SomaExtractor::writeUnshort2D(const char *fileName, UShortImageType2D::Pointer image)
{
	ushortImage2DWriter::Pointer writeImage = ushortImage2DWriter::New();
	writeImage->SetFileName(fileName);
	writeImage->SetInput(image);
	writeImage->Update();
}


SomaExtractor::ProbImageType2D::Pointer SomaExtractor::SetInputImageFloat2D(const char *fileName)
{
	probImage2DReader::Pointer reader = probImage2DReader::New();
	reader->SetFileName (fileName);	

	try
	{
		reader->Update();
	}
	catch ( itk::ExceptionObject &err)
	{
		std::cout << "Error in prob_image_writer: " << err << std::endl; 
	}

	ProbImageType2D::Pointer rtnImage = reader->GetOutput();
	return rtnImage;
}

SomaExtractor::ProbImageType::Pointer SomaExtractor::SetInputImageFloat(const char *fileName)
{
	probImageReader::Pointer reader = probImageReader::New();
	reader->SetFileName (fileName);	

	try
	{
		reader->Update();
	}
	catch ( itk::ExceptionObject &err)
	{
		std::cout << "Error in prob_image_writer: " << err << std::endl; 
	}

	ProbImageType::Pointer rtnImage = reader->GetOutput();
	return rtnImage;
}

SomaExtractor::UShortImageType::Pointer SomaExtractor::RescaleImage(ProbImageType::Pointer image, double globalMax, double intensityMax)
{
	typedef itk::MinimumMaximumImageCalculator <ProbImageType> ImageCalculatorFilterType;

	ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New();
	imageCalculatorFilter->SetImage(image);
	imageCalculatorFilter->Compute();
	double maxPixel = imageCalculatorFilter->GetMaximum();
	std::cout<<"Max Pixel: "<<maxPixel<<std::endl;

	RescaleFloatFilterType::Pointer rescalePointer = RescaleFloatFilterType::New();
	rescalePointer->SetInput(image);
	rescalePointer->SetOutputMaximum(maxPixel / globalMax * intensityMax);
	rescalePointer->SetOutputMinimum(0);

	typedef itk::CastImageFilter<ProbImageType, UShortImageType> CasterType2;
	CasterType2::Pointer castFilter2 = CasterType2::New();
	castFilter2->SetInput(rescalePointer->GetOutput());
	castFilter2->Update();
	UShortImageType::Pointer rtnImage = castFilter2->GetOutput();
	return rtnImage;
}

void SomaExtractor::writeImage(const char* writeFileName, UShortImageType::Pointer image)
{
	ushortImageWriter::Pointer writer = ushortImageWriter::New();
	writer->SetFileName(writeFileName);
	writer->SetInput(image);
	try
	{
		writer->Update();
	}
	catch ( itk::ExceptionObject &err)
	{
		std::cout << "Error in bin_image_writer: " << err << std::endl; 
	}
}

SomaExtractor::ProbImageType::Pointer SomaExtractor::OtsuThresholdImage(OutputImageType::Pointer image)
{
	OtsuThresholdImageFilterType::Pointer otsuFilter = OtsuThresholdImageFilterType::New();
	otsuFilter->SetInput(image);
	otsuFilter->Update();
	ProbImageType::Pointer probImage = otsuFilter->GetOutput();
	return probImage;
}