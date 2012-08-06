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
#include "itkInvertIntensityImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include <itkShapeLabelObject.h>
#include "itkCurvatureAnisotropicDiffusionImageFilter.h"

#include <itkSigmoidImageFilter.h>
#include <math.h>
#include "itkImageSliceIteratorWithIndex.h"
#ifdef _OPENMP
#include "omp.h"
#endif

#define N 11
std::string SomaInfo[N]={"ID", "centroid_x", "centroid_y", "centroid_z", "volume", "eccentricity", "elongation", "orientation", 
						"majorAxisLength", "minorAxisLength", "surface_area_volume_ratio"};


SomaExtractor::SomaExtractor()
{
}

SomaExtractor::~SomaExtractor()
{
}

SomaExtractor::ProbImageType::Pointer SomaExtractor::SetInputImage(const char * fileName)
{
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName (fileName);	

	typedef itk::CastImageFilter<OutputImageType, ProbImageType> CasterType;
    CasterType::Pointer caster = CasterType::New();
    caster->SetInput(reader->GetOutput());
	caster->Update();
	inputImage = caster->GetOutput();
	width = inputImage->GetLargestPossibleRegion().GetSize()[0];
	height = inputImage->GetLargestPossibleRegion().GetSize()[1];
	depth = inputImage->GetLargestPossibleRegion().GetSize()[2];
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

SomaExtractor::SegmentedImageType::Pointer SomaExtractor::SetInitalContourImage(const char * fileName)
{
	//somaImageReader::Pointer reader = somaImageReader::New();
	//reader->SetFileName(fileName);	
	//reader->Update();
	//SegmentedImageType::Pointer imagePtr = reader->GetOutput();
	OutputImageType::Pointer image = Read8BitImage(fileName);
	typedef itk::CastImageFilter<OutputImageType, SegmentedImageType> CasterType;
    CasterType::Pointer caster = CasterType::New();
	caster->SetInput(image);
	caster->Update();
	SegmentedImageType::Pointer imagePtr = caster->GetOutput();
	return imagePtr;
}

void SomaExtractor::SetInputImage( ProbImageType::Pointer probImage)
{
	inputImage = probImage;
}

void SomaExtractor::ReadSeedpoints(const char * fileName, std::vector< itk::Index<3> > &seedVec, bool bNucleusTable)
{
	std::cout << "ReadSeedpoints" << std::endl;
	seedVec.clear();
	int endX = startX + width;
	int endY = startY + height;
	int endZ = startZ + depth;
	std::cout<< "Bounding Area:"<<startX<<"\t"<<startY<<"\t"<<endX<<"\t"<<endY<<std::endl;
	if( false == bNucleusTable)
	{
		vtkSmartPointer<vtkTable> table = ftk::LoadXYZTable( std::string(fileName));
		for( vtkIdType i = 0; i < table->GetNumberOfRows(); i++)
		{
			itk::Index<3> index;
			index[0] = table->GetValue(i,0).ToUnsignedInt();
			index[1] = table->GetValue(i,1).ToUnsignedInt();
			index[2] = table->GetValue(i,2).ToUnsignedInt();
			if( index[0] >= startX && index[0] <= endX && index[1] >= startY && index[1] <= endY && index[2] >= startZ && index[2] <= endZ)
			{
				itk::Index<3> newIndex;
				newIndex[0] = index[0] - startX;
				newIndex[1] = index[1] - startY;
				newIndex[2] = index[2] - startZ;
				seedVec.push_back(newIndex);
			}
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
				if( index[0] >= startX && index[0] <= endX && index[1] >= startY && index[1] <= endY && index[2] >= startZ && index[2] <= endZ)
				{
					itk::Index<3> newIndex;
					newIndex[0] = index[0] - startX;
					newIndex[1] = index[1] - startY;
					newIndex[2] = index[2] - startZ;
					seedVec.push_back(newIndex);
				}
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

	SetParamValue<double>(opts, "-sigmoid_alfa", alfa, 30);
	SetParamValue<double>(opts, "-sigmoid_beta", beta, 3);
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

	SetParamValue<int>(opts, "-start_x", startX, 0);
	SetParamValue<int>(opts, "-start_y", startY, 0);
	SetParamValue<int>(opts, "-start_z", startZ, 0);
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
	typedef itk::CastImageFilter<ProbImageType, OutputImageType> CasterType;
    CasterType::Pointer caster = CasterType::New();
	WriterType::Pointer soma_image_writer = WriterType::New();

	if(bscale)
	{
		RescaleFloatFilterType::Pointer rescaleFilter = RescaleFloatFilterType::New();
		rescaleFilter->SetInput(image);
		rescaleFilter->SetOutputMinimum(0);
		rescaleFilter->SetOutputMaximum(255);
		rescaleFilter->Update();
		caster->SetInput(rescaleFilter->GetOutput());
	}
	else
	{
	caster->SetInput(image);
	}

	soma_image_writer->SetFileName(writeFileName);
	soma_image_writer->SetInput(caster->GetOutput());

	try
	{
		soma_image_writer->Update();
	}
	catch ( itk::ExceptionObject &err)
	{
		std::cout << "Error in bin_image_writer: " << err << std::endl; 
	}
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
	ProbImageSliceType::Pointer slice = extractFilter->GetOutput();
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

// Shape Detection Active Contour Without GVF:
//	double alfa;
//	double beta;
//	int timethreshold; 
//	double curvatureScaling; 
//	double rmsThres;
//	int holeSize;
//	int minObjSize;
SomaExtractor::SegmentedImageType::Pointer SomaExtractor::SegmentSoma( ProbImageType::Pointer input, std::vector< itk::Index<3> > &somaCentroids)
{
	//typedef itk::NearestNeighborInterpolateImageFunction< ProbImageType, float>  InterpolatorType;
	//InterpolatorType::Pointer I_Interpolator = InterpolatorType::New();
	//I_Interpolator->SetInputImage(input);

	int SM = input->GetLargestPossibleRegion().GetSize()[0];
    int SN = input->GetLargestPossibleRegion().GetSize()[1];
    int SZ = input->GetLargestPossibleRegion().GetSize()[2];
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

	std::cout << "RescaleIntensity "<<alfa<<"\t"<<beta<<std::endl;
    SigmoidImageFilterType::Pointer sigmoidFilter = SigmoidImageFilterType::New();
	sigmoidFilter->SetInput(input);
	sigmoidFilter->SetOutputMinimum(0);
	sigmoidFilter->SetOutputMaximum(1);
	sigmoidFilter->SetAlpha(alfa);
	sigmoidFilter->SetBeta(beta);
	sigmoidFilter->Update();
	writeImage("Sigmoid.tif",sigmoidFilter->GetOutput());

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
    fastMarching->SetOutputSize( input->GetBufferedRegion().GetSize() );
    fastMarching->SetStoppingValue(  timethreshold);
	fastMarching->SetSpeedConstant( 1.0 );
	fastMarching->Update();

	std::cout<< "Shape Detection " <<curvatureScaling<<"\t"<<rmsThres<<std::endl;
	
	ShapeDetectionFilterType::Pointer shapeDetection = ShapeDetectionFilterType::New();
	
	shapeDetection->SetPropagationScaling(1.0);
	shapeDetection->SetCurvatureScaling(curvatureScaling);
	shapeDetection->SetMaximumRMSError( rmsThres);
	shapeDetection->SetNumberOfIterations( maxIterations);

	shapeDetection->SetInput( fastMarching->GetOutput());
	shapeDetection->SetFeatureImage( sigmoidFilter->GetOutput());
	shapeDetection->Update();
	std::cout << "No. elpased iterations: " << shapeDetection->GetElapsedIterations() << std::endl;

	std::cout<< "Thresholding..."<<endl;
	
    BinaryThresholdingFilterType::Pointer thresholder = BinaryThresholdingFilterType::New();
	thresholder->SetLowerThreshold( -10000);
	thresholder->SetUpperThreshold(0.0);                                                           
	thresholder->SetOutsideValue( 0);
	thresholder->SetInsideValue(255);
	thresholder->SetInput( shapeDetection->GetOutput());

	/// morphological closing
	KernelType ball;
	KernelType::SizeType ballSize;
	ballSize[0] = holeSize;
	ballSize[1] = holeSize;
	ballSize[2] = 1;
	ball.SetRadius(ballSize);
	ball.CreateStructuringElement();
	CloseFilterType::Pointer closeFilter = CloseFilterType::New();
	closeFilter->SetInput( thresholder->GetOutput());
	closeFilter->SetKernel( ball );
	closeFilter->SetForegroundValue( 255);

	/// Label image
	LabelFilterType::Pointer label = LabelFilterType::New();
	label->SetInput(closeFilter->GetOutput());

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

/// Generate Inital Contours by Daniel Distance Map
SomaExtractor::ProbImageType::Pointer SomaExtractor::GetInitalContourByDanielssonDistanceMap(SegmentedImageType::Pointer labelImage, double outlierExpand)
{
	SegThresholdingFilterType::Pointer thresholder = SegThresholdingFilterType::New();
	thresholder->SetInput(labelImage);
	thresholder->SetLowerThreshold(1); 
	thresholder->SetUpperThreshold(100000); 
	thresholder->SetOutsideValue( 0);
	thresholder->SetInsideValue(255);

	///generate bounding box
	LabelFilterType::Pointer label = LabelFilterType::New();
	label->SetInput(thresholder->GetOutput());

	LabelGeometryImageFilterType::Pointer labelGeometryImageFilter = LabelGeometryImageFilterType::New();
	labelGeometryImageFilter->SetInput( label->GetOutput());
	labelGeometryImageFilter->Update();
	std::cout<<"Number of objects: "<< labelGeometryImageFilter->GetNumberOfObjects()<<std::endl;

	LabelGeometryImageFilterType::LabelsType allLabels = labelGeometryImageFilter->GetLabels();
	LabelGeometryImageFilterType::LabelsType::iterator allLabelsIt =  allLabels.begin();
	
	for( allLabelsIt++; allLabelsIt != allLabels.end(); allLabelsIt++ )
    {
		LabelGeometryImageFilterType::LabelPixelType labelValue = *allLabelsIt;
		LabelGeometryImageFilterType::BoundingBoxType boundingBox = labelGeometryImageFilter->GetBoundingBox(labelValue);

	}
	/// Generate Distance Map within the bounding box(not yet finished)
	DanielssonDistanceMapFilterType::Pointer distanceMapFilter = DanielssonDistanceMapFilterType::New();
	typedef itk::CastImageFilter<SegmentedImageType, ProbImageType> CasterType;
    CasterType::Pointer caster = CasterType::New();
	caster->SetInput(thresholder->GetOutput());
	distanceMapFilter->SetInput(caster->GetOutput());

	SubtractImageFilterType::Pointer substractImageFilter = SubtractImageFilterType::New();
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

	ProbImageType::Pointer initalContourImage = substractImageFilter->GetOutput();
	writeImage("InitialContour.tif",initalContourImage, true);
	return initalContourImage;
}

// GVF Active Contour with GVF and curvature constraint:
// edge image:
// double sigma
// GVF:
// unsigned int numberOfIterations; 
// double noiseLevel;
// double alfa;
// double beta;
// double outlierExpandValue to expand the initial contour
// activecontour:
// double curvatureScaling; 
// double advectScaling
// double rmsThres;
// int minObjSize;
SomaExtractor::SegmentedImageType::Pointer SomaExtractor::SegmentSoma2( ProbImageType::Pointer input, SegmentedImageType::Pointer initialContour, std::vector< itk::Index<3> > &somaCentroids) 																	  
{
	int SM = input->GetLargestPossibleRegion().GetSize()[0];
    int SN = input->GetLargestPossibleRegion().GetSize()[1];
    int SZ = input->GetLargestPossibleRegion().GetSize()[2];
	std::cout<<SM<<"\t"<<SN<<"\t"<<SZ<<std::endl;

	ProbImageType::Pointer edgeImage = GetEdgePotentialMap(input, sigma);

	std::cout<<"Speed Image throught sigmoidFilter: "<<alfa<<"\t"<<beta<<std::endl;
	SigmoidImageFilterType::Pointer sigmoidFilter = SigmoidImageFilterType::New();
	sigmoidFilter->SetInput(edgeImage);
	sigmoidFilter->SetOutputMinimum(0);
	sigmoidFilter->SetOutputMaximum(1);
	sigmoidFilter->SetAlpha(alfa);
	sigmoidFilter->SetBeta(beta);
	sigmoidFilter->Update();
	
	std::cout<< "Initial Contour..."<<std::endl;
	ProbImageType::Pointer initialContourByDistanceMap = GetInitalContourByDanielssonDistanceMap(initialContour, outlierExpandValue);

	std::cout<<"GVF: "<<noiseLevel<<"\t"<<numberOfIterations<<std::endl;
	GradientIFilterType::Pointer gradientFilter = GradientIFilterType::New();
	GVFFilterType::Pointer gvfFilter = GVFFilterType::New();
	gradientFilter->SetInput(edgeImage);
	gvfFilter->SetInput(gradientFilter->GetOutput());
	gvfFilter->SetNoiseLevel(noiseLevel);
	gvfFilter->SetIterationNum(numberOfIterations);
	gvfFilter->Update();

	CastFlowFilterType::Pointer castFlowFilter = CastFlowFilterType::New();
	castFlowFilter->SetInput(gvfFilter->GetOutput());
	castFlowFilter->Update();

	std::cout<<"GVF Snake: "<<curvatureScaling<<"\t"<<advectScaling<<"\t"<<rmsThres<<std::endl;
	GeodesicActiveContourFilterType::Pointer GVF_snake = GeodesicActiveContourFilterType::New();
	GVF_snake->SetAutoGenerateSpeedAdvection(false);
	GVF_snake->SetInput(initialContourByDistanceMap);
	GVF_snake->SetFeatureImage(sigmoidFilter->GetOutput());
	GVF_snake->GenerateSpeedImage();
	GVF_snake->SetAdvectionImage(castFlowFilter->GetOutput());
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

		double ratio = (double)boundaryPixSize[LtoIMap[labelValue] ] / volume;
		row->InsertNextValue(vtkVariant(ratio));
		table->InsertNextRow(row);
    }

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