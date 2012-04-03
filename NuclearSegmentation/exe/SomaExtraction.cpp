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
#include "itkLaplacianSharpeningImageFilter.h"
#include <itkGradientMagnitudeRecursiveGaussianImageFilter.h>
#include <itkSigmoidImageFilter.h>
#include <math.h>
#include "itkImageSliceIteratorWithIndex.h"

#define N 11
std::string SomaInfo[N]={"ID", "centroid_x", "centroid_y", "centroid_z", "volume", "eccentricity", "elongation", "orientation", 
						"majorAxisLength", "minorAxisLength", "surface_area_volume_ratio"};


SomaExtractor::SomaExtractor()
{
}

SomaExtractor::~SomaExtractor()
{
}

void SomaExtractor::SetInputImage(const char * fileName)
{
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName (fileName);	

	typedef itk::CastImageFilter<OutputImageType, ProbImageType> CasterType;
    CasterType::Pointer caster = CasterType::New();
    caster->SetInput(reader->GetOutput());
	caster->Update();
	inputImage = caster->GetOutput();
}

void SomaExtractor::SetInputImage( ProbImageType::Pointer probImage)
{
	inputImage = probImage;
}

SomaExtractor::ProbImageType::Pointer SomaExtractor::GetFloatInputImage()
{
	return inputImage;
}

void SomaExtractor::ReadSeedpoints(const char * fileName, std::vector< itk::Index<3> > &seedVec, bool bNucleusTable)
{
	std::cout << "ReadSeedpoints" << std::endl;
	seedVec.clear();

	if( false == bNucleusTable)
	{
		vtkSmartPointer<vtkTable> table = ftk::LoadXYZTable( std::string(fileName));
		for( vtkIdType i = 0; i < table->GetNumberOfRows(); i++)
		{
			itk::Index<3> index;
			index[0] = table->GetValue(i,0).ToUnsignedInt();
			index[1] = table->GetValue(i,1).ToUnsignedInt();
			index[2] = table->GetValue(i,2).ToUnsignedInt();
			seedVec.push_back(index);
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
				seedVec.push_back(index);
			}
		}
	}
}

void SomaExtractor::writeImage(const char* writeFileName, SegmentedImageType::Pointer image)
{
	typedef itk::CastImageFilter<SegmentedImageType, OutputImageType> CasterType;
    CasterType::Pointer caster = CasterType::New();
    caster->SetInput(image);

	WriterType::Pointer soma_image_writer = WriterType::New();
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

void SomaExtractor::writeImage(const char* writeFileName, ProbImageType::Pointer image)
{
	typedef itk::CastImageFilter<ProbImageType, OutputImageType> CasterType;
    CasterType::Pointer caster = CasterType::New();
    caster->SetInput(image);

	WriterType::Pointer soma_image_writer = WriterType::New();
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

void SomaExtractor::writeCentroids(const char* writeFileName, std::vector< itk::Index<3> > &seedVec)
{
	std::ofstream ofs(writeFileName);
	for( int i = 0; i< seedVec.size(); i++ )
	{
		ofs<< seedVec[i][0]<< "\t"<<seedVec[i][1]<<"\t"<<seedVec[i][2]<<std::endl;
	}
	ofs.close();
}

SomaExtractor::ProbImageType::Pointer SomaExtractor::EnhanceContrast( ProbImageType::Pointer inputImage, double alfa, double beta)
{
	//ProbImageType::IndexType start;
	//start[0] = 0;
	//start[1] = 0;
	//start[2] = centroid[2];
	//std::cout<< "Slide: "<<start[2]<<std::endl;
	//ProbImageType::SizeType size;
	//size[0] = inputImage->GetLargestPossibleRegion().GetSize()[0];
	//size[1] = inputImage->GetLargestPossibleRegion().GetSize()[1];
	//size[2] = 1;
	//ProbImageType::RegionType desiredRegion;
	//desiredRegion.SetSize(size);
	//desiredRegion.SetIndex(start);

	//RegionOfInterestFilter::Pointer regionFilter = RegionOfInterestFilter::New();
	//regionFilter->SetInput(inputImage);
	//regionFilter->SetRegionOfInterest(desiredRegion);

	//HuangThresholdFilter::Pointer huangThresholdFilter = HuangThresholdFilter::New();
	//huangThresholdFilter->SetInput(regionFilter->GetOutput());
	//huangThresholdFilter->SetNumberOfHistogramBins( 256);
	//huangThresholdFilter->Update();
	//double threshold = huangThresholdFilter->GetThreshold();

	//std::cout<< "HuangThreshold: "<< threshold<<std::endl;
	//BinaryProbThresholdingFilterType::Pointer thresholder = BinaryProbThresholdingFilterType::New();
	//thresholder->SetInput(inputImage);
	//thresholder->SetUpperThreshold(threshold);                                                           
	//thresholder->SetOutsideValue( 255);
	//thresholder->SetInsideValue( 0);
	//thresholder->Update();
	//ProbImageType::Pointer probImage = thresholder->GetOutput();
	
	SigmoidImageFilterType::Pointer sigmoidFilter = SigmoidImageFilterType::New();
	sigmoidFilter->SetInput(inputImage);
	sigmoidFilter->SetOutputMinimum(0);
	sigmoidFilter->SetOutputMaximum(255);
	sigmoidFilter->SetAlpha(alfa);
	sigmoidFilter->SetBeta(beta);
	sigmoidFilter->Update();
	ProbImageType::Pointer floatImage = sigmoidFilter->GetOutput();
	
	return floatImage;
}

SomaExtractor::SegmentedImageType::Pointer SomaExtractor::SegmentSoma( ProbImageType::Pointer input, std::vector< itk::Index<3> > &somaCentroids, 
													   double alfa, double beta, int timethreshold, double curvatureScaling, double rmsThres, int minObjSize)
{
	/// rescaled image as speed image
	//std::cout << "RescaleIntensity"<<endl;
    SigmoidImageFilterType::Pointer sigmoidFilter = SigmoidImageFilterType::New();
	sigmoidFilter->SetInput(input);
	sigmoidFilter->SetOutputMinimum(0);
	sigmoidFilter->SetOutputMaximum(1);
	sigmoidFilter->SetAlpha(alfa);
	sigmoidFilter->SetBeta(beta);
	sigmoidFilter->Update();
	ProbImageType::Pointer floatImage = sigmoidFilter->GetOutput();

	/// generate distance map of the seeds using Fastmarching method
	//std::cout<< "Generating Distance Map..." <<endl;

	//clock_t SomaExtraction_start_time = clock();
	
	FastMarchingFilterType::Pointer fastMarching = FastMarchingFilterType::New();
    NodeContainer::Pointer seeds = NodeContainer::New();
    seeds->Initialize();

	std::cout<< somaCentroids.size()<<std::endl;
	for( int i = 0; i< somaCentroids.size(); i++ )
	{
		ProbImageType::IndexType  seedPosition;
		seedPosition[0] = somaCentroids[i][0];
		seedPosition[1] = somaCentroids[i][1];
		seedPosition[2] = somaCentroids[i][2];
		NodeType node;
		const double seedValue = -3;
		node.SetValue( seedValue );
		node.SetIndex( seedPosition );
		seeds->InsertElement( i, node );
	}

	fastMarching->SetTrialPoints(  seeds);
    fastMarching->SetOutputSize( input->GetBufferedRegion().GetSize() );
    fastMarching->SetStoppingValue(  timethreshold);
	fastMarching->SetSpeedConstant( 1.0 );
	fastMarching->Update();
	//std::cout<< "Total time for Distance Map is: " << (clock() - SomaExtraction_start_time) / (float) CLOCKS_PER_SEC << std::endl;

	/// Shape Detection
	//SomaExtraction_start_time = clock();
	//std::cout<< "Shape Detection..." <<std::endl;
	
	ShapeDetectionFilterType::Pointer shapeDetection = ShapeDetectionFilterType::New();
	
	shapeDetection->SetPropagationScaling(1.0);
	shapeDetection->SetCurvatureScaling(curvatureScaling);
	shapeDetection->SetMaximumRMSError( rmsThres);
	shapeDetection->SetNumberOfIterations( 500);

	shapeDetection->SetInput( fastMarching->GetOutput());
	shapeDetection->SetFeatureImage( sigmoidFilter->GetOutput());
	shapeDetection->Update();
	//std::cout << "No. elpased iterations: " << shapeDetection->GetElapsedIterations() << std::endl;
	//std::cout << "RMS change: " << shapeDetection->GetRMSChange() << std::endl;
	//std::cout<< "Total time for Shape Detection is: " << (clock() - SomaExtraction_start_time) / (float) CLOCKS_PER_SEC << std::endl;

	/// Get binarized image by thresholding
	//std::cout<< "Thresholding..."<<endl;
	
    BinaryThresholdingFilterType::Pointer thresholder = BinaryThresholdingFilterType::New();
	thresholder->SetLowerThreshold( -10000);
	thresholder->SetUpperThreshold(0.0);                                                           
	thresholder->SetOutsideValue( 0);
	thresholder->SetInsideValue(255);
	thresholder->SetInput( shapeDetection->GetOutput());

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
