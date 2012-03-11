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
#include "itkBinaryThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include <vtkTable.h>
#include <vtkSmartPointer.h>

SomaExtractor::SomaExtractor()
{
}

void SomaExtractor::SetInputImage(char * fileName)
{
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName (fileName);	

	typedef itk::CastImageFilter<OutputImageType, ProbImageType> CasterType;
    CasterType::Pointer caster = CasterType::New();
    caster->SetInput(reader->GetOutput());
	caster->Update();
	binImage = caster->GetOutput();

	SM = binImage->GetLargestPossibleRegion().GetSize()[0];
    SN = binImage->GetLargestPossibleRegion().GetSize()[1];
    SZ = binImage->GetLargestPossibleRegion().GetSize()[2];
}

SomaExtractor::ProbImageType::Pointer SomaExtractor::GetFloatInputImage()
{
	return binImage;
}

void SomaExtractor::ReadSeedpoints(char * fileName, std::vector< itk::Index<3> > &seedVec)
{
	seedVec.clear();
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

void SomaExtractor::writeSomaImage(char* writeFileName)
{
	typedef itk::CastImageFilter<SegmentedImageType, OutputImageType> CasterType;
    CasterType::Pointer caster = CasterType::New();
    caster->SetInput(somaImage);

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

SomaExtractor::SegmentedImageType::Pointer SomaExtractor::SegmentSoma( ProbImageType::Pointer inputImage, std::vector< itk::Index<3> > &somaCentroids, 
													    int timethreshold, double curvatureScaling, double rmsThres)
{
	typedef itk::NearestNeighborInterpolateImageFunction< ProbImageType, float>  InterpolatorType;
	InterpolatorType::Pointer I_Interpolator = InterpolatorType::New();
	I_Interpolator->SetInputImage(inputImage);

	SM = inputImage->GetLargestPossibleRegion().GetSize()[0];
    SN = inputImage->GetLargestPossibleRegion().GetSize()[1];
    SZ = inputImage->GetLargestPossibleRegion().GetSize()[2];

	//move the seed points along z axis
    for( int i = 0; i < somaCentroids.size(); i++ )
	{
		SegmentedImageType::IndexType index;
		SegmentedImageType::IndexType index1;
		index[0] = index1[0] = somaCentroids[i][0];
		index[1] = index1[1] = somaCentroids[i][1];
		index[2] = somaCentroids[i][2];
		for( int j = 0; j < SZ; j++ )
		{
			index1[2] = j;
			if( I_Interpolator->EvaluateAtIndex(index1) > I_Interpolator->EvaluateAtIndex(index) )
			{
				somaCentroids[i][2] = j;
				index[2] = j;
			}
		}
	}

	/// rescaled image as speed image
	std::cout << "RescaleIntensity"<<endl;
	typedef itk::RescaleIntensityImageFilter< ProbImageType, ProbImageType> RescaleFilterType;
    RescaleFilterType::Pointer rescale = RescaleFilterType::New();
	rescale->SetInput( inputImage); 
    rescale->SetOutputMinimum( 0 );
    rescale->SetOutputMaximum( 1 );
    rescale->Update();

	/// generate distance map of the seeds using Fastmarching method
	std::cout<< "Generating Distance Map..." <<endl;

	clock_t SomaExtraction_start_time = clock();
	typedef  itk::FastMarchingImageFilter< ProbImageType, ProbImageType >    FastMarchingFilterType;
	FastMarchingFilterType::Pointer  fastMarching = FastMarchingFilterType::New();
	typedef FastMarchingFilterType::NodeContainer           NodeContainer;
    typedef FastMarchingFilterType::NodeType                NodeType;
    NodeContainer::Pointer seeds = NodeContainer::New();
    seeds->Initialize();

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
    fastMarching->SetOutputSize( inputImage->GetBufferedRegion().GetSize() );
    fastMarching->SetStoppingValue(  timethreshold);
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
	shapeDetection->SetMaximumRMSError( rmsThres);
	shapeDetection->SetNumberOfIterations( 800);

	shapeDetection->SetInput( fastMarching->GetOutput());
	shapeDetection->SetFeatureImage( rescale->GetOutput());
	shapeDetection->Update();
	std::cout << "No. elpased iterations: " << shapeDetection->GetElapsedIterations() << std::endl;
	std::cout << "RMS change: " << shapeDetection->GetRMSChange() << std::endl;
	std::cout<< "Total time for Shape Detection is: " << (clock() - SomaExtraction_start_time) / (float) CLOCKS_PER_SEC << std::endl;

	/// Get binarized image by thresholding
	std::cout<< "Thresholding..."<<endl;
	typedef itk::BinaryThresholdImageFilter< ProbImageType, SegmentedImageType> ShapeThresholdingFilterType;
    ShapeThresholdingFilterType::Pointer thresholder = ShapeThresholdingFilterType::New();
	thresholder->SetLowerThreshold( -10000);
	thresholder->SetUpperThreshold(0.0);
	thresholder->SetOutsideValue( 0);
	thresholder->SetInsideValue(255);
	thresholder->SetInput( shapeDetection->GetOutput());

	/// Label image
	typedef itk::ConnectedComponentImageFilter< SegmentedImageType, SegmentedImageType, SegmentedImageType> LabelFilterType;
	LabelFilterType::Pointer label = LabelFilterType::New();
	label->SetInput(thresholder->GetOutput());
	label->Update();
	somaImage = label->GetOutput();
	return somaImage;
}