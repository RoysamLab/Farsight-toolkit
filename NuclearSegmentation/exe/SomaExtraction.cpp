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
#include "itkLabelGeometryImageFilter.h"
#include <itkShapeLabelObject.h>
#include "itkCurvatureAnisotropicDiffusionImageFilter.h"

#define N 11
std::string SomaInfo[N]={"ID", "centroid_x", "centroid_y", "centroid_z", "volume", "eccentricity", "elongation", "orientation", 
						"majorAxisLength", "minorAxisLength", "surface_area_volume_ratio"};


SomaExtractor::SomaExtractor()
{
	somaSeedQueue = NULL;
}

void SomaExtractor::SetInputImage(char * fileName)
{
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName (fileName);	

	typedef itk::CastImageFilter<OutputImageType, ProbImageType> CasterType;
    CasterType::Pointer caster = CasterType::New();
    caster->SetInput(reader->GetOutput());
	caster->Update();
	inputImage = caster->GetOutput();
}

void SomaExtractor::SmoothByDiffusionImageFilter()
{
	std::cout<< "Smoothing the image with edges preserved."<<std::endl;
	typedef itk::CurvatureAnisotropicDiffusionImageFilter<ProbImageType, ProbImageType > SmoothFilter;
	SmoothFilter::Pointer filter = SmoothFilter::New();
	filter->SetInput(inputImage);
	filter->SetNumberOfIterations(4);
	filter->SetTimeStep(0.0625);
	filter->SetConductanceParameter(3);
	inputImage = filter->GetOutput();
	writeImage("SmoothImage.tif",inputImage);
}

SomaExtractor::ProbImageType::Pointer SomaExtractor::GetFloatInputImage()
{
	return inputImage;
}

void SomaExtractor::ReadSeedpoints(char * fileName, std::vector< itk::Index<3> > &seedVec, bool bNucleusTable)
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

void SomaExtractor::writeImage(char* writeFileName, ProbImageType::Pointer image)
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

SomaExtractor::SegmentedImageType::Pointer SomaExtractor::SegmentSoma( ProbImageType::Pointer inputImage, std::vector< itk::Index<3> > &somaCentroids, 
													    int timethreshold, double curvatureScaling, double rmsThres)
{
	typedef itk::NearestNeighborInterpolateImageFunction< ProbImageType, float>  InterpolatorType;
	InterpolatorType::Pointer I_Interpolator = InterpolatorType::New();
	I_Interpolator->SetInputImage(inputImage);

	int SM = inputImage->GetLargestPossibleRegion().GetSize()[0];
    int SN = inputImage->GetLargestPossibleRegion().GetSize()[1];
    int SZ = inputImage->GetLargestPossibleRegion().GetSize()[2];

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

SomaExtractor::SegmentedImageType::Pointer SomaExtractor::SegmentSoma( ProbImageType::Pointer floatImage, std::vector< itk::Index<3> > &somaCentroids, 
													    int timethreshold, double curvatureScaling, double rmsThres, int minObjSize, unsigned int volThres, bool bnucleusTable)
{
	typedef itk::NearestNeighborInterpolateImageFunction< ProbImageType, float>  InterpolatorType;
	InterpolatorType::Pointer I_Interpolator = InterpolatorType::New();
	I_Interpolator->SetInputImage(floatImage);

	int SM = floatImage->GetLargestPossibleRegion().GetSize()[0];
	int SN = floatImage->GetLargestPossibleRegion().GetSize()[1];
    int SZ = floatImage->GetLargestPossibleRegion().GetSize()[2];

	//move the seed points along z axis
	PointList3D centroidsList;
    for( int i = 0; i < somaCentroids.size(); i++ )
	{
		Point3D point(somaCentroids[i][0], somaCentroids[i][1], somaCentroids[i][2], 0, i);
		point.check_out_of_range_3D(SM, SN, SZ);

		ProbImageType::IndexType index;
		ProbImageType::IndexType index1;
		index[0] = index1[0] = point.x;
		index[1] = index1[1] = point.y;
		index[2] = point.z;

		for( int j = 0; j < SZ; j++ )
		{
			index1[2] = j;
			if( I_Interpolator->EvaluateAtIndex(index1) > I_Interpolator->EvaluateAtIndex(index) )
			{
				somaCentroids[i][2] = j;
				index[2] = j;
				point.z = j;
			}
		}
		centroidsList.AddPt(point);
	}

	/// rescaled image as speed image
	std::cout << "RescaleIntensity"<<endl;
	typedef itk::RescaleIntensityImageFilter< ProbImageType, ProbImageType> RescaleFilterType;
    RescaleFilterType::Pointer rescale = RescaleFilterType::New();
	rescale->SetInput( floatImage); 
    rescale->SetOutputMinimum( 0 );
    rescale->SetOutputMaximum( 1 );
    rescale->Update();
	floatImage = rescale->GetOutput();

	typedef itk::ConnectedComponentImageFilter< SegmentedImageType, SegmentedImageType, SegmentedImageType> LabelFilterType;
	LabelFilterType::Pointer label = LabelFilterType::New();
	SegmentedImageType::Pointer binImage;

	if(bnucleusTable)
	{
		// eliminate the irrelevant seeds
		binImage = FastMarchingShapeDectectSoma(floatImage, centroidsList, timethreshold, 0.3, rmsThres);    // rough active contours

		std::cout<< "Remove small objects and build centroids queue "<<std::endl;
		relabelBinaryImage(binImage, centroidsList, minObjSize);
		std::cout<< std::endl;

		if(somaSeedQueue)
		{
			bool bfirst = true;
			GetCentroids(centroidsList);  /// get First Centroids from the queue
			while( centroidsList.GetSize() > 0)
			{
				std::cout<< "Trial Seed Points Number:"<< centroidsList.GetSize()<<std::endl;
				binImage = FastMarchingShapeDectectSoma(floatImage, centroidsList, timethreshold, curvatureScaling, rmsThres);

				/// label the image
				label->SetInput(binImage);
				label->Update();
				binImage = label->GetOutput();

				if( bfirst)
				{
					somaFeatureTable = ComputeSomaFeatures(binImage, centroidsList, true);   // build a new table
					bfirst = false;
				}
				else
				{
					somaFeatureTable = ComputeSomaFeatures(somaFeatureTable, binImage, centroidsList);    // attach the features to the table
				}
				GetCentroids(centroidsList);
				std::cout<< std::endl;
			}
			GetSomaCentroids( volThres, somaCentroidsList);
			binImage = FastMarchingShapeDectectSoma(floatImage, somaCentroidsList, timethreshold, curvatureScaling, rmsThres);
			somaImage = binImage;

			label->SetInput(binImage);
			label->Update();
			binImage = label->GetOutput();
			somaFeatureTable = ComputeSomaFeatures(binImage, somaCentroidsList, false); 
		}
		else
		{
			std::cout<< "Seed Queue hasn't been built successfully"<<endl;
		}
	}
	else 
	{
		binImage = FastMarchingShapeDectectSoma(floatImage, centroidsList, timethreshold, curvatureScaling, rmsThres);    // directly active contours based on the seeds 
		somaImage = binImage;

		label->SetInput(binImage);
		label->Update();
		binImage = label->GetOutput();
		somaFeatureTable = ComputeSomaFeatures(binImage, centroidsList, false); 
	}
	return somaImage;   // unlabeled
}

SomaExtractor::SegmentedImageType::Pointer SomaExtractor::FastMarchingShapeDectectSoma(ProbImageType::Pointer inputImage, PointList3D &seg_seeds, int timeThreshold, double curvatureScaling, double rmsError)
{
	/// generate distance map of the seeds using Fastmarching method
	std::cout<< "Generating Distance Map..." <<endl;
	clock_t SomaExtraction_start_time = clock();
	typedef  itk::FastMarchingImageFilter< ProbImageType, ProbImageType >    FastMarchingFilterType;
	FastMarchingFilterType::Pointer  fastMarching = FastMarchingFilterType::New();
	typedef FastMarchingFilterType::NodeContainer           NodeContainer;
    typedef FastMarchingFilterType::NodeType                NodeType;
    NodeContainer::Pointer seeds = NodeContainer::New();
    seeds->Initialize();

	for( int i = 0; i < seg_seeds.GetSize(); i++)
	{
		ProbImageType::IndexType  seedPosition;
		seedPosition[0] = seg_seeds.GetPt(i).x;
		seedPosition[1] = seg_seeds.GetPt(i).y;
		seedPosition[2] = seg_seeds.GetPt(i).z;
		NodeType node;
		const double seedValue = -3;
		node.SetValue( seedValue );
		node.SetIndex( seedPosition );
		seeds->InsertElement( i, node );
	}

	fastMarching->SetTrialPoints(  seeds);
    fastMarching->SetOutputSize( inputImage->GetBufferedRegion().GetSize());
    fastMarching->SetStoppingValue(  timeThreshold);
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
	shapeDetection->SetFeatureImage( inputImage);
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
	thresholder->Update();
	SegmentedImageType::Pointer binImage = thresholder->GetOutput();
	return binImage;
}

void SomaExtractor::relabelBinaryImage(SegmentedImageType::Pointer binImage, PointList3D &seg_seeds, int minObjSize)
{
	typedef itk::ConnectedComponentImageFilter< SegmentedImageType, SegmentedImageType, SegmentedImageType> LabelFilterType;
	typedef itk::RelabelComponentImageFilter< SegmentedImageType, SegmentedImageType > RelabelFilterType;

	LabelFilterType::Pointer label = LabelFilterType::New();
	label->SetInput(binImage);

	RelabelFilterType::Pointer relabel = RelabelFilterType::New();
	relabel->SetInput( label->GetOutput());
	relabel->SetMinimumObjectSize( minObjSize );  

	//Calculate labels
	try
    {
		relabel->Update();
		binImage = relabel->GetOutput();
    }
    catch( itk::ExceptionObject & excep )
    {
		std::cerr << "Relabel: soma extraction exception caught !" << std::endl;
		std::cerr << excep << std::endl;
    }
	
	std::cout << "Originally there were " << relabel->GetOriginalNumberOfObjects()<< " objects" << std::endl;
	std::cout << "After relabel there are now " << relabel->GetNumberOfObjects() << " objects" << std::endl;

	BuildCentroidQueue(binImage, seg_seeds);  /// build up centroid queues
}

void SomaExtractor::BuildCentroidQueue(SegmentedImageType::Pointer binImage, PointList3D &seg_seeds)
{
	if( somaSeedQueue == NULL)
	{
		somaSeedQueue = new PointList3D();
	}
	else
	{
		somaSeedQueue->RemoveAllPts();
	}

	for( int i = 0; i < seg_seeds.GetSize(); i++)
	{
		Point3D point = seg_seeds.GetPt(i);
		SegmentedImageType::IndexType index;
		index[0] = point.x;
		index[1] = point.y;
		index[2] = point.z;
		SegmentedImageType::PixelType pixel = binImage->GetPixel(index);
		if( pixel > 0)
		{
			point.tag = pixel;
			somaSeedQueue->AddPt(point);
		}
	}
	somaSeedQueue->BuildNeighbourList();
}

void SomaExtractor::GetCentroids(PointList3D &pointList)
{
	pointList.RemoveAllPts();
	if( somaSeedQueue)
	{
		somaSeedQueue->GetFirstPointOfEachQueue(pointList);
	}
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

vtkSmartPointer<vtkTable> SomaExtractor::ComputeSomaFeatures(SegmentedImageType::Pointer inputImage, PointList3D &seg_seeds, bool bTag)
{
	typedef itk::LabelGeometryImageFilter< SegmentedImageType> LabelGeometryImageFilterType;
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

	/// mapping from label to point tag
	std::map< TLPixel, int> LtoTagMap;
	if( bTag)
	{
		for( int i = 0; i < seg_seeds.GetSize(); i++)
		{
			Point3D point = seg_seeds.GetPt(i);
			SegmentedImageType::IndexType index;
			index[0] = point.x;
			index[1] = point.y;
			index[2] = point.z;
			SegmentedImageType::PixelType pixel = inputImage->GetPixel(index);
			if( pixel > 0)
			{
				LtoTagMap[pixel] = point.tag;
			}
		}
	}

	/// Calculating geometry features and build a vtkTable
	LabelGeometryImageFilterType::LabelsType::iterator allLabelsIt =  allLabels.begin();
	std::cout << "Calculating geometry labels, number of labels: " << labelGeometryImageFilter->GetNumberOfLabels() << std::endl;
	
	if( bTag)
	{
		for( int i = 0; i < volumeVec.size(); i++)
		{
			volumeVec[i].clear();
		}
		volumeVec.clear();
		volumeVec.resize(labelGeometryImageFilter->GetNumberOfLabels() - 1);
	}


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

		if( bTag)
		{
			int tag = LtoTagMap[labelValue];
			row->InsertNextValue(vtkVariant(tag));
			volumePoint pt(point[0], point[1], point[2], volume);
			volumeVec[ tag - 1].push_back( pt);
		}
		else
		{
			row->InsertNextValue(vtkVariant(count++));
		}
		
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

vtkSmartPointer<vtkTable> SomaExtractor::ComputeSomaFeatures(vtkSmartPointer<vtkTable> table, SegmentedImageType::Pointer inputImage, PointList3D &seg_seeds)
{
	typedef itk::LabelGeometryImageFilter< SegmentedImageType> LabelGeometryImageFilterType;
	LabelGeometryImageFilterType::Pointer labelGeometryImageFilter = LabelGeometryImageFilterType::New();
	labelGeometryImageFilter->SetInput( inputImage);

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

	/// mapping from label to point tag
	std::map< TLPixel, int> LtoTagMap;
	for( int i = 0; i < seg_seeds.GetSize(); i++)
	{
		Point3D point = seg_seeds.GetPt(i);
		SegmentedImageType::IndexType index;
		index[0] = point.x;
		index[1] = point.y;
		index[2] = point.z;
		SegmentedImageType::PixelType pixel = inputImage->GetPixel(index);
		if( pixel > 0)
		{
			LtoTagMap[pixel] = point.tag;
		}
	}

	/// Calculating geometry features and add to the vtkTable
	LabelGeometryImageFilterType::LabelsType::iterator allLabelsIt =  allLabels.begin();
	std::cout << "Calculating geometry labels, number of labels: " << labelGeometryImageFilter->GetNumberOfLabels() << std::endl;
 
	// remove the first labeled object which is the background
	for( allLabelsIt++; allLabelsIt != allLabels.end(); allLabelsIt++ )
    {
		LabelGeometryImageFilterType::LabelPixelType labelValue = *allLabelsIt;
		vtkSmartPointer<vtkVariantArray> row = vtkSmartPointer<vtkVariantArray>::New();
		int tag = LtoTagMap[labelValue];
		row->InsertNextValue(vtkVariant(tag));
		
		LabelGeometryImageFilterType::LabelPointType point = labelGeometryImageFilter->GetCentroid(labelValue);
		row->InsertNextValue(vtkVariant(point[0]));
		row->InsertNextValue(vtkVariant(point[1]));
		row->InsertNextValue(vtkVariant(point[2]));

		unsigned int volume = labelGeometryImageFilter->GetVolume(labelValue);

		volumePoint pt(point[0], point[1], point[2], volume);
		volumeVec[ tag - 1].push_back( pt);

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

vtkSmartPointer<vtkTable> SomaExtractor::GetSomaFeatureTable()
{
	return somaFeatureTable;
}

void SomaExtractor::GetSomaCentroids(unsigned int volumeThreshold, PointList3D &seg_seeds)
{
	seg_seeds.RemoveAllPts();
	for( int i = 0; i < volumeVec.size(); i++)
	{
		unsigned int volumeMax = 0;
		float x = 0;
		float y = 0;
		float z = 0;
		bool bflag = false;
		for( int j = 0; j < volumeVec[i].size(); j++)
		{
			volumePoint pt = volumeVec[i][j];
			if( pt.vol > volumeMax && pt.vol > volumeThreshold)
			{
				volumeMax = pt.vol;
				x = pt.x;
				y = pt.y;
				z = pt.z;
				bflag = true;
			}
		}
		if( bflag)
		{
			seg_seeds.AddPt(x,y,z);
		}
	}
}

void SomaExtractor::WriteSomaSeedsIntoImage()
{
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
    size[0] = inputImage->GetLargestPossibleRegion().GetSize()[0];  // size along X
    size[1] = inputImage->GetLargestPossibleRegion().GetSize()[0];  // size along Y
	size[2] = inputImage->GetLargestPossibleRegion().GetSize()[0];  // size along Z
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
	for(int i=0; i < somaCentroidsList.GetSize(); i++)
	{	
		Point3D seedIndex = somaCentroidsList.GetPt(i);
		OutputImageType::IndexType pixelIndex;

		pixelIndex[0] = seedIndex.x;
		pixelIndex[1] = seedIndex.y;
		pixelIndex[2] = seedIndex.z;
			
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
}