#include "VolumeOfInterest.h"
VolumeOfInterest::VolumeOfInterest()
{
	this->VOIPolyData.clear();// =  std::vector<vtkPolyData*>;
	this->ROIPoints.clear();
}
int VolumeOfInterest::AddVOIPoint(double* newPT)
{
	this->ROIPoints.push_back(newPT);
	return (int) this->ROIPoints.size();
}
bool VolumeOfInterest::ExtrudeVOI()
{
	if (this->ROIPoints.size() < 3)
	{
		std::cout<< "not enough points\n";
		return false;
	}

	//// Add the points to a vtkPoints object
	vtkSmartPointer<vtkPoints> vtkROIpoints = vtkSmartPointer<vtkPoints>::New();
	
	std::vector<double*>::iterator ROIPoints_iter;

	int count = 0;
	for (ROIPoints_iter = this->ROIPoints.begin(); ROIPoints_iter != ROIPoints.end(); ROIPoints_iter++)
	{
		vtkROIpoints->InsertNextPoint(*ROIPoints_iter);
		//this->EditLogDisplay->append((QString("%1\t%2\t%3").arg((*ROIPoints_iter)[0]).arg((*ROIPoints_iter)[1]).arg((*ROIPoints_iter)[2])));
		//can do with string stream
		delete *ROIPoints_iter;
		*ROIPoints_iter = NULL;
		count++; //size of polygon vertex
	}

	// build a polygon 
	vtkSmartPointer<vtkPolygon> ROI_Poly = vtkSmartPointer<vtkPolygon>::New();
	ROI_Poly->GetPointIds()->SetNumberOfIds(count);
	for (int i =0; i< count; i++)
	{
		ROI_Poly->GetPointIds()->SetId(i,i);
	}
	
	//build cell array
	vtkSmartPointer<vtkCellArray> ROI_Poly_CellArray = vtkSmartPointer<vtkCellArray>::New();
	ROI_Poly_CellArray->InsertNextCell(ROI_Poly);

	// Create a 2d polydata to store outline in
	vtkSmartPointer<vtkPolyData> ROIpolydata = vtkSmartPointer<vtkPolyData>::New();
	ROIpolydata->SetPoints(vtkROIpoints);
	ROIpolydata->SetPolys(ROI_Poly_CellArray);
	
	//extrude the outline
	vtkSmartPointer<vtkLinearExtrusionFilter> extrude = vtkSmartPointer<vtkLinearExtrusionFilter>::New();
	extrude->SetInputData( ROIpolydata);
	extrude->SetExtrusionTypeToNormalExtrusion();
	extrude->SetScaleFactor (100); //adjust depending upon size of image stack
	extrude->Update();

	this->VOIPolyData.push_back( extrude->GetOutput());
	return true;
}
vtkSmartPointer<vtkQuadricLODActor> VolumeOfInterest::GetActor()
{
	vtkSmartPointer<vtkPolyDataMapper> VOImapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	VOImapper->SetInputData(this->VOIPolyData.back());
	vtkSmartPointer<vtkQuadricLODActor> VOIactor = vtkSmartPointer<vtkQuadricLODActor>::New();
	//VOIactor->GetProperty()->SetRepresentationToSurface();
	VOIactor->SetMapper(VOImapper);
	VOIactor->GetProperty()->SetInterpolationToFlat();
	return VOIactor;
}

void VolumeOfInterest::CalculateCellDistanceToVOI(CellTraceModel *CellModel)
{
	vtkSmartPointer<vtkCellLocator> cellLocator = vtkSmartPointer<vtkCellLocator>::New();
	cellLocator->SetDataSet(this->VOIPolyData.back());
	cellLocator->BuildLocator();

	std::map< int ,CellTrace*>::iterator cellCount = CellModel->GetCelliterator();
	for (; cellCount != CellModel->GetCelliteratorEnd(); cellCount++)
	{
		//double testPoint[3] = {500, 600, 50};
		double somaPoint[3];
		CellTrace* currCell = (*cellCount).second;
		currCell->getSomaCoord(somaPoint);
		//Find the closest points to TestPoint
		double closestPoint[3];//the coordinates of the closest point will be returned here
		double closestPointDist2; //the squared distance to the closest point will be returned here
		vtkIdType cellId; //the cell id of the cell containing the closest point will be returned here
		int subId; //this is rarely used (in triangle strips only, I believe)
		cellLocator->FindClosestPoint(somaPoint, closestPoint, cellId, subId, closestPointDist2);
		currCell->setDistanceToROI( std::sqrt(closestPointDist2), closestPoint[0], closestPoint[1], closestPoint[2]);
	}//end for cell count
}

float* VolumeOfInterest::CalculateCentroidDistanceToVOI(vtkSmartPointer<vtkTable> tbl)
{
	vtkSmartPointer<vtkCellLocator> cellLocator = vtkSmartPointer<vtkCellLocator>::New();
	cellLocator->SetDataSet(this->VOIPolyData.back());
	cellLocator->BuildLocator();

	float* dist_object = new float[(int)tbl->GetNumberOfRows()];
	for (int row=0; row < (int)tbl->GetNumberOfRows(); row++)
	{
		//double testPoint[3] = {500, 600, 50};
		double centroid[3];
		centroid[0] = tbl->GetValueByName(row,"centroid_x").ToDouble();
		centroid[1] = tbl->GetValueByName(row,"centroid_y").ToDouble();
		centroid[2] = tbl->GetValueByName(row,"centroid_z").ToDouble();
		//Find the closest points to TestPoint
		double closestPoint[3];//the coordinates of the closest point will be returned here
		double closestPointDist2; //the squared distance to the closest point will be returned here
		vtkIdType cellId; //the cell id of the cell containing the closest point will be returned here
		int subId; //this is rarely used (in triangle strips only, I believe)
		cellLocator->FindClosestPoint(centroid, closestPoint, cellId, subId, closestPointDist2);
		dist_object[row] = sqrt(closestPointDist2);		
	}//end for cell count

	return dist_object;
}

void VolumeOfInterest::ReadVesselDistanceMap(std::string fileName)
{
	ReaderType::Pointer vesselMaskReader = ReaderType::New();
	vesselMaskReader->SetFileName(fileName.c_str());	
	try
	{
		vesselMaskReader->Update();
	}
	catch( itk::ExceptionObject & exp )
	{
		std::cerr << "Exception thrown while reading the input file " << std::endl;
		std::cerr << exp << std::endl;
		//return EXIT_FAILURE;
	}
	vesselMaskImage = vesselMaskReader->GetOutput();
	
	this->vesselImageRegion = vesselMaskImage->GetLargestPossibleRegion();
	ImageType::SizeType size = this->vesselImageRegion.GetSize();
	std::cout << "Vessel image size: " << size << std::endl;
}

void VolumeOfInterest::ReadImageTypeFloat3D(std::string fileName, FloatImageType::Pointer& data_ptr){

	ReaderTypeFloat::Pointer image_reader = ReaderTypeFloat::New();
	image_reader->SetFileName(fileName);

	try{
		image_reader->Update();
	}
	catch( itk::ExceptionObject & exp ){
		std::cerr << "Exception thrown while reading the input file " << std::endl;
		std::cerr << exp << std::endl;
		//return EXIT_FAILURE;
	}
	
	data_ptr = image_reader->GetOutput();
}

FloatImageType::Pointer VolumeOfInterest::GetVesselMaskDistanceMap()
{
	typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> ThresholdFilterType;
	ThresholdFilterType::Pointer threshold_filter = ThresholdFilterType::New();
	threshold_filter->SetLowerThreshold(1);
	threshold_filter->SetInsideValue(255);
	threshold_filter->SetOutsideValue(0);
	threshold_filter->SetInput(this->vesselMaskImage);
	//threshold_filter->Update();

	SignedMaurerDistanceMapImageFilterType::Pointer MaurerFilter = SignedMaurerDistanceMapImageFilterType::New();
	MaurerFilter->SetInput(threshold_filter->GetOutput());
	MaurerFilter->SetSquaredDistance(false);
	MaurerFilter->SetUseImageSpacing(false);
	MaurerFilter->SetInsideIsPositive(false);
	MaurerFilter->Update();
   
	FloatImageType::Pointer distance_Map = MaurerFilter->GetOutput();
	return distance_Map;
}

void VolumeOfInterest::ReadBinaryVOI(std::string filename)
{
	ReaderType::Pointer contourReader = ReaderType::New();
	contourReader->SetFileName(filename.c_str());	
	try
	{
		contourReader->Update();
	}
	catch( itk::ExceptionObject & exp )
	{
		std::cerr << "Exception thrown while reading the input file " << std::endl;
		std::cerr << exp << std::endl;
		//return EXIT_FAILURE;
	}
	ConnectorType::Pointer connector = ConnectorType::New();
	connector->SetInput( contourReader->GetOutput() );

	vtkSmartPointer<vtkMarchingCubes> ContourFilter = vtkSmartPointer<vtkMarchingCubes>::New();
		ContourFilter->ComputeNormalsOff();
		ContourFilter->ComputeScalarsOff();
		ContourFilter->ComputeGradientsOff();
		ContourFilter->SetInputData(connector->GetOutput() );
		ContourFilter->SetNumberOfContours(1);
	ContourFilter->SetValue(0,1);

	vtkSmartPointer<vtkUnsignedCharArray> colorArray = vtkSmartPointer<vtkUnsignedCharArray>::New();
	colorArray->SetNumberOfComponents(3);
	unsigned char color[3] = {255,255,255};
	ContourFilter->Update();
	for (vtkIdType count = 0; count < ContourFilter->GetOutput()->GetNumberOfCells(); count++)
	{
		colorArray->InsertNextTupleValue(color);
	}
	ContourFilter->GetOutput()->GetCellData()->SetScalars(colorArray);
	this->VOIPolyData.push_back( ContourFilter->GetOutput());
}
void VolumeOfInterest::ReadVTPVOI(std::string filename)
{
	// Read volume data from VTK's .vtp file
	vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
	reader->SetFileName(filename.c_str());
	reader->Update();
	this->VOIPolyData.push_back( reader->GetOutput());
}
void VolumeOfInterest::ReadOBJVOI(std::string filename)
{
	// Read volume data from VTK's .vtp file
	vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
	reader->SetFileName(filename.c_str());
	reader->Update();
	this->VOIPolyData.push_back( reader->GetOutput());
}
void VolumeOfInterest::WriteVTPVOI(std::string filename)
{  
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName(filename.c_str());
	writer->SetInputData(this->VOIPolyData.back());
	writer->Write();
}

//Voronoi
void VolumeOfInterest::ReadNucleiLabelImage(std::string filename)
{
	//! Read the label image
	ReaderType::Pointer labelImageReader = ReaderType::New();
	labelImageReader->SetFileName(filename.c_str());	
	try
	{
		labelImageReader->Update();
	}
	catch( itk::ExceptionObject & exp )
	{
		std::cerr << "Exception thrown while reading the input file " << std::endl;
		std::cerr << exp << std::endl;
		//return EXIT_FAILURE;
	}
	nucleiLabelImage = labelImageReader->GetOutput();
}
void VolumeOfInterest::CalculateVoronoiLabelImage()
{
	//! Calculates the voronoi and its distance map
	VoronoiImageFilterType::Pointer voronoiFilter = VoronoiImageFilterType::New();
	voronoiFilter->SetInput(nucleiLabelImage);
	voronoiFilter->Update();
	voronoiImage = voronoiFilter->GetVoronoiMap();
	voronoiDistMapImage = voronoiFilter->GetDistanceMap();
}
void VolumeOfInterest::GetVoronoiBoundingBox()
{
	LabelGeometryImageFilterType::Pointer labelGeometryImageFilter = LabelGeometryImageFilterType::New();
	labelGeometryImageFilter->SetInput( voronoiImage );
	labelGeometryImageFilter->CalculateOrientedBoundingBoxOn();
	labelGeometryImageFilter->Update();

	LabelGeometryImageFilterType::LabelsType allLabels = labelGeometryImageFilter->GetLabels();
	LabelGeometryImageFilterType::LabelsType::iterator allLabelsIt;
	std::cout << "Number of labels: " << labelGeometryImageFilter->GetNumberOfLabels() << std::endl;
	std::cout << std::endl;

	for( allLabelsIt = allLabels.begin(); allLabelsIt != allLabels.end(); allLabelsIt++ )
	{
		LabelGeometryImageFilterType::LabelPixelType labelValue = *allLabelsIt;
		std::cout << "\tCentroid: " << labelGeometryImageFilter->GetOrientedBoundingBoxOrigin(labelValue) << std::endl;
		std::cout << "\tSize: " << labelGeometryImageFilter->GetOrientedBoundingBoxSize(labelValue) << std::endl;
		std::cout << "\tVertices 1: " << labelGeometryImageFilter->GetOrientedBoundingBoxVertices(labelValue)[0] << std::endl;
		std::cout << "\tVertices 2: " << labelGeometryImageFilter->GetOrientedBoundingBoxVertices(labelValue)[1] << std::endl;
		std::cout << "\tVertices 3: " << labelGeometryImageFilter->GetOrientedBoundingBoxVertices(labelValue)[2] << std::endl;
		std::cout << "\tVertices 4: " << labelGeometryImageFilter->GetOrientedBoundingBoxVertices(labelValue)[3] << std::endl;
		std::cout << "\tVolume: " << labelGeometryImageFilter->GetOrientedBoundingBoxVolume(labelValue) << std::endl;
		std::cout << "\tRegion: " << labelGeometryImageFilter->GetOrientedLabelImage(labelValue) << std::endl;

		//std::cout << "\tVolume: " << labelGeometryImageFilter->GetVolume(labelValue) << std::endl;
		//std::cout << "\tIntegrated Intensity: " << labelGeometryImageFilter->GetIntegratedIntensity(labelValue) << std::endl;
		//std::cout << "\tCentroid: " << labelGeometryImageFilter->GetCentroid(labelValue) << std::endl;
		//std::cout << "\tWeighted Centroid: " << labelGeometryImageFilter->GetWeightedCentroid(labelValue) << std::endl;
		//std::cout << "\tAxes Length: " << labelGeometryImageFilter->GetAxesLength(labelValue) << std::endl;
		//std::cout << "\tMajorAxisLength: " << labelGeometryImageFilter->GetMajorAxisLength(labelValue) << std::endl;
		//std::cout << "\tMinorAxisLength: " << labelGeometryImageFilter->GetMinorAxisLength(labelValue) << std::endl;
		//std::cout << "\tEccentricity: " << labelGeometryImageFilter->GetEccentricity(labelValue) << std::endl;
		//std::cout << "\tElongation: " << labelGeometryImageFilter->GetElongation(labelValue) << std::endl;
		//std::cout << "\tOrientation: " << labelGeometryImageFilter->GetOrientation(labelValue) << std::endl;
		//std::cout << "\tBounding box: " << labelGeometryImageFilter->GetBoundingBox(labelValue) << std::endl;

		std::cout << std::endl << std::endl;

	}
}
void VolumeOfInterest::WriteVoronoiLabelImage(std::string filename)
{
	//! Writes the voronoi label image and the distance map to file
	//Voronoi
	std::string saveVoronoiFile = filename.substr(0, filename.size()-4)+"_voronoi.tif";

	WriterType::Pointer voronoiImageWriter = WriterType::New();
	voronoiImageWriter->SetFileName(saveVoronoiFile);
	std::cout << "Set input..." << std::endl;
	voronoiImageWriter->SetInput( voronoiImage );
	
	std::cout << "Update writer" << std::endl;
	try {
		voronoiImageWriter->Update();
	}
	catch (itk::ExceptionObject & excp )
	{
		std::cerr << excp << std::endl;
	}

	//Distance Map
	//tif does not support float values
	std::string saveVoronoiDistMapFile = filename.substr(0, filename.size()-4)+"_voronoiDistMap.mhd";

	FloatWriterType::Pointer voronoiDistMapImageWriter = FloatWriterType::New();
	voronoiDistMapImageWriter->SetFileName(saveVoronoiDistMapFile);
	std::cout << "Set input..." << std::endl;
	voronoiDistMapImageWriter->SetInput( voronoiDistMapImage );
	
	std::cout << "Update writer" << std::endl;
	try {
		voronoiDistMapImageWriter->Update();
	}
	catch (itk::ExceptionObject & excp )
	{
		std::cerr << excp << std::endl;
	}
}