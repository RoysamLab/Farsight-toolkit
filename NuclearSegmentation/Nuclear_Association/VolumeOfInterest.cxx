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
	extrude->SetScaleFactor (100);
	extrude->Update();

	this->VOIPolyData.push_back( extrude->GetOutput());
	return true;
}
vtkSmartPointer<vtkActor> VolumeOfInterest::GetActor()
{
	vtkSmartPointer<vtkPolyDataMapper> VOImapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	VOImapper->SetInputData(this->VOIPolyData.back());
	vtkSmartPointer<vtkActor> VOIactor = vtkSmartPointer<vtkActor>::New();
	//VOIactor->GetProperty()->SetRepresentationToSurface();
	VOIactor->SetMapper(VOImapper);
	return VOIactor;
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
void VolumeOfInterest::WriteVTPVOI(std::string filename)
{  
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName(filename.c_str());
	writer->SetInputData(this->VOIPolyData.back());
	writer->Write();
}