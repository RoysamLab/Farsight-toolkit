#include "VolumeOfInterest.h"
VolumeOfInterest::VolumeOfInterest()
{
	this->VOIPolyData =  vtkSmartPointer<vtkPolyData>::New();
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
	extrude->SetInput( ROIpolydata);
	extrude->SetExtrusionTypeToNormalExtrusion();
	extrude->SetScaleFactor (100);
	extrude->Update();

	this->VOIPolyData = extrude->GetOutput();
	return true;
}
vtkSmartPointer<vtkActor> VolumeOfInterest::GetActor()
{
	vtkSmartPointer<vtkPolyDataMapper> VOImapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	VOImapper->SetInput(this->VOIPolyData);
	vtkSmartPointer<vtkActor> VOIactor = vtkSmartPointer<vtkActor>::New();
	VOIactor->SetMapper(VOImapper);
	return VOIactor;
}

void VolumeOfInterest::CalculateCellDistanceToVOI(CellTraceModel *CellModel)
{
	vtkSmartPointer<vtkCellLocator> cellLocator = vtkSmartPointer<vtkCellLocator>::New();
	cellLocator->SetDataSet(this->VOIPolyData);
	cellLocator->BuildLocator();

	unsigned int cellCount= CellModel->getCellCount();
	for (unsigned int i = 0; i < cellCount; i++)
	{
		//double testPoint[3] = {500, 600, 50};
		double somaPoint[3];
		CellTrace* currCell = CellModel->GetCellAtNoSelection( i);
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