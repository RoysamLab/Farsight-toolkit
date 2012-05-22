#include "ConvexHull3D.h"

ConvexHull3D::ConvexHull3D()
{
	this->convexHullArea = 0;
	this->convexHullVol = 0;
	refPt[0] = -1;
	refPt[1] = -1;
	refPt[2] = -1;
	cellCentroid[0] = -1;
	cellCentroid[1] = -1;
	cellCentroid[2] = -1;
}

void ConvexHull3D::setPoints(std::vector<TraceBit> &tips)
{
	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();

	vtkSmartPointer<vtkPoints> inputPoints = vtkSmartPointer<vtkPoints>::New();
	for(unsigned int counter=0; counter<tips.size(); counter++)
	{
		inputPoints->InsertNextPoint(tips[counter].x,tips[counter].y,tips[counter].z);
	}

	polydata->SetPoints(inputPoints);

	// Generate a tetrahedral mesh from the input points. By
	// default, the generated volume is the convex hull of the points.
	delaunay3D = vtkSmartPointer<vtkDelaunay3D>::New();
	delaunay3D->SetInput(polydata);
	delaunay3D->SetTolerance(0.01);
	delaunay3D->Update();

	//Surface filter is used to get the surfaces outside of the volume
	vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
	surfaceFilter->SetInputConnection(delaunay3D->GetOutputPort());
	surfaceFilter->Update();

	vtkSmartPointer<vtkDataSetMapper> delaunayMapper = vtkSmartPointer<vtkDataSetMapper>::New();
	delaunayMapper->SetInputConnection(surfaceFilter->GetOutputPort());

	delaunayActor = vtkSmartPointer<vtkActor>::New();
	delaunayActor->SetMapper(delaunayMapper);
	delaunayActor->GetProperty()->SetColor(1,0,0);

/*	vtkSmartPointer<vtkPoints> outputPoints = polyData->GetPoints();
	std::cout << "Num of output point: " << outputPoints->GetNumberOfPoints() << std::endl;*/

	//save output for area and volume calculation
	this->surfacePolyData = surfaceFilter->GetOutput();
}

void ConvexHull3D::setReferencePt(double point[3])
{
	this->refPt[0] = point[0];
	this->refPt[1] = point[1];
	this->refPt[2] = point[2];
	//std::cout << "Center: " << refPt[0] << " " << refPt[1] << " " << refPt[2] << std::endl;
}

bool ConvexHull3D::calculate()
{
	//area calculation
	vtkSmartPointer<vtkPoints> boundaryPoints = this->surfacePolyData->GetPoints();
	vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
	vtkSmartPointer<vtkTetra> tetra = vtkSmartPointer<vtkTetra>::New();

	//int count = 1;
	vtkCellArray * polys = this->surfacePolyData->GetPolys();
	vtkIdType nPts = 0;
	vtkIdType * ptIds = polys->GetPointer();
	for (polys->InitTraversal(); polys->GetNextCell(nPts,ptIds);)
	{
		//std::cout << "Count: " << count << std::endl;
		//count++;
		double vertex[3];
		double list [3][3]; // 3D triangle
		for (int j = 0; j < nPts; j++)
		{
			boundaryPoints->GetPoint(ptIds[j], vertex);
			list[j][0] = vertex[0];
			list[j][1] = vertex[1];
			list[j][2] = vertex[2];
			//std::cout << j << ": " << vertex[0] << " " << vertex[1] << " " << vertex[2] << std::endl;
		}
		//calculate area
		double point1[3],point2[3],point3[3];
		for (int i = 0; i<3; i++)
		{
			point1[i] = list[0][i];
			point2[i] = list[1][i];
			point3[i] = list[2][i];
		}
		this->convexHullArea += triangle->TriangleArea(point1,point2,point3);
		this->convexHullVol += abs(tetra->ComputeVolume(refPt,point1,point2,point3));
	}
		//std::cout << "Area: " << convexHullArea << std::endl;
		//std::cout << "Volume: " << convexHullVol << std::endl;

	////volume calculation
	//vtkSmartPointer<vtkTetra> tetra = vtkSmartPointer<vtkTetra>::New();

	//vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = delaunay3D->GetOutput();
	//vtkSmartPointer<vtkPoints> allPoints = unstructuredGrid->GetPoints();
	//vtkSmartPointer<vtkCellArray> cellTetras = unstructuredGrid->GetCells();

	//count = 1;
	//vtkIdType volNPts = 0;
	//vtkIdType * volPtIds = cellTetras->GetPointer();
	//for (cellTetras->InitTraversal(); cellTetras->GetNextCell(volNPts,volPtIds);)
	//{
	//	std::cout << "Count: " << count << std::endl;
	//	count++;
	//	double vertex[3];
	//	double list [4][3]; // tetrahedron
	//	for (int j = 0; j < volNPts; j++)
	//	{
	//		allPoints->GetPoint(volPtIds[j], vertex);
	//		list[j][0] = vertex[0];
	//		list[j][1] = vertex[1];
	//		list[j][2] = vertex[2];
	//		std::cout << j << ": " << vertex[0] << " " << vertex[1] << " " << vertex[2] << std::endl;
	//	}

	//	double point1[3],point2[3],point3[3],point4[3];
	//	for (int i = 0; i<3; i++)
	//	{
	//		point1[i] = list[0][i];
	//		point2[i] = list[1][i];
	//		point3[i] = list[2][i];
	//		point4[i] = list[3][i];
	//	}
	//	this->convexHullVol += tetra->ComputeVolume(point1,point2,point3,point4);
	//	std::cout << "Volume: " << convexHullVol << std::endl;
	//}

	//std::cout << "Number of surface points: " << surfacePolyData->GetNumberOfPoints() << std::endl;
	
	//calculate based on surface points only
	double totalX = 0; double totalY = 0; double totalZ = 0;
	for (vtkIdType i = 0; i < surfacePolyData->GetNumberOfPoints(); i++)
	{
		double point[3];
		surfacePolyData->GetPoint(i,point);
		totalX += point[0]-this->refPt[0];
		totalY += point[1]-this->refPt[1];
		totalZ += point[2]-this->refPt[2];
		//std::cout << "Point " << i << " : (" << point[0] << " " << point[1] << " " << point[2] << ")" << std::endl;
	}

	this->convexHullMagnitude = sqrt(pow(totalX,2)+pow(totalY,2)+pow(totalZ,2));
	this->convexHullAzimuth = atan2(totalY,totalX)*180/PI;
	double hypotenuse = sqrt(pow(totalX,2)+pow(totalY,2));
	this->convexHullElevation = atan2(totalZ,hypotenuse)*180/PI;

	int numOfPts = surfacePolyData->GetNumberOfPoints();
	this->cellCentroid[0] = totalX/numOfPts + this->refPt[0];
	this->cellCentroid[1] = totalY/numOfPts + this->refPt[1];
	this->cellCentroid[2] = totalZ/numOfPts + this->refPt[2];
	//std::cout << "Centroid: (" << cellCentroid[0] << " " << cellCentroid[1] << " " << cellCentroid[2] << ")" << std::endl;
	
	//add to table later

	return true;
}

vtkSmartPointer<vtkActor> ConvexHull3D::getActor()
{
	return delaunayActor;
}

//double ConvexHull3D::getVolume()
//{
//	
//	return this->convexHullVol;
//}