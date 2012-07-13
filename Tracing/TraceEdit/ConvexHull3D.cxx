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

void ConvexHull3D::setPoints(std::vector<TraceBit> &pts)
{
	/*!
	 * @author Audrey Cheong
	 * @param pts data points
	 */
	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();

	vtkSmartPointer<vtkPoints> inputPoints = vtkSmartPointer<vtkPoints>::New();
	for(unsigned int counter=0; counter<pts.size(); counter++)
	{
		inputPoints->InsertNextPoint(pts[counter].x,pts[counter].y,pts[counter].z);
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
	delaunayActor->GetProperty()->SetColor(1,0.0,0.0);

/*	vtkSmartPointer<vtkPoints> outputPoints = polyData->GetPoints();
	std::cout << "Num of output point: " << outputPoints->GetNumberOfPoints() << std::endl;*/

	//save output for area and volume calculation
	this->surfacePolyData = surfaceFilter->GetOutput();
}

void ConvexHull3D::setReferencePt(double point[3])
{
	/*!
	 * @author Audrey Cheong
	 * @param pts the reference point (e.g. soma center)
	 */
	this->refPt[0] = point[0];
	this->refPt[1] = point[1];
	this->refPt[2] = point[2];
	//std::cout << "Center: " << refPt[0] << " " << refPt[1] << " " << refPt[2] << std::endl;
}

bool ConvexHull3D::calculate()
{
	/*!
	 * Calculate area, volume, centroid, etc.
	 * @author Audrey Cheong
	 * @return check whether calculations is successful
	 */
	//area calculation
	vtkSmartPointer<vtkPoints> boundaryPoints = this->surfacePolyData->GetPoints();
	vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
	vtkSmartPointer<vtkTetra> tetra = vtkSmartPointer<vtkTetra>::New();

	vtkCellArray * polys = this->surfacePolyData->GetPolys();
	vtkIdType nPts = 0;
	vtkIdType * ptIds = polys->GetPointer();
	for (polys->InitTraversal(); polys->GetNextCell(nPts,ptIds);)
	{
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

	return true;
}

vtkSmartPointer<vtkActor> ConvexHull3D::getActor()
{
	/*!
	 * @author Audrey Cheong
	 * @return convex hull(vtkActor)
	 */
	return delaunayActor;
}

void ConvexHull3D::calculateEllipsoid()
{
	/*!
	 * Calculate best-fit 3D ellipse
	 * @author Audrey Cheong
	 * @return check whether calculations are successful
	 */
	//http://www.ahinson.com/algorithms/Sections/InterpolationRegression/EigenPlane.pdf
	// find best fit plane
	vnl_matrix<double> A(3,3, 0.0); //3x3 matrix, fill with zeroes

	double x_norm,y_norm,z_norm;
	double point[3];
	for (vtkIdType i = 0; i < surfacePolyData->GetNumberOfPoints(); i++)
	{
		//normalize points
		surfacePolyData->GetPoint(i,point);
		x_norm = point[0]-this->cellCentroid[0];
		y_norm = point[1]-this->cellCentroid[1];
		z_norm = point[2]-this->cellCentroid[2];

		//calculate covariance (symmetric matrix)
		A(0,0) += pow(x_norm,2);	A(0,1) += x_norm * y_norm;	A(0,2) +=  x_norm * z_norm;
									A(1,1) += pow(y_norm,2);	A(1,2) += y_norm * z_norm;
																A(2,2) += pow(z_norm,2);
	}
		A(1,0) = A(0,1);
		A(2,0) = A(0,2);
		A(2,1) = A(1,2);

	// get eigenvector corresponding to the smallest and largest eigenvalue (normal to plane)
	vnl_symmetric_eigensystem<double> eig(A);
	double eigenvalues[3];
	double eigenvalue_norm[3]; //length
	int num_of_points = (int) surfacePolyData->GetNumberOfPoints();

	//std::cout << "Number of points: " << num_of_points << std::endl;

	double min = eig.get_eigenvalue(0);
	double max = eig.get_eigenvalue(0);
	int min_index = 0;
	int max_index = 0;
	for (int i = 0; i < 3; i++)
	{
		eigenvalues[i] = eig.get_eigenvalue(i);
		eigenvalue_norm[i] = 4*sqrt(eigenvalues[i] / num_of_points); //normalize to get eigenvalue_norm
		if (eigenvalues[i] < min)
		{
			min = eigenvalues[i];
			min_index = i;
		}
		if (eigenvalues[i] >= max)
		{
			max = eigenvalues[i];
			max_index = i;
		}
	}

	int median_index = 3 - min_index - max_index;
	vnl_vector<double> eigenVector_normal = eig.get_eigenvector(min_index);		//normal axis (smallest)
	vnl_vector<double> eigenVector_minor = eig.get_eigenvector(median_index);	//minor axis
	vnl_vector<double> eigenVector_major = eig.get_eigenvector(max_index);		//major axis
	eigenVector_normal.normalize();
	eigenVector_minor.normalize();
	eigenVector_major.normalize();

	//std::cout << "EigenVector: ";

	//how to orientate correctly?
	vtkMatrix4x4 * matrix = vtkMatrix4x4::New();
	for (int i = 0; i < 3; i++)
	{
		vnl_vector<double> eigenVector = eig.get_eigenvector(i);

		for (int j = 0; j < 3; j++)
		{
			matrix->SetElement(i,j,eigenVector.get(j));
		}
		////rotation
		//matrix->SetElement(0,i,eigenVector_normal.get(i));
		//matrix->SetElement(1,i,eigenVector_minor.get(i));
		//matrix->SetElement(2,i,eigenVector_major.get(i));
		//position
		matrix->SetElement(i,3,cellCentroid[i]+200);

		//std::cout << eigenVector_normal.get(i) << " ";
	}

	//std::cout << std::endl;

	//double x_rotate = eigenVector_normal.get(0);
	//double y_rotate = eigenVector_normal.get(1);
	//double z_rotate = eigenVector_normal.get(2);


	//std::cout << "Eigenvalues: " << eigenvalue_norm[0] << " " << eigenvalue_norm[1] << " " << eigenvalue_norm[2] << std::endl;
	//still need to add to table!!

	////validate (draw ellipsoid)

	vtkSmartPointer<vtkParametricEllipsoid> parametricObject = vtkSmartPointer<vtkParametricEllipsoid>::New();
	//parametricObject->SetXRadius(eigenvalue_norm[min_index]);
	//parametricObject->SetYRadius(eigenvalue_norm[median_index]);
	//parametricObject->SetZRadius(eigenvalue_norm[max_index]);
	parametricObject->SetXRadius(eigenvalue_norm[0]);
	parametricObject->SetYRadius(eigenvalue_norm[1]);
	parametricObject->SetZRadius(eigenvalue_norm[2]);
	vtkSmartPointer<vtkParametricFunctionSource> parametricFunctionSource = vtkSmartPointer<vtkParametricFunctionSource>::New();
	parametricFunctionSource->SetParametricFunction(parametricObject);
	parametricFunctionSource->Update();

	// Setup mapper and actor
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(parametricFunctionSource->GetOutputPort());

	ellipsoidActor = vtkSmartPointer<vtkActor>::New();
	ellipsoidActor->SetMapper(mapper);
	ellipsoidActor->SetUserMatrix(matrix);
	//ellipsoidActor->SetPosition(cellCentroid);
	//ellipsoidActor->RotateX(180);
	//ellipsoidActor->RotateY(-90);
	//ellipsoidActor->RotateZ(-90);
}

//for validation purposes (unless desired for visualization)
vtkSmartPointer<vtkActor> ConvexHull3D::get3DEllipseActor()
{
	/*!
	 * @author Audrey Cheong
	 * @return ellipsoid(vtkActor)
	 */
	std::cout<< "Ellipsoid actor" << std::endl;
	return ellipsoidActor;
}
