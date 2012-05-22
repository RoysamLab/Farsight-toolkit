#ifndef CONVEXHULL3D_H
#define CONVEXHULL3D_H

#define PI 3.14159265
#include <cmath>

#include "TraceBit.h"

#include "vtkSmartPointer.h"
#include "vtkActor.h"
#include "vtkCellArray.h"
#include "vtkDataSetMapper.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkDelaunay3D.h"
#include "vtkIdList.h"
#include "vtkMath.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkProperty.h"
#include "vtkTetra.h"
#include "vtkTriangle.h"
#include "vtkUnstructuredGrid.h"

// by Audrey Cheong

class ConvexHull3D
{
public:
	ConvexHull3D();

	/**
	 * performs a convexhull in 3D
	 * @param	points the Vector3D cloud
	 * @return 
	 */
	void setPoints(std::vector<TraceBit> &points);
	void setReferencePt(double point[3]); //default is (0,0,0)
	bool calculate();
	vtkSmartPointer<vtkActor> getActor();
	//double[] getCellCentroid();
	//double getArea();
	//double getVolume();
	//static TraceBit getCentroid(std::vector<TraceBit> &points, int index, Face * face);

private:
	double convexHullMagnitude, convexHullAzimuth, convexHullElevation, convexHullArea, convexHullVol;
	double refPt[3];
	double cellCentroid[3];
	vtkSmartPointer<vtkActor> delaunayActor;
	vtkSmartPointer<vtkDelaunay3D> delaunay3D;
	vtkPolyData * surfacePolyData;
};
#endif