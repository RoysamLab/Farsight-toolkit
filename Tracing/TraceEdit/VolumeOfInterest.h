
#ifndef VOLUMEOFINTEREST_H_
#define VOLUMEOFINTEREST_H_

#include <stdio.h>
#include <string>

#include "vtkActor.h"
#include "vtkSmartPointer.h"

#include "vtkCellArray.h"
#include "vtkCellLocator.h"

#include "vtkVolume.h"
#include "vtkVolumeProperty.h"

#include "vtkPoints.h"
#include "vtkLinearExtrusionFilter.h"

#include "vtkPolygon.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"

//FTK Includes
#include "CellTrace.h"
#include "CellTraceModel.h"

class VolumeOfInterest
{
public:
	VolumeOfInterest();
	int AddVOIPoint(double* newPT);
	bool ExtrudeVOI();
	vtkSmartPointer<vtkActor> GetActor();
	void CalculateCellDistanceToVOI(CellTraceModel *CellModel);
private:
	std::vector<double*> ROIPoints;
	vtkSmartPointer<vtkPolyData> VOIPolyData;
};
#endif