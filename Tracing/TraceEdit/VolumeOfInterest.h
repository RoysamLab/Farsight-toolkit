
#ifndef VOLUMEOFINTEREST_H_
#define VOLUMEOFINTEREST_H_

#include <stdio.h>
#include <string>

#include "vtkActor.h"
#include "vtkQuadricLODActor.h"
#include "vtkSmartPointer.h"

#include "vtkCellArray.h"
#include "vtkCellLocator.h"
#include "vtkCellData.h"

#include "itkImageFileReader.h"
#include "itkImageToVTKImageFilter.h"

#include "vtkContourFilter.h"
#include "vtkMarchingCubes.h"

#include "vtkVolume.h"
#include "vtkVolumeProperty.h"

#include "vtkPoints.h"
#include "vtkLinearExtrusionFilter.h"

#include "vtkPolygon.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"

#include "vtkXMLPolyDataReader.h"
#include "vtkXMLPolyDataWriter.h"
#include <vtkOBJReader.h>

//FTK Includes
#include "CellTrace.h"
#include "CellTraceModel.h"

typedef itk::Image< unsigned char, 3 >   ImageType;

typedef itk::ImageFileReader< ImageType >    ReaderType;
typedef itk::ImageToVTKImageFilter<ImageType> ConnectorType;

class VolumeOfInterest
{
public:
	VolumeOfInterest();
	int AddVOIPoint(double* newPT);
	bool ExtrudeVOI();
	vtkSmartPointer<vtkQuadricLODActor> GetActor();
	void CalculateCellDistanceToVOI(CellTraceModel *CellModel);
	float* CalculateCentroidDistanceToVOI(vtkSmartPointer<vtkTable> tbl);
	void ReadBinaryVOI(std::string filename);
	void ReadVTPVOI(std::string filename);
	void ReadOBJVOI(std::string filename);
	void WriteVTPVOI(std::string filename);
private:
	std::vector< double * > ROIPoints;
	std::vector< vtkSmartPointer<vtkPolyData > > VOIPolyData;
};
#endif