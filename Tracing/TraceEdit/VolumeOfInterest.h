
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
#include "itkImageFileWriter.h"
#include "itkImageToVTKImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkLabelGeometryImageFilter.h"

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
typedef itk::Image< float, 3> FloatImageType;

typedef itk::ImageFileReader< ImageType >    ReaderType;
typedef itk::ImageFileReader<FloatImageType> ReaderTypeFloat;
typedef itk::ImageFileWriter< ImageType >    WriterType;
typedef itk::ImageFileWriter< FloatImageType > FloatWriterType;
typedef itk::ImageToVTKImageFilter<ImageType> ConnectorType;
typedef itk::SignedMaurerDistanceMapImageFilter<ImageType, FloatImageType>			SignedMaurerDistanceMapImageFilterType;
typedef itk::DanielssonDistanceMapImageFilter<ImageType, FloatImageType, ImageType> VoronoiImageFilterType;
typedef itk::LabelGeometryImageFilter< ImageType >	LabelGeometryImageFilterType;

class VolumeOfInterest
{
public:
	VolumeOfInterest();
	int AddVOIPoint(double* newPT);
	bool ExtrudeVOI();
	vtkSmartPointer<vtkQuadricLODActor> GetActor();
	void CalculateCellDistanceToVOI(CellTraceModel *CellModel);
	float* CalculateCentroidDistanceToVOI(vtkSmartPointer<vtkTable> tbl);
	FloatImageType::Pointer GetVesselMaskDistanceMap();
	void ReadVesselDistanceMap(std::string fileName);
	void ReadImageTypeFloat3D(std::string fileName, FloatImageType::Pointer& data_ptr);
	void ReadBinaryVOI(std::string filename);
	void ReadVTPVOI(std::string filename);
	void ReadOBJVOI(std::string filename);
	void WriteVTPVOI(std::string filename);
	void ReadNucleiLabelImage(std::string filename);
	void CalculateVoronoiLabelImage();
	void GetVoronoiBoundingBox();
	void WriteVoronoiLabelImage(std::string filename);

	ImageType::RegionType GetVesselImageRegion() {return vesselImageRegion;}
	
private:
	std::vector< double * > ROIPoints;
	std::vector< vtkSmartPointer<vtkPolyData > > VOIPolyData;
	ImageType::Pointer vesselMaskImage;
	ImageType::RegionType vesselImageRegion;

	ImageType::Pointer nucleiLabelImage;
	ImageType::Pointer voronoiImage;
	FloatImageType::Pointer voronoiDistMapImage;
};
#endif