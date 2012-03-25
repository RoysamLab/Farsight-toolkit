#ifndef _SOMA_EXTRACTION_H
#define _SOMA_EXTRACTION_H

#include "yousef_core/yousef_seg.h"
#include <itkImage.h>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include <vtkTable.h>
#include <vtkSmartPointer.h>
#include <iostream>
#include <time.h>
#include <vector>
#include <map>
#include "Tracing/TheTracingSystem/TracingCore/PointOperation.h"

class SomaExtractor
{

public: 
	
	static const int Dim = 3;
	typedef unsigned int TLPixel;
	typedef itk::Image< unsigned char, Dim > OutputImageType;
	typedef itk::Image< float, Dim > ProbImageType;
	typedef itk::Image< TLPixel, Dim > SegmentedImageType;

	typedef itk::ImageFileReader< OutputImageType > ReaderType;
	typedef itk::ImageFileWriter< OutputImageType > WriterType;
	typedef itk::ImageFileWriter< SegmentedImageType > somaImageWriter;

	//: constructor
	SomaExtractor();
 	//: destructor
	~SomaExtractor();

	void SetInputImage(char * fileName);
	void SmoothByDiffusionImageFilter();
	ProbImageType::Pointer GetFloatInputImage();
	void ReadSeedpoints(char * fileName, std::vector< itk::Index<3> > &seedVec, bool bNucleusTable);

	/// return labeled image for somas
	SegmentedImageType::Pointer SegmentSoma( ProbImageType::Pointer inputImage, std::vector< itk::Index<3> > &somaCentroids, int timethreshold, double curvatureScaling, double rmsThres);
	
	/// return unlabeled image for somas( with capacity to select seeds based on the volumes)
	SegmentedImageType::Pointer SegmentSoma( ProbImageType::Pointer inputImage, std::vector< itk::Index<3> > &somaCentroids, 
											int timethreshold, double curvatureScaling, double rmsThres, int minObjSize, unsigned int volThres, bool bnucleusTable);
	void writeSomaImage(char* writeFileName);
	void writeImage(char* writeFileName, ProbImageType::Pointer image);
	vtkSmartPointer<vtkTable> ComputeSomaFeatures(SegmentedImageType::Pointer inputImage, PointList3D &seg_seeds, bool bTag = false);
	vtkSmartPointer<vtkTable> ComputeSomaFeatures(vtkSmartPointer<vtkTable> table, SegmentedImageType::Pointer inputImage, PointList3D &seg_seeds);
	vtkSmartPointer<vtkTable> GetSomaFeatureTable();
	void WriteSomaSeedsIntoImage();

protected:
	struct volumePoint
	{
		volumePoint( float X, float Y, float Z, unsigned int Vol)
		{
			x = X;
			y = Y;
			z = Z;
			vol = Vol;
		};
		float x;
		float y;
		float z;
		unsigned int vol;
	};

	void SomaBoundaryScan(SegmentedImageType::Pointer labelImage, std::map< TLPixel, int> &LtoIMap, std::vector< int> &boundaryPixSize);
	SegmentedImageType::Pointer FastMarchingShapeDectectSoma(ProbImageType::Pointer inputImage, PointList3D &seg_seeds, 
															int timeThreshold, double curvatureScaling, double rmsError);
	void relabelBinaryImage(SegmentedImageType::Pointer binImage, PointList3D &seg_seeds, int minObjSize);
	void BuildCentroidQueue(SegmentedImageType::Pointer image, PointList3D &seg_seeds);
	void GetCentroids(PointList3D &pointList);
	void GetSomaCentroids(unsigned int volumeThreshold, PointList3D &seg_seeds);

private:
	ProbImageType::Pointer inputImage;
	SegmentedImageType::Pointer somaImage;
	yousef_nucleus_seg *NucleusSeg;

	vector<Seed> seeds;
	PointList3D *somaSeedQueue;
	vtkSmartPointer<vtkTable> somaFeatureTable;
	vector< vector< volumePoint> > volumeVec;
	PointList3D somaCentroidsList;
};

#endif
