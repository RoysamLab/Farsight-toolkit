#ifndef _SOMA_EXTRACTION_H
#define _SOMA_EXTRACTION_H

#include "yousef_core/yousef_seg.h"
#include <itkImage.h>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkShapeDetectionLevelSetImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkFastMarchingImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkLabelGeometryImageFilter.h"
#include <itkRelabelComponentImageFilter.h>
#include <vtkTable.h>
#include <vtkSmartPointer.h>
#include <iostream>
#include <time.h>
#include <vector>
#include <map>
#include "Tracing/TheTracingSystem/TracingCore/PointOperation.h"
#include <itkHuangThresholdImageFilter.h>
#include "itkRegionOfInterestImageFilter.h"
#include <itkSigmoidImageFilter.h>
#include <itkExtractImageFilter.h>

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

	typedef itk::RegionOfInterestImageFilter< ProbImageType, ProbImageType> RegionOfInterestFilter;
	typedef itk::RescaleIntensityImageFilter< ProbImageType, ProbImageType> RescaleFilterType;
	typedef itk::FastMarchingImageFilter< ProbImageType, ProbImageType >    FastMarchingFilterType;
	typedef FastMarchingFilterType::NodeContainer           NodeContainer;
    typedef FastMarchingFilterType::NodeType                NodeType;
	typedef itk::ConnectedComponentImageFilter< SegmentedImageType, SegmentedImageType, SegmentedImageType> LabelFilterType;
	typedef itk::RelabelComponentImageFilter< SegmentedImageType, SegmentedImageType > RelabelFilterType;
	typedef itk::LabelGeometryImageFilter< SegmentedImageType> LabelGeometryImageFilterType;
	typedef itk::ShapeDetectionLevelSetImageFilter< ProbImageType, ProbImageType> ShapeDetectionFilterType;
	typedef itk::BinaryThresholdImageFilter< ProbImageType, SegmentedImageType> BinaryThresholdingFilterType;
	typedef itk::BinaryThresholdImageFilter< ProbImageType, ProbImageType> BinaryProbThresholdingFilterType;
	typedef itk::HuangThresholdImageFilter< ProbImageType, ProbImageType> HuangThresholdFilter;
	typedef itk::SigmoidImageFilter <ProbImageType, ProbImageType> SigmoidImageFilterType;

	//: constructor
	SomaExtractor();
 	//: destructor
	virtual ~SomaExtractor();

	void SetInputImage(const char * fileName);
	void SetInputImage( ProbImageType::Pointer probImage);
	ProbImageType::Pointer GetFloatInputImage();
	void ReadSeedpoints(const char * fileName, std::vector< itk::Index<3> > &seedVec, bool bNucleusTable);
	ProbImageType::Pointer EnhanceContrast( ProbImageType::Pointer inputImage, double alfa, double beta);
	/// return labeled image for somas
	SegmentedImageType::Pointer SegmentSoma( ProbImageType::Pointer input, std::vector< itk::Index<3> > &somaCentroids, double alfa, double beta, int timethreshold, double curvatureScaling, double rmsThres, int minObjSize);
	
	void writeImage(const char* writeFileName, SegmentedImageType::Pointer image);
	void writeImage(const char* writeFileName, ProbImageType::Pointer image);
	void writeCentroids(const char* writeFileName, std::vector< itk::Index<3> > &seedVec);
	
	vtkSmartPointer<vtkTable> ComputeSomaFeatures(SegmentedImageType::Pointer inputImage);

protected:
	void SomaBoundaryScan(SegmentedImageType::Pointer labelImage, std::map< TLPixel, int> &LtoIMap, std::vector< int> &boundaryPixSize);

private:
	ProbImageType::Pointer inputImage;
};

#endif
