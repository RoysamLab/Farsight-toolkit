#ifndef _SOMA_EXTRACTION_H
#define _SOMA_EXTRACTION_H

#include "yousef_core/yousef_seg.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileWriter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkGrayscaleMorphologicalOpeningImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkLabelGeometryImageFilter.h"

//label object classes
#include "itkShapeLabelObject.h"
#include "itkShapeLabelMapFilter.h" 
#include "itkBinaryImageToShapeLabelMapFilter.h"
#include <tinyxml/tinyxml.h>

#include <iostream>
#include <time.h>
#include "Tracing/TheTracingSystem/TracingCore/PointOperation.h"

class SomaExtractor
{

public: 

	//: constructor
	SomaExtractor();
 	//: destructor
	~SomaExtractor();

	static const int Dim = 3;


	typedef itk::Image< unsigned char, Dim > OutputImageType;
	typedef itk::Image< float, Dim > ProbImageType;
	typedef std::vector<OutputImageType::IndexType> centroidVectorType;
	typedef itk::ImageFileReader< OutputImageType > ReaderType;
	typedef itk::ImageRegionConstIterator< OutputImageType > ConstIteratorType;
	typedef itk::Image< int, Dim > SegmentedImageType;
	typedef itk::ImageRegionIteratorWithIndex< SegmentedImageType > IteratorType;
	typedef itk::BinaryBallStructuringElement< unsigned int, Dim > KernelType;
	//typedef itk::FlatStructuringElement< Dim > KernelType;
	typedef itk::GrayscaleMorphologicalOpeningImageFilter< SegmentedImageType, SegmentedImageType, KernelType > morphOpenFilterType;
	typedef itk::RelabelComponentImageFilter< SegmentedImageType, SegmentedImageType > RelabelFilterType;
	typedef itk::BinaryThresholdImageFilter<SegmentedImageType, OutputImageType> BinaryThresholdImageType;
	typedef unsigned long LabelType;
	typedef itk::ShapeLabelObject< LabelType, Dim > LabelObjectType;
	typedef itk::LabelMap< LabelObjectType > LabelMapType;
	typedef itk::BinaryImageToShapeLabelMapFilter< OutputImageType, LabelMapType >	ConverterType;
	typedef itk::ImageFileWriter< OutputImageType > WriterType;
	typedef itk::ImageFileWriter< SegmentedImageType > somaImageWriter;


	void SetInputImage(char * fileName);
	void SetInputImage(OutputImageType::Pointer img);
	bool LoadSegParams(int kernel, int minObj);
	int binarizeImage( char* paramFile, unsigned short num_bins);
	OutputImageType::Pointer GetSomaBinaryImage() { return outputImage; };
	int relabelBinaryImage(void);
	centroidVectorType GetSomaCentroids();
	void writeSomaCentroids(char* writeFileName);
	void writeSomaImage(char* writeFileName);

	void SegmentSoma(char* centroidsFileName, const char * somafileName, int timethreshold, double curvatureScaling, double rmsThres);
	int GenerateSeedPoints(char* paramFile, unsigned short num_bins);
	void ImFastMarching_Soma(PointList3D seg_seeds, int timeThreshold, double curvatureScaling, double rmsError, const char *somaFileName);

private:
	OutputImageType::Pointer inputImage, outputImage;
	SegmentedImageType::Pointer binImage, somaImage;
	centroidVectorType Centroids;
	unsigned char *in_Image;
	//ConstIteratorType pix_buf;
	yousef_nucleus_seg * NucleusSeg;
	size_t size1, size2, size3;
	morphOpenFilterType::Pointer morphOpenFilter;
	int minObjSize, KernelSize;
	vector<Seed> seeds;
	int SM,SN,SZ;
};

#endif
