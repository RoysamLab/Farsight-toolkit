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
	ProbImageType::Pointer GetFloatInputImage();
	void ReadSeedpoints(char * fileName, std::vector< itk::Index<3> > &seedVec);
	SegmentedImageType::Pointer SegmentSoma( ProbImageType::Pointer inputImage, std::vector< itk::Index<3> > &somaCentroids, 
											int timethreshold, double curvatureScaling, double rmsThres);
	void writeSomaImage(char* writeFileName);
	vtkSmartPointer<vtkTable> ComputeSomaFeatures(SegmentedImageType::Pointer inputImage);

protected:
	void BuildLabelObjectMap(SegmentedImageType::Pointer labelImage);
	void SomaBoundaryScan(SegmentedImageType::Pointer labelImage);

private:
	ProbImageType::Pointer binImage;
	SegmentedImageType::Pointer somaImage;
	yousef_nucleus_seg *NucleusSeg;
	vector<Seed> seeds;
	int SM,SN,SZ;
	std::vector< int> boundaryPix;	//boundary pixels for each label
	std::map< TLPixel, int> LtoIMap;	
};

#endif
