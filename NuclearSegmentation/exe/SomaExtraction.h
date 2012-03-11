#ifndef _SOMA_EXTRACTION_H
#define _SOMA_EXTRACTION_H

#include "yousef_core/yousef_seg.h"
#include <itkImage.h>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include <iostream>
#include <time.h>

class SomaExtractor
{

public: 
	
	static const int Dim = 3;

	typedef itk::Image< unsigned char, Dim > OutputImageType;
	typedef itk::Image< float, Dim > ProbImageType;
	typedef itk::Image< unsigned int, Dim > SegmentedImageType;

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

private:
	ProbImageType::Pointer binImage;
	SegmentedImageType::Pointer somaImage;
	yousef_nucleus_seg *NucleusSeg;
	vector<Seed> seeds;
	int SM,SN,SZ;
};

#endif
