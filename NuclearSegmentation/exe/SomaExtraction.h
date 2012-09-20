#ifndef _SOMA_EXTRACTION_H
#define _SOMA_EXTRACTION_H

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
#include "itkVotingBinaryHoleFillingImageFilter.h"
#include "itkAdaptiveHistogramEqualizationImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"

#include "itkVector.h"
#include "itkLaplacianSharpeningImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkGradientImageFilter.h"
#include "itkGradientVectorFlowImageFilter.h"
#include "itkGeodesicActiveContourLevelSetImageFilter.h"
#include "itkSignedDanielssonDistanceMapImageFilter.h"
#include <itkSubtractImageFilter.h>
#include "itkMinimumMaximumImageCalculator.h"
#include "itkImageRegionIterator.h"
#include <itkSignedMaurerDistanceMapImageFilter.h>

class SomaExtractor
{

public: 
	
	static const int Dim = 3;
	typedef unsigned int TLPixel;
	typedef unsigned short UShortPixel;
	typedef itk::Image< unsigned char, 3 > OutputImageType;
	typedef itk::Image< float, 3 > ProbImageType;
	typedef itk::Image< TLPixel,3 > SegmentedImageType;
	typedef itk::Image< UShortPixel,3 > UShortImageType;

protected:
	typedef itk::ImageFileReader< OutputImageType > ReaderType;
	typedef itk::ImageFileWriter< OutputImageType > WriterType;
	typedef itk::ImageFileReader< SegmentedImageType > somaImageReader;
	typedef itk::ImageFileWriter< SegmentedImageType > somaImageWriter;
	typedef itk::ImageFileReader< UShortImageType > ushortImageReader;
	typedef itk::ImageFileWriter< UShortImageType > ushortImageWriter;
	typedef itk::ImageFileReader< ProbImageType > probImageReader;

	//typedef itk::RegionOfInterestImageFilter< ProbImageType, ProbImageType> RegionOfInterestFilter;
	typedef itk::RescaleIntensityImageFilter< ProbImageType, ProbImageType> RescaleFloatFilterType;
	typedef itk::RescaleIntensityImageFilter< OutputImageType, OutputImageType> RescaleUCharFilterType;
	typedef itk::RescaleIntensityImageFilter< UShortImageType, OutputImageType> RescaleUshortFilterType;
	typedef itk::FastMarchingImageFilter< ProbImageType, ProbImageType >    FastMarchingFilterType;
	typedef FastMarchingFilterType::NodeContainer           NodeContainer;
    typedef FastMarchingFilterType::NodeType                NodeType;
	typedef itk::ConnectedComponentImageFilter< SegmentedImageType, SegmentedImageType, SegmentedImageType> LabelFilterType;
	typedef itk::RelabelComponentImageFilter< SegmentedImageType, SegmentedImageType > RelabelFilterType;
	typedef itk::LabelGeometryImageFilter< SegmentedImageType> LabelGeometryImageFilterType;
	typedef itk::ShapeDetectionLevelSetImageFilter< ProbImageType, ProbImageType> ShapeDetectionFilterType;
	typedef itk::BinaryThresholdImageFilter< ProbImageType, SegmentedImageType> BinaryThresholdingFilterType;
	typedef itk::BinaryThresholdImageFilter< ProbImageType, ProbImageType> BinaryProbThresholdingFilterType;
	typedef itk::BinaryThresholdImageFilter< SegmentedImageType, SegmentedImageType> SegThresholdingFilterType;
	typedef itk::ExtractImageFilter< ProbImageType, ProbImageType> ExtractFilterType;
	typedef itk::HuangThresholdImageFilter< ProbImageType, ProbImageType> HuangThresholdFilter;
	typedef itk::SigmoidImageFilter < ProbImageType, ProbImageType> SigmoidImageFilterType;
	typedef itk::VotingBinaryHoleFillingImageFilter< SegmentedImageType, SegmentedImageType> HoleFillingFilterType;
	typedef itk::AdaptiveHistogramEqualizationImageFilter< ProbImageType> AdaptiveHistogramEqualizationImageFilterType;
	typedef itk::BinaryBallStructuringElement< SegmentedImageType::PixelType, Dim> KernelType;
	typedef itk::BinaryMorphologicalClosingImageFilter< SegmentedImageType, SegmentedImageType, KernelType > CloseFilterType;

	typedef itk::LaplacianSharpeningImageFilter<ProbImageType, ProbImageType >  LaplacianSharpeningImageFilterType;
	typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<ProbImageType, ProbImageType > GradientFilterType;
	typedef itk::CovariantVector<float, Dim> GradientType;
	typedef itk::Image<GradientType, Dim>   GradientImageType;
	typedef itk::ImageFileWriter< GradientImageType>  VectorImageWriterType;
	typedef itk::GradientImageFilter<ProbImageType, float, float> GradientIFilterType;
	typedef itk::GradientVectorFlowImageFilter<GradientImageType, GradientImageType> GVFFilterType;
	typedef itk::GeodesicActiveContourLevelSetImageFilter< ProbImageType, ProbImageType > GeodesicActiveContourFilterType;
	typedef GeodesicActiveContourFilterType::VectorImageType  AdvectionImageType;
	typedef itk::CastImageFilter<GradientImageType,AdvectionImageType> CastFlowFilterType;
	typedef itk::SignedDanielssonDistanceMapImageFilter<ProbImageType, ProbImageType> DanielssonDistanceMapFilterType;
	typedef itk::SignedMaurerDistanceMapImageFilter<ProbImageType, ProbImageType> MaurerDistanceMapFilterType;
	typedef itk::SubtractImageFilter <ProbImageType, ProbImageType, ProbImageType> SubtractImageFilterType;
	typedef itk::RegionOfInterestImageFilter< SegmentedImageType, SegmentedImageType> RegionOfInterestFilter;
	typedef itk::MinimumMaximumImageCalculator <ProbImageType> ImageCalculatorFilterType;
	typedef itk::ImageRegionConstIterator< ProbImageType > ProbConstIteratorType;
	typedef itk::ImageRegionIterator< ProbImageType>  ProbIteratorType;

public:
	//: constructor
	SomaExtractor();
 	//: destructor
	virtual ~SomaExtractor();

	ProbImageType::Pointer SetInputImage(const char * fileName);
	ProbImageType::Pointer SetInputImageByPortion(const char * fileName);
	SegmentedImageType::Pointer SetInitalContourImage(const char * fileName);
	SegmentedImageType::Pointer SetInitalContourImageByPortion(const char * fileName);
	OutputImageType::Pointer Read8BitImage(const char * fileName);
	OutputImageType::Pointer Read16BitImage(const char * fileName);
	void ReadSeedpoints(const char * fileName, std::vector< itk::Index<3> > &seedVec, bool bNucleusTable);
	void LoadOptions(const char* paramFileName);

	/// return labeled image for somas
	SegmentedImageType::Pointer SegmentSoma( OutputImageType::Pointer input, std::vector< itk::Index<3> > &somaCentroids, ProbImageType::Pointer binImagePtr);
	SegmentedImageType::Pointer SegmentSoma( ProbImageType::Pointer input, SegmentedImageType::Pointer initialContour, std::vector< itk::Index<3> > &somaCentroids);

	void writeImage(const char* writeFileName, SegmentedImageType::Pointer image);
	void writeImage(const char* writeFileName, OutputImageType::Pointer image);
	void writeImage(const char* writeFileName, ProbImageType::Pointer image, bool bscale = false);
	void writeImage(const char* writeFileName, GradientImageType::Pointer image);
	void writeCentroids(const char* writeFileName, std::vector< itk::Index<3> > &seedVec);
	
	vtkSmartPointer<vtkTable> ComputeSomaFeatures(SegmentedImageType::Pointer inputImage);
	void GetDebrisCentroids( OutputImageType::Pointer inputImage, std::vector< itk::Index<3> > &debrisSeeds);
	void AssociateDebris(OutputImageType::Pointer inputImage, std::vector< itk::Index<3> > &somaCentroids, std::vector< itk::Index<3> > &debrisSeeds);

	// generate soma seed points for the input image
	ProbImageType::Pointer GenerateSeedPoints(OutputImageType::Pointer inputImgPt, std::vector< itk::Index<3> > &somaCentroids);

protected:
	template <class T> bool SetParamValue(std::map<std::string,std::string> &opts, std::string str, T &value, T defVal);
	void SomaBoundaryScan(SegmentedImageType::Pointer labelImage, std::map< TLPixel, int> &LtoIMap, std::vector< int> &boundaryPixSize);
	ProbImageType::Pointer GetEdgePotentialMap(ProbImageType::Pointer inputImage, double sigma);
	ProbImageType::Pointer GetInitalContourByDistanceMap(SegmentedImageType::Pointer labelImage, double outlierExpand);
	void CheckBoundary(SegmentedImageType::IndexType &start, SegmentedImageType::IndexType &end, int SX, int SY, int SZ); 
	ProbImageType::Pointer EnhanceContrast( ProbImageType::Pointer inputImage, int sliceNum, double alfa, double beta, double &threshold);
	ProbImageType::Pointer EnhanceContrast( ProbImageType::Pointer inputImage, double alfa, double beta, double radius);

private:
	//ProbImageType::Pointer inputImage;
	int width;
	int height;
	int depth;

	// speed image
	double alfa;
	double beta;
	// for initial contour
	int timethreshold; 
	double seedValue;
	// active contour
	double curvatureScaling; 
	double advectScaling;
	double rmsThres;
	int maxIterations;
	int holeSize;
	int minObjSize;
	// GVF active contour
	double sigma;
	unsigned int numberOfIterations; 
	double noiseLevel;
	double outlierExpandValue;
	
	// read portion of image
	int startX;
	int startY;
	int startZ;
	int sizeX;
	int sizeY;
	int sizeZ;

	// for automatic soma seed detection
	ProbImageType::Pointer binProbImagePtr;

	int num_bins;
	int shift;
	double scaleMin;
	double scaleMax;
	double regionXY;
	double regionZ;
	int useDistMap;
	int sampling_ratio_XY_to_Z;
};

#endif
