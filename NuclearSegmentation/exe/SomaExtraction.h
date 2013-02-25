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
#include <itkDivideImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkStatisticsImageFilter.h>
#include <itkThresholdImageFilter.h>
#include "itkMedianImageFilter.h"
#include <itkOtsuThresholdImageFilter.h>

class SomaExtractor
{

public: 
	
	static const int Dim = 3;
	typedef unsigned int TLPixel;
	typedef unsigned short UShortPixel;
	typedef itk::Image< unsigned char, 3 > OutputImageType;
	typedef itk::Image< float, 3 > ProbImageType;
	typedef itk::Image< float, 2 > ProbImageType2D;
	typedef itk::Image< TLPixel,3 > SegmentedImageType;
	typedef itk::Image< UShortPixel,3 > UShortImageType;
	typedef itk::Image< UShortPixel,2 > UShortImageType2D;

protected:
	typedef itk::ImageFileReader< OutputImageType > ReaderType;
	typedef itk::ImageFileWriter< OutputImageType > WriterType;
	typedef itk::ImageFileReader< SegmentedImageType > somaImageReader;
	typedef itk::ImageFileWriter< SegmentedImageType > somaImageWriter;
	typedef itk::ImageFileReader< UShortImageType > ushortImageReader;
	typedef itk::ImageFileWriter< UShortImageType > ushortImageWriter;
	typedef itk::ImageFileReader< ProbImageType > probImageReader;
	typedef itk::ImageFileReader< ProbImageType2D > probImage2DReader;
	typedef itk::ImageFileWriter< ProbImageType2D > probImage2DWriter;
	typedef itk::ImageFileReader< UShortImageType2D > ushortImage2DReader;
	typedef itk::ImageFileWriter< UShortImageType2D > ushortImage2DWriter;

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
	typedef itk::HuangThresholdImageFilter< ProbImageType, OutputImageType> HuangThresholdFilter;
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
	typedef itk::SubtractImageFilter <ProbImageType2D, ProbImageType2D, ProbImageType2D> SubtractImageFilterType2D;
	typedef itk::AddImageFilter <ProbImageType, ProbImageType, ProbImageType> AddImageFilterType;
	typedef itk::AddImageFilter <ProbImageType2D, ProbImageType2D, ProbImageType2D> AddImageFilterType2D;
	typedef itk::DivideImageFilter <ProbImageType, ProbImageType, ProbImageType > DivideImageFilterType;
	typedef itk::DivideImageFilter <ProbImageType2D, ProbImageType2D, ProbImageType2D > DivideImageFilterType2D;
	typedef itk::RegionOfInterestImageFilter< SegmentedImageType, SegmentedImageType> RegionOfInterestFilter;
	typedef itk::MinimumMaximumImageCalculator <ProbImageType> ImageCalculatorFilterType;
	typedef itk::ImageRegionConstIterator< ProbImageType > ProbConstIteratorType;
	typedef itk::ImageRegionConstIterator< ProbImageType2D > Prob2DConstIteratorType;
	typedef itk::ImageRegionIterator< ProbImageType>  ProbIteratorType;
	typedef itk::ImageRegionIterator< OutputImageType>  UcharIteratorType;
	typedef itk::ImageRegionIterator< ProbImageType2D > Prob2DIteratorType;
	typedef itk::StatisticsImageFilter<ProbImageType> StatisticsImageFilterType;
	typedef itk::StatisticsImageFilter<ProbImageType2D> StatisticsImageFilterType2D;
	typedef itk::ThresholdImageFilter< ProbImageType2D> ThresholdImageFilterType;
	typedef itk::MedianImageFilter<ProbImageType2D, ProbImageType2D> MedianFilterType;
	typedef itk::MedianImageFilter<SegmentedImageType, SegmentedImageType> SegMedianFilterType;
	typedef itk::OtsuThresholdImageFilter <OutputImageType, ProbImageType> OtsuThresholdImageFilterType;

public:
	//: constructor
	explicit SomaExtractor();
 	//: destructor
	~SomaExtractor();


	ProbImageType::Pointer SetInputImage(const char * fileName);
	ProbImageType::Pointer SetInputImage8bit(const char * fileName);
	ProbImageType2D::Pointer SetInputImage2D(const char * fileName);
	ProbImageType::Pointer SetInputImageByPortion(const char * fileName);
	SegmentedImageType::Pointer SetInitalContourImage(const char * fileName);
	SegmentedImageType::Pointer SetInitalContourImageByPortion(const char * fileName);
	SegmentedImageType::Pointer SetInitalContourImage16bit(const char * fileName);
	OutputImageType::Pointer Read8BitImage(const char * fileName);
	OutputImageType::Pointer Read16BitImage(const char * fileName);
	void ReadSeedpoints(const char * fileName, std::vector< itk::Index<3> > &seedVec, bool bNucleusTable);
	void LoadOptions(const char* paramFileName);
	ProbImageType2D::Pointer SetInputImageFloat2D(const char *fileName);
	ProbImageType::Pointer SetInputImageFloat(const char *fileName);

	/// return labeled image for somas
	SegmentedImageType::Pointer SegmentSoma( std::vector< itk::Index<3> > &somaCentroids, ProbImageType::Pointer binImagePtr);
	SegmentedImageType::Pointer SegmentSomaUsingGradient( ProbImageType::Pointer input, SegmentedImageType::Pointer initialContour, std::vector< itk::Index<3> > &somaCentroids);
	SegmentedImageType::Pointer SegmentHeart(const char *imageName, const char *fileName, ProbImageType::Pointer inputImage, vnl_vector<int> &seperator, vnl_vector<double> &curvature);

	ProbImageType::Pointer OtsuThresholdImage(OutputImageType::Pointer image);
	void writeImage(const char* writeFileName, SegmentedImageType::Pointer image);
	void writeImage(const char* writeFileName, OutputImageType::Pointer image);
	void writeImage(const char* writeFileName, ProbImageType::Pointer image, bool bscale = false);
	void writeImage(const char* writeFileName, GradientImageType::Pointer image);
	void writeImage(const char* writeFileName, UShortImageType::Pointer image);
	void WriteFloat2DImage(const char* writeFileName, ProbImageType2D::Pointer image);
	void writeCentroids(const char* writeFileName, std::vector< itk::Index<3> > &seedVec);
	
	vtkSmartPointer<vtkTable> ComputeSomaFeatures(SegmentedImageType::Pointer inputImage);
	vtkSmartPointer<vtkTable> ComputeHeartFeatures(SegmentedImageType::Pointer inputImage);
	void GetDebrisCentroids( OutputImageType::Pointer inputImage, std::vector< itk::Index<3> > &debrisSeeds);
	void AssociateDebris(OutputImageType::Pointer inputImage, std::vector< itk::Index<3> > &somaCentroids, std::vector< itk::Index<3> > &debrisSeeds);

	// generate soma seed points for the input image
	ProbImageType::Pointer GenerateSeedPoints(OutputImageType::Pointer inputImgPt, std::vector< itk::Index<3> > &somaCentroids);
 	ProbImageType::Pointer GenerateSeedPoints( unsigned char* inputBuffer, int size1, int size2, int size3, std::vector< itk::Index<3> > &somaCentroids);
	
	ProbImageType2D::Pointer GetAverage(const char * channelName, int n, double sigma);
	ProbImageType2D::Pointer GetBackgroundImage(ProbImageType::Pointer image, double sigma);
	ProbImageType2D::Pointer GetBackgroundImageByFirstSlice(ProbImageType::Pointer image, double sigma);
	UShortImageType::Pointer DevideAndScaleToOriginalMean(ProbImageType::Pointer image, ProbImageType2D::Pointer backgroundImage, int border);
	UShortImageType::Pointer DevideAndScale(ProbImageType::Pointer oriImage, ProbImageType2D::Pointer backgroundImage, double median, double ratio_threshold = 0);
	UShortImageType::Pointer RescaleImage(ProbImageType::Pointer image, double globalMax, double intensityMax);
	ProbImageType2D::Pointer adjustMeanStd(ProbImageType2D::Pointer image, double globalMean, double globalStd);
	void NormalizeUsingBackgroundImage(ProbImageType2D::Pointer image, ProbImageType2D::Pointer backgroundimage, double sigma);
	void writeUnshort2D(const char *fileName, UShortImageType2D::Pointer image);
	void GetSeedpointsInRegion(std::vector< itk::Index<3> > &seedVec, std::vector< itk::Index<3> > &seedInRegion, int startX, int startY, int width, int height);
	void CaculateMeanStd(std::string fileName, ProbImageType::Pointer image);
	ProbImageType::Pointer diffusionsmoothing(ProbImageType::Pointer image);

protected:
	template <class T> bool SetParamValue(std::map<std::string,std::string> &opts, std::string str, T &value, T defVal);
	void SomaBoundaryScan(SegmentedImageType::Pointer labelImage, std::map< TLPixel, int> &LtoIMap, std::vector< int> &boundaryPixSize);
	ProbImageType::Pointer GetEdgePotentialMap(ProbImageType::Pointer inputImage, double sigma);
	ProbImageType::Pointer GetInitalContourByDistanceMap(SegmentedImageType::Pointer labelImage, double outlierExpand);
	void CheckBoundary(SegmentedImageType::IndexType &start, SegmentedImageType::IndexType &end, int SX, int SY, int SZ); 
	ProbImageType::Pointer EnhanceContrast( ProbImageType::Pointer inputImage, int sliceNum, double alfa, double beta, double &threshold);
	ProbImageType::Pointer EnhanceContrast( ProbImageType::Pointer inputImage, double alfa, double beta, double radius);
	void ShrinkPixel(ProbImageType2D::Pointer image, int border);
	ProbImageType2D::Pointer ExtractSlice(ProbImageType::Pointer image, int sliceId);
	ProbImageType::Pointer RemoveImageBorderByPixel(ProbImageType::Pointer image, int border);

private:
	//ProbImageType::Pointer inputImage;
	int width;
	int height;
	int depth;

	// speed image
	double alfa;
	double beta;
	double open_radius;
	double fill_radius;
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
	unsigned int smoothIteration;
	double conductance;
	// read portion of image
	int startX;
	int startY;
	int startZ;
	int sizeX;
	int sizeY;
	int sizeZ;

	int num_bins;
	int shift;
	double scaleMin;
	double scaleMax;
	double regionXY;
	double regionZ;
	int useDistMap;
	int sampling_ratio_XY_to_Z;
	int radius;
	int brerun;
};

#endif
