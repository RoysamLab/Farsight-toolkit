/* 
 * Copyright 2009 Rensselaer Polytechnic Institute
 * This program is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
 * for more details.
 * 
 * You should have received a copy of the GNU General Public License along 
 * with this program; if not, write to the Free Software Foundation, Inc., 
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */


#ifndef _ASTRO_TRACER_H_
#define _ASTRO_TRACER_H_

#include "itkTimeProbe.h"
#include "itkImage.h"
#include "itkArray.h"
#include "itkCovariantVector.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVector.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageLinearIteratorWithIndex.h"

#include "itkRescaleIntensityImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkShapedNeighborhoodIterator.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkSymmetricSecondRankTensor.h"
#include "itkSymmetricEigenAnalysis.h"
#include "itkMedianImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkApproximateSignedDistanceMapImageFilter.h"
#include "itkHessianToObjectnessMeasureImageFilter.h"
#include "itkMultiScaleHessianBasedMeasureImageFilter.h"
#include "itkHessian3DToVesselnessMeasureImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkMultiplyImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkPowImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkSquareImageFilter.h"
#include "itkLabelGeometryImageFilter.h"
#include "itkImageDuplicator.h"
//#include "itkCurvatureAnisotropicDiffusionImagefilter.h"

#include "vtkTable.h"
#include "vtkVariant.h"
#include "vtkVariantArray.h"
#include "vtkDoubleArray.h"
#include "vtkSmartPointer.h"
#include "vtkDelimitedTextWriter.h"

#include "vnl/vnl_math.h"

#include <ftkUtils.h>
#include "PatternAnalysis/activeLearning/mclr.h"
#include "Tracing/ftkVesselTracer/ftkVesselTracer.h"

#include <queue>
#include <algorithm>
#include <string>
#include <fstream>
#include <map>
#include <vector>
#include <time.h>
#include <sstream>
#include <numeric>
#include <math.h>

#ifdef _OPENMP
	#include <omp.h>
#endif

#define MAXVAL 100000.0f

typedef float PixelType;

class SWCNode_astro;
class HeapNode_astro;
class Comparison_astro;
//class AstroTracer;

class HeapNode_astro
{
public:
	
	itk::Index<3> ndx;
	PixelType KeyValue;
	
	HeapNode_astro(itk::Index<3> , PixelType);
	HeapNode_astro();
	bool operator==(const HeapNode_astro&);
};

class Comparison_astro
{
public:
	bool operator() (const HeapNode_astro* lhs, const HeapNode_astro* rhs) const  
	{
		return ((lhs->KeyValue) > (rhs->KeyValue));
	}
};

class RootPointFeatureVector{

public:
	HeapNode_astro node;	
	unsigned short int ID;
	double radius;
	float sphere_likelihood;
	double ballness;
	double plateness;
	double vesselness;
	PixelType intensity;
	PixelType meanIntensity;
	PixelType varianceIntensity;
	PixelType minIntensity;
	PixelType maxIntensity;
	double nucleusDistance;
	double radiusVBT;
	double likelihoodVBT;
	ODFFeatures odfFeatures;

	RootPointFeatureVector();
};

class CandidateRootPoint{

public:
	RootPointFeatureVector featureVector;
	
	int classValue;
	bool isRootPoint;
	double confidenceMeasure;

	CandidateRootPoint();
};

class IntrinsicFeatureVector{

public:

	HeapNode_astro centroid;
	//HeapNode_astro weightedCentroid;
	unsigned short int ID;
	double volume;
	double boundingBoxVolume;
	int integratedIntensity;
	double meanIntensity;
	double varianceIntensity;
	double minIntensity;
	double maxIntensity;
	double eccentricity;
	double elongation;
	double orientation;
	double meanSurfaceGradient;
	double interiorGradient;
	double surfaceIntensity;
	double interiorIntensity;
	double intensityRatio;
	double convexity;
	double radiusVariation;
	double shapeMeasure; // Ratio of surface voxels to total voxels
	double energy;
	double entropy;
	double inverseDiffMoment;
	double inertia;
	double clusterShade;
	double clusterProminence;
	double surfaceArea;
	double sharedBoundary;
	
	IntrinsicFeatureVector();
};

class AssociativeFeatureVector{

public:

	double astro_total;
	double astro_avg;
	double astro_surr;
	double micro_total;
	double micro_avg;
	double micro_surr;
	double neuro_total;
	double neuro_avg;
	double neuro_surr;
	double vessel_total;
	double vessel_avg;
	double vessel_surr;

	double minRootDist;
	double maxRootDist;
	double meanRootDist;
	double varRootDist;
	double nRoots;

	AssociativeFeatureVector();

};

class NucleiObject{

public:
	IntrinsicFeatureVector intrinsicFeatures;
	AssociativeFeatureVector associativeFeatures;
	
	// 1: Astrpcytes	2: Microglia	3:Neurons	4: Endotheliels
	int classValue; 
	double confidenceMeasure;

	double distanceToElectrode;

	NucleiObject();
};

class ObjectnessMeasures{
	
public:
	float alpha;
	float beta;
	float gamma;

	float sigma_min;
	float sigma_max;
	int sigma_intervals;
	int objectness_type; //0: Blobness, 1: Vesselness, 2: Plateness
	
	float noiseness;
	float ballness;
	float plateness;
	float vesselness;

	ObjectnessMeasures();
	ObjectnessMeasures(float alpha, float beta, float gamma);
	ObjectnessMeasures(float sigma_min, float sigma_max, float sigma_intervals, int obj_type);
};

class AstroTracer
{
public:
    
    typedef itk::SymmetricSecondRankTensor<double, 3> HessianPixelType;
    typedef itk::Image< HessianPixelType, 3 > HessianImageType3D;

	typedef itk::Index<3> IndexType;
	typedef itk::Offset<3> OffsetType;
	typedef itk::Image< PixelType, 3 >  ImageType3D;
	typedef itk::ImageFileReader<ImageType3D> ReaderType;
	typedef itk::RescaleIntensityImageFilter<ImageType3D, ImageType3D> RescalerType;
	typedef itk::MedianImageFilter<ImageType3D, ImageType3D> MedianFilterType;
	typedef itk::Image< unsigned char, 3 > CharImageType3D;
	typedef itk::Image< SWCNode_astro*, 3 > SWCImageType3D;
	typedef itk::Image< unsigned short, 3 > LabelImageType3D;
	typedef LabelImageType3D::PixelType * LabelArrayType;

	typedef itk::LaplacianRecursiveGaussianImageFilter< ImageType3D , ImageType3D> LoGFilterType;
	typedef itk::RegionOfInterestImageFilter<ImageType3D, ImageType3D> VolumeOfInterestFilterType;
	typedef itk::RegionOfInterestImageFilter<LabelImageType3D, LabelImageType3D> VolumeOfInterestFilterType_nuclei;
	typedef itk::StatisticsImageFilter<ImageType3D> StatisticsFilterType;
	typedef itk::SignedMaurerDistanceMapImageFilter<CharImageType3D, ImageType3D> SignedMaurerDistanceMapImageFilterType;
	
    typedef itk::MultiScaleHessianBasedMeasureImageFilter<ImageType3D, HessianImageType3D > MultiScaleHessianFilterType;
    typedef itk::HessianToObjectnessMeasureImageFilter<HessianImageType3D, ImageType3D > ObjectnessFilterType;

	typedef itk::MinimumMaximumImageCalculator<ImageType3D> MinMaxCalculatorType;
	typedef itk::MultiplyImageFilter<ImageType3D> MultiplyImageFilter;
	typedef itk::AddImageFilter<ImageType3D> AddImageFilter;
	typedef itk::PowImageFilter<ImageType3D> PowImageFilter;
	typedef itk::LabelGeometryImageFilter<LabelImageType3D, CharImageType3D> LabelGeometryFilterType;
	typedef itk::SubtractImageFilter<ImageType3D> SubtractImageFilterType;
	typedef itk::SquareImageFilter<ImageType3D, ImageType3D> SquareImageFilterType;
	typedef itk::ImageDuplicator<ImageType3D> DuplicatorType;
	//typedef itk::CurvatureAnisotropicDiffusionImageFilter<ImageType3D, ImageType3D> CurvatureAnisotropicDiffusionFilterType;
	
	//Constructor
	AstroTracer();
	//Destructor
	//~AstroTracer();

	void LoadParameters(const char* parametersFileName);

	void LoadCurvImageFromPath(std::string fname, unsigned int pad); 
	void LoadCurvImage(ImageType3D::Pointer &image, unsigned int pad);
	void ReadStartPointsFromPath(std::string fname, unsigned int padz);
	void ReadStartPoints(std::vector< itk::Index<3> > somaCentroids, unsigned int padz);
	void ReadStartPointsFromVTKTable(const std::vector<vtkSmartPointer<vtkTable> >);
	void SetCostThreshold(float thres){CostThreshold = thres;};
	void LoadSomaImage(std::string somaFileName);
	void RunTracing();
	void WriteMultipleSWCFiles(std::string fname, unsigned int );	
	void WriteSWCFile(std::string , unsigned int );
	void GenerateTestImage(); 
	//void LoadCurvImage_1(ImageType3D::Pointer &image, unsigned int pad);//
	//void ReadStartPoints_1(std::vector< itk::Index<3> > somaCentroids, unsigned int pad);//

	void ComputeAstroFeatures(std::string, std::string, unsigned int, const std::string);
	void ComputeAstroFeaturesPipeline(std::string, std::string, unsigned int, ImageType3D::RegionType, std::vector<vtkSmartPointer<vtkTable> >&, std::vector<LabelImageType3D::Pointer>&, const bool);
	bool PopulateLoGImages(std::vector<float> sigma_vec);
	void CallFeatureMainExternal();

	void SetNScales(int);
	void SetScaleRange(int, int);
	void SetInputDataPath(const std::string);

	void UseActiveLearningRootsModel(std::string);
	void ReadRootPointsExternal(std::string);
	void ReadRootPointsPipeline(const std::vector<vtkSmartPointer<vtkTable> >);
	void GetCentroidsForTracing(std::string, std::string);
	void GetCentroidsForTracingPipeline(std::string, std::string, unsigned int, ImageType3D::RegionType, std::vector<vtkSmartPointer<vtkTable> >&, std::vector<LabelImageType3D::Pointer>&, const bool);
	void ReadNucleiFeaturesExternal(std::string);
	void ReadNucleiFeaturesPipeline(const std::vector<vtkSmartPointer<vtkTable> >);
	void ComputeFeaturesFromCandidateRoots(void);
	void ComputeFeaturesFromCandidateRootsPipeline(const ImageType3D::RegionType, std::vector<vtkSmartPointer<vtkTable> >&, std::string, std::string); 
	void WriteNucleiFeatures(std::string);
	void ReadFinalNucleiTable(std::string);
	void ComputeObjectnessImage(ObjectnessMeasures obj_measures);
	void ComputeFTKObjectnessImage(void);
	void OptimizeCoverage(std::string coverageFileName, bool);
	void ReadStartPointsInternal(void);
	int optionsCreate(const char* optfile, std::map<std::string,std::string>& options);
	void Classification_Roots(std::vector< vtkSmartPointer<vtkTable> >&, std::vector< LabelImageType3D::Pointer >&, std::string, std::string, std::string, const bool, bool normalize_from_model = false);
	
	void Set_DistanceMapImage(ImageType3D::Pointer distance_map_image);
	void Set_NucleiLabelImage(LabelImageType3D::Pointer nuc_label_image);
	
	void DoPreprocessing(void);
	void LoadPreprocessingResults(void);
	void ComputeSomaDistanceMap(void);
	void ComputeRootPointFeatures();
	
	//external parameters
	float intensity_threshold;
	float contrast_threshold;
	int cost_threshold;
	int offshoot;
		
protected:
	void FeatureMainExternal();
	void FeatureMain();
	void GetFeatureExternal(float sigma, int scale_index); //void GetFeature( float );
	void GetFeature(float sigma); //void GetFeature( float );
	bool IsSeed(const itk::FixedArray<float, 3> & , unsigned int & );
	bool RegisterIndex(const float, itk::Index<3> &, itk::Size<3> &, long);
	SWCNode_astro* TBack(itk::Index<3> & ndx, std::vector<IndexType> &  );
	float GetCost(SWCNode_astro*, itk::Index<3>&);
	float GetCostLocal(SWCNode_astro* , itk::Index<3> & );
	float GetCostLocalLabel(SWCNode_astro* , itk::Index<3> & );
	float GetCostOld(SWCNode_astro* , itk::Index<3> &  );

	void ScanNeighbors( PixelType & a1,PixelType & a2,PixelType & a3, itk::Index<3> &);
	PixelType Update( PixelType a1,  PixelType a2,  PixelType a3,   PixelType P ) ;
	void Decimate();
	void Interpolate(PixelType);
	void RemoveIntraSomaNodes();	
	float getRadius(itk::Vector<float,3> & pos);
	float getRadiusAndLikelihood(itk::Vector<float,3> & pos, float& likelihood);
	void getHessianEigenFeatures(itk::Index<3> current_idx, double img_max_val, float& ballness, float& plateness, float& vesselness, float& noiseness);
	bool getLocalIntensityFeatures(itk::Index<3> current_idx, itk::Size<3> sz, float radius, float& max_intensity, float& min_intensity, float& mean_intensity, float& var_intensity);
	VBTNode getVBTFeatures(itk::Index<3> current_idx);
	void WriteImage3D(std::string , ImageType3D::Pointer );
	void BlackOut(itk::Index<3> &ndx );
	float GetCostLocal2(SWCNode_astro*, itk::Index<3>&);
	float GetCostLocal3(SWCNode_astro*, itk::Index<3>&);
	bool IsBall(const itk::FixedArray<float, 3>&, unsigned int&, double&);
	void GetHessianBasedObjectnessMeasures(itk::FixedArray<double, 3>&, ObjectnessMeasures&);
	

private:
	std::vector<SWCNode_astro*> SWCNode_astroContainer;
	//CharImageType3D::Pointer SomaImage;
	LabelImageType3D::Pointer SomaImage;
	PixelType CostThreshold;
	std::priority_queue < HeapNode_astro* , std::vector<HeapNode_astro*>,  Comparison_astro > PQ;
	ImageType3D::Pointer PaddedCurvImage, ConnImage, NDXImage, NDXImage2, NDXImage3;   //Input Image, EK image, CT image
	ImageType3D::Pointer LoGScaleImage;
	ImageType3D::Pointer ObjectnessImage;
	ImageType3D::Pointer ObjectnessHybridImage;
	ImageType3D::Pointer SomaDistanceMapImage;
	ImageType3D::Pointer AnisotropicDiffusedImage, gx, gy, gz, VesselnessImage;
	LabelImageType3D::Pointer IDImage, FinalRootsImage;	
	LabelImageType3D::Pointer RefinedRootImage;
	LabelImageType3D::Pointer NucleiLabelImage;
	SWCImageType3D::Pointer SWCImage; //swc label image
	itk::Size<3> size;
	std::vector<OffsetType> off;
	long CurrentID;
	std::vector<IndexType> StartPoints;
	unsigned int padz;

	vtkSmartPointer<vtkTable> features_table;
	vtkSmartPointer<vtkTable> roots_table;
	LabelImageType3D::Pointer roots_Image;

	int nScales;
	int startScale, endScale; // Perform LoG at these scales only
	bool isCoverageOptimized;

	std::vector<ImageType3D::Pointer> LoG_Vector;
	std::vector<std::vector<HeapNode_astro> > LoGPointsVector;
	std::vector<HeapNode_astro> AllLoGPointsVector;
	
	std::vector<CandidateRootPoint> CandidateRootPoints, AllRootPoints;
	std::vector<NucleiObject> NucleiObjects;

	std::vector<HeapNode_astro> CentroidListForTracing;
	std::vector<double> CentroidScales;
	std::vector<HeapNode_astro> DensityFilteredCentroidListForTracing;

	std::string InputDataPath;

	ftkVesselTracer *VBT;
};

///////////////////////////////////////////////////////////////

class SWCNode_astro 
{
public:
	long ID, PID, TreeID;
	itk::Index<3> ndx;
	itk::Vector<float,3> pos;
	bool IsLeaf, IsBranch, IsActive;
	SWCNode_astro *parent;
	std::vector<SWCNode_astro*> children;

	SWCNode_astro(); 
	SWCNode_astro(long, long, long, itk::Index<3> );
	SWCNode_astro(long, SWCNode_astro *, long, itk::Index<3> ); 
	
	//static bool IsIndexSame(itk::Index<3>, itk::Index<3>);

};


#endif

