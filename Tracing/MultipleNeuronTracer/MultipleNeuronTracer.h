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


#ifndef _MULTIPLE_NEURON_TRACER_H_
#define _MULTIPLE_NEURON_TRACER_H_

#include <stdlib.h>
#include "itkTimeProbe.h"
#include "PointOperation.h"
#include "ftkUtils.h"
#include "itkImage.h"
#include "itkArray.h"
#include "itkCovariantVector.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkStatisticsImageFilter.h"
//#include "itkHuangThresholdImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"

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
#include "itkHessianToObjectnessMeasureImageFilter.h"
#include "itkMultiScaleHessianBasedMeasureImageFilter.h"
#include "itkHessian3DToVesselnessMeasureImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkMultiplyImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkPowImageFilter.h"

#include "vnl/vnl_math.h"
#include "vtkSmartPointer.h"
#include "vtkTable.h"
#include "vtkVariantArray.h"
#include "vtkDoubleArray.h"



#include "Vesselness/itkMultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter.h"
#include "itkHessian3DToVesselnessMeasureImageFilter.h"
#include "itkMultiScaleHessianBasedMeasureImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkGradientVectorFlowImageFilter.h"
#include "itkGradientImageFilter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "itkCovariantVector.h"
#include "itkImageSliceConstIteratorWithIndex.h"
#include "itkImageSliceIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkImageAdaptor.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkLabelGeometryImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"

#include "itkRegionOfInterestImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkListSample.h"
#include <queue>
#include <algorithm>
#include <string>
#include <fstream>
#include <map>
#include <vector>
#include <time.h>

#define MAXVAL 100000.0f
#ifdef _OPENMP
    #include "omp.h"
#endif

typedef float PixelType;

typedef itk::CovariantVector< float, 3 >  GradientPixelType;
typedef itk::Image< GradientPixelType, 3 > GradientImageType;
typedef itk::Vector< float, 3 > SeedType;
typedef itk::Statistics::ListSample< SeedType > SampleType;
typedef SampleType::Pointer SamplePointer;

typedef GradientImageType::Pointer GradientImagePointer;
typedef itk::Image< PixelType, 3 >  ImageType3D;
typedef ImageType3D::Pointer ProbImagePointer;
typedef itk::Image<unsigned long int, 3> LabelImageTypeCCIF;//////////for connected component image filter in seed_centroid()
typedef LabelImageTypeCCIF::Pointer LabelImagePointerCCIF;

ProbImagePointer extract_one_component(int index, GradientImagePointer IG);
static void tred2(double V[3][3], double d[3], double e[3]);
static void tql2(double V[3][3], double d[3], double e[3]);
void eigen_decomposition(double A[3][3], double V[3][3], double d[3]);
void eigen_decomposition2D(double A[2][2], double V[2][2], double d[2]);
double MAX_1(double a, double b);
static double hypot2(double x, double y);
double absd(double val);



class SWCNode;
class DebrisNode;
class HeapNode;
class Comparison;
//class MultipleNeuronTracer;
class VectorPixelAccessor;
class HeapNode
{
public:
	itk::Index<3> ndx;
	PixelType KeyValue;
	HeapNode(itk::Index<3> , PixelType);
};

class Comparison
{
public:
	bool operator() (const HeapNode* lhs, const HeapNode* rhs) const  
	{
		return ((lhs->KeyValue) > (rhs->KeyValue));
	}
};

class ObjectnessMeasures_micro{
	
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
	//automatically selected threshold value

	ObjectnessMeasures_micro();
	ObjectnessMeasures_micro(float alpha, float beta, float gamma);
	ObjectnessMeasures_micro(float sigma_min, float sigma_max, float sigma_intervals, int obj_type);
};

class MultipleNeuronTracer
{
public:

	typedef itk::Index<3> IndexType;
	typedef itk::Offset<3> OffsetType;
	typedef itk::Image< PixelType, 3 >  ImageType3D;
	typedef itk::ImageFileReader<ImageType3D> ReaderType;
	typedef itk::RescaleIntensityImageFilter<ImageType3D, ImageType3D> RescalerType;
	typedef itk::MedianImageFilter<ImageType3D, ImageType3D> MedianFilterType;
	typedef itk::Image< unsigned char, 3 > CharImageType3D;
	typedef itk::RescaleIntensityImageFilter<CharImageType3D, CharImageType3D> CharRescalerType;
	typedef itk::Image< SWCNode*, 3 > SWCImageType3D;
	typedef itk::Image< unsigned int, 3 > LabelImageType3D;
	typedef LabelImageType3D::Pointer LabelImagePointer;
	typedef LabelImageType3D::PixelType * LabelArrayType;
	//typedef itk::HuangThresholdImageFilter<ImageType3D,ImageType3D> HuangThresholdFilterType;
	typedef itk::OtsuThresholdImageFilter<ImageType3D,ImageType3D>  OtsuThresholdImageFilterType;
	typedef itk::BinaryThresholdImageFilter <ImageType3D,ImageType3D>	BinaryThresholdImageFilterType;
	typedef itk::MinimumMaximumImageCalculator<ImageType3D> MinMaxImageCalculatorType;
	typedef itk::StatisticsImageFilter<ImageType3D> StatisticsFilterType;
	//typedef itk::HessianToObjectnessMeasureImageFilter<PixelType, 3> ObjectnessFilterType;
	//typedef itk::MultiScaleHessianBasedMeasureImageFilter<ImageType3D, ObjectnessFilterType> MultiScaleHessianFilterType;
	typedef itk::MinimumMaximumImageCalculator<ImageType3D> MinMaxCalculatorType;
	typedef itk::MultiplyImageFilter<ImageType3D> MultiplyImageFilter;
	typedef itk::AddImageFilter<ImageType3D> AddImageFilter;
	typedef itk::PowImageFilter<ImageType3D> PowImageFilter;




	//Constructor
	MultipleNeuronTracer();
	//Destructor
	~MultipleNeuronTracer();

	void LoadParameters(const char* parametersFileName, int _argc);
	void LoadParameters_1(const char* parametersFileName,float intensityThreshold,float contrastThreshold,int costThreshold);

	void LoadCurvImage(std::string fname, unsigned int pad); 
	void LoadCurvImage_1(ImageType3D::Pointer &image, unsigned int pad);
	void LoadCurvImage_2(ImageType3D::Pointer &image);
	void RunMask();
	
	void ReadStartPoints(std::string fname, unsigned int padz);
	void ReadStartPoints_1(std::vector< itk::Index<3> > somaCentroids, unsigned int padz);
	void ReadStartPoints_2(std::string fname, unsigned int pad,float startx,float starty,float startz, float widthx,float widthy,float widthz); //
	void SetCostThreshold(float thres){_CostThreshold = thres;};
	void LoadSomaImage(std::string somaFileName);
	void LoadSomaImage_1(LabelImageType3D::Pointer image){ _SomaImage = image; };
	void RunTracing();
	void RunGVFTracing(bool flagPreComputedGVFVessel);
	void WriteMultipleSWCFiles(std::string fname, unsigned int );	
	void WriteSWCFile(std::string , unsigned int );
	vtkSmartPointer< vtkTable > GetSWCTable(unsigned int);
	vtkSmartPointer< vtkTable > GetDebrisTable(unsigned int);
	void GenerateTestImage(); 

	int optionsCreate(const char* optfile, std::map<std::string,std::string>& options);


	//Additions after the pipeline
	void setLogScale( ImageType3D::Pointer, int scale );
	void setDiceSize( itk::Size<3> );
	void setDiceIndex( itk::Index<3> );
	void setFlagPipeline( bool flag ){_flagPipeline = flag;};
	void setFlagOutLog( bool flag ){_flagOutLog = flag;};
	void runNDX();
	ImageType3D::Pointer getNDX(){return _NDXImage;};
	CharImageType3D::Pointer getNDX2(){return _NDXImage2;};
	void setNDX( ImageType3D::Pointer image){_NDXImage = image;};

	//Parameters read from txt file
	float intensity_threshold;
	float contrast_threshold;
	int cost_threshold;
	float debris_threshold;
	int offshoot;
	int device;
	float mu;
	int tracing_type; // 1 for LOG; 2 for GVF - Default - GVF Tracing
	int noOfIteration;
	void RemoveSoma( LabelImageType3D::Pointer image2 );

	void OptimizeCoverage(std::string, bool);
	void ComputeObjectnessImage(ObjectnessMeasures_micro obj_measures);
	void Set_isCoverageOptimized(bool);
	void computeGVF(int mu, int ITER, int smoothing_scale);
	void computeGVF_2(ImageType3D::Pointer &image,int mu, int ITER, int smoothing_scale);
	void ComputeGVFVesselness();
	void ComputeGVFVesselness_2(ImageType3D::Pointer &image);

	GradientImageType::Pointer getGVFImage(){return _IGVF;};
	ImageType3D::Pointer getVessleness(){return _IVessel;};
	void setGVFImage(GradientImageType::Pointer &image){_IGVF = image;};
	void setVesselImage(ImageType3D::Pointer &image){ _IVessel = image;};
		
protected:
	void FeatureMain();
	void UpdateNDXImage_GVF(bool flagPreComputedGVFVessel);
	void ComputeMultiVesselness(double sigma_min, double sigma_max, int sigma_step);
	void SeedDetection(float th, int detection_method, int radius);
	void SeedAdjustment(int iter_num);
	bool SeedSparsify(SamplePointer seeds_candidate, SeedType query_point, int radius);
	void seed_centroid();
	void normalizeGVF();
	void threeLevelMinErrorThresh(unsigned char* im, float* Alpha1, float* Alpha2, float* Alpha3, float* P_I1, float* P_I2, size_t r, size_t c, size_t z);
	ProbImagePointer Upsampling(ProbImagePointer, int);
	LabelImagePointer UpsamplingLabel(LabelImagePointer,int);

	void GetFeature( float );
	
	bool IsPlate(const itk::FixedArray<float, 3> & , unsigned int & );
	bool IsPlate_control(const itk::FixedArray<float, 3> &, unsigned int & );
	bool IsDebris(const itk::FixedArray<float, 3> & , unsigned int &, float val );
	bool RegisterIndex(const float, itk::Index<3> &, itk::Size<3> &, long);
	bool RegisterIndexDebris(const float, itk::Index<3> &, itk::Size<3> &, long);
	SWCNode* TBack(itk::Index<3> & ndx, std::vector<IndexType> &  );
	float GetCost(SWCNode* , itk::Index<3> &  );
	float GetCostLocal(SWCNode* , itk::Index<3> & );
	float GetCostLocalLabel(SWCNode* , itk::Index<3> & );
	void ScanNeighbors( PixelType & a1,PixelType & a2,PixelType & a3, itk::Index<3> &);
	PixelType Update( PixelType a1,  PixelType a2,  PixelType a3,   PixelType P ) ;
	void Decimate();
	void Interpolate(PixelType);
	void RemoveIntraSomaNodes();	
	float getRadius(itk::Vector<float,3> & pos);
	void WriteImage3D(std::string , ImageType3D::Pointer );
	void BlackOut(itk::Index<3> &ndx );
	void GetFeature_2( float, int );
	//ProbImagePointer extract_one_component(int index, GradientImagePointer IG);

private:
	int SM;
	int SN;
	int SZ;
	vnl_vector<int> visit_label;

	std::vector<SWCNode*> _SWCNodeContainer;
	std::vector<DebrisNode*> _DebrisNodeContainer;
	//CharImageType3D::Pointer SomaImage;
	LabelImageType3D::Pointer _SomaImage;
	PixelType _CostThreshold;
	std::priority_queue < HeapNode* , std::vector<HeapNode*>,  Comparison > _PQ;
	ImageType3D::Pointer _PaddedCurvImage, _ConnImage, _NDXImage, _MaskedImage,_IVessel;   //Input Image, EK image, CT image
	GradientImagePointer _IGVF;
	PointList3D SeedPt;
	ImageType3D::Pointer ObjectnessImage;
	CharImageType3D::Pointer _NDXImage2;
	SWCImageType3D::Pointer _SWCImage; //swc label image
	itk::Size<3> _size;
	std::vector<OffsetType> _off;
	long _CurrentID;
	std::vector<IndexType> _StartPoints;
	unsigned int _padz;
	bool debug;
	//Additions after including in the pipeline
	itk::Size<3> _sizeDice;
	itk::Index<3> _indxDice;
	ImageType3D::Pointer _logScale_1; 
	ImageType3D::Pointer _logScale_2;
	ImageType3D::Pointer _logScale_3;
	ImageType3D::Pointer _logScale_4;
	ImageType3D::Pointer _logScale_5;
	ImageType3D::Pointer _logScale_6;
	
	bool _flagPipeline;
	bool _flagOutLog;
	bool _isCoverageOptimized;
	float _v_threshold; 

};

///////////////////////////////////////////////////////////////

class SWCNode 
{
public:
	long ID, PID, TreeID;
	itk::Index<3> ndx;
	itk::Vector<float,3> pos;
	bool IsLeaf, IsBranch, IsActive;
	SWCNode *parent;
	std::vector<SWCNode*> children;

	SWCNode(); 
	SWCNode(long id, long parent_id, long tree_id, itk::Index<3> index );
	SWCNode(long id, SWCNode *parent, long tree_id, itk::Index<3> index ); 
	
	//static bool IsIndexSame(itk::Index<3>, itk::Index<3>);

};

class DebrisNode 
{
public:
	itk::Index<3> ndx;

	DebrisNode(itk::Index<3> dndx); 

};


///////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////
class VectorPixelAccessor  
{
public:
  typedef itk::CovariantVector<float,3>   InternalType;
  typedef                      float      ExternalType;

  void operator=( const VectorPixelAccessor & vpa )
    {
      m_Index = vpa.m_Index;
    }
  ExternalType Get( const InternalType & input ) const 
    {
    return static_cast<ExternalType>( input[ m_Index ] );
    }
  void SetIndex( unsigned int index )
    {
    m_Index = index;
    }
private:
  unsigned int m_Index;
};

#endif

