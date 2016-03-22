#ifndef VOTINGSEG_H
#define VOTINGSEG_H

#include <time.h>
#include<iostream>
#include<vector>
#include <cmath>

#include "itkImage.h"
#include <itkIndex.h>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkLabelGeometryImageFilter.h"
#include "itkConfidenceConnectedImageFilter.h"
#include "itkCurvatureFlowImageFilter.h"
#include "itkLaplacianSharpeningImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkGradientImageFilter.h"
#include "itkGradientVectorFlowImageFilter.h"
#include "itkGeodesicActiveContourLevelSetImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkFastMarchingImageFilter.h"
#include "itkCurvatureAnisotropicDiffusionImageFilter.h"
#include "itkSigmoidImageFilter.h"

#include "ftkCommon/ftkUtils.h"
#include "ftkAdapSeedDetection/ftkVotingGlobal.h"
#include "ftkAdapSeedDetection/ftkVoting.h"
#include "ftkAdapSeedDetection/ftkVoting_3D.h"

#ifdef _OPENMP
#include "omp.h"
#endif

#include "GeneralIO.h"
//#include "seed_2D.h"

static const int Dimension_2D = 2;
typedef unsigned char UCPixelType;
typedef unsigned UIPixelType;
typedef double DPixelType;
typedef unsigned short USPixelType;
typedef itk::Image<UCPixelType, Dimension_2D> UCImageType_2D;
typedef itk::Image<UIPixelType, Dimension_2D> UIImageType_2D;
typedef itk::Image<DPixelType, Dimension_2D> DImageType_2D;
typedef itk::Image<USPixelType, Dimension_2D> USImageType_2D;
typedef itk::ImageFileWriter < USImageType_2D > USWriterType;
typedef itk::ImageRegionIteratorWithIndex< UCImageType_2D > UCIteratorType_2D;
typedef itk::ImageRegionIteratorWithIndex< DImageType_2D > DIteratorType_2D;

struct Paras 
{
	int voting_rmin;
	int voting_rmax;
	int voting_radius;
	double voting_gradthresh;
	double voting_scale;
	unsigned int smooth_Iteration;
	double smooth_conductance;
	double gradient_sigma;
	double sigmoid_alpha;
	double sigmoid_beta;
	int fastmarching_timethreshold;
	double ac_curvature;
	double ac_advect;
	int ac_iteration;
	double ac_propagation;
	double ac_maximumRMSE;
};

typedef struct { std::string name; double value; } Parameter;

class Voting_Seg
{

public:

	Voting_Seg();
	~Voting_Seg();
	void Initialize();
	void SetInputImage(DImageType_2D::Pointer im);
	void SetParameters(std::vector<Parameter>);
	void Update();
	void GetOutputImage();
	void WriteImage(DImageType_2D::Pointer image, const char* filename);

	DImageType_2D::Pointer InputImage;
	DImageType_2D::Pointer OutputImage;
	Paras parameters;

protected:

	void Binarize_2D();
	void SeedDetection_2D();
	void Watershed_2D();
	void LevelSet_2D();

private:

	DImageType_2D::Pointer VotingImage;
	UCImageType_2D::Pointer SeedImage;
	UCImageType_2D::Pointer SeedLabelImage;
	DImageType_2D::Pointer SegmentedImage;

	std::vector< itk::Index<Dimension_2D> > seedVector;

	DImageType_2D::SizeType imagesize;

	int seedNum;

	void DefultPara();
	void SeedExtraction();
};

#endif
