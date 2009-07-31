#ifndef __HELPERS_H
#define __HELPERS_H
#include<stdio.h>
#define MPICH_IGNORE_CXX_SEEK
//#include <mpi.h>
#include <stdio.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageSliceConstIteratorWithIndex.h>
#include <itkImageLinearIteratorWithIndex.h>
#include <itkConstNeighborhoodIterator.h>
#include <itkMedianImageFilter.h>
#include <vector>
#include <algorithm>
#include <itkBinaryThresholdImageFilter.h>
#include <itkBinaryMedianImageFilter.h>
#include <itkOtsuThresholdImageFilter.h>
#include <itkShiftScaleImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkScalarConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkBinaryErodeImageFilter.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkDanielssonDistanceMapImageFilter.h>
#include <itkVector.h>
#include <itkVTKImageExport.h>
#include <itkVTKImageImport.h>
#include "itkLabelGeometryImageFilter.h"
#include <iostream>
#include <sstream>
#include <string>
#include <iterator> 
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_hungarian_algorithm.h>

#include "itkLabelStatisticsImageFilter.h"
#include "itkScalarImageTextureCalculator.h"
#include "ftkLabelImageToFeatures.h"
#include "ftkIntrinsicFeatures.h"
#include "ftkTrackFeatures.h"
#include "itkBinaryThinningImageFilter3D.h"

//int npes,rank;
#define PRINTF(...) PRINTF1(__VA_ARGS__,1)
#define PRINTF1(s,...) {printf("%d/%d: "s,rank,npes,__VA_ARGS__);}//{printf("%d/%d: ",rank,npes);printf(__VA_ARGS__);}

#define STRIDE 1


// tasks

#define DO_UNMIX 0
#define DO_NOTHING 1

// MPI tags

#define IMAGE_UNMIX_INPUT 0
#define IMAGE_UNMIX_OUTPUT 1
#define MAX_CHANNELS 4

#define EPSILON 1e-6


#define MAX(a,b) (((a) > (b))?(a):(b))
#define MIN(a,b) (((a) < (b))?(a):(b))

//typedefs
namespace helpers{

typedef unsigned char InputPixelType;
typedef unsigned char OutputPixelType;

typedef itk::Vector<unsigned char, 3> VectorPixelType;
typedef itk::Image<VectorPixelType, 3> ColorImageType;
typedef itk::Image<VectorPixelType, 2> Color2DImageType;

typedef itk::Image<InputPixelType,3> InputImageType;
typedef itk::Image<OutputPixelType,3> OutputImageType;
typedef itk::Image<short int,3> LabelImageType;
typedef itk::Image<short int,2> Label2DImageType;
typedef itk::Image<float,3> FloatImageType;


typedef itk::Image<InputPixelType,2> Input2DImageType;
typedef itk::Image<InputPixelType,2> Output2DImageType;

typedef itk::ImageRegionConstIterator<InputImageType> ConstIteratorType;
typedef itk::ImageRegionIterator<InputImageType> IteratorType;

typedef itk::ImageRegionConstIterator<LabelImageType> ConstLabelIteratorType;
typedef itk::ImageRegionIterator<LabelImageType> LabelIteratorType;

typedef itk::ImageRegionConstIterator<Input2DImageType> Const2DIteratorType;
typedef itk::ImageRegionIterator<Input2DImageType> twoDIteratorType;

typedef itk::ImageRegionConstIterator<ColorImageType> ConstColorIteratorType;
typedef itk::ImageRegionIterator<ColorImageType> ColorIteratorType;

typedef itk::ImageLinearIteratorWithIndex< Input2DImageType > LinearIteratorType;
typedef itk::ImageSliceConstIteratorWithIndex< InputImageType > SliceIteratorType;

typedef itk::ImageLinearIteratorWithIndex< Color2DImageType > LinearColorIteratorType;
typedef itk::ImageSliceConstIteratorWithIndex< ColorImageType > SliceColorIteratorType;

typedef itk::ConstNeighborhoodIterator<InputImageType> NeighborhoodIteratorType;

typedef itk::MedianImageFilter<InputImageType,InputImageType> MedianFilterType;
typedef itk::BinaryMedianImageFilter<InputImageType,InputImageType> BinaryMedianFilterType;

typedef itk::Image<bool,3> BoolImageType;
typedef itk::BinaryThresholdImageFilter<InputImageType,OutputImageType> ThresholdFilterType;
typedef itk::OtsuThresholdImageFilter<InputImageType,OutputImageType> OtsuThresholdFilterType;

typedef itk::Image<short int,3> DistanceImageType;
typedef itk::DanielssonDistanceMapImageFilter<InputImageType,DistanceImageType> DistanceMapFilterType;
typedef DistanceMapFilterType::VectorImageType OffsetImageType;

typedef itk::ConnectedComponentImageFilter<InputImageType,LabelImageType> ConnectedFilterType;
typedef itk::RelabelComponentImageFilter<LabelImageType,LabelImageType> RelabelFilterType;

typedef ftk::IntrinsicFeatures FeaturesType;

typedef itk::VTKImageExport<InputImageType> ExportFilterType;

struct cubecoord{
	unsigned short sx,sy,sz;
	unsigned short ex,ey,ez;
};


struct PQdata{
	InputImageType::IndexType index;
	int depth;
};

struct Vec3f{
	double x,y,z;
	void Normalize()
	{
		double sum = sqrt(x*x+y*y+z*z);
		if(sum!=0)
		{
		x/=sum;
		y/=sum;
		z/=sum;
		}
	}
};


// all header declarations
//
//
double features_box_overlap(FeaturesType &f1, FeaturesType &f2);
double features_diff(FeaturesType &f1, FeaturesType &f2,bool overlap);

Input2DImageType::Pointer getProjection(InputImageType::Pointer im);
Color2DImageType::Pointer getColorProjection(ColorImageType::Pointer im);
Input2DImageType::Pointer getCollage(InputImageType::Pointer im[4]);
unsigned char getMedianValue(NeighborhoodIteratorType it,int size);
void unmix_median(InputImageType::Pointer im[4],InputImageType::Pointer om[4],InputImageType::Pointer assignment[4]);
void unmix_neighborhood(InputImageType::Pointer im[4],InputImageType::Pointer om[4]);
OutputImageType::Pointer getThresholded(InputImageType::Pointer im,int n);
OutputImageType::Pointer getOtsuThresholded(InputImageType::Pointer im);
OutputImageType::Pointer getBinaryMedianFiltered(InputImageType::Pointer im, InputImageType::SizeType radius);
OutputImageType::Pointer getScaledFromBool(BoolImageType::Pointer im);
InputImageType::Pointer getLargeComponents(InputImageType::Pointer im, int n);
LabelImageType::Pointer getLargeLabels(LabelImageType::Pointer im, int n);
LabelImageType::Pointer getFeatureVectors(LabelImageType::Pointer im, InputImageType::Pointer in_image,std::vector<FeaturesType> &feature_vector,int time,int tag);
void getFeatureVectorsOld(InputImageType::Pointer im,std::vector<FeaturesType> &feature_vector,int time,int tag);
InputImageType::Pointer getDilated(InputImageType::Pointer im, int n);
InputImageType::Pointer getEroded(InputImageType::Pointer im, int n);
DistanceMapFilterType::Pointer getDistanceMap(InputImageType::Pointer im);
InputImageType::Pointer getOldSegmented(InputImageType::Pointer im_input,int threshold, int min_component_size, int morph_opening_depth);
LabelImageType::Pointer getLabelled(InputImageType::Pointer im_input,int threshold, int min_component_size, int morph_opening_depth);
ColorImageType::Pointer getColorCompositeImage(InputImageType::Pointer im[4],VectorPixelType colors[4]);
ColorImageType::Pointer getColorComposite(InputImageType::Pointer im[],int n, VectorPixelType colors[]);
InputImageType::Pointer getImageFromNPTS(char *filename_npts,int imagesize[]);
void getClassified(DistanceImageType::Pointer dist, InputImageType::Pointer micro, InputImageType::Pointer &p1, InputImageType::Pointer &p2);
void writeNPTS(InputImageType::Pointer im, char*filename);
//void MSA_classifymicroglia(char*filename_vesselnpts,char*filename_microglianpts);
InputImageType::Pointer getPreprocessed(InputImageType::Pointer im);
ColorImageType::Pointer getColorCompositeImageFromLabelled(LabelImageType::Pointer im,VectorPixelType color);
Color2DImageType::Pointer getColor2DImage(LabelImageType::Pointer labelled,int channel);
Color2DImageType::Pointer getColorFromGrayScale(Input2DImageType::Pointer);
Input2DImageType::Pointer get2DBoundary(LabelImageType::Pointer label);
InputImageType::Pointer getEmpty(int,int,int);
Input2DImageType::Pointer get2DEmpty(int, int);
Color2DImageType::Pointer getColorBoundaryImage(LabelImageType::Pointer, InputImageType::Pointer,int);
LabelImageType::Pointer getLabelsMapped(LabelImageType::Pointer label, std::vector<FeaturesType> &fvec, unsigned int * indices);
void getFeatureVectorsFarsight(LabelImageType::Pointer im, InputImageType::Pointer in_image, std::vector<FeaturesType> & feature_vector, int time, int tag);
InputImageType::Pointer getLabelToBinary(LabelImageType::Pointer l);
InputImageType::Pointer getMaxImage(InputImageType::Pointer,InputImageType::Pointer);
std::vector<float> traverseCenterline(itk::ImageRegionIteratorWithIndex<InputImageType> iter,InputImageType::Pointer im,char neighbors[26][3],int n);
void AnalyzeTimeFeatures(std::vector<ftk::TrackFeatures> &tfs);
void AnalyzeVesselCenterlines(InputImageType::Pointer cline, std::vector<ftk::TrackFeatures> &tfs);
FloatImageType::IndexType searchNearestVesselDirection(FloatImageType::Pointer dir_image[3],FloatImageType::IndexType index,InputImageType::Pointer vesselim);
void AnalyzeDCContact(LabelImageType::Pointer segmented[][4], std::vector<ftk::TrackFeatures> &tfs, int c, int num_t);
}
#endif
