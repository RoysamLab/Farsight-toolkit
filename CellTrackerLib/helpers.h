#ifndef __HELPERS_H
#define __HELPERS_H
#include<stdio.h>
#include <iostream>
#define MPICH_IGNORE_CXX_SEEK
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageSliceConstIteratorWithIndex.h>
#include <itkImageLinearIteratorWithIndex.h>
#include <itkConstNeighborhoodIterator.h>
#include <itkLineIterator.h>
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
#include <itkGrayscaleDilateImageFilter.h>
#include <itkGrayscaleErodeImageFilter.h>
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
#include "itkScalarImageToTextureFeaturesFilter.h"
#include "ftkTrackFeatures.h"
#include "itkBinaryThinningImageFilter3D.h"

#ifdef _OPENMP
#include "omp.h"
#endif


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



#define MAX(a,b) (((a) > (b))?(a):(b))
#define MIN(a,b) (((a) < (b))?(a):(b))
#define SIGN(a)  (((a)>=0)?1:-1)

//typedefs
namespace helpers{

  typedef unsigned char InputPixelType;
  typedef unsigned char OutputPixelType;
  //typedef short int LabelPixelType;
  typedef unsigned short LabelPixelType;

  typedef itk::Vector<unsigned char, 3> VectorPixelType;
  typedef itk::Image<VectorPixelType, 3> ColorImageType;
  typedef itk::Image<VectorPixelType, 2> Color2DImageType;

  typedef itk::Image<InputPixelType,3> InputImageType;
  typedef itk::Image<OutputPixelType,3> OutputImageType;
  typedef itk::Image<LabelPixelType,3> LabelImageType;
  typedef itk::Image<LabelPixelType,2> Label2DImageType;
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

  typedef ftk::IntrinsicFeatures FeatureType;

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

  struct FeatureVariances{
    enum ScalarTypes
    {
      VOLUME, INTEGRATED_INTENSITY, ECCENTRICITY, ELONGATION, ORIENTATION, BBOX_VOLUME, \
        SUM, MEAN, MEDIAN, MINIMUM, MAXIMUM, SIGMA, VARIANCE, \
        SURFACE_GRADIENT, INTERIOR_GRADIENT, SURFACE_INTENSITY, INTERIOR_INTENSITY, \
        INTENSITY_RATIO, CONVEXITY,RADIUS_VARIATION, SURFACE_AREA, SHAPE, SHARED_BOUNDARY, \
        SKEW, ENERGY, ENTROPY,
      T_ENERGY, T_ENTROPY, INVERSE_DIFFERENCE_MOMENT, INERTIA, CLUSTER_SHADE, CLUSTER_PROMINENCE,

    };	//FEATURES WILL GET ASSIGNED INT 0,1,...N-1

    static const int N = CLUSTER_PROMINENCE + 1;
    float variances[N];
    float means[N];
    int Dimensions;

    float BoundingBox[6]; // start_x end_x start_y end_y start_z end_z
    float distVariance;
    float distMean;
    float boundDistVariance;
    float boundDistMean;
    float spacing[3];
    float timeVariance;
    float timeMean;
    int time_last;
    float overlapVariance;
    float overlapMean;

    float MS_prior;
    float AD_prior;
    float T_prior;

    bool FlagRoi; // ?? FIXME what is this ?
  };
  // all header declarations
  //
  //

  double features_box_overlap(FeatureType &f1, FeatureType &f2);
  double features_diff(FeatureType &f1, FeatureType &f2,bool overlap);

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
  LabelImageType::Pointer getFeatureVectors(LabelImageType::Pointer im, InputImageType::Pointer in_image,std::vector<FeatureType> &feature_vector,int time,int tag);
  void getFeatureVectorsOld(InputImageType::Pointer im,std::vector<FeatureType> &feature_vector,int time,int tag);
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
  LabelImageType::Pointer getLabelsMapped(LabelImageType::Pointer label, std::vector<FeatureType> &fvec, unsigned int * indices);
  void getFeatureVectorsFarsight(LabelImageType::Pointer im, InputImageType::Pointer in_image, std::vector<FeatureType> & feature_vector, int time, int tag);
  InputImageType::Pointer getLabelToBinary(LabelImageType::Pointer l);
  InputImageType::Pointer getMaxImage(InputImageType::Pointer,InputImageType::Pointer);
  std::vector<float> traverseCenterline(itk::ImageRegionIteratorWithIndex<InputImageType> iter,InputImageType::Pointer im,char neighbors[26][3],int n);
  void AnalyzeTimeFeatures(std::vector<ftk::TrackFeatures> &tfs, float spacing[3]);
  void AnalyzeTimeFeatures(std::vector<ftk::TrackFeatures> &tfs);
  void PrintTrackFeatures(std::vector<ftk::TrackFeatures> &tfs, std::string path);

  void AnalyzeVesselCenterlines(InputImageType::Pointer cline, std::vector<ftk::TrackFeatures> &tfs, float spacing[3]);
  FloatImageType::IndexType searchNearestVesselDirection(FloatImageType::Pointer dir_image[3],FloatImageType::IndexType index,InputImageType::Pointer vesselim);
  void AnalyzeDCContact(LabelImageType::Pointer segmented[][4], std::vector<ftk::TrackFeatures> &tfs, int c, int num_t, float spacing[3]);
  LabelImageType::Pointer extract_label_image(int label, float bbox[6],LabelImageType::Pointer l);
  InputImageType::Pointer extract_raw_image(float bbox[6],InputImageType::Pointer r);
  void annotateImage(Color2DImageType::Pointer number,Color2DImageType::Pointer orig, int n, int x, int y);
  ColorImageType::Pointer getColorImageFromColor2DImages(std::vector<Color2DImageType::Pointer> input);
  void drawLine(ColorImageType::Pointer input, VectorPixelType color1, VectorPixelType color2, int x1, int y1, int z1, int x2, int y2, int z2);
  std::vector<FeatureType> get_all_connected_components(LabelImageType::Pointer,FeatureType);
  void SplitCell(LabelImageType::Pointer lin, InputImageType::Pointer imin,FeatureType fin, FeatureVariances fvar,std::vector<LabelImageType::Pointer> &lout,std::vector<InputImageType::Pointer> &rout,std::vector<FeatureType> &fvecout);
  void MergeCells(std::vector<LabelImageType::Pointer> lin, std::vector<InputImageType::Pointer> imin, std::vector<FeatureType> fin, FeatureVariances fvar, LabelImageType::Pointer &lout, InputImageType::Pointer &rout, FeatureType &fout);
  LabelImageType::Pointer fillHoles(LabelImageType::Pointer im, int n);

  int relabelWells(std::vector<LabelImageType::Pointer> & tracked_images, int maxPreviousLabel );

}
#endif
