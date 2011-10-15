#ifndef _MICROGLIA_PROCESS_TRACER_H_
#define _MICROGLIA_PROCESS_TRACER_H_

#include "itkImage.h"
#include "itkArray.h"
#include "itkImageFileReader.h"
//debug purposes only
#include "itkImageFileWriter.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

#include "itkBresenhamLine.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkRobustAutomaticThresholdImageFilter.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "itkSymmetricSecondRankTensor.h"
#include "itkMaskNegatedImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "vnl/vnl_math.h"

#include <fstream>
#include <list>
#include <limits>
#include <map>
#include <string>
#include <vector>
#include <time.h>

typedef float PixelType;
typedef itk::Index<3> IndexType;

struct Node
{
  Node();
  
  long ID;
  IndexType index;
  bool IsRoot;
  bool IsOpen;
  Node *parent;
  std::vector< Node * > children; //needed for multiple .swc file write mode.
};

//this is only safe for images with extents that can be measured with 5 digits in each dimension...
struct CompareIndices
{
  bool operator() (const itk::Index<3>& lhs, const itk::Index<3>& rhs) const
    {
    unsigned long long lval = lhs[0] + (100000 * lhs[1]) + (10000000000 * lhs[2]); 
    unsigned long long rval = rhs[0] + (100000 * rhs[1]) + (10000000000 * rhs[2]); 
    return (lval < rval);
    }
};

class MicrogliaProcessTracer
{
public:

  typedef itk::Image< PixelType, 3 >  ImageType3D;
  typedef itk::Image< unsigned char, 3 > CharImageType3D;

  typedef itk::ImageFileReader<ImageType3D> ReaderType;
  typedef itk::RescaleIntensityImageFilter<ImageType3D, ImageType3D> RescalerType;
  typedef itk::MaskNegatedImageFilter<ImageType3D, CharImageType3D, ImageType3D> MaskFilterType;

  MicrogliaProcessTracer();
  ~MicrogliaProcessTracer();

  void LoadInputImage(std::string fname);
  void LoadSeedPoints(std::string fname);
  void LoadNodes(std::string fname);
  void LoadSomaImage(std::string somaFileName);
  void SetPadding(unsigned int i) { this->Padding = i; }
  void SetMaxDistance(double d) { this->MaxDistance = d; }

  void RunTracing();
  void WriteSWC( std::string fname );
    
protected:
  void LoadInputImage(ImageType3D::Pointer &image);
  void CalculateCriticalPoints();
  void CalculateCriticalPointsAtScale( float );
  void ComputeAdjacencies( std::vector< Node * > );
  unsigned int ShapeAnalysis(const itk::FixedArray<float, 3> & );
  bool RegisterIndex(const float, itk::Index<3> &, itk::Size<3> &, long);
  //float GetRadius(itk::Vector<float,3> & pos);
  float GetRadius(itk::Index<3> pos);
  void BuildTrees();
  std::vector< Node * > ReadListOfPoints(std::string fname);
  std::pair< Node *, Node * > FindClosestOpenNode();
  void MaskAwaySomas();
  double GetDistanceBetweenPoints(itk::Index<3> start, itk::Index<3> end);

private:
  CharImageType3D::Pointer SomaImage;
  ImageType3D::Pointer InputImage, PaddedInputImage, NDXImage, ThresholdedImage;
  //std::vector<IndexType> StartPoints;
  itk::BresenhamLine<3> Line;

  unsigned int Padding;
  std::vector< Node * > Open;
  std::vector< Node * > Closed;
  std::map< itk::Index<3>, Node *, CompareIndices > IndexToNodeMap;
  long NodeCounter;
  double MaxDistance;

  std::map< Node *, std::list< std::pair< double, Node *> > > AdjacencyMap;
};

#endif

