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
#include "itkImage.h"
#include "itkArray.h"
#include "itkCovariantVector.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkStatisticsImageFilter.h"
#include "itkHuangThresholdImageFilter.h"
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
#include "vnl/vnl_math.h"
#include "vtkSmartPointer.h"
#include "vtkTable.h"
#include "vtkVariantArray.h"
#include "vtkDoubleArray.h"


#include "itkRegionOfInterestImageFilter.h"
#include "itkImageDuplicator.h"

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

class SWCNode;
class DebrisNode;
class HeapNode;
class Comparison;
//class MultipleNeuronTracer;

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
	typedef LabelImageType3D::PixelType * LabelArrayType;
	typedef itk::HuangThresholdImageFilter<ImageType3D,ImageType3D> HuangThresholdFilterType;
	typedef itk::OtsuThresholdImageFilter<ImageType3D,ImageType3D>  OtsuThresholdImageFilterType;
	typedef itk::BinaryThresholdImageFilter <ImageType3D,ImageType3D>	BinaryThresholdImageFilterType;
	typedef itk::MinimumMaximumImageCalculator<ImageType3D> MinMaxImageCalculatorType;

	//Constructor
	MultipleNeuronTracer();
	//Destructor
	~MultipleNeuronTracer();

	void LoadParameters(const char* parametersFileName, int _argc);
	void LoadParameters_1(const char* parametersFileName,float intensityThreshold,float contrastThreshold,int costThreshold);

	void LoadCurvImage(std::string fname, unsigned int pad); 
	void LoadCurvImage_1(ImageType3D::Pointer &image, unsigned int pad);
	void LoadCurvImage_2(ImageType3D::Pointer &image);
	
	void ReadStartPoints(std::string fname, unsigned int padz);
	void ReadStartPoints_1(std::vector< itk::Index<3> > somaCentroids, unsigned int padz);
	void ReadStartPoints_2(std::string fname, unsigned int pad,float startx,float starty,float startz, float widthx,float widthy,float widthz); //
	void SetCostThreshold(float thres){_CostThreshold = thres;};
	void LoadSomaImage(std::string somaFileName);
	void LoadSomaImage_1(LabelImageType3D::Pointer image){ _SomaImage = image; };
	void RunTracing();
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
	void RemoveSoma( LabelImageType3D::Pointer image2 );

		
protected:
	void FeatureMain();
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

private:
	std::vector<SWCNode*> _SWCNodeContainer;
	std::vector<DebrisNode*> _DebrisNodeContainer;
	//CharImageType3D::Pointer SomaImage;
	LabelImageType3D::Pointer _SomaImage;
	PixelType _CostThreshold;
	std::priority_queue < HeapNode* , std::vector<HeapNode*>,  Comparison > _PQ;
	ImageType3D::Pointer _PaddedCurvImage, _ConnImage, _NDXImage, _MaskedImage;   //Input Image, EK image, CT image
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

#endif

