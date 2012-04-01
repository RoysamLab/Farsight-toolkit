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

#include "omp.h"

typedef float PixelType;

class SWCNode;
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
	typedef itk::Image< SWCNode*, 3 > SWCImageType3D;
	typedef itk::Image< unsigned int, 3 > LabelImageType3D;
	typedef LabelImageType3D::PixelType * LabelArrayType;

	//Constructor
	MultipleNeuronTracer();
	//Destructor
	~MultipleNeuronTracer();


	void LoadCurvImage(std::string fname, unsigned int pad); 
	void LoadCurvImage_1(ImageType3D::Pointer &image, unsigned int pad);
	void LoadCurvImage_2(ImageType3D::Pointer &image);
	
	void ReadStartPoints(std::string fname, unsigned int padz);
	void ReadStartPoints_1(std::vector< itk::Index<3> > somaCentroids, unsigned int padz);
	void SetCostThreshold(float thres){_CostThreshold = thres;};
	void LoadSomaImage(std::string somaFileName);
	void LoadSomaImage_1(LabelImageType3D::Pointer image){ _SomaImage = image; };
	void RunTracing();
	void WriteMultipleSWCFiles(std::string fname, unsigned int );	
	void WriteSWCFile(std::string , unsigned int );
	vtkSmartPointer< vtkTable > GetSWCTable(unsigned int);
	void GenerateTestImage(); 
	
	//Additions after the pipeline
	void setLogScale( ImageType3D::Pointer, int scale );
	void setDiceSize( itk::Size<3> );
	void setDiceIndex( itk::Index<3> );
	void setFlagPipeline( bool flag ){_flagPipeline = flag;};
	void setFlagOutLog( bool flag ){_flagOutLog = flag;};
	void runNDX();
	ImageType3D::Pointer getNDX(){return _NDXImage;};
	void setNDX( ImageType3D::Pointer image){_NDXImage = image;};
		
protected:
	void FeatureMain();
	void GetFeature( float );
	
	bool IsPlate(const itk::FixedArray<float, 3> & , unsigned int & );
	bool RegisterIndex(const float, itk::Index<3> &, itk::Size<3> &, long);
	SWCNode* TBack(itk::Index<3> & ndx, std::vector<IndexType> &  );
	float GetCost(SWCNode* , itk::Index<3> &  );
	float GetCostLocal(SWCNode* , itk::Index<3> & );
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
	//CharImageType3D::Pointer SomaImage;
	LabelImageType3D::Pointer _SomaImage;
	PixelType _CostThreshold;
	std::priority_queue < HeapNode* , std::vector<HeapNode*>,  Comparison > _PQ;
	ImageType3D::Pointer _PaddedCurvImage, _ConnImage, _NDXImage;   //Input Image, EK image, CT image
	SWCImageType3D::Pointer _SWCImage; //swc label image
	itk::Size<3> _size;
	std::vector<OffsetType> _off;
	long _CurrentID;
	std::vector<IndexType> _StartPoints;
	unsigned int _padz;

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

///////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////

#endif

