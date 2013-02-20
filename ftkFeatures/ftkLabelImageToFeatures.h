/*=========================================================================
Copyright 2009 Rensselaer Polytechnic Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 
=========================================================================*/

/*=========================================================================

  Program:   Farsight Biological Image Segmentation and Visualization Toolkit
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision: 0.00 $

=========================================================================*/
//
// FOR CLASS DETAILS SEE ftkLabelImageToFeatures.txx
//
#ifndef __ftkLabelImageToFeatures_h
#define __ftkLabelImageToFeatures_h

#include <itkLightObject.h>
#include <itkObjectFactory.h>

#include <itkLabelGeometryImageFilter.h>
#include <itkLabelStatisticsImageFilter.h>
#include <itkScalarImageToTextureFeaturesFilter.h>
#include <itkGradientMagnitudeImageFilter.h>
#include <itkSignedDanielssonDistanceMapImageFilter.h>
#include <itkLineIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>

#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkVariantArray.h>
#include <vtkTable.h>


#include "ftkFeatures/Qhull/libqhull.h"
#include "ftkFeatures/Qhull/user.h"
#include "ftkFeatures/Qhull/stat.h"
#include "ftkFeatures/Qhull/io.h"
#include "ftkFeatures/Qhull/merge.h"
#include "ftkFeatures/Qhull/geom.h"
#include "ftkFeatures/Qhull/poly.h"
#include "ftkFeatures/Qhull/qset.h"
#include "ftkFeatures/Qhull/mem.h"
#include "ftkFeatures/Qhull/random.h"
#include "ftkFeatures/Qhull/qhull_a.h"
#include "ftkFeatures/Qhull/qhull.h"

#include "ftkIntrinsicFeatures.h"
#include "ftkObject.h"
#include "ftkImage/ftkImage.h"

#ifdef ZERNIKE
#include "Zernike/zernike.h"
#endif

#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <set>

namespace ftk
{

class IntrinsicFeatureCalculator
{
public:
	typedef unsigned char IPixelT;
	typedef unsigned short LPixelT;

	IntrinsicFeatureCalculator();
	bool SetInputImages(ftk::Image::Pointer intImg, ftk::Image::Pointer labImg, int intChannel=0, int labChannel=0, bool CytoImage = false);
	void SetFeaturesOn(void);							//Turn on all features
	void SetFeaturesOn(std::set<int> onFeats);			//Turn on only these features
	void SetFeatureOn(int feat, bool v = true);			//Set the value of this feature
	void SetFeaturePrefix(std::string prefix);			//Set Prefix for feature names
	void SetRegion(int x1, int y1, int z1, int x2, int y2, int z2);	//Compute features for objects in this region
	void SetIDs(std::set<LPixelT> ids);					//Only update these ids
	void ClearRegion(void){ useRegion = false; };		//Clear the Region;
	void ClearIDs(void){ useIDs = false; IDs.clear(); };//Clear the IDs;

	vtkSmartPointer<vtkTable> Compute(void);			//Compute features that are ON and return table with values (for all objects)
	//void Update(vtkSmartPointer<vtkTable> table);		//Update the features in this table whose names match (sets doFeat)
//	void Update(vtkSmartPointer<vtkTable> table, std::map<int, ftk::Object::Point> * cc = NULL, std::map<int, ftk::Object::Box> * bbox = NULL, vtkSmartPointer<vtkTable> NucAdjTable = NULL);
	void GetObjectCentroids(vtkSmartPointer<vtkTable> table, int time=0);
	void Update(vtkSmartPointer<vtkTable> table, std::map<int, ftk::Object::Point> * cc = NULL, std::map<int, ftk::Object::Box> * bbox = NULL, vtkSmartPointer<vtkTable> NucAdjTable = NULL, int currtime = 0);
	void Append(vtkSmartPointer<vtkTable> table);		//Compute features that are ON and append them to the existing table (makes more columns)
	void GetDistanceToSurfaceMeasures(vtkSmartPointer<vtkTable> table, std::vector< ftk::Object::Point > surfacePoints);

private:
	ftk::Image::Pointer intensityImage;
	int intensityChannel;
	ftk::Image::Pointer labelImage;
	int labelChannel;
	std::string fPrefix;
	bool doFeat[IntrinsicFeatures::N];

	bool cyto_image;
	bool useRegion;
	LPixelT regionIndex[3];
	LPixelT regionSize[3];
	bool useIDs;
	std::set<LPixelT> IDs;


	int getMaxFeatureTurnedOn(void);
	bool needTextures(void);
	bool needHistogram(void);
	int needLevel(void);
};


//********************************************************************************************************
//THIS IS THE CLASS THAT DOES THE CALCULATIONS (THE ENGINE)
//********************************************************************************************************
template< typename TIPixel = unsigned char, typename TLPixel = unsigned short, unsigned int VImageDimension = 2> 
class LabelImageToFeatures : public itk::LightObject
{
public:

	typedef LabelImageToFeatures Self;
	typedef itk::LightObject Superclass;
	typedef itk::SmartPointer< Self > Pointer;
	typedef itk::SmartPointer< const Self > ConstPointer;

	typedef TIPixel IntensityPixelType;
	typedef TLPixel LabelPixelType;
	typedef itk::Image< IntensityPixelType, VImageDimension > IntensityImageType;
	typedef itk::Image< LabelPixelType, VImageDimension > LabelImageType;
	typedef itk::ImageRegionIterator< IntensityImageType> intIteratorType;
	typedef itk::ImageRegionIterator< LabelImageType> labIteratorType;
	typedef itk::Image< TIPixel,2 > Intensity2DType;
	typedef itk::Image< TLPixel,2 > Label2DType;
	typedef typename IntensityImageType::Pointer IntensityImagePointer;
	typedef typename LabelImageType::Pointer LabelImagePointer;
	typedef itk::LabelGeometryImageFilter< LabelImageType, IntensityImageType > LabelGeometryType;
	typedef typename LabelGeometryType::Pointer LabelGeometryPointer;

	//typedef zernike::ImageType zernikeImageType;
	typedef itk::Image< unsigned char,2> zernikeImageType;
	typedef itk::ImageRegionIterator< zernikeImageType >  zerIteratorType;

	itkNewMacro( Self );
	
	itkTypeMacro(LabelImageToFeatures, LightObject);

	bool SetImageInputs( IntensityImagePointer intImgIn, LabelImagePointer lblImgIn, bool CytoImage = false );
	bool SetImageInputs( IntensityImagePointer intImgIn, LabelImagePointer lblImgIn, TLPixel index[VImageDimension], TLPixel size[VImageDimension], bool CytoImage = false );
	bool SetCompleteImageInputs( IntensityImagePointer intImgIn, LabelImagePointer lblImgIn, bool CytoImage = false );

	void Update();
	LabelPixelType GetMaxLabel();
	float GetPercentSharedBoundary(TLPixel focusLabel, TLPixel neighborLabel);
	std::vector<TLPixel> GetContactNeighbors(TLPixel label);
	IntrinsicFeatures * GetFeatures( LabelPixelType label );
	std::vector< std::vector<double> > GetZernikeMoments( LabelPixelType label ){ return IDtoZernikeMap[label]; };
	std::vector< LabelPixelType > GetLabels() { return this->labels; };

	void ComputeSurfaceOn();
	void ComputeSurfaceOff(){ this->surfacecomputation = false; };
	void ComputeHistogramOn();
	void ComputeHistogramOff(){ this->computeHistogram = false; };
	void ComputeTexturesOn();
	void ComputeTexturesOff(){ this->computeTextures = false; };
	void SetLevel(short int newLevel);
	void SetZernikeOrder(int order){ this->zernikeOrder = order; };
	short int GetLevel(){ return computationLevel; };
	void GetAdjacency();

protected:
	LabelImageToFeatures();
	~LabelImageToFeatures(){};

private:
	LabelImageToFeatures(const Self&);  //purposely not implemented
	void operator=(const Self&);		//purposely not implemented

	//Internal Functions:
	bool RunLabelGeometryFilter();
	bool RunLabelStatisticsFilter();
	bool RunTextureFilter();
	bool RunTextureFilter2D();
	void RunZernikeFilter();
	void LabelImageScan();
	void CalculateScanFeatures();
	void SetHistogramParameters(int* numBins, int* lowerBound, int* upperBound);
	void RunSurfaceFeature(LabelGeometryPointer);

	//Internal Variables:
	IntensityImagePointer intensityImage;	//Input intensity image;
	LabelImagePointer labelImage;			//Input label image;
	int zernikeOrder;
	bool cyto_image;

	std::vector< std::vector< typename LabelImageType::IndexType > > boundaryPix;	//boundary pixels for each label
	std::vector< std::vector< typename LabelImageType::IndexType > > interiorPix;	//interior pixels for each label
	std::vector< std::map<TLPixel, int> > sharePix;		//number of edges shared between boundary pairs
								//the map will connect neighbors to number of edges shared.

	typedef std::map<TLPixel, IntrinsicFeatures> FeatureMapType;
	FeatureMapType featureVals;					//Holds all Features that have been calculated (including 0)
	std::map<TLPixel, int> LtoIMap;				//Map the label to the index in vectors that it is stored
	std::map<TLPixel, std::vector< std::vector<double> > > IDtoZernikeMap;
	std::vector< LabelPixelType > labels;		//Holds all of the Labels that have been found (including 0)
	
	//OPTIONS
	short int computationLevel;					//We have 3 levels of computation
	bool surfacecomputation;                    //Surface features needed
	bool computeHistogram;						//Requires Level 2
	bool computeTextures;						//Requires Level 3

};

}  // end namespace ftk

#include "ftkLabelImageToFeatures.txx"

#endif	// end __ftkLabelImageToFeatures_h

