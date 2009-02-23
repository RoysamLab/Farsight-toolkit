/*=========================================================================

  Program:   Farsight Biological Image Segmentation and Visualization Toolkit
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision: 0.00 $

=========================================================================*/
#ifndef __ftkLabelImageToFeatures_h
#define __ftkLabelImageToFeatures_h

#include <itkLightObject.h>
#include <itkObjectFactory.h>

#include "itkLabelGeometryImageFilter.h"
#include <itkLabelStatisticsImageFilter.h>
#include <itkGradientMagnitudeImageFilter.h>

#include "ftkFeatures.h"

#include <iostream>
#include <map>

namespace ftk
{

//********************************************************************************************************
//THIS IS THE CLASS THAT DOES THE CALCULATIONS (THE ENGINE)
//********************************************************************************************************
template< typename TIPixel = unsigned char, typename TLPixel = unsigned short, unsigned int VImageDimension = 2> 
class LabelImageToFeatures : public itk::LightObject
{
public:

	typedef LabelImageToFeatures Self;
	typedef LightObject Superclass;
	typedef itk::SmartPointer< Self > Pointer;
	typedef itk::SmartPointer< const Self > ConstPointer;

	typedef TIPixel IntensityPixelType;
	typedef TLPixel LabelPixelType;
	typedef float FloatPixelType;
	typedef itk::Image< IntensityPixelType, VImageDimension > IntensityImageType;
	typedef itk::Image< LabelPixelType, VImageDimension > LabelImageType;
	typedef itk::Image< FloatPixelType, VImageDimension > FloatImageType;
	typedef typename IntensityImageType::Pointer IntensityImagePointer;
	typedef typename LabelImageType::Pointer LabelImagePointer;
	typedef typename FloatImageType::Pointer FloatImagePointer;

	itkNewMacro( Self );
	
	itkTypeMacro(LabelImageToFeatures, LightObject);

	bool SetImageInputs( IntensityImagePointer intImgIn, LabelImagePointer lblImgIn );
	void Update();
	LabelPixelType GetMaxLabel();
	float GetPercentSharedBoundary(TLPixel focusLabel, TLPixel neighborLabel);
	std::vector<TLPixel> GetContactNeighbors(TLPixel label);
	std::vector<float> GetCentroid(TLPixel label);
	std::vector<float> GetWeightedCentroid(TLPixel label);
	std::vector<float> GetAxisLengths(TLPixel label);
	std::vector<int> GetBoundingBox(TLPixel label);
	LabelImageFeatureValueMapType GetFeatures( LabelPixelType label );
	std::vector< LabelPixelType > GetLabels() { return this->labels; };
	std::vector< std::string > GetAvailableFeatureNames(void);
	LabelImageFeatureInfoMapType GetFeatureInfo(void) { return this->featureInfo; };

	void ComputeHistogramOn();
	void ComputeHistogramOff(){ computeHistogram = false; };
	void ComputeAdvancedOn();
	void ComputeAdvancedOff(){ computeAdvanced = false; };
	void SetLevel(short int newLevel);
	short int GetLevel(){ return computationLevel; };

protected:
	LabelImageToFeatures();
	~LabelImageToFeatures(){};

private:
	LabelImageToFeatures(const Self&);  //purposely not implemented
	void operator=(const Self&);		//purposely not implemented

	//Internal Functions:
	bool CreateGradientMagnitudeImage();
	bool RunLabelGeometryFilter();
	bool RunLabelStatisticsFilter();
	void LabelImageScan();
	void ReadLabelGeometryFeatures();
	void ReadLabelStatisticsFeatures();
	void CalculateScanFeatures();
	void CalculateHistogramFeatures();
	void SetHistogramParameters(int* numBins, int* lowerBound, int* upperBound);

	//Internal types:
	typedef itk::LabelGeometryImageFilter< LabelImageType, IntensityImageType > LabelGeometryType;
	typedef typename LabelGeometryType::Pointer LabelGeometryPointer;
	typedef itk::LabelStatisticsImageFilter< IntensityImageType , LabelImageType > LabelStatisticsType;
	typedef typename LabelStatisticsType::Pointer LabelStatisticsPointer;
	

	//Internal Variables:
	IntensityImagePointer intensityImage;	//Input intensity image;
	LabelImagePointer labelImage;			//Input label image;
	FloatImagePointer gmImage;				//Calculated gradient magnitude image;

	LabelGeometryPointer labelGeometryFilter;		//Dirk's Filter;
	LabelStatisticsPointer labelStatisticsFilter;	//ITK label statistics filter

	std::vector< std::vector< typename LabelImageType::IndexType > > boundaryPix;	//boundary pixels for each label
	std::vector< std::vector< typename LabelImageType::IndexType > > interiorPix;	//interior pixels for each label
	std::vector< std::vector< LabelPixelType > > sharePix;				//number of edges shared between boundary pairs
																		//Values stored once with greater label first (outer array)
																		//example: sharePix[5][3] is correct, sharePix[3][5] does not exist

	typedef std::map<TLPixel, LabelImageFeatureValueMapType> FeatureMapType;
	std::vector< LabelPixelType > labels;		//Holds all of the Labels that have been found (including 0)
	FeatureMapType featureVals;					//Holds all Features that have been calculated (including 0)
	LabelImageFeatureInfoMapType featureInfo;	//Holds the Feature Info for features calculated here!!


	//OPTIONS
	short int computationLevel;					//We have 3 levels of computation
	bool computeHistogram;						//Requires Level 2
	bool computeAdvanced;						//Requires Level 3

};

}  // end namespace ftk

#include "ftkLabelImageToFeatures.txx"

#endif	// end __ftkLabelImageToFeatures_h

