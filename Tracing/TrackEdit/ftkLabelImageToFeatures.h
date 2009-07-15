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
#ifndef __ftkLabelImageToFeatures_h
#define __ftkLabelImageToFeatures_h

#include <itkLightObject.h>
#include <itkObjectFactory.h>

#include "itkLabelGeometryImageFilter.h"
#include <itkLabelStatisticsImageFilter.h>
#include <itkGradientMagnitudeImageFilter.h>

namespace ftk
{

typedef struct
{
	//Dirk's itkLabelGeometryImageFilter Features: (stored, others used in other calculations (level 1)
	unsigned int volume;
	unsigned int integratedintensity;
	std::vector<float> centroid;
	std::vector<float> weightedcentroid;
	std::vector<float> axislengths;
	float eccentricity; 
	float elongation;
	float orientation;
	std::vector<int> boundingbox;		// [min(X), max(X), min(Y), max(Y), ... ]
	float bbVolume;

	//itkLabelStatisticsImageFilter Features: (level 2)
	float sum;
	float mean;
	float median;	//advanced
	float minimum;
	float maximum;
	float sigma;
	float variance;

	//Others:	(level 3)
	float solidity;
	float texture;
	float radiusvariation;
	float skew;
	float energy;
	float entropy;
	float surfacegradient;
	float interiorgradient;
	float interiorintensity;
	float surfaceintensity;
	float intensityratio;
	float percentsharedboundary;
	float surfacearea;
	float shape;
	float averagedistance;
	float distancevariation;

	//User set
	int num;
	int time;
	int tag;
	std::vector<float> spacing;
	static int packed_size()
	{
		return 49;
	}
	void pack(float *p)
	{
		*(p++) = volume;
		*(p++) = integratedintensity;
		*(p++) = eccentricity;
		*(p++) = elongation;
		*(p++) = orientation;
		*(p++) = bbVolume;
		*(p++) = sum;
		*(p++) = mean;
		*(p++) = median;
		*(p++) = minimum;
		*(p++) = maximum;
		*(p++) = sigma;
		*(p++) = variance;
		*(p++) = solidity;
		*(p++) = texture;
		*(p++) = radiusvariation;
		*(p++) = skew;
		*(p++) = energy;
		*(p++) = entropy;
		*(p++) = surfacegradient;
		*(p++) = interiorintensity;
		*(p++) = surfaceintensity;
		*(p++) = intensityratio;
		*(p++) = percentsharedboundary;
		*(p++) = surfacearea;
		*(p++) = shape;
		*(p++) = averagedistance;
		*(p++) = distancevariation;
		*(p++) = spacing[0];//FIXME: assumed 3-D
		*(p++) = spacing[1];
		*(p++) = spacing[2];
		*(p++) = centroid[0];
		*(p++) = centroid[1];
		*(p++) = centroid[2];
		*(p++) = weightedcentroid[0];
		*(p++) = weightedcentroid[1];
		*(p++) = weightedcentroid[2];
		*(p++) = axislengths[0];
		*(p++) = axislengths[1];
		*(p++) = axislengths[2];
		*(p++) = boundingbox[0];
		*(p++) = boundingbox[1];
		*(p++) = boundingbox[2];
		*(p++) = boundingbox[3];
		*(p++) = boundingbox[4];
		*(p++) = boundingbox[5];
		*(p++) = num;
		*(p++) = time;
		*(p++) = tag;
	}
	void unpack(float *p)
	{
		
  		volume 			= (unsigned int)(*(p++)+0.5);
  		integratedintensity	= (unsigned int)(*(p++)+0.5);
  		eccentricity 		= *(p++);
  		elongation 		= *(p++);
  		orientation 		= *(p++);
  		bbVolume 		= *(p++);
  		sum 			= *(p++);
  		mean 			= *(p++);
  		median 			= *(p++);
  		minimum 		= *(p++);
  		maximum 		= *(p++);
  		sigma 			= *(p++);
  		variance 		= *(p++);
  		solidity 		= *(p++);
  		texture 		= *(p++);
  		radiusvariation 	= *(p++);
  		skew 			= *(p++);
  		energy 			= *(p++);
  		entropy 		= *(p++);
  		surfacegradient 	= *(p++);
  		interiorintensity	= *(p++);
  		surfaceintensity 	= *(p++);
  		intensityratio 		= *(p++);
  		percentsharedboundary   = *(p++);
  		surfacearea 		= *(p++);
  		shape 			= *(p++);
  		averagedistance 	= *(p++);
  		distancevariation	= *(p++);
		spacing.clear();
  		spacing.push_back(*(p++));
  		spacing.push_back(*(p++));
		spacing.push_back(*(p++));
		centroid.clear();
  		centroid.push_back(*(p++));
  		centroid.push_back(*(p++));
  		centroid.push_back(*(p++));
		weightedcentroid.clear();
  		weightedcentroid.push_back(*(p++));
  		weightedcentroid.push_back(*(p++));
  		weightedcentroid.push_back(*(p++));
		axislengths.clear();
  		axislengths.push_back(*(p++));
  		axislengths.push_back(*(p++));
  		axislengths.push_back(*(p++));
		boundingbox.clear();
  		boundingbox.push_back((int)(*(p++)+0.5));
  		boundingbox.push_back((int)(*(p++)+0.5));
  		boundingbox.push_back((int)(*(p++)+0.5));
  		boundingbox.push_back((int)(*(p++)+0.5));
  		boundingbox.push_back((int)(*(p++)+0.5));
  		boundingbox.push_back((int)(*(p++)+0.5));
  		num 			= (int)(*(p++)+0.5);
  		time 			= (int)(*(p++)+0.5);
  		tag 			= (int)(*(p++)+0.5);
		//delete [] p; //FIXME
	}
	int Sprintf(char *s,int level =1)
	{
		int bn;
		if(level >=1)
		{
			//ORDER:
			//num tag dims time spacing centroid axislengths eccentricity elongation boundingbox
			bn = sprintf(s," %d %d %d %d",num,tag,centroid.size(),time);
			for (int counter=0; counter< spacing.size(); counter++)
			{
				bn += sprintf(&s[bn]," %f", spacing[counter]);
			}
		        for (int counter=0; counter< centroid.size(); counter++)
			{
				bn += sprintf(&s[bn]," %0.3f",centroid[counter]);
			}
			for(int counter=0; counter< axislengths.size(); counter++)
			{
				bn += sprintf(&s[bn]," %0.3f",axislengths[counter]);
			}
			for(int counter=0; counter< boundingbox.size(); counter++)
			{
				bn += sprintf(&s[bn], " %d", boundingbox[counter]);
			}
			bn += sprintf(&s[bn]," %0.3f",bbVolume);
		}
		if(level >=2)
		{
			//ORDER:
			//mean, sigma
			bn += sprintf(&s[bn]," %0.3f %0.3f",mean,sigma);
		}
	}
	void Fprintf(FILE *fp = stdout,int level=1)
	{
		if(level >=1)
		{
			//ORDER:
			//num tag dims time spacing centroid axislengths eccentricity elongation boundingbox
			fprintf(fp,"%d %d %d %d",num,tag,centroid.size(),time);
			for (int counter=0; counter< spacing.size(); counter++)
			{
				fprintf(fp," %f", spacing[counter]);
			}
		        for (int counter=0; counter< centroid.size(); counter++)
			{
				fprintf(fp," %0.3f",centroid[counter]);
			}
			for(int counter=0; counter< axislengths.size(); counter++)
			{
				fprintf(fp," %0.3f",axislengths[counter]);
			}
			for(int counter=0; counter< boundingbox.size(); counter++)
			{
				fprintf(fp, " %d", boundingbox[counter]);
			}
			fprintf(fp," %0.3f",bbVolume);
		}
		if(level >=2)
		{
			//ORDER:
			//mean, sigma
			fprintf(fp," %0.3f %0.3f",mean,sigma);
		}
		fprintf(fp,"\n");
	}

/* FIXME */
/*	int Sscanf(char* s,int level=1)
	{
		int ret;
		if(level >=1)
		{
			int dims,ival;
			float fval;
			sscanf(s,"%d %d %d %d", &num, &tag, &dims, &time);
			for(int counter=0; counter < dims; counter++)
			{
				sscanf(s,"%f",&fval);
				spacing.push_back(fval);
			}
			for(int counter=0; counter < dims; counter++)
			{
				sscanf(s,"%f", &fval);
				centroid.push_back(fval);
			}
			for(int counter=0; counter < dims; counter++)
			{
				ret = sscanf(s,"%f", &fval);
				axislengths.push_back(fval);
			}
			for(int counter=0; counter < 2*dims; counter++)
			{
				sscanf(s,"%d",&ival);
				boundingbox.push_back(ival);
			}
			sscanf(s,"%f",bbVolume);
		}
		if(level >=2)
		{
			ret = sscanf(s," %f %f", &mean, &sigma);
		}
		return ret;
	} */

	int Fscanf(FILE *fp = stdin,int level=1)
	{
		int ret;
		if(level >=1)
		{
			int dims,ival;
			float fval;
			fscanf(fp,"%d %d %d %d", &num, &tag, &dims, &time);
			for(int counter=0; counter < dims; counter++)
			{
				fscanf(fp,"%f",&fval);
				spacing.push_back(fval);
			}
			for(int counter=0; counter < dims; counter++)
			{
				fscanf(fp,"%f", &fval);
				centroid.push_back(fval);
			}
			for(int counter=0; counter < dims; counter++)
			{
				ret = fscanf(fp,"%f", &fval);
				axislengths.push_back(fval);
			}
			for(int counter=0; counter < 2*dims; counter++)
			{
				fscanf(fp,"%d",&ival);
				boundingbox.push_back(ival);
			}
			fscanf(fp,"%f",bbVolume);
		}
		if(level >=2)
		{
			ret = fscanf(fp," %f %f", &mean, &sigma);
		}
		return ret;
	}

} LabelImageFeatures;

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
		LabelImageFeatures GetFeatures( LabelPixelType label );
		std::vector< LabelPixelType > GetLabels() { return this->labels; };

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
		void InitFeatureMap();
		LabelImageFeatures GetEmptyFeatures();
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

		typedef std::map<TLPixel, LabelImageFeatures> FeatureMapType;

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

		std::vector< LabelPixelType > labels;		//Holds all of the Labels that have been found (including 0)
		FeatureMapType allFeatures;					//Holds all Features that have been calculated (including 0)

		//OPTIONS
		short int computationLevel;					//We have 3 levels of computation
		bool computeHistogram;						//Requires Level 2
		bool computeAdvanced;						//Requires Level 3

};

}  // end namespace ftk

#include "ftkLabelImageToFeatures.txx"

#endif	// end __ftkLabelImageToFeatures_h

