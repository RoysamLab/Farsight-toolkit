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

#ifndef _YOUSEF_SEG_H_
#define _YOUSEF_SEG_H_

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

//STD INCLUDES
#include <iostream>
#include <list>
#include <vector>

//FTK INCLUDES
#include <ftkImage/ftkImage.h>
#include <ftkObject.h>

//LOCAL INCLUDES
#include "cell_binarization/cell_binarization.h"
#include "cell_binarization/adaptive_binarization.h"
#include "seed_detection/seedsdetection.h"
#include "seed_detection_3D/seedsdetection_3D.h"
#include "clustering_3D/local_max_clust_3D.h"
#include "clustering_2D/Local_Max_Clust_2D.h"
#include "graphColLearn_3D/Multi_Color_Graph_Learning_3D.h"
#include "alpha_expansion/alpha_expansion.h"
#include "EM_GMM/EM_Project_3_cpp_3D_comp.h"
//#include "gradient_weighted_distance_2D/gradientweighteddistancemap.h" //Still under testing
#include "preprocessing/GrAnisDiff.h" //Still under testing

//ITK INCLUDES
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkScalarConnectedComponentImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"

// FOR GMM ESTIMATION
#include "itkVector.h"
#include "itkListSample.h"
#include "itkGaussianMixtureModelComponent.h"
#include "itkExpectationMaximizationMixtureModelEstimator.h"
#include "itkNormalVariateGenerator.h"

//#include "tmp_itk_rev/itkBinaryMorphologicalClosingImageFilter.h"


void ucharToFloat(unsigned char* fromLoc, float* toLoc,int r, int c, int z, char invert);
void ucharToUShort(unsigned char* fromLoc, unsigned short* toLoc,int r, int c, int z, char invert);
void ushortToFloat(unsigned short* fromLoc, float* toLoc,int r, int c, int z, char invert);

unsigned char *** TriplePtr(int z, int r, int c);
  
class Seed;	//At bottom of file

typedef struct _ParamsEntry
{
	char m_pName[128];
	int m_pValue;
} TParamsEntry;

struct ConnComp
{
	int x1;
	int x2;
	int y1;
	int y2;
	int z1;
	int z2;
};

struct GMMParams
{
	double background_mean;
	double background_stdev;
	double foreground_mean;
	double foreground_stdev;
	double background_prior;
	double foreground_prior;
};

class yousef_nucleus_seg 
{
public:
	//: constructor
	yousef_nucleus_seg();
  
	//: destructor
	~yousef_nucleus_seg();

	void setDataImage(unsigned char* imgPtr, int x, int y, int z, const char* filename);	//The image is loaded elsewhere and passed here.  I do not delete the data.
	void setDataImage(unsigned short* imgPtr, int x, int y, int z, const char* filename);	//The image is loaded elsewhere and passed here.  I do not delete the data.
	void setInitialBinaryImage(unsigned short* imgPtr);	//The image is loaded elsewhere and passed here.  I do not delete the data.
	void setParams(int *params);		//All parameters passed as integers, set the parameters accordingly
  	void setParamsForSeedDetection(int highsensitivity, double sMin, double sMax, double rXY,  double rZ, int usedistMap, int samplingRatio, int minSize = 50);
	unsigned char* getDataImagePtr(){ return dataImagePtr; };								
	unsigned short* getBinImage(){ return binImagePtr; }; 
	void setBinImage(unsigned short* ptr); 
	unsigned short* getSeedImage(){ return seedImagePtr; };
	float* getLogImage(){ return logImagePtr; };
	unsigned short* getClustImage(){ return clustImagePtr; };
	unsigned short* getSegImage(){ return segImagePtr; };
	std::vector<int> getImageSize();		// Returns in form [0]numStacks, [1]numRows, [2]numColumns	

	//Return a list of seeds detected during Seeds Detection Function
	std::vector<Seed> getSeeds(){ return mySeeds; };
	void outputSeeds(void);

	//sub-modules that can be executed
	void runBinarization(unsigned short number_of_bins = 128);
	bool runBinarization16();
	void runSeedDetection();
	void runSeedDetection16();
	void runSeedDetection(int minScale,int maxScale); // Added by Raghav
	void runClustering();
	void runClustering16();
	void runAlphaExpansion();
	void runAlphaExpansion16();
	void runAlphaExpansion3D();
	void runAlphaExpansion2D();
	void runAlphaExpansion2D16();
	//added by Yousef on 9/14/2009
	//void runGradWeightedDistance();
	//added by Yousef on 11/11/2009
	void runGradAnisDiffSmoothing();
	bool EstimateGMMParameters();
	bool RunKmeansClustering();


	void readParametersFromFile(const char* pFName);					//This function reads the parameters based on Yousef's parameter format
	void writeParametersToFile();

	//Come from parameters
	int getShift(){ return shift; };
	int isAdaptiveBinEnabled(){ return adaptive_bin; };
	int getSigma(){ return sigma; };
	int getScaleMin(){ return scaleMin; };
	int getScaleMax(){ return scaleMax; };
	int getRegionXY(){ return regionXY; };
	int getRegionZ(){ return regionZ; };
	int isSegmentationFinEnabled() { return finalizeSegmentation;} ;
	int getSamplingRatio() { return sampling_ratio_XY_to_Z;} ;
	int isUseDapEnabled() { return useDistMap;} ;
	int getRefineRange(){ return refineRange; };
	int getMinObjSize(){ return minObjSize; };

	//Convert the label image to the previous IDL format
	int saveIntoIDLFormat(std::string imageName);
	//Read the output from the previous IDL format
	int readFromIDLFormat(std::string fileName);
	//return a list of seeds
	std::vector<Seed> getSeedsList() { return mySeeds; }

	//
	std::vector< int > SplitInit(ftk::Object::Point P1, ftk::Object::Point P2);
	ftk::Object::Point MergeInit(ftk::Object::Point P1, ftk::Object::Point P2, int *newID);
	bool DeleteInit(ftk::Object::Point P1);
	int getMaxID(int);
	//added by Yousef on 9/11/2009
	int AddObject(unsigned char* inImage, unsigned short* lbImage, std::vector<int> P1, std::vector<int> P2,
									std::vector<itk::SizeValueType> imSZ, int maxID);
	//added by Yousef on 9/26/2009
	int AddObject2D(unsigned char* inImage, unsigned short* lbImage, std::vector<int> P1, std::vector<int> P2,
									std::vector<itk::SizeValueType> imSZ, int maxID);
	

private:	
	void ExtractSeeds();
	void getConnCompInfo3D();
	int getConnCompImage(unsigned short* IM, int connectivity, int minSize, size_t r, size_t c, size_t z,int runConnComp);
	int getRelabeledImage(unsigned short* IM, int connectivity, int minSize, int r, int c, int z,int runConnComp);
	
	void clearBinImagePtr();
	void clearSeedImagePtr();
	void clearLogImagePtr();
	void clearSegImagePtr();
	void clearClustImagePtr();
	void clearMyConnComp();

	std::vector< ftk::Object::Point > getObjectBoundingBox(int id, int Int_Fin);

	//added by Yousef on 09/06/2009
	void fitMixGaussians();	

	//Internal Image information
	unsigned char* dataImagePtr;	//Created outside yousef_seg
	unsigned short* dataImagePtr16;	//Created outside yousef_seg
	unsigned short* initial_binaryImage;	//Created outside yousef_seg
	std::string dataFilename;
	unsigned short* binImagePtr;				//Created in yousef_seg
	unsigned short* seedImagePtr;				//Created in yousef_seg
	float* logImagePtr;				//Created in yousef_seg
	unsigned short* clustImagePtr;				//Created in yousef_seg
	unsigned short* segImagePtr;				//Created in yousef_seg
	
	//Size of all images above
	size_t numStacks;
	size_t numRows;
	size_t numColumns;

	std::vector<Seed> mySeeds;
	ConnComp* myConnComp;	//added by Yousef on 05-21-2008
	int numConnComp;		//added by Yousef on 05-21-2008
	int minLoGImg;			//minimum value from Laplacian of Gaussian image
	int numObjects;

	TParamsEntry* m_pData;
	GMMParams gmm_params_;
	
	//parameters
	int shift;
	int adaptive_bin;
	//added by vinay on 02/09/2012
	int sigma;
	double scaleMin;
	double scaleMax;
	double regionXY;
	double regionZ;
	int useDistMap;
	int finalizeSegmentation;
	int sampling_ratio_XY_to_Z;
	//added by yousef on 11/4/2008
	int refineRange;
	//added by yousef on 12/5/2008
	int minObjSize;
	//added by yousef on 09/02/2009
	bool autoParamEstimation;
};

class Seed
{
public:
	Seed(){
		xVal = 0;
		yVal = 0;
		zVal = 0;
		idVal = 0;
		ccVal = 0;
	}
	Seed( int x, int y, int z, int ID, int cc) {
		xVal = x;
		yVal = y;
		zVal = z;
		idVal = ID;
		ccVal = cc;
	}
	
	void setX(int x ){ xVal = x; }
	void setY(int y ){ yVal = y; }
	void setZ(int z ){ zVal = z; }
	void setID(int i){ idVal = i;}
	void setCC(int cc) { ccVal = cc;}

	int x() const { return xVal; }
	int y() const { return yVal; }
	int z() const { return zVal; }
	int ID() const {return idVal;}
	int CC() const {return ccVal;}

private:
	int xVal;
	int yVal;
	int zVal;
	int idVal;
	int ccVal;
};

#endif

