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
#include "seed_detection/seedsdetection.h"
#include "seed_detection_3D/seedsdetection_3D.h"
#include "clustering_3D/local_max_clust_3D.h"
#include "clustering_2D/Local_Max_Clust_2D.h"
#include "graphColLearn_3D/Multi_Color_Graph_Learning_3D.h"
#include "alpha_expansion/alpha_expansion.h"

//ITK INCLUDES
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkScalarConnectedComponentImageFilter.h"

using namespace std;


void ucharToFloat(unsigned char* fromLoc, float* toLoc,int r, int c, int z, char invert);
void ucharToUShort(unsigned char* fromLoc, unsigned short* toLoc,int r, int c, int z, char invert);

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

class yousef_nucleus_seg 
{
public:
	//: constructor
	yousef_nucleus_seg();
  
	//: destructor
	~yousef_nucleus_seg();

	void setDataImage(unsigned char* imgPtr, int x, int y, int z, const char* filename);	//The image is loaded elsewhere and passed here.  I do not delete the data.
	void setParams(int *params);													//All parameters passed as integers, set the parameters accordingly
  
	unsigned char* getDataImagePtr(){ return dataImagePtr; };								
	unsigned short* getBinImage(){ return binImagePtr; }; 
	unsigned short* getSeedImage(){ return seedImagePtr; };
	float* getLogImage(){ return logImagePtr; };
	unsigned short* getClustImage(){ return clustImagePtr; };
	unsigned short* getSegImage(){ return segImagePtr; };
	std::vector<int> getImageSize();		// Returns in form [0]numStacks, [1]numRows, [2]numColumns

	//Return a list of seeds detected during Seeds Detection Function
	vector<Seed> getSeeds(){ return mySeeds; };
	void outputSeeds(void);

	//sub-modules that can be executed
	void runBinarization();
	void runSeedDetection();
	void runClustering();
	void runAlphaExpansion();
	void runAlphaExpansion3D();
	void runAlphaExpansion2D();

	void readParametersFromFile(const char* pFName);					//This function reads the parameters based on Yousef's parameter format

	//Come from parameters
	int isSegmentationFinEnabled() { return finalizeSegmentation;} ;
	int getSamplingRatio() { return sampling_ratio_XY_to_Z;} ;
	int isUseDapEnabled() { return useDistMap;} ;

	//Convert the label image to the previous IDL format
	int saveIntoIDLFormat(std::string imageName);
	//Read the output from the previous IDL format
	int readFromIDLFormat(std::string fileName);
	//return a list of seeds
	vector<Seed> getSeedsList() { return mySeeds; }

	//
	std::vector< int > SplitInit(ftk::Object::Point P1, ftk::Object::Point P2);
	ftk::Object::Point MergeInit(ftk::Object::Point P1, ftk::Object::Point P2, int *newID);
	bool DeleteInit(ftk::Object::Point P1);
	int getMaxID(int);

private:	
	void ExtractSeeds();
	void getConnCompInfo3D();
	int getConnCompImage(unsigned short* IM, int connectivity, int minSize, int r, int c, int z,int runConnComp);
	int getRelabeledImage(unsigned short* IM, int connectivity, int minSize, int r, int c, int z,int runConnComp);

	void clearBinImagePtr();
	void clearSeedImagePtr();
	void clearLogImagePtr();
	void clearSegImagePtr();
	void clearClustImagePtr();
	void clearMyConnComp();

	std::vector< ftk::Object::Point > getObjectBoundingBox(int id, int Int_Fin);

	//Internal Image information
	unsigned char* dataImagePtr;	//Created outside yousef_seg
	string dataFilename;
	unsigned short* binImagePtr;				//Created in yousef_seg
	unsigned short* seedImagePtr;				//Created in yousef_seg
	float* logImagePtr;				//Created in yousef_seg
	unsigned short* clustImagePtr;				//Created in yousef_seg
	unsigned short* segImagePtr;				//Created in yousef_seg
	
	//Size of all images above
	int numStacks;
	int numRows;
	int numColumns;

	vector<Seed> mySeeds;
	ConnComp* myConnComp;	//added by Yousef on 05-21-2008
	int numConnComp;		//added by Yousef on 05-21-2008
	int minLoGImg;			//minimum value from Laplacian of Gaussian image
	int numObjects;

	TParamsEntry* m_pData;
	
	//parameters
	int shift;
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

