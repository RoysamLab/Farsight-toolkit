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

//#pragma warning(disable:4786)

#ifndef main_h
#define main_h

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <ctime>
#include <list>

//#include "itkVector.h"
//#include "itkListSample.h"
//#include "itkKdTree.h"
//#include "itkWeightedCentroidKdTreeGenerator.h"
//#include "itkKdTreeBasedKmeansEstimator.h"
//#include "itkMinimumDecisionRule.h"
//#include "itkEuclideanDistance.h"
//#include "itkSampleClassifier.h"
//#include "itkMeanCalculator.h"
//#include "itkCovarianceCalculator.h"

#include "timer.h"
#include "CONSTANTS.h"
#include "Mytypes.h"
#include "Config.h"
#include "Cpoints.h"
#include "Dllist.h"
#include "Cimage.h"    // a simple image class
#include "Cvessel.h"
#include "C3dimage.h"
#include "Cvector.h"
#include "Template.h"  // a template class
#include "Global.h"    // all global variables
#include "Ctree.h"
#include "Soma.h"
#include "StrEle.h"

extern void Draw_BranchPoints();
extern C3DImage* WobblyLineImage(int diameter);

extern CSomas *gTheSomas;

int giSomaVolume = 0;
CPoint *gaSomaPoints = NULL;
int giNumOfLeadingZeros = 0;

///////////////////////////
// Function Declarations //
///////////////////////////

extern double compute_Q(C3DImage * VesselnessImage, double alpha);

void LocateSomas();
void LocateSomas2();
void LocateSomas3();
//By Yousef: Try this ////////////
void LocateSomas3_v2();
void LocateSomas3_v3();
void LocateSomas3_v4();
//////////////////////////////////
extern void FreeResources();
extern void ProcessCommandLine(int argc, char *argv[]);
extern C3DImage * Close(C3DImage * source, StrEle *sDilate, StrEle *sErode=NULL);

extern void Initialize();
extern void DisplayVectors(char *fName);
extern void DisplayTemplates(char *fName);
//extern int FindSeedPoints(C3DImage &anImage, 
//									CQueue<CPoint> &aQueue, int sliceNum);

extern int FindSeedPoints2(CImage &anImage);
//extern void FindSeedPoints2(CImage& anImage, std::vector<CPoint> & seed_candidates);
extern void EstimateZCoord(std::vector<CPoint> & seed_candidates);
extern void EstimateXCoord(std::vector<CPoint> & seed_candidates);
extern void EstimateYCoord(std::vector<CPoint> & seed_candidates);
extern int ReadSeedPointsFromAFile(std::string filename);


extern void TraceA3DImage(C3DImage &);
extern int ConvertFromHex(unsigned char aChar, int Order = 0);
extern void DetectSomas(unsigned char **, unsigned char **);
extern void DetectSomas(CImage *, CImage *, int);
extern void DetectSomas(CImage *, CImage *, CPoint *, int);
extern void MarkSomas(C3DImage &);

extern bool VerifySeedPoint3D2(CPoint *aPoint);
extern void VerifyAndSortAllSeedPoints();
//extern void VerifyAndSortAllSeedPoints(const std::vector<CPoint> & seed_candidates);

extern void EstimateMedianAndStdDev(int *Hist, int &median, float &stdDev);
extern void DrawBall(C3DImage &anImage, CPoint &center, int radius, unsigned char color = 255);
extern void DrawFilledCircle(CImage &anImage, CPoint &center, int radius, int color = 255);
extern void	ConstructStructElem(int Radius, CPoint *aDiskPoints, int &iNumOfPoints);

extern void SetConfigParameters();
extern void SetImageVariables(char *);
extern void PrintUsageMessage();

// functions declared in draw.cpp
extern void Draw_SeedCandidates (std::vector<CPoint> & seed_candidates);
extern void Draw_Projections();
extern void Draw_CenterlineOnProjections();
extern void Draw_BorderlineOnProjections(char * plane);
extern void Draw_SeedPointsOnProjections(char * gridlines = NULL);
extern void Draw_PointsOnProjections(list<CPoint> center, list<CPoint> left,
							  list<CPoint> right, char * name = NULL);
extern void Draw_PointsOnProjections(list<CPoint> points, char *name = NULL);
extern void Draw_SeedPointDirections();
extern void Draw_3DVessels();
extern void Draw_Binary3DImage();
extern void Draw_Centerline3D();
extern void Draw_Residual3DImage();

CPoint **gapSortedArrayOfSeedPoints= NULL;
int giNumOfSeedPoints = 0;
extern void SortSeedPoints();

extern int giSeedPointsHistogram[256];
// for debugging purposes
//#define DEBUG 1
#define DEBUG 0
//#define Display_DEBUG 1
#define Display_DEBUG 0
//#define ReadNegateWrite 1
#define ReadNegateWrite 0


int giHitSomaFlag = 0;

ofstream gLogFile;

// the following two lines were added for CANCER images
float gfHWidth = 0.0;
float gfVWidth = 0.0;

#endif
