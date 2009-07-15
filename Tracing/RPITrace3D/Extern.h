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

#ifndef extern_h
#define extern_h

//	extern void Draw_XYCenterline(CVessel *current);
#include "Branches.h"

extern list<deque<CPoint> > Seeds;
extern bool*** TracedImage;

extern int HALFGRID;
extern ofstream gLogFile;

/////////////////
// operation parameters

extern int giMaxExtensionLength;
extern int giParam_Trace_MinSegmentLength;

// configuration items and their corresponding default values
extern char cfg_output_image_format[4];
extern int cfg_output_image_bgcolor;
extern int cfg_output_image_fgcolor;
extern int cfg_output_projections;
extern int cfg_output_soma_draw_trees;
extern int cfg_mode_debug;
extern int cfg_trace_nonparallel;
extern bool cfg_verify_amri;
extern bool cfg_seedsfromfile;
extern std::string cfg_seedsfilename;
// a copy of projection images for drawing purposes
extern CImage* CanvasXY;
extern CImage* CanvasXZ;
extern CImage* CanvasYZ;

///////////////////////////////////
// Image constants and variables //
///////////////////////////////////
//	extern char gachBasicFileName[];
extern char gachPath[];
extern char gachType[];
extern char gachFileOfSlices[];

extern unsigned short heated_object_colormap[];

extern int* giMedians;
extern float* gfStdDevs;
extern int giMedian;
extern float gfStdDev;
extern int giForegroundMedian;
extern float gfForegroundStdDev;
extern int giBackgroundMedian;
extern float gfBackgroundStdDev;

extern int giFirstSlice;
extern int giLastSlice;
extern int giSliceNumber;
// a margine from each side of the image. This is necessary because when
// considering a center point, we need space to shift the template.
extern int giMARGIN;

extern int giMinSeedPixelValue;

// a configuration object
extern CConfig gConfig;

extern list<CPoint> CenterPoints;
extern list<CPoint> HRPoints;
extern list<CPoint> HLPoints;
extern list<CPoint> VRPoints;
extern list<CPoint> VLPoints;
extern list<CPoint> TracedPoints;
extern list<CPoint> UnverifiedSeeds;
extern list<CPoint> VerifiedSeedsCenter;
extern list<CPoint> VerifiedSeedsRight;
extern list<CPoint> VerifiedSeedsLeft;
extern list<CPoint> AmriVerifiedSeedsCenter;
extern list<CPoint> TracedSeeds;

//extern int giXYDistBetweenSeeds;

// a variable holding the number of the image being currently
// processed
extern int giContrastThreshold;
extern long glSumOfSeedPointValues;
extern long glNumOfSeedPoints;
extern int gaiForegroundHistogram[256];
extern int gaiBackgroundHistogram[256];
extern int giSmallestResponseAccepted;
extern float gfOrthogonalResponsePercentage;

//
// the image that is being traked
extern CImage* TheTargetImage;
extern CImage* EdgeImage;
//	extern CImage *TrackImage;
extern CImage* ProjectionImage;
extern CImage* GlobalImage;
extern CImage* CenterlinesImage;
extern CImage* BoundariesImage;
extern CImage* BranchingPointsImage;

extern CImage* TrackImageXY;
extern CImage* TrackImageXZ;
extern CImage* TrackImageYZ;

// a 3D image that will hold the stack of images to be tracked. The tracking results
// will also be displayed on this image.
extern C3DImage* The3DImage;

//By Yousef
//This image will hold the 3D soma image. The image will be initialized to zeros, and will be filled with oned
//at the locations of soma
extern C3DImage* Soma3DImage;
//By YOusef
//Here is another image to hold the soma labels
extern C3DImage* SomaLabelsImage;
//By Yousef
//Use this object to detect branch points
extern CBranches* BranchPoints;

// an empty 3D image that is used to register the tracked vessels and seed points.
// this image is needed to avoid tracking a vessel more than one time.
extern C3DImage* Empty3DImage;

///////////
// the following two variables are used to make sure that we do
// not overstep the 3D-image boundaries
extern unsigned char* gpuchThe3DImageData;
extern long glThe3DImageLength;


//////////////////////////////////////
// Template constants and variables //
//////////////////////////////////////

// define 4-2D arrays of templates. left/right-horizontal and left/right-vertical
// The direction of each template is defined by the angles it makes with the
// positive x-axis in the xy plane (Hdir) and the angle it makes with the positive
// z-axis (Vdir)
extern CTemplate* gHLeftTemplatesArray[NumOfDirections][NumOfDirections];
extern CTemplate* gVLeftTemplatesArray[NumOfDirections][NumOfDirections];
extern CTemplate* gHRightTemplatesArray[NumOfDirections][NumOfDirections];
extern CTemplate* gVRightTemplatesArray[NumOfDirections][NumOfDirections];

//////////////////////////////////////
//  Vector constants and variables  //
//////////////////////////////////////

// an array of shift vectors, one for each direction
extern CVector* gVectorsArray[NumOfDirections][NumOfDirections];

////////////////////////////////////////
//  Tracking constants and variables  //
////////////////////////////////////////

//	extern CQueue<CPoint> *InitialPoints;   		  // an array of empty queue of points
// a 2D array of pointers to seed points
extern CPoint*** gapArrayOfSeedPoints; 
// a 2D array of filling points
//extern CDLList<CFillingPoint> ***ArrayOfFillings;


//////////////////////////////////////
//  	   Neurolucida Stuff		//
//////////////////////////////////////

extern char gachSpaces[1000];
extern char gachLeadingSpaces[1000];
extern int giNumOfLeadingSpaces;

//////////////////////////////////////
//  		   Functions			//
//////////////////////////////////////

extern float minimum(float a, float b);
extern int ReverseTheta(int dir);
extern int Round(double value);
extern bool VerifySeedPoint3D2(CPoint* aPoint);
extern int DirectionMinus(int d, int v);
extern int DirDistance(int d1, int d2);
extern int DirectionPlus(int d, int v);
extern void PrepareDirectionsArray(int* directions, int size, int dir);
// draw a line from point-a to point-b in the image
extern void draw_line(CPoint& a, CPoint& b, CImage& image, unsigned char color);
extern void draw_line(CPoint& a, CPoint& b, C3DImage& image, unsigned char color,
	int SliceNum);
extern void Add_line(CPoint& a, CPoint& b, CImage& image);
extern void DrawThickLine(CPoint* a, CPoint* b, CImage& image, unsigned char color);
extern void DrawThickLine(CPoint* a, CPoint* b, C3DImage& image, unsigned char color,
	int sliceNum);
extern void FindAndAddMoreSeedPoints(CPoint& aPoint);
extern int ConvertFromHex(unsigned char aChar, int Order = 0);
extern void EstimateMedianAndStdDev(int* Hist, int& median, float& stdDev);
extern int median(priority_queue<short> responses);
extern int median(list<int>& responses);
extern float Median(list<float>& responses);
extern void DrawBall(C3DImage& anImage, CPoint& center, int radius,
	unsigned char color = 255);
extern void DetectSomas(unsigned char** inImageData,
	unsigned char** outImageData);
extern void ConstructStructElem(int Radius, CPoint* aDiskPoints,
	int& iNumOfPoints);
extern void Construct3DLine(CPoint& from, CPoint& to, CDLList<CPoint>& result);
extern void DetectSomas(CImage* InImage, CImage* OutImage,
	CPoint* aDiskPoints, int iNumOfPoints);

extern void GetCenterLocation(CPoint& center_point, const CPoint& HR,
	const CPoint& HL, const CPoint& VR, const CPoint& VL);
extern void GetCenterDirection(CPoint& center_point, const CPoint& HR,
	const CPoint& HL, const CPoint& VR, const CPoint& VL);

// perform matrix multiplication
void MatrixMultiply(double [][4], double [][4], double [][4], int, int);


//Yousef
extern void Construct3DStructElem(int Radius, CPoint* aSphrPoints, int &iNumOfPoints);
extern void Detect3DSomas(C3DImage* InImage, C3DImage* OutImage, CPoint* sSphrPoints, int iNumOfPoints, int StructElemSize);


//////////////////////////////	
//	  Response Functions	//
//////////////////////////////
extern int ProcessResponseArrays(CPoint** pHR, CPoint** pHL, CPoint** pVR,
	CPoint** pVL, int& TemplateLength);
extern int GetAverageResponse(const CPoint& aPoint, CPoint& pHR, CPoint& pHL,
	CPoint& pVR, CPoint& pVL, int& TemplateLength);
extern int GetMedianResponse(const CPoint& aPoint, CPoint& pHR, CPoint& pHL,
	CPoint& pVR, CPoint& pVL, int& TemplateLength);
extern bool IndividualAverageResponse(int prev_shift, CPoint& point,
	int& TemplateLength, short responses[MaxShiftDistance][MaxTemplateLength],
	CPoint response_locations[MaxShiftDistance][MaxTemplateLength]);

extern bool IndividualMedianResponse(int prev_shift, CPoint& point,
	int& TemplateLength, short responses[MaxShiftDistance][MaxTemplateLength],
	CPoint response_locations[MaxShiftDistance][MaxTemplateLength]);

extern CVessels gTheVessels;

extern CPoint** gapSortedArrayOfSeedPoints;
extern int giNumOfSeedPoints;

// the sum of vessel widths at valid seed points
extern float gfWidthSum;
// the number of points participating with the sum
extern int giNumOfWidthSumMembers;

extern CIntersectionPoints gIntersectionPoints;

extern int giNumOfSomas;

extern int giReturnFlag ;

extern float gfContrast;

extern int giHitSomaFlag;

extern int giSomaVolume;
extern CPoint* gaSomaPoints;

extern int gi3DMedian;
extern float gf3DStdDev ;

extern int gaRGBColors[20][3];


// The following two lines were added for CANCER images
extern float gfHWidth;
extern float gfVWidth;

#endif
