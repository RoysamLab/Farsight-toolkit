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

#ifndef global_h
#define global_h

#include <deque>
#include "Branches.h"

list<deque<CPoint> > Seeds;

/////////////////
// operation parameters
// if set use fixed length templates and stricter stopping criteria
bool*** TracedImage;
int giParam_Trace_MinSegmentLength = 10;

// configuration items and their corresponding default values
char cfg_output_image_format[4] = "tif";
int cfg_output_image_bgcolor = 0;
int cfg_output_image_fgcolor = 1;
int cfg_output_projections = 1;
int cfg_output_soma_draw_trees = 0;
int cfg_mode_debug = 1;
int cfg_trace_nonparallel = 0;
bool cfg_verify_amri = false;
bool cfg_seedsfromfile = false;
std::string cfg_seedsfilename;

// a copy of projection images for drawing purposes
CImage* CanvasXY;
CImage* CanvasXZ;
CImage* CanvasYZ;

int giFirstSlice;
int giLastSlice;
int giSliceNumber;

// a margine from each side of the image. This is necessary because when
// considering a center point, we need space to shift the template.
int giMARGIN = 3;

int giMinSeedPixelValue;
int giMaxExtensionLength;
// the maximum shift away from the center. this value is more than enough
// since the images we are considering contains no structures wider than
// 6 pixels in diameter

int giSmallestResponseAccepted;

// if at any point the orthogonal response exceeds this value times the 
// point's response, then stop tracking. This parameter is used by the 
// EndVessel function in the file track.cpp
float gfOrthogonalResponsePercentage;

// a configuration object
CConfig gConfig;

list<CPoint> CenterPoints;
list<CPoint> HRPoints;
list<CPoint> HLPoints;
list<CPoint> VRPoints;
list<CPoint> VLPoints;
list<CPoint> UnverifiedSeeds;
list<CPoint> VerifiedSeedsCenter;
list<CPoint> VerifiedSeedsRight;
list<CPoint> VerifiedSeedsLeft;
list<CPoint> AmriVerifiedSeedsCenter;
list<CPoint> TracedPoints;
list<CPoint> TracedSeeds;

///////////////////////////////////
// Image constants and variables //
///////////////////////////////////

int* giMedians;
float* gfStdDevs;

int giMedian = - 1;
float gfStdDev = - 1.0;
int giForegroundMedian = - 1;
int giBackgroundMedian = - 1;
float gfForegroundStdDev = - 1.0;
float gfBackgroundStdDev = - 1.0;

int giXZMedian = - 1;
int giYZMedian = - 1;
float gfXZStdDev = - 1.0;
float gfYZStdDev = - 1.0;

// the image that is being traked
CImage* TheTargetImage;
CImage* EdgeImage;
//	CImage *TrackImage;
CImage* ProjectionImage;
CImage* GlobalImage;
CImage* CenterlinesImage;
CImage* BoundariesImage;
CImage* BranchingPointsImage;

CImage* TrackImageXY;
CImage* TrackImageXZ;
CImage* TrackImageYZ;

// a 3D image that will hold the stack of images to be tracked. The tracking results
// will also be displayed on this image.
C3DImage* The3DImage;
//By Yousef
//This image will hold the 3D soma image. The image will be initialized to zeros, and will be filled with oned
//at the locations of soma
C3DImage* Soma3DImage;
//By YOusef
//Here is another image to hold the soma labels
C3DImage* SomaLabelsImage;
//By Yousef
//Use this object to detect branch points
CBranches* BranchPoints;

// an empty 3D image that is used to register the tracked vessels and seed points.
// this image is needed to avoid tracking a vessel more than one time.
C3DImage* Empty3DImage;

// file name of the input image, its type, and the path to it
//	char gachBasicFileName[256];
char gachFileOfSlices[256];
char gachPath[128];
char gachType[10];

long glSumOfSeedPointValues = 0;
long glNumOfSeedPoints = 0;
int giContrastThreshold = 0;

int gaiForegroundHistogram[256];
int gaiBackgroundHistogram[256];


///////////
// the following two variables are used to make sure that we do
// not overstep the 3D-image boundaries
unsigned char* gpuchThe3DImageData = 0;
long glThe3DImageLength = 0;

//////////////////////////////////////
// Template constants and variables //
//////////////////////////////////////

// define 4-2D arrays of templates. left/right-horizontal and left/right-vertical
// The direction of each template is defined by the angles it makes with the
// positive x-axis in the xy plane (Hdir) and the angle it makes with the positive
// z-axis (Vdir)
CTemplate* gHLeftTemplatesArray[NumOfDirections][NumOfDirections];
CTemplate* gVLeftTemplatesArray[NumOfDirections][NumOfDirections];
CTemplate* gHRightTemplatesArray[NumOfDirections][NumOfDirections];
CTemplate* gVRightTemplatesArray[NumOfDirections][NumOfDirections];

//////////////////////////////////////
//  Vector constants and variables  //
//////////////////////////////////////

// an array of shift vectors, one for each direction
CVector* gVectorsArray[NumOfDirections][NumOfDirections];

////////////////////////////////////////
//  Tracking constants and variables  //
////////////////////////////////////////

//	CQueue<CPoint> *InitialPoints = NULL;   		  // an array of queues of points

// a 2D array of pointers to seed points
CPoint*** gapArrayOfSeedPoints = NULL; 

// a 2D array of link lists of filling points
//CDLList<CFillingPoint> ***ArrayOfFillings = NULL;


// the sum of vessel widths at valid seed points
float gfWidthSum = 0;
// the number of points participating with the sum
int giNumOfWidthSumMembers = 0;

CIntersectionPoints gIntersectionPoints;

CVessels gTheVessels;

int giNumOfSomas = 0;

// check this later
int giReturnFlag = 0;

int gi3DMedian = 0;
float gf3DStdDev = 0.0;

float gfContrast = - 1;


//////////////////////////////////////
//  	   Neurolucida Stuff		//
//////////////////////////////////////

char gachSpaces[1000];
char gachLeadingSpaces[1000];
int giNumOfLeadingSpaces = 0;

// the colors for drawing trees in Neurolucida format
int gaRGBColors[20][3] =
{
	{192, 192, 0}, {0, 192, 0}, {255, 0, 0}, {0, 0, 255}, {0, 0, 255},
	{5, 5, 5}, {250, 252, 100}, {0, 255, 255}, {255, 192, 255},
	{128, 128, 255}, {255, 192, 192}, {255, 128, 128}, {255, 0, 255},
	{253, 161, 248}, {192, 255, 255}, {251, 146, 43}, {87, 255, 190},
	{188, 255, 129}, {249, 252, 12}, {251, 9, 140}
};

unsigned short heated_object_colormap[] =
{
	0, 0, 0, 35, 0, 0, 52, 0, 0, 60, 0, 0, 63, 1, 0, 64, 2, 0, 68, 5, 0, 69,
	6, 0, 72, 8, 0, 74, 10, 0, 77, 12, 0, 78, 14, 0, 81, 16, 0, 83, 17, 0, 85,
	19, 0, 86, 20, 0, 89, 22, 0, 91, 24, 0, 92, 25, 0, 94, 26, 0, 95, 28, 0,
	98, 30, 0, 100, 31, 0, 102, 33, 0, 103, 34, 0, 105, 35, 0, 106, 36, 0,
	108, 38, 0, 109, 39, 0, 111, 40, 0, 112, 42, 0, 114, 43, 0, 115, 44, 0,
	117, 45, 0, 119, 47, 0, 119, 47, 0, 120, 48, 0, 122, 49, 0, 123, 51, 0,
	125, 52, 0, 125, 52, 0, 126, 53, 0, 128, 54, 0, 129, 56, 0, 129, 56, 0,
	131, 57, 0, 132, 58, 0, 134, 59, 0, 134, 59, 0, 136, 61, 0, 137, 62, 0,
	137, 62, 0, 139, 63, 0, 139, 63, 0, 140, 65, 0, 142, 66, 0, 142, 66, 0,
	143, 67, 0, 143, 67, 0, 145, 68, 0, 145, 68, 0, 146, 70, 0, 146, 70, 0,
	148, 71, 0, 148, 71, 0, 149, 72, 0, 149, 72, 0, 151, 73, 0, 151, 73, 0,
	153, 75, 0, 153, 75, 0, 154, 76, 0, 154, 76, 0, 154, 76, 0, 156, 77, 0,
	156, 77, 0, 157, 79, 0, 157, 79, 0, 159, 80, 0, 159, 80, 0, 159, 80, 0,
	160, 81, 0, 160, 81, 0, 162, 82, 0, 162, 82, 0, 163, 84, 0, 163, 84, 0,
	165, 85, 0, 165, 85, 0, 166, 86, 0, 166, 86, 0, 166, 86, 0, 168, 87, 0,
	168, 87, 0, 170, 89, 0, 170, 89, 0, 171, 90, 0, 171, 90, 0, 173, 91, 0,
	173, 91, 0, 174, 93, 0, 174, 93, 0, 176, 94, 0, 176, 94, 0, 177, 95, 0,
	177, 95, 0, 179, 96, 0, 179, 96, 0, 180, 98, 0, 182, 99, 0, 182, 99, 0,
	183, 100, 0, 183, 100, 0, 185, 102, 0, 185, 102, 0, 187, 103, 0, 187, 103,
	0, 188, 104, 0, 188, 104, 0, 190, 105, 0, 191, 107, 0, 191, 107, 0, 193,
	108, 0, 193, 108, 0, 194, 109, 0, 196, 110, 0, 196, 110, 0, 197, 112, 0,
	197, 112, 0, 199, 113, 0, 200, 114, 0, 200, 114, 0, 202, 116, 0, 202, 116,
	0, 204, 117, 0, 205, 118, 0, 205, 118, 0, 207, 119, 0, 208, 121, 0, 208,
	121, 0, 210, 122, 0, 211, 123, 0, 211, 123, 0, 213, 124, 0, 214, 126, 0,
	214, 126, 0, 216, 127, 0, 217, 128, 0, 217, 128, 0, 219, 130, 0, 221, 131,
	0, 221, 131, 0, 222, 132, 0, 224, 133, 0, 224, 133, 0, 225, 135, 0, 227,
	136, 0, 227, 136, 0, 228, 137, 0, 230, 138, 0, 230, 138, 0, 231, 140, 0,
	233, 141, 0, 233, 141, 0, 234, 142, 0, 236, 144, 0, 236, 144, 0, 238, 145,
	0, 239, 146, 0, 241, 147, 0, 241, 147, 0, 242, 149, 0, 244, 150, 0, 244,
	150, 0, 245, 151, 0, 247, 153, 0, 247, 153, 0, 248, 154, 0, 250, 155, 0,
	251, 156, 0, 251, 156, 0, 253, 158, 0, 255, 159, 0, 255, 159, 0, 255, 160,
	0, 255, 161, 0, 255, 163, 0, 255, 163, 0, 255, 164, 0, 255, 165, 0, 255,
	167, 0, 255, 167, 0, 255, 168, 0, 255, 169, 0, 255, 169, 0, 255, 170, 0,
	255, 172, 0, 255, 173, 0, 255, 173, 0, 255, 174, 0, 255, 175, 0, 255, 177,
	0, 255, 178, 0, 255, 179, 0, 255, 181, 0, 255, 181, 0, 255, 182, 0, 255,
	183, 0, 255, 184, 0, 255, 187, 7, 255, 188, 10, 255, 189, 14, 255, 191,
	18, 255, 192, 21, 255, 193, 25, 255, 195, 29, 255, 197, 36, 255, 198, 40,
	255, 200, 43, 255, 202, 51, 255, 204, 54, 255, 206, 61, 255, 207, 65, 255,
	210, 72, 255, 211, 76, 255, 214, 83, 255, 216, 91, 255, 219, 98, 255, 221,
	105, 255, 223, 109, 255, 225, 116, 255, 228, 123, 255, 232, 134, 255, 234,
	142, 255, 237, 149, 255, 239, 156, 255, 240, 160, 255, 243, 167, 255, 246,
	174, 255, 248, 182, 255, 249, 185, 255, 252, 193, 255, 253, 196, 255, 255,
	204, 255, 255, 207, 255, 255, 211, 255, 255, 218, 255, 255, 222, 255, 255,
	225, 255, 255, 229, 255, 255, 233, 255, 255, 236, 255, 255, 240, 255, 255,
	244, 255, 255, 247, 255, 255, 255
};

#endif
