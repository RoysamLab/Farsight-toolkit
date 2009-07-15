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

//////////////////////////////////////////////////////////////////////////
// FILE: track.cpp
//
// This file contains all tracking functions
//
//#pragma warning(disable:4786)  // disable STL-related warnings
//#pragma warning(disable:4710)  // compiler didn't inline function that got selected for inline expansion
//#pragma warning(disable:4514)  // The optimizer removed an inline function that is not called
//#pragma warning(disable:4702)  // unreachable STLport code

#include <iostream>
#include <fstream>
#include <cmath>
#include <queue>
#include <list>

#include "CONSTANTS.h"
#include "Mytypes.h"
#include "Config.h"
#include "Cpoints.h"
#include "Dllist.h"
#include "Cimage.h"
#include "C3dimage.h"
#include "Cvessel.h"
#include "Cvector.h"
#include "Template.h"
#include "Extern.h"
#include "StrEle.h"

using namespace std;

int giReversedFlag = 0;
extern ofstream gLogFile;
extern int giContrastThresholdMultiplier;

bool WithinTracedVessel(const CPoint& pSeedPoint);

// num of consecutive zero response of horizontal templates
int giHLViolationsCounter = 0;
int giHRViolationsCounter = 0;
int giVLViolationsCounter = 0;
int giVRViolationsCounter = 0;

//By Yousef (07-23-2207)
int giAnotherVesselViolation = 0;

// arrays to hold the location of template responses
CPoint gaHLeft[MaxShiftDistance][MaxTemplateLength];
CPoint gaHRight[MaxShiftDistance][MaxTemplateLength];
CPoint gaVLeft[MaxShiftDistance][MaxTemplateLength];
CPoint gaVRight[MaxShiftDistance][MaxTemplateLength];
// arrays to hold the values of template responses. Notice that the above
// arrays do contain such information, but these are much faster to reset
// using "memset", rather than two nested loops!!!
short gaiHRResp[MaxShiftDistance][MaxTemplateLength];
short gaiHLResp[MaxShiftDistance][MaxTemplateLength];
short gaiVRResp[MaxShiftDistance][MaxTemplateLength];
short gaiVLResp[MaxShiftDistance][MaxTemplateLength];
//
int MaxNumberOfResponses = MaxShiftDistance* MaxTemplateLength;

int Hdirections[32];
int Vdirections[32];

int iBestTemplateLength = 0;
int iBestTemplateLengthAtThisDirection = 0;

///////////////////////////////////////////////////////////////////////////////
// Function: Within Image Boundaries
// 
// Input:  a 3D point
// Output: bool Yes/NO
// Logic:
//    return yes if the given point is located within image boundaries
// 
//inline bool WithinImageBoundaries(const CPoint& aPoint)
//{
//	return (!(aPoint.m_iX <giMARGIN ||
//		aPoint.m_iX> iCols - giMARGIN ||
//		aPoint.m_iY <giMARGIN ||
//		aPoint.m_iY> iRows - giMARGIN ||
//		aPoint.m_iZ <giPadding ||
//		aPoint.m_iZ> iSlices - giPadding));
//}
//inline bool OutsideImageBoundaries(const CPoint& aPoint)
//{
//	return ((aPoint.m_iX <giMARGIN ||
//		aPoint.m_iX> iCols - giMARGIN ||
//		aPoint.m_iY <giMARGIN ||
//		aPoint.m_iY> iRows - giMARGIN ||
//		aPoint.m_iZ <giPadding ||
//		aPoint.m_iZ> iSlices - giPadding));
//}

bool ReversePoint(CPoint& aPoint)
{
	// amri: We just simply reverse the direction, and step forward
	aPoint.m_iHDir = static_cast<unsigned char>(ReverseTheta(aPoint.m_iHDir));

	// it turns out that VDir is relative to the HDir, so reversing is not as straight forward
	// we reverse the vertical direction so that it is relative to the hdir
	if (aPoint.m_iVDir > 0)
	{
		aPoint.m_iVDir = static_cast<unsigned char>(NumOfDirections - aPoint.m_iVDir);
	}

	aPoint.m_iX += gVectorsArray[aPoint.m_iHDir][aPoint.m_iVDir]->m_pIndices[3].m_iX;
	aPoint.m_iY += gVectorsArray[aPoint.m_iHDir][aPoint.m_iVDir]->m_pIndices[3].m_iY;
	aPoint.m_iZ += gVectorsArray[aPoint.m_iHDir][aPoint.m_iVDir]->m_pIndices[3].m_iZ;

	return true;
}

/////////////////////////////////////////////////////////////////////////
// Function: EndVessel
//
// Input: 
//   lPoint: left boundary point
//   rPoint: right boundary point
//   lResp : response of left template at "lPoint"
//   rResp : response of right template at "rPoint"
//
// Logic:
//   If the response of the left (right) template at lPoint (rPoint) in 
//   a direction perpendicular to that of lPoint.dir (rPoint.dir) is GREATER 
//   than half lResp (rResp), then these two point do not lie on a vessel's
//   boundary.
bool EndVessel(CPoint& aPoint, int, CPoint& HLPoint,
			   CPoint& HRPoint, CPoint& VLPoint, CPoint& VRPoint, CVessel&)
{
	// 1. any of the boundary points or the centerpoints is within another traced vessel
	//By Yousef: (07-23-2007)
	if (WithinTracedVessel(aPoint))
	{
		giAnotherVesselViolation++;
		//return true;
	}
	//count the number of concecutive errors
	if(giAnotherVesselViolation >1)
	{
		giAnotherVesselViolation = 0;
		return true;
	}
	///////////////////////////////

	int boundary_jumping = 0;
	if (The3DImage->WithinImagePadding(HRPoint, giMARGIN))
	{
		if (WithinTracedVessel(HRPoint))
		{
			boundary_jumping++;
		}
	}

	if (The3DImage->WithinImagePadding(HLPoint, giMARGIN))
	{
		if (WithinTracedVessel(HLPoint))
		{
			boundary_jumping++;
		}
	}

	if (The3DImage->WithinImagePadding(VRPoint, giMARGIN))
	{
		if (WithinTracedVessel(VRPoint))
		{
			boundary_jumping++;
		}
	}

	if (The3DImage->WithinImagePadding(VLPoint, giMARGIN))
	{
		if (WithinTracedVessel(VLPoint))
		{
			boundary_jumping++;
		}
	}
	if (boundary_jumping > 0)
		return true;

	// 2. If the response is below the threshold of 3 * contrast 
	//    Terminate after giParam_Trace_MaxNumOfConsecutiveZeroResponses
	int ContrastThresholdMultiplier = gConfig.GetContrastThresholdMultiplier();
	int Threshold = static_cast<int>(gfContrast *ContrastThresholdMultiplier);// * giUsedTemplateLength;

	int HLResp = HLPoint.m_iValue;
	int HRResp = HRPoint.m_iValue;
	int VLResp = VLPoint.m_iValue;
	int VRResp = VRPoint.m_iValue;

	// count the number of consecutive zero responses of the templates
	if (HLResp <= Threshold)
		giHLViolationsCounter++;
	else
		giHLViolationsCounter = 0;

	if (HRResp < Threshold)
		giHRViolationsCounter++;
	else
		giHRViolationsCounter = 0;

	if (VLResp <= Threshold)
		giVLViolationsCounter++;
	else
		giVLViolationsCounter = 0;

	if (VRResp < Threshold)
		giVRViolationsCounter++;
	else
		giVRViolationsCounter = 0;

	int iMaxAllowedStoppingViolations = gConfig.GetMaximumAllowedStoppingViolations();
	if (giHRViolationsCounter > iMaxAllowedStoppingViolations ||
		giHLViolationsCounter > iMaxAllowedStoppingViolations ||
		giVRViolationsCounter > iMaxAllowedStoppingViolations ||
		giVLViolationsCounter > iMaxAllowedStoppingViolations)
	{
		return true;
	}

	//// valid edges must have left and right boundaries
	//// allow only one template with low resonse
	//if (HLResp + HRResp + VLResp + VRResp < Threshold * 4)
	//{
	//	//	cout <<"\nEndPoint::SumLessThanThreshold" << endl;
	//	return true;
	//}

	return false;
}

void Display2DResponseArray(short aiResponses[][MaxTemplateLength])
{
	ofstream out("2DArray.txt");
	for (int j = 0; j < MaxShiftDistance; j++)
	{
		for (int i = 0; i < MaxTemplateLength; i++)
		{
			out << aiResponses[j][i] << "\t";
		}
		out << endl;
	}
}


int GetNonParallelResponse(const CPoint& aPoint, CPoint& HR, CPoint& HL,
						   CPoint& VR, CPoint& VL)
{
	register int i, j;
	register int Hdir, Vdir;

	int iDirectionFreedom = gConfig.GetDirectionalDegreeOfFreedom();

	PrepareDirectionsArray(Hdirections, iDirectionFreedom, aPoint.m_iHDir);
	PrepareDirectionsArray(Vdirections, iDirectionFreedom, aPoint.m_iVDir);
	CPoint HR_curr;
	CPoint HL_curr;
	CPoint VR_curr;
	CPoint VL_curr;

	// reset the best template responses at the current point
	int iBestResponse = 0;	

	int HR_length, HL_length, VR_length, VL_length;
	int HR_maxresponse, HL_maxresponse, VR_maxresponse, VL_maxresponse;
	HR_maxresponse = HL_maxresponse = VR_maxresponse = VL_maxresponse = 0;
	int HR_bestlength = 0, HL_bestlength = 0, VR_bestlength = 0, VL_bestlength = 0;


	// for each of the specified directions calculate the Hleft, HRight, Vleft
	// and Vright template responses.
	for (i = 0; i < iDirectionFreedom; i++)
	{
		Hdir = Hdirections[i];
		for (j = 0; j < iDirectionFreedom; j++)
		{
			// initialize the response arrays
			memset(gaiHRResp, 0, sizeof(short) * MaxNumberOfResponses);
			memset(gaiHLResp, 0, sizeof(short) * MaxNumberOfResponses);
			memset(gaiVRResp, 0, sizeof(short) * MaxNumberOfResponses);
			memset(gaiVLResp, 0, sizeof(short) * MaxNumberOfResponses);

			giHitSomaFlag = 0;

			Vdir = Vdirections[j];

			gHRightTemplatesArray[Hdir][Vdir]->
				CalculateKernelResponses(aPoint,
				Round(aPoint.m_fHWidth),
				gaHRight,
				gaiHRResp);

			if (giHitSomaFlag)
				return false;

			gHLeftTemplatesArray[Hdir][Vdir]->
				CalculateKernelResponses(aPoint,
				Round(aPoint.m_fHWidth),
				gaHLeft,
				gaiHLResp);

			//	Display2DResponseArray(gaiHLResp);

			if (giHitSomaFlag)
				return false;

			gVRightTemplatesArray[Hdir][Vdir]->
				CalculateKernelResponses(aPoint,
				Round(aPoint.m_fVWidth),
				gaVRight,
				gaiVRResp);

			if (giHitSomaFlag)
				return false;

			gVLeftTemplatesArray[Hdir][Vdir]->
				CalculateKernelResponses(aPoint,
				Round(aPoint.m_fVWidth),
				gaVLeft,
				gaiVLResp);

			if (giHitSomaFlag)
				return false;

			// find the best combination of responses from the templates lengths at
			// the current direction. Notice that the function must also select the
			// best template length.
			if (gConfig.GetMedianTracing())
			{
				if (IndividualMedianResponse(Round(aPoint.m_fHWidth),
					HR_curr,
					HR_length,
					gaiHRResp,
					gaHRight))
				{
					if (HR_curr.m_iValue > HR_maxresponse)
					{
						HR_maxresponse = HR_curr.m_iValue;
						HR = HR_curr;
						HR_bestlength = HR_length;
					}
				}
				if (IndividualMedianResponse(Round(aPoint.m_fHWidth),
					HL_curr,
					HL_length,
					gaiHLResp,
					gaHLeft))
				{
					if (HL_curr.m_iValue > HL_maxresponse)
					{
						HL_maxresponse = HL_curr.m_iValue;
						HL = HL_curr;
						HL_bestlength = HL_length;
					}
				}
				if (IndividualMedianResponse(Round(aPoint.m_fVWidth),
					VR_curr,
					VR_length,
					gaiVRResp,
					gaVRight))
				{
					if (VR_curr.m_iValue > VR_maxresponse)
					{
						VR_maxresponse = VR_curr.m_iValue;
						VR = VR_curr;
						VR_bestlength = VR_length;
					}
				}
				if (IndividualMedianResponse(Round(aPoint.m_fVWidth),
					VL_curr,
					VL_length,
					gaiVLResp,
					gaVLeft))
				{
					if (VL_curr.m_iValue > VL_maxresponse)
					{
						VL_maxresponse = VL_curr.m_iValue;
						VL = VL_curr;
						VL_bestlength = VL_length;
					}
				}
			}
			else
			{
				if (IndividualAverageResponse(Round(aPoint.m_fHWidth),
					HR_curr,
					HR_length,
					gaiHRResp,
					gaHRight))
				{
					if (HR_curr.m_iValue > HR_maxresponse)
					{
						HR_maxresponse = HR_curr.m_iValue;
						HR = HR_curr;
						HR_bestlength = HR_length;
					}
				}
				if (IndividualAverageResponse(Round(aPoint.m_fHWidth),
					HL_curr,
					HL_length,
					gaiHLResp,
					gaHLeft))
				{
					if (HL_curr.m_iValue > HL_maxresponse)
					{
						HL_maxresponse = HL_curr.m_iValue;
						HL = HL_curr;
						HL_bestlength = HL_length;
					}
				}
				if (IndividualAverageResponse(Round(aPoint.m_fVWidth),
					VR_curr,
					VR_length,
					gaiVRResp,
					gaVRight))
				{
					if (VR_curr.m_iValue > VR_maxresponse)
					{
						VR_maxresponse = VR_curr.m_iValue;
						VR = VR_curr;
						VR_bestlength = VR_length;
					}
				}
				if (IndividualAverageResponse(Round(aPoint.m_fVWidth),
					VL_curr,
					VL_length,
					gaiVLResp,
					gaVLeft))
				{
					if (VL_curr.m_iValue > VL_maxresponse)
					{
						VL_maxresponse = VL_curr.m_iValue;
						VL = VL_curr;
						VL_bestlength = VL_length;
					}
				}
			}
		}  // end loop over VDirections	
	}	// end loop over HDirections	
	list<int> lengths;
	if (HR_bestlength > 0)
		lengths.push_back(HR_bestlength);
	if (HL_bestlength > 0)
		lengths.push_back(HL_bestlength);
	if (VR_bestlength > 0)
		lengths.push_back(VR_bestlength);
	if (VL_bestlength > 0)
		lengths.push_back(VL_bestlength);

	lengths.sort();

	if (lengths.front() > 0)
	{
		iBestTemplateLength = lengths.front();
	}
	else
	{
		return false;
	}
	iBestResponse = HR_maxresponse + HL_maxresponse + VR_maxresponse + VL_maxresponse;
	return iBestResponse;
}


int GetParallelResponse(const CPoint& aPoint, CPoint& HR, CPoint& HL,
						CPoint& VR, CPoint& VL)
{
	int iDirectionFreedom = gConfig.GetDirectionalDegreeOfFreedom();
	CPoint HR_curr;
	CPoint HL_curr;
	CPoint VR_curr;
	CPoint VL_curr;
	register int i, j;
	register int Hdir, Vdir;
	PrepareDirectionsArray(Hdirections, iDirectionFreedom, aPoint.m_iHDir);
	PrepareDirectionsArray(Vdirections, iDirectionFreedom, aPoint.m_iVDir);

	// reset the best template responses at the current point
	int iBestResponse = 0;	

	// for each of the specified directions calculate the Hleft, HRight, Vleft
	// and Vright template responses.
	int iBestResponseAtThisDirection;
	for (i = 0; i < iDirectionFreedom; i++)
	{
		Hdir = Hdirections[i];
		for (j = 0; j < iDirectionFreedom; j++)
		{
			// initialize the response arrays
			memset(gaiHRResp, 0, sizeof(short) * MaxNumberOfResponses);
			memset(gaiHLResp, 0, sizeof(short) * MaxNumberOfResponses);
			memset(gaiVRResp, 0, sizeof(short) * MaxNumberOfResponses);
			memset(gaiVLResp, 0, sizeof(short) * MaxNumberOfResponses);
			iBestResponseAtThisDirection = 0;

			giHitSomaFlag = 0;

			Vdir = Vdirections[j];

			gHRightTemplatesArray[Hdir][Vdir]->CalculateKernelResponses(aPoint,Round(aPoint.m_fHWidth),gaHRight,gaiHRResp);

			if (giHitSomaFlag)
				return false;

			gHLeftTemplatesArray[Hdir][Vdir]->CalculateKernelResponses(aPoint,Round(aPoint.m_fHWidth),gaHLeft,gaiHLResp);

			if (giHitSomaFlag)
				return false;

			gVRightTemplatesArray[Hdir][Vdir]->CalculateKernelResponses(aPoint,Round(aPoint.m_fVWidth),gaVRight,gaiVRResp);

			if (giHitSomaFlag)
				return false;

			gVLeftTemplatesArray[Hdir][Vdir]->CalculateKernelResponses(aPoint,Round(aPoint.m_fVWidth),gaVLeft,gaiVLResp);

			if (giHitSomaFlag)
				return false;

			// find the best combination of responses from the templates lengths at
			// the current direction. Notice that the function must also select the
			// best template length.

			if (gConfig.GetMedianTracing())
			{
				iBestResponseAtThisDirection = GetMedianResponse(aPoint,HR_curr,HL_curr,VR_curr,VL_curr,
					iBestTemplateLengthAtThisDirection);
			}
			else
			{
				iBestResponseAtThisDirection = GetAverageResponse(aPoint,HR_curr,HL_curr,VR_curr,VL_curr,
					iBestTemplateLengthAtThisDirection);
			}

			if (iBestResponseAtThisDirection > iBestResponse)
			{
				iBestResponse = iBestResponseAtThisDirection;
				iBestTemplateLength = iBestTemplateLengthAtThisDirection;	
				HR = HR_curr;
				HL = HL_curr;
				VR = VR_curr;
				VL = VL_curr;
			}
		}  // end loop over VDirections	
	} // end of loop over HDirections
	return iBestResponse;
}

bool TraceAPoint(CPoint& center, CVessel& TheVessel,
				 TopEndMiddleNeither dirFlag)
{
	CPoint HR;
	CPoint HL;
	CPoint VR;
	CPoint VL;
	CPoint original_center = center;
	int iBestResponse;

	if (cfg_trace_nonparallel)
	{
		if ((iBestResponse = GetNonParallelResponse(center, HR, HL, VR, VL)) == 0)
			return false;
	}
	else
	{
		if ((iBestResponse = GetParallelResponse(center, HR, HL, VR, VL)) == 0)
			return false;
	}
	center.m_iValue = iBestResponse;

	// check for boundary points out of image
	if( !The3DImage->WithinImagePadding(center, giMARGIN)
		|| !The3DImage->WithinImagePadding(HR, giMARGIN)
		|| !The3DImage->WithinImagePadding(HL, giMARGIN))
		return false;

	int iStepSize = (int) (iBestTemplateLength / 4.0 + 0.5) ;

	int iMaxStepSize = gConfig.GetMaximumStepSize();
	if (iStepSize < 3)
		iStepSize = 3;
	if (iStepSize > iMaxStepSize)
		iStepSize = iMaxStepSize;


	// 3.1 Update the location and direction of the current center point.
	//     - The center point is located at midpoint between the left and right
	//  	 boundary points.
	//     - The center point follows the direction of the stronger edge.
	GetCenterLocation(center, HR, HL, VR, VL);

	//Yousef
	if(center.m_fHWidth > 7)
		return false;

	if (! The3DImage->WithinImagePadding(center, giMARGIN))
		return false;

	// Test if we reached the end of a vessel
	if (EndVessel(center, iBestResponse, HL, HR, VL, VR, TheVessel))
		return false;	

	GetCenterDirection(center, HR, HL, VR, VL);

	HRPoints.push_back(HR);
	HLPoints.push_back(HL);
	CenterPoints.push_back(center);

	// The following lines of code were added to accumulate 
	// stats for CANCER images

	float h_width = static_cast<float>(HL.FindDistance(&HR));
	float v_width = static_cast<float>(VL.FindDistance(&VR));
	
	int MinAllowableWidth = gConfig.GetMinimumShiftDistance();
	int MaxAllowableWidth = Round(static_cast<double>(gConfig.GetMaximumShiftDistance()) * 2.0 * sqrt(2.0));
	if(h_width < MinAllowableWidth || h_width > MaxAllowableWidth
		|| v_width < MinAllowableWidth || v_width > MaxAllowableWidth)
		return false;

	TheVessel.m_iNumOfPoints++;

	if (!TheVessel.ExtendVessel(&center, & HL, & HR, & VL, & VR, dirFlag))
		return false;

	// Estimate the location of the next center point.
	CPoint pPoint = gVectorsArray[center.m_iHDir][center.m_iVDir]->m_pIndices[iStepSize];

	// From Here 7-26-99
	center.m_iX += pPoint.m_iX;
	center.m_iY += pPoint.m_iY;
	center.m_iZ += pPoint.m_iZ;
	center.m_fHWidth = h_width;
	center.m_fVWidth = v_width;

	// to avoid getting into an infinite loop
	if (center == original_center)
		return false;

	giReversedFlag = 0;

	if (The3DImage->WithinImagePadding(center, giMARGIN))
		return true;
	else
		return false;
}

/////////////////////////////////////////////////////////////////////////
// Function: InitializeGlobalTrackingVariables
void InitializeGlobalTrackingVariables(C3DImage& anImage)
{
	//	register int i;
	//char* pchStr = NULL;

	//pchStr = gConfig.GetStringValue("Param.Trace.MinSegmentLength");
	//if (pchStr != NULL)
	//	giParam_Trace_MinSegmentLength = atoi(pchStr);
	//else
	//{
	//	cout << "The tracing parameter MinSegmentLength was not found\n"
	//		<< "Default value of 10 is used" << endl;
	//	giParam_Trace_MinSegmentLength = 10;
	//}

	// experimental
	//	gfContrast = giParam_Trace_MaxContrastValue;	

	int iMaxContrastValue = 15;
	if (gfContrast < 0.0)
	{
		gfContrast = static_cast<float>((gf3DStdDev + gfStdDev) / 2.0);// - 0.5;

		// if the foreground median is less than the median of the projection
		// then we have too much noise ...
		if (giForegroundMedian < 0.5 * giMedian)
			gfContrast += 1;//gf3DStdDev;	
		if (gfContrast < 0.5)
			gfContrast = 0.5;
		if (gfContrast > iMaxContrastValue)
			gfContrast = static_cast<float>(iMaxContrastValue);
	}

	gLogFile << "Contrast Used: " << gfContrast << endl;

	// to track an image, we need to associate the templates with that image
	register int Hdir, Vdir;
	for (Hdir = 0; Hdir < NumOfDirections; Hdir++)
	{
		for (Vdir = 0; Vdir < NumOfDirections; Vdir++)
		{
			gHLeftTemplatesArray[Hdir][Vdir]->AssociateWithImage(anImage);
			gHRightTemplatesArray[Hdir][Vdir]->AssociateWithImage(anImage);
			gVLeftTemplatesArray[Hdir][Vdir]->AssociateWithImage(anImage);
			gVRightTemplatesArray[Hdir][Vdir]->AssociateWithImage(anImage);
		}
	}

	giHRViolationsCounter = 0;
	giHLViolationsCounter = 0;
	giVRViolationsCounter = 0;
	giVLViolationsCounter = 0;
}

bool WithinTracedVessel(const CPoint& SeedPoint)
{
	return(TracedImage[SeedPoint.m_iZ][SeedPoint.m_iY][SeedPoint.m_iX]);
}

///////////////////////////////////////////////////////////////////////////////////
// Function: TraceA3DImage
// 
// Input:
//  	anImage  : A 3D  image to be tracked
// Logic:
//		For all valid seed points, track all vessels reached from that seed point
//		Upon encountering intersection/branching points, follow the stronger branch.
void TraceA3DImage(C3DImage& anImage)
{
	// Working variables
	CPoint seed;
	CPoint forward_point;
	CPoint reverse_point;
	CPoint pre_reverse_point;
	CVessel* aVessel = NULL;
	int vesselID = 1;
	int i;

	// initialize global variables
	InitializeGlobalTrackingVariables(anImage);
	int count = 0;

	list<deque<CPoint> >::reverse_iterator j;

	for (j = Seeds.rbegin(); j != Seeds.rend(); j++)
	{
		cout << "\r\tTracing seed " << ++count << " out of " << Seeds.size()
			<< " ... " << flush;

		// if the seed point or any of its boundaries is within another traced vessel, ignore it
		bool to_continue = false;
		for (i = 0; i < 5; i++)
		{
			if(!The3DImage->WithinImagePadding((*j)[i],giMARGIN) || WithinTracedVessel(((*j)[i])))
			{
				to_continue = true;
				break;
			}
		}
		if(to_continue)
			continue;

		seed = ((*j)[0]);		

		// amri 6-16-05 BUGFIX: ignore seeds with zero width
		if( Round( seed.m_fHWidth ) == 0 || Round( seed.m_fVWidth ) == 0 )
			continue;

		forward_point = seed;
	
		if (WithinTracedVessel(forward_point))
			continue;	

		// we'll go in the reverse direction relative to the seed point
		reverse_point = forward_point;

		TracedPoints.push_back(forward_point);

		// create a new vessel and assign it an ID
		aVessel = new CVessel(vesselID);

		giReversedFlag = 0;

		// trace the vessel in forward direction
		while (TraceAPoint(forward_point, * aVessel, OnTop));
		//By Yousef (07-23-2007)
		giAnotherVesselViolation = 0;
		///////////////////////////////

		// get the tail of the traced vessel if it exists
		if (aVessel->m_Center.tail && aVessel->m_Center.tail->data)
		{
			reverse_point = *aVessel->m_Center.tail->data;
		}

		// reverse the direction
		if (ReversePoint(reverse_point))
		{
			if (The3DImage->WithinImagePadding(reverse_point,giMARGIN) &&
				!WithinTracedVessel(reverse_point))
			{
				// trace in the reverse direction
				pre_reverse_point = reverse_point;
				// check if the refined point is too far from the tail of the vessel 
				if (pre_reverse_point.FindDistance(&reverse_point) < 6)
				{
					if (!WithinTracedVessel(reverse_point))
						while (TraceAPoint(reverse_point, * aVessel, OnEnd));
					//By Yousef (07-23-2007)
					giAnotherVesselViolation = 0;
					///////////////////////////////
				}
			}
		}


		// revise the widths using gaussian smooth
		if (aVessel->m_iLength > 0)
		{
			aVessel->ReviseWidths();

			// ignore vessels shorter than 3 times its width
			//				if ((aVessel->m_iLength >= 2 * Round(aVessel->m_fHWidth) &&
			//					aVessel->m_iLength >= 2 * Round(aVessel->m_fVWidth)) &&
			//					!((aVessel->m_fHWidth == 1) && (aVessel->m_fVWidth == 1)))
			if ((aVessel->m_iLength >= 1.0 * aVessel->m_fHWidth &&
				aVessel->m_iLength >= 1.0 * aVessel->m_fVWidth) &&
				!((aVessel->m_fHWidth == 1) && (aVessel->m_fVWidth == 1)))

				//&& aVessel->m_iLength >= giParam_Trace_MinSegmentLength)
			{
				gTheVessels.AddVessel(aVessel);
				vesselID++;
				

				// to prevent re-tracking the same vessel again, we mark it in
				// another image using spheres at each center point. 

				aVessel->MarkVessel();
				TracedSeeds.push_back(seed);
			}
			else
				delete aVessel;
		}
		else
			delete aVessel;
	}

	cout << "Done" << endl;
}

