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
#include <list>
#include <queue>
#include <string>
#include <iomanip>

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

using namespace std;

extern int giReversedFlag;
extern ofstream gLogFile;

// Variables declared in track.cpp
extern int giParam_Trace_MaxNumOfConsecutiveZeroResponses;
extern int giParam_Trace_MaxNumOfConsecutiveViolations;
extern int giParam_Trace_MaxContrastValue;
extern int giParam_Trace_MinSegmentLength;
extern CPoint gaHLeft[MaxShiftDistance][MaxTemplateLength];
extern CPoint gaHRight[MaxShiftDistance][MaxTemplateLength];
extern CPoint gaVLeft[MaxShiftDistance][MaxTemplateLength];
extern CPoint gaVRight[MaxShiftDistance][MaxTemplateLength];
extern short gaiHRResp[MaxShiftDistance][MaxTemplateLength];
extern short gaiHLResp[MaxShiftDistance][MaxTemplateLength];
extern short gaiVRResp[MaxShiftDistance][MaxTemplateLength];
extern short gaiVLResp[MaxShiftDistance][MaxTemplateLength];

extern short median__(short a, short b, short c);

void MedianFilterResponses(short responses[MaxTemplateLength])
{
	int iToLength = gConfig.GetMaximumTemplateLength();

	short smoothed_responses[MaxTemplateLength];
	memset(smoothed_responses, 0, sizeof(short) * MaxTemplateLength);

	for (int i = 0; i <= iToLength; i++)
	{
		if (i != 0 && i != iToLength)
		{
			smoothed_responses[i] = median__(responses[i], 
										responses[i + 1],
										responses[i - 1]);
		}
		else
		{
			smoothed_responses[i] = responses[i];
		}
	}
	memcpy(responses, smoothed_responses, sizeof(short) * MaxTemplateLength);
}

void GaussianFilterResponses(short responses[MaxTemplateLength])
{
	short smoothed_responses[MaxTemplateLength];
	memset(smoothed_responses, 0, sizeof(short) * MaxTemplateLength);

	int iToLength = gConfig.GetMaximumTemplateLength();

	for (int i = 0; i <= iToLength; i++)
	{
		if (i == 0)
		{
			smoothed_responses[i] = static_cast<short>(Round(static_cast<float> (responses[i] + responses[i + 1]) / 2.0));
		}
		else
		{
			if (i == iToLength)
			{
				smoothed_responses[i] = static_cast<short>(Round(static_cast<float>(responses[i] + responses[i - 1]) / 2.0));
			}
			else
			{
				smoothed_responses[i] = static_cast<short>(Round(0.5 * static_cast<float>(responses[i]) + 0.25 * static_cast<float> (responses[i + 1]) +
											0.25 * static_cast<float> (responses[i - 1])));
			}
		}
	}
	memcpy(responses, smoothed_responses, sizeof(short) * MaxTemplateLength);
}

bool IndividualAverageResponse(int prev_shift, CPoint& point,
	int& TemplateLength, short responses[MaxShiftDistance][MaxTemplateLength],
	CPoint response_locations[MaxShiftDistance][MaxTemplateLength])
{
	int iMaxShiftDistance = gConfig.GetMaximumShiftDistance();
	CPoint* best_point = NULL;

	// From each response array, select the boundary point resulting in the 
	// highest response density. 
	register unsigned char i, j;
	register float fRespDensity = 0;
	register float fMaxRespDensity = 0;
	register float fMaxRespDensityAtThisLength = 0;
	register float fMaxDensity = 0;

	register short sum[MaxShiftDistance];
	memset(sum, 0, sizeof(short) * MaxShiftDistance);

	int iRelativeShiftDistance = gConfig.GetRelativeShiftDistance();
	char min_shift, max_shift;
	int radius = Round((float) prev_shift / 2.0);

	int iFromLength = gConfig.GetMinimumTemplateLength();
	int iToLength = gConfig.GetMaximumTemplateLength();

	if (prev_shift == 0)
	{
		min_shift = 0;
		max_shift = static_cast<char>(iMaxShiftDistance);
	}
	else
	{
		min_shift = static_cast<char>(radius - iRelativeShiftDistance);
		max_shift = static_cast<char>(radius + iRelativeShiftDistance);
		if (min_shift < 0)
			min_shift = 0;
		if (min_shift > iMaxShiftDistance)
			min_shift = static_cast<char>(iMaxShiftDistance);
		if (max_shift > iMaxShiftDistance)
			max_shift = static_cast<char>(iMaxShiftDistance);
		if (min_shift > max_shift)
		{
			return false;
		}
	}

	for (i = min_shift; i <= max_shift; i++)
	{
		for (j = 0; j < iFromLength; j++)
		{
			sum[i] = static_cast<short>(sum[i] + responses[i][j]);
		}
	}

	float min_allow = 0.0;

	bool stop_shift = false;

	for (j = static_cast<unsigned char>(iFromLength); j <= iToLength; j++)
	{
		if (!stop_shift)
		{
			for (i = min_shift; i <= max_shift; i++)
			{
				min_allow = static_cast<float>(-0.8 * fMaxDensity);

				sum[i] = static_cast<short>(sum[i] + responses[i][j]);
				// is this length/shift the best density for the Horizontal Right Point?
				if (sum[i] > min_allow)
				{
					fRespDensity = sum[i] / (float) j;
					if (fRespDensity > fMaxDensity)
					{
						fMaxDensity = fRespDensity;
						if (best_point == NULL)
						{
							best_point = new CPoint(response_locations[i][j]);
						}
						else
							*best_point = response_locations[i][j];
						best_point->m_iValue = static_cast<int>(fRespDensity);
					}
				}
				else
				{
					if (fMaxDensity > 0)
					{
						stop_shift = true;
						break;
					}
				}
			}
		}

		// At this template length, this is the best density
		fMaxRespDensityAtThisLength = fMaxDensity;

		if (fMaxRespDensityAtThisLength >= fMaxRespDensity)
		{
			fMaxRespDensity = fMaxRespDensityAtThisLength;
			TemplateLength = j;
		}
	}

	if (best_point)
	{
		point = *best_point;
		return true;
	}
	else
		return false;
}

int GetAverageResponse(const CPoint& aPoint, CPoint& pHR, CPoint& pHL,
	CPoint& pVR, CPoint& pVL, int& TemplateLength)
{
	int HR_length, HL_length, VR_length, VL_length;
	if (!IndividualAverageResponse(Round(aPoint.m_fHWidth),
		 	pHR,
		 	HR_length,
		 	gaiHRResp,
		 	gaHRight) ||
		!IndividualAverageResponse(Round(aPoint.m_fHWidth),
						  	pHL,
						  	HL_length,
						  	gaiHLResp,
						  	gaHLeft) ||
		!IndividualAverageResponse(Round(aPoint.m_fVWidth),
										 	pVR,
										 	VR_length,
										 	gaiVRResp,
										 	gaVRight) ||
		!IndividualAverageResponse(Round(aPoint.m_fVWidth),
														  	pVL,
														  	VL_length,
														  	gaiVLResp,
														  	gaVLeft))
		return 0;

	list<int> lengths;
	if (HR_length > 0)
		lengths.push_back(HR_length);
	if (HL_length > 0)
		lengths.push_back(HL_length);
	if (VR_length > 0)
		lengths.push_back(VR_length);
	if (VL_length > 0)
		lengths.push_back(VL_length);

	lengths.sort();

	if (lengths.front() > 0)
	{
		TemplateLength = lengths.front();
	}
	else
	{
		return 0;
	}

	return (pHR.m_iValue + pHL.m_iValue + pVR.m_iValue + pVL.m_iValue);
}

#define rel_7_0 1

#ifdef rel_7_0
//bool IndividualMedianResponse(int prev_shift, CPoint& point,
//	int& TemplateLength, short responses[MaxShiftDistance][MaxTemplateLength],
//	CPoint response_locations[MaxShiftDistance][MaxTemplateLength])
//{
//	CPoint* best_point = NULL;
//
//	register unsigned char i, j;
//	register float fRespDensity = 0;
//	register float fMaxRespDensity = 0;
//	register float fMaxRespDensityAtThisLength = 0;
//	register float fMaxDensity = 0;
//
//	register short sum[MaxShiftDistance];
//	priority_queue<short> responses_queue[MaxShiftDistance];
//	memset(sum, 0, sizeof(short) * MaxShiftDistance);
//
//	int shift_freedom = cfg_trace_shift_freedom;
//	char min_shift, max_shift;
//	int radius = Round((float) prev_shift / 2.0);
//
//	if (prev_shift == 0)
//	{
//		min_shift = 0;
//		max_shift = giShiftDistance;
//	}
//	else
//	{
//		min_shift = radius - shift_freedom;
//		max_shift = radius + shift_freedom;
//		if (min_shift < 0)
//			min_shift = 0;
//		if (min_shift > giShiftDistance)
//			min_shift = giShiftDistance;
//		if (max_shift > giShiftDistance)
//			max_shift = giShiftDistance;
//		if (min_shift > max_shift)
//		{
//			return 0;
//		}
//	}
//
//	for (i = min_shift; i <= max_shift; i++)
//	{
//		for (j = 0; j < giFromLength; j++)
//		{
//			sum[i] += responses[i][j];
//			responses_queue[i].push(responses[i][j]);
//		}
//	}
//
//	float min_allow = 0.0;
//	bool stop_shift = false;
//
//	for (j = giFromLength; j <= giToLength; j++)
//	{
//		if (!stop_shift)
//		{
//			for (i = min_shift; i <= max_shift; i++)
//			{
//				sum[i] += responses[i][j];
//				responses_queue[i].push(responses[i][j]);
//
//				min_allow = -0.8 * fMaxDensity;
//
//				// is this length/shift the best density for the Horizontal Right Point?
//				if (sum[i] > min_allow)
//				{
//					fRespDensity = median(responses_queue[i]);
//					if (fRespDensity >= fMaxDensity)
//					{
//						fMaxDensity = fRespDensity;
//						if (best_point == NULL)
//						{
//							best_point = new CPoint(response_locations[i][j]);
//						}
//						else
//							*best_point = response_locations[i][j];
//						best_point->m_iValue = fRespDensity;
//					}
//				}
//				else
//				{
//					if (fMaxDensity > 0)
//					{
//						stop_shift = true;
//						break;
//					}
//				}
//			}
//		}
//
//		// At this template length, this is the best density
//		fMaxRespDensityAtThisLength = fMaxDensity;
//
//		if (fMaxRespDensityAtThisLength > fMaxRespDensity)
//		{
//			fMaxRespDensity = fMaxRespDensityAtThisLength;
//			TemplateLength = j;
//		}
//	}
//
//	if (best_point)
//	{
//		point = *best_point;
//		return true;
//	}
//	else
//		return false;
//}

int SortShortAscending(const void *pvArg1, const void *pvArg2)
{
	short iV1 = * static_cast<short *>(const_cast<void *>(pvArg1));
	short iV2 = * static_cast<short *>(const_cast<void *>(pvArg2));

	return iV1-iV2;
}

// new median calculation - inspired by Bob MBF
// Do quicksort for initial template length, then do insertion sort
bool IndividualMedianResponse(int prev_shift, CPoint& point,
	int& TemplateLength, short responses[MaxShiftDistance][MaxTemplateLength],
	CPoint response_locations[MaxShiftDistance][MaxTemplateLength])
{
	int iMaxShiftDistance = gConfig.GetMaximumShiftDistance();
	CPoint* best_point = NULL;

	register unsigned char i, j;
	register float fRespDensity = 0;
	register float fMaxRespDensity = 0;
	register float fMaxRespDensityAtThisLength = 0;
	register float fMaxDensity = 0;
	//static int diff = 0;

	register short sum[MaxShiftDistance];
	memset(sum, 0, sizeof(short) * MaxShiftDistance);

	int iRelativeShiftDistance = gConfig.GetRelativeShiftDistance();

	char min_shift, max_shift;
	int radius = Round((float) prev_shift / 2.0);

	int iFromLength = gConfig.GetMinimumTemplateLength();
	int iToLength = gConfig.GetMaximumTemplateLength();

	if (prev_shift == 0)
	{
		min_shift = 0;
		max_shift = static_cast<char>(iMaxShiftDistance);
	}
	else
	{
		min_shift = static_cast<char>(radius - iRelativeShiftDistance);
		max_shift = static_cast<char>(radius + iRelativeShiftDistance);
		if (min_shift < 0)
			min_shift = 0;
		if (min_shift > iMaxShiftDistance)
			min_shift = static_cast<char>(iMaxShiftDistance);
		if (max_shift > iMaxShiftDistance)
			max_shift = static_cast<char>(iMaxShiftDistance);
		if (min_shift > max_shift)
		{
			return 0;
		}
	}

	short **paiResponsesCopy = new short*[MaxShiftDistance];
	for(i = 0; i < MaxShiftDistance; i++) {
		paiResponsesCopy[i] = new short[MaxTemplateLength];
	}

	for (i = min_shift; i <= max_shift; i++)
	{
		for (j = 0; j < iFromLength; j++)
		{
			sum[i] = static_cast<short>(sum[i] + responses[i][j]);
		}
		memcpy(paiResponsesCopy[i],responses[i], sizeof (short)*iFromLength);
		qsort(paiResponsesCopy[i],iFromLength, sizeof(short), SortShortAscending);
	}

	float min_allow = 0.0;
	bool stop_shift = false;
	int d;

	for (j = static_cast<unsigned char>(iFromLength); j <= iToLength; j++)
	{
		if (!stop_shift)
		{
			for (i = min_shift; i <= max_shift; i++)
			{
				sum[i] = static_cast<short>(sum[i] + responses[i][j]);

				min_allow = static_cast<float>(-0.8 * fMaxDensity);
				
				int iToInsert = responses[i][j];
					for (d=j-1; d >= 0; d--)
						if ( paiResponsesCopy[i][d] > iToInsert )
						{
							paiResponsesCopy[i][d+1] = paiResponsesCopy[i][d];
						}
						else
						{
							paiResponsesCopy[i][d+1] = static_cast<short>(iToInsert);
							break;
						}
						if ( d < 0 )
							paiResponsesCopy[i][0] = static_cast<short>(iToInsert);

				// is this length/shift the best density for the Horizontal Right Point?
				if (sum[i] > min_allow)
				{
					int iMiddle = (j+1)/2;
					fRespDensity = static_cast<float>(((j+1)%2)? paiResponsesCopy[i][iMiddle] : (paiResponsesCopy[i][iMiddle-1]+paiResponsesCopy[i][iMiddle])/2);

					if (fRespDensity >= fMaxDensity)
					{
						fMaxDensity = fRespDensity;
						if (best_point == NULL)
						{
							best_point = new CPoint(response_locations[i][j]);
						}
						else
							*best_point = response_locations[i][j];
						best_point->m_iValue = static_cast<int>(fRespDensity);
					}
				}
				else
				{
					if (fMaxDensity > 0)
					{
						stop_shift = true;
						break;
					}
				}
			}
		}

		// At this template length, this is the best density
		fMaxRespDensityAtThisLength = fMaxDensity;

		if (fMaxRespDensityAtThisLength >= fMaxRespDensity)
		{
			fMaxRespDensity = fMaxRespDensityAtThisLength;
			TemplateLength = j;
		}
	}

	for(i = 0; i < MaxShiftDistance; i++) {
		delete paiResponsesCopy[i];
	}
	delete []paiResponsesCopy;

	if (best_point)
	{
		point = *best_point;
		return true;
	}
	else
		return false;
}
#endif

void Print2DArray(short responses[MaxShiftDistance][MaxTemplateLength])
{
	ofstream out("array.txt");
	for(int i = 0; i < MaxShiftDistance; i++) {
		for(int j = 0; j < MaxTemplateLength; j++) {
			out << responses[i][j] << "\t" << flush;
		}
		out << endl;
	}
	//	exit(0);
}

#ifndef rel_7_0
bool IndividualMedianResponse(int prev_shift, CPoint& point,
	int& TemplateLength, short responses[MaxShiftDistance][MaxTemplateLength],
	CPoint response_locations[MaxShiftDistance][MaxTemplateLength])
{
	//	Print2DArray(responses);
	CPoint* best_point = NULL;

	register unsigned char i, j;
	register float fRespDensity = 0;
	register float fMaxRespDensity = 0;
	register float fMaxRespDensityAtThisLength = 0;
	register float fMaxDensity = 0;

	register short sum[MaxShiftDistance];
	priority_queue<short> responses_queue[MaxShiftDistance];
	memset(sum, 0, sizeof(short) * MaxShiftDistance);

	int shift_freedom = cfg_trace_shift_freedom;
	char min_shift, max_shift;
	int radius = Round((float) prev_shift / 2.0);

	if (prev_shift == 0)
	{
		min_shift = 0;
		max_shift = giShiftDistance;
	}
	else
	{
		min_shift = radius - shift_freedom;
		max_shift = radius + shift_freedom;
		if (min_shift < 0)
			min_shift = 0;
		if (min_shift > giShiftDistance)
			min_shift = giShiftDistance;
		if (max_shift > giShiftDistance)
			max_shift = giShiftDistance;
		if (min_shift > max_shift)
		{
			return 0;
		}
	}
	
	register int shift, length;

	for (i = min_shift; i <= max_shift; i++)
	{
		for (j = 0; j < giFromLength; j++)
		{
			sum[i] += responses[i][j];
			responses_queue[i].push(responses[i][j]);
		}
	}

	float min_allow = 0.0;
	
	float max_resp_at_this_shift = 0.0;
	int best_length_at_this_shift = 0;
	float max_resp = 0.0;

	// shift from min_shift to max_shift
	for(shift = min_shift; shift <= max_shift; shift++) 
	{
		max_resp_at_this_shift = -9999;
		best_length_at_this_shift = 0;
		for(length = giFromLength; length <= giToLength; length++) 
		{
			sum[shift] += responses[shift][length];
			responses_queue[shift].push(responses[shift][length]);

			fRespDensity = median(responses_queue[shift]);
			
			if (fRespDensity >= max_resp_at_this_shift)
			{
				max_resp_at_this_shift = fRespDensity;
				best_length_at_this_shift = length;				
			}
		}
		
		if(max_resp_at_this_shift > max_resp) {
			max_resp = max_resp_at_this_shift;
			TemplateLength = best_length_at_this_shift;
			if (best_point == NULL)
			{
				best_point = new CPoint(response_locations[shift][best_length_at_this_shift]);
			}
			else
				*best_point = response_locations[shift][best_length_at_this_shift];
			best_point->m_iValue = max_resp;
		}
		min_allow = -0.8 * max_resp;
		if(max_resp_at_this_shift < min_allow && max_resp != 0.0)
			break;
		
	}

	if (best_point)
	{
		point = *best_point;
		return true;
	}
	else
		return false;
}
#endif

int GetMedianResponse(const CPoint& aPoint, CPoint& pHR, CPoint& pHL,
	CPoint& pVR, CPoint& pVL, int& TemplateLength)
{
	int HR_length, HL_length, VR_length, VL_length;

	if (!IndividualMedianResponse(Round(aPoint.m_fHWidth),
		 	pHR,
		 	HR_length,
		 	gaiHRResp,
		 	gaHRight) ||
		!IndividualMedianResponse(Round(aPoint.m_fHWidth),
						  	pHL,
						  	HL_length,
						  	gaiHLResp,
						  	gaHLeft) ||
		!IndividualMedianResponse(Round(aPoint.m_fVWidth),
										 	pVR,
										 	VR_length,
										 	gaiVRResp,
										 	gaVRight) ||
		!IndividualMedianResponse(Round(aPoint.m_fVWidth),
														  	pVL,
														  	VL_length,
														  	gaiVLResp,
														  	gaVLeft))
		return 0;

	list<int> lengths;
	if (HR_length > 0)
		lengths.push_back(HR_length);
	if (HL_length > 0)
		lengths.push_back(HL_length);
	if (VR_length > 0)
		lengths.push_back(VR_length);
	if (VL_length > 0)
		lengths.push_back(VL_length);

	lengths.sort();

	if (lengths.front() > 0)
	{
		TemplateLength = lengths.front();
	}
	else
	{
		return 0;
	}
	return (pHR.m_iValue + pHL.m_iValue + pVR.m_iValue + pVL.m_iValue);
}

