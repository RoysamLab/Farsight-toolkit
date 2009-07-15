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
// File: Vessel.cpp
// 
// This file contains the declaration of a vessel class. 
//
// - A vessel contains three linked list of CPoint instances: left, right, and center. 
//
// - Since these three linked lists might be of different lengths, there is no
//   association information at creation time. However, such informatino can be create
//   for each centerline point, by finding the corresponding left and right points in 
//   directions perpendicular to that of the centerline points.
//
// - A vessel keeps a list (an Array) of intersection points. Each intersection point 
//   contains location and other vessel ID information. 
//
// - A vessel is usually constructed from an instant of DLList<CProfile>. After
//   construction, the constructor checks if at the vessel intersects with any 
//   other vessels at each of its end. If it does, the appropriate points are added
//   to the corresponding vessel ends, and the intersection poiint(s) is(are) noted
//   and registered in the intersection points array. 
//
// - The length of the vessel is the lenfth of its center line, and its response at
//   any centerline point is the sum of the corresponding left and right points.
//
//#pragma warning(disable:4786)  // disable STL-related warnings
//#pragma warning(disable:4710)  // compiler didn't inline function that got selected for inline expansion
//#pragma warning(disable:4514)  // The optimizer removed an inline function that is not called
//#pragma warning(disable:4702)  // unreachable STLport code

#include <iostream>
#include <iomanip>
#include <fstream>
#include <list>
#include <queue>
#include <cmath>
#include <vector>
#include <map>

#include "CONSTANTS.h"
#include "Mytypes.h"  
#include "Config.h" 
#include "Cpoints.h"
#include "Dllist.h"
#include "Cimage.h"    // a simple image class
#include "C3dimage.h"
#include "Cvessel.h"
#include "Cvector.h"
#include "Template.h"  // a template class
#include "Extern.h"
#include "Numbers.h"
#include "Ctree.h"
#include "Soma.h"
#include "StrEle.h"

using namespace std;

#define LargeNumber 1000000
#define _INFINITY_ 9999999
int MaxNumOfAveragedPoints = 5;
int MaxExtensionLength = 20;
int LengthThreshold = 40; // don't write an ID for any vessel shorter than thes
int giSomaDistThreshold = 0; 

CLNode<CPoint>* FindClosestPoint(CDLList<CPoint>* aList, CPoint* aPoint);

extern CSomas* gTheSomas;
extern int giReversedFlag;

int inline ABS(int x)
{
	return ((x > 0) ? x : (0 - x));
}

/////////////////////////////
// Function: FindClosestIntersectionPoint
//
// Find and return the closest Intersection point
// this function is used to avoid selecting adjacent intersection points
// which will require merging later.
int FindClosestIntersectionPoint(CPoint& point, int& MinDistance)
{
	//int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;
	int pixelValue = 0;
	int distance = 99999;
	int IntersectionPointID = 0;
	int i, j, k;

	int iFrom = - 1 * giMaxExtensionLength;
	int iTo = giMaxExtensionLength;
	MinDistance = 99999;

	for (i = iFrom; i <= iTo; i++)
	{
		for (j = iFrom; j <= iTo; j++)
		{
			for (k = iFrom; k <= iTo; k++)
			{
				if ((i + point.m_iX) >= 0 &&
					(i + point.m_iX) < iCols &&
					(j + point.m_iY) >= 0 &&
					(j + point.m_iY) < iRows &&
					(j + point.m_iZ) >= 0 &&
					(j + point.m_iZ) < iRows)
				{
					pixelValue = Empty3DImage->data[j + point.m_iZ][j + point.m_iY][i + point.m_iX];
					if (pixelValue > StartingIntersectionColor &&
						pixelValue < StartingSomaColor)
					{
						distance = static_cast<int>( std::sqrt( static_cast<double>(i * i + j * j + k * k) ) );
						if (distance < MinDistance)
						{
							MinDistance = distance;
							IntersectionPointID = pixelValue -
								StartingIntersectionColor;
						}
					}
				}
			}
		}
	}

	return IntersectionPointID;
}

/////////////////////////////////////////////////////////////////////////////////
// METHOD: ExtendVessel
//
// Purpose: 
//	To extend the current vessel to include all points upto and including
//	the points in 'aProfile'.
// 
// Logic:
//	Use the draw line algorithm to find all such intervening points. However, 
//  rather than actually drawing the line, add its points as points to the link
//  list
//
// Called By: 
//	TrackAPoint() (see Track.cpp)
//
// Calls:
//	ExtendVessel (see below)
bool CVessel::ExtendVessel(CPoint* C, CPoint* HL, CPoint* HR, CPoint* VL,
	CPoint* VR, TopEndMiddleNeither DirFlag)
{
	if (giReversedFlag == 0)
	{
		// mark these points as original points (i.e. not filling)
		C->m_lUserFlag = 1;
		HL->m_lUserFlag = 1;
		HR->m_lUserFlag = 1;
		VL->m_lUserFlag = 1;
		VR->m_lUserFlag = 1;
	}

	// Furthermore, the flag at the centerline point means the width of the segment
	double xdiff = (HL->m_iX - HR->m_iX);   xdiff = xdiff * xdiff;
	double ydiff = (HL->m_iY - HR->m_iY);   ydiff = ydiff * ydiff;
	double zdiff = (HL->m_iZ - HR->m_iZ);   zdiff = zdiff * zdiff;

	C->m_lUserFlag = static_cast<long>(std::sqrt(xdiff + ydiff + zdiff));
	if (C->m_lUserFlag == 0)
		C->m_lUserFlag = 1;

	//	int prev_length = this->m_iLength;
	// if the point already exists, don't to the list
	CLNode<CPoint>* temp = m_Center.head;
	while (temp)
	{
		if (*C == *temp->data)
			return false;
		temp = temp->after;
	}

	if (ExtendVessel(m_Center, C, DirFlag))
	{
		ExtendVessel(m_HLeft, HL, DirFlag);
		ExtendVessel(m_HRight, HR, DirFlag);
		ExtendVessel(m_VLeft, VL, DirFlag);
		ExtendVessel(m_VRight, VR, DirFlag);
		UpdateLength();
		return true;
	}
	else
		return false;
}

/////////////////////////////////////////////////////////////////////////////////
// METHOD: ExtendVessel
//
// Purpose: 
//	To extend the current vessel to include all points upto and including 'aPoint',
//  see the 'ExtendVessel' method above. 
//
// Logic:
//	Use the drawing program to draw a line between two points, but rather than
//	drawing actual points, add them to the link list. The Direction flag "dirFlag"
//	specifies whether to add the elements on top or on the end of the link list
// 
// Notice that this method is essentially a line drawing function, just like those
// in the file "tools.cpp"
bool CVessel::ExtendVessel(CDLList<CPoint>& aList, CPoint* toPoint,
	TopEndMiddleNeither dirFlag)
{
	// a working CPoint instant that will be added to the link list
	int iToLength = gConfig.GetMaximumTemplateLength();
	CPoint fromPoint;
	if (dirFlag == OnTop)
	{
		if (aList.head)
		{
			fromPoint = *aList.head->data;
			if (fromPoint.FindDistance(toPoint) > iToLength)
				return false;
		}
		else
		{
			aList.AddElementOnTop(toPoint);
			return true;
		}
	}
	else
	{
		if (aList.tail)
		{
			fromPoint = *aList.tail->data;
			if (fromPoint.FindDistance(toPoint) > iToLength)
				return false;
		}
		else
		{
			aList.AddElementOnEnd(toPoint);
			return true;
		}
	}

	// make sure that the points are different
	if (fromPoint == *toPoint)
		return false;

	// don't extend when the distance is greater than min template length
	if (fromPoint.FindDistance(toPoint) > iToLength)
		return false;

	// a place holder for the resulting line
	CDLList<CPoint> aLine;
	// construct the line by linear interpolation (see the file tools.cpp)
	Construct3DLine(fromPoint, * toPoint, aLine);

	// the resulting line extends from "fromPoint" to "toPoint" inclusive
	// we need not add the point "fromPoint" to our list since it is
	// already there

	CLNode<CPoint>* tempNode = aLine.head->after;
	while (tempNode)
	{
		if (dirFlag == OnTop)
			aList.AddElementOnTop(tempNode->data);
		else
			aList.AddElementOnEnd(tempNode->data);

		tempNode = tempNode->after;
	}
	return true;
}

/////////////////////////////////////////////////////////////////////////////////
// METHOD: ConnectWithOtherVessels
//
// Purpose: 
//  To connect both ends of this vessel with other vessels if it intersects with
//  with any.

// Logic:
//	Check the regions surrounding the end of the vessel to see if the vessel
//	intersect with other vessels. Only those vessels lying along a line in 
//	in the direction of my vessel, and within LOOK_AHEAD distance from it are 
//	considered intersecting vessels.  

void CVessel::ConnectWithOtherVessels(TopEndMiddleNeither DirFlag)
{
	int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;
	int pixelValue = 0;
	int PathSum = 0;
	int BestLength = giMaxExtensionLength;
	CPoint* FromPoint = NULL;
	CPoint* Indices = NULL;
	CPoint FoundPoint;
	int newX, newY, newZ;
	int Hindex = 1;
	int Vindex = 1;
	int Hdirections[17];
	int Vdirections[17]; // NumOfDirections/2 + 1
	if (DirFlag == OnTop)
		FromPoint = m_Center.head->data;
	else
		FromPoint = m_Center.tail->data;

	int Hdir = FromPoint->m_iHDir;
	int Vdir = FromPoint->m_iVDir;

	Hdirections[0] = Hdir;
	Vdirections[0] = Vdir;

	int FoundVesselID = 0;

	register int i, j, k;

	for (i = 1; i < 9; i++)
	{
		Hdirections[Hindex] = DirectionPlus(Hdir, i);
		Hindex++;
		Hdirections[Hindex] = DirectionMinus(Hdir, i);
		Hindex++;

		Vdirections[Vindex] = DirectionPlus(Vdir, i);
		Vindex++;
		Vdirections[Vindex] = DirectionMinus(Vdir, i);
		Vindex++;
	}

	for (i = 0; i < 17; i++)
	{
		for (j = 0; j < 17; j++)
		{
			Indices = gVectorsArray[Hdirections[i]][Vdirections[j]]->m_pIndices;

			// the cost of this path is given by the sum of its gray level values
			// over the original image
			PathSum = 0;
			// we loop until the maximum vector length because in case of
			// wide vessels the centerline of such vessels might be more
			// than LOOK_AHEAD pixels away. 
			for (k = 1; k < giMaxExtensionLength; k++)
			{
				newX = FromPoint->m_iX + Indices[k].m_iX;
				newY = FromPoint->m_iY + Indices[k].m_iY;
				newZ = FromPoint->m_iZ + Indices[k].m_iZ;

				if (newX < 0 ||
					newX >= iCols ||
					newY < 0 ||
					newY >= iRows ||
					newZ < 0 ||
					newZ >= iSlices)
					break;

				pixelValue = Empty3DImage->data[newZ][newY][newX];	
				// at any point if the vessels intersects with itself, break
				// i.e. consider the next direction
				if (pixelValue == m_iID)
					break;

				PathSum += The3DImage->data[newZ][newY][newX];

				if (pixelValue)
				{
					// Consider this intersection only if it is better than the optimal
					if (k < BestLength)
					{
						BestLength = k;
						FoundPoint.m_iX = newX; 
						FoundPoint.m_iY = newY;
						FoundPoint.m_iZ = newZ;
						FoundPoint.m_iHDir = static_cast<unsigned char>(Hdirections[i]);
						FoundPoint.m_iVDir = static_cast<unsigned char>(Vdirections[j]);
						FoundVesselID = pixelValue;

						break;
					}
					else
						break;
				}
			}
		}
	}
	// this vessel can be extended to intersect with another vessel
	if (FoundVesselID)
	{
		/////////////////
		// 1. If the found point is near the edge of the found vessel
		//    or near an already existing intersection point, see if we
		//    want to snap to either.
		//    To do this, we find the distance from the found point, and the
		//    found vessel edges, say DistFromEdge. 
		//		Also, we find the distance between the found point and the
		//    closest intersection point, say DistFromClosestIntPoint
		//   1.1
		//   => If DistFromClosestIntPoint < DistFromEdge &&
		//  	   DistFromClosestIntPoint < threshold
		//  	   => extend this vessel to the already existing intersection
		//  		  point.
		//   1.2
		//   => Else If DistFromEdge < DistFromClosestIntPoint && 
		//  			DistFromEdge < threshold
		//				  => Take the found vessel edge as the location of the new 
		//  			   intersection point.
		//   1.3
		//   => Else 
		//  		extend this vessel to the found location.
		//

		int DistFromClosestIntPoint = 99999;
		int ClosestIntPointID = FindClosestIntersectionPoint(*FromPoint,
									DistFromClosestIntPoint);
		int DistFromEdge = 0;		
		CPoint* pEdgePoint = NULL;
		CVessel* pFoundVessel = NULL;
		TopEndMiddleNeither location = Neither;

		int d1 = 0;
		int d2 = 0;
		if ((pFoundVessel = gTheVessels.GetVessel(FoundVesselID)) == NULL)
			return ;

		d1 = static_cast<int>(FromPoint->FindDistance(pFoundVessel->m_Center.head->data));
		d2 = static_cast<int>(FromPoint->FindDistance(pFoundVessel->m_Center.tail->data));

		if (d1 < d2)
		{
			pEdgePoint = pFoundVessel->m_Center.head->data;
			DistFromEdge = d1;
			location = OnTop;
		}
		else
		{
			pEdgePoint = pFoundVessel->m_Center.tail->data;
			DistFromEdge = d2;
			location = OnEnd;
		}

		// 1.1, see above
		if (ClosestIntPointID &&
			DistFromClosestIntPoint <= giMaxExtensionLength &&
			DistFromClosestIntPoint <= DistFromEdge)
		{
			// add the vessel ID to the found intersection point
			gIntersectionPoints.m_apData[ClosestIntPointID - 1]->AddVessel(m_iID,
																 	DirFlag);
			// add the intersection point ID to the vessel
			AddIntersectionPoint(ClosestIntPointID);
			// extend the vessel
			ExtendVesselCenter(&(gIntersectionPoints.m_apData[ClosestIntPointID - 1]->m_Point),
				DirFlag);
		}
		// 1.2, see comments above
		else if (DistFromEdge < DistFromClosestIntPoint &&
			DistFromEdge < giMaxExtensionLength)
		{
			// create a new intersection point and add the vessel info to it
			CIntersectionPoint IntPoint(*pEdgePoint);
			IntPoint.AddVessel(m_iID, DirFlag);
			IntPoint.AddVessel(FoundVesselID, location);
			// add this point to the intersection point collection
			int pointID = gIntersectionPoints.Add(IntPoint);

			// add the intersection point to each of the vessels
			AddIntersectionPoint(pointID);
			pFoundVessel->AddIntersectionPoint(pointID);

			// extend this vessel
			ExtendVesselCenter(pEdgePoint, DirFlag);

			// mark the new intersection point in the image
			Empty3DImage->data[pEdgePoint->m_iZ][pEdgePoint->m_iY][pEdgePoint->m_iX] = 
				static_cast<unsigned char>(pointID + StartingIntersectionColor);
		}
		// 1.3, see comments above
		else
		{
			// create a new intersection point and add the vessel info to it
			CIntersectionPoint IntPoint(FoundPoint);
			IntPoint.AddVessel(m_iID, DirFlag);
			IntPoint.AddVessel(FoundVesselID, Middle);
			// add this point to the intersection points collection
			int pointID = gIntersectionPoints.Add(IntPoint);

			// add the intersection point to both vessels
			AddIntersectionPoint(pointID);
			pFoundVessel->AddIntersectionPoint(pointID);

			// extend this vessel
			ExtendVesselCenter(&FoundPoint, DirFlag);

			// mark the new intersection point in the image
			Empty3DImage->data[FoundPoint.m_iZ][FoundPoint.m_iY][FoundPoint.m_iX] = 
				static_cast<unsigned char>(pointID + StartingIntersectionColor);
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
// METHOD: ConnectWithSomaAndOtherVessels
//
// Purpose: 
//  To extend this vessel at either of its ends to connect it with other vessels or soma
//
// Logic:
//  First check if it does not hit another soma/vessel already (8-connected)
//  if it is not, try to extend the vessel to the nearest soma/vessel. 
//  If we can extend the vessel either to hit a vessel or a soma, we select the
//  whichever is closer
void CVessel::ConnectWithSomasAndOtherVessels()
{
	int aPointID = 0;  // the ID of an intersection point. 
	int FoundSomaID = 0;

	CPoint tempPoint;

	// find the possible directions of extension for either end, and extend 
	// the vessel in such directions
	// an extendable vessel is not 8-connected with a soma or an other vessel already
	if (Extendable(OnTop))
	{
		FoundSomaID = FindSomaDistance(OnTop, tempPoint);
		if (FoundSomaID)
		{
			// add this intersection point to the collection of intersection
			// points. And mark such point as a vessel-soma one.
			ExtendVessel(m_Center, & tempPoint, OnTop);
			CIntersectionPoint IntPoint(tempPoint);
			IntPoint.AddSoma(FoundSomaID);
			IntPoint.AddVessel(m_iID, OnTop);
			aPointID = gIntersectionPoints.Add(IntPoint);
			AddIntersectionPoint(aPointID);
			gTheSomas->m_aData[FoundSomaID - 1].AddIntersectionPoint(aPointID);
		}
		else
			ConnectWithOtherVessels(OnTop);
	}

	if (Extendable(OnEnd))
	{
		FoundSomaID = FindSomaDistance(OnEnd, tempPoint);

		if (FoundSomaID)
		{
			ExtendVessel(m_Center, & tempPoint, OnEnd);

			// add this intersection point to the collection of intersection
			// points. And mark such point as a vessel-soma one.
			CIntersectionPoint IntPoint(tempPoint);
			IntPoint.AddSoma(FoundSomaID);
			IntPoint.AddVessel(m_iID, OnTop);

			aPointID = gIntersectionPoints.Add(IntPoint);
			AddIntersectionPoint(aPointID);
			gTheSomas->m_aData[FoundSomaID - 1].AddIntersectionPoint(aPointID);
		}
		else
			ConnectWithOtherVessels(OnEnd);
	}
}

void CVessel::ConnectWithSomasAndOtherVessels2()
{
	// try to connect with the soma first from each end
	if (! ConnectWithSomas(OnTop))
		ConnectWithOtherVessels2(OnTop);

	if (! ConnectWithSomas(OnEnd))
		ConnectWithOtherVessels2(OnEnd);
}

vector<bool> CVessel::Connect()
{
	vector<bool> result;
	result.push_back(false);
	result.push_back(false);

		list<int> responses_Topthis; 
	list<int> responses_Endthis; 
	list<int> responses_pTopClosestVessel;
	list<int> responses_pEndClosestVessel;
	CLNode<CPoint>* tempNode = NULL;
	
	int sum_thistop = 0, count_thistop = 0;
	int sum_thisend = 0, count_thisend = 0;
	int sum_top = 0, count_top = 0;
	int sum_end = 0, count_end = 0;
	
	tempNode = this->m_Center.head;
	while (tempNode)
	{
		if(tempNode->data->m_iValue)
		{
			responses_Topthis.push_back(tempNode->data->m_iValue);
			sum_thistop += tempNode->data->m_iValue;
			count_thistop++;
		}
		if(count_thistop > 7)
			break;
		tempNode = tempNode->after;
	}

	tempNode = this->m_Center.tail;
	while (tempNode)
	{
		if(tempNode->data->m_iValue)
		{
			responses_Endthis.push_back(tempNode->data->m_iValue);
			sum_thisend += tempNode->data->m_iValue;
			count_thisend++;
		}
		if(count_thisend > 7)
			break;
		tempNode = tempNode->before;
	}

	tempNode = m_pTopClosestVessel->m_Center.head;
	while (tempNode)
	{
		if(tempNode->data->m_iValue)
		{
			responses_pTopClosestVessel.push_back(tempNode->data->m_iValue);
			sum_top += tempNode->data->m_iValue;
			count_top++;
		}
		tempNode = tempNode->after;
	}
	tempNode = m_pEndClosestVessel->m_Center.head;
	while (tempNode)
	{
		if(tempNode->data->m_iValue)
		{
			responses_pTopClosestVessel.push_back(tempNode->data->m_iValue);
			sum_end += tempNode->data->m_iValue;
			count_end++;
		}
		tempNode = tempNode->after;
	}

	int median_Topthis = median(responses_Topthis);
	int median_Endthis = median(responses_Endthis);

	int median_pTopClosestVessel = median(responses_pTopClosestVessel);
	int median_pEndClosestVessel = median(responses_pEndClosestVessel);
	sum_thistop /= count_thistop;
	sum_thisend /= count_thisend;
	sum_top /= count_top;
	sum_end /= count_end;

	bool top_intensity = false;
	bool end_intensity = false;

	float thresh = 0.6f;
	if(min(median_Topthis, median_pTopClosestVessel) > thresh * max(median_Topthis, median_pTopClosestVessel)
		&& min(sum_thistop, sum_top) > thresh * max(sum_thistop, sum_top) )
		top_intensity = true;
	if(min(median_Endthis, median_pEndClosestVessel) > thresh * max(median_Endthis, median_pEndClosestVessel)
		&& min(sum_thisend, sum_end) > thresh * max(sum_thisend, sum_end) )
		end_intensity = true;

	// check the width of end points
	//float iTopDistToVesselThreshold = max(this->m_Center.head->data->m_fHWidth * 1.33f, m_pTopClosestNode->data->m_fHWidth * 1.33f);
	// check the directions of the two closest point

	// if end-to-end
	CLNode<CPoint>* pTempNode1 = NULL;
	CLNode<CPoint>* pTempNode2 = NULL;
	bool bTopEndToEnd = (m_pTopClosestNode->data == m_pTopClosestVessel->m_Center.head->data 
		|| m_pTopClosestNode->data == m_pTopClosestVessel->m_Center.tail->data);
	if(bTopEndToEnd)
	{
		this->FindClosestCenterlinePoint(*m_pTopClosestVessel->m_Center.head->data, & pTempNode1);
		this->FindClosestCenterlinePoint(*m_pTopClosestVessel->m_Center.tail->data, & pTempNode2);
		bTopEndToEnd = (bTopEndToEnd && ((*pTempNode1->data == *this->m_Center.head->data) || (*pTempNode2->data == *this->m_Center.head->data)));
	}
	int TopDir = this->m_Center.head->data->m_iHDir;
	if(TopDir > NumOfDirections/2)
		TopDir -= NumOfDirections/2;
	int TopClosestDir = m_pTopClosestNode->data->m_iHDir;
	if(TopClosestDir > NumOfDirections/2)
		TopClosestDir -= NumOfDirections/2;
	bool bTopNotParallel = (DirDistance(TopDir, TopClosestDir) >= NumOfDirections/8 - 1
		&& DirDistance(TopDir, TopClosestDir) <= NumOfDirections/8 * 3 + 1);
//
//	// check if the intensity values of vessels to be connected is close enough
//	// get the median response of both vessels (do_connect)

	bool bTopConnect = (!top_intensity && !bTopNotParallel && bTopEndToEnd) 
		|| (!top_intensity && bTopNotParallel && bTopEndToEnd) 
		|| (top_intensity && !bTopNotParallel && bTopEndToEnd) 
		|| (top_intensity && bTopNotParallel && !bTopEndToEnd)  
		|| (top_intensity && bTopNotParallel && bTopEndToEnd);

	pTempNode1 = NULL;
	pTempNode2 = NULL;
	//float iEndDistToVesselThreshold = max(this->m_Center.tail->data->m_fHWidth * 1.33f, m_pEndClosestNode->data->m_fHWidth * 1.33f);
	// check the directions of the two closest points
	bool bEndEndToEnd = (m_pEndClosestNode->data == m_pEndClosestVessel->m_Center.head->data 
		|| m_pEndClosestNode->data == m_pEndClosestVessel->m_Center.tail->data);
	if(bEndEndToEnd)
	{
		this->FindClosestCenterlinePoint(*m_pEndClosestVessel->m_Center.head->data, & pTempNode1);
		this->FindClosestCenterlinePoint(*m_pEndClosestVessel->m_Center.tail->data, & pTempNode2);
		bEndEndToEnd = (bEndEndToEnd && ((*pTempNode1->data == *this->m_Center.tail->data) || (*pTempNode2->data == *this->m_Center.tail->data)));
	}
	int EndDir = this->m_Center.tail->data->m_iHDir;
	if(EndDir > NumOfDirections / 2)
		EndDir -= NumOfDirections/2;
	int ClosestEndDir = m_pEndClosestNode->data->m_iHDir;
	if(ClosestEndDir> NumOfDirections / 2)
		ClosestEndDir -= NumOfDirections/2;
	bool bEndNotParallel = (DirDistance(EndDir, ClosestEndDir) >= NumOfDirections/8 - 1
		&& DirDistance(EndDir, ClosestEndDir) <= NumOfDirections/8 * 3 + 1);
	
	bool bEndConnect = (!end_intensity && !bEndNotParallel && bEndEndToEnd) 
		|| (!end_intensity && bEndNotParallel && bEndEndToEnd) 
		|| (end_intensity && !bEndNotParallel && bEndEndToEnd) 
		|| (end_intensity && bEndNotParallel && !bEndEndToEnd)  
		|| (end_intensity && bEndNotParallel && bEndEndToEnd);

	result[0] = bTopConnect;
	result[1] = bEndConnect;

	return result;
}

void CVessel::ConnectWithSomasAndOtherVessels3()
{
	// amri: use the width info at vessel tips as distance threshold
	//	int iVesselDistThreshold = iFromLength + 2;
	int iTopDistToSoma = _INFINITY_;
	int iEndDistToSoma = _INFINITY_;
	int iTopDistToVessel = _INFINITY_;
	int iEndDistToVessel = _INFINITY_;
	int iTopSomaID = -1, iEndSomaID = -1;	
	CPoint TopSomaPoint, EndSomaPoint;	
	CIntersectionPoint* pIntersectionPointOnTop = NULL;
	CIntersectionPoint* pIntersectionPointOnEnd = NULL;

	// flags
	int iTopConnectedToSoma = 0;
	int iEndConnectedToSoma = 0;
	int id;
	register int i;

	m_iConnectOnTopFlag = 0;
	m_iConnectOnEndFlag = 0;
	m_pTopClosestNode = NULL;
	m_pTopClosestVessel = NULL;
	m_pEndClosestNode = NULL;
	m_pEndClosestVessel = NULL;

	// if the vessel is not connected with anything on top, find
	// the distances and locations of the closest somas and other vessels
	// from the top.
	if (Float(OnTop))
	{
		iTopDistToSoma = FindSomaDistance3(OnTop, TopSomaPoint, iTopSomaID);
		iTopDistToVessel = gTheVessels.FindClosestVessel(m_iID, m_Center.head->data, & m_pTopClosestVessel, & m_pTopClosestNode);
	}
	else
	{
		// get a pointer to the intersection point on top
		for (i = 0; i < m_iNumOfIntersectionPoints; i++)
		{
			id = m_aiMyIntersectionPoints[i];
			if (gIntersectionPoints.m_apData[id - 1]->GetLocation(m_iID) == OnTop)
			{
				pIntersectionPointOnTop = gIntersectionPoints.m_apData[id - 1];
				break;
			}
		}
	}

	// if the vessel is not connected with anything on End, find
	// the distances and locations of the closest somas and other vessels
	// from the End.
	if (Float(OnEnd))
	{
		iEndDistToSoma = FindSomaDistance3(OnEnd, EndSomaPoint, iEndSomaID);
		iEndDistToVessel = gTheVessels.FindClosestVessel(m_iID,
									   	m_Center.tail->data,
									   	& m_pEndClosestVessel,
									   	& m_pEndClosestNode);
	}
	else
	{
		// get a pointer to the intersection point on end
		for (i = 0; i < m_iNumOfIntersectionPoints; i++)
		{
			id = m_aiMyIntersectionPoints[i];
			if (gIntersectionPoints.m_apData[id - 1]->GetLocation(m_iID) ==
				OnEnd)
			{
				pIntersectionPointOnEnd = gIntersectionPoints.m_apData[id - 1];
				break;
			}
		}
	}

//	if(iTopSomaID == -1)
//	{
//		cerr << "local variable 'iTopSomaID' may be used without having been initialized" << endl;
//		exit(0);
//	}
//
//	if(iEndSomaID == -1)
//	{
//		cerr << "local variable 'iEndSomaID' may be used without having been initialized" << endl;
//		exit(0);
//	}


	// Find the shortest soma distance from either end
	if (iTopDistToSoma < giSomaDistThreshold &&
		iTopDistToSoma <= iEndDistToSoma &&
		iTopDistToSoma <= 3 * iTopDistToVessel)
	{
		// check that the other end is not connected to the same soma
		if (pIntersectionPointOnEnd)
		{
			if (pIntersectionPointOnEnd->m_iSomaID != iTopSomaID)
			{
				ConnectWithSoma(iTopSomaID, TopSomaPoint, OnTop);
				iTopConnectedToSoma = 1;
			}
		}
		else
		{
			ConnectWithSoma(iTopSomaID, TopSomaPoint, OnTop);
			iTopConnectedToSoma = 1;
		}
	}
	else if (iEndDistToSoma < giSomaDistThreshold &&
		iEndDistToSoma < iTopDistToSoma &&
		iEndDistToSoma <= 2.5 * iEndDistToVessel)
	{
		// check that the other end is not connected to the same soma
		if (pIntersectionPointOnTop)
		{
			if (pIntersectionPointOnTop->m_iSomaID != iEndSomaID)
			{
				ConnectWithSoma(iEndSomaID, EndSomaPoint, OnEnd);
				iEndConnectedToSoma = 1;
			}
		}
		else
		{
			ConnectWithSoma(iEndSomaID, EndSomaPoint, OnEnd);
			iEndConnectedToSoma = 1;
		}
	}

	// if the top end is free, try to connect to other vessels
	// as long as the other end is not connected to the same vessel, or
	// it should be.
	vector<bool> connect = Connect();

	float iTopDistToVesselThreshold = max(this->m_Center.head->data->m_fHWidth * 1.33f, m_pTopClosestNode->data->m_fHWidth * 1.33f);
	if (iTopConnectedToSoma == 0 && iTopDistToVessel < iTopDistToVesselThreshold && connect[0])
	{
		// also means that this end is free

		// if the other end has an intersection point, make sure
		// that is not connected with the same vessel I want to connect
		// with
		gLogFile << "OnTop: Vessel " << this->m_iID << " with Vessel " << m_pTopClosestVessel->m_iID << endl;
		if (pIntersectionPointOnEnd && m_pEndClosestVessel)
		{
			if (pIntersectionPointOnEnd->GetLocation(m_pEndClosestVessel->m_iID) == Neither)
				m_iConnectOnTopFlag = 1;
		}
		// the other end has no intersection point
		else
		{
			// if both ends wants to connect with the same vessel
			// connect the closer one
			if (m_pEndClosestVessel && m_pTopClosestVessel->m_iID == m_pEndClosestVessel->m_iID)
			{
				if (iTopDistToVessel <= iEndDistToVessel)
					m_iConnectOnTopFlag = 1;
			}
			else
				m_iConnectOnTopFlag = 1;
		}
	}

	float iEndDistToVesselThreshold = max(this->m_Center.tail->data->m_fHWidth * 1.33f, m_pEndClosestNode->data->m_fHWidth * 1.33f);
	if (iEndConnectedToSoma == 0 && iEndDistToVessel < iEndDistToVesselThreshold && connect[1])
	{
		// if the other end has an intersection point, make sure
		// that is not connected with the same vessel I want to connect
		// with
		gLogFile << "OnEnd: Vessel " << this->m_iID << " with Vessel " << m_pEndClosestVessel->m_iID << endl;
		if (pIntersectionPointOnTop && m_pTopClosestVessel)
		{
			if (pIntersectionPointOnTop->GetLocation(m_pTopClosestVessel->m_iID) ==
				Neither)
				m_iConnectOnEndFlag = 1;
		}
		// the other end has no intersection point
		else
		{
			// if both ends wants to connect with the same vessel
			// connect the closer one
			if (m_pTopClosestVessel &&
				m_pTopClosestVessel->m_iID == m_pEndClosestVessel->m_iID)
			{
				if (iEndDistToVessel < iTopDistToVessel)
					m_iConnectOnEndFlag = 1;
			}
			else
				m_iConnectOnEndFlag = 1;
		}
	}
}

// connect the vessel with the given soma point at the specified end
void CVessel::ConnectWithSoma(int iSomaID, CPoint& SomaPoint,
	TopEndMiddleNeither DirFlag)
{
	// add this intersection point to the collection of intersection
	// points. And mark such point as a vessel-soma one.
	// extend the vessel
	ExtendVesselCenter(&SomaPoint, DirFlag);
	CIntersectionPoint IntPoint(SomaPoint);
	IntPoint.AddSoma(iSomaID);
	IntPoint.AddVessel(m_iID, DirFlag);

	int iIntPointID = gIntersectionPoints.Add(IntPoint);
	AddIntersectionPoint(iIntPointID);
	gTheSomas->m_aData[iSomaID - 1].AddIntersectionPoint(iIntPointID);
}

int CVessel::ConnectWithSomas(TopEndMiddleNeither DirFlag)
{
	int result = 0;
	giMaxExtensionLength = 20; //iFromLength+2;
	CPoint SomaPoint;
	int iClosestSomaID = 0;
	// we need to prevent situations where the vessels is connected with the
	// soma from both ends

	if (iClosestSomaID)
	{
		// add this intersection point to the collection of intersection
		// points. And mark such point as a vessel-soma one.
		// extend the vessel
		ExtendVesselCenter(&SomaPoint, DirFlag);

		CIntersectionPoint IntPoint(SomaPoint);
		IntPoint.AddSoma(iClosestSomaID);
		IntPoint.AddVessel(m_iID, DirFlag);

		int iIntPointID = gIntersectionPoints.Add(IntPoint);
		AddIntersectionPoint(iIntPointID);
		gTheSomas->m_aData[iClosestSomaID - 1].AddIntersectionPoint(iIntPointID);
		result = 1;
	}

	return result;
}

void CVessel::ConnectWithOtherVessel(TopEndMiddleNeither DirFlag, int ,
	CVessel* pClosestVessel, CLNode<CPoint>* pClosestNode,
	CIntersectionPoint* pIntPoint)
{

	//int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;
	int iFromLength = gConfig.GetMinimumTemplateLength();

	CPoint* FromPoint = m_Center.head->data;
	if (DirFlag == OnEnd)
		FromPoint = m_Center.tail->data;

	int iDistFromClosestIntPoint = iCols* iRows;
	CIntersectionPoint* pClosestIntPoint = NULL;

	if (pIntPoint)
		iDistFromClosestIntPoint = gIntersectionPoints.FindClosestIntersectionPointExcluding(*FromPoint,
													   	& pClosestIntPoint,
													   	pIntPoint->m_iID);
	else
		// 1. Find the closest already existing intersection point
		iDistFromClosestIntPoint = gIntersectionPoints.FindClosestIntersectionPoint(*FromPoint,
													   	& pClosestIntPoint);

	// 1.1 If the point is closet than threshold, go for it
	if (iDistFromClosestIntPoint <= iFromLength)
	{
		// add the vessel ID to the found intersection point
		pClosestIntPoint->AddVessel(m_iID, DirFlag);
		// add the intersection point ID to the vessel
		AddIntersectionPoint(pClosestIntPoint->m_iID);
		// extend the vessel
		ExtendVesselCenter(&(pClosestIntPoint->m_Point), DirFlag);
	}
	else
	{
		// 2.1.1 Find the distance from the closest node and
		// its vessel edges
		int d1 = static_cast<int>(pClosestNode->data->FindDistance(pClosestVessel->m_Center.head->data));
		int d2 = static_cast<int>(pClosestNode->data->FindDistance(pClosestVessel->m_Center.tail->data));

		TopEndMiddleNeither Location = Middle;

		if (d1 < d2 && d1 <= iFromLength / 2)
		{
			pClosestNode = pClosestVessel->m_Center.head;
			Location = OnTop;
		}
		else if (d2 < d1 && d2 <= iFromLength / 2)
		{
			pClosestNode = pClosestVessel->m_Center.tail;
			Location = OnEnd;
		}

		// connect with other vessel at the specified point

		// create a new intersection point and add the vessel info to it
		CIntersectionPoint IntPoint(*(pClosestNode->data));
		IntPoint.AddVessel(m_iID, DirFlag);
		IntPoint.AddVessel(pClosestVessel->m_iID, Location);
		// add this point to the intersection points collection
		int pointID = gIntersectionPoints.Add(IntPoint);

		// add the intersection point to both vessels
		AddIntersectionPoint(pointID);
		pClosestVessel->AddIntersectionPoint(pointID);

		// extend this vessel
		ExtendVesselCenter(pClosestNode->data, DirFlag);
	}
}

void CVessel::ConnectWithOtherVessels2(TopEndMiddleNeither DirFlag)
{
	//int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;
	int iFromLength = gConfig.GetMinimumTemplateLength();
	CPoint* FromPoint = m_Center.head->data;
	if (DirFlag == OnEnd)
		FromPoint = m_Center.tail->data;

	// otherwise try to connect with other vessels
	int iDistFromClosestIntPoint = iCols* iRows;
	CIntersectionPoint* pClosestIntPoint = NULL;
	CLNode<CPoint>* pClosestNode = NULL;
	CVessel* pClosestVessel = NULL;
	int iDistFromClosestVessel = 0;

	// 1. Find the closest already existing intersection point	
	iDistFromClosestIntPoint = gIntersectionPoints.FindClosestIntersectionPoint(*FromPoint,
												   	& pClosestIntPoint);

	// 1.1 If the point is closet than threshold, go for it
	if (iDistFromClosestIntPoint < iFromLength + 2)
	{
		// add the vessel ID to the found intersection point
		pClosestIntPoint->AddVessel(m_iID, DirFlag);
		// add the intersection point ID to the vessel
		AddIntersectionPoint(pClosestIntPoint->m_iID);
		// extend the vessel
		ExtendVesselCenter(&(pClosestIntPoint->m_Point), DirFlag);

		return;
	}

	// 2. Find the closest point on any other vessel
	iDistFromClosestVessel = gTheVessels.FindClosestVessel(m_iID,
										 	FromPoint,
										 	& pClosestVessel,
										 	& pClosestNode);

	// 2.1
	if (iDistFromClosestVessel < iFromLength + 2)
	{
		// 2.1.1 Find the distance from the closest node and
		// its vessel edges
		int d1 = static_cast<int>(pClosestNode->data->FindDistance(pClosestVessel->m_Center.head->data));
		int d2 = static_cast<int>(pClosestNode->data->FindDistance(pClosestVessel->m_Center.tail->data));

		TopEndMiddleNeither Location = Middle;

		if (d1 < d2 && d1 < iFromLength)
		{
			pClosestNode = pClosestVessel->m_Center.head;
			Location = OnTop;
		}
		else if (d2 < d1 && d2 < iFromLength)
		{
			pClosestNode = pClosestVessel->m_Center.tail;
			Location = OnEnd;
		}

		// connect with other vessel at the specified point

		// create a new intersection point and add the vessel info to it
		CIntersectionPoint IntPoint(*(pClosestNode->data));
		IntPoint.AddVessel(m_iID, DirFlag);
		IntPoint.AddVessel(pClosestVessel->m_iID, Location);
		// add this point to the intersection points collection
		int pointID = gIntersectionPoints.Add(IntPoint);

		// add the intersection point to both vessels
		AddIntersectionPoint(pointID);
		pClosestVessel->AddIntersectionPoint(pointID);

		// extend this vessel
		ExtendVesselCenter(pClosestNode->data, DirFlag);
	}
}

// The possible connection scenarios
//
// 1. This wants to connect from both ends
//    1.1 With the same Vessel, V1
//  	  1.1.A. V1 wants to connect with This from both ends
//  	  1.1.B. V1 wants to connect with This from one end
//  	  1.1.C. V1 Does NOT want to connect with THIS.
//
//    1.2 With two Different Vessels V1, V2
//  	  1.2.A. V1/V2 Wants to connect with THIS from both ends
//  	  1.2.B. V1/V2 Wants to connect with THIS from one end only
//  	  1.2.C. V1/V2 Does NOT want to connect with THIS.
//
// 2. This wants to connect from one end only with V1
//    2.1. V1 wants to connect with THIS from both ends
//    2.2. V1 wants to connect with THIS from one end only
//    2.3. V1 Does NOT want to connect with THIS
//
void CVessel::ConnectWithOtherVessels3(CIntersectionPoint* pIntPoint)
{
	int d1 = _INFINITY_;
	int d2 = _INFINITY_;
	int d3 = _INFINITY_;
	int d4 = _INFINITY_;
	int d5 = _INFINITY_;
	int d6 = _INFINITY_;

	CVessel* V1 = NULL;
	CVessel* V2 = NULL;

	//1.  This wants to connect from both ends
	if (m_iConnectOnTopFlag && m_iConnectOnEndFlag)
	{
		m_iConnectOnTopFlag = 0;
		m_iConnectOnEndFlag = 0;

		//1.1 This wants to connect from both ends with the same vessel, V1
		if (m_pTopClosestVessel == m_pEndClosestVessel)
		{
			V1 = m_pTopClosestVessel;

			// Measure the distance from This to V1 and From V1 to This			
			d1 = static_cast<int>(m_Center.head->data->FindDistance(m_pTopClosestNode->data));
			d2 = static_cast<int>(m_Center.tail->data->FindDistance(m_pEndClosestNode->data));

			// 1.1.A. V1 wants to connect with This from both ends as well
			if (V1->m_pTopClosestVessel == V1->m_pEndClosestVessel &&
				V1->m_pTopClosestVessel->m_iID == m_iID)
			{
				// Measure the distance from V1 to This
				d3 = static_cast<int>(V1->m_Center.head->data->FindDistance(V1->m_pTopClosestNode->data));
				d4 = static_cast<int>(V1->m_Center.tail->data->FindDistance(V1->m_pEndClosestNode->data));

				V1->m_iConnectOnTopFlag = 0;
				V1->m_iConnectOnEndFlag = 0;
			}
			// 1.1.B V1 wants to connect with This from one end (On Top)
			else if (V1->m_pTopClosestVessel == this)
			{
				// Measure the distance from V1's Top to This
				d3 = static_cast<int>(V1->m_Center.head->data->FindDistance(V1->m_pTopClosestNode->data));
				V1->m_iConnectOnTopFlag = 0;
			}
			// 1.1.B  wants to connect with This from one end (On End)
			else if (V1->m_pEndClosestVessel == this)
			{
				// Measure the distance from V1's End to This
				d4 = static_cast<int>(V1->m_Center.tail->data->FindDistance(V1->m_pEndClosestNode->data));
				V1->m_iConnectOnEndFlag = 0;
			}
			// 1.1.C V1 Does not want to connect with This at all, do nothig

			// Connect only one end, the closest
			if (d1 <= d2 && d1 <= d3 && d1 <= d4)
				ConnectWithOtherVessel(OnTop,
					d1,
					V1,
					m_pTopClosestNode,
					pIntPoint);
			else if (d2 <= d1 && d2 <= d3 && d2 <= d4)
				ConnectWithOtherVessel(OnEnd,
					d2,
					V1,
					m_pEndClosestNode,
					pIntPoint);
			else if (d3 <= d1 && d3 <= d2 && d3 <= d4)
				V1->ConnectWithOtherVessel(OnTop,
						d3,
						this,
						V1->m_pTopClosestNode,
						pIntPoint);
			else if (d4 <= d1 && d4 <= d2 && d4 <= d3)
				V1->ConnectWithOtherVessel(OnEnd,
						d4,
						this,
						V1->m_pEndClosestNode,
						pIntPoint);
		}
		// 1.2 This wants to connect with two Different vessels, V1, and V2
		else
		{
			V1 = m_pTopClosestVessel;
			V2 = m_pEndClosestVessel;

			// Measure the distance from This to V1 and V2
			d1 = static_cast<int>(m_Center.head->data->FindDistance(m_pTopClosestNode->data));
			d2 = static_cast<int>(m_Center.tail->data->FindDistance(m_pEndClosestNode->data));

			// 1.2.A V1 wants to connect with This from both ends
			if (V1->m_pTopClosestVessel == V1->m_pEndClosestVessel &&
				V1->m_pTopClosestVessel == this)
			{
				// Measure the distances from V1 to This
				d3 = static_cast<int>(V1->m_Center.head->data->FindDistance(V1->m_pTopClosestNode->data));
				d4 = static_cast<int>(V1->m_Center.tail->data->FindDistance(V1->m_pEndClosestNode->data));
				V1->m_iConnectOnTopFlag = 0;
				V1->m_iConnectOnEndFlag = 0;
			}
			// 1.2.B V1 wants to connect with This from one end (onTop)
			else if (V1->m_pTopClosestVessel == this)
			{
				// Measure the distance from V1's Top to This
				d3 = static_cast<int>(V1->m_Center.head->data->FindDistance(V1->m_pTopClosestNode->data));
				V1->m_iConnectOnTopFlag = 0;
			}
			// 1.2.B V1 wants to connect with This from one end (OnEnd)
			else if (V1->m_pEndClosestVessel == this)
			{
				// Measure the distance from V1's End to This
				d4 = static_cast<int>(V1->m_Center.tail->data->FindDistance(V1->m_pEndClosestNode->data));
				V1->m_iConnectOnEndFlag = 0;
			}

			// 1.2.A V2 wants to connect with this From both ends
			if (V2->m_pTopClosestVessel == V2->m_pEndClosestVessel &&
				V2->m_pTopClosestVessel == this)
			{
				// Measure the distances from V2 to This
				d5 = static_cast<int>(V2->m_Center.head->data->FindDistance(V2->m_pTopClosestNode->data));
				d6 = static_cast<int>(V2->m_Center.tail->data->FindDistance(V2->m_pEndClosestNode->data));
				V2->m_iConnectOnTopFlag = 0;
				V2->m_iConnectOnEndFlag = 0;
			}
			// 1.2.B V1 wants to connect with This from one end (onTop)
			else if (V2->m_pTopClosestVessel == this)
			{
				// Measure the distance from V2's Top to This
				d5 = static_cast<int>(V2->m_Center.head->data->FindDistance(V2->m_pTopClosestNode->data));
				V2->m_iConnectOnTopFlag = 0;
			}
			// 1.2.B V1 wants to connect with This from one end (OnEnd)
			else if (V2->m_pEndClosestVessel == this)
			{
				// Measure the distance from V2's End to This
				d6 = static_cast<int>(V2->m_Center.tail->data->FindDistance(V2->m_pEndClosestNode->data));
				V2->m_iConnectOnEndFlag = 0;
			}

			// connect The Top of This using the best connection
			if (d1 <= d3 && d1 <= d4)
				ConnectWithOtherVessel(OnTop,
					d1,
					V1,
					m_pTopClosestNode,
					pIntPoint);
			else if (d3 <= d1 && d3 <= d4)
				V1->ConnectWithOtherVessel(OnTop,
						d3,
						this,
						V1->m_pTopClosestNode,
						pIntPoint);
			else if (d4 <= d1 && d4 <= d3)
				V1->ConnectWithOtherVessel(OnEnd,
						d4,
						this,
						V1->m_pEndClosestNode,
						pIntPoint);

			// Connect the End of This using the best connection
			if (d2 <= d5 && d2 <= d6)
				ConnectWithOtherVessel(OnEnd,
					d2,
					V2,
					m_pEndClosestNode,
					pIntPoint);
			else if (d5 <= d2 && d5 <= d6)
				V2->ConnectWithOtherVessel(OnTop,
						d5,
						this,
						V2->m_pTopClosestNode,
						pIntPoint);
			else if (d6 <= d2 && d6 <= d5)
				V2->ConnectWithOtherVessel(OnEnd,
						d6,
						this,
						V2->m_pEndClosestNode,
						pIntPoint);
		}
	}
	// 2. This wants to connect from one end only with V1
	else if (m_iConnectOnTopFlag || m_iConnectOnEndFlag)
	{
		if (m_iConnectOnTopFlag)
		{
			V1 = m_pTopClosestVessel;
			m_iConnectOnTopFlag = 0;
			// Measure the distance from This to V1
			d1 = static_cast<int>(m_Center.head->data->FindDistance(m_pTopClosestNode->data));
		}
		else if (m_iConnectOnEndFlag)
		{
			V1 = m_pEndClosestVessel;
			m_iConnectOnEndFlag = 0;
			d2 = static_cast<int>(m_Center.tail->data->FindDistance(m_pEndClosestNode->data));
		}

		// 2.A V1 wants to connect with This from Both Ends
		if (V1->m_pTopClosestVessel == V1->m_pEndClosestVessel &&
			V1->m_pTopClosestVessel == this)
		{
			// Measure the distance from V1 to This
			d3 = static_cast<int>(V1->m_Center.head->data->FindDistance(V1->m_pTopClosestNode->data));
			d4 = static_cast<int>(V1->m_Center.tail->data->FindDistance(V1->m_pEndClosestNode->data));
			V1->m_iConnectOnTopFlag = 0;
			V1->m_iConnectOnEndFlag = 0;
		}
		// 2.B V1 wants to connect with This from one end only (OnTop)
		else if (V1->m_pTopClosestVessel == this)
		{
			// Measure the distance from V1's Top to This
			d3 = static_cast<int>(V1->m_Center.head->data->FindDistance(V1->m_pTopClosestNode->data));
			V1->m_iConnectOnTopFlag = 0;
		}
		// 2.B V1 wants to connect with This from one end only (OnEnd)
		else if (V1->m_pEndClosestVessel == this)
		{
			// Measure the distance from V1's End to This
			d4 = static_cast<int>(V1->m_Center.tail->data->FindDistance(V1->m_pEndClosestNode->data));
			V1->m_iConnectOnEndFlag = 0;
		}
		// else V1 does not want to connect with This, do nothing

		// Connect the Top of this using the best connection
		if (d1 <= d2 && d1 <= d3 && d1 <= d4)
			ConnectWithOtherVessel(OnTop, d1, V1, m_pTopClosestNode, pIntPoint);
		else if (d2 <= d1 && d2 <= d3 && d2 <= d4)
			ConnectWithOtherVessel(OnEnd, d2, V1, m_pEndClosestNode, pIntPoint);
		else if (d3 <= d1 && d3 <= d2 && d3 <= d4)
			V1->ConnectWithOtherVessel(OnTop,
					d3,
					this,
					V1->m_pTopClosestNode,
					pIntPoint);
		else if (d4 <= d1 && d4 <= d2 && d4 <= d3)
			V1->ConnectWithOtherVessel(OnEnd,
					d4,
					this,
					V1->m_pEndClosestNode,
					pIntPoint);
	}
}

// find the cost to the nearst soma and the best way to connect with it
int CVessel::FindSomaDistance(TopEndMiddleNeither dirFlag, CPoint& FoundPoint)
{
	int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;
	int SomaID = 0;
	int pixelValue = 0;
	int PathSum = 0;
	int BestLength = giMaxExtensionLength;
	CPoint* Indices = NULL;
	int newX, newY, newZ;
	int Hindex = 1;
	int Vindex = 1;
	int Hdirections[17];
	int Vdirections[17]; // NumOfDirections/2 + 1
	CPoint* FromPoint = NULL;

	if (dirFlag == OnTop)
		FromPoint = m_Center.head->data;
	else
		FromPoint = m_Center.tail->data;

	int Hdir = FromPoint->m_iHDir;
	int Vdir = FromPoint->m_iVDir;

	Hdirections[0] = Hdir;
	Vdirections[0] = Vdir;

	register int i, j, k;

	for (i = 1; i < 9; i++)
	{
		Hdirections[Hindex] = DirectionPlus(Hdir, i);
		Hindex++;
		Hdirections[Hindex] = DirectionMinus(Hdir, i);
		Hindex++;

		Vdirections[Vindex] = DirectionPlus(Vdir, i);
		Vindex++;
		Vdirections[Vindex] = DirectionMinus(Vdir, i);
		Vindex++;
	}

	for (i = 0; i < 17; i++)
	{
		for (j = 0; j < 17; j++)
		{
			Indices = gVectorsArray[Hdirections[i]][Vdirections[j]]->m_pIndices;

			// the cost of this path is given by the sum of its gray level values
			// over the original image
			PathSum = 0;
			// we loop until the maximum vector length because in case of
			// wide vessels the centerline of such vessels might be more
			// than LOOK_AHEAD pixels away. 
			for (k = 1; k < giMaxExtensionLength; k++)
			{
				newX = FromPoint->m_iX + Indices[k].m_iX;
				newY = FromPoint->m_iY + Indices[k].m_iY;
				newZ = FromPoint->m_iZ + Indices[k].m_iZ;

				if (newX < 0 ||
					newX >= iCols ||
					newY < 0 ||
					newY >= iRows ||
					newZ < 0 ||
					newZ >= iSlices)
					break;

				pixelValue = Empty3DImage->data[newZ][newY][newX];

				// do not step over an other vessel to reach the soma
				if (pixelValue && pixelValue < StartingSomaColor)
					break;

				PathSum += The3DImage->data[newZ][newY][newX];

				// keep going until you hit a soma pixel
				if (pixelValue > StartingSomaColor)
				{
					// Consider this intersection only if it is better than the optimal
					if (k < BestLength)
					{
						BestLength = k;
						FoundPoint.m_iX = newX; 
						FoundPoint.m_iY = newY;
						FoundPoint.m_iZ = newZ;
						FoundPoint.m_iHDir = static_cast<unsigned char>(Hdirections[i]);
						FoundPoint.m_iVDir = static_cast<unsigned char>(Vdirections[j]);
						SomaID = pixelValue - StartingSomaColor;
						break;
					}
					else
						break;
				}
			}
		}
	}

	return SomaID;
}

// find the cost to the nearst soma and the best way to connect with it
int CVessel::FindSomaDistance2(TopEndMiddleNeither dirFlag, CPoint& FoundPoint)
{
	//int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;
	CPoint* FromPoint = NULL;

	if (dirFlag == OnTop)
		FromPoint = m_Center.head->data;
	else
		FromPoint = m_Center.tail->data;

	register int i;
	register int x = FromPoint->m_iX;
	register int y = FromPoint->m_iY;
	register int z = FromPoint->m_iZ;
	register int dx, dy, dz;
	register int distance;
	register int minDistance = iCols * iRows;

	register int result = 0;

	for (i = 0; i < giSomaVolume; i++)
	{
		dx = gaSomaPoints[i].m_iX - x;
		dy = gaSomaPoints[i].m_iY - y;
		dz = gaSomaPoints[i].m_iZ - z;
		if ((ABS(dx) + ABS(dy) + ABS(dz)) < 70)
		{
			distance = static_cast<int>(std::sqrt( static_cast<double>(dx * dx + dy * dy + dz * dz) )) ;
			if (distance <= 20 && distance < minDistance)
			{
				minDistance = distance;
				FoundPoint = gaSomaPoints[i];
				result = 1;
			}
		}
	}
	return result;
}

// find the cost to the nearst soma and the best way to connect with it
int CVessel::FindSomaDistance2(TopEndMiddleNeither dirFlag,
	CPoint& FoundPoint, int& FoundSomaID)
{
	//int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;
	CPoint* FromPoint = NULL;

	if (dirFlag == OnTop)
		FromPoint = m_Center.head->data;
	else
		FromPoint = m_Center.tail->data;

	register int i;
	register int x = FromPoint->m_iX;
	register int y = FromPoint->m_iY;
	register int z = FromPoint->m_iZ;
	register int dx, dy, dz;
	register int distance;
	register int minDistance = iCols * iRows;

	register int result = 0;

	for (i = 0; i < giSomaVolume; i++)
	{
		dx = gaSomaPoints[i].m_iX - x;
		dy = gaSomaPoints[i].m_iY - y;
		dz = gaSomaPoints[i].m_iZ - z;
		if ((ABS(dx) + ABS(dy) + ABS(dz)) < 70)
		{
			distance = static_cast<int>( std::sqrt(static_cast<double>(dx * dx + dy * dy + dz * dz)) ) ;
			if (distance <= 20 && distance < minDistance)
			{
				minDistance = distance;
				FoundPoint = gaSomaPoints[i];
				FoundSomaID = 1;				// For now this function assumes
				result = 1;							// only one soma
			}
		}
	}
	return result;
}


// find the cost to the nearst soma and the best way to connect with it
int CVessel::FindSomaDistance3(TopEndMiddleNeither dirFlag,
	CPoint& FoundPoint, int& FoundSomaID)
{
	//int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;
	CPoint* FromPoint = NULL;

	if (dirFlag == OnTop)
		FromPoint = m_Center.head->data;
	else
		FromPoint = m_Center.tail->data;

	register int i;
	register int x = FromPoint->m_iX;
	register int y = FromPoint->m_iY;
	register int z = FromPoint->m_iZ;
	register int dx, dy, dz;
	register int distance;
	register int minDistance = iCols * iRows;

	for (i = 0; i < giSomaVolume; i++)
	{
		dx = gaSomaPoints[i].m_iX - x;
		dy = gaSomaPoints[i].m_iY - y;
		dz = gaSomaPoints[i].m_iZ - z;
		if ((ABS(dx) + ABS(dy) + ABS(dz)) < 70)
		{
			distance = static_cast<int>(std::sqrt( static_cast<double>(dx * dx + dy * dy + dz * dz) )) ;
			if (distance < minDistance)
			{
				minDistance = distance;
				FoundPoint = gaSomaPoints[i];
				FoundSomaID = 1;				// For now this function assumes
				//Yousef: This has been changed
				//FoundSomaID=  SomaLabelsImage->data[gaSomaPoints[i].m_iZ][gaSomaPoints[i].m_iY][gaSomaPoints[i].m_iX];
			}
		}
	}
	return minDistance;
}

///////////////////////////////////////////////////////////////////////////////
// METHOD: Extendable
//
//  return 1 If the vessel can be extended to intersect with other vessels, and 
//  0 otherwise. A vessel is extendable if it is not already intersecting with another vessel
//  and some vessel lies in the neighborhood of the given end point.
//  Make sure that the pixel does not belong to a soma
int CVessel::Extendable(TopEndMiddleNeither DirFlag)
{
	int iSlices = The3DImage->m_iSlices;
	//int iRows = The3DImage->m_iRows;
	//int iCols = The3DImage->m_iCols;
	int result = 1;
	CPoint* aPoint = NULL;
	int i;
					
					// if this vessel already belongs to an intersection point at this end
					// then it is not extendable
					for(i = 0; i < m_iNumOfIntersectionPoints; i++)
	{
		if (gIntersectionPoints.m_apData[m_aiMyIntersectionPoints[i] - 1]->GetLocation(m_iID) ==
			DirFlag)
			return 0;
	}

	if (DirFlag == OnTop)
		aPoint = m_Center.head->data;
	else
		aPoint = m_Center.tail->data;
	if (aPoint)
	{
		int x = aPoint->m_iX;
		int y = aPoint->m_iY;
		int z = aPoint->m_iZ;

		register int i, j, k;
		int pixelValue;
		// check that none of the 8-neighboring pixels belong to other vessels
		for (i = -1; i <= 1; i++)
		{
			for (j = -1; j <= 1; j++)
			{
				for (k = -1; k <= 1; k++)
				{
					if (z + k <0 || z + k> iSlices)
						continue;

					pixelValue = Empty3DImage->data[z + k][y + i][x + j];
					// if there exist a point belonging to another vessel then add
					// the point to the list of intersection points for each of the vessels.
					// In addition, add such point to the list of my intersection points.
					if (pixelValue &&
						pixelValue != m_iID &&
						pixelValue < StartingIntersectionColor)
					{
						CPoint tempPoint(x + j,
							y + i,
							z + k,
							aPoint->m_iHDir,
							aPoint->m_iVDir);
						CIntersectionPoint IntPoint(tempPoint);
						IntPoint.AddVessel(m_iID, DirFlag);
						// find the location of the point with respect to the vessel
						// Find the location of the point with respect to the vessel
						TopEndMiddleNeither location = gTheVessels.m_apData[pixelValue - 1]->GetPositionOfPoint(&tempPoint);
						if (location != Neither)
							IntPoint.AddVessel(pixelValue, location);
						else
						{
							cout << "CVessel::Extendable.. Error2" << endl;
							exit(0);
						}

						return 0;
					}
					else if (pixelValue > StartingSomaColor)
					{
						// the point belongs to a soma, add it to
						// the list of intersection points for the vessel and the global 
						// array of intersection points.

						CPoint tempPoint(x + j,
							y + i,
							z + k,
							aPoint->m_iHDir,
							aPoint->m_iVDir);
						CIntersectionPoint IntPoint(tempPoint);
						IntPoint.AddVessel(m_iID, DirFlag);
						IntPoint.AddSoma(pixelValue - StartingSomaColor);
						IntPoint.m_Type = SOMA;
						int aPointID = gIntersectionPoints.Add(IntPoint);
						AddIntersectionPoint(aPointID);
						gTheSomas->m_aData[pixelValue - 1].AddIntersectionPoint(aPointID);

						return 0;
					}
				}
			}
		}
		/*
										int dir1 = aPoint->dir;
										int dir2 = DirectionMinus(aPoint->dir, 1);
										int dir3 = DirectionPlus(aPoint->dir, 1);
										unsigned char *aData = &(CenterlinesImage->data[aPoint->y][aPoint->x]);
										// to prevent a vessel from being wrapped around itself
										for(i = 1; i < 4; i++)
										{
											if( (*(aData + VectorsArray[dir1]->data[i]) == ID ) ||
												 (*(aData + VectorsArray[dir2]->data[i]) == ID ) ||
												 (*(aData + VectorsArray[dir3]->data[i]) == ID )   )
											{
												result = 0;
												break;
											}
										}
										*/
	}
	else
		result = 0;

	return result;
}

// return if the vessel is ankered by an intersection point ont
// the give edge
int CVessel::Float(TopEndMiddleNeither DirFlag)
{
	int result = 1;
	int id;
	for (register int i = 0; i < m_iNumOfIntersectionPoints; i++)
	{
		id = m_aiMyIntersectionPoints[i];
		if (gIntersectionPoints.m_apData[id - 1]->GetLocation(m_iID) ==	DirFlag)
		{
			result = 0;
			break;
		}
	}
	return result;
}
/////////////////////////////////////////////////////////////////////////////////
// METHOD: FindIntersectionPoint
//
// Purpose: 
//  To connect one of this vessel's ends with other intersecting vessels. See
//  ConnectWithOtherVessels(). 
//  with any.
//
// Logic:
//	Check the regions surrounding the end of the vessel to see if the vessel
//	intersect with other vessels. Only those vessels lying along a line in 
//	in the direction of my vessel, and within LOOK_AHEAD distance from it are 
//	considered intersecting vessels.  
int CVessel::FindIntersectionPoint(TopEndMiddleNeither DirFlag,
	CPoint& FoundPoint, int& FoundID, int& Cost)
{
	int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;
	int result = 0;
	int iFromLength = gConfig.GetMinimumTemplateLength();
	// this is the cost of the most expensive vessel extension. A vessel extension
	// must be less than this for it to be valid
	int MaxExtensionLength = iFromLength + 2;

	int bestDistance = MaxExtensionLength;
	int pathSum = 0;
	int optimalSum = 0;

	int Hdirections[17];
	int Vdirections[17];

	CPoint* FromPoint = NULL;

	if (DirFlag == OnTop)
		FromPoint = m_Center.head->data;
	else
		FromPoint = m_Center.tail->data;

	int Hdir = FromPoint->m_iHDir;
	int Vdir = FromPoint->m_iVDir;
	Hdirections[0] = Hdir;
	Vdirections[0] = Vdir;
	int Hindex = 1;
	int Vindex = 1;

	register int i, j, k;

	for (i = 1; i < 9; i++)
	{
		Hdirections[Hindex] = DirectionPlus(Hdir, i);
		Hindex++;
		Hdirections[Hindex] = DirectionMinus(Hdir, i);
		Hindex++;

		Vdirections[Vindex] = DirectionPlus(Vdir, i);
		Vindex++;
		Vdirections[Vindex] = DirectionMinus(Vdir, i);
		Vindex++;
	}

	CPoint* Indices = NULL;
	int newX, newY, newZ;

	for (i = 0; i < 17; i++)
	{
		for (j = 0; j < 17; j++)
		{
			Indices = gVectorsArray[Hdirections[i]][Vdirections[j]]->m_pIndices;

			// the cost of this path is given by the sum of its gray level values
			// over the original image
			pathSum = 0;
			// we loop until the maximum vector length because in case of
			// wide vessels the centerline of such vessels might be more
			// than LOOK_AHEAD pixels away. 
			for (k = 1; k < MaxExtensionLength; k++)
			{
				newX = FromPoint->m_iX + Indices[k].m_iX;
				newY = FromPoint->m_iY + Indices[k].m_iY;
				newZ = FromPoint->m_iZ + Indices[k].m_iZ;

				if (newX < 0 ||
					newX >= iCols ||
					newY < 0 ||
					newY >= iRows ||
					newZ < 0 ||
					newZ >= iSlices)
					break;

				pathSum += The3DImage->data[newZ][newY][newX];

				int pixelValue = Empty3DImage->data[newZ][newY][newX];

				// at any point if the vessels intersects with itself, break
				// i.e. consider the next direction
				if (pixelValue == m_iID)
					break;

				// keep going until you hit the centerline of the vessel
				// Remember that the centerline image is set to 255 originally
				if (pixelValue && pixelValue < StartingSomaColor)
				{
					int Threshold = (k - 1) * giMedian;
					// Consider this intersection only if it is better than the optimal
					if (pathSum >= Threshold && k < bestDistance)
					{
						optimalSum = pathSum;
						bestDistance = k;
						FoundPoint.m_iX = newX; 
						FoundPoint.m_iY = newY;
						FoundPoint.m_iZ = newZ;
						FoundPoint.m_iHDir = static_cast<unsigned char>(Hdirections[i]);
						FoundPoint.m_iVDir = static_cast<unsigned char>(Vdirections[j]);
						FoundID = pixelValue;
						result = 1;

						break;
					}
					else
						break;
				}
			}
		}
	}

	Cost = optimalSum;

	return result;
}

/////////////////////////////////////////////////////////////////////////////////
// METHOD: CheckAndConnectVesselEnd
//
// Purpose: 
//  To connect one of this vessel's ends with other intersecting vessels. See
//  ConnectWithOtherVessels(). 
//  with any.
//
// Logic:
//	Check the regions surrounding the end of the vessel to see if the vessel
//	intersect with other vessels. Only those vessels lying along a line in 
//	in the direction of my vessel, and within LOOK_AHEAD distance from it are 
//	considered intersecting vessels. 
/* 
void CVessel::CheckAndConnectVesselEnd(CDLList<CPoint> *aList, TopEnd dirFlag)
{
	CVector *tempVector  = 0;
	CPoint *ExtremePoint = 0;
	int intersectionFlag = 0;
	if(dirFlag == OnTop)
		ExtremePoint = aList->head->data;
	else
		ExtremePoint = aList->tail->data;

	if(ExtremePoint)
	{
		CVector *tempVector = VectorsArray[ExtremePoint->dir];
		int x = ExtremePoint->m_iX;
		int y = ExtremePoint->m_iY;
		int z = ExtremePoint->m_iZ;
		for(register int i = 2; i < MaxVectorLength; i++)
		{
			int newX = x + tempVector->IndexC[i];
			int newY = y + tempVector->IndexR[i];

			// make sure that we do not exceed the boundaries of the image
			if(newX >= iCols || newY >= iRows)
				break;

			if(i < 2*LOOK_AHEAD )
			{	
				if(The3DImage->data[z][newY][newX] == 255)
				{
					CPoint tempPoint(newX, newY, z, ExtremePoint->dir);
					ExtendVessel(aList, &tempPoint, dirFlag);
				}
			}
			// Reached LOOK_AHEAD without intersection with a vessel boundary
			// done
			//else
			//	break;
		}
	}
}
 
void CVessel::CheckAndConnectVesselEnd2(CDLList<CPoint> *aList, TopEnd dirFlag)
{
	CVector *tempVector  = 0;
	CPoint *ExtremePoint = 0;
	int intersectionFlag = 0;
	if(dirFlag == OnTop)
		ExtremePoint = aList->head->data;
	else
		ExtremePoint = aList->tail->data;

	if(ExtremePoint)
	{
		int x = ExtremePoint->m_iX;
		int y = ExtremePoint->m_iY;
		int z = ExtremePoint->m_iZ;
		int dir = DirectionPlus(ExtremePoint->dir, NumOfDirections/2);

		CVector *tempVector = VectorsArray[dir];
		for(register int i = 2; i < MaxVectorLength; i++)
		{
			int newX = x + tempVector->IndexC[i];
			int newY = y + tempVector->IndexR[i];

			// make sure that we do not exceed the boundaries of the image
			if(newX >= iCols || newY >= iRows)
				break;

			if(i < 2*LOOK_AHEAD )
			{
				if(The3DImage->data[z][newY][newX] == 255)
				{
					CPoint tempPoint(newX, newY, z, dir);
					ExtendVessel(aList, &tempPoint, dirFlag);
				}
			}
		}
	}
}
*/

/////////////////////////////////////////////////////////////////////////////////
// METHOD: DrawCenterline
//
// Purpose: 
//	To draw the instant vessel into the provided image.
//  The color of the vessel is its own ID; 
//
// Logic:
// Draw the given vessel in the given Image
void CVessel::DrawCenterline(CImage& anImage, unsigned char color)
{
	unsigned char * *data = anImage.data;

	if (color == 0)
		color = static_cast<unsigned char>(m_iID);

	CLNode<CPoint>* tempNode = NULL;

	// centerline
	tempNode = m_Center.head;
	while (tempNode)
	{
		data[tempNode->data->m_iY][tempNode->data->m_iX] = color;
		//anImage.MarkCrosshairXY(tempNode->data, color);
		tempNode = tempNode->after;
	}
}
void CVessel::DrawCenterlineXZ(CImage& anImage, unsigned char  color)
{
	unsigned char * *data = anImage.data;
	CLNode<CPoint>* tempNode = NULL;
	if (color == 0)
		color = static_cast<unsigned char>(m_iID);

	tempNode = m_Center.head;
	while (tempNode)
	{
		data[tempNode->data->m_iZ][tempNode->data->m_iX] = color;
		tempNode = tempNode->after;
	}
}

void CVessel::DrawBoundariesXZ(CImage& anImage, unsigned char  color)
{
	unsigned char * *data = anImage.data;
	CLNode<CPoint>* tempNode = NULL;
	if (color == 0)
		color = static_cast<unsigned char>(m_iID);

	tempNode = m_HRight.head;
	while (tempNode)
	{
		data[tempNode->data->m_iZ][tempNode->data->m_iX] = 1;
		tempNode = tempNode->after;
	}
	tempNode = m_HLeft.head;
	while (tempNode)
	{
		data[tempNode->data->m_iZ][tempNode->data->m_iX] = 2;
		tempNode = tempNode->after;
	}
	tempNode = m_VRight.head;
	while (tempNode)
	{
		data[tempNode->data->m_iZ][tempNode->data->m_iX] = 3;
		tempNode = tempNode->after;
	}
	tempNode = m_VLeft.head;
	while (tempNode)
	{
		data[tempNode->data->m_iZ][tempNode->data->m_iX] = 4;
		tempNode = tempNode->after;
	}
}
void CVessel::DrawCenterlineYZ(CImage& anImage, unsigned char  color)
{
	if (color == 0)
		color = static_cast<unsigned char>(m_iID);

	unsigned char * *data = anImage.data;
	CLNode<CPoint>* tempNode = NULL;
	tempNode = m_Center.head;
	while (tempNode)
	{
		data[tempNode->data->m_iZ][tempNode->data->m_iY] = color; //ID;
		tempNode = tempNode->after;
	}
}

void CVessel::DrawBoundariesYZ(CImage& anImage, unsigned char  color)
{
	if (color == 0)
		color = static_cast<unsigned char>(m_iID);

	unsigned char * *data = anImage.data;
	CLNode<CPoint>* tempNode = NULL;

	tempNode = m_HRight.head;
	while (tempNode)
	{
		data[tempNode->data->m_iZ][tempNode->data->m_iY] = 1; //ID;
		tempNode = tempNode->after;
	}
	tempNode = m_HLeft.head;
	while (tempNode)
	{
		data[tempNode->data->m_iZ][tempNode->data->m_iY] = 2;
		tempNode = tempNode->after;
	}
	tempNode = m_VRight.head;
	while (tempNode)
	{
		data[tempNode->data->m_iZ][tempNode->data->m_iY] = 3;
		tempNode = tempNode->after;
	}
	tempNode = m_VLeft.head;
	while (tempNode)
	{
		data[tempNode->data->m_iZ][tempNode->data->m_iY] = 4;
		tempNode = tempNode->after;
	}
}


void CVessel::DrawBoundaries(CImage& anImage, unsigned char  aColor)
{
	int color = static_cast<unsigned char>(m_iID);
	if (aColor)
		color = aColor;

	unsigned char * *data = anImage.data;

	CLNode<CPoint>* temp = m_HRight.head;
	while (temp)
	{
		data[temp->data->m_iY][temp->data->m_iX] = 1; //ID;
		temp = temp->after;
	}

	temp = m_HLeft.head;
	while (temp)
	{
		data[temp->data->m_iY][temp->data->m_iX] = 2; //ID;
		temp = temp->after;
	}

	temp = m_VRight.head;
	while (temp)
	{
		data[temp->data->m_iY][temp->data->m_iX] = 3; //ID;
		temp = temp->after;
	}

	temp = m_VLeft.head;
	while (temp)
	{
		data[temp->data->m_iY][temp->data->m_iX] = 4; //ID;
		temp = temp->after;
	}
}

void CVessel::DrawBoundariesHXY(CImage& anImage, unsigned char  aColor)
{
	int color = static_cast<unsigned char>(m_iID);
	if (aColor)
		color = aColor;

	unsigned char * *data = anImage.data;

	CLNode<CPoint>* temp = m_HRight.head;
	while (temp)
	{
		// saturate the borderline to the border
		if (temp->data->m_iY >= anImage.m_iRows)
			temp->data->m_iY = anImage.m_iRows - 1;
		if (temp->data->m_iY < 0)
			temp->data->m_iY = 0;
		if (temp->data->m_iX >= anImage.m_iCols)
			temp->data->m_iX = anImage.m_iCols - 1;
		if (temp->data->m_iX < 0)
			temp->data->m_iX = 0;

		data[temp->data->m_iY][temp->data->m_iX] = 1; //ID;
		temp = temp->after;
	}

	temp = m_HLeft.head;
	while (temp)
	{
		if (temp->data->m_iY >= anImage.m_iRows)
			temp->data->m_iY = anImage.m_iRows - 1;
		if (temp->data->m_iY < 0)
			temp->data->m_iY = 0;
		if (temp->data->m_iX >= anImage.m_iCols)
			temp->data->m_iX = anImage.m_iCols - 1;
		if (temp->data->m_iX < 0)
			temp->data->m_iX = 0;
		data[temp->data->m_iY][temp->data->m_iX] = 2; //ID;
		temp = temp->after;
	}
}

void CVessel::DrawBoundariesHXZ(CImage& anImage, unsigned char  aColor)
{
	int color = static_cast<unsigned char>(m_iID);
	if (aColor)
		color = aColor;

	unsigned char * *data = anImage.data;
	CLNode<CPoint>* tempNode = NULL;
	tempNode = m_HRight.head;
	while (tempNode)
	{
		if (tempNode->data->m_iZ >= anImage.m_iRows)
			tempNode->data->m_iZ = anImage.m_iRows - 1;
		if (tempNode->data->m_iZ < 0)
			tempNode->data->m_iZ = 0;
		if (tempNode->data->m_iX >= anImage.m_iCols)
			tempNode->data->m_iX = anImage.m_iCols - 1;
		if (tempNode->data->m_iX < 0)
			tempNode->data->m_iX = 0;
		data[tempNode->data->m_iZ][tempNode->data->m_iX] = 1;
		tempNode = tempNode->after;
	}
	tempNode = m_HLeft.head;
	while (tempNode)
	{
		if (tempNode->data->m_iZ >= anImage.m_iRows)
			tempNode->data->m_iZ = anImage.m_iRows - 1;
		if (tempNode->data->m_iZ < 0)
			tempNode->data->m_iZ = 0;
		if (tempNode->data->m_iX >= anImage.m_iCols)
			tempNode->data->m_iX = anImage.m_iCols - 1;
		if (tempNode->data->m_iX < 0)
			tempNode->data->m_iX = 0;
		data[tempNode->data->m_iZ][tempNode->data->m_iX] = 2;
		tempNode = tempNode->after;
	}
}

void CVessel::DrawBoundariesHYZ(CImage& anImage, unsigned char  aColor)
{
	int color = static_cast<unsigned char>(m_iID);
	if (aColor)
		color = aColor;

	unsigned char * *data = anImage.data;
	CLNode<CPoint>* tempNode = NULL;
	tempNode = m_HRight.head;
	while (tempNode)
	{
		if (tempNode->data->m_iZ >= anImage.m_iRows)
			tempNode->data->m_iZ = anImage.m_iRows - 1;
		if (tempNode->data->m_iZ < 0)
			tempNode->data->m_iZ = 0;
		if (tempNode->data->m_iY >= anImage.m_iCols)
			tempNode->data->m_iY = anImage.m_iCols - 1;
		if (tempNode->data->m_iY < 0)
			tempNode->data->m_iY = 0;
		data[tempNode->data->m_iZ][tempNode->data->m_iY] = 1; //ID;
		tempNode = tempNode->after;
	}
	tempNode = m_HLeft.head;
	while (tempNode)
	{
		if (tempNode->data->m_iZ >= anImage.m_iRows)
			tempNode->data->m_iZ = anImage.m_iRows - 1;
		if (tempNode->data->m_iZ < 0)
			tempNode->data->m_iZ = 0;
		if (tempNode->data->m_iY >= anImage.m_iCols)
			tempNode->data->m_iY = anImage.m_iCols - 1;
		if (tempNode->data->m_iY < 0)
			tempNode->data->m_iY = 0;
		data[tempNode->data->m_iZ][tempNode->data->m_iY] = 2;
		tempNode = tempNode->after;
	}
}

void CVessel::DrawBoundariesVXY(CImage& anImage, unsigned char aColor)
{
	int color = static_cast<unsigned char>(m_iID);
	if (aColor)
		color = aColor;

	unsigned char * *data = anImage.data;
	CLNode<CPoint>* temp = NULL;
	temp = m_VRight.head;
	while (temp)
	{
		// saturate the borderline to the border
		if (temp->data->m_iY >= anImage.m_iRows)
			temp->data->m_iY = anImage.m_iRows - 1;
		if (temp->data->m_iY < 0)
			temp->data->m_iY = 0;
		if (temp->data->m_iX >= anImage.m_iCols)
			temp->data->m_iX = anImage.m_iCols - 1;
		if (temp->data->m_iX < 0)
			temp->data->m_iX = 0;
		data[temp->data->m_iY][temp->data->m_iX] = 1; //ID;
		temp = temp->after;
	}

	temp = m_VLeft.head;
	while (temp)
	{
		// saturate the borderline to the border
		if (temp->data->m_iY >= anImage.m_iRows)
			temp->data->m_iY = anImage.m_iRows - 1;
		if (temp->data->m_iY < 0)
			temp->data->m_iY = 0;
		if (temp->data->m_iX >= anImage.m_iCols)
			temp->data->m_iX = anImage.m_iCols - 1;
		if (temp->data->m_iX < 0)
			temp->data->m_iX = 0;
		data[temp->data->m_iY][temp->data->m_iX] = 2; //ID;
		temp = temp->after;
	}
}

void CVessel::DrawBoundariesVXZ(CImage& anImage, unsigned char  aColor)
{
	int color = static_cast<unsigned char>(m_iID);
	if (aColor)
		color = aColor;

	unsigned char * *data = anImage.data;
	CLNode<CPoint>* tempNode = NULL;
	tempNode = m_VRight.head;
	while (tempNode)
	{
		// saturate the borderline to the border
		if (tempNode->data->m_iZ >= anImage.m_iRows)
			tempNode->data->m_iZ = anImage.m_iRows - 1;
		if (tempNode->data->m_iZ < 0)
			tempNode->data->m_iZ = 0;
		if (tempNode->data->m_iX >= anImage.m_iCols)
			tempNode->data->m_iX = anImage.m_iCols - 1;
		if (tempNode->data->m_iX < 0)
			tempNode->data->m_iX = 0;

		data[tempNode->data->m_iZ][tempNode->data->m_iX] = 1;
		tempNode = tempNode->after;
	}
	tempNode = m_VLeft.head;
	while (tempNode)
	{
		if (tempNode->data->m_iZ >= anImage.m_iRows)
			tempNode->data->m_iZ = anImage.m_iRows - 1;
		if (tempNode->data->m_iZ < 0)
			tempNode->data->m_iZ = 0;
		if (tempNode->data->m_iX >= anImage.m_iCols)
			tempNode->data->m_iX = anImage.m_iCols - 1;
		if (tempNode->data->m_iX < 0)
			tempNode->data->m_iX = 0;
		data[tempNode->data->m_iZ][tempNode->data->m_iX] = 2;
		tempNode = tempNode->after;
	}
}

void CVessel::DrawBoundariesVYZ(CImage& anImage, unsigned char  aColor)
{
	int color = static_cast<unsigned char>(m_iID);
	if (aColor)
		color = aColor;

	unsigned char * *data = anImage.data;
	CLNode<CPoint>* tempNode = NULL;
	tempNode = m_VRight.head;
	while (tempNode)
	{
		if (tempNode->data->m_iZ >= anImage.m_iRows)
			tempNode->data->m_iZ = anImage.m_iRows - 1;
		if (tempNode->data->m_iZ < 0)
			tempNode->data->m_iZ = 0;
		if (tempNode->data->m_iY >= anImage.m_iCols)
			tempNode->data->m_iY = anImage.m_iCols - 1;
		if (tempNode->data->m_iY < 0)
			tempNode->data->m_iY = 0;
		data[tempNode->data->m_iZ][tempNode->data->m_iY] = 1;
		tempNode = tempNode->after;
	}
	tempNode = m_VLeft.head;
	while (tempNode)
	{
		if (tempNode->data->m_iZ >= anImage.m_iRows)
			tempNode->data->m_iZ = anImage.m_iRows - 1;
		if (tempNode->data->m_iZ < 0)
			tempNode->data->m_iZ = 0;
		if (tempNode->data->m_iY >= anImage.m_iCols)
			tempNode->data->m_iY = anImage.m_iCols - 1;
		if (tempNode->data->m_iY < 0)
			tempNode->data->m_iY = 0;
		data[tempNode->data->m_iZ][tempNode->data->m_iY] = 2;
		tempNode = tempNode->after;
	}
}

void CVessel::DrawCenterline(C3DImage& anImage, unsigned char  color)
{
	if (color == 0)
		color = static_cast<unsigned char>(m_iID);

	unsigned char * **data = anImage.data;
	CLNode<CPoint>* temp = m_Center.head;
	while (temp)
	{
		data[temp->data->m_iZ][temp->data->m_iY][temp->data->m_iX] = static_cast<unsigned char>(m_iID);
		temp = temp->after;
	}
}
void CVessel::DrawBoundaries(C3DImage& anImage, unsigned char  aColor)
{
	unsigned char color = static_cast<unsigned char>(m_iID);
	if (aColor)
		color = aColor;

	unsigned char * **data = anImage.data;
	CLNode<CPoint>* temp = m_HLeft.head;
	while (temp)
	{
		data[temp->data->m_iZ][temp->data->m_iY][temp->data->m_iX] = color; //ID;
		temp = temp->after;
	}
	temp = m_HRight.head;
	while (temp)
	{
		data[temp->data->m_iZ][temp->data->m_iY][temp->data->m_iX] = static_cast<unsigned char>(color - 1); //ID;
		temp = temp->after;
	}
	temp = m_VLeft.head;
	while (temp)
	{
		data[temp->data->m_iZ][temp->data->m_iY][temp->data->m_iX] = static_cast<unsigned char>(color - 1); //ID;
		temp = temp->after;
	}

	temp = m_VRight.head;
	while (temp)
	{
		data[temp->data->m_iZ][temp->data->m_iY][temp->data->m_iX] = static_cast<unsigned char>(color - 1); //ID;
		temp = temp->after;
	}
}

void CVessel::DrawThickLine(CImage& , unsigned char )
{
	/*
	if(color == 0)
		color = m_iID;
	CLNode<CPoint> *from = m_Center.head;
	CLNode<CPoint> *to = m_Center.tail;
	while(to != NULL)
	{
		from = to;
		to = to->after;
		if(to == NULL)
		{
			anImage.MarkPoint(from->data, 255);
			break;
		}
		else
			::DrawThickLine(from->data, to->data, anImage, 255);
	}
	*/
}

void CVessel::DrawThickLine(C3DImage& anImage, unsigned char color)
{
	int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;
	if (color == 0)
		color = static_cast<unsigned char>(m_iID);

	//	register int i, j, k;

	unsigned char * **data = anImage.data;
	CLNode<CPoint>* temp = m_Center.head;
	while (temp)
	{
		// since all points are garanteed to be within the margin, we don't need to
		// check that we stepped out of boundaries
		data[temp->data->m_iZ][temp->data->m_iY][temp->data->m_iX] = color;
		if (temp->data->m_iZ > 0)
			data[temp->data->m_iZ - 1][temp->data->m_iY][temp->data->m_iX] = color;
		if (temp->data->m_iZ < iSlices - 1)
			data[temp->data->m_iZ + 1][temp->data->m_iY][temp->data->m_iX] = color;		
		if (temp->data->m_iY > 0)
			data[temp->data->m_iZ][temp->data->m_iY - 1][temp->data->m_iX] = color;
		if (temp->data->m_iY < iRows - 1)
			data[temp->data->m_iZ][temp->data->m_iY + 1][temp->data->m_iX] = color;
		if (temp->data->m_iX > 0)
			data[temp->data->m_iZ][temp->data->m_iY][temp->data->m_iX - 1] = color;
		if (temp->data->m_iX < iCols)
			data[temp->data->m_iZ][temp->data->m_iY][temp->data->m_iX + 1] = color;

		temp = temp->after;
	}
}




void CVessel::Print(ostream& out)
{
	int x1 = m_Center.head->data->m_iX; 
	int x2 = m_Center.tail->data->m_iX;
	int y1 = m_Center.head->data->m_iY;
	int y2 = m_Center.tail->data->m_iY;
	int z1 = m_Center.head->data->m_iZ;
	int z2 = m_Center.tail->data->m_iZ;

	out << setw(3) << m_iID << "\t" << setw(3) << x1 << "," << setw(3) << y1
		<< "," << setw(3) << z1 << "\t" << setw(3) << x2 << "," << setw(3)
		<< y2 << "," << setw(3) << z2 << "\t" << setw(3) << m_Center.length
		<< "\t" << setprecision(3) << m_fHWidth << "\t" << setprecision(3)
		<< m_fVWidth << endl;
}


void CVessel::MergeVessels()
{
	/*
	int DistanceThreshold = (1+giLOOK_AHEAD)*(1+giLOOK_AHEAD);
	int DirectionThr = 3;
	// for each of my intersection points
	for(register int i = 0; i < m_m_iNumOfIntersectionPoints; i++)
	{
		int PointID = m_aim_aiMyIntersectionPoints[i]->PointID;
		// if the point lies on one of the ends
		if(EndPoint(gTheVessels.m_aim_aiMyIntersectionPoints[PointID]->thePoint))
		{
			// for each (different) intersection point from "TheVessels" 
			// that is closer than threshold
			for(register int j = 0; j < gTheVessels.m_m_iNumOfIntersectionPoints; j++)
			{
			   if(j != PointID && gTheVessels.m_apDistanceArray[PointID][j] < DistanceThreshold)
				{
					CPoint *tempPoint = gTheVessels.m_aiMyIntersectionPoints[j]->thePoint;
					for(register int k = 0; k < gTheVessels.m_iNumOfElements; k++)
					{
						// skip this vessel
						if(k+1 == ID)
							continue;
						// if such intersection point is also an end point for
						// this vessel, 
						if(gTheVessels.m_apData[k]->EndPoint(tempPoint))
						{
							// if such vessel and "this" are of similar orientation
							// then merge them.
							int dirDist = DirDistance(tempPoint->dir,
								gTheVessels.m_aiMyIntersectionPoints[PointID]->thePoint->dir);
							int dirDist2 = DirDistance(DirectionPlus(tempPoint->dir, NumOfDirections/2),
								gTheVessels.m_aiMyIntersectionPoints[PointID]->thePoint->dir);
							
							if(dirDist < DirectionThr || dirDist2 < DirectionThr)
							{
								cout << "\n\n\n The Two vessels: " << ID << ", and "
									<< k+1 << " are mergable \n" << endl;
							}
						}
					}
				}
			}
		}
	}
	*/
}
///////////////////////////////////////////////////////////////////////////
void CVessel::ExtendByAnotherVessel(CVessel* rhs, TopEndMiddleNeither pos1,
	TopEndMiddleNeither pos2)
{
	if (pos1 == OnTop && pos2 == OnTop)
	{
		m_Center.ExtendListOnTopFromTop(rhs->m_Center);
		m_iLength += rhs->m_iLength;
	}
	else if (pos1 == OnTop && pos2 == OnEnd)
	{
		m_Center.ExtendListOnTopFromEnd(rhs->m_Center);
		m_iLength += rhs->m_iLength;
	}
	else if (pos1 == OnEnd && pos2 == OnTop)
	{
		m_Center.ExtendListOnEndFromTop(rhs->m_Center);
		m_iLength += rhs->m_iLength;
	}
	else if (pos1 == OnEnd && pos2 == OnEnd)
	{
		m_Center.ExtendListOnEndFromEnd(rhs->m_Center);
		m_iLength += rhs->m_iLength;
	}
	else
	{
		// this should never fire
		cout << "CVessel::ExtendByAnotherVessel.... ERROR" << endl;
		exit(0);
	}
}
//////////////////////////////////////////////////////////////////////////////
// METHOD: WriteID
//
// Purpose: Write the vessels ID at each end of the vessel. Ignore all vessels
// that are shorter than 20 pixels long. For all vessels that are > 20 pixels
// but are < 40 pixels only write an ID on one end of the vessel where it does
// have an intersection point
void CVessel::WriteIDXY(CImage& anImage, unsigned char color)
{
	if (m_Center.length > LengthThreshold)
	{
		CPoint tempPoint;

		FindLocationToWriteIDXY(tempPoint);
		WriteID(anImage, tempPoint, color);
	}
}
void CVessel::WriteIDXZ(CImage& anImage, unsigned char color)
{
	if (m_Center.length > LengthThreshold)
	{
		CPoint tempPoint;

		FindLocationToWriteIDXZ(tempPoint);
		WriteID(anImage, tempPoint, color);
	}
}
void CVessel::WriteIDYZ(CImage& anImage, unsigned char color)
{
	if (m_Center.length > LengthThreshold)
	{
		CPoint tempPoint;

		FindLocationToWriteIDYZ(tempPoint);
		WriteID(anImage, tempPoint, color);
	}
}

//////////////////////////////////////////////////////////////////////////////
// METHOD: WriteID
//
void CVessel::WriteID(CImage& anImage, CPoint& aPoint, unsigned char color)
{
	//int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;
	int HFrom = 0 - DigitWidth / 2;  // 0 - 9/2 = -4
	int HTo = DigitWidth / 2;      // 4;
	int VFrom = 0 - DigitHight / 2;  // -4
	int VTo = DigitHight / 2;      // 2;
	int X = aPoint.m_iX;
	int Y = aPoint.m_iY;
	unsigned char * *data = anImage.data;

	int Digits[2] =
	{
		-1, -1
	};

	if (m_iID < 10)
		Digits[0] = m_iID;
	else
	{
		Digits[0] = m_iID / 10;
		Digits[1] = m_iID % 10;
	}

	for (int d = 0; d < 2; d++)
	{
		if (Digits[d] == -1)
			break;
		int Hindex = 0;
		int Vindex = 0;
		if (d == 1)
			X += DigitWidth - 2;

		for (register int i = VFrom; i <= VTo; i++)
		{
			Hindex = 0;
			for (register int j = HFrom; j <= HTo; j++)
			{
				int Xprime = j + X;
				int Yprime = i + Y;

				if (Xprime >= 0 &&
					Xprime < iCols &&
					Yprime >= 0 &&
					Yprime < iRows)
				{
					if (Numbers[Digits[d]][Vindex][Hindex])
						data[Yprime][Xprime] = color;
				}
				//	else
				//		data[i+Y][j+X] = 0;

				Hindex++;
			}
			Vindex++;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// METHOD: FindLocationToWriteID
//
int CVessel::FindLocationToWriteIDXY(CPoint& atPoint)
{
	//int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;
	int HFrom = 0 - (DigitWidth / 2 + 1);  // 0 - 9/2 = -4 - 1 = -5
	int HTo = DigitWidth / 2 + 1;      // 4 + 1;
	int VFrom = 0 - (DigitHight / 2 + 1);  // -4 -1
	int VTo = DigitHight / 2 + 1;      // 2;

	// if the vessels ID is a two digit number, then we should cosider more area
	if (m_iID > 9)
		HTo += (DigitWidth - 2);
	int Sum = 0;
	int MinSum = 10000;
	int margin = 10;
	int pixelValue = 0;


	CLNode<CPoint>* pNodeFrom = m_Center.head;
	CLNode<CPoint>* pNodeTo = NULL;
	int iCounter = 0;
	int iFromCount = (m_Center.length / 2) - 5;
	int iToCount = (m_Center.length / 2) + 5;

	while (iCounter < iFromCount)
	{
		if (pNodeFrom == NULL)
			return 0;

		iCounter++;
		pNodeFrom = pNodeFrom->after;
	}

	pNodeTo = pNodeFrom;
	while (iCounter < iToCount)
	{
		if (pNodeTo == NULL)
			return 0;

		iCounter++;
		pNodeTo = pNodeTo->after;
	}

	CLNode<CPoint>* pNode = pNodeFrom;

	register int i, j, k;
	int X, Y;
	int iIncrement = - 1;
	iCounter = 5;

	// for each possible node on the center line of the 
	// vessel, select the best one
	while (pNode != pNodeTo)
	{
		iCounter += iIncrement;
		if (iCounter == 0)
			iIncrement = 1;

		for (i = 1; i < NumOfDirections; i++)
		{
			X = pNode->data->m_iX + gVectorsArray[i][0]->m_pIndices[5].m_iX;
			Y = pNode->data->m_iY + gVectorsArray[i][0]->m_pIndices[5].m_iY;

			if (Y < margin ||
				Y >= iRows - margin ||
				X < margin ||
				X >= iCols - margin)
				continue;

			Sum = 0;

			// count how many other vessel's pixels are around this point
			// and select the point with the fewest such points
			for (j = VFrom; j <= VTo; j++)
			{
				for (k = HFrom; k <= HTo; k++)
				{
					pixelValue = TrackImageXY->data[Y + j][X + k];

					// penalty for writing over other numbers
					if (pixelValue == IDColor)
						Sum += 20;

					// prefer darker areas for writing the number
					Sum += (pixelValue / 50);
				}
			}

			Sum += iCounter;
			if (Sum < MinSum)
			{
				MinSum = Sum;
				atPoint.m_iX = X;
				atPoint.m_iY = Y;
			}
		}

		pNode = pNode->after;
	}

	return MinSum;
}
int CVessel::FindLocationToWriteIDXZ(CPoint& atPoint)
{
	int iSlices = The3DImage->m_iSlices;
	//int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;
	int HFrom = 0 - (DigitWidth / 2 + 1);  // 0 - 9/2 = -4 - 1 = -5
	int HTo = DigitWidth / 2 + 1;      // 4 + 1;
	int VFrom = 0 - (DigitHight / 2 + 1);  // -4 -1
	int VTo = DigitHight / 2 + 1;      // 2;

	// if the vessels ID is a two digit number, then we should cosider more area
	if (m_iID > 9)
		HTo += (DigitWidth - 2);
	int Sum = 0;
	int MinSum = 10000;
	int margin = 10;
	int pixelValue = 0;


	CLNode<CPoint>* pNodeFrom = m_Center.head;
	CLNode<CPoint>* pNodeTo = NULL;
	int iCounter = 0;
	int iFromCount = (m_Center.length / 2) - 5;
	int iToCount = (m_Center.length / 2) + 5;

	while (iCounter < iFromCount)
	{
		if (pNodeFrom == NULL)
			return 0;

		iCounter++;
		pNodeFrom = pNodeFrom->after;
	}

	pNodeTo = pNodeFrom;
	while (iCounter < iToCount)
	{
		if (pNodeTo == NULL)
			return 0;

		iCounter++;
		pNodeTo = pNodeTo->after;
	}

	CLNode<CPoint>* pNode = pNodeFrom;

	register int i, j, k;
	int X, Z;
	int iIncrement = - 1;
	iCounter = 5;

	// for each possible node on the center line of the 
	// vessel, select the best one
	while (pNode != pNodeTo)
	{
		iCounter += iIncrement;
		if (iCounter == 0)
			iIncrement = 1;

		for (i = 1; i < NumOfDirections; i++)
		{
			X = pNode->data->m_iX + gVectorsArray[0][i]->m_pIndices[5].m_iX;
			Z = pNode->data->m_iZ + gVectorsArray[0][i]->m_pIndices[5].m_iZ;

			if (Z < margin ||
				Z >= iSlices - margin ||
				X < margin ||
				X >= iCols - margin)
				continue;

			Sum = 0;

			// count how many other vessel's pixels are around this point
			// and select the point with the fewest such points
			for (j = VFrom; j <= VTo; j++)
			{
				for (k = HFrom; k <= HTo; k++)
				{
					pixelValue = TrackImageXZ->data[Z + j][X + k];

					// penalty for writing over other numbers
					if (pixelValue == IDColor)
						Sum += 20;

					// prefer darker areas for writing the number
					Sum += (pixelValue / 50);
				}
			}

			Sum += iCounter;
			if (Sum < MinSum)
			{
				MinSum = Sum;
				atPoint.m_iX = X;
				atPoint.m_iY = Z;
			}
		}

		pNode = pNode->after;
	}

	return MinSum;
}
int CVessel::FindLocationToWriteIDYZ(CPoint& atPoint)
{
	int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	//int iCols = The3DImage->m_iCols;
	int HFrom = 0 - (DigitWidth / 2 + 1);  // 0 - 9/2 = -4 - 1 = -5
	int HTo = DigitWidth / 2 + 1;      // 4 + 1;
	int VFrom = 0 - (DigitHight / 2 + 1);  // -4 -1
	int VTo = DigitHight / 2 + 1;      // 2;

	// if the vessels ID is a two digit number, then we should cosider more area
	if (m_iID > 9)
		HTo += (DigitWidth - 2);
	int Sum = 0;
	int MinSum = 10000;
	int margin = 10;
	int pixelValue = 0;


	CLNode<CPoint>* pNodeFrom = m_Center.head;
	CLNode<CPoint>* pNodeTo = NULL;
	int iCounter = 0;
	int iFromCount = (m_Center.length / 2) - 5;
	int iToCount = (m_Center.length / 2) + 5;

	while (iCounter < iFromCount)
	{
		if (pNodeFrom == NULL)
			return 0;

		iCounter++;
		pNodeFrom = pNodeFrom->after;
	}

	pNodeTo = pNodeFrom;
	while (iCounter < iToCount)
	{
		if (pNodeTo == NULL)
			return 0;

		iCounter++;
		pNodeTo = pNodeTo->after;
	}

	CLNode<CPoint>* pNode = pNodeFrom;

	register int i, j, k;
	int Z, Y;
	int iIncrement = - 1;
	iCounter = 5;

	// for each possible node on the center line of the 
	// vessel, select the best one
	while (pNode != pNodeTo)
	{
		iCounter += iIncrement;
		if (iCounter == 0)
			iIncrement = 1;

		for (i = 1; i < NumOfDirections; i++)
		{
			Y = pNode->data->m_iY + gVectorsArray[0][i]->m_pIndices[5].m_iY;
			Z = pNode->data->m_iZ + gVectorsArray[0][i]->m_pIndices[5].m_iZ;

			if (Y < margin ||
				Y >= iRows - margin ||
				Z < margin ||
				Z >= iSlices - margin)
				continue;

			Sum = 0;

			// count how many other vessel's pixels are around this point
			// and select the point with the fewest such points
			for (j = VFrom; j <= VTo; j++)
			{
				for (k = HFrom; k <= HTo; k++)
				{
					pixelValue = TrackImageYZ->data[Z + j][Y + k];

					// penalty for writing over other numbers
					if (pixelValue == IDColor)
						Sum += 20;

					// prefer darker areas for writing the number
					Sum += (pixelValue / 50);
				}
			}

			Sum += iCounter;
			if (Sum < MinSum)
			{
				MinSum = Sum;
				atPoint.m_iX = Y;
				atPoint.m_iY = Z;
			}
		}

		pNode = pNode->after;
	}

	return MinSum;
}

//////////////////////////////////////////////////////////////////////
// METHOD: 	RemoveIntersectionPoint
//
// Remove the intersection point with the given ID from my array
/*
void CVessel::RemoveIntersectionPoint(int pointID)
{
	int index = -1;
	for(register int i = 0; i < m_iNumOfIntersectionPoints; i++)
	{
		if(m_aiMyIntersectionPoints[i]->PointID == pointID)
		{
			index = i;
			break;
		}
	}
	for(i = index; i < m_iNumOfIntersectionPoints - 1; i++)
		m_aiMyIntersectionPoints[i] = m_aiMyIntersectionPoints[i+1];

	m_iNumOfIntersectionPoints--;
	
}
*/
///////////////////////////////////////////////////////////////////////////
// Method: MarkInImage
//
// Mark the vessel in the given image. For each original centerline point 
// draw a ball. The radius of the ball is determined by the corresponding
// boundary points.
void CVessel::MarkInImage(C3DImage& anImage, unsigned char color)
{
	CLNode<CPoint>* C = m_Center.head;
	CLNode<CPoint>* HL = m_HLeft.head;
	CLNode<CPoint>* HR = m_HRight.head;
	CLNode<CPoint>* VL = m_VLeft.head;
	CLNode<CPoint>* VR = m_VRight.head;

	int HX, HY, HZ, VX, VY, VZ;
	int radius, radius2;
	while (C)
	{
		HX = (HL->data->m_iX - HR->data->m_iX);
		HX = HX * HX;
		HY = (HL->data->m_iY - HR->data->m_iY);
		HY = HY * HY;
		HZ = (HL->data->m_iZ - HR->data->m_iZ);
		HZ = HZ * HZ;

		VX = (VL->data->m_iX - VR->data->m_iX);
		VX = VX * VX;
		VY = (VL->data->m_iY - VR->data->m_iY);
		VY = VY * VY;
		VZ = (VL->data->m_iZ - VR->data->m_iZ);
		VZ = VZ * VZ;
		radius = static_cast<int>(std::sqrt( static_cast<double>(HX + HY + HZ) ));
		radius2 = static_cast<int>(std::sqrt( static_cast<double>(VX + VY + VZ) ));
		if (radius2 < radius)
			radius = radius2;

		if (radius > 5)
			radius = 5;

		DrawBall(anImage, * (C->data), radius, color);

		// get the next points
		C = C->after;
		HL = HL->after;
		HR = HR->after;
		VL = VL->after;
		VR = VR->after;

		// get the next points
		while (C && C->data->m_lUserFlag == 0)
			C = C->after;
		while (HL && HL->data->m_lUserFlag == 0)
			HL = HL->after;
		while (HR && HR->data->m_lUserFlag == 0)
			HR = HR->after;
		while (VL && VL->data->m_lUserFlag == 0)
			VL = VL->after;
		while (VR && VR->data->m_lUserFlag == 0)
			VR = VR->after;

		if (! HL || ! HR || ! VL || ! VR)
			break;
	}
}

/////////////////////////////////// 
// Method:GetPositionOfFloatingEnd
//
// If the vessel has a floating end, return a pointer to its position
//
CPoint* CVessel::GetPositionOfFloatingEnd()
{
	CPoint* result = NULL;

	if (m_iNumOfIntersectionPoints == 1)
	{
		CIntersectionPoint* temp = gIntersectionPoints.GetPoint(m_aiMyIntersectionPoints[0]);
		if (temp)
		{
			if (*m_Center.head->data == temp->m_Point)
				result = m_Center.tail->data;
			else if (*m_Center.tail->data == temp->m_Point)
				result = m_Center.head->data;
		}
	}

	return result;
}

/////////////////////////////////////////
// Method: GetIdOfOtherIntersectionPoint
//
// return the id of the other intersection point.  After breaking on intersection
// points, each vessel can have a maximum of two such points at each end. This
// function returns the id of the point with the id different from the
// argument
int CVessel::GetIdOfOtherIntersectionPoint(int id)
{
	int result = - 1;
	// if I have more than two intersection points, error
	if (m_iNumOfIntersectionPoints > 2)
	{
		cout << "\n CVessel::GetIdOfOtherIntersectionPoint. Error1" << endl;	
		//exit(0);
		return result;
	}
	if (m_iNumOfIntersectionPoints == 2)
	{
		if (m_aiMyIntersectionPoints[0] == m_aiMyIntersectionPoints[1])
		{
			cout << "\n CVessel::GetIdOfOtherIntersectionPoint. Error2"
				<< endl;
			//exit(0);
			return result;
		}
		if (m_aiMyIntersectionPoints[0] == id)
			result = m_aiMyIntersectionPoints[1];
		else
			result = m_aiMyIntersectionPoints[0];
	}
	return result;
}

// calculate the vessel's width at the given point. Notice that
// if the point does not lie on the vessel, the function returns -1
int CVessel::GetWidthAtPoint(CPoint* pPoint)
{
	int result = - 1;
	CLNode<CPoint>* pcPoint = m_Center.head;
	CLNode<CPoint>* pcPoint2 = NULL;
	while (pcPoint)
	{
		if (*pcPoint->data == *pPoint)
		{
			// if the found point is original, then it has the width
			// stored int "m_lUserFlag", else find the closes 
			// original point, and use its width
			if (pcPoint->data->m_lUserFlag)
				result = pcPoint->data->m_lUserFlag;
			else
			{
				pcPoint2 = pcPoint->after;
				while (pcPoint2)
				{
					if (pcPoint2->data->m_lUserFlag)
					{
						result = pcPoint2->data->m_lUserFlag;
						break;
					}
					pcPoint2 = pcPoint2->after;
				}
				// not found in the after direction, try the before direction
				if (pcPoint2 == NULL)
				{
					pcPoint2 = pcPoint->before;

					while (pcPoint2)
					{
						if (pcPoint2->data->m_lUserFlag)
						{
							result = pcPoint2->data->m_lUserFlag;
							break;
						}
						pcPoint2 = pcPoint2->before;
					}
				}
			}

			break;
		}
		pcPoint = pcPoint->after;
	}

	return result;
}

/////////////////////////////////////////
// Method: GetPositionOfPoint
//
// Determine whether the given point is onTop, onEnd, in the Middle, or not
// on the vessel
TopEndMiddleNeither CVessel::GetPositionOfPoint(CPoint* aPoint)
{
	if (! aPoint)
		return Neither;

	if (*aPoint == *(m_Center.head->data))
		return OnTop;

	if (*aPoint == *(m_Center.tail->data))
		return OnEnd;

	CLNode<CPoint>* tempNode = m_Center.head;
	while (tempNode)
	{
		if (*(tempNode->data) == *aPoint)
			return Middle;

		tempNode = tempNode->after;
	}

	return Neither;
}

//////////////////////////////////
// Method: AddIntersectionPoint
//
// Add the given intersection point ID to my collection
void CVessel::AddIntersectionPoint(int pointID)
{
	// if the point is already there, don't add it
	for (register int i = 0; i < m_iNumOfIntersectionPoints; i++)
		if (m_aiMyIntersectionPoints[i] == pointID)
			return;

	if (m_iNumOfIntersectionPoints > m_iSizeOfIntersectionPointsArray)
	{
		int* tempPtr = m_aiMyIntersectionPoints;
		m_aiMyIntersectionPoints = new int [m_iNumOfIntersectionPoints + BlockSize];
		memcpy(m_aiMyIntersectionPoints,
			tempPtr,
			sizeof(int) * m_iNumOfIntersectionPoints);

		m_aiMyIntersectionPoints[m_iNumOfIntersectionPoints] = pointID;
		m_iNumOfIntersectionPoints++;
		m_iSizeOfIntersectionPointsArray += BlockSize;

		if (tempPtr)
			delete [] tempPtr;
	}
	else
	{
		m_aiMyIntersectionPoints[m_iNumOfIntersectionPoints] = pointID;
		m_iNumOfIntersectionPoints++;
	}
}
// see below
void CVessel::BreakOnIntersectionPoint(CPoint* pPoint)
{
	if (EndPoint(pPoint))
		return;

	// find the break points for each of the linked lists in this vessel
	CLNode<CPoint>* centerBreakNode = FindClosestPoint(&m_Center, pPoint);

	int iFlag1 = 0;
	int iFlag2 = 0;
	CVessel* newVessel = NULL;
	// if we found such points, break the vessel into two
	if (centerBreakNode)
	{
		newVessel = new CVessel();
		newVessel->m_Center = m_Center.BreakIntoTwo(centerBreakNode);

		// make sure the linked lists were broken successfully
		if (newVessel->m_Center.head == NULL)
			delete newVessel;
		else
		{
			newVessel->m_iLength = newVessel->m_Center.length;
			m_iLength = m_Center.length;

			// The following three lines were added for CANCER Images
			newVessel->m_fHWidth = m_fHWidth;
			newVessel->m_fVWidth = m_fVWidth;
			newVessel->m_iNumOfPoints = m_iNumOfPoints;

			// add the vessel to the collection of vessels
			gTheVessels.AddVessel(newVessel);					

			// determine which of my intersection points belong to me and
			// which belong to the new vessel
			CIntersectionPoint* pIntPoint = NULL;
			TopEndMiddleNeither location = Neither;
			int IntPointId = 0;
			for (register int i = 0; i < m_iNumOfIntersectionPoints; i++)
			{
				IntPointId = m_aiMyIntersectionPoints[i];
				pIntPoint = gIntersectionPoints.GetPoint(IntPointId);

				if (pIntPoint)
				{
					// if this interseciton point lies on this vessel, add it
					// else remove it from its list of intersection points
					location = GetPositionOfPoint(&pIntPoint->m_Point);
					if (location != Neither)
						pIntPoint->AddVessel(m_iID, location);
					else
					{
						pIntPoint->RemoveVessel(m_iID);
						if (m_aiMyIntersectionPoints[i])
						{
							m_aiMyIntersectionPoints[i] = LargeNumber;
							iFlag1 = 1;
						}
					}

					// if this intersection point lies on the new vessel
					// add it to its list of intersection point, 
					// ELSE, remove it
					location = newVessel->GetPositionOfPoint(&pIntPoint->m_Point);
					if (location != Neither)
					{
						pIntPoint->AddVessel(newVessel->m_iID, location);
						newVessel->AddIntersectionPoint(IntPointId);
					}
					else
					{
						pIntPoint->RemoveVessel(newVessel->m_iID);
						if (newVessel->m_aiMyIntersectionPoints[i])
						{
							newVessel->m_aiMyIntersectionPoints[i] = LargeNumber;
							iFlag2 = 1;
						}
					}
				}
			}
		}
	}
	else
	{
		cout << "Error! Can't find closest point to intersection Point ";
		pPoint->Print();
	}

	// now we need to shrink m_aiMyIntersectionPoints array to exclude those with
	// LargeNumber entries (i.e. have been deleted in the above loop). 

	if (iFlag1)
		ShrinkMyIntersectionPoints();

	if (iFlag2)
		newVessel->ShrinkMyIntersectionPoints();
}
/* This is the original function.
 Above we have one that works for the centerline only. 
// break the vessel at the intersection point with the index "index"
void CVessel::BreakOnIntersectionPoint(CPoint *pPoint)
{
	if(EndPoint(pPoint))
		return;

	// find the break points for each of the linked lists in this vessel
	CLNode<CPoint> *centerBreakNode  = FindClosestPoint(&m_Center, pPoint);
	CLNode<CPoint> *HleftBreakNode   = FindClosestPoint(&m_HLeft, pPoint);
	CLNode<CPoint> *HrightBreakNode  = FindClosestPoint(&m_HRight, pPoint);
	CLNode<CPoint> *VleftBreakNode   = FindClosestPoint(&m_VLeft, pPoint);
	CLNode<CPoint> *VrightBreakNode  = FindClosestPoint(&m_VRight, pPoint);

	int iFlag1 = 0;
	int iFlag2 = 0;
	CVessel *newVessel = NULL;
	// if we found such points, break the vessel into two
	if(centerBreakNode && HleftBreakNode && HrightBreakNode &&
		VleftBreakNode && VrightBreakNode )
	{
	
		newVessel = new CVessel();

		newVessel->m_HLeft   = m_HLeft.BreakIntoTwo(HleftBreakNode);
		newVessel->m_HRight  = m_HRight.BreakIntoTwo(HrightBreakNode);
		newVessel->m_VLeft   = m_VLeft.BreakIntoTwo(VleftBreakNode);
		newVessel->m_VRight  = m_VRight.BreakIntoTwo(VrightBreakNode);
		newVessel->m_Center  = m_Center.BreakIntoTwo(centerBreakNode);
		
		// make sure the linked lists were broken successfully
		if(! newVessel->m_HLeft.head || ! newVessel->m_HRight.head ||
			! newVessel->m_VLeft.head || ! newVessel->m_VRight.head ||
			! newVessel->m_Center.head)
			delete newVessel;
		else
		{
			newVessel->m_iLength = newVessel->m_Center.length;
			m_iLength = m_Center.length;

			// add the vessel to the collection of vessels
			gTheVessels.AddVessel(newVessel);					

			// determine which of my intersection points belong to me and
			// which belong to the new vessel
			CIntersectionPoint *pIntPoint = NULL;
			TopEndMiddleNeither location = Neither;
			int IntPointId = 0;
			for(register int i = 0; i < m_iNumOfIntersectionPoints; i++)
			{
				IntPointId = m_aiMyIntersectionPoints[i];
				pIntPoint = gIntersectionPoints.GetPoint(IntPointId);

				if(pIntPoint)
				{
					// if this interseciton point lies on this vessel, add it
					// else remove it from its list of intersection points
					location = GetPositionOfPoint(&pIntPoint->m_Point);
					if(location != Neither)
						pIntPoint->AddVessel(m_iID, location);
					else
					{
						pIntPoint->RemoveVessel(m_iID);
						if(m_aiMyIntersectionPoints[i])
						{
							m_aiMyIntersectionPoints[i] = LargeNumber;
							iFlag1 = 1;
						}
					}
					
					// if this intersection point lies on the new vessel
					// add it to its list of intersection point, 
					// ELSE, remove it
					location = newVessel->GetPositionOfPoint(&pIntPoint->m_Point);
					if(location != Neither)
					{	
						pIntPoint->AddVessel(newVessel->m_iID, location);
						newVessel->AddIntersectionPoint(IntPointId);
					}
					else
					{
						pIntPoint->RemoveVessel(newVessel->m_iID);
						if(newVessel->m_aiMyIntersectionPoints[i])
						{
							newVessel->m_aiMyIntersectionPoints[i] = LargeNumber;
							iFlag2 = 1;
						}
					}
				}
			}
		}
	}
	else
	{
		cout << "Error! Can't find closest point to intersection Point ";
		pPoint->Print();
	}

	// now we need to shrink m_aiMyIntersectionPoints array to exclude those with
	// LargeNumber entries (i.e. have been deleted in the above loop). 
	
	int iIndex = 0;
	int *pTemp = NULL;
	if(iFlag1)
		ShrinkMyIntersectionPoints();

	if(iFlag2)
		newVessel->ShrinkMyIntersectionPoints();
}
*/
////////////////////////////////////
// Method: Shrinkm_aiMyIntersectionPoints
//
// Purpose: When adding/replacing intersection point IDs, the possibility exist
// that the intersection points array contains duplicates. In such situation, 
// remove such duplicates and shrink the intersection points array accordingly
void CVessel::ShrinkMyIntersectionPoints()
{
	register int i, j;
	for (i = 0; i < m_iNumOfIntersectionPoints - 1; i++)
	{
		for (j = i + 1; j < m_iNumOfIntersectionPoints; j++)
		{
			if (m_aiMyIntersectionPoints[i] == m_aiMyIntersectionPoints[j])
				m_aiMyIntersectionPoints[j] = LargeNumber;
		}
	}

	int iIndex = 0;
	int* pTemp = new int[m_iSizeOfIntersectionPointsArray];
	for (i = 0; i < m_iNumOfIntersectionPoints; i++)
	{
		if (m_aiMyIntersectionPoints[i] != LargeNumber)
		{
			pTemp[iIndex] = m_aiMyIntersectionPoints[i];
			iIndex++;
		}
	}
	delete [] m_aiMyIntersectionPoints;
	m_iNumOfIntersectionPoints = iIndex;
	m_aiMyIntersectionPoints = pTemp;
}


///////////////////////////////////
// Method: ExtendVesselCenter
//
// 
int CVessel::ExtendVesselCenter(CPoint* pPoint, TopEndMiddleNeither dirFlag)
{
	// extend the vessel's center but since the flag "gReturnFlag" is 
	// reset to 0,  this forces the extension regardless of the
	// existence of other vessels. 
	giReturnFlag = 0;
	ExtendVessel(m_Center, pPoint, dirFlag);
	giReturnFlag = 1;

	return 1;
}


//////////////////////
// Method: ReplaceIntersectionPoint
// 
// if I am using the oldID, replace it with the newID
void CVessel::ReplaceIntersectionPoint(int oldID, int newID)
{
	int iFlag = 0;
	for (register int i = 0; i < m_iNumOfIntersectionPoints; i++)
	{
		if (m_aiMyIntersectionPoints[i] == oldID)
		{
			m_aiMyIntersectionPoints[i] = newID;
			iFlag = 1;
			break;
		}
	}
	// the old point ID was not found, add the new ID
	if (iFlag == 0)
		AddIntersectionPoint(newID);

	// in case this operation produced duplicates, remove them
	ShrinkMyIntersectionPoints();
}

///////////////////////////////////////////////////////////////////////////////
// UpdateLength
int CVessel::UpdateLength()
{
	m_iLength = m_Center.length;
	return m_iLength;
}

///////////////////////////////////////////////////////////////////////////////
// Compute the distance from the point "FromPoint" to the list "List"
// and point to the closest Node in the list by "ClosestNode"
//
int CVessel::ComputeDistanceToList(CPoint& FromPoint, CLNode<CPoint>* List,
	CLNode<CPoint>** pClosestNode)
{
	int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;
	int iMinDistance = iCols* iRows* iSlices;
	int iDistance = 0;
	int dx, dy, dz;

	CLNode<CPoint>* temp = List;
	while (temp)
	{
		if (*(temp->data) == FromPoint)
		{
			*pClosestNode = temp;
			return 0;
		}

		dx = temp->data->m_iX - FromPoint.m_iX;
		dy = temp->data->m_iY - FromPoint.m_iY;
		dz = temp->data->m_iZ - FromPoint.m_iZ;

		iDistance = dx * dx + dy * dy + dz * dz;
		if (iDistance < iMinDistance)
		{
			iMinDistance = iDistance;
			*pClosestNode = temp;
		}

		temp = temp->after;
	}

	return static_cast<int>( (std::sqrt(static_cast<double>(iMinDistance)) ) );
}


///////////////////////////////
// METHOD: FindClosestCenterlinePoint
//
int CVessel::FindClosestCenterlinePoint(CPoint& aPoint,
	CLNode<CPoint>** pResult)
{
	return ComputeDistanceToList(aPoint, m_Center.head, pResult);
}

///////////////////////////////
// METHOD: FindClosestHLeftPoint
//
int CVessel::FindClosestHLeftPoint(CPoint& aPoint, CLNode<CPoint>** pResult)
{
	return ComputeDistanceToList(aPoint, m_HLeft.head, pResult);
}

///////////////////////////////
// METHOD: FindClosestHRightPoint
//
int CVessel::FindClosestHRightPoint(CPoint& aPoint, CLNode<CPoint>** pResult)
{
	return ComputeDistanceToList(aPoint, m_HRight.head, pResult);
}

///////////////////////////////
// METHOD: FindClosestVLeftPoint
//
int CVessel::FindClosestVLeftPoint(CPoint& aPoint, CLNode<CPoint>** pResult)
{
	return ComputeDistanceToList(aPoint, m_VLeft.head, pResult);
}

///////////////////////////////
// METHOD: FindClosestVRightPoint
//
int CVessel::FindClosestVRightPoint(CPoint& aPoint, CLNode<CPoint>** pResult)
{
	return ComputeDistanceToList(aPoint, m_VRight.head, pResult);
}


//////////////////////////////////////////////////////
// Method: RemoveIntersectionPoint
//
void CVessel::RemoveIntersectionPoint(int iID)
{
	int iIndex = - 1;
	int i;
	for (i = 0; i < m_iNumOfIntersectionPoints; i++)
	{
		if (m_aiMyIntersectionPoints[i] == iID)
		{
			iIndex = i;
			break;
		}
	}

	if (iIndex != -1)
	{
		for (i = iIndex; i < m_iNumOfIntersectionPoints - 1; i++)
		{
			m_aiMyIntersectionPoints[i] = m_aiMyIntersectionPoints[i + 1];
		}
		m_iNumOfIntersectionPoints--;
	}
}

//////////////////////////////////////////////////////
// METHOD: PrintNeurolucidaFormat
// print the vessel's centerline in Neurolucida format to the given file
void CVessel::PrintNeurolucidaFormat(ostream& outFile,
	TopEndMiddleNeither DirFlag)
{
	CLNode<CPoint>* temp = m_Center.head;

	if (DirFlag == OnEnd)
		temp = m_Center.tail;

	while (temp)
	{
		outFile << "(" << temp->data->m_iX << " " << temp->data->m_iY << " "
			<< temp->data->m_iZ << ")\n";

		if (DirFlag == OnTop)
			temp = temp->after;
		else
			temp = temp->before;
	}
}

extern float median__(float a, float b, float c);

void SmoothWidths(vector<float> &widths) 
{	

	size_t i = 0;
	size_t j = 0;
	vector<float> result;
	for(i = 0; i < widths.size(); i++) {
		if( i == 0) {
			if(widths.size() <= 1)
				result.push_back(widths[i]);
			else
				result.push_back((widths[i] + widths[i+1])/2);
		}
		else {
			if(i != widths.size() - 1) {
				result.push_back(median__(widths[i], widths[i-1], widths[i+1]));
			}
			else
				result.push_back((widths[i] + widths[i-1])/2);
		}
	}

	// get the median width and saturate the width at 0.5 and 2.0 times the median

	list<float> values;
	for(j = 0; j < result.size(); j++)
		values.push_back(result[j]);
	float median = Median(values);
	//	cout << median << endl;

	for(j = 0; j < result.size(); j++) {
		if(result[j] > 2.0f * median)
			result[j] = 2.0f * median;
		if(result[j] < 0.5f * median) 
			result[j] = max(0.5f * median, 1.0f);
	}

	if(result.size() != widths.size()) {
		cerr<< "SmoothWidths:ERROR" << endl;
		return;
	}
	else {
		for(j = 0; j < result.size(); j++)
			widths[j] = result[j];
	}

}

// amri 6-16-05 BUGFIX: major revision
// Revise vessel widths only on a complete vessel
void CVessel::ReviseWidths()
{
	int iMaxShiftDistance = gConfig.GetMaximumShiftDistance();
	CLNode<CPoint>* temp = m_Center.head;		
	float h_width = temp->data->m_fHWidth;
	float v_width = temp->data->m_fVWidth;
	float hwidth_sum = 0.0;
	float vwidth_sum = 0.0;
		
	vector<float> hwidth;
	vector<float> vwidth;

	// collect the widths for smoothing
	temp = m_Center.head;
	while (temp)
	{
		if( Round(temp->data->m_fHWidth) != 0 || Round(temp->data->m_fVWidth != 0) ) {
			hwidth.push_back(temp->data->m_fHWidth);
			vwidth.push_back(temp->data->m_fVWidth);
		}
		temp = temp->after;
	}

	// smooth the widths
	SmoothWidths(hwidth);
	SmoothWidths(vwidth);
	
	// copy back to smoothed widths
	temp = m_Center.head;
	int index = 0;
	while (temp)
	{
		if( Round(temp->data->m_fHWidth) != 0 || Round(temp->data->m_fVWidth != 0) ) {
			temp->data->m_fHWidth = hwidth[index];
			temp->data->m_fVWidth = vwidth[index];
			index++;
		}
		temp = temp->after;
	}
	
	// for points without width info, extrapolate from nearest point with width info
	temp = m_Center.head;
	while (temp)
	{
		if (temp->data->m_fHWidth < iMaxShiftDistance &&
			temp->data->m_fVWidth < iMaxShiftDistance)
		{
			if (temp->data->m_fHWidth)
				h_width = temp->data->m_fHWidth;
			else
				temp->data->m_fHWidth = h_width;
			if (temp->data->m_fVWidth)
				v_width = temp->data->m_fVWidth;
			else
				temp->data->m_fVWidth = v_width;

			// if the width is still zero
			if (temp->data->m_fHWidth < 1.0)
				temp->data->m_fHWidth = 1.0;
			if (temp->data->m_fVWidth < 1.0)
				temp->data->m_fVWidth = 1.0;
		}
		else
		{
			if (temp->data->m_fHWidth >= iMaxShiftDistance)
				temp->data->m_fHWidth = static_cast<float>(iMaxShiftDistance);
			if (temp->data->m_fVWidth >= iMaxShiftDistance)
				temp->data->m_fVWidth = static_cast<float>(iMaxShiftDistance);
		}
		temp = temp->after;
	}

	vector<float> h_widths;
	vector<float> v_widths;

	// do a gaussian smooth on the widths
	// width_i = 1/4 width_i-1 + 1/2 width_i + 1/4 width_i+1
	temp = m_Center.head;
	while (temp)
	{
		if (temp->before && temp->after)
		{
			h_widths.push_back(0.25f * temp->before->data->m_fHWidth +
				0.5f * temp->data->m_fHWidth +
				0.25f * temp->after->data->m_fHWidth);
			v_widths.push_back(0.25f * temp->before->data->m_fVWidth +
				0.5f * temp->data->m_fVWidth +
				0.25f * temp->after->data->m_fVWidth);
		}
		else
		{
			h_widths.push_back(temp->data->m_fHWidth);
			v_widths.push_back(temp->data->m_fVWidth);
		}
		temp = temp->after;
	}

	temp = m_Center.head;
	for (size_t i = 0; i < h_widths.size(); i++)
	{
		if (temp->before && temp->after)
		{
			temp->data->m_fHWidth = h_widths[i];
			temp->data->m_fVWidth = v_widths[i];
		}
		hwidth_sum += h_widths[i];
		vwidth_sum += v_widths[i];
		temp = temp->after;
	}

	this->m_fHWidth = hwidth_sum / (float)this->m_iLength;
	this->m_fVWidth = vwidth_sum / (float)this->m_iLength;
}

void CVessels::ReviseWidths()
{
	for (int i = 0; i < m_iNumOfElements; i++)
	{
		if (m_apData[i])
		{
			if (m_apData[i]->m_iLength > 0)
				m_apData[i]->ReviseWidths();
		}
	}
}

// Vessel Marking Using Deformed Spheres
// amri 7-19-02
void CVessel::MarkVessel()
{
	CLNode<CPoint>* temp = m_Center.head;	
	
	// prebuild the structuring elements
	map< pair<int,int>, StrEle *> ellipsoids;	
	while(temp) {
		pair<int,int> dimension;
		dimension.first = Round(temp->data->m_fHWidth);
		dimension.second = Round(temp->data->m_fVWidth);
		if(!ellipsoids.count(dimension))
			ellipsoids[dimension] = new StrEle(The3DImage, "sphere", dimension.first, dimension.second);
		temp = temp->after;
	}

	temp = m_Center.head;

	int i;
	for(i = 0; i < Round(temp->data->m_fHWidth)/2;i++)
		temp = temp->after;

	CLNode<CPoint>* probe;

	while (temp)
	{
		pair<int,int> size;
		size.first = Round(temp->data->m_fHWidth);
		size.second = Round(temp->data->m_fVWidth);
		StrEle * sphere = ellipsoids[size];
		for (int j = 0; j < sphere->GetSize(); j++)
		{
			// if the center point plus the offsets is within the image boundary, then mark it
			if(The3DImage->WithinImagePadding(temp->data->m_iZ + sphere->GetOffsets()[j][0],
				temp->data->m_iY + sphere->GetOffsets()[j][1],
				temp->data->m_iX + sphere->GetOffsets()[j][2], giMARGIN))
			{
				TracedImage[temp->data->m_iZ + sphere->GetOffsets()[j][0]] [temp->data->m_iY + sphere->GetOffsets()[j][1]]
				[temp->data->m_iX + sphere->GetOffsets()[j][2]] = true;
			}
		}
		probe = temp;
		i = 0;
		while(i < Round(temp->data->m_fHWidth)/2)
		{
			if(!(probe = probe->after)) break;
			i++;
		}
		if(i == Round(temp->data->m_fHWidth)/2)
			temp = temp->after;
		else
			break;
	}
}

void CVessel::Mark_In_3D_Image(C3DImage& image, unsigned char color)
{
	CLNode<CPoint>* temp = m_Center.head;
	
	// prebuild the structuring elements
	map< pair<int,int>, StrEle *> ellipsoids;	
	while(temp) {
		pair<int,int> dimension;
		dimension.first = Round(temp->data->m_fHWidth);
		dimension.second = Round(temp->data->m_fVWidth);
		if(!ellipsoids.count(dimension))
			ellipsoids[dimension] = new StrEle(&image, "sphere", dimension.first, dimension.second);
		temp = temp->after;
	}

	temp = m_Center.head;
	while (temp)
	{
		pair<int,int> size;
		size.first = Round(temp->data->m_fHWidth);
		size.second = Round(temp->data->m_fVWidth);
		StrEle * sphere = ellipsoids[size];
		sphere->SetCenter(temp->data->m_iZ, temp->data->m_iY, temp->data->m_iX);
		sphere->Fill(color);
		temp = temp->after;
	}
}

void CVessels::Mark_In_3D_Image(C3DImage& image, unsigned char  color)
{
	for (int i = 0; i < m_iNumOfElements; i++)
	{
		if (m_apData[i])
		{
			m_apData[i]->Mark_In_3D_Image(image, color);
		}
	}
}

// amri tried to add this
//void CVessel::PrintNeurolucidaFormat(ostream &outFile)
//{
//	CLNode<CPoint> *temp = m_Center.head;
//
//	outFile << "("
//
//	while(temp) {
//		outFile << "(" 
//			<< temp->data->m_iX << " "
//			<< temp->data->m_iY << " " 
//			<< temp->data->m_iZ 
//			<< ")\n";
//			temp = temp->after;
//	}
//}
//
////////////////////////////////////////////////////////
//// METHOD: PrintNeurolucidaFormat
//// print the vessel's centerline in Neurolucida format to the given file
//void CVessels::PrintNeurolucidaFormat(ostream &outFile)
//{
//	for(register int i = 0; i < m_iNumOfElements; i++) {
//		m_apData[i]->PrintNeurolucidaFormat(outFile);
//	}
//}


///////////////////////////////////////////////////////////////////////////////
//									CLASS CVessels
///////////////////////////////////////////////////////////////////////////////

void CVessels::AddVessel(CVessel* aVessel)
{
	if (aVessel)
	{
		if (m_iNumOfElements < m_iSize)
		{
			m_apData[m_iNumOfElements] = aVessel;
			m_iNumOfElements++;
		}
		else
		{
			CVessel** temp = m_apData;
			m_apData = new CVessel * [m_iSize + BlockSize];
			memcpy(m_apData, temp, sizeof(CVessel *) * m_iNumOfElements);
			m_apData[m_iNumOfElements] = aVessel;
			m_iNumOfElements++;
			m_iSize += BlockSize;
			delete [] temp;
		}
		aVessel->m_iID = m_iNextID;
		m_iNextID++;
		aVessel->m_iLength = aVessel->m_Center.length;
	}
}

void CVessels::DeleteShortNetwork(int min_length)
{
	bool* to_be_deleted = new bool[m_iNumOfElements];
	memset(to_be_deleted, 0, sizeof(bool) * m_iNumOfElements);

	int i;
	for (i = 0; i < m_iNumOfElements; i++)
	{
		if (m_apData[i] != 0 && m_apData[i]->m_iDrawFlag == 0)
		{
			if (m_apData[i]->m_iNumOfIntersectionPoints > 0 ||
				m_apData[i]->m_iLength > min_length)
			{
			}
			else
			{
				to_be_deleted[i] = true;
			}
		}
		else
			to_be_deleted[i] = true;
	}

	CVessel** temp = new CVessel*[m_iSize];

	int vessels_kept = 0;
	for (i = 0; i < m_iNumOfElements; i++)
	{
		if (!to_be_deleted[i])
		{
			temp[vessels_kept] = m_apData[i];
			vessels_kept++;
		}
		//	else
			//	Draw_XYCenterline(m_apData[i]);
	}
	delete [] m_apData;
	m_apData = temp;
	m_iNumOfElements = vessels_kept;
	//	Draw_XYCenterline(0);
}

void CVessels::DeleteNarrowVessels(double min_width)
{
	bool* to_be_deleted = new bool[m_iNumOfElements];
	memset(to_be_deleted, 0, sizeof(bool) * m_iNumOfElements);
	int i;

	for (i = 0; i < m_iNumOfElements; i++)
	{
		if (m_apData[i] != 0)
		{
			if (m_apData[i]->m_fHWidth < min_width)
				to_be_deleted[i] = true;
		}
		else
			to_be_deleted[i] = true;
	}

	CVessel** temp = new CVessel*[m_iSize];

	int vessels_kept = 0;
	for (i = 0; i < m_iNumOfElements; i++)
	{
		if (!to_be_deleted[i])
		{
			temp[vessels_kept] = m_apData[i];
			vessels_kept++;
		}
	}
	delete [] m_apData;
	m_apData = temp;
	m_iNumOfElements = vessels_kept;
}

////////////////////////////////////////////////////////////////////
// Method: ExtendVessels
//
// Purpose: extend vessel end points to connect with other vessels
void CVessels::ExtendVessels()
{
	int iFromLength = gConfig.GetMinimumTemplateLength();
	int iToLength = gConfig.GetMaximumTemplateLength();

	giMaxExtensionLength = iFromLength + 2;

	giSomaDistThreshold = iToLength + 2;

	int i;
	// connect vessels with somas, and determine intervessel connectivity
	for (i = 0; i < m_iNumOfElements; i++)
	{
		if (m_apData[i])
		{
			m_apData[i]->ConnectWithSomasAndOtherVessels3();
		}
	}
	
	// do the actual connections.
	for (i = 0; i < m_iNumOfElements; i++)
	{
		if (m_apData[i])
			m_apData[i]->ConnectWithOtherVessels3();
	}

	return;

	int iDistToSoma = _INFINITY_;
	int iTopDistToVessel = _INFINITY_;
	int iEndDistToVessel = _INFINITY_;
	int iSomaID = 0;
	int iFlag = 0;
	register int  j;
	TopEndMiddleNeither pos = Neither;
	CVessel* pVessel = NULL;
	CPoint SomaPoint;
	CIntersectionPoint* pIntPoint = NULL;

	// Finally, for each intersection point that does not have a path to a soma
	// and inolving two vessels connected at their tips, Extend one of them to the soma
	for (i = 0; i < gIntersectionPoints.m_iNumOfElements; i++)
	{
		pIntPoint = gIntersectionPoints.m_apData[i];
		if (pIntPoint->m_Type == VESSEL)
		{
			iFlag = 0;
			// The point must be at the tips of each of the intersecting vessels
			for (j = 0; j < pIntPoint->m_iIndex; j++)
			{
//				BUG? vijay & amri - 11-25.02	
//				if (pIntPoint->m_aPosition[i] == Middle)
//				{
//					iFlag = 1;
//					break;
//				}
				if (pIntPoint->m_aPosition[j] == Middle)
				{
					iFlag = 1;
					break;
				}
			}

			if (iFlag)
				continue; 

			// if there is no path to a soma
			if (! gIntersectionPoints.IsThereAPathToASoma(pIntPoint->m_iID))
			{
				// Extend one of the vessels to the soma
				pVessel = GetVessel(pIntPoint->m_aVesselIDs[0]);
				pos = pIntPoint->m_aPosition[0];

				iDistToSoma = pVessel->FindSomaDistance3(pos,
									   	SomaPoint,
									   	iSomaID);
				if (iDistToSoma < giSomaDistThreshold)
					pVessel->ConnectWithSoma(iSomaID, SomaPoint, pos);
				else
				{
					//	continue;
					// see if you can connect it to a vessel

					iTopDistToVessel = FindClosestVesselExcluding(pVessel->m_iID,
									   	pVessel->m_Center.head->data,
									   	& pVessel->m_pTopClosestVessel,
									   	& pVessel->m_pTopClosestNode,
									   	pIntPoint);
					iEndDistToVessel = FindClosestVesselExcluding(pVessel->m_iID,
									   	pVessel->m_Center.tail->data,
									   	& pVessel->m_pEndClosestVessel,
									   	& pVessel->m_pEndClosestNode,
									   	pIntPoint);

					if (iTopDistToVessel <= iEndDistToVessel &&
						iTopDistToVessel <= iFromLength + 2)
					{
						pVessel->m_iConnectOnTopFlag = 1;
						pVessel->ConnectWithOtherVessels3(pIntPoint);
					}
					else if (iEndDistToVessel <= iTopDistToVessel &&
						iEndDistToVessel <= iFromLength + 2)
					{
						pVessel->m_iConnectOnEndFlag = 1;
						pVessel->ConnectWithOtherVessels3(pIntPoint);
					}
				}
			}
		}
	}
}

void CVessels::DrawVessels(CImage& anImage, unsigned char color)
{
	for (register int i = 0; i < m_iNumOfElements; i++)
	{
		if (m_apData[i] != 0 && m_apData[i]->m_iDrawFlag == 0)
		{
			//			if (m_apData[i]->m_iNumOfIntersectionPoints > 0 ||
			//				m_apData[i]->m_iLength > giParam_Trace_MinSegmentLength)
			//			{

			///////////////////////////////////////////////////
			//By Yousef, 11/09/2006
			//Remove the vessels shorter than 20 
			int l = m_apData[i]->GetLength();
			if(l<12)
				continue;
			///////////////////////////////////////////////////

			m_apData[i]->DrawCenterline(anImage, color);
			//m_apData[i]->DrawBoundaries(anImage, color);
			//if (m_apData[i]->m_iID < 100)
				//m_apData[i]->WriteIDXY(anImage, IDColor);
			//	}
		}
	}
}

void CVessels::UpdateIntersectionPoints()
{
	//	ofstream out("updateIntersection.txt");
	int i, j, k;
	for (i = 0; i < gIntersectionPoints.m_iNumOfElements; i++)
	{
		//	out << i << " Intersection " << gIntersectionPoints.m_apData[i]->m_iID << " "; 
		//	gIntersectionPoints.m_apData[i]->m_Point.Print(out);
		for (j = 0; j < gIntersectionPoints.m_apData[i]->m_iIndex; j++)
		{
			//			out << "connected vessel ID = " << gIntersectionPoints.m_apData[i]->m_aVesselIDs[j] << endl;
			//			cout << "babi " << endl;
			for (k = 0; k < this->m_iNumOfElements; k++)
			{
				if (this->m_apData[k]->m_iID ==
					gIntersectionPoints.m_apData[i]->m_aVesselIDs[j])
				{
					if (this->m_apData[k])
					{
						this->m_apData[k]->AddIntersectionPoint(gIntersectionPoints.m_apData[i]->m_iID);
						//	this->m_apData[k]->Print(out);
					}
				}
			}
		}
		//	out << endl << endl;
	}
}


void CVessel::Draw3DCenterline(C3DImage& anImage, unsigned char color)
{
	CLNode<CPoint>* tempNode = NULL;
	tempNode = m_Center.head;
	while (tempNode)
	{
		anImage.MarkSinglePoint(tempNode->data, color);
		tempNode = tempNode->after;
	}
}

void CVessels::Draw3DCenterline(C3DImage& anImage, unsigned char color)
{
	for (register int i = 0; i < m_iNumOfElements; i++)
	{
		if (m_apData[i] != 0 && m_apData[i]->m_iDrawFlag == 0)
		{
			m_apData[i]->Draw3DCenterline(anImage, color);
		}
	}
}



void CVessels::DrawVesselsXZ(CImage& anImage, unsigned char color)
{
	for (register int i = 0; i < m_iNumOfElements; i++)
	{
		if (m_apData[i] != 0 && m_apData[i]->m_iDrawFlag == 0)
		{
			m_apData[i]->DrawCenterlineXZ(anImage, color);
			//m_apData[i]->DrawBoundariesXZ(anImage, color);
			if (m_apData[i]->m_iID < 100)
				m_apData[i]->WriteIDXZ(anImage, IDColor);
		}
	}
}

void CVessels::DrawVesselsYZ(CImage& anImage, unsigned char color)
{
	for (register int i = 0; i < m_iNumOfElements; i++)
	{
		if (m_apData[i] != 0 && m_apData[i]->m_iDrawFlag == 0)
		{
			m_apData[i]->DrawCenterlineYZ(anImage, color);
			//m_apData[i]->DrawBoundariesYZ(anImage, color);
			if (m_apData[i]->m_iID < 100)
				m_apData[i]->WriteIDYZ(anImage, IDColor);
		}
	}
}
void CVessels::DrawVessels(C3DImage& anImage, unsigned char color)
{
	for (register int i = 0; i < m_iNumOfElements; i++)
	{
		if (m_apData[i] != 0)
		{
			m_apData[i]->DrawCenterline(anImage, color);
			//data[i]->WriteID(anImage, 254);
		}
	}
}

/*
void CVessels::RemoveAllDuplicateIntersectionPoints()
{
	for(register int i = 0; i < m_iNumOfIntersectionPoints; i++)
	{
		int v1 = m_aiMyIntersectionPoints[i]->VesselID1;
		int v2 = muIntersectionPoints[i]->VesselID2;

		for(register int j = i+1; j < m_iNumOfIntersectionPoints; j++)
		{
			int v3 = 
		}
	}
}
*/
//////////////////////////////////////////////////////////////////////////////
// Method: CombineVesselFragments
//
// Purpose: to combine two vessel fragments that are seperated by another intersecting
//  		vessel.
void CVessels::MergeVessels()
{
	register int i;
	int iPosFlag = 0;
	CIntersectionPoint* pIntPoint = NULL;
	TopEndMiddleNeither pos1, pos2;
	CVessel* pVessel1,* pVessel2;

	// for each intersection point, if it ivolves two vessels only connected
	// at their tips, merge them	
	cout << "Number of points: %d" << gIntersectionPoints.m_iNumOfElements;
	for (i = 0; i < gIntersectionPoints.m_iNumOfElements; i++)
	{		
		cout << "%d" << i;
		pIntPoint = gIntersectionPoints.m_apData[i];

		pos1 = pIntPoint->m_aPosition[0];
		pos2 = pIntPoint->m_aPosition[1];

		iPosFlag = 0;

		if (pIntPoint->m_iIndex == 2 &&
			pIntPoint->m_Type == VESSEL &&
			pos1 != Middle &&
			pos2 != Middle)
		{
			// Merge the two vessels.
			pVessel1 = GetVessel(pIntPoint->m_aVesselIDs[0]);
			pVessel2 = GetVessel(pIntPoint->m_aVesselIDs[1]);
			pVessel1->ExtendByAnotherVessel(pVessel2, pos1, pos2);
			pVessel2->m_iMergedFlag = 1;

			// Inform all other intersection points.
			if (pos1 == pos2)
				iPosFlag = 1;
			gIntersectionPoints.UpdateVesselID(pVessel2->m_iID,
									pVessel1->m_iID,
									iPosFlag);

			// remove the intersection point from the vessel's points
			pVessel1->RemoveIntersectionPoint(pIntPoint->m_iID);

			// Mark this intersection point for removal.
			pIntPoint->m_ProcessedFlag = MERGED;
		}
	}

	// for each vessel that is marked as merged, remove it from the collection
	CVessel** temp = new CVessel*[m_iSize];
	int iIndex = 0;
	for (i = 0; i < m_iNumOfElements; i++)
	{
		if (m_apData[i]->m_iMergedFlag == 0)
		{
			temp[iIndex] = m_apData[i];
			iIndex++;
		}
	}
	delete [] m_apData;
	m_apData = temp;
	m_iNumOfElements = iIndex;

	// do the same for the intersection points
	gIntersectionPoints.RemovePointsBetweenMergedVessels();
}

////////////////////////////////////////////////////////////////////////////////
// METHOD: MergeTwoVessels
//
// Purpose: The given two vessels intersect with a third vessel at the given 
//  		intersection points (respectively). 
void CVessels::MergeTwoVessels(int , int , int , int , int )
{
	/*
	cout << "\n The segments: " << v1+1 << " and " << v2+1 
		  << " have been merged into one." << endl;
	if(!(data[v1] && data[v2] && data[v3]))
		return;
	if(data[v1]->center->length < data[v2]->center->length)
	{
		int temp = v1;
		v1 = v2;
		v2 = temp;
		temp = p1;
		p1 = p2;
		p2 = temp;
	}
	data[v1]->ExtendByAnotherVessel(data[v2], m_aiMyIntersectionPoints[p1]->thePoint,
			m_aiMyIntersectionPoints[p2]->thePoint);
	data[v3]->RemoveIntersectionPoint(m_aiMyIntersectionPoints[p2]->PointID);
	m_aiMyIntersectionPoints[p1]->Type = Intersection;
	m_aiMyIntersectionPoints[p2] = 0;
	data[v2] = 0;
	*/
}
/*
int CVessels::AddIntersectionPoint(IntersectionPoint *aIntPoint)
{
	int result = 0;
	if(aIntPoint)
	{
		aIntPoint->PointID = m_iNumOfIntersectionPoints;
		result = m_iNumOfIntersectionPoints;


		if(m_iNumOfIntersectionPoints < iSizeOfIntersectionPointsArray)
		{
			m_aiMyIntersectionPoints[m_iNumOfIntersectionPoints] = aIntPoint;
			m_iNumOfIntersectionPoints++;
		}
		else
		{
			IntersectionPoint **temp = m_aiMyIntersectionPoints;
			m_aiMyIntersectionPoints = new IntersectionPoint* [iSizeOfIntersectionPointsArray + BlockSize];
			memcpy(m_aiMyIntersectionPoints, temp, sizeof(IntersectionPoint *) * m_iNumOfIntersectionPoints);
			m_aiMyIntersectionPoints[m_iNumOfIntersectionPoints] = aIntPoint;
			m_iNumOfIntersectionPoints++;
			iSizeOfIntersectionPointsArray += BlockSize;
			delete [] temp;
		}
	}

	return result;
}
*/


//////////////////////////////////////////////////////////////////////////////
// METHOD: ResolveBranches
//
// Figure out which vessel is a branch of which other vessel. 
// Logic:
//  
//    For each branching point, assign the vessel that ends at the other
//   as a branch of it. If neither ends at the other, designate the shorter
//   vessel as a branch of the longer intersecting vessel. 
//   If a vessel is designated as a branch of two or more vessels, 
//   the longer vessel wins.
//
/* 
void CVessels::ResolveBranches()
{

	int Branch = -1;
	int Parent = -1;
	// for each branching point
	for(register int i = 0; i < m_iNumOfIntersectionPoints; i++)
	{
		if(m_aiMyIntersectionPoints[i]->Type == Branching)
		{
			int v1 = m_aiMyIntersectionPoints[i]->VesselID1 - 1;
			int v2 = m_aiMyIntersectionPoints[i]->VesselID2 - 1;
			CPoint *aPoint = m_aiMyIntersectionPoints[i]->thePoint;

			if( !(data[v1] && data[v2]) )
				continue;

			// assign the shorter one as a branch to the longer one
			if(data[v1]->EndPoint(aPoint) && !data[v2]->EndPoint(aPoint))
			{
				Branch = v1;
				Parent = v2;
			}
			else if(!data[v1]->EndPoint(aPoint) && data[v2]->EndPoint(aPoint))
			{
				Branch = v2;
				Parent = v1;
			}
			else if(data[v1]->center->length > data[v2]->center->length)
			{
				Branch = v2;
				Parent = v1;
			}
			else
			{
				Branch = v1;
				Parent = v2;
			}
			// if the branch already have a parent
			if(data[Branch]->Parent != -1)
			{
				if(data[Parent]->center->length > 
					data[data[Branch]->Parent]->center->length)
					data[Branch]->Parent = Parent;
			}
			else
			{
				data[Branch]->Parent = Parent;
			}
			data[Parent]->AddBranch(Branch);
		}
	}

	for(i = 0; i < iNumOfElements; i++)
	{
		if(data[i] == 0)
		{
			for(register int j = i; j < iNumOfElements - 1; j++)
				data[j] = data[j+1];

			iNumOfElements--;
			i--;
		}
	}

}

*/

// for each vessel, break it at its intersection points
void CVessels::BreakOnIntersectionPoints()
{
	CVessel* pVessel = NULL;
	CIntersectionPoint* pIntPoint = NULL;
	for (register int i = 0; i < gIntersectionPoints.m_iNumOfElements; i++)
	{
		pIntPoint = gIntersectionPoints.m_apData[i];

		// except those indicating soma intersection
		if (pIntPoint->m_Type == SOMA)
			continue;

		for (register int j = 0; j < pIntPoint->m_iIndex; j++)
		{
			if ((pVessel = GetVessel(pIntPoint->m_aVesselIDs[j])) != 0)
				pVessel->BreakOnIntersectionPoint(&(pIntPoint->m_Point));
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// METHOD: UpdateVesselLengths
//
// Update the length of each vessel to include its branches
//
// Logic: Start with all vessels that have only one branch, then two and so on. At
//  	  each level update the length of the parent vessel
void CVessels::UpdateVesselLengths()
{
	for (register int i = 0; i < m_iNumOfElements; i++)
	{
		m_apData[i]->UpdateLength();
	}
}


int compare3(const void* , const void* )
{
	return 0;
	/*
	CVessel **v1 = ((CVessel **) arg1);
	CVessel **v2 = ((CVessel **) arg2);
	int Sum1 = (*v1)->Center.m_iLength + (*v1)->BranchesSum;
	int Sum2 = (*v2)->Center.m_iLength + (*v2)->BranchesSum;
	if(Sum1 < Sum2)
		return 1;
	else if(Sum1 > Sum2)
		return -1;
	else
		return 0;
		*/
}


void CVessels::Print(ostream& out)
{
	out <<
		"ID    Starts      Ends     Length   BranchesSum    Parent   Brnaches\n"
		<<
		"========================================================================="
		<<
		endl;
	for (register int i = 0; i < m_iNumOfElements; i++)
	{
		if (m_apData[i])
			m_apData[i]->Print(out);
	}
}

void CVessels::Print(char* fName)
{
	register int i;
	register int iLength = 0;
	register float fHWidth = 0;
	register float fVWidth = 0;

	ofstream outFile(fName);

	for (i = 0; i < m_iNumOfElements; i++)
	{
		if (m_apData[i])
		{
			fHWidth += (m_apData[i]->m_fHWidth * m_apData[i]->m_Center.length);
			fVWidth += (m_apData[i]->m_fVWidth * m_apData[i]->m_Center.length);
			iLength += (m_apData[i]->m_Center.length);
		}
	}

	outFile <<
		"A total of " <<
		giNumOfSeedPoints
		<<
		" Seed Points were Found.\n" <<
		"Based on the seed points, \n"
		<<
		"\tThe Average HWidth is: " <<
		gfHWidth /
		giNumOfSeedPoints <<
		"\n"
		<<
		"\tThe Average VWidth is: " <<
		gfVWidth /
		giNumOfSeedPoints <<
		endl;

	outFile << "\nA total of " << m_iNumOfElements
		<< " vessel segments were found.\n"
		<< "The total length of all vessels is: " << iLength << "\n"
		<< "Based on all segments, \n" << "\tThe Average HWidth is: "
		<< fHWidth / iLength << "\n" << "\tThe Average VWidth is: "
		<< fVWidth / iLength << "\n\nThe segments are:\n\n" << endl;



	//outFile << "ID\tStart\tEnd\tLength\tActualStart\tActualEnd\tNoOfOffCenterPoints" << endl;
	outFile <<
		" ID\t   Start\t    End\t       Length\tHWidth\tVWidth" <<
		endl;

	for (i = 0; i < m_iNumOfElements; i++)
	{
		if (m_apData[i])
			m_apData[i]->Print(outFile);
	}
}


int CVessels::FindClosestVessel(int vesselID, CPoint* FromPoint,
	CVessel** pClosestVessel, CLNode<CPoint>** pFoundNode)
{
	int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;
	int iMinDistance = iCols* iRows* iSlices;
	int iDistance = 0;
	register int i;

	CLNode<CPoint>* pTempNode = NULL;
	// otherwise try to connect with other vessels
	// for each vessel other than this, find the closest centerline point
	for (i = 0; i < m_iNumOfElements; i++)
	{
		if (m_apData[i])
		{
			if (m_apData[i]->m_iID == vesselID)
				continue;

			iDistance = m_apData[i]->FindClosestCenterlinePoint(*FromPoint,
									 	& pTempNode);

			if (iDistance < iMinDistance)
			{
				iMinDistance = iDistance;
				*pClosestVessel = m_apData[i];
				*pFoundNode = pTempNode;
			}
		}
	}

	return iMinDistance;
}

//////////////////////////////////////////////
// Same like above but only consider those vessels not belonging to the
// give intersection point
int CVessels::FindClosestVesselExcluding(int , CPoint* FromPoint,
	CVessel** pClosestVessel, CLNode<CPoint>** pFoundNode,
	CIntersectionPoint* pIntPoint)
{
	int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;
	int iMinDistance = iCols* iRows* iSlices;
	int iDistance = 0;
	register int i;

	CLNode<CPoint>* pTempNode = NULL;

	for (i = 0; i < m_iNumOfElements; i++)
	{
		if (m_apData[i])
		{
			if (pIntPoint->GetLocation(m_apData[i]->m_iID) != Neither)
				continue;

			iDistance = m_apData[i]->FindClosestCenterlinePoint(*FromPoint,
									 	& pTempNode);

			if (iDistance < iMinDistance)
			{
				iMinDistance = iDistance;
				*pClosestVessel = m_apData[i];
				*pFoundNode = pTempNode;
			}
		}
	}

	return iMinDistance;
}
// find the node in the list aList that is closest to the
// point aPoint. 
CLNode<CPoint>* FindClosestPoint(CDLList<CPoint>* aList, CPoint* aPoint)
{
	if (! (aList && aPoint))
		return NULL;

	CLNode<CPoint>* result = NULL;
	CLNode<CPoint>* temp = aList->head;
	double distance = 99999;
	double minDistance = 9999;

	while (temp)
	{
		distance = temp->data->FindDistance(aPoint);

		if (distance == 0)
		{
			result = temp;
			break;
		}
		else if (distance < minDistance)
		{
			minDistance = distance;
			result = temp;
		}

		temp = temp->after;
	}

	return result;
}

