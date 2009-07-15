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

//#pragma warning(disable:4786)  // disable STL-related warnings
//#pragma warning(disable:4710)  // compiler didn't inline function that got selected for inline expansion
//#pragma warning(disable:4514)  // The optimizer removed an inline function that is not called
//#pragma warning(disable:4702)  // unreachable code
//#pragma warning(disable:4127)  // to allow while(1)

#include <iostream>
#include <fstream>
#include <list>
#include <queue>
#include <string>
#include <cmath>

#include "Mytypes.h"	// all global variables
#include "Config.h"
#include "Cpoints.h"
#include "CONSTANTS.h"
#include "Dllist.h"
#include "Cimage.h"    // a simple image class
#include "Cvessel.h"
#include "Cvector.h"
#include "Template.h"  // a template class
#include "Extern.h"    // all global variables
#include "Ctree.h"
#include "Soma.h"

using namespace std;

extern CSomas* gTheSomas;
///////////////////////////////////////////////////////////////////////////////
//  					 Class CIntersectionPoint									  //
///////////////////////////////////////////////////////////////////////////////

// CTORS
CIntersectionPoint::CIntersectionPoint()
{
	Init();
}

CIntersectionPoint::CIntersectionPoint(CPoint& aPoint)
{
	Init();

	m_Point.m_iX = aPoint.m_iX;
	m_Point.m_iY = aPoint.m_iY;
	m_Point.m_iZ = aPoint.m_iZ;
}

CIntersectionPoint::CIntersectionPoint(CIntersectionPoint& aPoint)
{
	m_Point.m_iX = aPoint.m_Point.m_iX;
	m_Point.m_iY = aPoint.m_Point.m_iY;
	m_Point.m_iZ = aPoint.m_Point.m_iZ;
	m_iIndex = aPoint.m_iIndex;
	m_iID = 0;
	m_iSomaID = aPoint.m_iSomaID;
	m_Type = aPoint.m_Type;
	m_ProcessedFlag = aPoint.m_ProcessedFlag;

	memcpy(m_aVesselIDs, aPoint.m_aVesselIDs, sizeof(int) * MaxNumOfVesselIDs);
	memcpy(m_aPosition,
		aPoint.m_aPosition,
		sizeof(TopEndMiddleNeither) * MaxNumOfVesselIDs);
}
///////////////////
// Method: Init
//
// Initialization
void CIntersectionPoint::Init()
{
	m_Point.m_iX = m_Point.m_iY = m_Point.m_iZ = -1;
	m_iIndex = 0;
	m_iID = 0;
	m_iSomaID = 0;
	m_Type = VESSEL;
	m_ProcessedFlag = ORIGINAL;
	memset(m_aVesselIDs, 0, sizeof(int) * MaxNumOfVesselIDs);
	memset(m_aPosition, 0, sizeof(int) * MaxNumOfVesselIDs);
}

bool CIntersectionPoint::operator==(CIntersectionPoint& rhs)
{
	if (m_Point.m_iX != rhs.m_Point.m_iX ||
		m_Point.m_iY != rhs.m_Point.m_iY ||
		m_Point.m_iZ != rhs.m_Point.m_iZ)
	{
		return false;
	}
	else
	{
		int similar = 0;
		for (int i = 0; i < m_iIndex; i++)
		{
			if (m_aVesselIDs[i] == rhs.m_aVesselIDs[i])
			{
				similar++;
			}
		}
		if (similar == m_iIndex)
			return true;
		else
			return false;
	}
}

/////////////////
// Method: Update
// 
// If one of my particpating vessels is "oldId" change it to "newId"
//
void CIntersectionPoint::UpdateVesselID(int oldId, int newId, int iPosFlag)
{
	for (register int i = 0; i < m_iIndex; i++)
	{
		if (m_aVesselIDs[i] == oldId)
		{
			m_aVesselIDs[i] = newId;

			// check if we have to change the position as well
			if (iPosFlag)
			{
				if (m_aPosition[i] == OnTop)
					m_aPosition[i] = OnEnd;
				if (m_aPosition[i] == OnEnd)
					m_aPosition[i] = OnTop;
			}
		}
	}
}

///////////////////////
// Method: GetLocation
// return the location of the vessel with the given ID 
TopEndMiddleNeither CIntersectionPoint::GetLocation(int iID)
{
	for (register int i = 0; i < m_iIndex; i++)
	{
		if (m_aVesselIDs[i] == iID)
			return m_aPosition[i];
	}

	return Neither;
}
////////////////////
// Method: AddVessel
//
// Add the given vesselID to the list of my participating vessels
int CIntersectionPoint::AddVessel(int vesselID, TopEndMiddleNeither position)
{
	int result = 1;
	if (m_iIndex == MaxNumOfVesselIDs)
		result = 0;
	else
	{
		// make sure the vessel is not already added
		for (register int j = 0; j < m_iIndex; j++)
		{
			if (m_aVesselIDs[j] == vesselID)
				return 0;
		}

		m_aVesselIDs[m_iIndex] = vesselID;
		m_aPosition[m_iIndex] = position;
		m_iIndex++;
	}

	return result;
}
/////////////////////////////////
// Method: IsThereACommonVessel
// 
// return yes if this and the point at indx share a vessel, 0 otherwise
int CIntersectionPoint::IsThereACommonVessel(CIntersectionPoint* rhs)
{
	if (! rhs)
		return 0;

	for (register int i = 0; i < m_iIndex; i++)
	{
		for (register int j = 0; j < rhs->m_iIndex; j++)
		{
			if (m_aVesselIDs[i] == rhs->m_aVesselIDs[j])
				return 1;
		}
	}
	return 0;
}

///////////////////////////////////////
// Method: AddVessels
//
// Add the vessels of the given interseciton point to my list
void CIntersectionPoint::AddVessels(CIntersectionPoint* rhs)
{
	for (register int i = 0; i < rhs->m_iIndex; i++)
		AddVessel(rhs->m_aVesselIDs[i], rhs->m_aPosition[i]);
}

////////////////////////////
// Method: Print
//
// The argument is defaulted to cout
void CIntersectionPoint::Print(ostream& out)
{
	if (m_Type == VESSEL)
	{
		out << "Intersection " << m_iID << ": " << "at (" << m_Point.m_iX
			<< ", " << m_Point.m_iY << ", " << m_Point.m_iZ << ") "
			<< "between the segments: ";
		for (register int i = 0; i < m_iIndex; i++)
			out << m_aVesselIDs[i] << ", ";

		out << endl;
	}
	else
	{
		out << "Intersection " << m_iID << ": " << "at (" << m_Point.m_iX
			<< ", " << m_Point.m_iY << ", " << m_Point.m_iZ << ") "
			<< "Between Soma: " << gTheSomas->m_aData[m_iSomaID - 1].m_achID
			<< ", and segment " << m_aVesselIDs[0] << endl;
	}
}

///////////////////////////////////////////
// Method: RemoveVessel
//
// Remove the vessel with the give ID from my list. Return 1 upon sucess,
// 0, otherwise
int CIntersectionPoint::RemoveVessel(int ID)
{
	int result = 0;
	for (register int i = 0; i < m_iIndex; i++)
	{
		if (m_aVesselIDs[i] == ID)
		{
			if (i < m_iIndex - 1)
			{
				for (register int j = i; j < m_iIndex - 1; j++)
				{
					m_aVesselIDs[j] = m_aVesselIDs[j + 1];
					m_aPosition[j] = m_aPosition[j + 1];
				}
			}
			else
			{
				m_aVesselIDs[i] = 0;
				m_aPosition[i] = Neither;
			}
			result = 1;
			break;
		}
	}
	if (result)
		m_iIndex -= 1;

	return result;
}

////////////////////////////////////////
// Method: AdjustLocation
//
// IF the given point is located near the end of any of its vessels, snap it
// to such end
void CIntersectionPoint::AdjustLocation()
{
	if (m_Type == SOMA)
		return;

	CVessel* pVessel = NULL;
	CPoint* pPoint1 = NULL;
	CPoint* pPoint2 = NULL;
	int distance1 = 0;
	int distance2 = 0;
	int distThresh = 0;

	for (register int i = 0; i < m_iIndex; i++)
	{
		if (m_aPosition[i] == OnTop || m_aPosition[i] == OnEnd)
			continue;

		// the point is not on any of the ends.... is it close enough?
		pVessel = gTheVessels.GetVessel(m_aVesselIDs[i]);
		pPoint1 = pVessel->GetFirstPoint();
		distance1 = static_cast<int>(m_Point.FindDistance(pPoint1));
		pPoint2 = pVessel->GetLastPoint();
		distance2 = static_cast<int>(m_Point.FindDistance(pPoint2));
		distThresh = pVessel->GetWidthAtPoint(&m_Point);

		if (distance1 < distance2)
		{
			if (distance1 < distThresh)
			{
				m_Point = *pPoint1;
				m_aPosition[i] = OnTop;
			}
		}
		else if (distance2 < distThresh)
		{
			m_Point = *pPoint2;
			m_aPosition[i] = OnEnd;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// Method: IsThereAPath
//
// Find a path from This to the point with the given ID
int CIntersectionPoint::IsThereAPath(int iID)
{
	int iResult = 0;
	register int i, j;
	CVessel* pVessel = NULL;

	// if I am that point, then DONE
	if (m_iID == iID)
		return 1;

	// if I was processed before, return
	if (m_ProcessedFlag == YES)
		return 0;

	m_ProcessedFlag = YES;

	// find a path with each of the intersection points connected to my vessels
	for (i = 0; i < m_iIndex; i++)
	{
		// for each of my vessels
		if (m_aPosition[i] != Middle)
		{
			pVessel = gTheVessels.GetVessel(m_aVesselIDs[i]);

			// for each of the vessels intersection points
			for (j = 0; j < pVessel->m_iNumOfIntersectionPoints; j++)
			{
				iResult = gIntersectionPoints.m_apData[pVessel->m_aiMyIntersectionPoints[j] - 1]->
				IsThereAPath(iID);

				if (iResult)
					break;
			}
		}
		if (iResult)
			break;
	}

	return iResult;
}

///////////////////////////////////////////////////////////////////////////////
//  						 Class CIntersectionPoints							  //
///////////////////////////////////////////////////////////////////////////////


/////////////////////////
// Method: Init (private)
//
// Initialize data members
void CIntersectionPoints::Init()
{
	m_apData = new CIntersectionPoint * [BlockSize];
	memset(m_apData, 0, sizeof(CIntersectionPoint *) * BlockSize);
	m_iNumOfElements = 0;
	m_iSize = BlockSize;
	m_aDistanceArray = 0;
}

// DTOR
CIntersectionPoints::~CIntersectionPoints()
{
	for (register int i = 0; i < m_iNumOfElements; i++)
		delete m_apData[i];

	delete [] m_apData;
}

///////////////////
// Method:: Add
//
// Add a new point to my collection
// This operation requires the following:
// - check if the intersection points image already has a point in the given 
//   location.
// - If the point is adjacent to an end of the involved vessels, snap it to 
//   its end.
int CIntersectionPoints::Add(CIntersectionPoint* pElement)
{
	if (! pElement)
		return -1;

	// see if the point has been added earlier
	for (int i = 0; i < m_iNumOfElements; i++)
	{
		if (*m_apData[i] == *pElement)
			return -1;
	}

	// adjust the location of the point according to whether it lies
	// near the end of any of its vessels.
	if (pElement->m_Type != SOMA)
	{
		//		pElement->AdjustLocation();

		/*
								int x = pElement->m_Point.m_iX;
								int y = pElement->m_Point.m_iY;
								int pointID = 0;
								// if the intersection image already has a point in the given location
								if( (pointID = IntersectionPointsImage->data[y][x]) != 0)
								{
									// merge the information from both points into the old one
									m_apData[pointID-1]->AddVessels(pElement);
									return pointID;
								}
								*/
	}


	if (m_iNumOfElements < m_iSize)
	{
		m_apData[m_iNumOfElements] = new CIntersectionPoint(*pElement);
	}
	else
	{
		CIntersectionPoint** apTemp = new CIntersectionPoint*[m_iSize + BlockSize];
		m_iSize += BlockSize;
		memcpy(apTemp,
			m_apData,
			sizeof(CIntersectionPoint *) * m_iNumOfElements);
		apTemp[m_iNumOfElements] = new CIntersectionPoint(*pElement);

		delete [] m_apData;
		m_apData = apTemp;
	}
	m_iNumOfElements++;
	m_apData[m_iNumOfElements - 1]->m_iID = m_iNumOfElements;
	/*
		// mark the point in the image
		IntersectionPointsImage->
				data[pElement->m_Point.m_iY][pElement->m_Point.m_iX] = m_iNumOfElements;
	*/
	return m_iNumOfElements;
}

///////////////////
// Method:: Add
//
// Add a new element to my collection
int CIntersectionPoints::Add(CIntersectionPoint& Element)
{
	return Add(&Element);
}

//////////////////////////////////////////////////////////////////////
// METHOD: 	RemoveIntersectionPoint
//
// Remove the intersection point with the given ID from my collection
void CIntersectionPoints::Remove(int pointID)
{
	int index = - 1;
	int i;
	for (i = 0; i < m_iNumOfElements; i++)
	{
		if (m_apData[i]->m_iID == pointID)
		{
			index = i;
			break;
		}
	}

	delete m_apData[index];

	for (i = index; i < m_iNumOfElements - 1; i++)
		m_apData[i] = m_apData[i + 1];

	m_iNumOfElements--;
}

//////////////////////////////////////////
// Method:: UpdataIntersectionInformation
//
// For each point in my collection, if the point contains the vessel with the ID
// "oldID" as a participant, change it to "newID"
//
void CIntersectionPoints::UpdateVesselID(int oldID, int newID, int iPosFlag)
{
	for (register int i = 0; i < m_iNumOfElements; i++)
	{
		m_apData[i]->UpdateVesselID(oldID, newID, iPosFlag);
	}
}

///////////////////////////////////
// Method: CompleteVesselInfo
// Called By: MergeIntersectionPoints
//
// for each of the intersection points, find all of its adjacent vessels.
// In addition, for each of the vessels, determine if the point is 
// OnTop, OnEnd, or Middle point
void CIntersectionPoints::CompleteVesselInfo()
{
	TopEndMiddleNeither result = Neither;

	for (register int i = 0; i < m_iNumOfElements; i++)
	{
		if (m_apData[i]->m_Type == SOMA)
			continue;

		for (register int j = 0; j < gTheVessels.m_iNumOfElements; j++)
		{
			result = gTheVessels.m_apData[j]->GetPositionOfPoint(&m_apData[i]->m_Point);
			if (result == OnTop)
				m_apData[i]->AddVessel(j + 1, OnTop);
			else if (result == OnEnd)
				m_apData[i]->AddVessel(j + 1, OnEnd);
			else if (result == Middle)
				m_apData[i]->AddVessel(j + 1, Middle);
		}
	}
}

//////////////////////////////////////
// Method: MoreThanOneMiddleIntersection
//
// return 1 if the points with the given intersect more than one vessel in its middle
// and 0 otherwise
int CIntersectionPoints::MoreThanOneMiddleIntersection(int index1, int index2)
{
	int vesselIDs[MaxNumOfVesselIDs];
	CIntersectionPoint* pPoint1 = m_apData[index1];
	CIntersectionPoint* pPoint2 = m_apData[index2];
	int j = 0;
	int i;
	for (i = 0; i < pPoint1->m_iIndex; i++)
	{
		if (pPoint1->m_aPosition[i] == Middle)
			vesselIDs[j++] = pPoint1->m_aVesselIDs[i];
	}

	for (i = 0; i < pPoint2->m_iIndex; i++)
	{
		if (pPoint2->m_aPosition[i] == Middle)
			vesselIDs[j++] = pPoint2->m_aVesselIDs[i];
	}

	if (j > 1)
	{
		for (i = 0; i < j - 1; i++)
		{
			if (vesselIDs[i] != vesselIDs[i + 1])
				return 1;
		}
	}
	return 0;
}

////////////////////////////////////
// Method: FindLocationOfNewPoint
//
// Given two intersection points, find the location of a new intersection point
// representing the merging of both. Notice that the new intersection point
// must be a middle point to all vessels having either of the old intersection
// points as middle point.
int CIntersectionPoints::FindLocationOfNewPoint(int index1, int index2,
	CPoint& newPoint)
{
	newPoint.m_iX = m_apData[index1]->m_Point.m_iX +
		m_apData[index2]->m_Point.m_iX;
	newPoint.m_iY = m_apData[index1]->m_Point.m_iY +
		m_apData[index2]->m_Point.m_iY;
	newPoint.m_iX = (int) ((float) newPoint.m_iX * 0.5 + 0.5);
	newPoint.m_iY = (int) ((float) newPoint.m_iY * 0.5 + 0.5);

	int problemFlag = 0;
	int changeIndex = 0;
	CPoint ChangeArray[8];
	int i,j;
	for (i = -1; i <= 1; i++)
	{
		for (j = -1; j <= 1; j++)
		{
			if (i == 0 && j == 0)
				continue;

			ChangeArray[changeIndex].m_iX = i;
			ChangeArray[changeIndex].m_iY = j;
			changeIndex++;
		}
	}

	changeIndex = 0;
	while (1)
	{
		problemFlag = 0;
		for (i = 0; i < m_apData[index1]->m_iIndex; i++)
		{
			if (m_apData[index1]->m_aPosition[i] == Middle)
			{
				if (gTheVessels.m_apData[m_apData[index1]->m_aVesselIDs[i] - 1]
					->IsMemberOfCenterline(newPoint))
					continue;
				else // problem, modify the location of the new point
				{
					newPoint.m_iX += ChangeArray[changeIndex].m_iX;
					newPoint.m_iY += ChangeArray[changeIndex].m_iY;
					changeIndex++;
					problemFlag = 1;
				}
			}
		}

		if (changeIndex == 8)
			break;

		for (i = 0; i < m_apData[index2]->m_iIndex; i++)
		{
			if (m_apData[index2]->m_aPosition[i] == Middle)
			{
				if (gTheVessels.m_apData[m_apData[index2]->m_aVesselIDs[i] - 1]->IsMemberOfCenterline(newPoint))
					continue;
				else // problem, modify the location of the new point
				{
					newPoint.m_iX += ChangeArray[changeIndex].m_iX;
					newPoint.m_iY += ChangeArray[changeIndex].m_iY;
					changeIndex++;
					problemFlag = 1;
				}
			}
		}
		if (changeIndex == 8)
			break;

		if (! problemFlag)
			break;
	}
	if (problemFlag)
		return 0;
	else
		return 1;
}

///////////////////////////////////
// Method: MergeIntersectionPoints
// 
// For each two intersection points that are closer than the any of the 
// participating segment widths, merge them into one point.
//
void CIntersectionPoints::MergeIntersectionPoints()
{
	// before we start merging the intersection points, we must
	// find all of the vessels involved in every intersection point
	CompleteVesselInfo();

	int iFlag = 0;
	int i;
	for (i = 0; i < m_iNumOfElements; i++)
	{
		iFlag = 0;
		// To avoid processing the same point more than once, a flag is 
		// used. Also don't try to merge points that indicate
		// intersection with a soma
		if (m_apData[i]->m_ProcessedFlag == MERGED ||
			m_apData[i]->m_Type == SOMA)
			continue;

		int pointIndex = - 1;
		int Distance = FindClosestIntersectionPoint(i, pointIndex);
		if (pointIndex == -1)
			continue;

		if (pointIndex != -1)
		{
			int WidestSegment = 0;
			int Width1 = FindWidestSegment(i);
			int Width2 = FindWidestSegment(pointIndex);
			if (Width1 > Width2)
				WidestSegment = Width1;
			else
				WidestSegment = Width2;

			// the distance between the two points must be less than the widest
			// vessel for them to be merged
			if (2 * WidestSegment >= Distance)
			{
				// the two points must have vessels in common
				if (m_apData[i]->IsThereACommonVessel(m_apData[pointIndex]))
				{
					// if the two points intersect more than one vessel in the
					// middle, they are unmergable
					if (MoreThanOneMiddleIntersection(i, pointIndex))
						continue;

					// Find the location of the new point
					// notice that the new point must lie on all common vessels 
					// with either of the intersection points as middle points
					CPoint tempPoint;
					if (FindLocationOfNewPoint(i, pointIndex, tempPoint))
					{
						iFlag = 1;

						CIntersectionPoint tempIntPoint(tempPoint);
						tempIntPoint.m_ProcessedFlag = NEW;
						tempIntPoint.AddVessels(m_apData[i]);
						tempIntPoint.AddVessels(m_apData[pointIndex]);

						int newPointID = Add(tempIntPoint);	
						tempPoint = tempIntPoint.m_Point;

						// for each of the vessels, extend and/or update
						for (register int j = 0;
							j < tempIntPoint.m_iIndex;
							j++)
						{
							if (tempIntPoint.m_aPosition[j] == Middle)
							{
								gTheVessels.m_apData[tempIntPoint.m_aVesselIDs[j] - 1]->
								ReplaceIntersectionPoint(m_apData[pointIndex]->m_iID,
									newPointID);

								gTheVessels.m_apData[tempIntPoint.m_aVesselIDs[j] - 1]->
								ReplaceIntersectionPoint(m_apData[i]->m_iID,
									newPointID);
							}
							else if (tempIntPoint.m_aPosition[j] == OnTop)
							{
								gTheVessels.m_apData[tempIntPoint.m_aVesselIDs[j] - 1]->
								ExtendVesselCenter(&tempPoint, OnTop);

								gTheVessels.m_apData[tempIntPoint.m_aVesselIDs[j] - 1]->
								ReplaceIntersectionPoint(m_apData[pointIndex]->m_iID,
									newPointID);

								gTheVessels.m_apData[tempIntPoint.m_aVesselIDs[j] - 1]->
								ReplaceIntersectionPoint(m_apData[i]->m_iID,
									newPointID);
							}
							else if (tempIntPoint.m_aPosition[j] == OnEnd)
							{
								gTheVessels.m_apData[tempIntPoint.m_aVesselIDs[j] - 1]->
								ExtendVesselCenter(&tempPoint, OnEnd);

								gTheVessels.m_apData[tempIntPoint.m_aVesselIDs[j] - 1]->
								ReplaceIntersectionPoint(m_apData[pointIndex]->m_iID,
									newPointID);

								gTheVessels.m_apData[tempIntPoint.m_aVesselIDs[j] - 1]->
								ReplaceIntersectionPoint(m_apData[i]->m_iID,
									newPointID);
							}
							else
							{
								cout <<
									"\nCIntersectionPoints::MergeIntersectionPoints."
									<<
									"ERROR1. " <<
									endl;
								exit(0);
							}
						}

						// raise the flags for the old intersection points.
						m_apData[i]->m_ProcessedFlag = MERGED;
						m_apData[pointIndex]->m_ProcessedFlag = MERGED;
					} // ::FindTheLocationOfNewPoint
				} // :: IsThereACommonVessel
			}// :: if(WidestSegment >= Distance)
		}
	} // For loop

	// for all those points that were not mergable with any other point
	// extend their vessels to touch them
	TopEndMiddleNeither location = Neither;
	int vIndex = 0;
	for (i = 0; i < m_iNumOfElements; i++)
	{
		if (m_apData[i]->m_Type == VESSEL &&
			m_apData[i]->m_ProcessedFlag == ORIGINAL)
		{
			// for each vessel of this intersection point
			for (register int j = 0; j < m_apData[i]->m_iIndex; j++)
			{
				location = m_apData[i]->m_aPosition[j];
				vIndex = m_apData[i]->m_aVesselIDs[j] - 1;

				if (location == OnTop || location == OnEnd)
					gTheVessels.m_apData[vIndex]-> ExtendVesselCenter(&m_apData[i]->m_Point,
												   	location);

				if (location != Neither)
					gTheVessels.m_apData[vIndex]->AddIntersectionPoint(m_apData[i]->m_iID);
			}
		}
	}

	// remove all points with raised flags
	while (1)
	{
		for (i = 0; i < m_iNumOfElements; i++)
		{
			if (m_apData[i]->m_ProcessedFlag == MERGED)
			{
				Remove(m_apData[i]->m_iID);

				// Remove this ID from any vessel that is using it

				// since removing such point will change the value of
				// m_iNumOfElements, we must start all over again
				break;
			}
		}

		if (i == m_iNumOfElements)
			break;
	}
}

/////////////////////////////
// Method: FillDistanceArray
//
// Calculate the distance between each pair of points in my collection
//
void CIntersectionPoints::FillDistanceArray()
{
	// form a 2D Array of intersection Points, such that the entry[i][j] represents
	// the distance between the intersection points [i][j]
	m_aDistanceArray = new int * [m_iNumOfElements];
	int i;
	for (i = 0; i < m_iNumOfElements; i++)
		m_aDistanceArray[i] = new int[m_iNumOfElements];


	// calculate the square distance between any pair of intersection points
	for (i = 0; i < m_iNumOfElements; i++)
	{
		for (register int j = 0; j < m_iNumOfElements; j++)
		{
			if (i == j)
				m_aDistanceArray[i][j] = 0;
			else
			{
				int dx = m_apData[i]->m_Point.m_iX -
					m_apData[j]->m_Point.m_iX;
				int dy = m_apData[i]->m_Point.m_iY -
					m_apData[j]->m_Point.m_iY;

				m_aDistanceArray[i][j] = (dx * dx + dy * dy);
			}
		}
	}
}

////////////////////////////////////////
// Method:: FindClosestIntersectionPoint
//	
// Find the closest intersection point from the point with id1. 
// Return the distance, and set id2 to the found point id.
int CIntersectionPoints::FindClosestIntersectionPoint(int id1, int& id2)
{
	int minDistance = 999999;
	for (register int i = 0; i < m_iNumOfElements; i++)
	{
		if (i == id1)
			continue;

		if (m_apData[i]->m_ProcessedFlag == MERGED ||
			m_apData[i]->m_Type == SOMA)
			continue;

		int dx = m_apData[id1]->m_Point.m_iX - m_apData[i]->m_Point.m_iX;
		int dy = m_apData[id1]->m_Point.m_iY - m_apData[i]->m_Point.m_iY;

		int distance = dx* dx + dy* dy;
		if (distance < minDistance)
		{
			minDistance = distance;
			id2 = i;
		}
	}
	return (int) (sqrt((float) minDistance) + 0.5);
}

////////////////////////////////////////
// Method:: FindClosestIntersectionPoint
//	
// Find the closest intersection point from the given point
int CIntersectionPoints::FindClosestIntersectionPoint(CPoint& FromPoint,
	CIntersectionPoint** pClosestPoint)
{
	int minDistance = 999999;
	int dx, dy, dz;
	for (register int i = 0; i < m_iNumOfElements; i++)
	{
		dx = m_apData[i]->m_Point.m_iX - FromPoint.m_iX;
		dy = m_apData[i]->m_Point.m_iY - FromPoint.m_iY;
		dz = m_apData[i]->m_Point.m_iZ - FromPoint.m_iZ;

		int distance = dx* dx + dy* dy + dz* dz;
		if (distance < minDistance)
		{
			minDistance = distance;
			*pClosestPoint = m_apData[i];
		}
	}
	return (int) (sqrt((float) minDistance) + 0.5);
}

////////////////////////////////////////
// Method:: FindClosestIntersectionPoint
//	
// Find the closest intersection point from the given point
int CIntersectionPoints::FindClosestIntersectionPointExcluding(CPoint& FromPoint,
	CIntersectionPoint** pClosestPoint, int iID)
{
	int minDistance = 999999;
	int dx, dy, dz;
	for (register int i = 0; i < m_iNumOfElements; i++)
	{
		if (m_apData[i]->m_iID == iID)
			continue;

		dx = m_apData[i]->m_Point.m_iX - FromPoint.m_iX;
		dy = m_apData[i]->m_Point.m_iY - FromPoint.m_iY;
		dz = m_apData[i]->m_Point.m_iZ - FromPoint.m_iZ;

		int distance = dx* dx + dy* dy + dz* dz;
		if (distance < minDistance)
		{
			minDistance = distance;
			*pClosestPoint = m_apData[i];
		}
	}
	return (int) (sqrt((float) minDistance) + 0.5);
}

/////////////////////////////
// Method: GetCommonVesselID
//
// given the ids of two intersection points, id1 & id2, Find the id of the
// vessel common to both of these intersection points.
int CIntersectionPoints::GetCommonVesselID(int id1, int id2)
{
	if (m_apData[id1]->m_ProcessedFlag == MERGED &&
		m_apData[id2]->m_ProcessedFlag == MERGED)
		return -1;

	for (register int i = 0; i < m_apData[id1]->m_iIndex; i++)
	{
		for (register int j = 0; j < m_apData[id2]->m_iIndex; j++)
		{
			if (m_apData[id1]->m_aVesselIDs[i] ==
				m_apData[id2]->m_aVesselIDs[j])
				return m_apData[id1]->m_aVesselIDs[i];
		}
	}

	return -1;
}


//////////////////////////////////////////
// Method: FindWidestSegment
//
// Among the vessels touching the given intersection point, find the widest
// one. Notice that since the intersection point might be an extension, the 
// the given vessel might not have left and right boundary points at this
int CIntersectionPoints::FindWidestSegment(int id)
{
	int iWidth;
	iWidth= - 1;
	int iMaxWidth = - 1;
	for (register int i = 0; i < m_apData[id]->m_iIndex; i++)
	{
		int iVesselID = m_apData[id]->m_aVesselIDs[i];

		int iWidth = gTheVessels.m_apData[iVesselID - 1]->GetWidthAtPoint(&m_apData[id]->m_Point);
		if (iWidth != -1 && iWidth > iMaxWidth)
			iMaxWidth = iWidth;
	}
	return iMaxWidth;
}

////////////////////////////
// Method: Print
//
void CIntersectionPoints::Print(ostream& out)
{
	for (register int i = 0; i < m_iNumOfElements; i++)
		m_apData[i]->Print(out);
}


////////////////////////////
// Method: GetPoint
//
// return a pointer to the point with the given ID
CIntersectionPoint* CIntersectionPoints::GetPoint(int id)
{
	for (register int i = 0; i < m_iNumOfElements; i++)
	{
		if (m_apData[i]->m_iID == id)
			return m_apData[i];
	}
	return NULL;
}

//////////////////////////
// Method: ResetFlags
//
// For each intersection point, reset processed flag
void CIntersectionPoints::ResetFlags()
{
	for (register int i = 0; i < m_iNumOfElements; i++)
		m_apData[i]->ResetFlag();
}

///////////////////////////
//
// Remove all points that are located between vessels that has been already merged
void CIntersectionPoints::RemovePointsBetweenMergedVessels()
{
	int iIndex = 0;
	CIntersectionPoint** temp = new CIntersectionPoint*[m_iSize];

	for (register int i = 0; i < m_iNumOfElements; i++)
	{
		if (m_apData[i]->m_ProcessedFlag != MERGED)
		{
			temp[iIndex] = m_apData[i];
			iIndex++;
		}
	}
	delete [] m_apData;
	m_apData = temp;
	m_iNumOfElements = iIndex;
}


//////////////////////////////////////////////////////////
// is there a path from this intersection point to a soma
int CIntersectionPoints::IsThereAPathToASoma(int iID)
{
	int result = 0;
	register int i = 0;

	// Reset all Flags
	ResetFlags();

	for (i = 0; i < m_iNumOfElements; i++)
	{
		if (i == iID)
			continue;

		if (m_apData[i]->m_Type == SOMA)
		{
			if (m_apData[i]->IsThereAPath(iID))
			{
				result = 1;
				break;
			}
		}
	}

	ResetFlags();

	return result;
}
