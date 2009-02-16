#pragma warning(disable:4786)

#include <iostream>
#include <fstream>
#include <list>
#include <queue>
#include <cmath>

#include "CONSTANTS.h"
#include "Mytypes.h"
#include "Config.h"
#include "Cpoints.h"
#include "Dllist.h"
#include "Cimage.h"
#include "C3dimage.h"
#include "Cvector.h"
#include "Cvessel.h"
#include "Template.h"
#include "Extern.h"

#define DISTANCE_THRESHOLD 10

/// used below in GetNumOfOppositPairs
int giMinDirDistance = 3;

//////////////////////////////////////////////////////////////////////////
// Function: FindTwoMaxima
// type: support function for seed point verification
// 
// given the array of boundary points, find two maxima such and point to them
// using the given arguments in the specified order
void FindTwoMaxima(CPoint aPoints[], CPoint** pPoint1, CPoint** pPoint2)
{
	register int i = 0;
	int iIndex1, iMax1 = 0;
				int iIndex2, iMax2 = 0;
	int iBefore, iAfter;

	*pPoint1 = 0;
	*pPoint2 = 0;

	for (i = 0; i < NumOfDirections; i++)
	{
		iBefore = i - 1;
		iAfter = i + 1;
		if (i == 0)
			iBefore = NumOfDirections - 1;
		if (i == NumOfDirections - 1)
			iAfter = 0;
		if (aPoints[i].m_iValue > aPoints[iAfter].m_iValue &&
			aPoints[i].m_iValue > aPoints[iBefore].m_iValue)
		{
			if (i < NumOfDirections / 2)
			{
				if (aPoints[i].m_iValue > iMax1)
				{
					iMax1 = aPoints[i].m_iValue;
					iIndex1 = i;
				}
			}
			else
			{
				if (aPoints[i].m_iValue > iMax2)
				{
					iMax2 = aPoints[i].m_iValue;
					iIndex2 = i;
				}
			}
		}
	}

	if (iMax1 && iMax2)
	{
		*pPoint1 = &aPoints[iIndex1];
		*pPoint2 = &aPoints[iIndex2];
	}
}

/////////////////////////////////////////////////////////////////////////
// Function: IsParallel
// type: Support function for seed point verification
//
bool IsParallel(int dir1, int dir2)
{
	if (DirDistance(dir1, dir2) < 3)
		return true;
	else
		return false;
}

/////////////////////////////////////////////////////////////////////////
// Function: IsOpposite
// type: Support function for seed point verification
//
bool IsOpposite(int dir1, int dir2)
{
	int diff = DirDistance(dir1, dir2);
	if (diff >= (NumOfDirections / 4 - 2) && diff <= (NumOfDirections / 4 + 2))
		return true;
	else
		return false;
}

/////////////////////////////////////////////////////////////////////////
// Function: IsClose
// type: Support function for seed point verification
//
bool IsClose(CPoint* pPoint1, CPoint* pPoint2)
{
	int XDiff = pPoint1->m_iX - pPoint1->m_iX;
	int YDiff = pPoint1->m_iY - pPoint1->m_iY;
	int ZDiff = pPoint1->m_iZ - pPoint1->m_iZ;

	int distance = XDiff* XDiff + YDiff* YDiff + ZDiff* ZDiff;

	if (distance < DISTANCE_THRESHOLD)
		return true;
	else
		return false;
}

// Given an array of directions, calculate the number of "almost" opposite
// pairs in the array.
int GetNumOfOppositePairs(int anArray[], int arraySize)
{
	int count = 0;
	int index = - 1;
	int dirDist;
	int minDirDist;
	int tempArray[NumOfDirections];

	memcpy(tempArray, anArray, NumOfDirections * sizeof(int));

	for (register int i = 0; i < arraySize - 1; i++)
	{
		minDirDist = NumOfDirections;

		for (register int j = i + 1; j < arraySize; j++)
		{
			if (tempArray[i] != -1 && tempArray[j] != -1)
			{
				dirDist = abs(DirDistance(anArray[i], anArray[j]) -
						  	NumOfDirections /
						  	2);
				if (dirDist < minDirDist)
				{
					minDirDist = dirDist;
					index = j;
				}
			}
		}

		// if we found one, mark these points as a pair (by deleting them)
		// and increment the count
		if (minDirDist < giMinDirDistance)
		{
			tempArray[i] = -1;
			tempArray[index] = -1;
			count++;
		}
	}
	return count;
}

////////////////////////////////////////////////////////////////////
// Function: ValidSeedCandidate
// CalledBy: FindSeedPoints
// Calls   : -
// Purpose :
//   To avoid having many initial that are adjacent to each other
//
// Logic:
//   Use the 2D array of seed point pointers to check if we already have another
//   seed point in the immediate neighborhood of the given point. If we do, 
//   keep the best one among all of such points
bool AddSeedCandidate(int x, int y, int dir, int value)
{
	int WindowSize = 2;
	bool result = true;

	int giROWS = The3DImage->m_iRows;
	int giCOLS = The3DImage->m_iCols;

	for (register int i = y - WindowSize; i <= y + WindowSize; i++)
	{
		for (register int j = x - WindowSize; j <= x + WindowSize; j++)
		{
			if (i < 0 || i >= giROWS || j < 0 || j >= giCOLS)
				continue;

			// if we alread have a point in this window,
			if (gapArrayOfSeedPoints[i][j])
			{
				if (gapArrayOfSeedPoints[i][j]->m_iValue > value)
					result = false;
				else
				{
					delete gapArrayOfSeedPoints[i][j];
					gapArrayOfSeedPoints[i][j] = NULL;
				}
			}
		}
	}

	// if the current point is the best, add it to the array
	if (result)
	{
		gapArrayOfSeedPoints[y][x] = new CPoint(x, y, 0, dir, 0, value);
	}
	return result;
}
////////////////////////////////////////////////////////////
// Function: FindSeedPoints
//
// Find seed points for tracking the given image and store them into
// the given queue. Notice that the seed points need to pass through
// a verification process before they can be considered as valid seed
// points. We do not do this here, since many of such 
// points will not be used, thus avoiding unnecessary computations.
// 
// Logic: seed points are defined as local minima points along 
// horisontal and vertical profiles. To reduce the large number of
// resulting points, only the smallest minima from an interval 
// (i.e. grid_size) is selected.
//
// return the number of points found
//int FindSeedPoints(C3DImage &anImage, CQueue<CPoint> &aQueue, int SliceNum)
int FindSeedPoints2(CImage& anImage)
{
	int giMARGIN2 = 3;
	register int x, y, i;
	register int max, temp, result;
	CPoint aPoint;
	int count = 0;

	unsigned char * *ImageData = anImage.data;

	int giROWS = The3DImage->m_iRows;
	int giCOLS = The3DImage->m_iCols;
	int giGRIDSIZE = gConfig.GetGridSpacing();

	// search along vertical profiles
	// from the left boundary to the right boundary
	for (x = giMARGIN2; x < giCOLS - giMARGIN2; x += giGRIDSIZE)
	{
		// from the top boundary to the bottome boundary
		for (y = giMARGIN2; y < giROWS - giMARGIN2; y++)
		{
			result = 1;
			// the average value of the point (without division)
			max = ImageData[y][x] + ImageData[y][x - 1] + ImageData[y][x + 1];

			// Is this value (i.e. min) the minimum for "HALFGRID" pixels
			// ahead? 
			for (i = y + 1 ; i < (y + giGRIDSIZE); i++)
			{
				if (i >= giROWS - giMARGIN2)
					break;

				// the average of another point along the profile
				temp = ImageData[i][x] +
					ImageData[i][x - 1] +
					ImageData[i][x + 1];

				if (max < temp)
				{
					result = 0;
					break;
				}
			}

			// YES. 
			// add it to the queue and advance x to be HALFGRID points
			// ahead. 
			if (result && max > 3.0 * giMinSeedPixelValue)
			{
				// add the seed point if and only if it satisfies the following
				// test
				if (AddSeedCandidate(x, y, NumOfDirections / 4, max))
					count++;

				// update the average inside pixel and outside pixel
				glSumOfSeedPointValues += ImageData[y][x];
				glNumOfSeedPoints++;

				y += giGRIDSIZE;
			}
			// NO.
			// consider the next point in line
			else
			{
				y = i;
			}
		}
	}

	// search along horizontal profiles.
	// from top boundary to lower one
	for (y = giMARGIN2; y < giROWS - giMARGIN2; y += giGRIDSIZE)
	{
		//////////////////////////////////////////
		// from the left boundary to the right one

		// if the point is a global minimum in the specified interval
		// then add it to the collection of seed points
		for (x = giMARGIN; x < giCOLS - giMARGIN2; x++)
		{
			// if one of the neighboring pixels is set to 255, then it belongs to a
			// a vessel, skip it and a grid size around it
			if (ImageData[y][x] == 255 ||
				ImageData[y - 1][x] == 255 ||
				ImageData[y + 1][x] == 255)
			{
				x += giGRIDSIZE;
				continue;
			}

			result = 1;
			// the average value of the point (without division)
			max = ImageData[y][x] + ImageData[y - 1][x] + ImageData[y + 1][x];

			// Is this value (i.e. min) the minimum for "HALFGRID" pixels
			// to come? 
			for (i = x + 1 ; i < x + giGRIDSIZE; i++)
			{
				if (i >= giCOLS - giMARGIN2)
					break;

				// the average of another point along the profile
				temp = ImageData[y][i] +
					ImageData[y - 1][i] +
					ImageData[y + 1][i];

				if (max < temp)
				{
					result = 0;
					break;
				}
			}

			// YES. 
			// add it to the queue and advance x to be HALFGRID points
			// ahead.
			if (result && max > 3.0 * giMinSeedPixelValue)
			{
				// add the seed point if and only if it satisfies the following
				// test
				if (AddSeedCandidate(x, y, 0, max))
					count++;
				/*
																		CPoint temp(x, y, SliceNum, 0);
																		aQueue.Push(temp);				
																		count++;
														*/				
				// update the average inside pixel and outside pixel
				glSumOfSeedPointValues += ImageData[y][x];
				glNumOfSeedPoints++;

				x += giGRIDSIZE;
			}
			// NO.
			// consider the next point in line
			else
			{
				x = i;
			}
		}
	}

	return count;
}

// these are defined here because MSC++ does not support local functions!
int Average[NumOfDirections];
int compare1(const void* arg1, const void* arg2)
{
	int result = 0;
	if (Average[*(int *) arg1] > Average[*(int *) arg2])
		result = 1;
	else if (Average[*(int *) arg1] < Average[*(int *) arg2])
		result = -1;
	return result;
}
int compare2(const void* arg1, const void* arg2)
{
	int result = 0;
	if (Average[*(int *) arg1] < Average[*(int *) arg2])
		result = 1;
	else if (Average[*(int *) arg1] > Average[*(int *) arg2])
		result = -1;
	return result;
}
//////////////////////////////////////////////////////////////////////////
// Function: VerifySeedPoint
// Input:  A seed point
// Output: Yes/No depending if the point qualifies as an seed point for
//   the tracking process.
//
// Logic:
// 
// A valid seed point is one that:
// 1. Has not been already covered by another tracking cycle.
// 2. Produces two maxima and two minima template responses such that:
//    2.1 the maxima are (almost) opposite to each other, and
//    2.2 the minima are (almost) opposite to the maxima directions.
//		2.3 the maximas and minima are almost perpendicular to each other
// 3. Satisfy minimum response requirements.
//



/////////////////////////////////////////////////////////////////////////////
// Verify that the array of responses (a response for each direction) 
// represents a valid response array according to the old 2D verification rules
bool VerifySeedResponses2D(int* Response, int* dirVotes)
{
	bool result = true;

	//	return true;

	int min[NumOfDirections] =
	{
		0
	};
	int max[NumOfDirections] =
	{
		0
	};

	// variables to be used by subsequent processes
	int dirBefore, dirAfter;

	int minCount = 0, maxCount = 0;

	// weighted average of the responses. This results in responses that
	// are less sensitive to noise. 
	// Notice that we are using the global "Average" array
	int dir, hDir, vDir, hDirBefore, hDirAfter, vDirBefore, vDirAfter;
	int lastIndex = NumOfDirections* NumOfDirections - 1;
	for (register int i = 0; i < NumOfDirections; i++)
	{
		for (register int j = 0; j < NumOfDirections; j++)
		{
			if (i == 0)
				dirBefore = lastIndex;
			if (i == NumOfDirections * NumOfDirections - 1)
				Average[i] = Response[i];

			hDir = i / NumOfDirections;
			vDir = i % NumOfDirections;

			hDirBefore = DirectionMinus(hDir, 1);
			hDirAfter = DirectionPlus(hDir, 1);
			vDirBefore = DirectionMinus(vDir, 1);
			vDirAfter = DirectionPlus(vDir, 1);

			dir = hDir * NumOfDirections + vDir;
			dirAfter = hDirAfter * NumOfDirections + vDirAfter;
			dirBefore = hDirAfter * NumOfDirections + vDirAfter;

			Average[dir] = Response[dirBefore] +
				Response[dirAfter] +
				Response[dir] +
				Response[dir];
		}
	}

	// calculate the number of maxima and minima in "Average"
	for (dir = 0; dir < NumOfDirections; dir++)
	{
		dirBefore = DirectionMinus(dir, 1);
		dirAfter = DirectionPlus(dir, 1);

		if (Average[dir] > Average[dirBefore] &&
			Average[dir] >= Average[dirAfter])
		{
			max[maxCount] = dir;
			++maxCount;
		}
		else if (Average[dir] < Average[dirBefore] &&
			Average[dir] <= Average[dirAfter])
		{
			min[minCount] = dir;
			++minCount;
		}
	}

	////////////////////////////////////////
	// 2. Do we have at least two maxima and two minima?
	if (maxCount < 2)//|| minCount < 2)
	{
		result = false;
	}
	else
	{
		// sort them so that the valid responses are listed first
		qsort((void *) max, maxCount, sizeof(int), compare2);
		qsort((void *) min, minCount, sizeof(int), compare1);

		int max0 = Average[max[0]];
		int max1 = Average[max[1]];
		int min0 = 2 * Average[min[0]];
		int min1 = 2 * Average[min[1]];

		// the threshold is set based on the minimum detectable edge.
		// for a template of length 12
		if (min0 > max0 || min0 > max1 || min1 > max0 || min1 > max1)
			result = false;


		/////////////////////////////////////////////////////////////
		// 3. all the above conditions were satisfied. perform the last
		//    threshold test
		if (result)
		{
			// if the number of maximas (minimas) is not even, only use those
			// the strongest (weakest) responses.
			if (maxCount % 2 != 0)
				maxCount--;

			if (minCount % 2 != 0)
				minCount--;

			// 2.1 If we have two maximas, then we expect them to be almost opposite
			//     If we have four maximas, then we expect to have two pairs of almost
			//     opposite points, so perform such check
			//
			// count how many opposite pairs we have in each array
			int NumOfOppositePairs = GetNumOfOppositePairs(max, maxCount);
			if (NumOfOppositePairs < 1)//(maxCount/2))
				result = false;

			//
			// vote for the desired initial direction
			if (result)
				dirVotes[max[0]]++;
		}
	}
	return result;
}

deque<CPoint> boundaries;

/////////////////////////////////////////////////////////////
// Function: VerifySeedPoint3D2
//
// Verify that the given point is a valid seed point.
bool VerifySeedPoint3D2(CPoint* aPoint)
{
	bool result = true;

	/////////////////////
	// Notice that the given point does not have its z-component assigned yet
	// becuase it was obtained from the projection image. So find the plane
	// with the highest average and assign it as the possible z-value
	int y = aPoint->m_iY;
	int x = aPoint->m_iX;
	int z = 0;
	int sum, Max = 0;
	
	int giSLICES = The3DImage->m_iSlices;
	register int i;
	for (i = 0; i < giSLICES; i++)
	{
		sum = The3DImage->data[i][y - 1][x - 1] +
			The3DImage->data[i][y - 1][x] +
			The3DImage->data[i][y - 1][x + 1] +
			The3DImage->data[i][y][x - 1] +
			The3DImage->data[i][y][x] +
			The3DImage->data[i][y][x + 1] +
			The3DImage->data[i][y + 1][x - 1] +
			The3DImage->data[i][y + 1][x] +
			The3DImage->data[i][y + 1][x + 1];

		if (sum > Max)
		{
			Max = sum;
			z = i;
		}
	}

	aPoint->m_iZ = z;
	
	

	// an array of votes for each direction
	//int HDirVotes[NumOfDirections];
	int* DirVotes = new int [NumOfDirections* NumOfDirections];
	int* Response = new int [NumOfDirections* NumOfDirections];
	//memset(HDirVotes, 0, sizeof(int)*NumOfDirections);
	memset(DirVotes, 0, sizeof(int) * NumOfDirections * NumOfDirections);
	memset(Response, 0, sizeof(int) * NumOfDirections * NumOfDirections);

	// since we call this function before we perform any tracking, 
	// we don't have to see if the candidate point lies on an already 
	// tracked vessel

	CPoint tempPoint(*aPoint);
	CPoint HLeftPoint, HRightPoint;
	CPoint VLeftPoint, VRightPoint;
	CPoint BestHLeftPoint, BestVLeftPoint;
	CPoint BestHRightPoint, BestVRightPoint;

	int BestResponse = 0;
	int index = 0;

	int BestForegroundEstimate = 0;
	int BestBackgroundEstimate = 0;

	int VDirections[9];
	PrepareDirectionsArray(VDirections, 9, 0);
	register int Hdir, Vdir, v;

	// notice that since we limit starting points to those withing +/- 45
	// vertical degrees. This usually produces more stable tracking results.

	// for each possible direction, record the maximum left and right
	// template responses in that direction
	for (Hdir = 0; Hdir < NumOfDirections; Hdir++)
	{
		for (v = 0; v < 9; v++)
		{
			Vdir = VDirections[v];
			Response[index] = gHLeftTemplatesArray[Hdir][Vdir]->CalculateMaxResponse(&tempPoint,
																	& HLeftPoint) +
				gHRightTemplatesArray[Hdir][Vdir]->CalculateMaxResponse(&tempPoint,
												   	& HRightPoint) +
				gVLeftTemplatesArray[Hdir][Vdir]->CalculateMaxResponse(&tempPoint,
												  	& VLeftPoint) +
				gVRightTemplatesArray[Hdir][Vdir]->CalculateMaxResponse(&tempPoint,
												   	& VRightPoint);

			if (Response[index] > BestResponse)
			{
				BestResponse = Response[index];
				BestHLeftPoint = HLeftPoint;
				BestHRightPoint = HRightPoint;
				BestVLeftPoint = VLeftPoint;
				BestVRightPoint = VRightPoint;

				BestForegroundEstimate = gHLeftTemplatesArray[Hdir][Vdir]->m_iForegroundPixelValue +
					gHRightTemplatesArray[Hdir][Vdir]->m_iForegroundPixelValue +
					gVLeftTemplatesArray[Hdir][Vdir]->m_iForegroundPixelValue +
					gVRightTemplatesArray[Hdir][Vdir]->m_iForegroundPixelValue;

				BestBackgroundEstimate = gHLeftTemplatesArray[Hdir][Vdir]->m_iBackgroundPixelValue +
					gHRightTemplatesArray[Hdir][Vdir]->m_iBackgroundPixelValue +
					gVLeftTemplatesArray[Hdir][Vdir]->m_iBackgroundPixelValue +
					gVRightTemplatesArray[Hdir][Vdir]->m_iBackgroundPixelValue;
			}
			index++;
		}
	}

	BestForegroundEstimate /= 4;
	BestBackgroundEstimate /= 4;

	int giUsedTemplateLength = gConfig.GetMinimumTemplateLength();
	int threshold = 0;
	if (gf3DStdDev > 3.0)
		threshold = 3 * 3 * giUsedTemplateLength;
	else if (gf3DStdDev > 0.0)
		threshold = gf3DStdDev * 3.0 * giUsedTemplateLength;
	else
		threshold = 3.0 * giUsedTemplateLength;

	// we assume starting points must have at least a contrast of 3
	if (BestResponse < threshold)
		result = false;

	if (BestHLeftPoint.m_iValue < threshold ||
		BestHRightPoint.m_iValue < threshold ||
		BestVRightPoint.m_iValue < threshold ||
		BestVLeftPoint.m_iValue < threshold)
		result = false;

	if (result && VerifySeedResponses2D(Response, DirVotes))
	{
		aPoint->m_iValue = BestResponse / 4;

		aPoint->m_iHDir = BestHRightPoint.m_iHDir;
		if (BestHLeftPoint.m_iValue > BestHRightPoint.m_iValue)
			aPoint->m_iHDir = BestHLeftPoint.m_iHDir;

		aPoint->m_iVDir = BestVRightPoint.m_iVDir;
		if (BestVLeftPoint.m_iValue > BestVRightPoint.m_iValue)
			aPoint->m_iVDir = BestVLeftPoint.m_iVDir;

		GetCenterLocation(*aPoint,
			BestHRightPoint,
			BestHLeftPoint,
			BestVRightPoint,
			BestVLeftPoint);

		aPoint->m_iValue = BestResponse;

		double Xdiff = BestHLeftPoint.m_iX - BestHRightPoint.m_iX;
		double Ydiff = BestHLeftPoint.m_iY - BestHRightPoint.m_iY;
		double Zdiff = BestHLeftPoint.m_iZ - BestHRightPoint.m_iZ;

		double width = std::sqrt(Xdiff* Xdiff + Ydiff* Ydiff + Zdiff* Zdiff);

		// the following two lines were added for CANCER images
		float h_width = BestHLeftPoint.FindDistance(&BestHRightPoint);
		float v_width = BestVLeftPoint.FindDistance(&BestVRightPoint);
		gfHWidth += h_width;
		gfVWidth += v_width;
		aPoint->m_fHWidth = h_width;
		aPoint->m_fVWidth = v_width;


		gfWidthSum += width;
		giNumOfWidthSumMembers++;
	}
	else
		result = false;

	// 1. if the point is a valid point, update the foreground and background histograms
	// 2. add the point to the list of seed points
	if (result == true)
	{
		deque<CPoint> points; 
		points.push_back(*aPoint);
		points.push_back(BestHRightPoint);
		points.push_back(BestHLeftPoint);
		points.push_back(BestVRightPoint);
		points.push_back(BestVLeftPoint);
		Seeds.push_back(points);

		gaiForegroundHistogram[BestForegroundEstimate]++;
		gaiBackgroundHistogram[BestBackgroundEstimate]++;
	}

	delete [] Response;
	delete [] DirVotes;

	return result;
}

bool VerifySeedPoint3D3(CPoint* aPoint)
{
	bool result = true;

	/////////////////////
	// Notice that the given point does not have its z-component assigned yet
	// becuase it was obtained from the projection image. So find the plane
	// with the highest average and assign it as the possible z-value
	int y = aPoint->m_iY;
	int x = aPoint->m_iX;
	int z = 0;
	int sum, Max = 0;
	
	int giSLICES = The3DImage->m_iSlices;
	register int i;
	for (i = 0; i < giSLICES; i++)
	{
		sum = The3DImage->data[i][y - 1][x - 1] +
			The3DImage->data[i][y - 1][x] +
			The3DImage->data[i][y - 1][x + 1] +
			The3DImage->data[i][y][x - 1] +
			The3DImage->data[i][y][x] +
			The3DImage->data[i][y][x + 1] +
			The3DImage->data[i][y + 1][x - 1] +
			The3DImage->data[i][y + 1][x] +
			The3DImage->data[i][y + 1][x + 1];

		if (sum > Max)
		{
			Max = sum;
			z = i;
		}
	}

	aPoint->m_iZ = z;	

	int Response;

	CPoint tempPoint(*aPoint);
	CPoint HLeftPoint, HRightPoint;
	CPoint VLeftPoint, VRightPoint;
	CPoint BestHLeftPoint, BestVLeftPoint;
	CPoint BestHRightPoint, BestVRightPoint;

	BestHLeftPoint.m_iValue = BestHRightPoint.m_iValue = 0;
	BestVLeftPoint.m_iValue = BestVRightPoint.m_iValue = 0;
	int BestResponse = 0;

	int giShiftDistance = gConfig.GetMaximumShiftDistance();
	int ShiftDistance = giShiftDistance;
	
	int BestForegroundEstimate = 0;
	int BestBackgroundEstimate = 0;

	int VDirections[9];
	PrepareDirectionsArray(VDirections, 9, 0);
	register int Hdir, Vdir, v;

	// notice that since we limit starting points to those withing +/- 45
	// vertical degrees. This usually produces more stable tracking results.

	// for each possible direction, record the maximum left and right
	// template responses in that direction
	for (Hdir = 0; Hdir < NumOfDirections; Hdir++)
	{
		for (v = 0; v < 9; v++)
		{
			Vdir = VDirections[v];
			Response = gHLeftTemplatesArray[Hdir][Vdir]->CalculateMaxResponse(&tempPoint,
														 	& HLeftPoint) +
				gHRightTemplatesArray[Hdir][Vdir]->CalculateMaxResponse(&tempPoint,
												   	& HRightPoint) +
				gVLeftTemplatesArray[Hdir][Vdir]->CalculateMaxResponse(&tempPoint,
												  	& VLeftPoint) +
				gVRightTemplatesArray[Hdir][Vdir]->CalculateMaxResponse(&tempPoint,
												   	& VRightPoint);

			if (Response > BestResponse)
			{
				BestResponse = Response;
				BestHLeftPoint = HLeftPoint;
				BestHRightPoint = HRightPoint;
				BestVLeftPoint = VLeftPoint;
				BestVRightPoint = VRightPoint;

				BestForegroundEstimate = gHLeftTemplatesArray[Hdir][Vdir]->m_iForegroundPixelValue +
					gHRightTemplatesArray[Hdir][Vdir]->m_iForegroundPixelValue +
					gVLeftTemplatesArray[Hdir][Vdir]->m_iForegroundPixelValue +
					gVRightTemplatesArray[Hdir][Vdir]->m_iForegroundPixelValue;

				BestBackgroundEstimate = gHLeftTemplatesArray[Hdir][Vdir]->m_iBackgroundPixelValue +
					gHRightTemplatesArray[Hdir][Vdir]->m_iBackgroundPixelValue +
					gVLeftTemplatesArray[Hdir][Vdir]->m_iBackgroundPixelValue +
					gVRightTemplatesArray[Hdir][Vdir]->m_iBackgroundPixelValue;
			}
		}
	}

	BestForegroundEstimate /= 4;
	BestBackgroundEstimate /= 4;

	int giUsedTemplateLength = gConfig.GetMinimumTemplateLength();
	int threshold = 0;
	if (gf3DStdDev > 3.0)
		threshold = 3 * 3 * giUsedTemplateLength;
	else if (gf3DStdDev > 0.0)
		threshold = gf3DStdDev * 3.0 * giUsedTemplateLength;
	else
		threshold = 3.0 * giUsedTemplateLength;

	// we assume starting points must have at least a contrast of 3
	if (BestResponse < threshold)
		result = false;

	if (BestHLeftPoint.m_iValue < threshold ||
		BestHRightPoint.m_iValue < threshold ||
		BestVRightPoint.m_iValue < threshold ||
		BestVLeftPoint.m_iValue < threshold)
		result = false;

	// revert to original value
	giShiftDistance = ShiftDistance;

	if (result)
	{
		aPoint->m_iValue = BestResponse / 4;
		aPoint->m_iHDir = BestHRightPoint.m_iHDir;
		if (BestHLeftPoint.m_iValue > BestHRightPoint.m_iValue)
			aPoint->m_iHDir = BestHLeftPoint.m_iHDir;

		aPoint->m_iVDir = BestVRightPoint.m_iVDir;
		if (BestVLeftPoint.m_iValue > BestVRightPoint.m_iValue)
			aPoint->m_iVDir = BestVLeftPoint.m_iVDir;

		// fix the position of the starting point
		aPoint->m_iX = (int)
			((BestHLeftPoint.m_iX + BestHRightPoint.m_iX) / 2.0 + 0.5);
		aPoint->m_iY = (int)
			((BestHLeftPoint.m_iY + BestHRightPoint.m_iY) / 2.0 + 0.5);
		aPoint->m_iZ = (int)
			((BestVLeftPoint.m_iZ + BestVRightPoint.m_iZ) / 2.0 + 0.5);

		int Xdiff = BestHLeftPoint.m_iX - BestHRightPoint.m_iX;
		int Ydiff = BestHLeftPoint.m_iY - BestHRightPoint.m_iY;
		int Zdiff = BestHLeftPoint.m_iZ - BestHRightPoint.m_iZ;

		float width = sqrt(static_cast<float>(Xdiff* Xdiff + Ydiff* Ydiff + Zdiff* Zdiff));

		gfWidthSum += width;
		giNumOfWidthSumMembers++;
	}
	else
		result = false;

	// if the point is a valid point, update the 
	// foreground and background histograms
	if (result == true)
	{
		gaiForegroundHistogram[BestForegroundEstimate]++;
		gaiBackgroundHistogram[BestBackgroundEstimate]++;
	}

	return result;
}

bool SearchForSeedAlongA_Line(CPoint& aPoint, CVector* aVector,
	CPoint& minPoint)
{
	cout <<
		"This function : SearchForSeedAlongA_Line Should've not been called "
		<<
		endl;
	exit(0);
	return false;
	/* 3-3-1999
	bool result = false;
	int min, temp;
	// determine that the whole line is within the boundaries of the image 
	// being tracked (i.e. TheTargetImage) 
	int x1 = aPoint.m_iX + aVector->IndexC[0];
	int x2 = aPoint.m_iX + aVector->IndexC[MaxVectorLength-1];
	int y1 = aPoint.m_iY + aVector->IndexR[0];
	int y2 = aPoint.m_iY + aVector->IndexR[MaxVectorLength-1];
	if( (x1 > giMARGIN && x1 < giCOLS-giMARGIN) &&
		 (x2 > giMARGIN && x2 < giCOLS-giMARGIN) &&
		 (y1 > giMARGIN && y1 < giROWS-giMARGIN) &&
		 (y2 > giMARGIN && y2 < giROWS-giMARGIN) )
	{
		// the orthogonal directions
		int dir1 = DirectionPlus(aPoint.dir, NumOfDirections/4);
		int dir2 = DirectionMinus(aPoint.dir, NumOfDirections/4);
		// a pointer to the image pixel where the vector is starting
		unsigned char *baseData  = &(The3DImage->data[aPoint.m_iZ][y1][x1]);
		unsigned char *baseData1 = baseData + VectorsArray[dir1]->data[1];
		unsigned char *baseData2 = baseData + VectorsArray[dir2]->data[1];
		for(register int i = 0; i < MaxVectorLength; i++)
		{
			result = true;
			min    = *(baseData  + aVector->data[i]) + 
					  *(baseData1 + aVector->data[i]) +
					  *(baseData2 + aVector->data[i]);
			// Is this value (i.e. min) the minimum until the end of the line
			// ahead? 
			for(register int j = i + 1 ; j < MaxVectorLength; j++)
			{
				// the average of another point along the profile
				temp = *(baseData  + aVector->data[j]) + 
						*(baseData1 + aVector->data[j]) +
						*(baseData2 + aVector->data[j]);
				if(min > temp)
				{
					result = false;
					break;
				}
			}
				
			// YES. 
			// add it to the queue and advance x to be HALFGRID points
			// ahead. 
			if(result)
			{
				minPoint.m_iX = aPoint.m_iX + aVector->IndexC[i];
				minPoint.m_iY = aPoint.m_iY + aVector->IndexR[i];
				minPoint.m_iZ = aPoint.m_iZ;
				minPoint.dir = aPoint.dir;
				break;
			}
			// NO.
			// consider the next point in line
			else
			{
				i = j;
			}
		}
	}
	return result;
	*/
}

///////////////////////////////////////////////////////////////////////////////
// Function: FindAndAddMoreSeedPoints
//
// Purpose: 
//		To locate potential seed points that are encountered while tracking
// Logic:
//		Starting from the given point, search along the two lines extending in 
//  	the opposite directions for one local minima to be considered later as a
//		seed point candidate. Such point must pass through the verification 
//		process before it can be used as a seed in the tracking algorithm
//
/*
void FindAndAddMoreSeedPoints(CPoint &aPoint)
{
	// two vectors, one for each direction
	CVector *Vector1 = VectorsArray[aPoint.dir];
	CVector *Vector2 = VectorsArray[DirectionPlus(aPoint.dir, NumOfDirections/2)];
	
	// a point to hold the seed candidate
	CPoint result(0,0,0,0);
	// Search for the seed in one direction. The result is stored in (result)
	if(SearchForSeedAlongA_Line(aPoint, Vector1, result))
	{
		InitialPoints.Push(result);
		cout << "\n Initial Point Pushed: ";
		BranchingPointsImage->MarkPoint(&result, 254);
		//result.Print();
		//WorkingImage->MarkPoint(&result, 250);
	}
	// do the same for the other direcion
	if(SearchForSeedAlongA_Line(aPoint, Vector2, result))
	{
		InitialPoints.Push(result);
		cout << "\n Initial Point Pushed: ";
		BranchingPointsImage->MarkPoint(&result, 254);
		//result.Print();
		//WorkingImage->MarkPoint(&result, 250);
	}
}
*/

int compare5(const void* arg1, const void* arg2)
{
	int result = 0;
	int value1 = (*(CPoint**) arg1)->m_iValue;
	int value2 = (*(CPoint**) arg2)->m_iValue; 
	if (value1 < value2)
		result = 1;
	else if (value1 > value2)
		result = -1;
	return result;
}
void SortSeedPoints()
{
	register int Index = 0;
	int giROWS = The3DImage->m_iRows;
	int giCOLS = The3DImage->m_iCols;

	// fill the sorted array of pointers
	for (register int i = 0; i < giROWS; i++)
	{
		for (register int j = 0; j < giCOLS; j++)
		{
			if (gapArrayOfSeedPoints[i][j])
			{
				gapSortedArrayOfSeedPoints[Index] = gapArrayOfSeedPoints[i][j];
				Index++;
			}
		}
	}

	giNumOfSeedPoints = Index;
	qsort((void *) gapSortedArrayOfSeedPoints,
		giNumOfSeedPoints,
		sizeof(CPoint *),
		compare5);
}

void VerifyAndSortAllSeedPoints()
{
	register int HDir, VDir;

	for (HDir = 0; HDir < NumOfDirections; HDir++)
	{
		for (VDir = 0; VDir < NumOfDirections; VDir++)
		{
			gHLeftTemplatesArray[HDir][VDir]->AssociateWithImage(*The3DImage);
			gHRightTemplatesArray[HDir][VDir]->AssociateWithImage(*The3DImage);
			gVLeftTemplatesArray[HDir][VDir]->AssociateWithImage(*The3DImage);
			gVRightTemplatesArray[HDir][VDir]->AssociateWithImage(*The3DImage);
		}
	}

	CImage tempImage(*ProjectionImage);

	register int i, j;
	register int iIndex = 0;
	register CPoint * pPoint = NULL;
	register double dArea = 0;
	list<deque<CPoint> >::iterator k;

	int giROWS = The3DImage->m_iRows;
	int giCOLS = The3DImage->m_iCols;
	for (i = 0; i < giROWS; i++)
	{
		for (j = 0; j < giCOLS; j++)
		{
			pPoint = gapArrayOfSeedPoints[i][j];
			if (pPoint)
			{
				UnverifiedSeeds.push_back(*pPoint);
				if (VerifySeedPoint3D2(pPoint))
				{
					tempImage.data[pPoint->m_iY][pPoint->m_iX] = 254;
					pPoint->m_iValue = pPoint->m_iValue * ((int)
						(The3DImage->data[pPoint->m_iZ][pPoint->m_iY][pPoint->m_iX]));

					gapSortedArrayOfSeedPoints[iIndex] = pPoint;
					iIndex++;
					dArea += pPoint->m_iValue;
					VerifiedSeedsCenter.push_back(*pPoint);
				}
				else
				{
					tempImage.data[pPoint->m_iY][pPoint->m_iX] = 255;
					delete gapArrayOfSeedPoints[i][j];
					gapArrayOfSeedPoints[i][j] = NULL;
				}
			}
		}
	}

	giNumOfSeedPoints = iIndex;
	qsort((void *) gapSortedArrayOfSeedPoints,
		giNumOfSeedPoints,
		sizeof(CPoint *),
		compare5);
	Seeds.sort();

	VerifiedSeedsCenter.clear();
	list<deque<CPoint> >::reverse_iterator tp;
	
	CPoint Seed;
	for (tp = Seeds.rbegin(); tp != Seeds.rend(); tp++) {
		Seed = (*tp)[0];
		VerifiedSeedsCenter.push_back(Seed);
	}

	/*
		// output the points
		ofstream outFile("Points.txt");
		for(i = 0; i < giNumOfSeedPoints; i++)
			outFile << gapSortedArrayOfSeedPoints[i]->m_iValue << endl;
		register double dArea2 = 0;
		for( i = 0; i < giNumOfSeedPoints; i++) {
			
			dArea2 += gapSortedArrayOfSeedPoints[i]->m_iValue;
			// include all seeds points comprising 90% of the area
			if((dArea2/dArea) > 0.95) {
				cout << "\nThrshold: " << gapSortedArrayOfSeedPoints[i]->m_iValue << endl;
				break;
			}
			tempImage.data[gapSortedArrayOfSeedPoints[i]->m_iY][gapSortedArrayOfSeedPoints[i]->m_iX] = 254;
		}
		for(j = i; j < giNumOfSeedPoints; j++) {
			tempImage.data[gapSortedArrayOfSeedPoints[j]->m_iY][gapSortedArrayOfSeedPoints[j]->m_iX] = 255;
			giNumOfSeedPoints--;
			delete gapSortedArrayOfSeedPoints[j];
			gapSortedArrayOfSeedPoints[j] = NULL;
		}
	*/	
	//	cout << "\n\n A total of : " << giNumOfSeedPoints << " were found.\n"
	//		  << endl;

	//	tempImage.Write("AllSeeds.pgm");
	//exit(0);
}


