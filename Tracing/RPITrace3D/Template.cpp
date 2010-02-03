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

////////////////////////////////////////////////////////////////////////////
// File: template.cpp
// 3-24-98
//
// definition of the template class
//
// NOTE:
//   In all vector and template rotations, we first rotate around the y-axis
// then around the z-axix.
//#pragma warning(disable:4786)

#include <iostream>
#include <cmath>
#include <string>
#include <list>
#include <queue>

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

using namespace std;

// the first dimension stands for 45, 135, 225, and 315
// the second dimension stands for minus2, minus1, plus2, and plus1
extern int FourtyFiveDegreeVectorsLeft[4][4];
extern int FourtyFiveDegreeVectorsRight[4][4];

//////////////////////////////////////////////////////////////////////
////////////////// Class CTemplate  //////////////////////////////////
//////////////////////////////////////////////////////////////////////

// CTOR for the base class
CTemplate::CTemplate(int row, int col, int slice, float Hdir, float Vdir,
	int len) : m_fHdir(Hdir), m_fVdir(Vdir), m_iLength(len), m_iRows(row),
	m_iCols(col), m_iSlices(slice), m_puchBase(0), m_piPlus2(0), m_piPlus1(0),
	m_piMinus1(0), m_piMinus2(0), m_iPlus2Resp(0), m_iPlus1Resp(0),
	m_iMinus1Resp(0), m_iMinus2Resp(0), m_pMyImage(0), m_pMy3DImage(0)
{
	// create space for my rows
	m_piPlus2 = new int[m_iLength];
	m_piPlus1 = new int[m_iLength];
	m_piMinus1 = new int[m_iLength];
	m_piMinus2 = new int[m_iLength];
}

// DTOR
CTemplate::~CTemplate()
{
	if (m_piPlus2)
	{
		delete [] m_piPlus2;
		delete [] m_piPlus1;
		delete [] m_piMinus1;
		delete [] m_piMinus2;
	}
}

void CTemplate::ConstructTemplate()
{
	// the indices of the vectors parallel to the template
	int hIndex = (int) (m_fHdir / DirectionStep);
	int vIndex = (int) (m_fVdir / DirectionStep);

	// the vector that will shift the template away from the 
	// dendrite center
	CVector* minusVector = this->GetShiftFromCenterVector();

	int minus2Base = minusVector->m_piData[1];
	int minus1Base = minusVector->m_piData[2];
	int index = 3;
	while (minus1Base == minus2Base)
	{
		minus1Base = minusVector->m_piData[index];
		index++;
	}

	// the vector that will bring the template to the dendrite centerline
	CVector* plusVector = this->GetShiftToCenterVector();

	int plus2Base = plusVector->m_piData[1];
	int plus1Base = plusVector->m_piData[2];
	index = 3;
	while (plus1Base == plus2Base)
	{
		plus1Base = plusVector->m_piData[index];
		index++;
	}

	// fill the template rows
	CVector* myVector = gVectorsArray[hIndex][vIndex];
	for (register int i = 0; i < m_iLength; i++)
	{
		m_piMinus2[i] = minus2Base + myVector->m_piData[i];
		m_piMinus1[i] = minus1Base + myVector->m_piData[i];
		m_piPlus2[i] = plus2Base + myVector->m_piData[i];
		m_piPlus1[i] = plus1Base + myVector->m_piData[i];
	}
}
///////////////////////////////////////////
// Method: Position
//
// Position the template at the given point
void CTemplate::Position(int x, int y, int z)
{
	if (z != -1)
		m_puchBase = &(m_pMy3DImage->data[z][y][x]);
	else
		m_puchBase = &(m_pMyImage->data[y][x]);
}

///////////////////////////////////////////
// Method: Position
//
// Position the template at the given point
void CTemplate::Position(CPoint& aPoint)
{
	if (m_pMy3DImage)
	{
		// check that the point is valid
		if (m_pMy3DImage->ValidPoint(aPoint))
			m_puchBase = & (m_pMy3DImage->data[aPoint.m_iZ][aPoint.m_iY][aPoint.m_iX]);
	}
	else
	{
		// check that the point is valid
		if (m_pMyImage->ValidPoint(aPoint))
			m_puchBase = &(m_pMyImage->data[aPoint.m_iY][aPoint.m_iX]);
	}
}

// draw the template in the given image
void CTemplate::DisplayTemplate(C3DImage& aImage)
{
	for (register int i = 0; i < m_iLength; i++)
	{
		*(m_puchBase + m_piPlus1[i]) = (unsigned char) 255;
		*(m_puchBase + m_piPlus2[i]) = (unsigned char) 254;
		*(m_puchBase + m_piMinus2[i]) = (unsigned char) 253;
		*(m_puchBase + m_piMinus1[i]) = (unsigned char) 252;
	}
}

////////////////////////////////////////////////////////////////
// Method: Apply
//
// Apply the template at the given Point and return the response
int CTemplate::Apply(CPoint* aPoint, int usedWidth)
{
	// make sure that we are not at the image boundary
//	if (abs(aPoint->m_iX - 0) <
//		giMARGIN ||
//		abs(aPoint->m_iX - giCOLS) <
//		giMARGIN ||
//		abs(aPoint->m_iX -
//			0) <
//		giMARGIN ||
//		abs(aPoint->m_iY -
//			giROWS) <
//		giMARGIN ||
//		abs(aPoint->m_iZ -
//			0) <
//		giMARGIN ||
//		abs(aPoint->m_iZ -
//			giSLICES) <
//		giMARGIN)
//		return 0;
	if(!The3DImage->WithinImageMargin(*aPoint,giMARGIN))
		return 0;


	// get a pointer to the orthogonal shift vector
	//CVector* orthogonalVector = GetShiftFromCenterVector();

	// position my base at the given point
	Position(aPoint->m_iX, aPoint->m_iY, aPoint->m_iZ);

	// initialize
	register int minus1resp = 0;
	register int minus2resp = 0;
	register int plus1resp = 0;
	register int plus2resp = 0;
	register int i , pixelValue1, pixelValue2, pixelValue3, pixelValue4;

	long diff;
	//long pixelSum = 0;
	for (i = 0; i < usedWidth; i++)
	{
		if ((diff = (m_puchBase - gpuchThe3DImageData + m_piMinus1[i])) < 0 ||
			diff >= glThe3DImageLength ||
			(diff = (m_puchBase - gpuchThe3DImageData + m_piMinus2[i])) < 0 ||
			diff >= glThe3DImageLength ||
			(diff = (m_puchBase - gpuchThe3DImageData + m_piPlus1[i])) < 0 ||
			diff >= glThe3DImageLength ||
			(diff = (m_puchBase - gpuchThe3DImageData + m_piPlus2[i])) < 0 ||
			diff >= glThe3DImageLength)
			return 0;

		pixelValue1 = *(m_puchBase + m_piMinus1[i]);
		pixelValue2 = *(m_puchBase + m_piMinus2[i]);
		pixelValue3 = *(m_puchBase + m_piPlus1[i]);
		pixelValue4 = *(m_puchBase + m_piPlus2[i]);

		if (pixelValue1 < HighestReservedColor ||
			pixelValue2 < HighestReservedColor ||
			pixelValue3 < HighestReservedColor ||
			pixelValue4 < HighestReservedColor)
			return 0;
		else
		{
			minus1resp += pixelValue1;
			minus2resp += pixelValue2;
			plus1resp += pixelValue3;
			plus2resp += pixelValue4;
		}
	}

	int maxResult = plus1resp +
		plus2resp +
		plus2resp -
		minus1resp -
		minus2resp -
		minus2resp;

	return maxResult;
}

//////////////////////////////////////////////////////////////////
// METHOD: Apply
// Apply the templage at the given location and return the response
// this function is the same like the above but it takes a pointer to
// the image data at the base of the template.
int CTemplate::Apply(unsigned char* dataPtr)
{
	return 0;
}

// Shift the template by the given amount in a direction perpendicular
// to its direction
int CTemplate::ShiftAndApply(int index)
{
	return 0;
}


////////////////////////////////////////////////////////////////////////////////
// METHOD: CalculateMaxResponse
//
// Calculate the template responses at all points along the
// perpendicular vector originating from the given point and return
// the maximum of such responses. This function is shared by derived classes.
// the only difference between the behavior of the left and right templates 
// is due to the fact that the perpendicular direction to the left template is
// its direction + 4 (mod NumOfDirections) while the perpendicular to the right
// one is given by its direction -4 (mod NumOfDirections). This is taken care of
// by the pure virual function "GetOrthogonalDir()"
//
int CTemplate::CalculateMaxResponse(CPoint* aPoint, CPoint* atPoint,
	CVessel* aVessel, CPoint* pPoints)
{
	// make sure that we are not at the image boundary
	if(!The3DImage->WithinImageMargin(*aPoint,giMARGIN))
		return 0;

	// get a pointer to the orthogonal shift vector
	CVector* orthogonalVector = GetShiftFromCenterVector();

	// position my base at the given point
	Position(aPoint->m_iX, aPoint->m_iY, aPoint->m_iZ);

	// initialize
	register int minus1resp = 0;
	register int minus2resp = 0;
	register int plus1resp = 0;
	register int plus2resp = 0;
	register int i , j, pixelValue1, pixelValue2, pixelValue3, pixelValue4;
	//register int oldCount = 0;
	//register int newCount = 0;
	long diff;
	long pixelSum = 0;
	int HDir = (int) (m_fHdir / DirectionStep);
	int VDir = (int) (m_fVdir / DirectionStep);

	int giUsedTemplateLength = gConfig.GetMinimumTemplateLength();
	double SmallestAcceptedResponse = gfContrast * 3.0 * giUsedTemplateLength;

	for (i = 0; i < giUsedTemplateLength; i++)
	{
		if ((diff = (m_puchBase - gpuchThe3DImageData + m_piMinus1[i])) < 0 ||
			diff >= glThe3DImageLength ||
			(diff = (m_puchBase - gpuchThe3DImageData + m_piMinus2[i])) < 0 ||
			diff >= glThe3DImageLength ||
			(diff = (m_puchBase - gpuchThe3DImageData + m_piPlus1[i])) < 0 ||
			diff >= glThe3DImageLength ||
			(diff = (m_puchBase - gpuchThe3DImageData + m_piPlus2[i])) < 0 ||
			diff >= glThe3DImageLength)
			break;

		pixelValue1 = *(m_puchBase + m_piMinus1[i]);
		pixelValue2 = *(m_puchBase + m_piMinus2[i]);
		pixelValue3 = *(m_puchBase + m_piPlus1[i]);
		pixelValue4 = *(m_puchBase + m_piPlus2[i]);

		if (pixelValue1 < HighestReservedColor ||
			pixelValue2 < HighestReservedColor ||
			pixelValue3 < HighestReservedColor ||
			pixelValue4 < HighestReservedColor)
		{
			break;
		}
		else
		{
			minus1resp += pixelValue1;
			minus2resp += pixelValue2;
			plus1resp += pixelValue3;
			plus2resp += pixelValue4;
		}
	}

	// store the forground and background pixel values
	m_iForegroundPixelValue = (plus1resp + plus2resp) / giUsedTemplateLength;
	m_iBackgroundPixelValue = (minus1resp + minus1resp) /
		giUsedTemplateLength;

	pixelSum = (plus1resp + plus2resp) / 2;

	int maxResult = plus1resp +
		plus2resp +
		plus2resp -
		minus1resp -
		minus2resp -
		minus2resp;

	if (pPoints)
	{
		pPoints[0].m_iValue = maxResult;
		pPoints[0].m_iX = aPoint->m_iX;
		pPoints[0].m_iY = aPoint->m_iY;
		pPoints[0].m_iZ = aPoint->m_iZ;
		pPoints[0].m_iHDir = HDir;
		pPoints[0].m_iVDir = VDir;
		pPoints[0].m_lUserFlag = plus1resp + plus2resp;
	}

	int index = 0;
	i = 0;

	int giShiftDistance = gConfig.GetMaximumShiftDistance();

	// shift and apply
	for (j = 1; j < giShiftDistance; j++)
	{
		unsigned char * ShiftedBase = m_puchBase +
			orthogonalVector->m_piData[j];
		//pShiftedLinkList = pBaseLinkList + orthogonalVector->data[j];

		plus1resp = 0;
		plus2resp = 0;
		minus1resp = 0;
		minus2resp = 0;

		for (i = 0; i < giUsedTemplateLength; i++)
		{
			// make sure that we are still within the image boundaries
			if ((diff = (ShiftedBase - gpuchThe3DImageData + m_piMinus1[i])) <
				0 ||
				diff >= glThe3DImageLength ||
				(diff = (ShiftedBase - gpuchThe3DImageData + m_piMinus2[i])) <
				0 ||
				diff >= glThe3DImageLength ||
				(diff = (ShiftedBase - gpuchThe3DImageData + m_piPlus1[i])) <
				0 ||
				diff >= glThe3DImageLength ||
				(diff = (ShiftedBase - gpuchThe3DImageData + m_piPlus2[i])) <
				0 ||
				diff >= glThe3DImageLength)
				break;

			pixelValue1 = *(ShiftedBase + m_piMinus1[i]);
			pixelValue2 = *(ShiftedBase + m_piMinus2[i]);
			pixelValue3 = *(ShiftedBase + m_piPlus1[i]);
			pixelValue4 = *(ShiftedBase + m_piPlus2[i]);

			if (pixelValue1 < HighestReservedColor ||
				pixelValue2 < HighestReservedColor ||
				pixelValue3 < HighestReservedColor ||
				pixelValue4 < HighestReservedColor)
			{
				break;
			}
			else
			{
				minus1resp += pixelValue1;
				minus2resp += pixelValue2;
				plus1resp += pixelValue3;
				plus2resp += pixelValue4;
			}
		}

		int result = plus1resp +
			plus2resp +
			plus2resp -
			minus1resp -
			minus2resp -
			minus2resp;

		if (pPoints)
		{
			pPoints[j].m_iValue = result;
			pPoints[j].m_iX = aPoint->m_iX +
				orthogonalVector->m_pIndices[j].m_iX;
			pPoints[j].m_iY = aPoint->m_iY +
				orthogonalVector->m_pIndices[j].m_iY;
			pPoints[j].m_iZ = aPoint->m_iZ +
				orthogonalVector->m_pIndices[j].m_iZ;
			pPoints[j].m_iHDir = HDir;
			pPoints[j].m_iVDir = VDir;

			pPoints[j].m_lUserFlag += (pPoints[j - 1].m_lUserFlag + plus2resp);
		}

		//////////////////////////////////////////////////////
		// calculate the width difference penalty
		// according to the average width over the past 3 tracking steps
		// 
		if (result > maxResult)
		{
			maxResult = result;
			index = j;

			// store the forground and background pixel values
			m_iForegroundPixelValue = (plus1resp + plus2resp) /
				giUsedTemplateLength;
			m_iBackgroundPixelValue = (minus1resp + minus2resp) /
				giUsedTemplateLength;
		}


		////////
		// At any point, if we obtain a negative response after obtaining a
		// response that is larger than the "SmallestResponseAccepted" then stop
		// shifting. This is because such response is most likely to have been 
		// obtained at the boundary of an adjacent vessel.
		//if(maxResult > SmallestResponseAccepted && result < -100)
		if ((result < 0.2 * maxResult) && maxResult > SmallestAcceptedResponse)
			break;

		if (plus1resp < 0.4 * pixelSum)
			break;

		if (j >= 3 && maxResult < 0)
			break;
	}

	if (pPoints)
	{
		// calculate the pixel density of segment within each of the boundary points
		for (i = 0; i < j; i++)
		{
			pPoints[i].m_lUserFlag = pPoints[i].m_lUserFlag /
				((i + 1) * giUsedTemplateLength);
		}
	}

	int giROWS = The3DImage->m_iRows;
	int giCOLS = The3DImage->m_iCols;
	// if the caller needs the location of the maximum response point, 
	// return it
	if (maxResult > 0 && atPoint)
	{
		atPoint->m_iHDir = (int) (m_fHdir / DirectionStep);
		atPoint->m_iVDir = (int) (m_fVdir / DirectionStep);

		int x = aPoint->m_iX + orthogonalVector->m_pIndices[index].m_iX;
		int y = aPoint->m_iY + orthogonalVector->m_pIndices[index].m_iY;
		int z = aPoint->m_iZ + orthogonalVector->m_pIndices[index].m_iZ;

		if (x >= 0 && x < giCOLS && y >= 0 && y < giROWS)
		{
			atPoint->m_iX = x;
			atPoint->m_iY = y;
			atPoint->m_iZ = z;
			atPoint->m_iValue = maxResult;
		}
		else
		{
			maxResult = 0;
			atPoint->m_iValue = 0;
		}
	}
	else
		atPoint->m_iValue = maxResult;

	m_iForegroundPixelValue /= 2;
	m_iBackgroundPixelValue /= 2;

	return maxResult;
}

/////////////////////////
/// yet another version of the calculate response function
int CTemplate::CalculateResponses(CPoint* aPoint, CVessel* aVessel,
	CPoint aPoints[][MaxTemplateLength], int aiResponses[][MaxTemplateLength])
{
	// make sure that we are not at the image boundary
	if(!The3DImage->WithinImageMargin(*aPoint,giMARGIN))
		return 0;

	// get a pointer to the orthogonal shift vector
	CVector* orthogonalVector = GetShiftFromCenterVector();

	// position my base at the given point
	Position(aPoint->m_iX, aPoint->m_iY, aPoint->m_iZ);

	// initialize

	int giFromLength = gConfig.GetMinimumTemplateLength();
	int giToLength = gConfig.GetMaximumTemplateLength();

	register int minus1resp = 0;
	register int minus2resp = 0;
	register int plus1resp = 0;
	register int plus2resp = 0;
	register int i , j, pixelValue1, pixelValue2, pixelValue3, pixelValue4;
	//register int oldCount = 0;
	//register int newCount = 0;
	long diff;
	long pixelSum = 0;
	int HDir = (int) (m_fHdir / DirectionStep);
	int VDir = (int) (m_fVdir / DirectionStep);
	int response = 0, maxResponse = 0;
	int maxResponseAtThisShift = 0;
	//int iGoodResponseThreshold = 3.0 * gfContrast* giFromLength;

	//int giUsedTemplateLength = gConfig.GetMinimumTemplateLength();
	//int SmallestAcceptedResponse = gfContrast*3.0 * giUsedTemplateLength;

	for (i = 0; i <= giToLength; i++)
	{
		if ((diff = (m_puchBase - gpuchThe3DImageData + m_piMinus1[i])) < 0 ||
			diff >= glThe3DImageLength ||
			(diff = (m_puchBase - gpuchThe3DImageData + m_piMinus2[i])) < 0 ||
			diff >= glThe3DImageLength ||
			(diff = (m_puchBase - gpuchThe3DImageData + m_piPlus1[i])) < 0 ||
			diff >= glThe3DImageLength ||
			(diff = (m_puchBase - gpuchThe3DImageData + m_piPlus2[i])) < 0 ||
			diff >= glThe3DImageLength)
			break;

		pixelValue1 = *(m_puchBase + m_piMinus1[i]);
		pixelValue2 = *(m_puchBase + m_piMinus2[i]);
		pixelValue3 = *(m_puchBase + m_piPlus1[i]);
		pixelValue4 = *(m_puchBase + m_piPlus2[i]);

		if (pixelValue1 < HighestReservedColor ||
			pixelValue2 < HighestReservedColor ||
			pixelValue3 < HighestReservedColor ||
			pixelValue4 < HighestReservedColor)
		{
			// if you hit the soma and the template is short set a flag
			// and return. The flag will be used to stop tracing
			if ((pixelValue1 == SomaColor ||
				pixelValue2 == SomaColor ||
				pixelValue3 == SomaColor ||
				pixelValue4 == SomaColor) &&
				i <
				giFromLength)
			{
				giHitSomaFlag = 1;
				return 0;
			}
			else
				break;
		}
		else
		{
			minus1resp += pixelValue1;
			minus2resp += pixelValue2;
			plus1resp += pixelValue3;
			plus2resp += pixelValue4;
		}	

		if (i >= giFromLength)
		{
			maxResponse = plus1resp +
				plus2resp +
				plus2resp -
				minus1resp -
				minus2resp -
				minus2resp;

			aiResponses[0][i] = maxResponse;

			aPoints[0][i].m_iValue = maxResponse;
			aPoints[0][i].m_iX = aPoint->m_iX;
			aPoints[0][i].m_iY = aPoint->m_iY;
			aPoints[0][i].m_iZ = aPoint->m_iZ;
			aPoints[0][i].m_iHDir = HDir;
			aPoints[0][i].m_iVDir = VDir;
			aPoints[0][i].m_lUserFlag = plus1resp + plus2resp;
		}
	}
	pixelSum = (plus1resp + plus2resp) / 2;

	int giShiftDistance = gConfig.GetMaximumShiftDistance();
	//int index = 0;
	i = 0;
	// shift and apply
	for (j = 1; j < giShiftDistance; j++)
	{
		unsigned char * ShiftedBase = m_puchBase +
			orthogonalVector->m_piData[j];
		//pShiftedLinkList = pBaseLinkList + orthogonalVector->data[j];		

		plus1resp = 0;
		plus2resp = 0;
		minus1resp = 0;
		minus2resp = 0; 
		for (i = 0; i <= giToLength; i++)
		{
			// make sure that we are still within the image boundaries

			if ((diff = (ShiftedBase - gpuchThe3DImageData + m_piMinus1[i])) <
				0 ||
				diff >= glThe3DImageLength ||
				(diff = (ShiftedBase - gpuchThe3DImageData + m_piMinus2[i])) <
				0 ||
				diff >= glThe3DImageLength ||
				(diff = (ShiftedBase - gpuchThe3DImageData + m_piPlus1[i])) <
				0 ||
				diff >= glThe3DImageLength ||
				(diff = (ShiftedBase - gpuchThe3DImageData + m_piPlus2[i])) <
				0 ||
				diff >= glThe3DImageLength)
				break;

			pixelValue1 = *(ShiftedBase + m_piMinus1[i]);
			pixelValue2 = *(ShiftedBase + m_piPlus2[i]);			
			pixelValue3 = *(ShiftedBase + m_piMinus2[i]);
			pixelValue4 = *(ShiftedBase + m_piPlus1[i]);

			if (pixelValue1 < HighestReservedColor ||
				pixelValue2 < HighestReservedColor ||
				pixelValue3 < HighestReservedColor ||
				pixelValue4 < HighestReservedColor)
			{
				// if you hit the soma and the template is short set a flag
				// and return. The flag will be used to stop tracing
				if ((pixelValue1 == SomaColor || pixelValue2 == SomaColor) &&
					(i < giFromLength && j < 4))
				{
					giHitSomaFlag = 1;
					return 0;
				}
				else
					break;
			}
			else
			{
				minus1resp += pixelValue1;
				plus2resp += pixelValue2;
				minus2resp += pixelValue3;
				plus1resp += pixelValue4;
			}

			maxResponseAtThisShift = 0;
			if (i >= giFromLength)
			{
				response = plus1resp +
					plus2resp +
					plus2resp -
					minus1resp -
					minus2resp -
					minus2resp;

				aiResponses[j][i] = response;

				aPoints[j][i].m_iValue = response;
				aPoints[j][i].m_iX = aPoint->m_iX +
					orthogonalVector->m_pIndices[j].m_iX;
				aPoints[j][i].m_iY = aPoint->m_iY +
					orthogonalVector->m_pIndices[j].m_iY;
				aPoints[j][i].m_iZ = aPoint->m_iZ +
					orthogonalVector->m_pIndices[j].m_iZ;
				aPoints[j][i].m_iHDir = HDir;
				aPoints[j][i].m_iVDir = VDir;

				aPoints[j][i].m_lUserFlag += (aPoints[j - 1][i].m_lUserFlag +
					plus2resp);

				if (response > maxResponseAtThisShift)
					maxResponseAtThisShift = response;
			}
		}

		if (maxResponseAtThisShift > maxResponse)
			maxResponseAtThisShift = maxResponse;
	}

	return maxResponse;
}

/////////////////////////
/// yet another version of the calculate response function
int CTemplate::CalculateKernelResponses(const CPoint& aPoint, int prev_shift,
	CPoint aPoints[][MaxTemplateLength],
	short aiResponses[][MaxTemplateLength])
{
	// make sure that we are not at the image boundary
	if(!The3DImage->WithinImageMargin(aPoint, giMARGIN))
		return 0;

	// get a pointer to the orthogonal shift vector
	CVector* orthogonalVector = GetShiftFromCenterVector();

	// position my base at the given point
	Position(aPoint.m_iX, aPoint.m_iY, aPoint.m_iZ);

	// initialize
	register unsigned char minus1resp = 0;
	register unsigned char minus2resp = 0;
	register unsigned char plus1resp = 0;
	register unsigned char plus2resp = 0;
	register unsigned char i , j , t;
	register long diff_minus1, diff_minus2, diff_plus1, diff_plus2;
	register unsigned char HDir = (int) (m_fHdir / DirectionStep);
	register unsigned char VDir = (int) (m_fVdir / DirectionStep);
	
	int giFromLength = gConfig.GetMinimumTemplateLength();
	int giToLength = gConfig.GetMaximumTemplateLength();

	register unsigned char to_length = giToLength;
	register unsigned char from_length = giFromLength;
	register long image_length = glThe3DImageLength;

	i = 0;
	// shift and apply
	int shift_freedom = gConfig.GetRelativeShiftDistance();

	char min_shift, max_shift;
	int radius = Round((float) prev_shift / 2.0);

	int giShiftDistance = gConfig.GetMaximumShiftDistance();
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

	for (j = min_shift; j <= max_shift; j++)
	{
		unsigned char * ShiftedBase = m_puchBase +
			orthogonalVector->m_piData[j];

		diff_minus1 = ShiftedBase - gpuchThe3DImageData + m_piMinus1[i];
		diff_minus2 = ShiftedBase - gpuchThe3DImageData + m_piMinus2[i];
		diff_plus1 = ShiftedBase - gpuchThe3DImageData + m_piPlus1[i];
		diff_plus2 = ShiftedBase - gpuchThe3DImageData + m_piPlus2[i];

		if (diff_minus1 < 0 ||
			diff_minus2 < 0 ||
			diff_plus1 < 0 ||
			diff_plus2 < 0 ||
			diff_minus1 >= image_length ||
			diff_minus2 >= image_length ||
			diff_plus1 >= image_length ||
			diff_plus2 >= image_length)
			continue;

		for (i = 0; i <= to_length; i++)
		{
			// make sure that we are still within the image boundaries

			diff_minus1 = ShiftedBase - gpuchThe3DImageData + m_piMinus1[i];
			diff_minus2 = ShiftedBase - gpuchThe3DImageData + m_piMinus2[i];
			diff_plus1 = ShiftedBase - gpuchThe3DImageData + m_piPlus1[i];
			diff_plus2 = ShiftedBase - gpuchThe3DImageData + m_piPlus2[i];

			if (diff_minus1 < 0 ||
				diff_minus2 < 0 ||
				diff_plus1 < 0 ||
				diff_plus2 < 0 ||
				diff_minus1 >= image_length ||
				diff_minus2 >= image_length ||
				diff_plus1 >= image_length ||
				diff_plus2 >= image_length)
			{
				if (i <= from_length)
				{
					for (t = 0; t <= i; t++)
					{
						aiResponses[j][t] = 0;
					}
				}
				break;
			}	

			minus1resp = *(ShiftedBase + m_piMinus1[i]);				
			minus2resp = *(ShiftedBase + m_piMinus2[i]);
			plus1resp = *(ShiftedBase + m_piPlus1[i]);
			plus2resp = *(ShiftedBase + m_piPlus2[i]);

			if (minus1resp < HighestReservedColor ||
				minus2resp < HighestReservedColor ||
				plus1resp < HighestReservedColor ||
				plus2resp < HighestReservedColor)
			{
				// if you hit the soma and the template is short set a flag
				// and return. The flag will be used to stop tracing
				if ((minus1resp == SomaColor ||
					minus2resp == SomaColor ||
					plus1resp == SomaColor ||
					plus2resp == SomaColor) &&
					i <
					from_length &&
					j <
					4)
				{
					giHitSomaFlag = 1;
					if (i <= from_length)
					{
						for (t = 0; t <= i; t++)
						{
							aiResponses[j][t] = 0;
						}
					}
					return 0;
				}
				else
				{
					if (i <= from_length)
					{
						for (t = 0; t <= i; t++)
						{
							aiResponses[j][t] = 0;
						}
					}
					break;
				}
			}

			aiResponses[j][i] = plus1resp +
				plus2resp +
				plus2resp -
				minus1resp -
				minus2resp -
				minus2resp;

			aPoints[j][i].m_iValue = aiResponses[j][i];
			aPoints[j][i].m_iX = aPoint.m_iX +
				orthogonalVector->m_pIndices[j].m_iX;
			aPoints[j][i].m_iY = aPoint.m_iY +
				orthogonalVector->m_pIndices[j].m_iY;
			aPoints[j][i].m_iZ = aPoint.m_iZ +
				orthogonalVector->m_pIndices[j].m_iZ;
			aPoints[j][i].m_iHDir = HDir;
			aPoints[j][i].m_iVDir = VDir;

			//	aPoints[j][i].m_lUserFlag += (aPoints[j-1][i].m_lUserFlag + acc_plus2resp);
		}
	}

	return 0;
}

///////////////////////////////////////////////////////////////////////////
// MEHTOD: ComputeShifts
// 
// A utility method used by the constructors of derived classes. The constructor
// passes the array (i.e. minus2, minus2, plus2, or plus2) which is by definition
// is parallel to the x-axis and passed through the point (0,y). 
// 
// Note: The CTORs calls this method for angles other than 45, 135, 225, or 315.
void CTemplate::ComputeShifts(int* row, int y)
{
	/* 2-21-99
	// my directions in radians, notice the minus sign because the y-axis
	// (i.e. the row axis) in an image increasing downward
	double theta = (-2.0 * 3.1415927 * direction)/360.0;
	  // working variables
	  register int x, index;
	  double fcol = 0, frow = 0;
	  int icol, irow;
	  
		double cosTheta = cos(theta);
		double sinTheta = sin(theta);
		
		  // for each cell in the row, calculate its value
		  //for(x = 0; x < length; x++)
		  for(x = 0, index = 0; index < length; x++)
		  {
		  ///////////////////////////////////////////////////////////////////
		  // calculate the new coordinates of the row after being rotated by 
		  // the template's angle. 
		  fcol = x * cosTheta - y * sinTheta;
		  frow = x * sinTheta + y * cosTheta;
		  
			// round to the nearest integer
			if(fcol   > 0)   icol   = (int)(fcol   + 0.5); 
			else				icol   = (int)(fcol   - 0.5);
			if(frow   > 0)   irow   = (int)(frow   + 0.5); 
			else				 irow   = (int)(frow   - 0.5);
			
			  // calculate the shift from the base of the template to each of the
			  // minus2 entries
			  //row[x] = islice*sliceSize + irow*cols + icol; 
			  int temp = irow*cols + icol;
			  
				if(index-1 >= 0)
				{
				if(temp != row[index-1])
				{
				row[index] = temp;
				index++;
				}
				}
				else
				{
				row[index] = temp;
				index++;
				}
				}
		*/
}	

//////////////////////////////////////////////////////////////////////
////////////////// Class CHLeftTemplate  /////////////////////////////
//////////////////////////////////////////////////////////////////////

// CTOR. Needs the template's length, directions, and the dimensions of the 
// image it is to be applied on.
//
// the rows of the template can be thought of as vectors (see CVector.h) which are
// in the same direction as the template itself
// the bases of these vectors are alligned along another vector lying in the 
// xy-plane and perpendicular to the template's horizontal direction
CHLeftTemplate::CHLeftTemplate(int row, int col, int slice, float Hdir,
	float Vdir, int len) : CTemplate(row, col, slice, Hdir, Vdir, len)
{
	ConstructTemplate();
}

// DTOR
CHLeftTemplate::~CHLeftTemplate()
{
}

// return the direction orthogonal to my direction (i.e. direction of shift)
int CHLeftTemplate::GetInPlaneShiftDir()
{
	int x = (int) (m_fHdir* NumOfDirections / 360 + 0.5);
	return DirectionMinus(x, NumOfDirections / 4);
}

// return the direction orthogonal to my direction (i.e. direction of shift)
int CHLeftTemplate::GetOrthogonalShiftDir()
{
	int x = (int) (m_fVdir* NumOfDirections / 360 + 0.5);
	return DirectionPlus(x, NumOfDirections / 4);
}

/////////////////////////////////
// Method: GetShiftFromCenterVector
//
// return a pointer to a vector pointing in the direction that
// will shift this template a way from a dendrite's centerline
CVector* CHLeftTemplate::GetShiftFromCenterVector()
{
	//int x = (int) (m_fVdir* NumOfDirections / 360 + 0.5);
	int index = GetInPlaneShiftDir();
	return gVectorsArray[index][0];
}

/////////////////////////////////
// Method: GetShiftToCenterVector
//
// return a pointer to a vector pointing in the direction that
// will shift this template towards a dendrite's centerline
CVector* CHLeftTemplate::GetShiftToCenterVector()
{
	//int x = (int) (m_fVdir* NumOfDirections / 360 + 0.5);
	int index = GetInPlaneShiftDir();
	index += (NumOfDirections / 2); // make index point to the opposite direction
	index = (index % NumOfDirections);

	return gVectorsArray[index][0];
}

//////////////////////////////////////////////////////////////////////
////////////////// Class CHRightTemplate /////////////////////////////
//////////////////////////////////////////////////////////////////////

// CTOR. needs the template's length, direction, and the image it is
// created on.
CHRightTemplate::CHRightTemplate(int row, int col, int slice, float Hdir,
	float Vdir, int length) : CTemplate(row, col, slice, Hdir, Vdir, length)
{
	ConstructTemplate();
}


// DTOR
CHRightTemplate::~CHRightTemplate()
{
}

// return the direction orthogonal to my direction (i.e. direction of shift)
int CHRightTemplate::GetInPlaneShiftDir()
{
	int x = (int) (m_fHdir* NumOfDirections / 360 + 0.5);
	return DirectionPlus(x, NumOfDirections / 4);
}

// return the direction orthogonal to my direction (i.e. direction of shift)
int CHRightTemplate::GetOrthogonalShiftDir()
{
	int x = (int) (m_fVdir* NumOfDirections / 360 + 0.5);
	return DirectionPlus(x, NumOfDirections / 4);
}

/////////////////////////////////
// Method: GetShiftFromCenterVector
//
// return a pointer to a vector pointing in the direction that
// will shift this template a way from a dendrite's centerline
CVector* CHRightTemplate::GetShiftFromCenterVector()
{
	//int x = (int) (m_fVdir* NumOfDirections / 360 + 0.5);
	int index = GetInPlaneShiftDir();
	return gVectorsArray[index][0];
}

/////////////////////////////////
// Method: GetShiftToCenterVector
//
// return a pointer to a vector pointing in the direction that
// will shift this template towards a dendrite's centerline
CVector* CHRightTemplate::GetShiftToCenterVector()
{
	//int x = (int) (m_fVdir* NumOfDirections / 360 + 0.5);
	int index = GetInPlaneShiftDir();
	index += (NumOfDirections / 2); // make index point to the opposite direction
	index = (index % NumOfDirections);

	return gVectorsArray[index][0];
}

//////////////////////////////////////////////////////////////////////
////////////////// Class CVLeftTemplate  /////////////////////////////
//////////////////////////////////////////////////////////////////////

// CTOR. needs the template's length, direction, and the dimensions of the 
// image it is to be applied on.
CVLeftTemplate::CVLeftTemplate(int row, int col, int slice, float Hdir,
	float Vdir, int len) : CTemplate(row, col, slice, Hdir, Vdir, len)
{
	ConstructTemplate();
}

// DTOR
CVLeftTemplate::~CVLeftTemplate()
{
}

// return the direction orthogonal to my direction (i.e. direction of shift)
int CVLeftTemplate::GetInPlaneShiftDir()
{
	int x = (int) (m_fVdir* NumOfDirections / 360 + 0.5);
	return DirectionMinus(x, NumOfDirections / 4);
}

// return the direction orthogonal to my direction (i.e. direction of shift)
int CVLeftTemplate::GetOrthogonalShiftDir()
{
	int x = (int) (m_fHdir* NumOfDirections / 360 + 0.5);
	return DirectionMinus(x, NumOfDirections / 4);
}

/////////////////////////////////
// Method: GetShiftFromCenterVector
//
// return a pointer to a vector pointing in the direction that
// will shift this template a way from a dendrite's centerline
CVector* CVLeftTemplate::GetShiftFromCenterVector()
{
	int x = (int) (m_fHdir* NumOfDirections / 360 + 0.5);
	int index = GetInPlaneShiftDir();
	return gVectorsArray[x][index];
}

/////////////////////////////////
// Method: GetShiftToCenterVector
//
// return a pointer to a vector pointing in the direction that
// will shift this template towards a dendrite's centerline
CVector* CVLeftTemplate::GetShiftToCenterVector()
{
	int x = (int) (m_fHdir* NumOfDirections / 360 + 0.5);
	int index = GetInPlaneShiftDir();
	index += (NumOfDirections / 2); // make index point to the opposite direction
	index = (index % NumOfDirections);

	return gVectorsArray[x][index];
}
//////////////////////////////////////////////////////////////////////
////////////////// Class CVRightTemplate /////////////////////////////
//////////////////////////////////////////////////////////////////////

// CTOR. needs the template's length, direction, and the image it is
// created on.
CVRightTemplate::CVRightTemplate(int row, int col, int slice, float Hdir,
	float Vdir, int length) : CTemplate(row, col, slice, Hdir, Vdir, length)
{
	ConstructTemplate();
}

// DTOR
CVRightTemplate::~CVRightTemplate()
{
}

// return the direction orthogonal to my direction (i.e. direction of shift)
int CVRightTemplate::GetInPlaneShiftDir()
{
	int x = (int) (m_fVdir* NumOfDirections / 360 + 0.5);
	return DirectionPlus(x, NumOfDirections / 4);
}

// return the direction orthogonal to my direction (i.e. direction of shift)
int CVRightTemplate::GetOrthogonalShiftDir()
{
	int x = (int) (m_fHdir* NumOfDirections / 360 + 0.5);
	return DirectionPlus(x, NumOfDirections / 4);
}

/////////////////////////////////
// Method: GetShiftFromCenterVector
//
// return a pointer to a vector pointing in the direction that
// will shift this template a way from a dendrite's centerline
CVector* CVRightTemplate::GetShiftFromCenterVector()
{
	int x = (int) (m_fHdir* NumOfDirections / 360 + 0.5);
	int index = GetInPlaneShiftDir();
	return gVectorsArray[x][index];
}

/////////////////////////////////
// Method: GetShiftToCenterVector
//
// return a pointer to a vector pointing in the direction that
// will shift this template towards a dendrite's centerline
CVector* CVRightTemplate::GetShiftToCenterVector()
{
	int x = (int) (m_fHdir* NumOfDirections / 360 + 0.5);
	int index = GetInPlaneShiftDir();
	index += (NumOfDirections / 2); // make index point to the opposite direction
	index = (index % NumOfDirections);

	return gVectorsArray[x][index];
}

