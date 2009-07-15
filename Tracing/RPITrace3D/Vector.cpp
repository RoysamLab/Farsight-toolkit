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

///////////////////////////////////////////////////////////////////////////
// File: Vector.cpp
// By: Khalid Al-Kofahi
// Created: 2-21-99
//
// Description:
//		contains the definition of the class CVector, representing a 3D vector
// 

//#pragma warning(disable:4786)  // disable STL-related warnings
//#pragma warning(disable:4710)  // compiler didn't inline function that got selected for inline expansion
//#pragma warning(disable:4514)  // The optimizer removed an inline function that is not called
//#pragma warning(disable:4702)  // unreachable STLport code

#include <iostream>
#include <cmath>
#include <stdlib.h> /* exit() */

#include "CONSTANTS.h"
#include "Mytypes.h"
#include "Config.h"
#include "Cpoints.h"
#include "Dllist.h"
#include "Cimage.h"
#include "C3dimage.h"
#include "Cvector.h"

using namespace std;

extern void MatrixMultiply(double Left[][4], double Right[][4],
	double Result[][4], int x, int y);

#define LargeNegativeValue -99999
#define Epselon 0.3

typedef struct _junkPoint
{
	float X;
	float Y;
	float Z;
} JunkPoint;

int LessThanEpselon(double value1, double value2)
{
	double error = (value1 < value2) ? (value2 - value1) : (value1 - value2);
	return (error < Epselon);
}

extern void Construct3DLine(CPoint& from, CPoint& to, CDLList<CPoint>& result);
extern void MatrixMultiply(CPoint& aPoint, double Right[][4], CPoint& result);

/////////////////////////////////////
// Method: CTOR
//
// To construct a vector in the directions theta, and phi, we do the following
// 1. Rotate the point P1(m_iLength, 0, 0) by Theta around the z-axis,
//    resulting in P1'(x1, y1, 0). R_Theta
// 2. Rotate the point P2(0, m_iLength, 0) by Theta around the z-axis,
//		resulting in P2'(x2, y2, 0). R_Theta
// 3. Rotate the point P1' by Phi around the axis L from the origin to P2'
//    To do this:
//		- Aligh the axis L with the z-axis as follows
//			3.1 - Rotate P2' & P1' by 90 around the X-axis, 
//					resulting in P2'', and P1'' . R_X
//			3.2 - Rotate P2'' & P1'' by 90-Theta around the Y-axis,
//					resulting in P1'''. R_Y
//		- Rotate the point P1''' by Phi around the z-axis. R_Phi
//	4. Perform the inverse rotations to those in 3.1, and 3.2, R_YInv, R_XInv
//
// where	
//				 CosT	SinT	0		0			
//				-SinT	CosT	0		0			
// R_Theta = 0		0     1	   0	  
//			  0    0	 0		1			
//
//				 CosP SinP	0		0
//				-SinP CosP	0		0
// R_Phi   = 0		0		1		0
//				 0		0		0		1
//
//				 1		0		0		0					  1	0	0	0
//				 0		0		1		0					  0	0 -1	0
//	R_X 	= 0   -1		0		0		, R_XInv = 0	1  0  0
//				 0    0		0		1					  0   0  0  1
//
//				 CosT 0		SinT	0					  CosT	0	-SinT	0
//				 0		1		0		0					  0 	 1   0		0
// R_Y     =-SinT 0		CosT	0		, R_YInv = SinT	0   CosT 0
//				 0		0		0		1					  0 	 0   0    1
//
// and the final rotation is given by
// [R_Theta][R_X][R_Y][R_Phi][R_YInv][R-XInv]
//
//
//						 	 CosT*CosP	 SinT	  CosT*SinP
//	R_Theta * R_Phi = -SinT*CosP   CosT  -SinT*SinP	
//							-SinP			 0		  CosP
//
//
CVector::CVector(int , int rows, int cols, float Hdir, float Vdir,
	int len) : m_fHDir(Hdir), m_fVDir(Vdir), m_iLength(len), m_iRows(rows),
	m_iCols(cols), m_piData(0), m_iBase(0), m_pIndices(0)
{
	// my directions in radians, notice the minus sign because the x-axis
	// (i.e. the row axis) in an image increasing downward
	double theta = (2.0 * 3.1415927 * m_fHDir) / 360.0;
	double Phi = (-2.0 * 3.1415927 * m_fVDir) / 360.0;
	double CosT = cos(theta);
	double SinT = sin(theta);
	double CosP = cos(Phi);
	double SinP = sin(Phi);
	int sliceSize = m_iCols* m_iRows;

	// the end point of the vector before rotation (lies on the x-axis)
	CPoint oldPoint(m_iLength, 0, 0);
	// the end point after rotation
	CPoint newPoint;

	double R_Theta[4][4] =
	{
		{CosT, SinT, 0, 0}, { - SinT, CosT, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}
	};
	double R_Phi[4][4] =
	{
		{CosP, SinP, 0, 0}, { - SinP, CosP, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}
	};
	double R_X[4][4] =
	{
		{1, 0, 0, 0}, {0, 0, 1, 0}, {0, -1, 0, 0}, {0, 0, 0, 1}
	};
	double R_XInv[4][4] =
	{
		{1, 0, 0, 0}, {0, 0, -1, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}
	};
	double R_Y[4][4] =
	{
		{CosT, 0, -SinT, 0}, {0, 1, 0, 0}, {SinT, 0, CosT, 0}, {0, 0, 0, 1}
	};
	double R_YInv[4][4] =
	{
		{CosT, 0,SinT, 0}, {0, 1, 0, 0}, { - SinT, 0, CosT, 0}, {0, 0, 0, 1}
	};

	double Result1[4][4];
	double Result2[4][4];

	CPoint ResultPoint(m_iLength, 0, 0);
	CPoint point1(m_iLength, 0, 0);

	MatrixMultiply(point1, R_Theta, ResultPoint);
	MatrixMultiply(ResultPoint, R_X, point1);
	MatrixMultiply(point1, R_Y, ResultPoint);
	MatrixMultiply(ResultPoint, R_Phi, point1);
	MatrixMultiply(point1, R_YInv, ResultPoint);
	MatrixMultiply(ResultPoint, R_XInv, point1);

	// multiply
	MatrixMultiply(R_Theta, R_X, Result1, 4, 4);
	MatrixMultiply(Result1, R_Y, Result2, 4, 4);
	MatrixMultiply(Result2, R_Phi, Result1, 4, 4);
	MatrixMultiply(Result1, R_YInv, Result2, 4, 4);
	MatrixMultiply(Result2, R_XInv, Result1, 4, 4);

	double x = m_iLength;
	CPoint EndPoint;
	// Ideally we should use the following formula
	//EndPoint.m_iX = x*Result1[0][0] + y*Result1[1][0] + z*Result1[2][0];
	//EndPoint.m_iY = x*Result1[0][1] + y*Result1[1][1] + z*Result1[2][1];
	//EndPoint.m_iZ = x*Result1[0][2] + y*Result1[1][2] + z*Result1[2][2];

	double Xprime = x* Result1[0][0];
	double Yprime = x* Result1[0][1];
	double Zprime = x* Result1[0][2];
	if (Xprime >= 0.0)
		Xprime = Xprime + 0.5;
	else
		Xprime = Xprime - 0.5;
	if (Yprime >= 0.0)
		Yprime = Yprime + 0.5;
	else
		Yprime = Yprime - 0.5;
	if (Zprime >= 0.0)
		Zprime = Zprime + 0.5;
	else
		Zprime = Zprime - 0.5;

	EndPoint.m_iX = (int) Xprime;
	EndPoint.m_iY = (int) Yprime;
	EndPoint.m_iZ = (int) Zprime;

	CPoint Origin(0, 0, 0);
	CDLList<CPoint> thePoints;

	// now, draw the vector (i.e. origin -to- (newX, newY, newZ)) and store
	// the points in the given list
	Construct3DLine(Origin, EndPoint, thePoints);

	// the actual length of the vector is given by the number of points
	// in the constructed list
	m_iLength = thePoints.length;	
	m_piData = new int[m_iLength];
	m_pIndices = new CPoint[m_iLength];
	if (!m_piData || ! m_pIndices)
	{
		cout << "\n CVector::Out of Memory" << endl;
		exit(0);
	}
	memset(m_piData, 0, m_iLength * sizeof(int));
	memset(m_pIndices, 0, m_iLength * sizeof(CPoint));

	// store the point indices and the corresponding shifts from the base
	int i = 0;
	CLNode<CPoint>* temp = thePoints.head;
	while (temp)
	{
		m_pIndices[i] = *temp->data;
		m_piData[i] = m_pIndices[i].m_iZ * sliceSize +
			m_pIndices[i].m_iY * m_iCols +
			m_pIndices[i].m_iX;

		i++;
		temp = temp->after;
	}
}


// DTOR
CVector::~CVector()
{
	if (m_piData)
		delete [] m_piData;
	if (m_pIndices)
		delete [] m_pIndices;
}

// draw the template in the given image
void CVector::DisplayVector(C3DImage& aImage, int color)
{
	// a pointer to the base pixel in the image
	unsigned char * basePtr = &(aImage.data[0][0][0]) + m_iBase;
	for (register int i = 0; i < m_iLength; i++)
	{
		*(basePtr + m_piData[i]) = (unsigned char) color;
	}
}

//////////////////////////////////////////////
// Method: Position
// 
// position the base of the vector at the given point
void CVector::Position(int x, int y, int z)
{
	m_pIndices[0].m_iX = x;
	m_pIndices[0].m_iY = y; 
	m_pIndices[0].m_iZ = z;
	m_iBase = z * (m_iRows * m_iCols) + y * m_iRows + x;
}
