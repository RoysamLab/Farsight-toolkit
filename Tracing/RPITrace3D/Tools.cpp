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

//////////////////////////////////////////////////////////////////////////////
// FILE: tools.cpp
// 
// This file contains some usefull tools available for use of other functions
//
//#pragma warning(disable:4786)  // disable STL-related warnings
//#pragma warning(disable:4710)  // compiler didn't inline function that got selected for inline expansion
//#pragma warning(disable:4514)  // The optimizer removed an inline function that is not called
//#pragma warning(disable:4702)  // unreachable STLport code

#include <iostream>
#include <fstream>
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

//#define StrucElemSize 5
int StrucElemSize = 0;

int Round(double value)
{
	return static_cast<int>(((value > 0.0) ? (value + 0.5) : (value - 0.5)));
}

int DirectionMinus(int d, int v)
{
	register int x = d - v;
	return ((x < 0) ? (x + NumOfDirections) : x);
}

int DirectionPlus(int d, int v)
{
	register int x = d + v;
	return ((x >= NumOfDirections) ? (x - NumOfDirections) : x);
}

int DirDistance(int d1, int d2)
{
	int x = (d1 > d2) ? d1 - d2 : d2 - d1;

	return ((x > (NumOfDirections / 2)) ? NumOfDirections - x : x);
}

// calculate the search directions around my current direction
void PrepareDirectionsArray(int* directions, int size, int dir)
{
	directions[0] = dir;
	int index = 1;
	for (register int i = 1; i < size; i += 2)
	{
		directions[i] = DirectionMinus(dir, index);
		directions[i + 1] = DirectionPlus(dir, index);
		index++;
	}
}

// draw a line from point-a to point-b in the image
void draw_line(CPoint& a, CPoint& b, CImage& image, unsigned char color)
{
	int xstart = a.m_iX, ystart = a.m_iY;
	int xend = b.m_iX, yend = b.m_iY;
	int i,x,y,dx,dy,e;
	int signx, signy, change, temp;

	if (xstart == xend && ystart == yend)
	{
		image.data[ystart][xstart] = color;
		return;
	}
	if (xend < xstart)
	{
		signx = -1;
		dx = xstart - xend;
	}
	else
	{
		signx = 1;
		dx = xend - xstart;
	}
	if (yend < ystart)
	{
		signy = -1;
		dy = ystart - yend;
	}
	else
	{
		signy = 1;
		dy = yend - ystart;
	}

	if (dy <= dx)
		change = 0;
	else
	{
		change = 1;
		temp = dx; dx = dy; dy = temp;
	}
	e = (dy << 1) - dx;
	x = xstart; y = ystart;
	for (i = 0; i <= dx; i++)
	{
		image.data[y][x] = color;

		while (e >= 0)
		{
			if (change == 0)
				y += signy;
			else
				x += signx;
			e -= (dx << 1);
		}
		if (change == 0)
			x += signx;
		else
			y += signy;
		//e-=(dy<<1);
		e += (dy << 1);
	}
}   

// draw a line from point-a to point-b in the image
void draw_line(CPoint& a, CPoint& b, C3DImage& image, int SliceNum, unsigned char color)
{
	int xstart = a.m_iX, ystart = a.m_iY;
	int xend = b.m_iX, yend = b.m_iY;
	int i,x,y,dx,dy,e;
	int signx, signy, change, temp;

	if (xstart == xend && ystart == yend)
	{
		image.data[SliceNum][ystart][xstart] = color;
		return;
	}
	if (xend < xstart)
	{
		signx = -1;
		dx = xstart - xend;
	}
	else
	{
		signx = 1;
		dx = xend - xstart;
	}
	if (yend < ystart)
	{
		signy = -1;
		dy = ystart - yend;
	}
	else
	{
		signy = 1;
		dy = yend - ystart;
	}

	if (dy <= dx)
		change = 0;
	else
	{
		change = 1;
		temp = dx; dx = dy; dy = temp;
	}
	e = (dy << 1) - dx;
	x = xstart; y = ystart;
	for (i = 0; i <= dx; i++)
	{
		image.data[SliceNum][y][x] = color;

		while (e >= 0)
		{
			if (change == 0)
				y += signy;
			else
				x += signx;
			e -= (dx << 1);
		}
		if (change == 0)
			x += signx;
		else
			y += signy;
		//e-=(dy<<1);
		e += (dy << 1);
	}
}   
// draw a line from point-a to point-b in the image
void Add_line(CPoint& a, CPoint& b, CImage& image)
{
	int xstart = a.m_iX, ystart = a.m_iY;
	int xend = b.m_iX, yend = b.m_iY;
	int i,x,y,dx,dy,e;
	int signx, signy, change, temp;

	if (xstart == xend && ystart == yend)
	{
		image.data[ystart][xstart] += 1;
		//image.data[ystart][xstart+1] += 1;
		//image.data[ystart][xstart-1] += 1;
		//image.data[ystart+1][xstart] += 1;
		//image.data[ystart-1][xstart] += 1;
		return;
	}
	if (xend < xstart)
	{
		signx = -1;
		dx = xstart - xend;
	}
	else
	{
		signx = 1;
		dx = xend - xstart;
	}
	if (yend < ystart)
	{
		signy = -1;
		dy = ystart - yend;
	}
	else
	{
		signy = 1;
		dy = yend - ystart;
	}

	if (dy <= dx)
		change = 0;
	else
	{
		change = 1;
		temp = dx; dx = dy; dy = temp;
	}
	e = (dy << 1) - dx;
	x = xstart; y = ystart;
	for (i = 0; i <= dx; i++)
	{
		image.data[y][x] += 1;
		//image.data[y][x-1] += 1;
		//image.data[y][x+1] += 1;
		//image.data[y-1][x] += 1;
		//image.data[y+1][x] += 1;

		while (e >= 0)
		{
			if (change == 0)
				y += signy;
			else
				x += signx;
			e -= (dx << 1);
		}
		if (change == 0)
			x += signx;
		else
			y += signy;
		//e-=(dy<<1);
		e += (dy << 1);
	}
}   


// draw a line from point-a to point-b in the image
void DrawThickLine(CPoint* a, CPoint* b, CImage& image, unsigned char color)
{
	int xstart = a->m_iX, ystart = a->m_iY;
	int xend = b->m_iX, yend = b->m_iY;
	int i,x,y,dx,dy,e;
	int signx, signy, change, temp;

	if (xstart == xend && ystart == yend)
	{
		image.data[ystart][xstart] = color;
		return;
	}
	if (xend < xstart)
	{
		signx = -1;
		dx = xstart - xend;
	}
	else
	{
		signx = 1;
		dx = xend - xstart;
	}
	if (yend < ystart)
	{
		signy = -1;
		dy = ystart - yend;
	}
	else
	{
		signy = 1;
		dy = yend - ystart;
	}

	if (dy <= dx)
		change = 0;
	else
	{
		change = 1;
		temp = dx; dx = dy; dy = temp;
	}
	e = (dy << 1) - dx;
	x = xstart; y = ystart;
	for (i = 0; i <= dx; i++)
	{
		image.data[y][x] = color;
		image.data[y][x + 1] = color;
		image.data[y][x - 1] = color;
		image.data[y - 1][x] = color;
		image.data[y + 1][x] = color;

		while (e >= 0)
		{
			if (change == 0)
				y += signy;
			else
				x += signx;
			e -= (dx << 1);
		}
		if (change == 0)
			x += signx;
		else
			y += signy;
		//e-=(dy<<1);
		e += (dy << 1);
	}
}   

// draw a line from point-a to point-b in the image
void DrawThickLine(CPoint* a, CPoint* b, C3DImage& image, unsigned char sliceNumber,
	unsigned char color)
{
	//int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;

	int xstart = a->m_iX, ystart = a->m_iY;
	int xend = b->m_iX, yend = b->m_iY;
	int i,x,y,dx,dy,e;
	int signx, signy, change, temp;

	if (xstart == xend && ystart == yend)
	{
		image.data[sliceNumber][ystart][xstart] = color;
		return;
	}
	if (xend < xstart)
	{
		signx = -1;
		dx = xstart - xend;
	}
	else
	{
		signx = 1;
		dx = xend - xstart;
	}
	if (yend < ystart)
	{
		signy = -1;
		dy = ystart - yend;
	}
	else
	{
		signy = 1;
		dy = yend - ystart;
	}

	if (dy <= dx)
		change = 0;
	else
	{
		change = 1;
		temp = dx; dx = dy; dy = temp;
	}
	e = (dy << 1) - dx;
	x = xstart; y = ystart;
	for (i = 0; i <= dx; i++)
	{
		if (x > 0 && x <iCols - 1 && y> 0 && y < iRows - 1)
		{
			image.data[sliceNumber][y][x] = color;
			image.data[sliceNumber][y][x + 1] = color;
			image.data[sliceNumber][y][x - 1] = color;
			image.data[sliceNumber][y - 1][x] = color;
			image.data[sliceNumber][y + 1][x] = color;
		}

		while (e >= 0)
		{
			if (change == 0)
				y += signy;
			else
				x += signx;
			e -= (dx << 1);
		}
		if (change == 0)
			x += signx;
		else
			y += signy;
		//e-=(dy<<1);
		e += (dy << 1);
	}
}   

void Construct3DLine2(CLNode<CPoint>* head, CLNode<CPoint>* tail)
{
	// if the two points are adjacent, done
	if (abs(head->data->m_iX - tail->data->m_iX) <= 1 &&
		abs(head->data->m_iY -
														 	tail->data->m_iY) <= 1 &&
		abs(head->data->m_iZ -
			tail->data->m_iZ) <= 1)
		return;

	// find the middle point
	CPoint temp = *tail->data;
	float x = static_cast<float>((static_cast<float>(head->data->m_iX + tail->data->m_iX) / 2.0));
	float y = static_cast<float>((static_cast<float>(head->data->m_iY + tail->data->m_iY) / 2.0));
	float z = static_cast<float>((static_cast<float>(head->data->m_iZ + tail->data->m_iZ) / 2.0));

	temp.m_iX = static_cast<int>((x > 0.0) ? (x + 0.5) : (x - 0.5));
	temp.m_iY = static_cast<int>((y > 0.0) ? (y + 0.5) : (y - 0.5));
	temp.m_iZ = static_cast<int>((z > 0.0) ? (z + 0.5) : (z - 0.5));

	// add it to the list
	CLNode<CPoint>* tempNode = new CLNode<CPoint>(&temp);
	tempNode->after = tail;
	tempNode->before = head;
	head->after = tempNode;
	tail->before = tempNode;
	// call recursively on the first half and the second half
	Construct3DLine2(head, tempNode);
	Construct3DLine2(tempNode, tail);
}

void Construct3DLine(CPoint& from, CPoint& to, CDLList<CPoint>& result)
{
	result.AddElementOnEnd(&from);
	result.AddElementOnEnd(&to);


	//recursive call
	Construct3DLine2(result.head, result.tail);
	// update the length of the linked list
	int iCount = 0;
	CLNode<CPoint>* temp = result.head;
	while (temp)
	{
		iCount++;
		temp = temp->after;
	}
	result.length = iCount;
}


// convert the given byte from hex to int. The byte is organized as follows:
// HighHalf, LowHalf. The LowHalf has the order "Order"
#define LowMask 0x0F
#define HighMask 0xF0

int ConvertFromHex(unsigned char aChar, int Order)
{
	int LowValue = (int) (aChar& LowMask);
	int HighValue = (int) (aChar& HighMask);

	int Base = 1;
	for (register int i = 0; i < Order; i++)
		Base = 16 * Base;

	int result = LowValue* Base;

	result += Base * HighValue;

	return result;
}

// The structuring Element is defined global (for now) it is a 3x3 square with
// constant gray scale value of 0


/////////////////////////////////////////////////////////////////////////////
// FUNCTION: GrayScaleErosion
//
// Performs gray scale erosion. The following implementation although less readable
// than the original one (see ../2DTrack/Tools.cpp), but is much faster. For example,
// a 1000x1000 image can be eroded by a 11x11 structuring element in 1.8 seconds.
// compare this with 10 seconds for the original implementation.
void GrayScaleErosion(CImage* InImage, CImage* OutImage)
{
	int iRows = InImage->m_iRows;
	int iCols = InImage->m_iCols;
	int iMargin = StrucElemSize + 1;

	CImage tempImage(iRows, iCols); // to hold temp results
	tempImage.SetColor(0);

	register unsigned char * inData = 0; //inImage->data;
	register unsigned char * outData = 0; //tempImage.data;
	unsigned char * valuePtr;
	register int i, j;
	register int max; //, value;

	int From = - 1 * StrucElemSize;
	int To = StrucElemSize;

	unsigned char * ImageBase = &(InImage->data[0][0]);
	unsigned char * ImageEnd = &(InImage->data[iRows - 1][iCols - 1]);

	// perform erosion on the rows of the image and produce temp image
	for (i = 0; i < iRows; i++)
	{
		inData = &(InImage->data[i][iMargin]);
		outData = &(tempImage.data[i][iMargin]);
		for (j = iMargin; j < iCols - iMargin; j++)
		{
			max = 0;	
			for (register int x = From; x <= To; x++)
			{
				valuePtr = (inData + x);
				if (valuePtr <ImageBase || valuePtr> ImageEnd)
					continue;

				if (*valuePtr > max)
					max = *valuePtr;
			}

			*outData = static_cast<unsigned char>(max);
			outData++;
			inData++;
		}
	}

	ImageBase = &(tempImage.data[0][0]);
	ImageEnd = &(tempImage.data[iRows - 1][iCols - 1]);

	From = -1 * StrucElemSize * iCols;
	To = -1 * From;
	// Perform Erosion on the Columns of the image resulted from the above
	// step
	for (j = 0; j < iCols; j++)
	{
		inData = &(tempImage.data[iMargin][j]);
		outData = &(OutImage->data[iMargin][j]);
		for (i = iMargin; i < iRows - iMargin; i++)
		{
			max = 0;
			for (register int x = From; x <= To; x += iCols)
			{
				valuePtr = (inData + x); ;

				if (valuePtr <ImageBase || valuePtr> ImageEnd)
					continue;

				if (*valuePtr > max)
					max = *valuePtr;
			}

			*outData = static_cast<unsigned char>(max);
			outData += iCols;
			inData += iCols;
		}
	}
}
void GrayScaleErosion(CImage* InImage, CImage* OutImage, CPoint* aDiskPoints,
	int iNumOfPoints)
{
	int iRows = InImage->m_iRows;
	int iCols = InImage->m_iCols;
	int iMargin = static_cast<int>(std::sqrt(static_cast<double>(iNumOfPoints)) + 2);

	register unsigned char * *inData = InImage->data;
	register unsigned char * *outData = OutImage->data;
	register int i, j, k;
	register int max, value;


	// perform erosion on the rows of the image and produce temp image
	for (i = iMargin; i < iRows - iMargin; i++)
	{
		for (j = iMargin; j < iCols - iMargin; j++)
		{
			max = 0;	
			for (k = 0; k < iNumOfPoints; k++)
			{
				value = inData[i + aDiskPoints[k].m_iY][j + aDiskPoints[k].m_iX];
				if (value > max)
					max = value;
			}

			outData[i][j] = static_cast<unsigned char>(max);
		}
	}
}

void GrayScaleErosion(unsigned char** inImageData,
	unsigned char** outImageData)
{
	//int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;

	CImage tempImage(iRows, iCols); // to hold temp results
	tempImage.SetColor(0);

	register unsigned char * inData = 0; //inImage->data;
	register unsigned char * outData = 0; //tempImage.data;
	unsigned char * valuePtr;
	register int i, j;
	register int max;//, value;

	int From = - 1 * StrucElemSize;
	int To = StrucElemSize;

	long ImageLength = iCols* iRows;

	// perform erosion on the rows of the image and produce temp image
	for (i = 0; i < iRows; i++)
	{
		inData = &(inImageData[i][giMARGIN]);
		outData = &(tempImage.data[i][giMARGIN]);
		for (j = giMARGIN; j < iCols - giMARGIN; j++)
		{
			max = 0;	
			for (register int x = From; x <= To; x++)
			{
				valuePtr = (inData + x);
				if (valuePtr <*inImageData ||
					valuePtr> (*inImageData + ImageLength))
					continue;

				if (*valuePtr > max)
					max = *valuePtr;
			}

			*outData = static_cast<unsigned char>(max);
			outData++;
			inData++;
		}
	}

	unsigned char * tempImageBase = &(tempImage.data[0][0]);
	unsigned char * tempImageEnd = &(tempImage.data[iRows - 1][iCols - 1]);

	From = -1 * StrucElemSize * iCols;
	To = -1 * From;
	// Perform Erosion on the Columns of the image resulted from the above
	// step
	for (j = 0; j < iCols; j++)
	{
		inData = &(tempImage.data[giMARGIN][j]);
		outData = &(outImageData[giMARGIN][j]);
		for (i = giMARGIN; i < iRows - giMARGIN; i++)
		{
			max = 0;
			for (register int x = From; x <= To; x += iCols)
			{
				valuePtr = (inData + x); ;

				if (valuePtr <tempImageBase || valuePtr> tempImageEnd)
					continue;

				if (*valuePtr > max)
					max = *valuePtr;
			}

			*outData = static_cast<unsigned char>(max);
			outData += iCols;
			inData += iCols;
		}
	}
}

/////////////////////////////////////////////////////////////////////////////
// FUNCTION: GrayScaleDilation
//
// A fast implementation of the gray scale dilation procedure
void GrayScaleDilation(CImage* InImage, CImage* OutImage)
{
	int iCols = InImage->m_iCols;
	int iRows = InImage->m_iRows;
	int iMargin = StrucElemSize + 1;

	CImage tempImage(iRows, iCols); // to hold temp results
	tempImage.SetColor(0);

	register unsigned char * inData = 0; //inImage->data;
	register unsigned char * outData = 0; //tempImage.data;
	register int i, j;
	register int min, value;

	int From = - 1 * StrucElemSize;
	int To = StrucElemSize;
	// perform erosion on the rows of the image and produce temp image
	for (i = 0; i < iRows; i++)
	{
		inData = &(InImage->data[i][iMargin]);
		outData = &(tempImage.data[i][iMargin]);
		for (j = iMargin; j < iCols - iMargin; j++)
		{
			min = 300;	
			for (register int x = From; x <= To; x++)
			{
				value = *(inData + x);
				if (value < min)
					min = value;
			}

			*outData = static_cast<unsigned char>(min);

			outData++;
			inData++;
		}
	}

	From = -1 * StrucElemSize * iCols;
	To = -1 * From;

	// Perform Erosion on the Columns of the image resulted from the above
	// step
	for (j = 0; j < iCols; j++)
	{
		inData = &(tempImage.data[iMargin][j]);
		outData = &(OutImage->data[iMargin][j]);
		for (i = iMargin; i < iRows - iMargin; i++)
		{
			min = 300;
			for (register int x = From; x <= To; x += iCols)
			{
				value = *(inData + x); ;
				if (value < min)
					min = value;
			}

			*outData = static_cast<unsigned char>(min);
			outData += iCols;
			inData += iCols;
		}
	}
}


/////////////////////////////////////////////////////////////////////////////
// FUNCTION: GrayScaleDilation
//
// A fast implementation of the gray scale dilation procedure
void GrayScaleDilation(CImage* InImage, CImage* OutImage, CPoint* aDiskPoints,
	int iNumOfPoints)
{
	int iCols = InImage->m_iCols;
	int iRows = InImage->m_iRows;
	int iMargin = static_cast<int>(std::sqrt( static_cast<double>(iNumOfPoints) ) + 2); 

	register unsigned char * *inData = InImage->data;
	register unsigned char * *outData = OutImage->data;
	register int i, j, k;
	register int min, value;

	// perform erosion on the rows of the image and produce temp image
	for (i = iMargin; i < iRows - iMargin; i++)
	{
		for (j = iMargin; j < iCols - iMargin; j++)
		{
			min = 300;
			for (k = 0; k < iNumOfPoints; k++)
			{
				value = inData[i + aDiskPoints[k].m_iY][j + aDiskPoints[k].m_iX];
				if (value < min)
					min = value;
			}

			outData[i][j] = static_cast<unsigned char>(min);
		}
	}
}


/////////////////////////////////////////////////////////////////////////////
// FUNCTION: GrayScaleDilation
//
// A fast implementation of the gray scale dilation procedure
void GrayScaleDilation(unsigned char** inImageData,
	unsigned char** outImageData)
{
	//int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;
	CImage tempImage(iRows, iCols); // to hold temp results
	tempImage.SetColor(0);

	register unsigned char * inData = 0; //inImage->data;
	register unsigned char * outData = 0; //tempImage.data;
	register int i, j;
	register int min, value;

	int From = - 1 * StrucElemSize;
	int To = StrucElemSize;
	// perform erosion on the rows of the image and produce temp image
	for (i = 0; i < iRows; i++)
	{
		inData = &(inImageData[i][giMARGIN]);
		outData = &(tempImage.data[i][giMARGIN]);
		for (j = giMARGIN; j < iCols - giMARGIN; j++)
		{
			min = 300;	
			for (register int x = From; x <= To; x++)
			{
				value = *(inData + x);
				if (value < min)
					min = value;
			}

			*outData = static_cast<unsigned char>(min);

			outData++;
			inData++;
		}
	}

	From = -1 * StrucElemSize * iCols;
	To = -1 * From;
	// Perform Erosion on the Columns of the image resulted from the above
	// step
	for (j = 0; j < iCols; j++)
	{
		inData = &(tempImage.data[giMARGIN][j]);
		outData = &(outImageData[giMARGIN][j]);
		for (i = giMARGIN; i < iRows - giMARGIN; i++)
		{
			min = 300;
			for (register int x = From; x <= To; x += iCols)
			{
				value = *(inData + x); ;
				if (value < min)
					min = value;
			}

			*outData = static_cast<unsigned char>(min);
			outData += iCols;
			inData += iCols;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
// Function: GrayScaleOpen
//
// Perform gray scale opening with the specified structuring element
//
// Important:
//  When opening a 1000x1000 image by a structuring element of 11x11 using the 
//  combination "GrayScaleErosion + GrayScaleDilation", the operation took 21.5
//  seconds. Compare this with the combination
//  "GrayScaleErosion3 + GrayScaleDilation3", which took only 3.7 seconds. 
//  => a 5.8 folds decrease in execution time.
void GrayScaleOpen(CImage* inImage, CImage* outImage)
{
	//int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;
	CImage* tempImage = new CImage(iRows, iCols);
	tempImage->SetColor(250);

	GrayScaleErosion(inImage, tempImage);
	GrayScaleDilation(tempImage, outImage);

	//tempImage->Write("ErrodedImage.pgm");

	delete tempImage;
}
void GrayScaleOpen(CImage* inImage, CImage* outImage, CPoint* aDiskPoints,
	int iNumOfPoints)
{
	//int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;
	CImage* tempImage = new CImage(iRows, iCols);
	tempImage->SetColor(250);

	GrayScaleErosion(inImage, tempImage, aDiskPoints, iNumOfPoints);
	GrayScaleDilation(tempImage, outImage, aDiskPoints, iNumOfPoints);

	//tempImage->Write("ErrodedImage.pgm");

	delete tempImage;
}

void GrayScaleOpen(unsigned char** inImageData, unsigned char** outImageData)
{
	//int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;
	CImage* tempImage = new CImage(iRows, iCols);
	tempImage->SetColor(250);

	GrayScaleErosion(inImageData, tempImage->data);
	GrayScaleDilation(tempImage->data, outImageData);

	//tempImage->Write("ErrodedImage.pgm");

	delete tempImage;
}

void GrayScaleClose(unsigned char** inImageData, unsigned char** outImageData)
{
	//int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;
	CImage* tempImage = new CImage(iRows, iCols);
	tempImage->SetColor(0);

	GrayScaleDilation(inImageData, tempImage->data);
	//tempImage->Write("DilatedImage.pgm");
	GrayScaleErosion(tempImage->data, outImageData);

	delete tempImage;
}

void GrayScaleClose(CImage* InImage, CImage* OutImage)
{
	CImage* tempImage = new CImage(InImage->m_iRows, InImage->m_iCols);
	tempImage->SetColor(0);

	GrayScaleDilation(InImage, tempImage);
	//tempImage->Write("DilatedImage.pgm");
	GrayScaleErosion(tempImage, OutImage);

	delete tempImage;
}

void GrayScaleClose(CImage* InImage, CImage* OutImage, CPoint* aDiskPoints,
	int iNumOfElem)
{
	CImage* tempImage = new CImage(InImage->m_iRows, InImage->m_iCols);
	tempImage->SetColor(0);

	GrayScaleDilation(InImage, tempImage, aDiskPoints, iNumOfElem);
	//tempImage->Write("DilatedImage.pgm");
	GrayScaleErosion(tempImage, OutImage, aDiskPoints, iNumOfElem);

	delete tempImage;
}

/////////////////////////////////////////////////////////////////////////////////
// Function: DetectSomas
//
// Purpose:
//    To detect the somas in an input image and write them into the output image
// 
// Logic:
//    perform gray scale opening with a disk structuring elements
void DetectSomas(unsigned char** inImageData, unsigned char** outImageData)
{
	StrucElemSize = (int) (gfWidthSum / (float) giNumOfWidthSumMembers + 0.5);
	//StrucElemSize = StrucElemSize / 2.0;

	GrayScaleClose(inImageData, outImageData);
}

void DetectSomas(CImage* InImage, CImage* OutImage, int ElemSize)
{
	//StrucElemSize = (int) (gfWidthSum / (float) giNumOfWidthSumMembers + 0.5) + 2.0;
	StrucElemSize = ElemSize;

	GrayScaleClose(InImage, OutImage);
}


void DetectSomas(CImage* InImage, CImage* OutImage, CPoint* aDiskPoints,
	int iNumOfPoints)
{
	//StrucElemSize = (int) (gfWidthSum / (float) giNumOfWidthSumMembers + 0.5) + 2.0;

	GrayScaleClose(InImage, OutImage, aDiskPoints, iNumOfPoints);
}

// FUNCTION: ComputeMedianAndStdDev
void EstimateMedianAndStdDev(int* Hist, int& median, float& stdDev)
{
	int HistLow[256];
	int HistHigh[256];
	int Diff[256];
	int StdDevs[256];
	int i;

	memset(HistLow, 0, sizeof(int) * 256);
	memset(HistHigh, 0, sizeof(int) * 256);
	memset(Diff, 0, sizeof(int) * 256);
	memset(StdDevs, 0, sizeof(int) * 256);

	for (i  = 1; i < 256; i++)
		HistLow[i] = Hist[i - 1] + HistLow[i - 1];

	for (i = 254; i >= 0; i--)
		HistHigh[i] = Hist[i + 1] + HistHigh[i + 1];

	for (i = 0; i < 256; i++)
		Diff[i] = abs(HistHigh[i] - HistLow[i]);

	int minValue = Diff[0] + 1;
	for (i = 0; i < 256; i++)
	{
		if (Diff[i] < minValue)
		{
			minValue = Diff[i];
			median = i;
		}
	}

	memset(HistLow, 0, sizeof(int) * 256);
	memset(HistHigh, 0, sizeof(int) * 256);
	memset(Diff, 0, sizeof(int) * 256);

	int StdDev;
	StdDev = - 1;
	for (i = 0 ; i < 256; i++)
	{
		int diff = abs(i - median);
		StdDevs[diff] += Hist[i];
	}

	for (i = 1; i < 256; i++)
		HistLow[i] = StdDevs[i - 1] + HistLow[i - 1];

	for (i = 254; i >= 0; i--)
		HistHigh[i] = StdDevs[i + 1] + HistHigh[i + 1];

	for (i = 0; i < 256; i++)
		Diff[i] = abs(HistHigh[i] - HistLow[i]);

	minValue = Diff[0] + 1;
	for (i = 0; i < 256; i++)
	{
		if (Diff[i] < minValue)
		{
			minValue = Diff[i];
			stdDev = static_cast<float>(i);
		}
	}


	stdDev = static_cast<float>(1.4826 * stdDev);
}


//////////////////////////////////////
// Function: MatrixMultiply
//
// Perfomr matrix multiplication
void MatrixMultiply(double Left[][4], double Right[][4], double Result[][4],
	int x, int y)
{
	register int i, j, k;
	for (i = 0; i < x; i++)
	{
		for (j = 0; j < y; j++)
		{
			Result[i][j] = 0;
			for (k = 0; k < x; k++)
				Result[i][j] += (Left[i][k] * Right[k][j]);
		}
	}
}

void MatrixMultiply(CPoint& aPoint, double Right[][4], CPoint& result)
{
	result.m_iX = static_cast<int>(aPoint.m_iX * Right[0][0] + aPoint.m_iY * Right[1][0] + aPoint.m_iZ * Right[2][0]);
	result.m_iY = static_cast<int>(aPoint.m_iX * Right[0][1] + aPoint.m_iY * Right[1][1] + aPoint.m_iZ * Right[2][1]);
	result.m_iZ = static_cast<int>(aPoint.m_iX * Right[0][2] + aPoint.m_iY * Right[1][2] + aPoint.m_iZ * Right[2][2]);
}

///////////////////////////////////////////////////////////////
// Function: ConstructStructElem
// 
// construct a structuring elements of point indices
void ConstructStructElem(int Radius, CPoint* aDiskPoints, int& iNumOfPoints)
{
	register int i, j;
	register int index = 0;
	int RadiusSquared = Radius* Radius;

	for (i = -1 * Radius; i <= Radius; i++)
	{
		for (j = -1 * Radius; j <= Radius; j++)
		{
			if ((i * i + j * j) <= RadiusSquared)
			{
				aDiskPoints[index].m_iX = i;
				aDiskPoints[index].m_iY = j;
				index++;
			}
		}
	}
	iNumOfPoints = index;
}



/////////////////////////////////////////////////////////////////
// BY Yousef: The following functions are used for soma detection
//            From the 3D image directly
//1- Constructing a Structure Element
void Construct3DStructElem(int Radius, CPoint* aSphrPoints, int& iNumOfPoints)
{	
	register int i, j, k;
	register int index = 0;
	int RadiusSquare = Radius * Radius;

	for (i = -1 * Radius; i <= Radius; i++)
	{
		for (j = -1 * Radius; j <= Radius; j++)
		{
			for (k = -1 * Radius; k <= Radius; k++)
			{
				if ((i * i + j * j + k * k) <= RadiusSquare)
				{
					aSphrPoints[index].m_iX = i;
					aSphrPoints[index].m_iY = j;
					aSphrPoints[index].m_iZ = k;
					index++;
				}
			}			
		}
	}
	iNumOfPoints = index;
}

//2- Detect the somas
// Dilation Function
void GrayScaleDilation3D(C3DImage* InImage, C3DImage* OutImage, CPoint* aSphrPoints,
	int iNumOfPoints, int StructElemSize)
{
	int iCols = InImage->m_iCols;
	int iRows = InImage->m_iRows;
	int iSlcs = InImage->m_iSlices; 
	int iMargin = StructElemSize;
	register int i, j, l, k;
	register int max, value;
	

	// perform Dilation on the rows of the image and produce temp image
	for (i = iMargin; i < iRows - iMargin; i++)
	{
		for (j = iMargin; j < iCols - iMargin; j++)
		{
			for (l = iMargin; l < iSlcs - iMargin; l++)
			{
				max = 0;
				for (k = 0; k < iNumOfPoints; k++)
				{
					value = InImage->data[l+aSphrPoints[k].m_iZ][i + aSphrPoints[k].m_iY][j + aSphrPoints[k].m_iX];
					if(value == 255)
					{
						max = 255;
						break;
					}
					else if (value > max)
						max = value;
				}
				OutImage->data[l][i][j] = static_cast<unsigned char>(max);
			}									
		}
	}
}

// Erosion Function
void GrayScaleErosion3D(C3DImage* InImage, C3DImage* OutImage, CPoint* aSphrPoints,
	int iNumOfPoints, int StructElemSize)
{
	int iCols = InImage->m_iCols;
	int iRows = InImage->m_iRows;
	int iSlcs = InImage->m_iSlices; 
	int iMargin = StructElemSize;	
	register int i, j, l, k;
	register int min, value;


	// perform erosion on the rows of the image and produce temp image
	for (i = iMargin; i < iRows - iMargin; i++)
	{
		for (j = iMargin; j < iCols - iMargin; j++)
		{
			for (l = iMargin; l < iSlcs - iMargin; l++)
			{
				min = 255;
				for (k = 0; k < iNumOfPoints; k++)
				{
					value = InImage->data[l+aSphrPoints[k].m_iZ][i + aSphrPoints[k].m_iY][j + aSphrPoints[k].m_iX];
					if(value == 0)
					{
						min = 0;
						break;
					}
					else if (value < min)
						min = value;
				}
				OutImage->data[l][i][j] = static_cast<unsigned char>(min);
			}
		}
	}	
}
// Opening function 
void GrayScaleOpen3D(C3DImage* InImage, C3DImage* outImage, CPoint* aSphrPoints,
	int iNumOfPoints, int StructElemSize)
{
	C3DImage* tempImage = new C3DImage(InImage->m_iSlices, InImage->m_iRows, InImage->m_iCols);		
	//StructElemSize = (int) (StructElemSize/2);
	GrayScaleErosion3D(InImage, tempImage, aSphrPoints, iNumOfPoints, StructElemSize);
	GrayScaleDilation3D(tempImage, outImage, aSphrPoints, iNumOfPoints, StructElemSize);
	delete tempImage;
}
//Closing function
void GrayScaleClose3D(C3DImage* InImage, C3DImage* OutImage, CPoint* aSphrPoints,
	int iNumOfElem, int StructElemSize)
{	
	C3DImage* tempImage = new C3DImage(InImage->m_iSlices, InImage->m_iRows, InImage->m_iCols);	
	GrayScaleDilation3D(InImage, tempImage, aSphrPoints, iNumOfElem, StructElemSize);
	GrayScaleErosion3D(tempImage, OutImage, aSphrPoints, iNumOfElem, StructElemSize);	
	delete tempImage;
}
//Calling function
void Detect3DSomas(C3DImage* InImage, C3DImage* OutImage, CPoint* aSphrPoints,
	int iNumOfPoints, int StructElemSize)
{			
	GrayScaleClose3D(InImage, OutImage, aSphrPoints, iNumOfPoints, StructElemSize);	
	GrayScaleOpen3D(InImage, OutImage, aSphrPoints, iNumOfPoints, StructElemSize);	
}

/////////////////////////////////////////////////
// Function: DrawBall
// 
// Draw a ball (representing a soma).
void DrawBall(C3DImage& anImage, CPoint& center, int radius, unsigned char color)
{
	// for each pixel in the cube containing the ball, if the pixel is within
	// a distance r from the center, it belongs to the ball
	register int i, j, k;
	int sideLength = (int) (sqrt(2.0) * (float) radius + 0.5);
	int rSquared = radius* radius;
	int dx, dy, dz;
	for (i = center.m_iX - sideLength; i <= center.m_iX + sideLength; i++)
	{
		for (j = center.m_iY - sideLength; j <= center.m_iY + sideLength; j++)
		{
			for (k = center.m_iZ - sideLength;
				k <= center.m_iZ + sideLength;
				k++)
			{
				dx = (i - center.m_iX) * (i - center.m_iX);
				dy = (j - center.m_iY) * (j - center.m_iY);
				dz = (k - center.m_iZ) * (k - center.m_iZ);
				if ((dx + dy + dz) < rSquared)
					anImage.data[k][j][i] = color;
			}
		}
	}
}


////////////////////////////////
// Function: DisplayTemplates
//
// write all the templates to an image and display the image in the given file
void DisplayVectors(char* fName)
{
	int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;
	/////////////////////
	// working variables
	register int i, j;

	// the row position of each template. we want to print 8 rows of templates
	// and vectors each seprated by 60 pixels
	// bug?? amri: 11-30-02 unused variable
//	int rowPos[] =
//	{
//		40, 100, 160, 220, 280, 340, 400, 460
//	};
	// the col positions of each templage. We want to print 8 cols of templates
	// and vectors each separated by 80 pixels
	int colPos[] =
	{
		30, 70, 110, 150, 190, 230, 270, 310, 350, 390, 440, 490, 540, 590,
		640, 690
	};

	// create an empty image of the specified size
	C3DImage tempImage(iSlices, iRows, iCols);
	tempImage.ReadHeader("BioRadHeader.txt");
	// write each of the left and right templates into the image file at different
	// direction
	int index = 0;
	for (i = 9; i < 16; i++)
	{
		for (j = 8; j < 9; j++)
		{
			index = i;
			if (i > 8)
				index = i - 8;
			gVectorsArray[i][j]->Position(colPos[index], colPos[index], 50);			
			gVectorsArray[i][j]->DisplayVector(tempImage, 255);
		}
	}

	tempImage.Write(fName);
}
////////////////////////////////
// Function: DisplayTemplates
//
// write all the templates to an image and display the image in the given file
void DisplayTemplates(char* fName)
{
	int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;
	/////////////////////
	// working variables
	register int i, j;
	// the row position of each template. we want to print 8 rows of templates
	// and vectors each seprated by 60 pixels

	// bug?? amri: 11-30-02 unused variable
//	int rowPos[] =
//	{
//		40, 100, 160, 220, 280, 340, 400, 460
//	};
	// the col positions of each templage. We want to print 8 cols of templates
	// and vectors each separated by 80 pixels
	int colPos[] =
	{
		40, 80, 120, 160, 200, 240, 280, 320, 360, 390, 440, 490, 540, 590,
		640, 690
	};

	// create an empty image of the specified size
	C3DImage tempImage(iSlices, iRows, iCols);
	tempImage.ReadHeader("BioRadHeader.txt");

	int index = 0;
	// write each of the left and right templates into the image file at different
	// direction
	for (i = 0; i < NumOfDirections; i++)
	{
		for (j = 0; j < NumOfDirections; j++)
		{
			gHLeftTemplatesArray[i][j]->AssociateWithImage(tempImage);
			gVLeftTemplatesArray[i][j]->AssociateWithImage(tempImage);
			gHRightTemplatesArray[i][j]->AssociateWithImage(tempImage);
			gVRightTemplatesArray[i][j]->AssociateWithImage(tempImage);

			index = j;
			if (j > 8)
				index = j - 8;

			gHLeftTemplatesArray[i][j]->Position(colPos[index],
											colPos[index],
											50);
			gVLeftTemplatesArray[i][j]->Position(colPos[index],
											colPos[index],
											50);
			gHRightTemplatesArray[i][j]->Position(colPos[index],
										 	colPos[index],
										 	50);
			gVRightTemplatesArray[i][j]->Position(colPos[index],
										 	colPos[index],
										 	50);		

			// display the left and right templates in the given image
			gHLeftTemplatesArray[i][j]->DisplayTemplate(tempImage);
			gVLeftTemplatesArray[i][j]->DisplayTemplate(tempImage);
			gHRightTemplatesArray[i][j]->DisplayTemplate(tempImage);
			gVRightTemplatesArray[i][j]->DisplayTemplate(tempImage);
		}
	}

	tempImage.Write(fName);
}

//////////////////////////////////////////////////////////////////////////////////
// Function: strlower
//
// Purpose:  convert the given string to lower case
//
// Called By: Any
//
void strlower(char* pchString)
{
	if (pchString)
	{
		char* pch1 = pchString;
		while (*pch1 != '\0')
		{
			*pch1 = static_cast<char>(tolower(*pch1));
			pch1++;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
// Function: SetConfigParameters
// Purpose: sets the value of the program global parameters from the 
// configuration object

void SetConfigParameters()
{
	//char* pchValue = NULL;

	///////////////////////////////////
	//// Output configurations
	//if ((pchValue = gConfig.GetStringValue("Output.Image.Format")) != NULL)
	//{
	//	strcpy(cfg_output_image_format, pchValue);
	//}
	//else
	//{
	//	cout <<
	//		"ConfigError::Output.Image.Format not found. Defaulting to tif.\n";
	//}

	//if ((pchValue = gConfig.GetStringValue("Output.Projections")) != NULL)
	//{
	//	if (strcmp(pchValue, "YES") == 0)
	//		cfg_output_projections = 1;
	//	else if (strcmp(pchValue, "NO") == 0)
	//		cfg_output_projections = 0;
	//	else
	//	{
	//		cout <<
	//			"ConfigError::Output.Projections is set to invalid value\n"
	//			<<
	//			"Output.Projections is set to yes" <<
	//			endl;
	//		cfg_output_projections = 1;
	//	}
	//}
	//else
	//{
	//	cout << "ConfigError::Output.Projections not found.\n"
	//		<< "Output.Projections is set to yes" << endl;
	//	cfg_output_projections = 1;
	//}

	//if ((pchValue = gConfig.GetStringValue("Output.Image.FGColor")) != NULL)
	//{
	//	if (strcmp(pchValue, "YES") == 0)
	//		cfg_output_image_fgcolor = 1;
	//	else if (strcmp(pchValue, "NO") == 0)
	//		cfg_output_image_fgcolor = 0;
	//	else
	//	{
	//		cout <<
	//			"ConfigError::Output.Image.FGColor is set to invalid value\n"
	//			<<
	//			"Output.Image.FGColor is set to yes" <<
	//			endl;
	//		cfg_output_image_fgcolor = 1;
	//	}
	//}
	//else
	//{
	//	cout << "ConfigError::Output.Image.FGColor not found.\n"
	//		<< "Output.Image.FGColor is set to yes" << endl;
	//	cfg_output_image_fgcolor = 1;
	//}

	//if ((pchValue = gConfig.GetStringValue("Output.Image.BGColor")) != NULL)
	//{
	//	if (strcmp(pchValue, "YES") == 0)
	//		cfg_output_image_bgcolor = 1;
	//	else if (strcmp(pchValue, "NO") == 0)
	//		cfg_output_image_bgcolor = 0;
	//	else
	//	{
	//		cout <<
	//			"ConfigError::Output.Image.BGColor is set to invalid value\n"
	//			<<
	//			"Output.Image.BGColor is set to yes" <<
	//			endl;
	//		cfg_output_image_bgcolor = 1;
	//	}
	//}
	//else
	//{
	//	cout << "ConfigError::Output.Image.BGColor not found.\n"
	//		<< "Output.Image.BGColor is set to yes" << endl;
	//	cfg_output_image_bgcolor = 1;
	//}

	//if (giDetectSoma)
	//{
	//	if ((pchValue = gConfig.GetStringValue("Output.Soma.Draw.Trees")) !=
	//		NULL)
	//	{
	//		if (strcmp(pchValue, "YES") == 0)
	//			cfg_output_soma_draw_trees = 1;
	//		else if (strcmp(pchValue, "NO") == 0)
	//			cfg_output_soma_draw_trees = 0;
	//		else
	//		{
	//			cout <<
	//				"ConfigError::Output.Soma.Draw.Trees is set to invalid value\n"
	//				<<
	//				"Output.Soma.Draw.Trees is set to yes" <<
	//				endl;
	//			cfg_output_soma_draw_trees = 1;
	//		}
	//	}
	//	else
	//	{
	//		cout << "ConfigError::Output.Soma.Draw.Trees not found.\n"
	//			<< "Output.Soma.Draw.Trees is set to no" << endl;
	//		cfg_output_soma_draw_trees = 0;
	//	}
	//}

	//if ((pchValue = gConfig.GetStringValue("Mode.Debug")) != NULL)
	//{
	//	if (strcmp(pchValue, "YES") == 0)
	//	{
	//		cfg_mode_debug = 1;
	//		cout << "Mode Change: Operating in Debug Mode" << endl;
	//	}
	//	else if (strcmp(pchValue, "NO") == 0)
	//		cfg_mode_debug = 0;
	//	else
	//	{
	//		cout << "ConfigError::Mode.Debug is set to invalid value\n"
	//			<< "Mode.Debug is set to NO" << endl;
	//		cfg_mode_debug = 0;
	//	}
	//}
	//else
	//{
	//	cfg_mode_debug = 0;
	//}

	//if ((pchValue = gConfig.GetStringValue("Seeds.FromFile")) != NULL)
	//{
	//	cfg_seedsfromfile = 1;
	//	cout << "Config: Reading Seedpoints from file initiated" << endl;
	//	cfg_seedsfilename = pchValue;
	//}
	//else
	//{
	//	cout << "\n Seeds.FromFile option invalid: "
	//		<< cfg_seedsfilename << " not found." << endl;
	//}
}

void FreeResources()
{
	register int i, j;

	for (i = 0; i < NumOfDirections; i++)
	{
		for (j = 0; j < NumOfDirections; j++)
		{
			delete gHLeftTemplatesArray[i][j];
			delete gVLeftTemplatesArray[i][j];
			delete gHRightTemplatesArray[i][j];
			delete gVRightTemplatesArray[i][j];

			delete gVectorsArray[i][j];
		}
	}

	delete [] gapArrayOfSeedPoints[0];	
	delete [] gapArrayOfSeedPoints;
}

///////////////////////////////////////////////////////////////////////////////
// Function: ProcessCommandLine
// Called By: main
//
// Process command line and read the configuration file. See the comments 
// with the function "main" for the command line format
void ProcessCommandLine(int argc, char* argv[])
{
	// Set the program configuration parameters from the config object
	SetConfigParameters();
}

// Returns a median of a priority queue 
//	
// (only sorts the upper half of the queue)
// insertion time to the priority queue is O(log n)
// access time to the priority queue top element is constant time
int median(priority_queue<short> responses)
{
	int val = 0;
	int middle = responses.size() / 2;
	if (responses.size() == 1) 
	{
		return responses.top();
	}
	else
	{
		if (responses.size() % 2 == 0)
		{
			while (static_cast<int>(responses.size()) > middle)
			{
				val = responses.top();
				responses.pop();
			}
			return (responses.top() + val) / 2;
		}
		else
		{
			while (static_cast<int>(responses.size()) > middle + 1)
			{
				responses.pop();
			}
			return responses.top();
		}
	}
}

// Returns a median of a list
//	
// Although insertion time to a list is constant compared to log n for a priority_queue
// the sorting process involved slows it down to 50% compared to median calculator above 
// that uses priority queue
int median(list<int>& responses)
{
	int count;
	int middle = responses.size() / 2;
	float fmiddle = static_cast<float>(static_cast<float>(responses.size()) / 2.0);
	list<int>::iterator i;
	if (responses.size() == 1)
		return responses.front();
	else
	{
		responses.sort();
		if (responses.size() % 2 == 0)
		{
			int temp1, temp2;
			count = 1;
			i = responses.begin();
			while (count < fmiddle)
			{
				i++;
				count++;
			}
			temp1 = *i;
			i++;
			temp2 = *i;

			return (temp1 + temp2) / 2;
		}
		else
		{
			i = responses.begin();
			count = 0;
			while (count < middle)
			{
				i++;
				count++;
			}
			return *i;
		}
	}
}

float Median(list<float>& responses)
{
	int count;
	int middle = responses.size() / 2;
	float fmiddle = static_cast<float>(static_cast<float>(responses.size()) / 2.0);
	list<float>::iterator i;
	if (responses.size() == 1)
		return responses.front();
	else
	{
		responses.sort();
		if (responses.size() % 2 == 0)
		{
			int temp1, temp2;
			count = 1;
			i = responses.begin();
			while (count < fmiddle)
			{
				i++;
				count++;
			}
			temp1 = static_cast<int>(*i);
			i++;
			temp2 = static_cast<int>(*i);

			return static_cast<float>((temp1 + temp2) / 2.0);
		}
		else
		{
			i = responses.begin();
			count = 0;
			while (count < middle)
			{
				i++;
				count++;
			}
			return *i;
		}
	}
	return 0;
}

unsigned char __Median(list<unsigned char>& responses)
{
	int count;
	int middle = responses.size() / 2;
	float fmiddle = static_cast<float>(static_cast<float>(responses.size()) / 2.0);
	list<unsigned char>::iterator i;
	if (responses.size() == 1)
		return responses.front();
	else
	{
		responses.sort();
		if (responses.size() % 2 == 0)
		{
			int temp1, temp2;
			count = 1;
			i = responses.begin();
			while (count < fmiddle)
			{
				i++;
				count++;
			}
			temp1 = *i;
			i++;
			temp2 = *i;

			return static_cast<unsigned char>((temp1 + temp2) / 2.0);
		}
		else
		{
			i = responses.begin();
			count = 0;
			while (count < middle)
			{
				i++;
				count++;
			}
			return *i;
		}
	}
}

int ReverseTheta(int dir)
{
	if (dir == 0)
		return NumOfDirections / 2;
	if (dir == NumOfDirections / 2)
		return 0;
	if (dir > NumOfDirections / 2)
		return dir - NumOfDirections / 2;
	if (dir < NumOfDirections / 2)
		return dir + NumOfDirections / 2;
	return -1;
}

///////////////////////////////////////////////////////////////////////////////
// Function: FindAvgeDirection
//
// Find and return the avge of the two directions
int FindAvgeDirection(int dir1, int dir2)
{
	int result = 0;
	int dirDistance = DirDistance(dir1, dir2);
	if (((dir1 + dirDistance) % NumOfDirections) == dir2)
		result = static_cast<int>(dir1 + dirDistance / 2.0);
	else
		result = static_cast<int>(dir2 + dirDistance / 2.0);

	return ((result >= NumOfDirections) ? (result - NumOfDirections) : result);
}

// algorithm for centering seed points based on best template locations
// - amri july 20th 2002			
// we only center based on templates that does matter (that changes)
//
//								x best vertical left
//
//
//    best horizontal left x					 x best horizontal right
//  			
//  							x best vertical right
// 
// Say these are (x,y) coordinates
// - we calculate the middle y coordinate just using vertical points
//   because the the horizontal points shares the same y value
// - in this case, only vertical points are used for centering the y value,
//   and only horizontal points are used for centering the x value
void GetCenterLocation(CPoint& center_point, const CPoint& HR,
	const CPoint& HL, const CPoint& VR, const CPoint& VL)
{
	// previous centering method	
	//	CPoint tempPoint;
	//	int XYResponseSum = BestHLeftPoint.m_iValue + BestHRightPoint.m_iValue;
	//	int ZResponseSum  = BestVLeftPoint.m_iValue + BestVRightPoint.m_iValue;
	//		float AvgeX = (float) (BestHLeftPoint.m_iX  * BestHLeftPoint.m_iValue  + 
	//			BestHRightPoint.m_iX * BestHRightPoint.m_iValue) /
	//			XYResponseSum;
	//		
	//		float AvgeY = (float) (BestHLeftPoint.m_iY  * BestHLeftPoint.m_iValue  + 
	//			BestHRightPoint.m_iY * BestHRightPoint.m_iValue) /
	//			XYResponseSum;
	//		
	//		float AvgeZ = (float) (BestVLeftPoint.m_iZ  * BestVLeftPoint.m_iValue  + 
	//			BestVRightPoint.m_iZ * BestVRightPoint.m_iValue) /
	//			ZResponseSum;
	//		
	//		aPoint->m_iX   = (int) (( AvgeX + (float) aPoint->m_iX)/2.0 + 0.5);
	//		aPoint->m_iY   = (int) (( AvgeY + (float) aPoint->m_iY)/2.0 + 0.5);
	//		aPoint->m_iZ   = (int) (( AvgeZ + (float) aPoint->m_iZ)/2.0 + 0.5);

	bool V_dx, V_dy, V_dz, H_dx, H_dy, H_dz;
	int center_X, center_Y, center_Z;

	V_dx = !((VL.m_iX - VR.m_iX) == 0);
	V_dy = !((VL.m_iY - VR.m_iY) == 0);
	V_dz = !((VL.m_iZ - VR.m_iZ) == 0);
	H_dx = !((HL.m_iX - HR.m_iX) == 0);
	H_dy = !((HL.m_iY - HR.m_iY) == 0);
	H_dz = !((HL.m_iZ - HR.m_iZ) == 0);

	// get the center X
	if (V_dx && H_dx)
	{
		center_X = Round((float) (HL.m_iX + HR.m_iX + VL.m_iX + VR.m_iX) / 4.0);
	}
	else
	{
		if (V_dx && !H_dx)
		{
			center_X = Round((float) (VL.m_iX + VR.m_iX) / 2.0);
		}
		else
		{
			if (!V_dx && H_dx)
			{
				center_X = Round((float) (HL.m_iX + HR.m_iX) / 2.0);
			}
			else
			{
				center_X = Round((float) (VR.m_iX + HR.m_iX) / 2.0);
			}
		}
	}

	// get the center Y
	if (V_dy && H_dy)
	{
		center_Y = Round((float) (HL.m_iY + HR.m_iY + VL.m_iY + VR.m_iY) / 4.0);
	}
	else
	{
		if (V_dy && !H_dy)
		{
			center_Y = Round((float) (VL.m_iY + VR.m_iY) / 2.0);
		}
		else
		{
			if (!V_dy && H_dy)
			{
				center_Y = Round((float) (HL.m_iY + HR.m_iY) / 2.0);
			}
			else
			{
				center_Y = Round((float) (VR.m_iY + HR.m_iY) / 2.0);
			}
		}
	}

	// get the center Z
	if (V_dz && H_dz)
	{
		center_Z = Round((float) (HL.m_iZ + HR.m_iZ + VL.m_iZ + VR.m_iZ) / 4.0);
	}
	else
	{
		if (V_dz && !H_dz)
		{
			center_Z = Round((float) (VL.m_iZ + VR.m_iZ) / 2.0);
		}
		else
		{
			if (!V_dz && H_dz)
			{
				center_Z = Round((float) (HL.m_iZ + HR.m_iZ) / 2.0);
			}
			else
			{
				center_Z = Round((float) (VR.m_iZ + HR.m_iZ) / 2.0);
			}
		}
	}
	center_point.m_iX = center_X;
	center_point.m_iY = center_Y;
	center_point.m_iZ = center_Z;
}

void GetCenterDirection(CPoint& center_point, const CPoint& HR,
	const CPoint& HL, const CPoint& VR, const CPoint& VL)
{
	//	int hdir = HL.m_iHDir;
	//	if(HR.m_iValue > HL.m_iValue)
	//		hdir = HR.m_iHDir;

	//	int vdir = VL.m_iVDir;
	//	if(VR.m_iValue > VL.m_iValue)
	//		vdir = VR.m_iVDir;

	int hdir = FindAvgeDirection(HR.m_iHDir, HL.m_iHDir);
	int vdir = FindAvgeDirection(VR.m_iVDir, VL.m_iVDir);
	center_point.m_iHDir = static_cast<unsigned char>(hdir);
	center_point.m_iVDir = static_cast<unsigned char>(vdir);
}

template <typename _Tp>
inline const _Tp& median_(const _Tp& __a, const _Tp& __b, const _Tp& __c)
{
	// concept requirements
	if (__a < __b)
		if (__b < __c)
			return __b; // a b c
		else 
		{
			if (__a < __c)
				return __c; // a c b
			else
				return __a; // c a b
		}  
	else {
		if (__a < __c) 
			return __a;  // b a c
		else 
		{
			if (__b < __c)
				return __c; // b c a
			else
				return __b; // c b a 
		}
	}
}

int median__(int a, int b, int c)
{
	return median_(a,b,c);
}

float median__(float a, float b, float c)
{
	return median_(a,b,c);
}

short median__(short a, short b, short c)
{
	return median_(a,b,c);
}

