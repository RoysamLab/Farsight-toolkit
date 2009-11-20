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

// File: Cimage.h
// Last Updated: 1-8-97
//
// This file contains the declaration of a simple image class to read
// and write PGM images. It is done in a hurry so be warned
//
// Image Types:
//	P5: gray scale raw data. The header looks like this:
// "P5
//  # <Comments>
//  640 480
//  255\n\n"

#ifndef CImage_h
#define CImage_h

#include <string.h> /* memset() */

using namespace std;

// forward declration
class C3DImage;

enum PgmPic
{
	pgm		= 0,
	pic		= 1
};

class CImage
{
public:
	// default CTOR, do some initialization. 
	CImage();
	// copy CTOR, create a new image and intialize it to argument
	CImage(CImage& image);
	// CTOR from a file,
	CImage(const char* fName);
	// Create an image from a data matrix.
	CImage(unsigned char* temp, int x, int y);
	// create an empty image of the specified size
	CImage(int rows, int cols);
	// DTOR
	~CImage();
	void FillImage(unsigned char* matrix, int x, int y);

	// fill the area between the two given profiles. 
	//void Fill(CProfile *from, CProfile *to, int color);

	// Draw a vessel boundaries on the image. The boundaries are
	// given as a link list of profiles.
	//void DrawVessel(CDLList<CProfile> *aVessel, int color = 255);
	// for each vessel pixel, add one to my corresponding pixel
	//void SumVessel(CDLList<CProfile> *aVessel, int color);
	void ThresholdImage(int threshold, int flag = 0);

	// allocate space for the image
	unsigned char* AllocateSpace();

	// mark the given point in the image
	inline void MarkPoint(CPoint* aPoint, unsigned char color = 255)
	{
		data[aPoint->m_iY][aPoint->m_iX] = color;
	}

	inline void MarkPointXY(CPoint* aPoint, unsigned char color = 255)
	{
		data[aPoint->m_iY][aPoint->m_iX] = color;
	}

	inline void MarkPointXZ(CPoint* aPoint, unsigned char color = 255)
	{
		data[aPoint->m_iZ][aPoint->m_iX] = color;
	}

	inline void MarkPointYZ(CPoint* aPoint, unsigned char color = 255)
	{
		data[aPoint->m_iZ][aPoint->m_iY] = color;
	}

	inline void MarkCrosshairXY(CPoint* aPoint, unsigned char color = 255)
	{
		data[aPoint->m_iY][aPoint->m_iX] = color;
		data[aPoint->m_iY][aPoint->m_iX + 1] = color;
		data[aPoint->m_iY][aPoint->m_iX - 1] = color;
		data[aPoint->m_iY - 1][aPoint->m_iX] = color;
		data[aPoint->m_iY + 1][aPoint->m_iX] = color;
	}

	inline void MarkCrosshairXZ(CPoint* aPoint, unsigned char color = 255)
	{
		data[aPoint->m_iZ][aPoint->m_iX] = color;
		data[aPoint->m_iZ][aPoint->m_iX + 1] = color;
		data[aPoint->m_iZ][aPoint->m_iX - 1] = color;
		data[aPoint->m_iZ - 1][aPoint->m_iX] = color;
		data[aPoint->m_iZ + 1][aPoint->m_iX] = color;
	}

	inline void MarkCrosshairYZ(CPoint* aPoint, unsigned char color = 255)
	{
		data[aPoint->m_iZ][aPoint->m_iY] = color;
		data[aPoint->m_iZ][aPoint->m_iY + 1] = color;
		data[aPoint->m_iZ][aPoint->m_iY - 1] = color;
		data[aPoint->m_iZ - 1][aPoint->m_iY] = color;
		data[aPoint->m_iZ + 1][aPoint->m_iY] = color;
	}

	// return 1 if the point lies withing the image boundaries, 0 otherwise
	inline int ValidPoint(CPoint& aPoint)
	{
		return (aPoint.m_iX >= 0 &&
			aPoint.m_iX < m_iCols &&
			aPoint.m_iY >= 0 &&
			aPoint.m_iY < m_iRows);
	}

	// return a reference to the given point in the image
	unsigned char& operator[](CPoint& rhs);
	unsigned char& operator[](CPoint* rhs);

	// read an image from a file
	int Read(const string& fName);
	// write an image to afile
	int Write(const string& fName);
	// write an image into a TIFF file format
	int WriteTIFF(const string& fName);
	// write the data to a file in a format readabl by matlab
	void WriteDataMatLabFormat(const char* fName);
	// take the negative of the image
	void NegateImage();
	// replace all pixels with intensitied higher than the given pixel
	// value with the value.
	void RemovePixels(int pixelValue);

	// clear the contents of this image to all zero
	inline void Clear()
	{
		memset(data[0], 0, m_iRows * m_iCols);
	}
	inline void SetColor(int color)
	{
		memset(data[0], color, m_iRows * m_iCols);
	}

	void ProcessHeader(istream& inFile);
	void Histogram(int*);

	void ComputeStatistics(int&, float&);
	///////////////////////////////
	// Data.

	int m_iRows;				// number of rows in an image
	int m_iCols;				// number of cols in an image
	unsigned char** data;	// image data array 
	PgmPic type;
	unsigned char ImageHeader[100];

	int MinimaIndex;
};

#endif
