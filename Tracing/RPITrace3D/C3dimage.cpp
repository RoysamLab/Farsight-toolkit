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

// File: CImage.cpp
// Last Updated: 1-8-97
//
// This file contains the implementation of a simple image class to read
// and write PGM images. It is done in a hurry so be warned
//
//#pragma warning(disable:4786)  // disable STL-related warnings
//#pragma warning(disable:4710)  // compiler didn't inline function that got selected for inline expansion
//#pragma warning(disable:4514)  // The optimizer removed an inline function that is not called
//#pragma warning(disable:4702)  // unreachable STLport code

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <list>
#include <queue>
#include <sstream>
#include <stdlib.h> /* exit */
#include <string.h> /* memset() */

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

// default CTOR
C3DImage::C3DImage() : m_iSlices(0), m_iRows(0), m_iCols(0), m_iPadding(0), data(0),
	type(pgm)
{
	memset(m_aiHistogram, 0, sizeof(int) * 256);
}

// copy CTOR
C3DImage::C3DImage(const C3DImage& image) : m_iSlices(image.m_iSlices),
  m_iRows(image.m_iRows), m_iCols(image.m_iCols), m_iPadding(image.m_iPadding), data(0),
	type(image.type)
{
	if (m_iRows != 0 && m_iCols != 0 && m_iSlices != 0)
	{
		unsigned char * buffer = AllocateSpace();
		memcpy(buffer, & (image.data[0][0][0]), m_iSlices * m_iRows * m_iCols);
	}
	memset(m_aiHistogram, 0, sizeof(int) * 256);
	// copy the header as well
	memset(m_auchPicHeader, 0, sizeof(unsigned char) * 76);
	memcpy(m_auchPicHeader,
		& (image.m_auchPicHeader),
		sizeof(unsigned char) * 76);
}

// CTOR from a file
C3DImage::C3DImage(const string& fName, const string& atype) : data(0)
{
	m_iPadding = gConfig.GetImagePadding();
	if (atype == "pgm")
		type = pgm;
	else if (atype == "pic")
		type = pic;

	cout << "\tReading Input Image: " << fName << " .... " << flush;

	if (type == pic)
	{
		ifstream inFile;
		inFile.open(fName.c_str(), ios::binary);
		if (!inFile)
		{
			cout << "Could not open the file: " << fName << endl;
			exit(0);
		}
		unsigned char achTemp[77];
		inFile.read(reinterpret_cast<char*> (achTemp), 76);
		m_iCols = ConvertFromHex(achTemp[0], 0) + ConvertFromHex(achTemp[1], 2);
		m_iRows = ConvertFromHex(achTemp[2], 0) + ConvertFromHex(achTemp[3], 2);
		m_iSlices = ConvertFromHex(achTemp[4], 0) + ConvertFromHex(achTemp[5], 2);
		m_iSlices += (2 * m_iPadding);
	}
	
	AllocateSpace();
	// call 'ReadImage()' to do the work
	Read(fName);
	memset(m_aiHistogram, 0, sizeof(int) * 256);
}

// CTOR create an empty image of the specified size
C3DImage::C3DImage(int z, int y, int x) : m_iSlices(z), m_iRows(y), m_iCols(x), m_iPadding(0)
{
	unsigned char * buffer = AllocateSpace();
	memset(buffer, 0, m_iRows * m_iCols * m_iSlices);
	memset(m_aiHistogram, 0, sizeof(int) * 256);
	memset(m_auchPicHeader, 0, sizeof(unsigned char) * 76);

	int tempCols = m_iCols;
	int tempRows = m_iRows;
	int tempSlices = m_iSlices;
	m_auchPicHeader[0] = static_cast<char>(tempCols & 0x00FF);
	tempCols = m_iCols >> 8;
	m_auchPicHeader[1] = static_cast<char>(tempCols & 0x00FF);

	m_auchPicHeader[2] = static_cast<char>(tempRows & 0x00FF);
	tempRows = m_iRows >> 8;
	m_auchPicHeader[3] = static_cast<char>(tempRows & 0x00FF);

	m_auchPicHeader[4] = static_cast<char>(tempSlices & 0x00FF);
	tempSlices = m_iSlices >> 8;
	m_auchPicHeader[5] = static_cast<char>(tempSlices & 0x00FF);


	// image header added by Amri
	// If this is set to 1, then each pixel is 8-bits; otherwise pixels
	// are 16-bits
	m_auchPicHeader[14] = 0x0001;
	// valid .PIC file, header[54-55]=12345
	m_auchPicHeader[54] = 0x0039;
	m_auchPicHeader[55] = 0x0030;
}

// allocate the needed space and make the data member "data" point
// to the correct locations. This way, we can access the image as a single
// chunck of memory, as indivual slices, m_iRows, and/or cells.
unsigned char* C3DImage::AllocateSpace()
{
	unsigned char * buffer = 0;
	if (m_iRows != 0 && m_iCols != 0 && m_iSlices != 0)
	{
		int sliceSize = m_iRows* m_iCols;
		//int imageSize = sliceSize * m_iSlices;
		// allocate extra space for 5 slices before the first slice and another 5 after 
		// the last slice

		//m_iSlices += (2 * giPADDING);

		int imageSize = sliceSize* m_iSlices;

		if ((buffer = new unsigned char [imageSize]) == NULL)
		{
			cout << "\n C3DImage::C3DImage => Out of Memory" << endl;
			exit(0);
		}

		// fill buffer with zeros
		memset(buffer, 0, imageSize);
		//buffer[imageSize] = '\0';

		if ((data = new unsigned char * *[m_iSlices]) == NULL)
		{
			cout << "\n C3DImage::CTOR, Out of Memory " << endl;
			exit(1);
		}

		register int i, j;

		// assign the slice pointers
		for (i = 0; i < m_iSlices; i++)
		{
			if ((data[i] = new unsigned char * [m_iRows]) == NULL)
			{
				cout << "\n C3DImage::Copy_CTOR, Out of Memory " << endl;
				exit(1);
			}
		}

		for (i = 0; i < m_iSlices; i++)
		{
			for (j = 0; j < m_iRows; j++)
			{
				data[i][j] = buffer + (i * sliceSize + j * m_iCols);
			}
		}
	}

	return buffer;
}


// DTOR
C3DImage::~C3DImage()
{
	unsigned char * imageData = &(data[0][0][0]);

	if (m_iRows > 0 && m_iCols > 0 && m_iSlices > 0)
	{
		// delete the space holding the 3D image data
		delete [] imageData;

		// delete the pointers to the data
		for (register int i = 0; i < m_iSlices; i++)
		{
			if (data[i])
				delete [] data[i];
		}

		// and its pointers
		delete [] data;
	}
}

void C3DImage::RevertToFile(const char* fName, int)
{
	/*
	assert(fName);	
	char tempFName[100];
	int sliceSize = m_iRows*m_iCols;
	int imageSize = sliceSize*m_iSlices;
	// compose the image file name, open it and check
	sprintf(tempFName, "%s%s.pic", gachPath, fName);
	ifstream inFile(tempFName, ios::binary, ios::nocreate);
	if(!inFile)
	{
		cout << "\n Could not open the image file: " 
			  << fName << endl;
		exit(1);
	}
		
	// read the header and store it
	inFile.read(m_auchPicHeader, 76);
	/////////////////////////////////////////
	// read the file into the corresponding slice
	inFile.read(&data[0][0][0], imageSize);
	inFile.close();
	*/
	Read(fName);
}

// read an image from a series of files
int C3DImage::Read(const string& achFName)
{
	char tempFName[100];
	char tempLine[100];
	char Type[100];
	char* pchPath = NULL;
	int sliceSize = m_iRows* m_iCols;
	int imageSize = sliceSize* m_iSlices;

	if (type == pic)
	{
		//	ifstream inFile(tempFName, ios::binary | ios::nocreate);
		ifstream inFile(achFName.c_str(), ios::binary);
		if (!inFile)
		{
			cout << "\n Could not open the image file: " << achFName.c_str()
				<< endl;
			exit(1);
		}

		// read the header and store it
		inFile.read(m_auchPicHeader, 76);
		m_auchPicHeader[4] = static_cast<char>(m_iSlices);

		/////////////////////////////////////////
		// read the file into the corresponding slice
		inFile.read(reinterpret_cast<char*> (&data[m_iPadding][0][0]),
			   	imageSize);

		inFile.close();
	}
	else if (type == pgm)
	{
		// read image slices from a file	

		//register int i = 0;
		register int i = m_iPadding;
		//	char achFName.c_str()Part1[128];
		//	char achTemp[128];
		//	char achFName.c_str()[256];
		char achBuffer[256];
		ifstream FileOfSlices;
		ifstream SliceFile;

		// compose the image file name, open it and read the image
		if ((pchPath = gConfig.GetStringValue("Path.Input.Image")) == NULL)
			strcpy(tempFName, gachFileOfSlices);
		else
			sprintf(tempFName, "%s%s", pchPath, gachFileOfSlices);

		FileOfSlices.open(tempFName, ios::in);

		if (!FileOfSlices)
		{
			cout << "C3DImage:;Read Fatal Error: Could Not open the file "
				<< gachFileOfSlices << endl;
			exit(0);
		}	

		// loop through the file and read each slice
		while (FileOfSlices)
		{
			FileOfSlices.getline(achBuffer, 128, '\n');
			// skip empty lines
			if (achBuffer[0] == '\0')
				continue;

			if (pchPath)
				sprintf(tempFName, "%s%s", pchPath, achBuffer);
			else
				strcpy(tempFName, achBuffer);

			SliceFile.open(tempFName, ios::binary);
			if (!SliceFile)
			{
				cout << "\n Could not open the image file: " << tempFName
					<< endl;
				exit(1);
			}

			// read the header information
			SliceFile.getline(tempLine, 90, '\n');

      //the following condition is always true...
			//if (tempLine)
			//{
				strcpy(Type, tempLine);
				if (strcmp(Type, "P5") != 0)
				{
					cout << "\n C3DImage::Read improper image formate: "
						<< tempFName << endl;
					exit(0);
				}

				while (SliceFile)
				{
					SliceFile.getline(tempLine, 90, '\n');
					if (tempLine[0] != '#')
						break;
				}

				m_iCols = atoi(strtok(tempLine, " "));
				m_iRows = atoi(strtok(NULL, " "));
				// check if the values read are the same as the global ones
				if (m_iRows != m_iRows || m_iCols != m_iCols)
				{
					cout << "\n C3DImage::Read() imporper image size: "
						<< tempFName << endl;
					exit(0);
				}

				//////////////////////////////////////
				// read the maximum pixel value
				SliceFile.getline(tempLine, 90, '\n');
				if (atoi(tempLine) != 255)
				{
					cout <<
						"\n C3DImage::Read() image max pixel value is not 255 "
						<<
						endl;
					exit(0);
				}

				/////////////////////////////////////////
				// read the file into the corresponding slice
				SliceFile.read(reinterpret_cast<char*> (&data[i][0][0]),
						  	sliceSize);
				i++;
				SliceFile.close();
			//}
		}

		FileOfSlices.close();
	}

	cout << m_iSlices - (2 * m_iPadding) << " slices read " << endl;
	return 1;
}

////////////////////////////////////////////////////////
// Method: ReadHeader
//
// Read the standard bio_rad header from a file. but replace the dimensions given 
// in the file but the actual image sizes
void C3DImage::ReadHeader(const string& fName)
{
	ifstream inFile(fName.c_str(), ios::binary);
	memset(m_auchPicHeader, 0, sizeof(char) * 77);
	inFile.read(reinterpret_cast<char*> (m_auchPicHeader), 76);

	int tempCols = m_iCols;
	int tempRows = m_iRows;
	int tempSlices = m_iSlices;
	m_auchPicHeader[0] = static_cast<char>(tempCols & 0x00FF);
	tempCols = m_iCols >> 8;
	m_auchPicHeader[1] = static_cast<char>(tempCols & 0x00FF);

	m_auchPicHeader[2] = static_cast<char>(tempRows & 0x00FF);
	tempRows = m_iRows >> 8;
	m_auchPicHeader[3] = static_cast<char>(tempRows & 0x00FF);

	m_auchPicHeader[4] = static_cast<char>(tempSlices & 0x00FF);
	tempSlices = m_iSlices >> 8;
	m_auchPicHeader[5] = static_cast<char>(tempSlices & 0x00FF);
}

// write an image to a file
int C3DImage::WriteSlices(const string& fName)
{
	string tempFName;
	int sliceSize = m_iRows* m_iCols;
	for (register int i = 0; i < m_iSlices; i++)
	{
		// formulate the file name
		stringstream ss;
		ss << i;   // write the number into the string steam
		string strNumber;
		ss >> strNumber;  // extract the string from the string stream
		ss.str(""); // clear the contents of the stringstream
		// strNumber now contains the index i

		tempFName = fName + strNumber;

		ofstream outFile(tempFName.c_str(), ios::binary);

		// write header info
		outFile << "P5\n" << "# Created by C3DImage\n" << m_iCols << " "
			<< m_iRows << "\n255\n" << endl;

		// write the image data
		outFile.write(reinterpret_cast<char*> (data[i][0]), sliceSize);

		outFile.close();
	}
	return 1;
}

// write an image to a file
// the pic header consistes of 76 bytes, the first 5 stand for
// bytes #0,#1 = x size
// bytes #2,#3 = y size
// bytes #4,#5 = Num of slices
int C3DImage::Write(const string& fName)
{
	ofstream outFile;
	outFile.open(fName.c_str(), ios::binary);

	outFile.write(m_auchPicHeader, 76);
	outFile.write(reinterpret_cast<char*> (data[0][0]),
				m_iRows * m_iCols * m_iSlices);
	outFile.close();

	return 1;
}


// Perform a connected component operation on the image
void C3DImage::ConnectedComponent()
{
	for (register int z = 1; z < m_iSlices - 1; z++)
	{
		for (register int y = giMARGIN; y < m_iRows - giMARGIN; y++)
		{
			for (register int x = giMARGIN; x < m_iCols - giMARGIN; x++)
			{
				if (data[z][y][x])
				{
					int SumOfCenterlinePixels = CountNonZeroPixels(z - 1, y, x);
					SumOfCenterlinePixels += CountNonZeroPixels(z + 1, y, x);

					data[z][y][x] = static_cast<unsigned char>(SumOfCenterlinePixels + 1);
				}
			}
		}
	}
}

int C3DImage::CountNonZeroPixels(int z, int y, int x)
{
	int result = data[z][y - 1][x - 1] +
		data[z][y - 1][x] +
		data[z][y - 1][x + 1] +
		data[z][y][x - 1] +
		data[z][y][x] +
		data[z][y][x + 1] +
		data[z][y + 1][x - 1] +
		data[z][y + 1][x] +
		data[z][y + 1][x + 1];

	return result;
}
// take the negative of the image
void C3DImage::NegateImage()
{
	unsigned char * dataPtr = &data[0][0][0];
	int imageSize = m_iRows* m_iCols* m_iSlices;
	for (register int i = 0; i < imageSize; i++)
	{
		*dataPtr = static_cast<unsigned char>(255 - *dataPtr);
		dataPtr++;
	}
}

// Threshold the image with the corresponding threshold
void C3DImage::ThresholdImage(int threshold)
{
	unsigned char * dataPtr = &data[0][0][0];
	int imageSize = m_iRows* m_iCols* m_iSlices;
	for (register int i = 0; i < imageSize; i++)
	{
		if (*dataPtr > threshold)
			*dataPtr = 255;
		else
			*dataPtr = 0;

		dataPtr++;
	}
}


void C3DImage::ThresholdImageSlices(int iThreshold)
{
	register int k, j, i;
	//	register unsigned char * dataPtr = NULL;
	for (k = 0; k < m_iSlices; k++)
	{
		for (j = 0; j < m_iRows; j++)
		{
			for (i = 0; i < m_iCols; i++)
			{
				if (data[k][j][i] < iThreshold)
					data[k][j][i] = 0;
				else
					data[k][j][i] = 255;
			}
		}
	}
}


// generate the xy projection of the 3D image onto the given image
void C3DImage::ProjectXY(CImage& a2DImage)
{
	/*
	register int iMaxValue = 0;	
	register int i, j, k;
	register int SliceSize = m_iRows*m_iCols;
	register unsigned char *puchCurrent = NULL;
	for(i = 0; i < SliceSize; i++) {
		puchCurrent = data[i][0];
		for(j = 0; j < m_iSlices; j++) {
			
			if(*puchCurrent > iMaxValue)
				iMaxValue = *puchCurrent;
			puchCurrent += SliceSize;
		}
		a2DImage.data[i/m_iCols][i%m_iCols];
	}
	*/
	register int x, y, z;
	register int maxValue = 0, tempValue = 0;

	for (y = 0; y < m_iRows; y++)
	{
		for (x = 0; x < m_iCols; x++)
		{
			maxValue = 0;
			for (z = 0; z < m_iSlices; z++)
			{
				if ((tempValue = data[z][y][x]) > maxValue)
					maxValue = tempValue;
			}

			a2DImage.data[y][x] = static_cast<unsigned char>(maxValue);
		}
	}
}

// generate the xy projection of the 3D image onto the given image
void C3DImage::MinProjectXY(CImage& a2DImage)
{
	register int tempData;
	for (register int y = 0; y < m_iRows; y++)
	{
		for (register int x = 0; x < m_iCols; x++)
		{
			register int min = 255;
			for (register int z = 0; z < m_iSlices; z++)
			{
				tempData = data[z][y][x];
				if (tempData < min)
					min = tempData;
			}
			a2DImage.data[y][x] = static_cast<unsigned char>(min);
		}
	}
}
// generate the xy projection of the 3D image onto the given image
void C3DImage::MaxProjectXY(CImage& a2DImage, int z0, int z1)
{
	z1 = (z1 == -1) ? m_iSlices : z1;
	register unsigned char * dataPtr = NULL;
	register int sliceSize = m_iRows * m_iCols;
	for (register int y = 0; y < m_iRows; y++)
	{
		for (register int x = 0; x < m_iCols; x++)
		{
			register int max = 0;
			dataPtr = &(data[z0][y][x]);
			for (register int z = z0; z < z1; z++)
			{
				if (*dataPtr > max)
					max = *dataPtr;

				dataPtr += sliceSize;
			}

			a2DImage.data[y][x] = static_cast<unsigned char>(max);
		}
	}
}
// generate the XZ projection of the 3D image onto the given image
void C3DImage::ProjectXZ(CImage& a2DImage)
{
	register int maxValue = 0;
	register int tempValue = 0;
	for (register int z = 0; z < m_iSlices; z++)
	{
		for (register int x = 0; x < m_iCols; x++)
		{
			maxValue = 0;
			for (register int y = 0; y < m_iRows; y++)
			{
				if ((tempValue = data[z][y][x]) > maxValue)
					maxValue = tempValue;
			}

			a2DImage.data[z][x] = static_cast<unsigned char>(maxValue);
		}
	}
}

// generate the XZ projection of the 3D image onto the given image
void C3DImage::MaxProjectXZ(CImage& a2DImage, int y0, int y1)
{
	y1 = (y1 == -1) ? m_iRows : y1;
	register int tempData;
	for (register int z = 0; z < m_iSlices; z++)
	{
		for (register int x = 0; x < m_iCols; x++)
		{
			register int max = 0;
			for (register int y = y0; y < y1; y++)
			{
				tempData = data[z][y][x];
				if (tempData > max)
					max = tempData;
			}
			a2DImage.data[z][x] = static_cast<unsigned char>(max);
		}
	}
}

// generate the XZ projection of the 3D image onto the given image
void C3DImage::ProjectYZ(CImage& a2DImage)
{
	register int maxValue = 0;
	register int tempValue = 0;
	for (register int z = 0; z < m_iSlices; z++)
	{
		for (register int y = 0; y < m_iRows; y++)
		{
			maxValue = 0;
			for (register int x = 0; x < m_iCols; x++)
			{
				if ((tempValue = data[z][y][x]) > maxValue)
					maxValue = tempValue;
			}

			a2DImage.data[z][y] = static_cast<unsigned char>(maxValue);
		}
	}
}

// generate the XZ projection of the 3D image onto the given image
void C3DImage::MaxProjectYZ(CImage& a2DImage, int x0, int x1)
{
	x1 = (x1 == -1) ? m_iCols : x1;
	register int tempData;
	for (register int z = 0; z < m_iSlices; z++)
	{
		for (register int y = 0; y < m_iRows; y++)
		{
			register int max = 0;
			for (register int x = x0; x < x1; x++)
			{
				tempData = data[z][y][x];
				if (tempData > max)
					max = tempData;
			}
			a2DImage.data[z][y] = static_cast<unsigned char>(max);
		}
	}
}

// write each of my slices into a separate 2DImage
void C3DImage::GenerateSlices(char* )
{
	cout <<
		"\n 'C3DImage::GenerateSlices' Empty Implementation for now. Exiting "
		<<
		endl;
	exit(0);

	/*
	assert(fName);
	char OutFName[100];
	CImage *tempImage;
	for(register int i = 0; i < slices; i++)
	{
		sprintf(OutFName, "%s%i.pgm", fName, i);
		tempImage = new CImage(data[i][0], rows, cols);
		tempImage->Write(OutFName);
		delete tempImage;
	}
	*/
}
/////////////////////////////////////////////////////////////////////////////
// Method: Compute Slice Statistics
//
void C3DImage::ComputeSliceStatistics(int SliceNum, int& median, float& stdDev)
{
	if (SliceNum < 0 || SliceNum >= m_iSlices)
	{
		cout << "\n C3DImage::ComputeSliceStatistics => Error. "
			<< " \n Invalid Slice Number .... exiting" << endl;
		exit(0);
	}

	register int pixelValue = 0;

	int Hist[256];

	memset(Hist, 0, sizeof(int) * 256);

	for (register int j = 1; j < m_iRows - 1; j++)
	{
		for (register int i = 1; i < m_iCols - 1; i ++)
		{
			pixelValue = data[SliceNum][j][i];
			Hist[pixelValue]++;
		}
	}
	// call the function that estimates the median and the stddev.
	EstimateMedianAndStdDev(Hist, median, stdDev);
}

/*
//////////////////////////////////////////////////////////////////////////////
// METHOD: Fill3D
//
// Purpose: to fill a 3D region with the given color
// Input:
//   - aVessel: the vessel portion tracked so far
//   - to     : a profile across the tracked vessel at the current tracking step
//   - minZ   : the minimum slice number with vessel information
//   - maxZ   : the largest slice number with vessel information
//   - DirFlag: Indicates whether the region is adjacent to the top or the tail of
//  			the tracked vessel.
//   - color  : filling color.
void C3DImage::Fill3D(CVessel &aVessel, CProfile *to, 
							 int minZ, int maxZ, TopEnd DirFlag, int color)
{	
	// if the vessel contains data
	if(aVessel.center)
	{
		CPoint *fromCenter = 0;
		CPoint *fromLeft   = 0;
		CPoint *fromRight  = 0;

		if(DirFlag == OnTop)
		{
			fromCenter = aVessel.center->head->data;
			fromLeft   = aVessel.left->head->data;
			fromRight  = aVessel.right->head->data;
		}
		else
		{
			fromCenter = aVessel.center->tail->data;
			fromLeft   = aVessel.left->tail->data;
			fromRight  = aVessel.right->tail->data;
		}
		
		if(to)
		{
			
			//int z = (int) ( ((float)(minZ + maxZ )) / 2.0 + 0.5); 
			for(register int z = minZ; z <= maxZ; z++)
			{
				//draw_line(fromRight,  to->right,  *this, z, color);
				draw_line(*fromCenter,   to->m_Center,   *this, z, aVessel.ID);
			//	DrawThickLine(fromRight,  to->right,  *this, z, color);
			//	DrawThickLine(fromLeft,   to->left,   *this, z, color);
			//	DrawThickLine(fromCenter, to->center, *this, z, aVessel.ID);
			}
		}
	}
	// the vessel is empty, mark the current point only
	else if(to)
	{
		CPoint *temp = &to->m_Center;

		if(temp->m_iX > 0 && temp->m_iX < m_iCols-1 &&
		   temp->m_iY > 0 && temp->m_iY < m_iRows-1)
		{
			for(register int z = minZ; z <= maxZ; z++)
			{
				data[z][temp->m_iY][temp->m_iX]   = aVessel.ID;
				data[z][temp->m_iY][temp->m_iX-1] = aVessel.ID;
				data[z][temp->m_iY][temp->m_iX+1] = aVessel.ID;
				data[z][temp->m_iY-1][temp->m_iX] = aVessel.ID;
				data[z][temp->m_iY+1][temp->m_iX] = aVessel.ID;
			}
		}
	}	
}
*/
///////////////////////////////////////////////////////
// Method: RemovePixels
//
// replace all pixels with intensitied higher than the given pixel
// value with the value.
void C3DImage::RemovePixels(int Threshold)
{
	unsigned char * dataPtr = &(data[m_iPadding][0][0]);
	int imageSize = m_iRows* m_iCols*(m_iSlices - 2 * m_iPadding);
	int PixelValue = 0;
	for (register int i = 0; i < imageSize; i++)
	{
		m_aiHistogram[*dataPtr]++;

		PixelValue = *dataPtr + Threshold;
		if (PixelValue > 255)
			PixelValue = 255;
		*dataPtr = static_cast<unsigned char>(PixelValue);

		dataPtr++;
	}
}

void C3DImage::MarkImageBoundaries()
{
	register int i, j, k;
	for (k = 0; k < m_iSlices; k++)
	{
		for (j = 0; j < m_iRows; j++)
		{
			data[k][j][giMARGIN] = IMAGE_BOUNDARY;
			data[k][j][giMARGIN - 1] = IMAGE_BOUNDARY;
		}
		for (j = 0; j < m_iRows; j++)
		{
			data[k][j][m_iCols - giMARGIN] = IMAGE_BOUNDARY;
			data[k][j][m_iCols - giMARGIN + 1] = IMAGE_BOUNDARY;
		}

		for (i = 0; i < m_iCols; i++)
		{
			data[k][giMARGIN][i] = IMAGE_BOUNDARY;
			data[k][giMARGIN - 1][i] = IMAGE_BOUNDARY;
		}
		for (i = 0; i < m_iCols; i++)
		{
			data[k][m_iRows - giMARGIN][i] = IMAGE_BOUNDARY;
			data[k][m_iRows - giMARGIN + 1][i] = IMAGE_BOUNDARY;
		}
	}
}

int compare4(const void* arg1, const void* arg2)
{
	int result = 0;
	int value1 = *((int*) arg1);
	int value2 = *((int*) arg2);

	if (value1 > value2)
		result = 1;
	else if (value1 < value2)
		result = -1;
	return result;
}
void C3DImage::ComputeStatistics()
{
	int Hist[256];
	int HistLow[256];
	int HistHigh[256];
	int Diff[256];
	float StdDevs[256];

	float StdDev;
	register int i, k;
	register unsigned char * dataPtr = NULL;
	int median = -1, minValue;
	memset(Hist, 0, sizeof(int) * 256);
	memset(HistLow, 0, sizeof(int) * 256);
	memset(HistHigh, 0, sizeof(int) * 256);
	memset(Diff, 0, sizeof(int) * 256);
	memset(StdDevs, 0, sizeof(float) * 256);

	for (k = 0; k < m_iSlices; k++)
	{
		/*
								totalCount = 0;
								memset(Hist, 0, sizeof(int)*256);
								memset(HistLow, 0, sizeof(int)*256);
								memset(HistHigh, 0, sizeof(int)*256);
								memset(Diff, 0, sizeof(int)*256);
								memset(StdDevs, 0, sizeof(float)*256);
							*/
		// compute the histogram for this image
		dataPtr = &(data[k][0][0]);
		for (i = 0; i < m_iRows*m_iCols; i++)
		{
			Hist[*dataPtr]++;
			dataPtr++;
		}
	}
	// compute the lef to right integeration of histogram
	for (i = 1; i < 256; i++)
		HistLow[i] = Hist[i - 1] + HistLow[i - 1];

	// compute the right to left integeration of histogram
	for (i = 254; i >= 0; i--)
		HistHigh[i] = Hist[i + 1] + HistHigh[i + 1];

	// compute the difference of two histograms
	for (i = 0; i < 256; i++)
		Diff[i] = abs(HistHigh[i] - HistLow[i]);

	// find the minimum value as our estimate for the median
	minValue = Diff[0] + 1;
	for (i = 0; i < 256; i++)
	{
		if (Diff[i] < minValue)
		{
			minValue = Diff[i];
			median = i;
		}
	}

	// estimate the std deviation
	memset(HistLow, 0, sizeof(int) * 256);
	memset(HistHigh, 0, sizeof(int) * 256);
	memset(Diff, 0, sizeof(int) * 256);

	StdDev = -1;

	if(median == -1)
	{
		cerr << "local variable 'median' may be used without having been initialized" << endl;
		exit(0);
	}

	for (i = 0 ; i < 256; i++)
	{
		int diff = abs(i - median);
		StdDevs[diff] += Hist[i];
	}

	for (i = 1; i < 256; i++)
		HistLow[i] = static_cast<int>(StdDevs[i - 1]) + HistLow[i - 1];

	for (i = 254; i >= 0; i--)
		HistHigh[i] = static_cast<int>(StdDevs[i + 1]) + HistHigh[i + 1];

	for (i = 0; i < 256; i++)
		Diff[i] = abs(HistHigh[i] - HistLow[i]);

	minValue = Diff[0] + 1;
	for (i = 0; i < 256; i++)
	{
		if (Diff[i] < minValue)
		{
			minValue = Diff[i];
			StdDev = static_cast<float>(i);
		}
	}
	//	if(StdDev < 0.5)
	//		StdDev = 0.5;

	for (k = 0; k < m_iSlices; k++)
	{
		// store the median value
		giMedians[k] = median;
		gfStdDevs[k] = 1.4826f * StdDev;

		cout << "Slice Num: " << k << ", M: " << giMedians[k] << ", Std: "
			<< gfStdDevs[k] << endl;
	}
}

// Crop a 3D image
C3DImage* C3DImage::Crop(int x1, int x2, int y1, int y2, int z1, int z2)
{
	if (x1 >= x2 || y1 >= y2 || z1 >= z2)
	{
		cerr << "Invalid crop points" << endl;
		exit(0);
	}
	int rows, cols, slices;

	cols = x2 - x1 + 1;
	rows = y2 - y1 + 1;
	slices = z2 - z1 + 1;

	C3DImage* result = new C3DImage(slices, rows, cols);

	int s_source, r_source, c_source;
	int s_dest, r_dest, c_dest;

	s_dest = r_dest = c_dest = 0;

	for (s_source = z1; s_source <= z2; s_source++)
	{
		s_dest = s_source - z1;
		for (r_source = y1; r_source <= y2; r_source++)
		{
			r_dest = r_source - y1;
			for (c_source = x1; c_source <= x2; c_source++)
			{
				c_dest = c_source - x1;
				result->data[s_dest][r_dest][c_dest] = this->data[s_source][r_source][c_source];
			}
		}
	}
	return result;
}


/////////////////////////////////////////////////////
// Method: FillPaddingSlices
//
// Fill the extra slices that were added before the first slice
// and after the last slice.
// we do this by grapping random squares and writing them in these slices
void C3DImage::FillPaddingSlices(int iMedian)
{
	if (m_iPadding > 0)
	{
		cout << "\tFilling Padding Slices ....";
		//	int iSlice, iRow, iCol;
		int iNumOfBytes = m_iRows* m_iCols* m_iPadding;
		register unsigned char * pData1 = &(data[0][0][0]);
		register unsigned char * pData2 = &(data[m_iSlices - m_iPadding][0][0]);
		register int i;

		// Fill the extra slices with the median intensity value
		for (i = 0; i < iNumOfBytes; i++)
		{
			*pData1++ = static_cast<unsigned char>(iMedian);
			*pData2++ = static_cast<unsigned char>(iMedian);
		}
		cout << " Padded with " << m_iPadding << " slices" << endl;
	}
}

void C3DImage::Saturate(int low_end, int high_end) 
{
	if(low_end < 0 || high_end > 255) {
		cerr << "C3DImage::Saturate:ERROR: Invalid saturation values" << endl;
		return;
	}
	for (register int s = 0; s < m_iSlices; s++)
	{
		for (register int r = 0; r < m_iRows; r++)
		{
			for (register int c = 0; c < m_iCols; c++)
			{
				if(data[s][r][c] < low_end)
					data[s][r][c] = static_cast<unsigned char>(low_end);
				else {
					if(data[s][r][c] > high_end)
						data[s][r][c] = static_cast<unsigned char>(high_end);
				}
			}
		}
	}

}

void C3DImage::CopyImageData(const C3DImage& rhs)
{
	for (register int s = 0; s < m_iSlices; s++)
	{
		for (register int r = 0; r < m_iRows; r++)
		{
			for (register int c = 0; c < m_iRows; c++)
			{
				data[s][r][c] = rhs.data[s][r][c];
			}
		}
	}
}

//////////////////////////////////////
//     Morphological Operators  	//
//////////////////////////////////////

// operations below are interfaces for functions defined in StrEle.cpp

extern C3DImage* mtk_Dilate(C3DImage* source, StrEle* S);
extern C3DImage* mtk_Erode(C3DImage* source, StrEle* S); 
extern C3DImage* mtk_Median(C3DImage* source, StrEle* S);
extern C3DImage* mtk_Open(C3DImage* source, StrEle* sErode,
	StrEle* sDilate = NULL);
extern C3DImage* mtk_Close(C3DImage* source, StrEle* sDilate,
	StrEle* sErode = NULL);

void C3DImage::Dilate(StrEle* S)
{
	C3DImage* result = mtk_Dilate(this, S);
	CopyImageData(*result);
	delete result;
}

void C3DImage::Erode(StrEle* S)
{
	C3DImage* result = mtk_Erode(this, S);
	CopyImageData(*result);
	delete result;
}

void C3DImage::Median(StrEle* S)
{
	C3DImage* result = mtk_Median(this, S);
	CopyImageData(*result);
	delete result;
}

void C3DImage::Open(StrEle* S_erode, StrEle* S_dilate)
{
	if (S_dilate == NULL)
		S_dilate = S_erode;
	C3DImage* result = mtk_Open(this, S_erode, S_dilate);
	CopyImageData(*result);
	delete result;
}

void C3DImage::Close(StrEle* S_dilate, StrEle* S_erode)
{
	if (S_erode == NULL)
		S_erode = S_dilate;
	C3DImage* result = mtk_Close(this, S_dilate, S_erode);
	CopyImageData(*result);
	delete result;
}

void C3DImage::Histogram(string file_name)
{
	int pixel_count = this->m_iCols * this->m_iRows * (this->m_iSlices - 2 * this->m_iPadding);
	ofstream out(file_name.c_str());
	out << "i\th(i)" << endl;
	for(int i = 0; i < 256; i++) {
		out << i << "\t" << static_cast<float>(m_aiHistogram[i])/static_cast<float>(pixel_count) << endl;
	}
}

//By Yousef
void C3DImage::CreateHistogram()
{
	int pix_value;
	for (register int s = 0; s < m_iSlices; s++)
	{
		for (register int r = 0; r < m_iRows; r++)
		{
			for (register int c = 0; c < m_iRows; c++)
			{
				pix_value = data[s][r][c];
				m_aiHistogram[pix_value]++;
			}
		}
	}
}
