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
#include <list>
#include <queue>
#include <string>
#include <cmath>
#include <stdlib.h> /* exit */
//#include "tiffio.h"
#include "itk_tiff.h"
#include "CONSTANTS.h"
#include "Mytypes.h"
#include "Config.h"
#include "Cpoints.h"
#include "Dllist.h"
#include "Cimage.h"
#include "Cvessel.h"
#include "Soma.h"

class C3DImage;
class CTemplate;
class CVector;

using namespace std;

#include "Extern.h"

// default CTOR
CImage::CImage() : m_iRows(0), m_iCols(0), data(0)
{
}

// copy CTOR
CImage::CImage(CImage& image) : m_iRows(image.m_iRows),
	m_iCols(image.m_iCols), type(image.type)
{
	if (m_iRows != 0 && m_iCols != 0)
	{
		unsigned char * buffer = AllocateSpace();
		memcpy(buffer, (image.data[0]), m_iRows * m_iCols);
	}
}

// CTOR from a file
CImage::CImage(const char* fName) : data(0)
{
	if (strstr(fName, ".pgm") || strstr(fName, ".PGM"))
		type = pgm;
	else if (strstr(fName, ".PIC") || strstr(fName, ".pic"))
		type = pic;
	else
	{
		cout << "Imporper Image format: " << fName << endl;
		exit(0);
	}

	//AllocateSpace();
	// call 'ReadImage()' to do the work
	Read(fName);
}


// CTOR create an empty image of the specified size
CImage::CImage(int y, int x) : m_iRows(y), m_iCols(x), type(pgm)
{
	AllocateSpace();
}

// CTOR from a data array
CImage::CImage(unsigned char* aData, int y, int x) : m_iRows(y), m_iCols(x),
	type(pgm)
{
	if (m_iRows != 0 && m_iCols != 0)
	{
		unsigned char * buffer = AllocateSpace();
		memcpy(buffer, aData, m_iRows * m_iCols);
	}
}
// allocate the needed space and make the data member "data" point
// to the correct locations.
unsigned char* CImage::AllocateSpace()
{
	unsigned char * buffer = 0;
	if (m_iRows != 0 && m_iCols != 0)
	{
		if ((buffer = new unsigned char [m_iRows * m_iCols + 1]) == NULL)
		{
			cout << "\n CImage::CImage => Out of Memory" << endl;
			exit(0);
		}

		memset(buffer, 0, m_iRows * m_iCols);
		buffer[m_iRows * m_iCols] = '\0';

		if ((data = new unsigned char * [m_iRows]) == NULL)
		{
			cout << "\n CImage::CTOR, Out of Memory " << endl;
			exit(1);
		}

		register int i;

		for (i = 0; i < m_iRows; i++)
			data[i] = buffer + (i * m_iCols);
	}
	return buffer;
}




// read an image from a series of files
int CImage::Read(const string& fName)
{
	// open the file and check 
	ifstream inFile(fName.c_str(), ios::binary);
	if (!inFile.is_open())
	{
		cout << "\n Could not open the image file: " << fName << endl;
		exit(1);
	}

	ProcessHeader(inFile);
	AllocateSpace();

	/////////////////////////////////////////
	// read the file into the corresponding slice
	inFile.read(reinterpret_cast<char*> (data[0]), m_iRows * m_iCols);
	inFile.close();

	return 1;
}

// write an image to a file
int CImage::Write(const string& fName)
{
	ofstream outFile(fName.c_str(), ios::binary);

	// write header info
	outFile << "P5\n" << "# Created by CImage\n" << m_iCols << " " << m_iRows
		<< "\n255\n" << endl;

	// write the image data
	outFile.write(reinterpret_cast<char*> (data[0]), m_iRows * m_iCols);
	outFile.close();

	return 1;
}


//
// This member funciton write the date of the image in a format readable by
// matlab
void CImage::WriteDataMatLabFormat(const char*)
{
	/*
	assert(fName);
	ofstream outFile(fName);
	outFile << endl;
	for(register int i = 0; i < rows; i++)
	{
		for(register int j = 0; j < cols; j++)
		{
				outFile << (int) data[i][j] << ", ";
		}
		outFile << ";" << endl;
	}
	outFile.close();
	*/
}

// make the image data that of the arguments
void CImage::FillImage(unsigned char* matrix, int x, int y)
{
	// release resources first
	if (m_iRows != 0 || m_iCols != 0)
	{
		delete [] data[0];
		delete [] data;
	}

	m_iRows = x;
	m_iCols = y;

	int size = m_iRows* m_iCols;

	unsigned char * temp;
	if ((temp = new unsigned char [size + 1]) == NULL)
		cout << "\n ERROR. CImage::CTOR \n" << endl;

	memcpy(temp, matrix, size);

	if ((data = new unsigned char * [m_iRows]) == NULL)
		cout << "\n CImage::Read(), Out of Memory " << endl;

	for (register int i = 0; i < m_iRows; i++)
		data[i] = (temp + (i * m_iCols));
}

/*
// fill the area between the two given profiles. 
// notice that the second argument might be null. In such case, we only
// fill the line of the "from" profile itself
void CImage::Fill(CProfile *from, CProfile *to, int color)
{	
	CPoint *temp = &from->m_Center;
	data[temp->m_iY][temp->m_iX] = color;
	
	//data[temp->y][temp->x-1] = color;
	//data[temp->y][temp->x+1] = color;
	data[temp->m_iY-1][temp->m_iX] = color;
	data[temp->m_iY+1][temp->m_iX] = color;
	//DrawThickLine(from->right, from->left, *this, color );
	
	//draw_line(from->right, from->left, *this, color );

	if(to)
	{
		//draw_line(from->center, to->center, *this, color);
	   
		//DrawThickLine(to->right, to->left, *this, color);
		DrawThickLine(&from->m_Center, &to->m_Center, *this, color);
		//DrawThickLine(from->right, to->right, *this, color);
		//DrawThickLine(from->left, to->left, *this, color);
		
	}
}


// Draw a vessel boundaries on the image. The boundaries are
// given as a link list of profiles.
void CImage::DrawVessel(CDLList<CProfile> *aVessel, int color)
{
	CLNode<CProfile> *temp = aVessel->head;
	//draw_line(temp->data->left, temp->data->right, *this, color);
	CLNode<CProfile> *tempNext = temp->after;

	while(temp && tempNext)
	{
		//draw_line(temp->data->left, tempNext->data->left, *this, color-1);
		//draw_line(temp->data->right, tempNext->data->right, *this, color-2);
		draw_line(temp->data->m_Center, tempNext->data->m_Center, *this, color);

		temp = tempNext;
		tempNext = tempNext->after;
	}
}

// For each vessel pixel, add 1 to my pixel value
void CImage::SumVessel(CDLList<CProfile> *aVessel, int color)
{
	CLNode<CProfile> *temp = aVessel->head;
	//draw_line(temp->data->left, temp->data->right, *this, color);
	CLNode<CProfile> *tempNext = temp->after;

	while(temp && tempNext)
	{
		//draw_line(temp->data->left, tempNext->data->left, *this, color-1);
		//draw_line(temp->data->right, tempNext->data->right, *this, color-2);
		Add_line(temp->data->m_Center, tempNext->data->m_Center, *this);

		temp = tempNext;
		tempNext = tempNext->after;
	}
}
*/
// take the negative of the image
void CImage::NegateImage()
{
	unsigned char * dataPtr = &data[0][0];
	int imageSize = m_iRows* m_iCols;
	for (register int i = 0; i < imageSize; i++)
	{
		*dataPtr = static_cast<unsigned char>(255 - *dataPtr);
		dataPtr++;

		/*
						if(*dataPtr == 0)
							*dataPtr++ = 200;
						else if(*dataPtr > 200)
							*dataPtr++ = 0;
						else
							*dataPtr++ = 200 - *dataPtr;
						*/
	}
}

// Threshold the image with the corresponding threshold
void CImage::ThresholdImage(int threshold, int ThresholdFlag)
{
	unsigned char * dataPtr = &data[0][0];
	int imageSize = m_iRows* m_iCols;
	for (register int i = 0; i < imageSize; i++)
	{
		if (! ThresholdFlag)
		{
			if (*dataPtr < threshold)
				*dataPtr = 0;
			//else
			//	*dataPtr = 255;
		}
		else
		{
			if (*dataPtr < threshold)
				*dataPtr = 255;
			else
				*dataPtr = 0;
		}

		dataPtr++;
	}
}

// Change all pixels with intensities higher than the given value to it
void CImage::RemovePixels(int pixelValue)
{
	unsigned char * dataPtr = &data[0][0];
	int imageSize = m_iRows* m_iCols;
	for (register int i = 0; i < imageSize; i++)
	{
		//if(*dataPtr < pixelValue)
		if (*dataPtr > 5 && *dataPtr < pixelValue)
			*dataPtr = 5;
		//*dataPtr = pixelValue;

		dataPtr++;
	}
}

void CImage::ProcessHeader(istream& inFile)
{
	char tempLine[90];
	char Type[10];

	if (type == pgm)
	{
		// read the header information
		inFile.getline(tempLine, 90, '\n');
    //the following condition is always true...
		//if (tempLine)
		//{
			strncpy(Type, tempLine, 3);
			if (strcmp(Type, "P5") != 0)
			{
				cout << "\n CImage::Read improper image formate: " << endl;
				exit(0);
			}
			// read the coment line
			inFile.getline(tempLine, 90);
			///////////////////////////////////////
			// read the values for rows and columns
			inFile.getline(tempLine, 90);
			m_iCols = atoi(strtok(tempLine, " "));
			m_iRows = atoi(strtok(NULL, " "));
			// check if the values read are the same as the global ones

			//////////////////////////////////////
			// read the maximum pixel value
			inFile.getline(tempLine, 90);
			if (atoi(tempLine) != 255)
			{
				cout << "\n CImage::Read() image max pixel value is not 255 "
					<< endl;
				exit(0);
			}
	//	}
	}
	else if (type == pic)
	{
		inFile.getline(reinterpret_cast<char*> (ImageHeader), 77, '\n'); // read up to 76 characters
		m_iCols = ConvertFromHex(ImageHeader[0], 0);
		m_iCols += ConvertFromHex(ImageHeader[1], 2);

		int m_iRows = ConvertFromHex(ImageHeader[2], 0);
		m_iRows += ConvertFromHex(ImageHeader[3], 2);
	}
}


/////////////////////////////////////////////////////////////////////////////
// METHOD: Histogram
//
void CImage::Histogram(int* histogram)
{
	memset(histogram, 0, sizeof(int) * 256);

	unsigned char * dataPtr = data[0];
	for (register int i = 0; i < m_iCols*m_iRows; i++)
	{
		histogram[(int) (*dataPtr)] += 1;
		dataPtr++;
	}
}

void CImage::ComputeStatistics(int& median, float& stdDev)
{
	register int pixelValue = 0;

	int Hist[256];

	memset(Hist, 0, sizeof(int) * 256);

	for (register int j = 1; j < m_iRows - 1; j++)
	{
		for (register int i = 1; i < m_iCols - 1; i ++)
		{
			pixelValue = data[j][i];
			Hist[pixelValue]++;
		}
	}

	// call the function that estimates the median and the stddev.
	EstimateMedianAndStdDev(Hist, median, stdDev);

	//	cout << "\n Median is : " << giMedian << endl;
	//	cout << "\n StdDev: " << gfStdDev << endl;
}

/////////////////////////////////////////////////////////////////////
// Method: operator[]
//
// return a reference to the indicated point
unsigned char& CImage::operator[](CPoint& rhs)
{
	if (rhs.m_iX >= 0 &&
		rhs.m_iX < m_iCols &&
		rhs.m_iY >= 0 &&
		rhs.m_iY < m_iRows)
		return data[rhs.m_iY][rhs.m_iX];

	// otherwise, always return a reference to the origin and print a message
	cout << "CImage::operator[] Invalid Point" << rhs << endl;
	return data[0][0];
}

/////////////////////////////////////////////////////////////////////
// Method: operator[]
//
// return a reference to the indicated point
unsigned char& CImage::operator[](CPoint* rhs)
{
	if (rhs->m_iX >= 0 &&
		rhs->m_iX < m_iCols &&
		rhs->m_iY >= 0 &&
		rhs->m_iY < m_iRows)
		return data[rhs->m_iY][rhs->m_iX];

	// otherwise, always return a reference to the origin and print a message
	cout << "CImage::operator[] Invalid Point" << *rhs << endl;
	return data[0][0];
}

// forward declaration to use class CSomas instantiation
class CSomas;
/////////////////////////////////////////////////////////////////////
// Method: WriteTIFF
// 	
// write an 8-bit tiff file based on a heated object color palette 
int CImage::WriteTIFF(const string& fName)
{
	// create the colormap (default colormap is heatedobject256.txt)
	uint16 rmap[256];
	uint16 gmap[256];
	uint16 bmap[256];
	register int i;

	int index = 0;
	for (i = 0; i < 256; i ++)
	{
		if (cfg_output_image_bgcolor)
		{
			rmap[i] = heated_object_colormap[index]; index++;
			gmap[i] = heated_object_colormap[index]; index++;
			bmap[i] = heated_object_colormap[index]; index++;
		}
		else
		{
			rmap[i] = gmap[i] = bmap [i] = static_cast<uint16>(i);
		}
	}

	//overwrite the heated object colormap low-end for reserved colors
	rmap[0] = 255;	gmap[0] = 192;	bmap[0] = 192;
	rmap[1] = 0;	gmap[1] = 0;	bmap[1] = 255; 
	rmap[2] = 0;	gmap[2] = 255;	bmap[2] = 0; 
	rmap[3] = 255;	gmap[3] = 0;	bmap[3] = 0; 


	if (gConfig.GetDetectSoma() && cfg_output_soma_draw_trees)
	{
		extern CSomas* gTheSomas;
		if (gTheSomas)
		{
			int iNumOfTrees = gTheSomas->m_aData[0].m_iNumOfIntersectionPoints;
			for (i = 0; i < iNumOfTrees + 7; i++)
			{
				if (i != 5)
				{
					rmap[i] = static_cast<uint16>(abs(255 - (i * 4 + 1)));
					gmap[i] = static_cast<uint16>(abs(0 + (i * 32 + 1)));
					bmap[i] = static_cast<uint16>(abs(0 + (i * 8 + 1)));
				}
			}
		}
	}

	// scale up from 0-255(8-bit color) to 0-65535 (16-bit color)
	for (i = 0; i < 256; i ++)
	{
		rmap[i] = static_cast<uint16>(static_cast<float>(rmap[i]) / 255.0f * 65535.0f);
		gmap[i] = static_cast<uint16>(static_cast<float>(gmap[i]) / 255.0f * 65535.0f);
		bmap[i] = static_cast<uint16>(static_cast<float>(bmap[i]) / 255.0f * 65535.0f);
	}

	TIFF* out = TIFFOpen(fName.c_str(), "w");
	uint32 width = m_iCols;
	uint32 height = m_iRows;
	int sampleperpixel = 1;    // or 3 if it is RGB, 4 if it is RGBA 

	//	Now we need to set the tags in the new image file, and the essential ones are the following: 
	TIFFSetField(out, TIFFTAG_IMAGEWIDTH, width);					// set the width of the image
	TIFFSetField(out, TIFFTAG_IMAGELENGTH, height);					// set the height of the image
	TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, sampleperpixel);		// set number of channels per pixel
	TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 8);					// set the size of the channels
	TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);	// set the origin of the image.
	TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);	// store data contiguously
	TIFFSetField(out, TIFFTAG_PHOTOMETRIC, 3);						// use the heated object palette

	// We set the strip size of the file to be size of one row of pixels
	TIFFSetField(out,
		TIFFTAG_ROWSPERSTRIP,
		TIFFDefaultStripSize(out, width * sampleperpixel));
	TIFFSetField(out, TIFFTAG_COLORMAP, rmap, gmap, bmap);

	//set the resolution to 72dpi
	TIFFSetField(out, TIFFTAG_XRESOLUTION, (float) 72.0);
	TIFFSetField(out, TIFFTAG_YRESOLUTION, (float) 72.0);
	//Now writing image to the file one strip (row) at a time
	for (uint32 row = 0; row < height; row++)
	{
		//if write fails, we return
		if (TIFFWriteScanline(out, data[row], row, 0) < 0)
		{
			return -1;
		}
	}

	// close the output file 
	(void) TIFFClose(out); 

	return 0;
}

// DTOR
CImage::~CImage()
{
	if (m_iRows > 0 && m_iCols > 0)
	{
		// delete the image's buffer
		delete [] data[0];
		// and its pointers
		delete [] data;
		//		int i;
		//		for(i = 0; i < m_iRows; i++) {
		//			delete [] data[i];
		//		}
		//		delete [] data;
	}
}
