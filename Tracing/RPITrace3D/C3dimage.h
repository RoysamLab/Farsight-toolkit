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

#ifndef C3DImage_h
#define C3DImage_h

// forward declaration
class CVessel;
class StrEle;

class C3DImage
{
public:

	enum PgmPic
	{
		pgm		= 0,
		pic		= 1
	};

	// default CTOR, do some initialization. 
	C3DImage();
	// copy CTOR, create a new image and intialize it to argument
	C3DImage(const C3DImage& image);

	// CTOR from a file, create an image and read its data from
	// the given files. This CTOR assumes the existence of a sequence
	// of files given by "FNamefrom", ..., "FNamefrom+NumOfFiles-1"
	C3DImage(const string& fName, const string& atype);

	// create an empty image of the specified size
	C3DImage(int slices, int rows, int cols);

	// DTOR
	~C3DImage();

	// allocate space for the image
	unsigned char* AllocateSpace();

	void FillImage(unsigned char* matrix, int x, int y);
	C3DImage* Crop(int x1, int x2, int y1, int y2, int z1, int z2);

	void ThresholdImage(int threshold);
	void ThresholdImageSlices(int);
	void CopyImageData(const C3DImage& rhs);

	// mark the given point in the image
	inline void MarkPoint(CPoint* aPoint, int minZ, int maxZ,
		unsigned char color = FillingColor)
	{
		unsigned char * dataPtr = &data[minZ][aPoint->m_iY][aPoint->m_iX];
		register int sliceSize = m_iRows * m_iCols;
		for (register int z = minZ; z <= maxZ; z++)
		{
			*dataPtr = color;		// [z][y][x]
			//*(dataPtr-1) = FillingColor;  // [z][y][x-1]
			//*(dataPtr+1) = FillingColor;  // [z][y][x+1]
			//*(dataPtr+cols) = FillingColor;  // [z][y+1][x]
			//*(dataPtr-cols) = FillingColor;  // [z][y-1][x]
			//*(dataPtr+cols+1) = FillingColor;  // [z][y+1][x+1]
			//*(dataPtr+cols-1) = FillingColor;  // [z][y+1][x-1]
			//*(dataPtr-cols+1) = FillingColor;  // [z][y-1][x+1]
			//*(dataPtr-cols-1) = FillingColor;  // [z][y-1][x-1]

			dataPtr += sliceSize;	// z = z+1;
		}
	}

	inline void MarkCrosshair(CPoint* aPoint, unsigned char color = 255)
	{
		data[aPoint->m_iZ][aPoint->m_iY][aPoint->m_iX] = color;
		data[aPoint->m_iZ][aPoint->m_iY][aPoint->m_iX + 1] = color;
		data[aPoint->m_iZ][aPoint->m_iY][aPoint->m_iX - 1] = color;
		data[aPoint->m_iZ][aPoint->m_iY - 1][aPoint->m_iX] = color;
		data[aPoint->m_iZ][aPoint->m_iY + 1][aPoint->m_iX] = color;
		data[aPoint->m_iZ - 1][aPoint->m_iY][aPoint->m_iX] = color;
		data[aPoint->m_iZ + 1][aPoint->m_iY][aPoint->m_iX] = color;
	}

	inline void MarkSinglePoint(CPoint* aPoint, unsigned char color = 255)
	{
		data[aPoint->m_iZ][aPoint->m_iY][aPoint->m_iX] = color;
	}

	// return 1 if the point lies withing the image boundaries, 0 otherwise
	inline int ValidPoint(CPoint& aPoint)
	{
		return (aPoint.m_iX >= 0 &&
			aPoint.m_iX < m_iCols &&
			aPoint.m_iY >= 0 &&
			aPoint.m_iY < m_iRows &&
			aPoint.m_iZ >= 0 &&
			aPoint.m_iZ < m_iSlices);
	}

	void Saturate(int low_end, int high_end); 

	// fill the area between the two given profiles. 
	//void Fill3D(CVessel &, CProfile *, int , int , TopEnd, int color = FillingColor);

	// Perform a connected component operation on the image
	void ConnectedComponent();
	int CountNonZeroPixels(int z, int y, int x);
	// read a series of images from a file. The file names are given
	// by "fNamefrom" to "fNameFrom+slices-1"
	int Read(const string& fName);
	// read the header information from a file.
	// this method is needed when we generate synthetic images for
	// testing purposes.

	void RevertToFile(const char* fName, int from);

	void ReadHeader(const string& fName);

	// write an image to afile (pic format)
	int Write(const string& fName);
	// write slices
	int WriteSlices(const string& fName);
	// take the negative of the image
	void NegateImage();

	// generate the xy projection of the 3D image onto the given image
	void ProjectXY(CImage&);
	// generate the maximum xy projection of the 3D image onto the given image
	void MinProjectXY(CImage&);
	void MaxProjectXY(CImage&, int z0=0, int z1=-1);

	// generate the XZ projection of the 3D image onto the given image
	void ProjectXZ(CImage&);
	void MaxProjectXZ(CImage&, int y0=0, int y1=-1);

	// generate the XZ projection of the 3D image onto the given image
	void ProjectYZ(CImage&);
	void MaxProjectYZ(CImage&, int x0=0, int x1=-1);

	// write each of my slices into a separate 2DImage
	void GenerateSlices(char* fName);
	// compute image statistics for the given slice
	void ComputeSliceStatistics(int SliceNum, int& median, float& stdDev);

	// replace all pixels with intensitied higher than the given pixel
	// value with the value.
	void RemovePixels(int pixelValue);
	void MarkImageBoundaries();

	void ComputeStatistics();

	void FillPaddingSlices(int iMedian);

	bool WithinImageMargin(const int &z, const int &y, const int &x, const int &margin);
	bool WithinImageMargin(const CPoint & point, const int &margin);
	bool WithinImagePadding(const int &z, const int &y, const int &x, const int &margin);
	bool WithinImagePadding(const CPoint & point, const int &margin);

	void Histogram(string file_name);

	//Yousef
	void CreateHistogram();

	// Morphological Filters
	void Dilate(StrEle* S);
	void Erode(StrEle* S);
	void Median(StrEle* S);
	void Open(StrEle* S_erode, StrEle* S_dilate);
	void Close(StrEle* S_dilate, StrEle* S_erode);

	///////////////////////////////
	// Data.

	int m_iSlices;
	int m_iRows;				// number of rows in an image
	int m_iCols;				// number of cols in an image
	int m_iPadding;
	unsigned char*** data;	// image data array 
	PgmPic type;
	char m_auchPicHeader[76];

	int m_aiHistogram[256];
};

inline bool C3DImage::WithinImageMargin(const int &z, const int &y, const int &x, const int &margin) 
{
	return (!(x < margin || x > this->m_iCols - margin || y < margin || y > this->m_iRows - margin ||
		z < margin || z > this->m_iSlices - margin) );
}

inline bool C3DImage::WithinImageMargin(const CPoint & point, const int &margin) 
{
	return WithinImageMargin(point.m_iZ, point.m_iY, point.m_iX, margin);
}

inline bool C3DImage::WithinImagePadding(const int &z, const int &y, const int &x, const int &margin) 
{
	return (!(x < margin || x > this->m_iCols - margin || y < margin || y > this->m_iRows - margin ||
		z < this->m_iPadding || z > this->m_iSlices - this->m_iPadding));
}

inline bool C3DImage::WithinImagePadding(const CPoint & point, const int &margin) 
{
	return WithinImagePadding(point.m_iZ, point.m_iY, point.m_iX, margin);
}

#endif
