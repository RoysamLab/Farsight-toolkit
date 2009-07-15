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
//#pragma warning(disable:4702)  // unreachable STLport code

#include <iostream>
#include <fstream>
#include <list>
#include <queue>
#include <cmath>
#include <ctime>

#include "Mytypes.h"
#include "Config.h"
#include "Cpoints.h"
#include "CONSTANTS.h"
#include "Dllist.h"
#include "Cimage.h"
#include "C3dimage.h"
#include "Cvessel.h"
#include "Cvector.h"
#include "Template.h"
#include "Ctree.h"
#include "Soma.h"
#include "Letters.h"
#include "Extern.h"
#include "StrEle.h"

extern int iRows;
extern int iCols;

EquivalentColors* ColorsArray = NULL;

extern int giNumOfSomas;
CSomas* gTheSomas = 0;

EquivalentColors::EquivalentColors()
{
	label = 0;
	NumOfColors = 0;
	validColor = 0;
	memset(myColors, 0, sizeof(int) * MaxNumOfColors);
}

// Merge two color labels into one. Use the smaller label of the two
void EquivalentColors::Merge(int color, int recursiveFlag)
{
	if (color < label)
		label = color;

	validColor = 1;

	int flag = 0;
	for (register int i = 0; i < NumOfColors; i++)
	{
		if (myColors[i] == color)
		{
			flag = 1;
			break;
		}
	}
	if (flag)
		return;

	myColors[NumOfColors] = color;
	NumOfColors++;

	if (recursiveFlag)
	{
		for (register int i = 0; i < NumOfColors - 1; i++)
		{
			ColorsArray[myColors[i]].Merge(color, 0);
			ColorsArray[color].Merge(myColors[i], 0);
		}
	}
}


//////////////////////////////////
//
CSoma::CSoma() : m_Center(0, 0, 0, 0, 0)
{
	m_aTrees = NULL;
	m_iLabel = 0;
	m_iVolume = 0;
	m_iNumOfIntersectionPoints = 0;
	m_iNumOfSegments = 0;
	m_iNumOfSomas = 0;
	m_iInitialized = 0;
	m_iHitsImageBoundary = 0;
	m_iSumOfAllSegments = 0;
	m_achID[0] = '\0';

	memset(m_aiConnectedSomaIDs, 0, sizeof(int) * MaxNumOfSomas);
	memset(m_aiSegments, 0, sizeof(int) * MaxNumOfSegments);
	memset(m_aiIntersectionPoints, 0, sizeof(int) * MaxNumOfSegments);
}

///////////////////////////////////
// Method: AddIntersectionPoint
//
// Add the intersectionPoint with the given ID to my list
void CSoma::AddIntersectionPoint(int id)
{
	if (m_iNumOfIntersectionPoints >= MaxNumOfSegments - 1)
		return;

	for (register int i = 0; i < m_iNumOfIntersectionPoints; i++)
	{
		if (m_aiIntersectionPoints[i] == id)
			return;
	}

	m_aiIntersectionPoints[m_iNumOfIntersectionPoints++] = id;
}	

///////////////////////////////////
// Method: AddSegment
// 
// Add the given segment to my list of segments. 
// retrun 1 upon success, 0 otherwise
int CSoma::AddSegment(int id)
{
	if (m_iNumOfSegments >= MaxNumOfSegments - 1)
		return 0;

	for (register int i = 0; i < m_iNumOfSegments; i++)
		if (m_aiSegments[i] == id)
			return 0;

	m_aiSegments[m_iNumOfSegments++] = id;
	return 1;
}

///////////////////////////////////
// Method: AddConnectedSomaID
//
// Add the given soma id to my list of somas I am connected to
void CSoma::AddConnectedSomaID(int id)
{
	if (m_iNumOfSomas >= MaxNumOfSomas - 1)
		return;

	for (register int i = 0; i < m_iNumOfSomas; i++)
	{
		if (m_aiConnectedSomaIDs[i] == id)
			return;
	}

	m_aiConnectedSomaIDs[m_iNumOfSomas++] = id;
}
void CSoma::WriteID(CImage& anImage, unsigned char color)
{
	int VFrom = m_Center.m_iY - DigitWidth / 2;
	int VTo = m_Center.m_iY + DigitWidth / 2;
	int HFrom = m_Center.m_iX - DigitWidth / 2;
	int HTo = m_Center.m_iX + DigitWidth / 2;

	WriteLabel(anImage, VFrom, VTo, HFrom, HTo, color);
}

void CSoma::WriteIDXZ(CImage& anImage, unsigned char color)
{
	int VFrom = m_Center.m_iZ - DigitWidth / 2;
	int VTo = m_Center.m_iZ + DigitWidth / 2;
	int HFrom = m_Center.m_iX - DigitWidth / 2;
	int HTo = m_Center.m_iX + DigitWidth / 2;

	WriteLabel(anImage, VFrom, VTo, HFrom, HTo, color);
}

void CSoma::WriteIDYZ(CImage& anImage, unsigned char color)
{
	int VFrom = m_Center.m_iZ - DigitWidth / 2;
	int VTo = m_Center.m_iZ + DigitWidth / 2;
	int HFrom = m_Center.m_iY - DigitWidth / 2;
	int HTo = m_Center.m_iY + DigitWidth / 2;

	WriteLabel(anImage, VFrom, VTo, HFrom, HTo, color);
}

void CSoma::WriteLabel(CImage& anImage, int VFrom, int VTo, int HFrom,
	int HTo, unsigned char color)
{
	int iCols = anImage.m_iCols;
	int iRows = anImage.m_iRows;
	// to avoid crowding we only write labels to first
	// 26 somas
	if (m_iLabel > 26)
		return;

	int Vindex = 0;
	int Hindex = 0;
	for (register int i = VFrom; i <= VTo; i++)
	{
		Hindex = 0;
		for (register int j = HFrom; j <= HTo; j++)
		{
			if (j >= 0 && j < iCols && i >= 0 && i < iRows)
			{
				if (Letters[m_iLabel - 1][Vindex][Hindex])
					anImage.data[i][j] = color;
			}
			Hindex++;
		}
		Vindex++;
	}
}
// set the label and name of the soma
void CSoma::SetLabel(int i)
{
	m_iLabel = i;
	if (i <= 26)   // 26 the number of alphabet
	{
		m_achID[0] = (char) ('A' + (i - 1));
		m_achID[1] = '\0';
	}
	else
	{
		int index1 = m_iLabel / 26;   
		int index2 = m_iLabel % 26;
		m_achID[0] = (char) ('A' + index1);
		m_achID[1] = (char) ('A' + index2);
		m_achID[2] = '\0';
	}
}

/////////////////////////////////////////
// Mehtod: ConstructTrees
//
// Each of my segments (dendrite/axons) is thought as starting a tree of 
// dendrites. This method constructs such dendrites
void CSoma::ConstructTrees()
{
	// allocate space
	if (m_iNumOfIntersectionPoints == 0)
		return;
	m_aTrees = new CTree[m_iNumOfIntersectionPoints];

	CIntersectionPoint* pIntPoint = NULL;
	register int index = 0;
	for (register int i = 0; i < m_iNumOfIntersectionPoints; i++)
	{
		pIntPoint = gIntersectionPoints.GetPoint(m_aiIntersectionPoints[i]);
		m_aTrees[index].ConstructTree(m_iLabel, pIntPoint);
		if (m_aTrees[index].root != NULL)
		{
			index++;
			m_iSumOfAllSegments += m_aTrees[i].m_iSumOfAllBranches;
		}
	}
	m_iNumOfIntersectionPoints = index;
}

////////////////////////
// Method: Print
//
// Print the soma information
void CSoma::Print(ostream& outFile, int id)
{	
	outFile << "\n Soma " << id+1 << ":\n";
	outFile << "============\n";
	outFile << "\tCenter: " << gTheSomas->m_aData[id].m_Center.m_iX << ", " << gTheSomas->m_aData[id].m_Center.m_iY << ", " << gTheSomas->m_aData[id].m_Center.m_iZ << "\n";
	outFile << "\tVolume:   " << gTheSomas->m_aData[id].m_iVolume << "\n";

	/*if (m_iNumOfSomas)
	{
		if (m_iNumOfSomas == 1)
			outFile << "\t Connected to " << m_iNumOfSomas << " Soma: ";
		else
			outFile << "\t Connected to " << m_iNumOfSomas << " Somas: ";

		for (register int i = 0; i < m_iNumOfSomas; i++)
		{
			if (i == m_iNumOfSomas - 1)
				outFile << gTheSomas->m_aData[m_aiConnectedSomaIDs[i] - 1].m_achID;
			else
				outFile << gTheSomas->m_aData[m_aiConnectedSomaIDs[i] - 1].m_achID << ", ";
		}
	}
	else
		outFile << "\tIsolated From Other Somas.\n";

	if (m_iHitsImageBoundary)
		outFile << "\tHits Image Boundary? Yes\n";
	else
		outFile << "\tHits Image Boundary? No\n";

	outFile << "\tSum Of All Branches: " << m_iSumOfAllSegments << "\n";
	outFile << "\tTotal Of " << m_iNumOfIntersectionPoints << " Trees:\n";

	int number;
	for (register int i = 0; i < m_iNumOfIntersectionPoints; i++)
	{
		number = i+1;
		outFile << "\n\tTREE " << number; 
		m_aTrees[i].Print(outFile);
	}*/
}

////////////////////////
// Method: Print
//
// Similar to the previous one, but print tree information as well
void CSoma::PrintTrees(ostream& outFile)
{
	outFile << "\n Soma " << m_achID << ": \n";
	outFile << "\t center: " << m_Center.m_iX << ", " << m_Center.m_iY << ", " << m_Center.m_iZ << "\n";
	outFile << "\t Volume: " << m_iVolume << "\n";

	outFile << "The Trees: \n";
	for (register int i = 0; i < m_iNumOfIntersectionPoints; i++)
	{
		m_aTrees[i].Print(outFile);
	}
}


////////////////////////
// Method: DrawTrees
//
// Draw the trees hanging from this soma
void CSoma::DrawTrees()
{
	for (register int i = 0; i < m_iNumOfIntersectionPoints; i++)
		m_aTrees[i].Draw(static_cast<unsigned char>(m_iLabel + 5 + i));
}


/////////////////////////////////
// Method: PrintNeurolucidaFormat
//
void CSoma::PrintNeurolucidaFormat(ofstream& outFile)
{
	register int iIndex;
	register int i, j, k;
	int X, Y, Z;
	int iColor;

	// print the soma body
	outFile << "(\"Soma " << m_achID << "\"\n (Color DarkGreen)\n (CellBody)"
		<< endl;

	for (iIndex = 0; iIndex < gTheSomas->m_aData[0].m_iVolume; iIndex++)
	{
		X = gaSomaPoints[iIndex].m_iX;
		Y = gaSomaPoints[iIndex].m_iY;
		Z = gaSomaPoints[iIndex].m_iZ;

		for (i = -1; i <= 1; i++)
		{
			for (j = -1; j <= 1; j++)
			{
				for (k = -1; k <= 1; k++)
				{
					if (The3DImage->data[Z + k][Y + j][X + i] != SomaColor)
					{
						outFile << "( " << X << " " << Y << " " << Z << " )"
							<< endl;
						i = 2;
						j = 2;
						break;
					}
				}
			}
		}
	}
	outFile << "); End Cell Body" << endl;

	// output the trees
	for (i = 0; i < m_iNumOfIntersectionPoints; i++)
	{
		iColor = m_iLabel + 5 + i;
		if (iColor >= MaxNumOfColors)
			iColor = MaxNumOfColors - 1; 
		m_aTrees[i].PrintNeurolucidaFormat(outFile, iColor);
	}
}

/////////////////////////////////////////////////////////////////////////////
//									Class CSomas											   //
/////////////////////////////////////////////////////////////////////////////

// CTOR, default
CSomas::CSomas(int iCount) : m_iNumOfSomas(0), m_aData(0)
{
	if (iCount <= 0)
		return;

	m_iNumOfSomas = iCount;
	m_aData = new CSoma[m_iNumOfSomas];

	register int i;
	//register unsigned char *dataPtr = &SomaImage->data[0][0];

	for (i = 0; i < m_iNumOfSomas; i++)
	{
		m_aData[i].SetLabel(i + 1);
	}
	/*
	// calculate the center and area of each of the somas
	for(i = 0; i < iRows; i++)
	{
		for(j = 0; j < iCols; j++)
		{
			// the pixel belongs to a soma
			if( (pixelValue = *dataPtr) != 0)
			{
				m_aData[pixelValue -1].m_iCenterX += j;
				m_aData[pixelValue -1].m_iCenterY += i;
				m_aData[pixelValue - 1].m_iArea++;
			}
			dataPtr++;
		}
	}
	// for each soma, its center is given by
	for(i = 0; i < m_iNumOfSomas; i++)
	{
		if(m_aData[i].m_iArea)
		{
			m_aData[i].m_iCenterX /= m_aData[i].m_iArea;
			m_aData[i].m_iCenterY /= m_aData[i].m_iArea;
		}
	}
	// There is a bug in the code that determines the number of somas
	// (see ConnectedComponent). Sometimes it will overestimates the number
	// of somas available in an image. Such extra somas have no area. So remove
	// all such somas
	*/
}

void CSomas::WriteIDs(CImage& anImage, unsigned char color)
{
	m_iNumOfSomas = giNumOfSomas;
	/*for (register int i = 0; i < m_iNumOfSomas; i++)
	{
		if (i < 16)
		{
			m_aData[i].WriteID(anImage, color);
		}
	}*/
}

void CSomas::WriteIDsXZ(CImage& anImage, unsigned char color)
{
	m_iNumOfSomas = giNumOfSomas;
	for (register int i = 0; i < m_iNumOfSomas; i++)
	{
		if (i < 16)
		{
			m_aData[i].WriteIDXZ(anImage, color);
		}
	}
}
void CSomas::WriteIDsYZ(CImage& anImage, unsigned char color)
{
	m_iNumOfSomas = giNumOfSomas;
	for (register int i = 0; i < m_iNumOfSomas; i++)
	{
		if (i < 16)
		{
			m_aData[i].WriteIDYZ(anImage, color);
		}
	}
}

/////////////////////////////////////
// Mehtod: ConstructTrees
//
//
// For each of the somas, Construct all the trees rooted at each of its 
// dendrites/axon. 
void CSomas::ConstructTrees()
{
	for (register int i = 0; i < m_iNumOfSomas; i++)
		m_aData[i].ConstructTrees();
}

/////////////////////////////////////////
// Method: ConnectSomas
//
// inform the two given somas that they are connected
void CSomas::ConnectSomas(int id1, int id2)
{
	m_aData[id1 - 1].AddConnectedSomaID(id2);
	m_aData[id2 - 1].AddConnectedSomaID(id1);
}
/////////////////////////////
// Method: Print
// 
// Print the somas
void CSomas::Print(char* fName)
{
	std::string image_name = gConfig.GetImageName();
	if (fName)
	{
		ofstream outFile(fName);

		if (m_iNumOfSomas == 0)
			outFile << "\n There is no somas in the image " << image_name.c_str() << std::endl;
		else if (m_iNumOfSomas == 1)
			outFile << "\n There is one soma in the image " << image_name.c_str() << std::endl;
		else
			outFile << "\n \n There are a total of " << m_iNumOfSomas << ", in the image " << image_name.c_str() << std::endl;

		for (register int i = 0; i < m_iNumOfSomas; i++)
			m_aData[i].Print(outFile, i);

		/*
						fprintf(outFile,"\n==============> For Ali <===============\n");
						for(i = 0; i < gIntersectionPoints.m_iNumOfElements; i++)
						{
							CPoint *pPoint = &gIntersectionPoints.m_apData[i]->m_Point;
							fprintf(outFile, "%3i %3i %3i\n", pPoint->m_iX, pPoint->m_iY, pPoint->m_iZ);
						}
						fprintf(outFile, "-1 -1 -1\n");
						*/

		outFile.close();
	}
}

/////////////////////////////
// Method: Print
// 
// Print the somas
void CSomas::PrintTrees(char* fName)
{
	if (fName)
	{
		ofstream outFile(fName);
		for (register int i = 0; i < m_iNumOfSomas; i++)
			m_aData[i].PrintTrees(outFile);
		outFile.close();
	}
}

/////////////////////////////
// Method: DrawTrees
// 
// Draw the trees hanging from all the somas
void CSomas::DrawTrees()
{
	for (register int i = 0; i < m_iNumOfSomas; i++)
		m_aData[i].DrawTrees();
}

/////////////////////////////////
// Method: SomaHitsImageBoundary
//
// inform the some with the given id that it hits the image boundary
void CSomas::SomaHitsImageBoundary(int somaID)
{
	if (somaID > 0 && somaID <= m_iNumOfSomas)
		m_aData[somaID - 1].m_iHitsImageBoundary = 1;
}

/////////////////////////////////
// Method: AddSegment
//
// Add the following segment to the soma's list of ids
int CSomas::AddSegment(int somaID, int segmentID)
{
	if (somaID > 0 && somaID <= m_iNumOfSomas)
		return m_aData[somaID - 1].AddSegment(segmentID);

	return 0;
}

/////////////////////////////////
// Method: PrintNeurolucidaFormat
//
void CSomas::PrintNeurolucidaFormat(char* fName)
{
	ofstream outFile;
	outFile.open(fName);

	for (register int i = 0; i < m_iNumOfSomas; i++)
	{
		outFile << "\n; New Soma " << endl;

		m_aData[i].PrintNeurolucidaFormat(outFile);
	}
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


int compare(const void* arg1, const void* arg2)
{
	int v1 = (* (int*) arg1);
	int v2 = (* (int*) arg2);

	if (v1 > v2)
		return 1;
	else if (v1 < v2)
		return -1;
	else
		return 0;
}
inline void MergeTwoColors(int color1, int color2)
{
	ColorsArray[color1].Merge(color2, 1);
	ColorsArray[color2].Merge(color1, 1);
}

// perform connected component on the given image. The function returns the
// number of such components found. The argument NumOfColors denotes the number
// of different colors available in the image
int ConnectedComponent(CImage* anImage, int NumOfColors)
{
	int iResult = 0;
	register int i;
	register int j;
	register int x;
	register int y;

	// create and initialize the colors array
	ColorsArray = new EquivalentColors[NumOfColors + 1];
	for (i = 1; i <= NumOfColors; i++)
		ColorsArray[i].SetLabel(i);

	unsigned char * *Data = anImage->data;

	int iRows = anImage->m_iRows;
	int iCols = anImage->m_iCols;
	int pixelValue;
	int pixelValue2;
	for (i = 2; i < iRows - 2; i++)
	{
		for (j = 2; j < iCols - 2; j++)
		{
			pixelValue = Data[i][j];
			// if any of the neighbors is set to a different color, 
			// unify the colors
			if (pixelValue)
			{
				ColorsArray[pixelValue].validColor = 1;

				for (x = -1; x <= 1; x++)
				{
					for (y = -1; y <= 1; y++)
					{
						if (x == 0 && y == 0)
							continue;

						pixelValue2 = Data[i + x][j + y];
						if (pixelValue2 && pixelValue2 != pixelValue)
						{
							MergeTwoColors(pixelValue, pixelValue2);
						}
					} // for y
				} // for x
			} // if
		} // for j
	} // for i


	// make color labels consequtive
	int LabelsArray[256];
	memset(LabelsArray, 0, sizeof(int) * 256);
	int Index = 1;
	int NumOfLabels = 0;
	for (i = 1; i <= NumOfColors; i++)
	{
		if (! ColorsArray[i].validColor)
			continue;

		ColorsArray[i].label += 1000;
		LabelsArray[Index++] = ColorsArray[i].label;

		NumOfLabels++;
	}
	qsort(LabelsArray, NumOfLabels, sizeof(int), compare);

	int label = 1;
	// for the first label case
	for (i = 1; i < NumOfColors; i++)
	{
		if (ColorsArray[i].label == LabelsArray[1])
			ColorsArray[i].label = label;
	}
	label++;
	for (i = 2; i < NumOfColors; i++)
	{
		if (LabelsArray[i] != LabelsArray[i - 1])
		{
			for (j = 1; j <= NumOfColors; j++)
			{
				if (ColorsArray[j].label == LabelsArray[i])
					ColorsArray[j].label = label;
			}	
			label++;
		}
	}

	// This is how many connected-components (somas) we have in the image
	iResult = label - 1;

	int labelIndex = 0;
	// now do the coloring
	for (i = 0; i < iRows; i++)
	{
		for (j = 0; j < iCols; j++)
		{
			if ((pixelValue = Data[i][j]) != 0)
			{
				labelIndex = ColorsArray[pixelValue].label;
				Data[i][j] = static_cast<unsigned char>(labelIndex);
			}
		}
	}

	return iResult;
}


////////////////////////////////////////////////////////////////////////////
// Function: LocateSomas
//
// Find the location of the Somas
void LocateSomas2()
{
	int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;
	CImage tempImage(iRows, iCols);

	DetectSomas(ProjectionImage->data, tempImage.data);

	// use the histogram of seed points to estimate the median and the stddev of 
	// the foreground pixels

	int Threshold = static_cast<int>(giForegroundMedian + 2.0 * gfForegroundStdDev);

	ProjectionImage->ThresholdImage(Threshold, 0);
	tempImage.ThresholdImage(Threshold, 0);

	//tempImage.Write("tempImage.pgm");

	///////////////////////////////////////////////////////
	// NOTE: All what follows does not apply for images with multiple somas
	// change it in a manner similar to the 2dTrack programs to hadle multiple
	// somas
	// 3-4-1999

	int* PlaneSum = new int[iSlices];
	memset(PlaneSum, 0, sizeof(int) * iSlices);

	register int i, j, k;
	int minX, minY, maxX, maxY;
	minX = minY = 9999;
	maxX = maxY = 0;

	for (k = 0; k < iSlices; k++)
	{
		for (j = giMARGIN; j < iRows - giMARGIN; j++)
		{
			for (i = giMARGIN; i < iCols - giMARGIN; i++)
			{
				// if the point belongs to the soma, add its value
				// to the soma pixel sum in this image plane
				if (tempImage.data[j][i])
				{
					PlaneSum[k] += The3DImage->data[k][j][i];

					// figure out the boundaries of the region
					if (i < minX)
						minX = i;
					if (i > maxX)
						maxX = i;
					if (j < minY)
						minY = j;
					if (j > maxY)
						maxY = j;
				}
			}
		}
	}

	// find the plane with the higest sum
	int max = 0;
	int PlaneIndex = 0;
	for (i = 0; i < iSlices; i++)
	{
		if (PlaneSum[i] > max)
		{
			max = PlaneSum[i];
			PlaneIndex = i;
		}
	}
	int FromPlane = PlaneIndex;
	int ToPlane = PlaneIndex;
	int threshold = static_cast<int>(0.2 * max);
	while (PlaneSum[FromPlane] >= threshold)
	{
		if (FromPlane == 0)
			break;
		FromPlane -= 1;
	}
	while (PlaneSum[ToPlane] >= threshold)
	{
		if (ToPlane == iSlices - 1)
			break;
		ToPlane += 1;
	}

	// for now, we fill the region between FromPlane and ToPlane as the single
	// one soma in the image.	
	int SomaVolume = 0;
	int pixelValue = 0;
	int SomaID = 1;
	for (j = minY; j <= maxY; j++)
	{
		for (i = minX; i <= maxX; i++)
		{
			// if the point belong to the soma, mark it in the 3D image and 
			// update soma volume
			if (tempImage.data[j][i] != 0)
			{
				TrackImageXY->data[j][i] = SomaColor;
				for (k = FromPlane; k <= ToPlane; k++)
				{
					pixelValue = The3DImage->data[k][j][i];
					The3DImage->data[k][j][i] = SomaColor;
					Empty3DImage->data[k][j][i] = static_cast<unsigned char>(SomaID);

					if (pixelValue >= threshold)
					{
						TrackImageXZ->data[k][i] = SomaColor;
						TrackImageYZ->data[k][j] = SomaColor;

						SomaVolume++;
					}
				}
			}
		}
	}

	//TrackImageXY->Write("somaImage.pgm");

	//exit(0);
	CPoint centerPoint;
	int x = static_cast<int>(static_cast<float>(minX + maxX) / 2.0 + 0.5);
	int y = static_cast<int>(static_cast<float>(minY + maxX) / 2.0 + 0.5);
	int z = PlaneIndex;
	centerPoint.m_iX = x;
	centerPoint.m_iY = y;
	centerPoint.m_iZ = z;

	double radius = (maxX - minX) * (maxX - minX) +
		(maxY - minY) * (maxY - minY);
	radius = static_cast<int>(std::sqrt( radius) );
	if ((ToPlane - FromPlane) > radius)
		radius = ToPlane - FromPlane;

	radius = radius / 2;

	delete [] PlaneSum;

	// create the soma collection
	giNumOfSomas = 1;
	gTheSomas = new CSomas(giNumOfSomas);
	gTheSomas->m_aData[0].m_iVolume = SomaVolume;
	gTheSomas->m_aData[0].m_Center = centerPoint;
}


////////////////////////////////////////////////////////////////////////////
// Function: LocateSomas
//
// Find the location of the Somas
void LocateSomas3()
{
	int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;
	register int i,j,k;

	CImage tempImage(iRows, iCols);
	CImage tempImageXZ(iSlices, iCols);
	CImage tempImageYZ(iSlices, iRows);

	int iWidth = (int) (gfWidthSum / giNumOfWidthSumMembers + 0.5);

	int StructElemSize = iWidth + 7; // it was + 2 11-5-99

	CPoint* aDiskPoints = new CPoint[4 * StructElemSize* StructElemSize];
	int iNumOfPoints = 0;

	ConstructStructElem(StructElemSize, aDiskPoints, iNumOfPoints);

	// DetectSomas(ProjectionImage, &tempImage, StructElemSize);
	DetectSomas(ProjectionImage, & tempImage, aDiskPoints, iNumOfPoints);	

	// estimate the threshold from the brightest region in the image
	int Threshold = static_cast<int>(giForegroundMedian + gfForegroundStdDev);
	int aiHistogram[256];
	tempImage.Histogram(aiHistogram);

	for (i = 255; i > giForegroundMedian; i--)
	{
		if (aiHistogram[i] >= iNumOfPoints)
		{
			Threshold = static_cast<int>(0.6 * i);
			break;
		}
	}

	// use the histogram of seed points to estimate the median and the stddev of 
	// the foreground pixels

	if (Threshold > 250)
		Threshold = 250;


	ProjectionImage->ThresholdImage(Threshold, 0);

	tempImage.ThresholdImage(Threshold, 0);

	CImage tempImageXY(tempImage);

	// The soma is defined as the region of intersection betweent the three images

	int SomaVolume = 0;
	int minX = iCols + 1, minY = iRows + 1, minZ = iSlices + 1;
	int maxX = - 1;
	int maxY = - 1;
	int maxZ = - 1;

	int* aiPlaneSum = new int[iSlices];
	
	int iDenom = iRows;
	if (iCols < iRows)
		iDenom = iCols;
	if (2 * iSlices < iDenom)
		StructElemSize = static_cast<int>(StructElemSize * (2.0 * iSlices / iDenom));

	if (StructElemSize < iWidth)
		StructElemSize = iWidth + 1;

	ConstructStructElem(StructElemSize, aDiskPoints, iNumOfPoints);


	DetectSomas(TrackImageYZ, & tempImageYZ, aDiskPoints, iNumOfPoints);	
	tempImageYZ.Histogram(aiHistogram);
	for (i = 255; i > giForegroundMedian; i--)
	{
		if (aiHistogram[i] >= iNumOfPoints)
		{
			Threshold = static_cast<int>(0.6 * i);
			break;
		}
	}
	tempImageYZ.ThresholdImage(Threshold, 0);

	DetectSomas(TrackImageXZ, & tempImageXZ, aDiskPoints, iNumOfPoints);	
	tempImageXZ.Histogram(aiHistogram);
	for (i = 255; i > giForegroundMedian; i--)
	{
		if (aiHistogram[i] >= iNumOfPoints)
		{
			Threshold = static_cast<int>(0.6 * i);
			break;
		}
	}			

	delete [] aDiskPoints;
	delete [] aiPlaneSum;

	
	// count the number of pixels in the soma
	giSomaVolume = 0;
	for (k = 0; k < iSlices; k++)
	{
		for (i = 0; i < iRows; i++)
		{
			for (j = 0; j < iCols; j++)
			{
				if (tempImageXY.data[i][j] &&
					tempImageXZ.data[k][j] &&
					tempImageYZ.data[k][i])
				{
					giSomaVolume++;
				}
			}
		}
	}

	gaSomaPoints = new CPoint[giSomaVolume];

	register int index = 0;
	for (k = 0; k < iSlices; k++)
	{
		for (i = 0; i < iRows; i++)
		{
			for (j = 0; j < iCols; j++)
			{
				if (tempImageXY.data[i][j] &&
					tempImageXZ.data[k][j] &&
					tempImageYZ.data[k][i])
				{
					TrackImageXY->data[i][j] = SomaColor;
					TrackImageXZ->data[k][j] = SomaColor;
					TrackImageYZ->data[k][i] = SomaColor;
					The3DImage->data[k][i][j] = SomaColor;

					gaSomaPoints[index].m_iX = j;
					gaSomaPoints[index].m_iY = i;
					gaSomaPoints[index].m_iZ = k;
					index++;
					//						Empty3DImage->data[k][i][j] = StartingSomaColor+1;


					SomaVolume++;
					if (k < minZ)
						minZ = k;
					if (k > maxZ)
						maxZ = k;
					if (i < minY)
						minY = i;
					if (i > maxY)
						maxY = i;
					if (j < minX)
						minX = j;
					if (j > maxX)
						maxX = j;
				}
			}
		}
	}

	// temp stuff to write all soma points that are exterior to a file for 
	// vizualization with ASAD's program
	//	register int ii, jj, kk;
	//	register int iFlag;

	CPoint centerPoint;
	centerPoint.m_iX = static_cast<int>(static_cast<float>(minX + maxX) / 2.0 + 0.5);
	centerPoint.m_iY = static_cast<int>(static_cast<float>(minY + maxY) / 2.0 + 0.5);
	centerPoint.m_iZ = static_cast<int>(static_cast<float>(minZ + maxZ) / 2.0 + 0.5); 

	double radius = (maxX - minX) * (maxX - minX) +
		(maxY - minY) * (maxY - minY);
	radius = static_cast<int>(sqrt(radius) / 2.0 + 0.5);

	// create the soma collection
	giNumOfSomas = 1;
	gTheSomas = new CSomas(giNumOfSomas);
	gTheSomas->m_aData[0].m_iVolume = SomaVolume;
	gTheSomas->m_aData[0].m_Center = centerPoint;
}



//By Yousef
//Try this
void LocateSomas3_v2()
{
	cout << "\tDetecting Somas ... ";	


	int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;
	register int i,j,k;

	CImage tempImage(iRows, iCols);
	CImage tempImageXZ(iSlices, iCols);
	CImage tempImageYZ(iSlices, iRows);

	int iWidth = (int) (gfWidthSum / giNumOfWidthSumMembers + 0.5);

	int StructElemSize = iWidth + 7; // it was + 2 11-5-99

	CPoint* aDiskPoints = new CPoint[4 * StructElemSize* StructElemSize];
	int iNumOfPoints = 0;

	ConstructStructElem(StructElemSize, aDiskPoints, iNumOfPoints);

	// DetectSomas(ProjectionImage, &tempImage, StructElemSize);
	DetectSomas(ProjectionImage, & tempImage, aDiskPoints, iNumOfPoints);	

	// estimate the threshold from the brightest region in the image
	int Threshold = static_cast<int>(giForegroundMedian + gfForegroundStdDev);
	int aiHistogram[256];
	tempImage.Histogram(aiHistogram);

	for (i = 255; i > giForegroundMedian; i--)
	{
		if (aiHistogram[i] >= iNumOfPoints)
		{
			Threshold = static_cast<int>(0.6 * i);
			break;
		}
	}

	// use the histogram of seed points to estimate the median and the stddev of 
	// the foreground pixels

	if (Threshold > 250)
		Threshold = 250;


	ProjectionImage->ThresholdImage(Threshold, 0);

	tempImage.ThresholdImage(Threshold, 0);

	CImage tempImageXY(tempImage);
	//Yousef: Try this
	//tempImage.Write("somaXY.pgm");

	// The soma is defined as the region of intersection betweent the three images
	int SomaVolume = 0;
	
	int* aiPlaneSum = new int[iSlices];
	
	int iDenom = iRows;
	if (iCols < iRows)
		iDenom = iCols;

	StructElemSize = static_cast<int>(StructElemSize * (iSlices / iDenom));
	/*if (2 * iSlices < iDenom)
		StructElemSize = static_cast<int>(StructElemSize * (2.0 * iSlices / iDenom));

	if (StructElemSize < iWidth)
		StructElemSize = iWidth + 1;*/

	ConstructStructElem(StructElemSize, aDiskPoints, iNumOfPoints);


	DetectSomas(TrackImageYZ, & tempImageYZ, aDiskPoints, iNumOfPoints);	
	tempImageYZ.Histogram(aiHistogram);
	for (i = 255; i > giForegroundMedian; i--)
	{
		if (aiHistogram[i] >= iNumOfPoints)
		{
			Threshold = static_cast<int>(0.6 * i);
			break;
		}
	}
	tempImageYZ.ThresholdImage(Threshold, 0);

	//Yousef: Try this
	//tempImageYZ.Write("somaYZ.pgm");

	DetectSomas(TrackImageXZ, & tempImageXZ, aDiskPoints, iNumOfPoints);	
	tempImageXZ.Histogram(aiHistogram);
	for (i = 255; i > giForegroundMedian; i--)
	{
		if (aiHistogram[i] >= iNumOfPoints)
		{
			Threshold = static_cast<int>(0.6 * i);
			break;
		}
	}			

	delete [] aDiskPoints;
	delete [] aiPlaneSum;

	//Yousef: Try this
	//tempImageXZ.Write("somaXZ.pgm");

	// count the number of pixels in the soma
	giSomaVolume = 0;
	for (k = 0; k < iSlices; k++)
	{
		for (i = 0; i < iRows; i++)
		{
			for (j = 0; j < iCols; j++)
			{
				if (tempImageXY.data[i][j] &&
					tempImageXZ.data[k][j] &&
					tempImageYZ.data[k][i])
				{
					giSomaVolume++;
				}
			}
		}
	}
	

	gaSomaPoints = new CPoint[giSomaVolume];

	register int index = 0;
	for (k = 0; k < iSlices; k++)
	{
		for (i = 0; i < iRows; i++)
		{
			for (j = 0; j < iCols; j++)
			{
				//if (Soma3DImage->data[k][i][j])
				if (tempImageXY.data[i][j] &&
					tempImageXZ.data[k][j] &&
					tempImageYZ.data[k][i])
				{					
					gaSomaPoints[index].m_iX = j;
					gaSomaPoints[index].m_iY = i;
					gaSomaPoints[index].m_iZ = k;						
					TrackImageXY->data[i][j] = SomaColor;
					TrackImageXZ->data[k][j] = SomaColor;
					TrackImageYZ->data[k][i] = SomaColor;
					The3DImage->data[k][i][j] = SomaColor;
					SomaVolume++;
					index++;
				}
			}
		}
	}
	
	if(SomaVolume == 0)
	{
		int minX = iCols + 1, minY = iRows + 1, minZ = iSlices + 1;
		int maxX = - 1;
		int maxY = - 1;
		int maxZ = - 1;
		CPoint centerPoint;
		centerPoint.m_iX = static_cast<int>(static_cast<float>(minX + maxX) / 2.0 + 0.5);
		centerPoint.m_iY = static_cast<int>(static_cast<float>(minY + maxY) / 2.0 + 0.5);
		centerPoint.m_iZ = static_cast<int>(static_cast<float>(minZ + maxZ) / 2.0 + 0.5); 

		double radius = (maxX - minX) * (maxX - minX) +
			(maxY - minY) * (maxY - minY);
		radius = static_cast<int>(sqrt(radius) / 2.0 + 0.5);

		// create the soma collection
		giNumOfSomas = 1;
		gTheSomas = new CSomas(giNumOfSomas);
		gTheSomas->m_aData[0].m_iVolume = SomaVolume;
		gTheSomas->m_aData[0].m_Center = centerPoint;
		giNumOfSomas = 0;
		cout<<giNumOfSomas<<" somas detected"<<endl;
		return;
	}
	//Yousef: 01-26-2006
	//Apply connected components to seperate somas
	int x, y, z, x_min, y_min, z_min, x_max, y_max, z_max, More_Somas, Soma_Not_Full, soma_ID, br;
	More_Somas = soma_ID = br = 1;
	giNumOfSomas = 1;
	SomaLabelsImage->data[gaSomaPoints[0].m_iZ][gaSomaPoints[0].m_iY][gaSomaPoints[0].m_iX] = 1;

	
	while(More_Somas)
	{
		More_Somas = 0;
		Soma_Not_Full = 1;		
		while(Soma_Not_Full)
		{
			Soma_Not_Full = 0;
			for (i=0; i<index; i++)
			{
				if(SomaLabelsImage->data[gaSomaPoints[i].m_iZ][gaSomaPoints[i].m_iY][gaSomaPoints[i].m_iX] != 0)
					continue;
				x_min = max(0,gaSomaPoints[i].m_iX-1);
				x_max = min(iCols,gaSomaPoints[i].m_iX+1);
				y_min = max(0,gaSomaPoints[i].m_iY-1);
				y_max = min(iRows,gaSomaPoints[i].m_iY+1);
				z_min = max(0,gaSomaPoints[i].m_iZ-1);
				z_max = min(iSlices,gaSomaPoints[i].m_iZ+1);
				
				for(z=z_min; z<=z_max; z++)
				{
					br = 0;
					for(y=y_min; y<=y_max; y++)
					{
						for(x=x_min; x<=x_max; x++)
						{
							if(SomaLabelsImage->data[z][y][x] == soma_ID)
							{
								SomaLabelsImage->data[gaSomaPoints[i].m_iZ][gaSomaPoints[i].m_iY][gaSomaPoints[i].m_iX] = soma_ID;
								Soma_Not_Full = 1;
								br = 1;
								break;
							}
						}
						if(br == 1)
							break;
					}
					if(br == 1)
						break;
				}
			}			
		}
		//See if there is any soma point that is not assigned a label
		 for (i=0; i<index; i++)
		{
			if(SomaLabelsImage->data[gaSomaPoints[i].m_iZ][gaSomaPoints[i].m_iY][gaSomaPoints[i].m_iX] == 0)
			{
				soma_ID++;
				giNumOfSomas++;
				SomaLabelsImage->data[gaSomaPoints[i].m_iZ][gaSomaPoints[i].m_iY][gaSomaPoints[i].m_iX] = soma_ID;
				More_Somas = 1;
				break;
			}
		}
	}

	int* SomaVolumes  = new int[giNumOfSomas];
	int* SomaCentersX = new int[giNumOfSomas];
	int* SomaCentersY = new int[giNumOfSomas];
	int* SomaCentersZ = new int[giNumOfSomas];
	for (i=0; i<soma_ID; i++)
	{
		SomaVolumes[i] = SomaCentersX[i] = SomaCentersY[i] = SomaCentersZ[i] = 0;
	}
	for( i=0; i<index; i++)
	{
		SomaVolumes[SomaLabelsImage->data[gaSomaPoints[i].m_iZ][gaSomaPoints[i].m_iY][gaSomaPoints[i].m_iX]-1]++;
		SomaCentersX[SomaLabelsImage->data[gaSomaPoints[i].m_iZ][gaSomaPoints[i].m_iY][gaSomaPoints[i].m_iX]-1]+=gaSomaPoints[i].m_iX;
		SomaCentersY[SomaLabelsImage->data[gaSomaPoints[i].m_iZ][gaSomaPoints[i].m_iY][gaSomaPoints[i].m_iX]-1]+=gaSomaPoints[i].m_iY;
		SomaCentersZ[SomaLabelsImage->data[gaSomaPoints[i].m_iZ][gaSomaPoints[i].m_iY][gaSomaPoints[i].m_iX]-1]+=gaSomaPoints[i].m_iZ;
	}
	for (i=0; i<soma_ID; i++)
	{
		SomaCentersX[i] = (int) (SomaCentersX[i] / SomaVolumes[i]);
		SomaCentersY[i] = (int) (SomaCentersY[i] / SomaVolumes[i]);
		SomaCentersZ[i] = (int) (SomaCentersZ[i] / SomaVolumes[i]);
	}

	gTheSomas = new CSomas(giNumOfSomas);	
	for (i=0; i<soma_ID; i++)
	{
		gTheSomas->m_aData[i].m_iVolume = SomaVolumes[i];
		gTheSomas->m_aData[i].m_Center.m_iX = SomaCentersX[i];
		gTheSomas->m_aData[i].m_Center.m_iY = SomaCentersY[i];
		gTheSomas->m_aData[i].m_Center.m_iZ = SomaCentersZ[i];
	}


	
	cout<<giNumOfSomas<<" somas detected"<<endl;
	TrackImageXY->Write("somamm.pgm");	
}



//By Yousef
//Yet another version...
void LocateSomas3_v3()
{
	cout << "\tDetecting Somas ... ";	

	int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;
	register int i,j,k;
	int iWidth = (int) (gfWidthSum / giNumOfWidthSumMembers + 0.5);
	Soma3DImage = new C3DImage(iSlices,iRows,iCols);
	
	int StructElemSize = iWidth; // it was + 2 11-5-99

	CPoint* aSphrPoints = new CPoint[(2*StructElemSize+1)* (2*StructElemSize+1) * (2*StructElemSize+1)];
	int iNumOfPoints = 0;

	Construct3DStructElem(StructElemSize, aSphrPoints, iNumOfPoints);
	Detect3DSomas(The3DImage, Soma3DImage, aSphrPoints, iNumOfPoints, StructElemSize);	

	Soma3DImage->Write("Somas.pic");
	
	
	// estimate the threshold from the brightest region in the image
	int *mdn = new int[iSlices];
	float *std = new float[iSlices];
	int Total_mdn = 0;
	float Total_std = 0;

	for(i=0;i<iSlices;i++)
	{
		Soma3DImage->ComputeSliceStatistics(i,mdn[i],std[i]);
		Total_mdn += mdn[i];
		Total_std += std[i];
	}

	Total_mdn = (int) (Total_mdn / iSlices);
	Total_std = Total_std/ (float)iSlices;
	
	int Threshold = static_cast<int>(Total_mdn + Total_std);
	//int Threshold = static_cast<int>(giForegroundMedian + gfForegroundStdDev);		

	Soma3DImage->CreateHistogram();
	for (i = 255; i > giForegroundMedian; i--)
	{
		if (Soma3DImage->m_aiHistogram[i] >= iNumOfPoints)
		{
			Threshold = static_cast<int>(0.6 * i);
			break;
		}
	}

	// use the histogram of seed points to estimate the median and the stddev of 
	// the foreground pixels
	if (Threshold > 250)
		Threshold = 250;

	Soma3DImage->ThresholdImage(Threshold);
	Soma3DImage->Write("Somas.pic");
			
	delete [] aSphrPoints;


	giSomaVolume = 0;

	// count the number of pixels in the soma
	for (k = 0; k < iSlices; k++)
	{
		for (i = 0; i < iRows; i++)
		{
			for (j = 0; j < iCols; j++)
			{
				if (Soma3DImage->data[k][i][j])
				{
					giSomaVolume++;
				}
			}
		}
	}
		
	gaSomaPoints = new CPoint[giSomaVolume];
	int SomaVolume = 0;
	register int index = 0;
	for (k = 0; k < iSlices; k++)
	{
		for (i = 0; i < iRows; i++)
		{
			for (j = 0; j < iCols; j++)
			{
				if (Soma3DImage->data[k][i][j])				
				{					
					gaSomaPoints[index].m_iX = j;
					gaSomaPoints[index].m_iY = i;
					gaSomaPoints[index].m_iZ = k;						
					TrackImageXY->data[i][j] = SomaColor;
					TrackImageXZ->data[k][j] = SomaColor;
					TrackImageYZ->data[k][i] = SomaColor;
					The3DImage->data[k][i][j] = SomaColor;
					SomaVolume++;
					index++;
				}
			}
		}
	}


	if(SomaVolume == 0)
	{
		int minX = iCols + 1, minY = iRows + 1, minZ = iSlices + 1;
		int maxX = - 1;
		int maxY = - 1;
		int maxZ = - 1;
		CPoint centerPoint;
		centerPoint.m_iX = static_cast<int>(static_cast<float>(minX + maxX) / 2.0 + 0.5);
		centerPoint.m_iY = static_cast<int>(static_cast<float>(minY + maxY) / 2.0 + 0.5);
		centerPoint.m_iZ = static_cast<int>(static_cast<float>(minZ + maxZ) / 2.0 + 0.5); 

		double radius = (maxX - minX) * (maxX - minX) +
			(maxY - minY) * (maxY - minY);
		radius = static_cast<int>(sqrt(radius) / 2.0 + 0.5);

		// create the soma collection
		giNumOfSomas = 1;
		gTheSomas = new CSomas(giNumOfSomas);
		gTheSomas->m_aData[0].m_iVolume = SomaVolume;
		gTheSomas->m_aData[0].m_Center = centerPoint;
		giNumOfSomas = 0;
		cout<<giNumOfSomas<<" somas detected"<<endl;
		return;
	}
	//Yousef: 01-26-2006
	//Apply connected components to seperate somas
	int x, y, z, x_min, y_min, z_min, x_max, y_max, z_max, More_Somas, Soma_Not_Full, soma_ID, br;
	More_Somas = soma_ID = br = 1;
	giNumOfSomas = 1;
	SomaLabelsImage->data[gaSomaPoints[0].m_iZ][gaSomaPoints[0].m_iY][gaSomaPoints[0].m_iX] = 1;

	
	while(More_Somas)
	{
		More_Somas = 0;
		Soma_Not_Full = 1;		
		while(Soma_Not_Full)
		{
			Soma_Not_Full = 0;
			for (i=0; i<index; i++)
			{
				if(SomaLabelsImage->data[gaSomaPoints[i].m_iZ][gaSomaPoints[i].m_iY][gaSomaPoints[i].m_iX] != 0)
					continue;
				x_min = max(0,gaSomaPoints[i].m_iX-1);
				x_max = min(iCols,gaSomaPoints[i].m_iX+1);
				y_min = max(0,gaSomaPoints[i].m_iY-1);
				y_max = min(iRows,gaSomaPoints[i].m_iY+1);
				z_min = max(0,gaSomaPoints[i].m_iZ-1);
				z_max = min(iSlices,gaSomaPoints[i].m_iZ+1);
				
				for(z=z_min; z<=z_max; z++)
				{
					br = 0;
					for(y=y_min; y<=y_max; y++)
					{
						for(x=x_min; x<=x_max; x++)
						{
							if(SomaLabelsImage->data[z][y][x] == soma_ID)
							{
								SomaLabelsImage->data[gaSomaPoints[i].m_iZ][gaSomaPoints[i].m_iY][gaSomaPoints[i].m_iX] = soma_ID;
								Soma_Not_Full = 1;
								br = 1;
								break;
							}
						}
						if(br == 1)
							break;
					}
					if(br == 1)
						break;
				}
			}			
		}
		//See if there is any soma point that is not assigned a label
		 for (i=0; i<index; i++)
		{
			if(SomaLabelsImage->data[gaSomaPoints[i].m_iZ][gaSomaPoints[i].m_iY][gaSomaPoints[i].m_iX] == 0)
			{
				soma_ID++;
				giNumOfSomas++;
				SomaLabelsImage->data[gaSomaPoints[i].m_iZ][gaSomaPoints[i].m_iY][gaSomaPoints[i].m_iX] = soma_ID;
				More_Somas = 1;
				break;
			}
		}
	}

	int* SomaVolumes  = new int[giNumOfSomas];
	int* SomaCentersX = new int[giNumOfSomas];
	int* SomaCentersY = new int[giNumOfSomas];
	int* SomaCentersZ = new int[giNumOfSomas];
	for (i=0; i<soma_ID; i++)
	{
		SomaVolumes[i] = SomaCentersX[i] = SomaCentersY[i] = SomaCentersZ[i] = 0;
	}
	for( i=0; i<index; i++)
	{
		SomaVolumes[SomaLabelsImage->data[gaSomaPoints[i].m_iZ][gaSomaPoints[i].m_iY][gaSomaPoints[i].m_iX]-1]++;
		SomaCentersX[SomaLabelsImage->data[gaSomaPoints[i].m_iZ][gaSomaPoints[i].m_iY][gaSomaPoints[i].m_iX]-1]+=gaSomaPoints[i].m_iX;
		SomaCentersY[SomaLabelsImage->data[gaSomaPoints[i].m_iZ][gaSomaPoints[i].m_iY][gaSomaPoints[i].m_iX]-1]+=gaSomaPoints[i].m_iY;
		SomaCentersZ[SomaLabelsImage->data[gaSomaPoints[i].m_iZ][gaSomaPoints[i].m_iY][gaSomaPoints[i].m_iX]-1]+=gaSomaPoints[i].m_iZ;
	}
	for (i=0; i<soma_ID; i++)
	{
		SomaCentersX[i] = (int) (SomaCentersX[i] / SomaVolumes[i]);
		SomaCentersY[i] = (int) (SomaCentersY[i] / SomaVolumes[i]);
		SomaCentersZ[i] = (int) (SomaCentersZ[i] / SomaVolumes[i]);
	}

	gTheSomas = new CSomas(giNumOfSomas);	
	for (i=0; i<soma_ID; i++)
	{
		gTheSomas->m_aData[i].m_iVolume = SomaVolumes[i];
		gTheSomas->m_aData[i].m_Center.m_iX = SomaCentersX[i];
		gTheSomas->m_aData[i].m_Center.m_iY = SomaCentersY[i];
		gTheSomas->m_aData[i].m_Center.m_iZ = SomaCentersZ[i];
	}


	
	cout<<giNumOfSomas<<" somas detected"<<endl;
}


//By Yousef
//Last version... The soma image is provided
void LocateSomas3_v4()
{
	cout << "\tDetecting Somas ... ";	

	int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;
	register int i,j,k;
	giSomaVolume = 0;

	// count the number of pixels in the soma
	for (k = 0; k < iSlices; k++)
	{
		for (i = 0; i < iRows; i++)
		{
			for (j = 0; j < iCols; j++)
			{
				if (Soma3DImage->data[k][i][j])
				{
					giSomaVolume++;
				}
			}
		}
	}
		
	gaSomaPoints = new CPoint[giSomaVolume];
	int SomaVolume = 0;
	register int index = 0;
	for (k = 0; k < iSlices; k++)
	{
		for (i = 0; i < iRows; i++)
		{
			for (j = 0; j < iCols; j++)
			{
				if (Soma3DImage->data[k][i][j])				
				{					
					gaSomaPoints[index].m_iX = j;
					gaSomaPoints[index].m_iY = i;
					gaSomaPoints[index].m_iZ = k;						
					TrackImageXY->data[i][j] = SomaColor;
					TrackImageXZ->data[k][j] = SomaColor;
					TrackImageYZ->data[k][i] = SomaColor;
					The3DImage->data[k][i][j] = SomaColor;
					SomaVolume++;
					index++;
				}
			}
		}
	}


	//Yousef: Try this
	//TrackImageXY->Write("trXY.pgm");
	//TrackImageXZ->Write("trXZ.pgm");
	//TrackImageYZ->Write("trYZ.pgm");
	if(SomaVolume == 0)
	{
		int minX = iCols + 1, minY = iRows + 1, minZ = iSlices + 1;
		int maxX = - 1;
		int maxY = - 1;
		int maxZ = - 1;
		CPoint centerPoint;
		centerPoint.m_iX = static_cast<int>(static_cast<float>(minX + maxX) / 2.0 + 0.5);
		centerPoint.m_iY = static_cast<int>(static_cast<float>(minY + maxY) / 2.0 + 0.5);
		centerPoint.m_iZ = static_cast<int>(static_cast<float>(minZ + maxZ) / 2.0 + 0.5); 

		double radius = (maxX - minX) * (maxX - minX) +
			(maxY - minY) * (maxY - minY);
		radius = static_cast<int>(sqrt(radius) / 2.0 + 0.5);

		// create the soma collection
		giNumOfSomas = 1;
		gTheSomas = new CSomas(giNumOfSomas);
		gTheSomas->m_aData[0].m_iVolume = SomaVolume;
		gTheSomas->m_aData[0].m_Center = centerPoint;
		giNumOfSomas = 0;
		cout<<giNumOfSomas<<" somas detected"<<endl;
		return;
	}
	//Yousef: 01-26-2006
	//Apply connected components to seperate somas
	int x, y, z, x_min, y_min, z_min, x_max, y_max, z_max, More_Somas, Soma_Not_Full, soma_ID, br;
	More_Somas = soma_ID = br = 1;
	giNumOfSomas = 1;
	SomaLabelsImage->data[gaSomaPoints[0].m_iZ][gaSomaPoints[0].m_iY][gaSomaPoints[0].m_iX] = 1;

	
	while(More_Somas)
	{
		More_Somas = 0;
		Soma_Not_Full = 1;		
		while(Soma_Not_Full)
		{
			Soma_Not_Full = 0;
			for (i=0; i<index; i++)
			{
				if(SomaLabelsImage->data[gaSomaPoints[i].m_iZ][gaSomaPoints[i].m_iY][gaSomaPoints[i].m_iX] != 0)
					continue;
				x_min = max(0,gaSomaPoints[i].m_iX-1);
				x_max = min(iCols-1,gaSomaPoints[i].m_iX+1);
				y_min = max(0,gaSomaPoints[i].m_iY-1);
				y_max = min(iRows-1,gaSomaPoints[i].m_iY+1);
				z_min = max(0,gaSomaPoints[i].m_iZ-1);
				z_max = min(iSlices-1,gaSomaPoints[i].m_iZ+1);
				
				for(z=z_min; z<=z_max; z++)
				{
					br = 0;
					for(y=y_min; y<=y_max; y++)
					{
						for(x=x_min; x<=x_max; x++)
						{
							if(SomaLabelsImage->data[z][y][x] == soma_ID)
							{
								SomaLabelsImage->data[gaSomaPoints[i].m_iZ][gaSomaPoints[i].m_iY][gaSomaPoints[i].m_iX] = soma_ID;
								Soma_Not_Full = 1;
								br = 1;
								break;
							}
						}
						if(br == 1)
							break;
					}
					if(br == 1)
						break;
				}
			}			
		}
		//See if there is any soma point that is not assigned a label
		 for (i=0; i<index; i++)
		{
			if(SomaLabelsImage->data[gaSomaPoints[i].m_iZ][gaSomaPoints[i].m_iY][gaSomaPoints[i].m_iX] == 0)
			{
				soma_ID++;
				giNumOfSomas++;
				SomaLabelsImage->data[gaSomaPoints[i].m_iZ][gaSomaPoints[i].m_iY][gaSomaPoints[i].m_iX] = soma_ID;
				More_Somas = 1;
				break;
			}
		}
	}

	int* SomaVolumes  = new int[giNumOfSomas];
	int* SomaCentersX = new int[giNumOfSomas];
	int* SomaCentersY = new int[giNumOfSomas];
	int* SomaCentersZ = new int[giNumOfSomas];
	for (i=0; i<soma_ID; i++)
	{
		SomaVolumes[i] = SomaCentersX[i] = SomaCentersY[i] = SomaCentersZ[i] = 0;
	}
	for( i=0; i<index; i++)
	{
		SomaVolumes[SomaLabelsImage->data[gaSomaPoints[i].m_iZ][gaSomaPoints[i].m_iY][gaSomaPoints[i].m_iX]-1]++;
		SomaCentersX[SomaLabelsImage->data[gaSomaPoints[i].m_iZ][gaSomaPoints[i].m_iY][gaSomaPoints[i].m_iX]-1]+=gaSomaPoints[i].m_iX;
		SomaCentersY[SomaLabelsImage->data[gaSomaPoints[i].m_iZ][gaSomaPoints[i].m_iY][gaSomaPoints[i].m_iX]-1]+=gaSomaPoints[i].m_iY;
		SomaCentersZ[SomaLabelsImage->data[gaSomaPoints[i].m_iZ][gaSomaPoints[i].m_iY][gaSomaPoints[i].m_iX]-1]+=gaSomaPoints[i].m_iZ;
	}
	for (i=0; i<soma_ID; i++)
	{
		SomaCentersX[i] = (int) (SomaCentersX[i] / SomaVolumes[i]);
		SomaCentersY[i] = (int) (SomaCentersY[i] / SomaVolumes[i]);
		SomaCentersZ[i] = (int) (SomaCentersZ[i] / SomaVolumes[i]);
	}

	gTheSomas = new CSomas(giNumOfSomas);	
	for (i=0; i<soma_ID; i++)
	{
		gTheSomas->m_aData[i].m_iVolume = SomaVolumes[i];
		gTheSomas->m_aData[i].m_Center.m_iX = SomaCentersX[i];
		gTheSomas->m_aData[i].m_Center.m_iY = SomaCentersY[i];
		gTheSomas->m_aData[i].m_Center.m_iZ = SomaCentersZ[i];
	}


	
	cout<<giNumOfSomas<<" somas detected"<<endl;
}

