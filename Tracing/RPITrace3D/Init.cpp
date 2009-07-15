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

////////////////////////////////////////////////////////////////////
//  File: Init.cpp
//  Created: 4-3-98
//
//  Contains files necessary to initialize the global templates and 
//  shift vectors.
///////////////////////////////////////////////////////////////////
//#pragma warning(disable:4786)  // disable STL-related warnings
//#pragma warning(disable:4710)  // compiler didn't inline function that got selected for inline expansion
//#pragma warning(disable:4514)  // The optimizer removed an inline function that is not called

#include <fstream>
#include <iostream>
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

using namespace std;

class CVessel;
class CVessels;
#include "Template.h"
#include "Extern.h"

// the first dimension stands for 45, 135, 225, and 315
// the second dimension stands for minus2, minus1, plus2, and plus1
int FourtyFiveDegreeVectorsLeft[4][4];
int FourtyFiveDegreeVectorsRight[4][4];

void InitData()
{
	//int iSlices = The3DImage->m_iSlices;
  //int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;

	// the first dimension stands for 45, 135, 225, and 315
	// the second dimension stands for minus2, minus1, plus2, and plus1
	FourtyFiveDegreeVectorsLeft[0][0] = 1;
	FourtyFiveDegreeVectorsLeft[0][1] = iCols + 1;
	FourtyFiveDegreeVectorsLeft[0][2] = -1 * iCols;
	FourtyFiveDegreeVectorsLeft[0][3] = -1 * iCols - 1;

	FourtyFiveDegreeVectorsLeft[1][0] = -1 * iCols;
	FourtyFiveDegreeVectorsLeft[1][1] = -1 * iCols + 1;
	FourtyFiveDegreeVectorsLeft[1][2] = -1;
	FourtyFiveDegreeVectorsLeft[1][3] = iCols - 1;


	FourtyFiveDegreeVectorsLeft[2][0] = -1;
	FourtyFiveDegreeVectorsLeft[2][1] = -1 * iCols - 1;
	FourtyFiveDegreeVectorsLeft[2][2] = iCols;
	FourtyFiveDegreeVectorsLeft[2][3] = iCols + 1;


	FourtyFiveDegreeVectorsLeft[3][0] = iCols;
	FourtyFiveDegreeVectorsLeft[3][1] = iCols - 1;
	FourtyFiveDegreeVectorsLeft[3][2] = 1;
	FourtyFiveDegreeVectorsLeft[3][3] = -1 * iCols + 1;

	FourtyFiveDegreeVectorsRight[0][0] = -1 * iCols;
	FourtyFiveDegreeVectorsRight[0][1] = -1 * iCols - 1;
	FourtyFiveDegreeVectorsRight[0][2] = 1;
	FourtyFiveDegreeVectorsRight[0][3] = iCols + 1;

	FourtyFiveDegreeVectorsRight[1][0] = -1;
	FourtyFiveDegreeVectorsRight[1][1] = iCols - 1;
	FourtyFiveDegreeVectorsRight[1][2] = -1 * iCols;
	FourtyFiveDegreeVectorsRight[1][3] = -1 * iCols + 1;

	FourtyFiveDegreeVectorsRight[2][0] = iCols;
	FourtyFiveDegreeVectorsRight[2][1] = iCols + 1;
	FourtyFiveDegreeVectorsRight[2][2] = -1;
	FourtyFiveDegreeVectorsRight[2][3] = -1 * iCols - 1;

	FourtyFiveDegreeVectorsRight[3][0] = 1;
	FourtyFiveDegreeVectorsRight[3][1] = -1 * iCols + 1;
	FourtyFiveDegreeVectorsRight[3][2] = iCols;
	FourtyFiveDegreeVectorsRight[3][3] = iCols - 1;
}


////////////////////////
// Function: Initialize
//
// Create the shift vectors, left and right templates
// that are defined globally.

void Initialize()
{
	int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;
	
	cout << "\tInitializing variables ... ";
	// working variables
	register int i, j, k;
	double theta = DirectionStep;
	float Hdir, Vdir;

	InitData();
	// left and right templates and shift vectors
	for (i = 0; i < NumOfDirections; i++)
	{
		for (j = 0; j < NumOfDirections; j++)
		{
			Hdir = static_cast<float>(i * theta);
			Vdir = static_cast<float>(j * theta);
			gVectorsArray[i][j] = new CVector(iSlices,
									  	iRows,
									  	iCols,
									  	Hdir,
									  	Vdir,
									  	MaxVectorLength);
		}
	}

	for (i = 0; i < NumOfDirections; i++)
	{
		for (j = 0; j < NumOfDirections; j++)
		{
			Hdir = static_cast<float>(i * theta);
			Vdir = static_cast<float>(j * theta);

			gHLeftTemplatesArray[i][j] = new CHLeftTemplate(iSlices,
											 	iRows,
											 	iCols,
											 	Hdir,
											 	Vdir,
											 	MaxTemplateLength);


			gVLeftTemplatesArray[i][j] = new CVLeftTemplate(iSlices,
											 	iRows,
											 	iCols,
											 	Hdir,
											 	Vdir,
											 	MaxTemplateLength);

			gHRightTemplatesArray[i][j] = new CHRightTemplate(iSlices,
											  	iRows,
											  	iCols,
											  	Hdir,
											  	Vdir,
											  	MaxTemplateLength);

			gVRightTemplatesArray[i][j] = new CVRightTemplate(iSlices,
											  	iRows,
											  	iCols,
											  	Hdir,
											  	Vdir,
											  	MaxTemplateLength);
		}
	}

	// create and initialize the traced image.
	int size = iSlices* iRows* iCols;
	TracedImage = new bool * *[size];
	for (i = 0; i < iSlices; i++)
	{
		TracedImage[i] = new bool * [iRows];
		for (j = 0; j < iRows; j++)
		{
			TracedImage[i][j] = new bool[iCols];
			for (k = 0; k < iCols; k++)
				TracedImage[i][j][k] = false;
		}
	}

	// allocate space for the 2D array of seed points
	CPoint** tempSeedsArray = new CPoint *[iRows* iCols];
	memset(tempSeedsArray, 0, sizeof(CPoint *) * iCols * iRows);
	gapArrayOfSeedPoints = new CPoint * *[iRows];
	for (j = 0; j < iRows; j++)
	{
		gapArrayOfSeedPoints[j] = &(tempSeedsArray[j * iCols]);
	}

	giMedians = new int[iSlices];
	gfStdDevs = new float[iSlices];
	memset(giMedians, 0, sizeof(int) * iSlices);
	memset(gfStdDevs, 0, sizeof(float) * iSlices);

	// create an array of initial points queues, one for each image slice
	//	InitialPoints = new CQueue<CPoint> [iSlices];

	cout << "Done. " << endl;
}
