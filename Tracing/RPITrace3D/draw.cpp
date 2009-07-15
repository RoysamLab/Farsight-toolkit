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
// FILE: draw.cpp
// 
// This file contains some drawing tools available for use of other functions
//
//#pragma warning(disable:4786)  // disable STL-related warnings
//#pragma warning(disable:4710)  // compiler didn't inline function that got selected for inline expansion
//#pragma warning(disable:4514)  // The optimizer removed an inline function that is not called
//#pragma warning(disable:4702)  // unreachable STLport code

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
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
#include "Soma.h"

class CSomas;
using namespace std;

void Draw_Projections()
{
	string file_name;
	string output_path = gConfig.GetOutputPath();
	string image_name = gConfig.GetImageName();

	if (cfg_output_projections)
	{
		if (strcmp(cfg_output_image_format, "tif") == 0 ||
			cfg_output_image_fgcolor)
		{
			file_name = output_path + image_name + "XY.tif";
			CanvasXY->WriteTIFF(file_name);
			file_name = output_path + image_name + "XZ.tif";
			CanvasXZ->WriteTIFF(file_name);
			file_name = output_path + image_name + "YZ.tif";
			CanvasYZ->WriteTIFF(file_name);
		}
		else
		{
			if (strcmp(cfg_output_image_format, "pgm") == 0)
			{
				file_name = output_path + image_name + "XY.pgm";
				CanvasXY->Write(file_name);
				file_name = output_path + image_name + "XZ.pgm";
				CanvasXZ->Write(file_name);
				file_name = output_path + image_name + "YZ.pgm";
				CanvasYZ->Write(file_name);
			}
		}
	}
}

void Draw_Centerline3D()
{
	string output_path = gConfig.GetOutputPath();
	string image_name = gConfig.GetImageName();
	string file_name;
	C3DImage* VolumeImage;
	VolumeImage = new C3DImage(*The3DImage);
	VolumeImage->Saturate(0, 250);

	gTheVessels.Draw3DCenterline(*VolumeImage, 255);

	for (int i = 0; i < gIntersectionPoints.m_iNumOfElements; i++)
	{
		CPoint* pPoint =& gIntersectionPoints.m_apData[i]->m_Point;
		VolumeImage->MarkCrosshair(pPoint, 254);
	}
	file_name = output_path + image_name + "Trace.pic";
	VolumeImage->Write(file_name);
}

void Draw_CenterlineOnProjections()
{
	string output_path = gConfig.GetOutputPath();
	string image_name = gConfig.GetImageName();
	string file_name;
	int i;
	CImage* XY = NULL;
	CImage* XZ = NULL;
	CImage* YZ = NULL;

	if (gConfig.GetDetectSoma() && giNumOfSomas != 0)
	{
		extern CSomas* gTheSomas;

		if (cfg_output_soma_draw_trees)
		{
			gTheSomas->ConstructTrees();
			// remove all pixels interfering with tree coloring
			//int iNumOfTrees = gTheSomas->m_aData[0].m_iNumOfIntersectionPoints;
			//TrackImageXY->RemovePixels(iNumOfTrees + 5 + 2);
			//TrackImageXZ->RemovePixels(iNumOfTrees + 5 + 2);
			//TrackImageYZ->RemovePixels(iNumOfTrees + 5 + 2);
			gTheSomas->DrawTrees();
		}

		XY = new CImage(*TrackImageXY);
		XZ = new CImage(*TrackImageXZ);
		YZ = new CImage(*TrackImageYZ);
		gTheSomas->WriteIDs(*XY, LetterColor);
		gTheSomas->WriteIDsXZ(*XZ, LetterColor);
		gTheSomas->WriteIDsYZ(*YZ, LetterColor);
	}
	else
	{
		XY = new CImage(*CanvasXY);
		XZ = new CImage(*CanvasXZ);
		YZ = new CImage(*CanvasYZ);
	}

	//	XY = new CImage("M3B2.pgm");

	// draw vessels onto the projection images



	gTheVessels.DrawVessels(*XY, CenterlineColor);
	gTheVessels.DrawVesselsXZ(*XZ, CenterlineColor);
	gTheVessels.DrawVesselsYZ(*YZ, CenterlineColor);	

	// display the intersection points
	for (i = 0; i < gIntersectionPoints.m_iNumOfElements; i++)
	{
		CPoint* pPoint =& gIntersectionPoints.m_apData[i]->m_Point;

		XY->MarkCrosshairXY(pPoint,IntersectionPointColor);
		XZ->MarkCrosshairXZ(pPoint,IntersectionPointColor);
		YZ->MarkCrosshairYZ(pPoint,IntersectionPointColor);
	}

	// display seed points

	/*list<CPoint>::iterator k;
	for (k = TracedSeeds.begin(); k != TracedSeeds.end(); k++)
	{
		XY->MarkCrosshairXY(&(*k), 1);
		XZ->MarkCrosshairXZ(&(*k), 1);
		YZ->MarkCrosshairYZ(&(*k), 1);
	}*/

	cout << endl;
	if (gConfig.GetDetectSoma())
	{
		cout << "Soma and centerline images: " << endl;
	}
	else
	{
		cout << "Centerline Images: " << endl;
	}

	// we draw either in tif or pgm format
	if (strcmp(cfg_output_image_format, "tif") == 0 ||
		cfg_output_image_fgcolor)
	{
		// draw centerlines onto these images
		file_name = output_path + image_name + "TraceCenterLineXY.tif";
		XY->WriteTIFF(file_name);
		cout << "\t%OutputPath%/" << image_name << "TraceCenterLineXY.tif" << endl;
		file_name = output_path + image_name + "TraceCenterLineXZ.tif";
		XZ->WriteTIFF(file_name);
		cout << "\t%OutputPath%/" << image_name << "TraceCenterLineXZ.tif" << endl;
		file_name = output_path + image_name + "TraceCenterLineYZ.tif";
		YZ->WriteTIFF(file_name);
		cout << "\t%OutputPath%/" << image_name << "TraceCenterLineYZ.tif" << endl;
	}
	else
	{
		file_name = output_path + image_name + "TraceCenterLineXY.pgm";
		XY->Write(file_name);
		cout << "\t%OutputPath%/" << image_name << "TraceCenterLineXY.pgm" << endl;
		file_name = output_path + image_name + "TraceCenterLineXZ.pgm";
		XZ->Write(file_name);
		cout << "\t%OutputPath%/" << image_name << "TraceCenterLineXZ.pgm" << endl;
		file_name = output_path + image_name + "TraceCenterLineYZ.pgm";
		YZ->Write(file_name);
		cout << "\t%OutputPath%/" << image_name << "TraceCenterLineYZ.pgm" << endl;
	}
	cout << endl;
}

void Draw_BranchPoints()
{
	string output_path = gConfig.GetOutputPath();
	string image_name = gConfig.GetImageName();
	string file_name;
	int i, j, k;
	CImage* XY = NULL;

	XY = new CImage(*CanvasXY);

	// display the intersection points
	for (i = 0; i < gIntersectionPoints.m_iNumOfElements; i++)
	{
		CPoint* pPoint =& gIntersectionPoints.m_apData[i]->m_Point;
		
		list<CVessel> connected_vessels;
		for(j = 0; j < gTheVessels.m_iNumOfElements; j++)
		{
			if(gTheVessels.m_apData[j])
			{
				for(k = 0; k < gTheVessels.m_apData[j]->m_iNumOfIntersectionPoints; k++)
				{
					if(gTheVessels.m_apData[j]->m_aiMyIntersectionPoints[k] == gIntersectionPoints.m_apData[i]->m_iID)
					{
						
						// mark 5 pixels away from the point as well
						CLNode<CPoint>* tempNode = NULL;
						CLNode<CPoint>* point_loc = NULL;
						
						// centerline
						tempNode = gTheVessels.m_apData[j]->m_Center.head;
						while (tempNode)
						{
							if(*tempNode->data == *pPoint)
								point_loc = tempNode;
							tempNode = tempNode->after;
						}
						tempNode = point_loc;
						for(int a = 0; a < 10; a++)
						{
							if(tempNode)
								XY->MarkPoint(tempNode->data,CenterlineColor);
							else
								break;
							tempNode = tempNode->after;
						}
						tempNode = point_loc;
						for(int b = 0; b < 10; b++)
						{
							if(tempNode)
								XY->MarkPoint(tempNode->data,CenterlineColor);
							else
								break;
							tempNode = tempNode->before;
						}
						XY->MarkCrosshairXY(pPoint,IntersectionPointColor);
					}
				}
			}
		}
	}

	file_name = output_path + image_name + "BranchPoints.tif";
	XY->WriteTIFF(file_name);

}

void Draw_BorderlineOnProjections(char* plane)
{
	string output_path = gConfig.GetOutputPath();
	string image_name = gConfig.GetImageName();
	string file_name;
	int i;
	CImage* XY;
	CImage* XZ;
	CImage* YZ;

	XY = new CImage(*CanvasXY);
	XZ = new CImage(*CanvasXZ);
	YZ = new CImage(*CanvasYZ);

	if (strcmp(plane, "vertical") && strcmp(plane, "horizontal"))
	{
		cerr << "Error::Draw_BorderlineOnProjections: Invalid plane " << plane
			<< endl;
		exit(0);
	}

	if (strcmp(plane, "vertical") == 0)
	{
		for (i = 0; i < gTheVessels.m_iNumOfElements; i++)
		{
			if (gTheVessels.m_apData[i])
			{
				gTheVessels.m_apData[i]->DrawBoundariesVXY(*XY, 2);
				gTheVessels.m_apData[i]->DrawBoundariesVXZ(*XZ, 2);
				gTheVessels.m_apData[i]->DrawBoundariesVYZ(*YZ, 2);
			}
		}
	}

	if (strcmp(plane, "horizontal") == 0)
	{
		for (i = 0; i < gTheVessels.m_iNumOfElements; i++)
		{
			if (gTheVessels.m_apData[i])
			{
				gTheVessels.m_apData[i]->DrawBoundariesHXY(*XY, 2);
				gTheVessels.m_apData[i]->DrawBoundariesHXZ(*XZ, 2);
				gTheVessels.m_apData[i]->DrawBoundariesHYZ(*YZ, 2);
			}
		}
	}

	// display the intersection points
	for (i = 0; i < gIntersectionPoints.m_iNumOfElements; i++)
	{
		CPoint* pPoint =& gIntersectionPoints.m_apData[i]->m_Point;

		XY->data[pPoint->m_iY][pPoint->m_iX] = IntersectionPointColor;
		XY->data[pPoint->m_iY + 1][pPoint->m_iX] = IntersectionPointColor;
		XY->data[pPoint->m_iY - 1][pPoint->m_iX] = IntersectionPointColor;
		XY->data[pPoint->m_iY][pPoint->m_iX + 1] = IntersectionPointColor;
		XY->data[pPoint->m_iY][pPoint->m_iX - 1] = IntersectionPointColor;

		XZ->data[pPoint->m_iZ][pPoint->m_iX] = IntersectionPointColor;
		XZ->data[pPoint->m_iZ + 1][pPoint->m_iX] = IntersectionPointColor;
		XZ->data[pPoint->m_iZ - 1][pPoint->m_iX] = IntersectionPointColor;
		XZ->data[pPoint->m_iZ][pPoint->m_iX + 1] = IntersectionPointColor;
		XZ->data[pPoint->m_iZ][pPoint->m_iX - 1] = IntersectionPointColor;

		YZ->data[pPoint->m_iZ][pPoint->m_iY] = IntersectionPointColor;		
		YZ->data[pPoint->m_iZ + 1][pPoint->m_iY] = IntersectionPointColor;
		YZ->data[pPoint->m_iZ - 1][pPoint->m_iY] = IntersectionPointColor;
		YZ->data[pPoint->m_iZ][pPoint->m_iY + 1] = IntersectionPointColor;		
		YZ->data[pPoint->m_iZ][pPoint->m_iY - 1] = IntersectionPointColor;
	}


	// we draw either in tif or pgm format
	if (strcmp(cfg_output_image_format, "tif") == 0 ||
		cfg_output_image_fgcolor)
	{
		// draw centerlines onto these images
		file_name = output_path + image_name + "BorderLine_" + plane + "_XY.tif";
		XY->WriteTIFF(file_name);
		file_name = output_path + image_name + "BorderLine_" + plane + "_XZ.tif";
		XZ->WriteTIFF(file_name);
		file_name = output_path + image_name + "BorderLine_" + plane + "_YZ.tif";
		YZ->WriteTIFF(file_name);
	}
	else
	{
		file_name = output_path + image_name + "BorderLine_" + plane + "_XY.pgm";
		XY->Write(file_name);
		file_name = output_path + image_name + "BorderLine_" + plane + "_XZ.pgm";
		XZ->Write(file_name);
		file_name = output_path + image_name + "BorderLine_" + plane + "_YZ.pgm";
		YZ->Write(file_name);
	}
	delete XY;
	delete XZ;
	delete YZ;
}

void Draw_SeedPointsOnProjections(char* gridlines = NULL)
{
	//int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;
	string output_path = gConfig.GetOutputPath();
	string image_name = gConfig.GetImageName();
	int iGridSpacing = gConfig.GetGridSpacing();

  /* come up with a better way to do this...
	// assign default values
	if (gridlines == NULL)
	{
		gridlines = "no";
	}
  */

	string file_name;
	CImage* XY = NULL;
	CImage* XZ = NULL;
	CImage* YZ = NULL;
	CPoint* pPoint = NULL;
	int x,y,i;
	int giMARGIN2 = 3;

	XY = new CImage(*CanvasXY);
	XZ = new CImage(*CanvasXZ);
	YZ = new CImage(*CanvasYZ);

	XY->RemovePixels(5);
	XZ->RemovePixels(5);
	YZ->RemovePixels(5);

	//draw gridlines if requested
	if (strcmp(gridlines, "yes") == 0)
	{
		for (x = giMARGIN2; x < iCols - giMARGIN2; x += iGridSpacing)
		{
			// from the top boundary to the bottome boundary
			for (y = giMARGIN2; y < iRows - giMARGIN2; y++)
			{
				for (i = y + 1 ; i < (y + iGridSpacing); i++)
				{
					XY->data[y][x] = 255;
				}
			}
		}

		for (y = giMARGIN2; y < iRows - giMARGIN2; y += iGridSpacing)
		{
			for (x = giMARGIN; x < iCols - giMARGIN2; x++)
			{
				for (i = x + 1 ; i < x + iGridSpacing; i++)
				{
					XY->data[y][x] = 255;
				}
			}
		}
	}

	for (i = 0; i < giNumOfSeedPoints; i++)
	{
		pPoint = gapSortedArrayOfSeedPoints[i];
		XY->MarkCrosshairXY(pPoint, ValidSeedColor);
		XZ->MarkCrosshairXZ(pPoint, ValidSeedColor);
		YZ->MarkCrosshairYZ(pPoint, ValidSeedColor);
	}

	// we draw either in tif or pgm format
	if (strcmp(cfg_output_image_format, "tif") == 0 ||
		cfg_output_image_fgcolor)
	{
		file_name = output_path + image_name + "SeedPointsXY.tif";
		XY->WriteTIFF(file_name);
		file_name = output_path + image_name + "SeedPointsXZ.tif";
		XZ->WriteTIFF(file_name);
		file_name = output_path + image_name + "SeedPointsYZ.tif";
		YZ->WriteTIFF(file_name);
	}
	else
	{
		file_name = output_path + image_name + "SeedPointsXY.pgm";
		XY->Write(file_name);
		file_name = output_path + image_name + "SeedPointsXZ.pgm";
		XZ->Write(file_name);
		file_name = output_path + image_name + "SeedPointsYZ.pgm";
		YZ->Write(file_name);
	}
}

void Draw_SeedCandidates (std::vector<CPoint> & seed_candidates)
{
	string output_path = gConfig.GetOutputPath();
	string image_name = gConfig.GetImageName();
	string file_name;
	CImage* XY = NULL;

	CPoint* pPoint = NULL;
	//int giMARGIN2 = 3;

	XY = new CImage(*CanvasXY);

	XY->RemovePixels(5);

	int min_X, max_X, min_Y, max_Y;

	for (std::vector<CPoint>::iterator i = seed_candidates.begin (); i != seed_candidates.end (); i++)
	{	
		pPoint = &*i;
		if (!The3DImage->WithinImagePadding(*pPoint,0))
			continue;
		XY->MarkCrosshairXY(pPoint,2);
		int radius = Round( static_cast<double>(i->m_fHWidth) / 2.0 );
		if(i->m_iHDir != 0)
		{
			min_X = max(0,i->m_iX - radius);
			max_X = min(i->m_iX + radius, XY->m_iCols);
			pPoint->m_iY = i->m_iY;
			for (int x = min_X; x <= max_X; x++)
			{
				pPoint->m_iX = x;
				XY->MarkPoint(pPoint, 3);
			}
		}
		else
		{
			min_Y = max(0,i->m_iY - radius);
			max_Y = min(i->m_iY + radius, XY->m_iRows);
			pPoint->m_iX = i->m_iX;
			for (int y = min_Y; y <= max_Y; y++)
			{
				pPoint->m_iY = y;
				if(XY->ValidPoint(*pPoint))
					XY->MarkPoint(pPoint, 3);
			}
		}
	}

	file_name = output_path + image_name + "UnverifiedSeeds.tif";
	XY->WriteTIFF(file_name);

	delete XY; XY = NULL;
}

void Draw_PointsOnProjections(list<CPoint> center, list<CPoint> left,
	list<CPoint> right, char* name = NULL)
{
  /* come up with a better way to handle this case that doesn't violate
     char* = const char*
	if (name == NULL)
		name = "";
  */

	string output_path = gConfig.GetOutputPath();
	string image_name = gConfig.GetImageName();
	string file_name;
	CImage* XY = NULL;
	CImage* XZ = NULL;
	CImage* YZ = NULL;
	CPoint pPoint;

	XY = new CImage(*CanvasXY);
	XZ = new CImage(*CanvasXZ);
	YZ = new CImage(*CanvasYZ);

	XY->RemovePixels(5);
	XZ->RemovePixels(5);
	YZ->RemovePixels(5);

	list<CPoint>::iterator i;

	for (i = center.begin(); i != center.end(); i++)
	{
		pPoint = *i;
		XY->MarkPointXY(&pPoint, 0);
		XZ->MarkPointXZ(&pPoint, 0);
		YZ->MarkPointYZ(&pPoint, 0);
	}
	for (i = left.begin(); i != left.end(); i++)
	{
		pPoint = *i;
		if (The3DImage->WithinImagePadding(pPoint, giMARGIN))
		{
			XY->MarkPointXY(&pPoint, 1);
			XZ->MarkPointXZ(&pPoint, 1);
			YZ->MarkPointYZ(&pPoint, 1);
		}
	}
	for (i = right.begin(); i != right.end(); i++)
	{
		pPoint = *i;
		if (The3DImage->WithinImagePadding(pPoint, giMARGIN))
		{
			XY->MarkPointXY(&pPoint, 2);
			XZ->MarkPointXZ(&pPoint, 2);
			YZ->MarkPointYZ(&pPoint, 2);
		}
	}

	// we draw either in tif or pgm format
	if (strcmp(cfg_output_image_format, "tif") == 0 ||
		cfg_output_image_fgcolor)
	{
		file_name = output_path + image_name + name + "PointsXY.tif";
		XY->WriteTIFF(file_name);
		file_name = output_path + image_name + name + "PointsXZ.tif";
		XZ->WriteTIFF(file_name);
		file_name = output_path + image_name + name + "PointsYZ.tif";
		YZ->WriteTIFF(file_name);
	}
	else
	{
		file_name = output_path + image_name + name + "PointsXY.pgm";
		XY->Write(file_name);
		file_name = output_path + image_name + name + "PointsXZ.pgm";
		XZ->Write(file_name);
		file_name = output_path + image_name + name + "PointsYZ.pgm";
		YZ->Write(file_name);
	}
}

void Draw_PointsOnProjections(list<CPoint> points, char* name = NULL)
{
	//int iSlices = The3DImage->m_iSlices;
	//int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;
  /* come up with a better way to handle this case that doesn't violate
     char* = const char*
	if (name == NULL)
		name = "";
  */

	string output_path = gConfig.GetOutputPath();
	string image_name = gConfig.GetImageName();
	string file_name;
	CImage* XY = NULL;
	CPoint pPoint;

	XY = new CImage(*CanvasXY);
	XY->RemovePixels(5);

	list<CPoint>::iterator i;
	for (i = points.begin(); i != points.end(); i++)
	{
		pPoint = *i;
		if(pPoint.m_iX > iCols)
			pPoint.Print();
		else
			XY->MarkCrosshairXY(&pPoint, 2);
	}

	file_name = output_path + image_name + name + "PointsXY.tif";
	XY->WriteTIFF(file_name);
}

void Draw_SeedPointDirections()
{
	string output_path = gConfig.GetOutputPath();
	string image_name = gConfig.GetImageName();
	string file_name;
	CImage* XY = NULL;
	CPoint pPoint;

	XY = new CImage(*CanvasXY);
	XY->RemovePixels(5);

	list<CPoint>::iterator i;
	for (i = VerifiedSeedsCenter.begin(); i != VerifiedSeedsCenter.end(); i++)
	{
		pPoint = *i;
		XY->MarkCrosshairXY(&pPoint, 0);
		CPoint temp;
		for (int j = 0; j < 5; j++)
		{
			CPoint* dirPoint =& gVectorsArray[pPoint.m_iHDir][0]->
																m_pIndices[j];
			temp.m_iX = pPoint.m_iX + dirPoint->m_iX;
			temp.m_iY = pPoint.m_iY + dirPoint->m_iY;
			temp.m_iZ = pPoint.m_iZ;
			if (The3DImage->WithinImagePadding(temp, giMARGIN))
				XY->MarkPointXY(&temp, 1);
			else
				break;
		}
	}

	file_name = output_path + image_name + "SeedDirectionsXY.tif";
	XY->WriteTIFF(file_name);
	delete XY;
}

void Draw_AmriSeedPointDirections()
{
	string output_path = gConfig.GetOutputPath();
	string image_name = gConfig.GetImageName();
	string file_name;
	CImage* XY = NULL;
	CPoint pPoint;

	XY = new CImage(*CanvasXY);
	XY->RemovePixels(5);

	list<CPoint>::iterator i;
	for (i = AmriVerifiedSeedsCenter.begin();
		i != AmriVerifiedSeedsCenter.end();
		i++)
	{
		pPoint = *i;
		XY->MarkCrosshairXY(&pPoint, 0);
		CPoint temp;
		for (int j = 0; j < 5; j++)
		{
			CPoint* dirPoint =& gVectorsArray[pPoint.m_iHDir][0]->
																m_pIndices[j];
			temp.m_iX = pPoint.m_iX + dirPoint->m_iX;
			temp.m_iY = pPoint.m_iY + dirPoint->m_iY;
			temp.m_iZ = pPoint.m_iZ;
			if (The3DImage->WithinImagePadding(temp, giMARGIN))
				XY->MarkPointXY(&temp, 1);
			else
				break;
		}
	}

	file_name = output_path + image_name + "AmriSeedDirectionsXY.tif";
	XY->WriteTIFF(file_name);
}

void Draw_Binary3DImage()
{
	int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;
	string output_path = gConfig.GetOutputPath();
	string image_name = gConfig.GetImageName();
	C3DImage* newimage;
	newimage = new C3DImage(iSlices, iRows, iCols);

	int s,r,c;
	for (s = 0; s < iSlices; s++)
	{
		for (r = 0; r < iRows; r++)
		{
			for (c = 0; c < iCols; c++)
			{
				if (TracedImage[s][r][c])
					newimage->data[s][r][c] = 255;
				else
					newimage->data[s][r][c] = 0;
			}
		}
	}
	string file_name;
	file_name = output_path + image_name + "SegmentedVessels.pic";
	newimage->Write(file_name);
}

void Draw_Residual3DImage()
{
	int iPadding = The3DImage->m_iPadding;
	int iSlices = The3DImage->m_iSlices - 2 * iPadding;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;
	string output_path = gConfig.GetOutputPath();
	string image_name = gConfig.GetImageName();
	C3DImage* newimage;
	newimage = new C3DImage(iSlices, iRows, iCols);

	int s,r,c;
	int s_padded;
	for (s = 0; s < iSlices; s++)
	{
		s_padded = s + iPadding;
		for (r = 0; r < iRows; r++)
		{
			for (c = 0; c < iCols; c++)
			{
				if (!TracedImage[s_padded][r][c]) // if it is not traced, then keep it
					newimage->data[s][r][c] = The3DImage->data[s_padded][r][c];
				else
					newimage->data[s][r][c] = 0;
			}
		}
	}
	string file_name;
	file_name = output_path + image_name + "_ResidualImage.pic";
	newimage->Write(file_name);
}

void Draw_3DVessels()
{
	int iSlices = The3DImage->m_iSlices;
	int iRows = The3DImage->m_iRows;
	int iCols = The3DImage->m_iCols;
	string output_path = gConfig.GetOutputPath();
	string image_name = gConfig.GetImageName();
	C3DImage* newimage;
	newimage = new C3DImage(iSlices, iRows, iCols);

	gTheVessels.ReviseWidths();
	gTheVessels.Mark_In_3D_Image(*newimage, 255);	
	int i;

	// generate a histogram of intensity values at marked voxels

	long histogram_fg[256];
	long histogram_bg[256];
	long histogram_CR[256];		// CR = contrast reduced
	double dhistogram_fg[256];
	double dhistogram_bg[256];
	double dhistogram_CR[256];

	for(i = 0; i < 256; i++) {
		histogram_fg[i] = 0;
		histogram_bg[i] = 0;
		histogram_CR[i] = 0;
		dhistogram_fg[i] = 0;
		dhistogram_bg[i] = 0;
		dhistogram_CR[i] = 0;
	}

	long pixel_count_fg = 0;
	long pixel_count_bg = 0;
	long pixel_count_CR = 0;

	for(int s = The3DImage->m_iPadding; s < The3DImage->m_iSlices - The3DImage->m_iPadding; s++) 
	{
		for(int r = 0; r < newimage->m_iRows; r++) 
		{
			for(int c = 0; c < newimage->m_iCols; c++) 
			{
				if(newimage->data[s][r][c] == 255) {
					histogram_fg[The3DImage->data[s][r][c]]++;
					pixel_count_fg++;
				}
				else {
					histogram_bg[The3DImage->data[s][r][c]]++;
					pixel_count_bg++;
				}
				histogram_CR[The3DImage->data[s][r][c]]++;
				pixel_count_CR++;
			}
		}
	}

	for(i = 0; i < 256; i++) {
		dhistogram_fg[i] = static_cast<double>(histogram_fg[i]) / static_cast<double>(pixel_count_fg);
		dhistogram_bg[i] = static_cast<double>(histogram_bg[i]) / static_cast<double>(pixel_count_bg);
		dhistogram_CR[i] = static_cast<double>(histogram_CR[i]) / static_cast<double>(pixel_count_CR);
	}
	
	ofstream out((output_path + image_name + "_3DVessels_Histogram.txt").c_str());
	out << "i\tfg\tbg\tCR" << endl;
	for(i = 0; i < 256; i++) {
		out << i << "\t" << dhistogram_fg[i] << "\t" << dhistogram_bg[i] << "\t" << dhistogram_CR[i] << endl;
	}
	out << "Unnormalized" << endl;
	out << "i\tfg\tbg\tCR" << endl;
	for(i = 0; i < 256; i++) {
		out << i << "\t" << histogram_fg[i] << "\t" << histogram_bg[i] << "\t" << histogram_CR[i] << endl;
	}
	
	string file_name;
	file_name = output_path + image_name + "_3DVessels.pic";
	newimage->Write(file_name);
}
