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

//////////////////////////////////////////////////////////////////////////
// File: Vessel.h
// 
// This file contains the declaration of a vessel class. 
//
// - A vessel contains three linked list of CPoint instances: left, right, and center. 
//
// - Since these three linked lists might be of different lengths, there is no
//   association information at creation time. However, such informatino can be create
//   for each centerline point, by finding the corresponding left and right points in 
//   directions perpendicular to that of the centerline points.
//
// - A vessel keeps a list (an Array) of intersection points. Each intersection point 
//   contains location and other vessel ID information. 
//
// - A vessel is usually constructed from an instant of DLList<CProfile>. After
//   construction, the constructor checks if at the vessel intersects with any 
//   other vessels at each of its end. If it does, the appropriate points are added
//   to the corresponding vessel ends, and the intersection poiint(s) is(are) noted
//   and registered in the intersection points array. 
//
// - The length of the vessel is the lenfth of its center line, and its response at
//   any centerline point is the sum of the corresponding left and right points.
//

#include <vector>

#ifndef CVessel_h
#define CVessel_h

#include <stdlib.h> /* exit */

using namespace std;

#define BlockSize 50

// forward declaration
class CImage;
class C3DImage;

class CVessel
{
public:
	// CTOR
	CVessel(int id = 0) : m_iID(id), m_iLength(0), m_iSomaConnected(0),
		m_iNumOfIntersectionPoints(0), m_iHitsImageBoundary(0),
		m_iDrawFlag(0), m_iMergedFlag(0), m_fHWidth(0.0), m_fVWidth(0.0),
		m_iNumOfPoints(0), m_iParentID(-1), m_iParentLocation(0)
	{
		m_aiMyIntersectionPoints = new int[BlockSize];
		memset(m_aiMyIntersectionPoints, 0, sizeof(int) * BlockSize);
		m_iSizeOfIntersectionPointsArray = BlockSize;
	}

	// COPY CTOR. This function was added to avoid 
	// illegal use of default copy CTOR
	CVessel(CVessel&)
	{
		cout << "CVessel::CVessel(CVessel&) ... NOT ALLOWED " << endl;
		exit(0);
	}

	// DTOR
	~CVessel()
	{
		if(m_aiMyIntersectionPoints)
			delete [] m_aiMyIntersectionPoints;
	}

	void Draw3DCenterline(C3DImage& anImage, unsigned char  color);
	void MarkVessel();
	void ReviseWidths();

	/////////////////////////////////////////////////////////////////////////////////
	// Extend the current vessel upto and including the given profile. The 
	/// Add an element to the list (always add on on the head)
	// always make your own copy
	bool ExtendVessel(CPoint*, CPoint*, CPoint*, CPoint*, CPoint*,
		TopEndMiddleNeither DirFlag);

	void MergeVessels();


	// Check the regions surrounding the given point in the specified
	// direction to see if the vessel intersects other vessels or somas
	void ConnectWithSomasAndOtherVessels();
	void ConnectWithSomasAndOtherVessels2();
	void ConnectWithSomasAndOtherVessels3();

	// Connect the vessel with adjacent somas if it lies close-by along its end direction
	int FindSomaDistance(TopEndMiddleNeither, CPoint&);
	int FindSomaDistance2(TopEndMiddleNeither, CPoint&);
	int FindSomaDistance2(TopEndMiddleNeither, CPoint&, int&);
	int FindSomaDistance3(TopEndMiddleNeither, CPoint&, int&);

	int FindIntersectionPoint(TopEndMiddleNeither, CPoint&, int&, int&);

	// return 1 if the current vessel can be extended in at the specified end
	int Extendable(TopEndMiddleNeither DirFlag);

	void ConnectWithOtherVessel(TopEndMiddleNeither DirFlag, int iDist,
		CVessel* pClosestVessel, CLNode<CPoint>* pClosestNode,
		CIntersectionPoint* p = NULL);

	void ConnectWithOtherVessels(TopEndMiddleNeither DirFlag);
	void ConnectWithOtherVessels2(TopEndMiddleNeither DirFlag);

	void ConnectWithOtherVessels3(CIntersectionPoint* p = NULL);


	int ConnectWithSomas(TopEndMiddleNeither DirFlag);
	void ConnectWithSoma(int iSomaID, CPoint&, TopEndMiddleNeither DirFlag);
	
	void Mark_In_3D_Image(C3DImage & image, unsigned char  color);

	//Draw the vessel in a 2D image
	inline void Draw(CImage& anImage)
	{
		DrawCenterline(anImage);
		DrawBoundaries(anImage);
	}
	// Draw the centerline of the vessel in the given image
	void DrawCenterline(CImage& anImage, unsigned char  color = 0);
	void DrawCenterlineXZ(CImage& anImage, unsigned char  color = 0);
	void DrawCenterlineYZ(CImage& anImage, unsigned char  color = 0);
	// Draw the boundaries of the vesse in the given Image
	void DrawBoundaries(CImage& anImage, unsigned char  color = BoundaryColor);
	void DrawBoundariesXZ(CImage& anImage, unsigned char  color = BoundaryColor);
	void DrawBoundariesYZ(CImage& anImage, unsigned char  color = BoundaryColor);
	// Draw the vessel in a 3D image
	inline void Draw(C3DImage& anImage)
	{
		DrawBoundaries(anImage);
		DrawCenterline(anImage);
	}
	// Draw in 3D
	void DrawCenterline(C3DImage& anImage, unsigned char  color = CenterlineColor);
	void DrawBoundaries(C3DImage& anImage, unsigned char  color = BoundaryColor);

	// Draw in projections
	void DrawBoundariesHXY(CImage& anImage, unsigned char  aColor);
	void DrawBoundariesHXZ(CImage& anImage, unsigned char  aColor);
	void DrawBoundariesHYZ(CImage& anImage, unsigned char  aColor);
	void DrawBoundariesVXY(CImage& anImage, unsigned char  aColor);
	void DrawBoundariesVXZ(CImage& anImage, unsigned char  aColor);
	void DrawBoundariesVYZ(CImage& anImage, unsigned char  aColor);

	// return yes if this vessel is longer than rhs
	inline int operator >(CVessel& rhs)
	{
		return (m_Center.length > rhs.m_Center.length);
	}

	void DrawThickLine(CImage& anImage, unsigned char color = 0);
	void DrawThickLine(C3DImage& anImage, unsigned char color = 0);

	void AddIntersectionPoint(int pointID);
	void ShrinkMyIntersectionPoints();

	//////////////////////
	// Replace the old intersection point ID with the new one (if any)
	void ReplaceIntersectionPoint(int oldID, int newID);

	inline int EndPoint(CPoint* aPoint)
	{
		if (*m_Center.head->data == *aPoint || *m_Center.tail->data == *aPoint)
			return 1;
		else
			return 0;
	}

	int Float(TopEndMiddleNeither);

	void WriteID(CImage& anImage, CPoint& aPoint, unsigned char color = 255);
	void WriteIDXY(CImage& anImage, unsigned char color = 255);
	void WriteIDXZ(CImage& anImage, unsigned char color = 255);
	void WriteIDYZ(CImage& anImage, unsigned char color = 255);

	int FindLocationToWriteIDXY(CPoint& atPoint);
	int FindLocationToWriteIDXZ(CPoint& atPoint);
	int FindLocationToWriteIDYZ(CPoint& atPoint);

	void ExtendByAnotherVessel(CVessel* rhs, TopEndMiddleNeither,
		TopEndMiddleNeither);

	// Mark the vessel in the given image. This is more than drawing the
	// vessels centerline and boundaries.
	void MarkInImage(C3DImage&, unsigned char color = FillingColor);

	int UpdateLength();

	////////////////////////////////////////////
	// Print vessel Info
	void Print(ostream& out = cout);

	void PrintCenterLinePoints(ostream& out = cout)
	{
		m_Center.Print(out);
	}

	void PrintAdjacency(ostream& out = cout);

	// print the vessel's centerline in Neurolucida format
	void PrintNeurolucidaFormat(ostream&, TopEndMiddleNeither);

	// break the vessel at the intersection point with the index "index"
	void BreakOnIntersectionPoint(CPoint*);

	// calculate the vessel's width at the given point. Notice that
	// if the point does not lie on the vessel, the function returns -1
	int GetWidthAtPoint(CPoint*);

	inline CPoint* GetFirstPoint()
	{
		if (m_Center.head)
			return m_Center.head->data;

		return NULL;
	}
	inline CPoint* GetLastPoint()
	{
		if (m_Center.tail)
			return m_Center.tail->data;
		return NULL;
	}

	int ExtendVesselCenter(CPoint* pPoint, TopEndMiddleNeither dirFlag);

	// determine whether the given point is on Top, on End, in the 
	// middle, or not on the vessel
	TopEndMiddleNeither GetPositionOfPoint(CPoint*);

	int IsMemberOfCenterline(CPoint& newPoint)
	{
		return m_Center.IsMember(newPoint);
	}

	// get the longest path of vessel segments including this
	int GetLongestPath(CDLList<int>& Path);

	// A vessel can have a maximum of two intersectoin points.
	// this method return the ID of the point different from the
	// argument. 
	int GetIdOfOtherIntersectionPoint(int id);

	// return the location of the floating end of this vessel
	CPoint* GetPositionOfFloatingEnd();

	// does the vessel hit the image boundar?
	inline int HitsImageBoundary()
	{
		return m_iHitsImageBoundary;
	}

	// get length
	inline int GetLength()
	{
		return m_iLength;
	}

	// set the connected to soma data member to the soma label
	inline void SetSomaConnected(int somaLabel)
	{
		m_iSomaConnected = somaLabel;
	}

	int FindClosestCenterlinePoint(CPoint& aPoint, CLNode<CPoint>** pResult);

	int FindClosestHLeftPoint(CPoint& aPoint, CLNode<CPoint>** pResult);
	int FindClosestHRightPoint(CPoint& aPoint, CLNode<CPoint>** pResult);
	int FindClosestVLeftPoint(CPoint& aPoint, CLNode<CPoint>** pResult);
	int FindClosestVRightPoint(CPoint& aPoint, CLNode<CPoint>** pResult);

	int ComputeDistanceToList(CPoint& FromPoint, CLNode<CPoint>* List,
		CLNode<CPoint>** pClosestNode);


	void RemoveIntersectionPoint(int iID);

	//added by YOusef on 3-24-2009
	//set/get the parent ID
	void SetParentID(int ID) { m_iParentID = ID; }
	int GetParentID()
	{
		return m_iParentID;
	}
	void SetParentLocation(int Loc) { m_iParentLocation = Loc; }
	int GetParentLocation()
	{
		return m_iParentLocation;
	}
	/////////////////
	// data members//
	/////////////////

	// a vessel is described by 5 linked lists of points
	// the center, horizontal left and right, and vertical left and right
	CDLList<CPoint> m_Center;
	CDLList<CPoint> m_HLeft;
	CDLList<CPoint> m_HRight;
	CDLList<CPoint> m_VLeft;
	CDLList<CPoint> m_VRight;

	// the vessel's id
	int m_iID;
	
  // vessel's length equals that of the centerline length. 
	int m_iLength;

	// an array of intersection point IDs
	int* m_aiMyIntersectionPoints;
	// the size of the allocated array
	int m_iSizeOfIntersectionPointsArray;

	int m_iSomaConnected;

	// index indicating where to store the next intersection point
	int m_iNumOfIntersectionPoints;

	int m_iHitsImageBoundary;
	int m_iDrawFlag;
	int m_iMergedFlag;

	// these are used to enable vessel extensions	
	int m_iConnectOnTopFlag;
	int m_iConnectOnEndFlag;
	CVessel* m_pTopClosestVessel;
	CVessel* m_pEndClosestVessel;
	CLNode<CPoint>* m_pTopClosestNode;	
	CLNode<CPoint>* m_pEndClosestNode;


	// added to estimate vessel width (CANCER images)
	float m_fHWidth;
	float m_fVWidth;
	int m_iNumOfPoints;

	//added by Yousef (3-24-2009): the vessel (segment) parent ID
	int m_iParentID;
	int m_iParentLocation; //default is 0 and means begining (if any), else it is at the end
	//added by Yousef on 3-25-2009
	CPoint* ParentBranchPoint;

private:

	// extend the link list of points given by aList to include all those points
	// and including aPoint. See the definition in "CVessel.cpp" for more details.
	// notice that it is made private so that it will not be called from the 
	// outside world. This method is invoked only by the "ExtendVessel" method 
	// declared above
	bool ExtendVessel(CDLList<CPoint>&, CPoint*, TopEndMiddleNeither dirFlag);

	// check one of the vessel's boundaries to see if it intersects with other vessels
	// and connect it if it does.
	void CheckAndConnectVesselEnd(CDLList<CPoint>* aList,
		TopEndMiddleNeither dirFlag);
	void CheckAndConnectVesselEnd2(CDLList<CPoint>* aList,
		TopEndMiddleNeither dirFlag);
	vector<bool> Connect();
};


// a collection of vessels. In addition to being a collection it contains
// extra functionality needed for vessel manipulations, hence we could not
// use a general purpose collection class.
class CVessels
{
public:
	CVessels() : m_iNumOfElements(0), m_iSize(BlockSize), m_apData(0),
		m_iNextID(1)
	{
		m_apData = new CVessel * [BlockSize];
	}
	~CVessels()
	{
		for (register int i = 0; i < m_iNumOfElements; i++)
		{
			if (m_apData[i])
				delete m_apData[i];
		}

		delete m_apData;
	}
	void DeleteNarrowVessels(double min_width);
	void DeleteShortNetwork(int min_length);
	void AddVessel(CVessel* aVessel);
	void Draw3DCenterline(C3DImage& anImage, unsigned char  color);
	void UpdateIntersectionPoints();
	void ExtendVessels();

	void Mark_In_3D_Image(C3DImage & image, unsigned char  color);
	void ReviseWidths();

	void MergeTwoVessels(int v1, int v2, int p1, int p2, int v3);
	void MergeVessels();
	void Print(ostream& out);
	void Print(char* FileName);

	// print all those vessels not connected to a soma (directly or indirectly)
	void PrintIsolatedVessels(char* fName);

	void PrintAdjacency(char* fName = 0);

	void DrawCenterlines(CImage& anImage)
	{
		for (register int i = 0; i < m_iNumOfElements; i++)
		{
			if (m_apData[i])
				m_apData[i]->DrawCenterline(anImage, 255);
		}
	}

	// draw the vessels and their IDs on the given image
	void DrawVessels(CImage& anImage, unsigned char  color = CenterlineColor);
	void DrawVesselsXZ(CImage& anImage, unsigned char  color = CenterlineColor);
	void DrawVesselsYZ(CImage& anImage, unsigned char  color = CenterlineColor);
	//
	void DrawVessels(C3DImage& a3DImage, unsigned char  color = 0);
	// draw the vessel boundaries
	void DrawBoundaries(CImage& anImage, unsigned char  color = BoundaryColor)
	{
		for (register int i = 0; i < m_iNumOfElements; i++)
		{
			if (m_apData[i])
				m_apData[i]->DrawBoundaries(anImage, color);
		}
	}

	// Update the lengths of all vessels to include the length of their
	// branches
	void UpdateVesselLengths();

	void PrintIntersectionPoints(char* fName = NULL); 
	CVessel* GetVessel(int id)
	{
		for (register int i = 0; i < m_iNumOfElements; i++)
		{
			if (m_apData[i])
			{
				if (m_apData[i]->m_iID == id)
					return m_apData[i];
			}
		}
		return 0;
	}
	// for each vessel, break it at its intersection points
	void BreakOnIntersectionPoints();

	void GetLongestPath();

	// does the vessel with the given ID hit the image boundary?
	inline int HitsImageBoundary(int vesselID)
	{
		if (vesselID > 0 && vesselID <= m_iNumOfElements)
			return m_apData[vesselID - 1]->HitsImageBoundary();
		return 0;
	}
	// return the length of the vessel with the given ID
	inline int GetVesselLength(int vesselID)
	{
		if (vesselID > 0 && vesselID <= m_iNumOfElements)
			return m_apData[vesselID - 1]->GetLength();
		return 0;
	}
	// set connected to soma data member of the vessel with the
	// given ID
	inline void SetSomaConnected(int vesselID, int somaLabel)
	{
		if (vesselID > 0 && vesselID <= m_iNumOfElements)
			m_apData[vesselID - 1]->SetSomaConnected(somaLabel);
	}

	int FindClosestVessel(int vesselID, CPoint* FromPoint,
		CVessel** pClosestVessel, CLNode<CPoint>** pFoundNode);

	int FindClosestVesselExcluding(int vesselID, CPoint* FromPoint,
		CVessel** pClosestVessel, CLNode<CPoint>** pFoundNode,
		CIntersectionPoint*);



	int m_iNumOfElements;
	int m_iSize;
	//private:
	// an array of vessel pointers.
	CVessel** m_apData;

	int** m_apDistanceArray;

	// the largest index in the array. 
	// Notice that this is not necessarly
	// the same as iNumOfElements. This is due to 
	// the merging operation, where two vessels are merged
	// into one.
	int m_iNextID;


private:
};
#endif
