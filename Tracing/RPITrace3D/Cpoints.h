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

#ifndef Points_H
#define Points_H

#include <cmath>

using namespace std;

//typedef enum boolean {FALSE = 0, TRUE};
enum VesselOrSoma
{
	VESSEL		= 1,
	SOMA		= 2
};
enum PointType
{
	ORIGINAL	= 1,
	MERGED		= 2,
	NEW			= 3,
	YES			= 4,
	NO			= 5
};

#define MaxNumOfVesselIDs 10

class CPoint
{
public:
	// CTOR with the default direction and pixel value
	CPoint(int a, int b, int c, unsigned char Hd = 0, unsigned char Vd = 0, int v = 0) : m_iX(a),
		m_iY(b), m_iZ(c), m_iHDir(Hd), m_iVDir(Vd), m_iValue(v),
		m_lUserFlag(0), m_iVerifiable(0), m_fVWidth(0), m_fHWidth(0), m_iID(0), m_iParID(0)
	{
	}

	// default CTOR
	CPoint() : m_iX(0), m_iY(0), m_iZ(0), m_iHDir(0), m_iVDir(0), m_iValue(0),
		m_lUserFlag(0), m_iVerifiable(0), m_iPixelValue(0), m_fVWidth(0),
    m_fHWidth(0), m_iID(0), m_iParID(0)
	{
	}

	// copy CTOR
	CPoint(const CPoint& aPoint) : m_iX(aPoint.m_iX), m_iY(aPoint.m_iY),
		m_iZ(aPoint.m_iZ), m_iHDir(aPoint.m_iHDir), m_iVDir(aPoint.m_iVDir),
		m_iValue(aPoint.m_iValue), m_lUserFlag(aPoint.m_lUserFlag),
		m_iVerifiable(aPoint.m_iVerifiable),
		m_iPixelValue(aPoint.m_iPixelValue), m_fVWidth(aPoint.m_fVWidth),
		m_fHWidth(aPoint.m_fHWidth), m_iID(aPoint.m_iID), m_iParID(aPoint.m_iParID)
	{
	}

	// DTOR (NOP)
	~CPoint()
	{
	}

	// asignment operator
	CPoint& operator=(const CPoint& rhs)
	{
		m_iX = rhs.m_iX;
		m_iY = rhs.m_iY;
		m_iZ = rhs.m_iZ;
		m_iHDir = rhs.m_iHDir;  
		m_iVDir = rhs.m_iVDir;  
		m_iValue = rhs.m_iValue;
		m_fHWidth = rhs.m_fHWidth;
		m_fVWidth = rhs.m_fVWidth;

		m_lUserFlag = rhs.m_lUserFlag;
		m_iVerifiable = rhs.m_iVerifiable;
		m_iID = rhs.m_iID;
		m_iParID = rhs.m_iParID;

		return *this;
	}
	// are the two points equal
	bool operator==(const CPoint& rhs) const
	{
		return (m_iX == rhs.m_iX && m_iY == rhs.m_iY && m_iZ == rhs.m_iZ);
	}

	// operator that always returns no. Added for the sake of
	// the linked list
	inline int operator >(CPoint&)
	{
		return 0;
	}

	bool operator >(const CPoint& right) const
	{
		return (m_iValue > right.m_iValue);
	}

	bool operator <(const CPoint& right) const
	{
		return (m_iValue < right.m_iValue);
	}


	// find the distance between the given point and myself
	inline double FindDistance(CPoint* aPoint)
	{
		if (aPoint)
		{
			double distance = (aPoint->m_iX - m_iX) * (aPoint->m_iX - m_iX) +
				   	(aPoint->m_iY - m_iY) * (aPoint->m_iY - m_iY) +
				   	(aPoint->m_iZ - m_iZ) * (aPoint->m_iZ - m_iZ);
			return std::sqrt(distance);
		}
		else
			return 0;
	}
	CPoint operator-(CPoint& rhs)
	{
		CPoint temp;
		temp.m_iX = m_iX - rhs.m_iX;
		temp.m_iY = m_iY - rhs.m_iY;
		temp.m_iZ = m_iZ - rhs.m_iZ;

		return temp;
	}
	// Print the point
	inline void Print(ostream& out = cout)
	{
		out << "Point(" << m_iX << ", " << m_iY << ", " << m_iZ << ", "
			<< (int) m_iHDir << ", " << (int) m_iVDir << ", " << m_lUserFlag
			<< ", " << (int) m_fHWidth << ", " << (int) m_fVWidth << ")= "
			<< m_iValue << endl;
	}

	friend ostream& operator <<(ostream& out, CPoint& rhs)
	{
//		out << "(" << rhs.m_iX << ", " << rhs.m_iY << ", " << rhs.m_iZ << ", "
//			<< rhs.m_lUserFlag << ")" << flush;
		rhs.Print(out);

		return out;
	}

	//added by yousef on 3-25-2009
	void SetPointID(int ID) { m_iID = ID; }
    void SetParentID(int PID) { m_iParID = PID; }
	int GetPointID() { return m_iID; }
	int GetParentID() { return m_iParID; }

	// data
	int m_iX;
	int m_iY;
	int m_iZ;
	unsigned char m_iHDir;  // the direction number (0-NumOfDirections)};
	unsigned char m_iVDir;

	int m_iValue; // this can represent the pixel value or the template's
	// response at this location

	// general purpose user flag
	long m_lUserFlag;
	int m_iVerifiable;
	
	int m_iPixelValue;

  // we keep the width of the point here
	float m_fVWidth;
	float m_fHWidth;

	//added by Yousef on 3-25-2009
	int m_iID;
	int m_iParID;

};


class CIntersectionPoint
{
public:

	// CTORs
	CIntersectionPoint();
	CIntersectionPoint(CPoint& aPoint);
	CIntersectionPoint(CIntersectionPoint& aPoint);

	// DTOR
	~CIntersectionPoint()
	{
	}

	void UpdateVesselID(int oldId, int newId, int iPosFlag);

	TopEndMiddleNeither GetLocation(int iID);

	int AddVessel(int vesselID, TopEndMiddleNeither position);
	// add the vessels of the given point to me
	void AddVessels(CIntersectionPoint*);

	void AdjustLocation();

	inline void AddSoma(int id)
	{
		m_iSomaID = id; m_Type = SOMA;
	}

	// are the two points equal
	bool operator==(CIntersectionPoint& rhs);

	int RemoveVessel(int ID);

	// is there a common vessel between this and the point at the given index.
	int IsThereACommonVessel(CIntersectionPoint*);

	// is there a path from this intersection point to the point with the give ID
	int IsThereAPath(int iID);
	//
	void Print(ostream& out = cout);

	void Init();

	// reset the processed flag
	inline void ResetFlag()
	{
		m_ProcessedFlag = NO;
	}


	//data
	CPoint m_Point;
	int m_iID;
	int m_iSomaID;
	int m_aVesselIDs[MaxNumOfVesselIDs];
	TopEndMiddleNeither m_aPosition[MaxNumOfVesselIDs];
	int m_iIndex;
	PointType m_ProcessedFlag;
	VesselOrSoma m_Type;
};

// a collection of intersection points with some specific functionality
class CIntersectionPoints
{
public:

	CIntersectionPoints()
	{
		Init();
	}
	~CIntersectionPoints();

	// better than using "#define"
	enum Constants
	{
		BlockSize = 20
	};

	int Add(CIntersectionPoint* pElement);

	int Add(CIntersectionPoint& Element);

	// Remove the intersection point with the given ID from my array
	void Remove(int pointID);

	void UpdateVesselID(int oldID, int newID, int iPosFlag);

	void CompleteVesselInfo();
	// For each two intersection points that are closer than the any of the 
	// participating segment widths, merge them into one point.
	void MergeIntersectionPoints();

	void FillDistanceArray();

	int MoreThanOneMiddleIntersection(int index1, int index2);

	// Find the closest intersection point from the point id1. 
	// Return the distance, and set id2 to the found point id.
	int FindClosestIntersectionPoint(int id1, int& id2);
	int FindClosestIntersectionPoint(CPoint& FromPoint,
		CIntersectionPoint** pClosestPoint);
	int FindClosestIntersectionPointExcluding(CPoint& FromPoint,
		CIntersectionPoint** pClosestPoint, int iID);

	// given the ids of two intersection points, id1 & id2, Find the id of the
	// vessel common to both of these intersection points.
	int GetCommonVesselID(int id1, int id2);

	int FindWidestSegment(int);
	void Print(ostream& out = cout);
	int FindLocationOfNewPoint(int, int, CPoint&);

	// return a pointer to the point with the given ID
	CIntersectionPoint* GetPoint(int id);

	// For each intersection points, ResetFlags 
	void ResetFlags();

	// remove all points that were removed because they are located between
	// two vessels that has been merged.
	void RemovePointsBetweenMergedVessels();

	int IsThereAPathToASoma(int);

	// data members
	CIntersectionPoint** m_apData;
	int m_iNumOfElements;
	int m_iSize;
	int** m_aDistanceArray;

private:
	void Init();
};

#endif
