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

/////////////////////////////////////////////////////////////////////////////
// File: Profile.h
//
// declaration file for a Class representing the line a cross a blood vessel. 
// The profile is defined by its left and right boundary points, and its 
// center point. In addition, the profile contains information about the 
// left and right template responses.
//
#ifndef profile_h
#define profile_h

/////////////////////////////////////////////////////////////////////
// Class: CProfile
// a class that represents a line across a blood vessel as described
// by three points; center, left, and right points. It also contains
// information regarding the response of the left and right templates
// at the specified points.
class CProfile
{
public:
	// default CTOR
	CProfile()
		: m_iWidth(0), m_iRespL(0), m_iRespR(0)
	{
	}

	// copy CTOR
	CProfile(CProfile &aProfile)
		: m_Center(aProfile.m_Center), m_Left(aProfile.m_Left),
		  m_Right(aProfile.m_Right), m_iWidth(aProfile.m_iWidth),
		  m_iRespL(aProfile.m_iRespL), m_iRespR(aProfile.m_iRespR)
	{
	}

	// CTOR
	CProfile(CPoint *c, CPoint *l, CPoint *r)
	{
		m_Center = *c;
		m_Left   = *l;
		m_Right  = *r;

		int xdiff = m_Left.m_iX - m_Right.m_iX;
		int ydiff = m_Left.m_iY - m_Right.m_iY;		
		m_iWidth = (int) (sqrt( (double)(xdiff*xdiff + ydiff*ydiff)) + 0.5);
	}
	// DTOR, NOP
	~CProfile(){}

		// operator that always returns no added for the sake of
	// the linked list
	inline int operator >(CProfile &rhs) { return 0; }

	// Print
	void Print(ostream &out)
	{
		out << "\n Profile: ";
		out << "Left: "		<< m_Left;
		out << ", Center: "	<< m_Center;
		out << ", Right: "	<< m_Right;

		out << " width: " << m_iWidth 
			 << ", resp(L, R): " << m_iRespL 
			 << ", " << m_iRespR 
			 << endl;
	}

	// data
	CPoint m_Center;
	CPoint m_Left;
	CPoint m_Right;
	int m_iWidth;
	int m_iRespL;
	int m_iRespR;
};

#endif
