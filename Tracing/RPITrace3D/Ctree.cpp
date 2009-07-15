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

///////////////////////////////////////////////////////////////////////
// File: CTree.cpp
// By:   Khalid Al-Kofahi
// Created: 2-12-99
//
// This file contains the declaration of a template tree class and the required
// supporting classes

//#pragma warning(disable:4786)  // disable STL-related warnings
//#pragma warning(disable:4710)  // compiler didn't inline function that got selected for inline expansion
//#pragma warning(disable:4514)  // The optimizer removed an inline function that is not called


#include <iostream>
#include <fstream>
#include <list>
#include <queue>
#include <iomanip>
#include <string>
#include <cmath>

#include "Mytypes.h"	// all global variables
#include "Config.h"
#include "Cpoints.h"
#include "CONSTANTS.h"
#include "Dllist.h"
#include "Cimage.h"    // a simple image class
#include "Ctree.h"
#include "Soma.h"
#include "Cvessel.h"
#include "Cvector.h"
#include "Template.h"  // a template class
#include "Extern.h"

using namespace std;

extern CSomas* gTheSomas;
int giFlag = 1;

//////////////////////////////////////////////////////////////////////////////
//									Class CTreeNode											 //
//////////////////////////////////////////////////////////////////////////////

//////////////////
// CTOR.
// A node represents an intersection point. the argument "SegmentID" denotes
// the id of the segment leading to this node from the parent node
CTreeNode::CTreeNode(int ParentSomaID, int SegmentID,
	CIntersectionPoint* aIntPoint)
{
	CVessel* pVessel = NULL;

	// initialization
	Init();

	// the ID ot the segment to this node
	m_iSegmentID = SegmentID;

	if (SegmentID)
	{
		// if the segment hits the image boundary, record it
		if (gTheVessels.HitsImageBoundary(SegmentID))
			gTheSomas->SomaHitsImageBoundary(ParentSomaID);

		// record that the segment is connected to this soma
		gTheVessels.SetSomaConnected(SegmentID, ParentSomaID);
	}

	// if there is no segment coming to this node, then I am a root node
	if (m_iSegmentID == 0)
		m_Type = ROOT;

	// if the point is hitting another soma, create a leaf node and stop
	else if (SegmentID && aIntPoint && aIntPoint->m_Type == SOMA)
	{
		m_Type = LEAF;
		m_Point = aIntPoint->m_Point;

		// inform both this soma and my parent soma
		// that they hit each other
		gTheSomas->ConnectSomas(aIntPoint->m_iSomaID, ParentSomaID);

		// set the flag to avoid cycling
		aIntPoint->m_ProcessedFlag = YES;

		return;
	}
	// if there is no intersection point at this node, then I am a leaf node
	else if (aIntPoint == NULL)
	{
		 	m_Type = LEAF;
		 	// still we need to determine the location of this node
		 	pVessel = gTheVessels.GetVessel(SegmentID);
		CPoint* temp = pVessel->GetPositionOfFloatingEnd();


		 	//CPoint *temp = gTheVessels.m_apData[SegmentID - 1]->GetPositionOfFloatingEnd();
		if (temp)
			m_Point = *temp;
		 	//	else
		 	//		cout << "\nCTreeNode::CTreeNode. Error1" << endl;

		 	return;
	}
		// if the given intersection point, has already been processed, then we have
		// a cycle in this tree. Create a leaf node and return
	else if (aIntPoint->m_ProcessedFlag == YES)
	{
		m_Type = LEAF;
		m_Point = aIntPoint->m_Point;
		return;
	}

		// fill some of my data members
	m_iNodeID = aIntPoint->m_iID;
	// the location of this node
	m_Point = aIntPoint->m_Point;

	int ChildID = 0;
	int IntPoint2ID = 0;
	CIntersectionPoint* pIntPoint2 = NULL;

	aIntPoint->m_ProcessedFlag = YES;

	// for each of this intersection point segment, create a new child node
	for (register int i = 0; i < aIntPoint->m_iIndex; i++)
	{
		// the id of the child segment
		ChildID = aIntPoint->m_aVesselIDs[i];

		// if the segment already exists in my tree, ignore it
		if (! gTheSomas->AddSegment(ParentSomaID, ChildID))
			continue;

		// get the second intersection point adjacent to this vessel. Notice that
		// each segment is bounded by at most two intersection points, one at each
		// end. One of such points is the argument to this
		pVessel = gTheVessels.GetVessel(ChildID);
		IntPoint2ID = pVessel->GetIdOfOtherIntersectionPoint(m_iNodeID);

		// get a pointer to the other intersection point
		pIntPoint2 = gIntersectionPoints.GetPoint(IntPoint2ID);

		// Create a child node for such segment
		// Notice that we don't test to see if such intersection point exists.
		// This because if pIntPoint2 is NULL, then we will create a leaf node
		m_aChildren[m_iNumOfChildren] = new CTreeNode(ParentSomaID,
												ChildID,
												pIntPoint2);

		// the length sum of my branches, is given by the length sum of each
		// of my children nodes, and the segments going to them.
		m_iSumOfAllBranches += gTheVessels.GetVesselLength(ChildID);
		m_iSumOfAllBranches += m_aChildren[m_iNumOfChildren]->m_iSumOfAllBranches;

		m_iNumOfChildren++;

		if (m_iNumOfChildren >= MaxNumOfChildren)
			break;
	}
}

// DTOR
CTreeNode::~CTreeNode()
{
	for (register int i = 0; i < m_iNumOfChildren; i++)
	{
		if (m_aChildren[i])
			delete m_aChildren[i];
	}
}

///////////////////////////
// Method: GetLongestPath
//
// determine the longest path from this node on..
int CTreeNode::FindLongestPath(CDLList<int>& CurrentPath)
{
	if (m_Type == ROOT)
		return m_aChildren[0]->FindLongestPath(CurrentPath);

	// add the current segment to the path, and update its length
	// accordingly
	CurrentPath.AddElementOnEnd(&m_iSegmentID);
	CurrentPath.m_lUserFlag += gTheVessels.GetVesselLength(m_iSegmentID);

	// working paths up to the current path
	CDLList<int> tempPath(CurrentPath);
	CDLList<int> MaxPath;

	int tempLength = 0;
	int MaxLength = CurrentPath.m_lUserFlag;

	// for each of my children
	for (register int i = 0; i < m_iNumOfChildren; i++)
	{
		tempLength = m_aChildren[i]->FindLongestPath(tempPath);
		if (tempLength > MaxLength)
		{
			MaxLength = tempLength;
			MaxPath = tempPath;
		}
		// reset the temporary path to the current path
		tempPath = CurrentPath;
	}

	// if the maximum path is longer than the current path, update it
	if (MaxLength > CurrentPath.m_lUserFlag)
	{
		CurrentPath = MaxPath;
		return MaxLength;
	}

	return CurrentPath.m_lUserFlag;
}

/////////////////
// Method: Init
// 
// Initialization. It is declared private to limit access
void CTreeNode::Init()
{
	m_iSegmentID = 0;
	m_iNodeID = 0;
	m_iProcessedFlag = 0;
	m_iNumOfChildren = 0;
	m_iSumOfAllBranches = 0;
	m_Type = NORMAL;
	memset(m_aChildren, 0, sizeof(CTreeNode *) * MaxNumOfChildren);
}

/////////////////////
// Method: Print
//
// Print the contents of the tree
// for some reason the stream manipulation operators setw, etc does not work
// when displaying results in notepad. so use space padding instead
void CTreeNode::Print(ostream& outFile)
{
	if (m_Type != ROOT)
	{
		CVessel* aVessel = gTheVessels.GetVessel(m_iSegmentID);
		if (aVessel)
		{
			CPoint* aPoint1 = aVessel->GetFirstPoint();
			CPoint* aPoint2 = aVessel->GetLastPoint();

			if (*aPoint1 == m_Point)
			{
				outFile << "\t\t" << m_iSegmentID << "\t\t\t(" 
					<< aPoint1->m_iX << ", " << aPoint1->m_iY << ", " << aPoint1->m_iZ << ")-(" 
					<< aPoint2->m_iX << ", " << aPoint2->m_iY << ", " << aPoint2->m_iZ << ")";
			}
			else
			{
				outFile << "\t\t" << m_iSegmentID << "\t\t\t(" 
					<< aPoint1->m_iX << ", " << aPoint1->m_iY << ", " << aPoint1->m_iZ << ")-(" 
					<< aPoint2->m_iX << ", " << aPoint2->m_iY << ", " << aPoint2->m_iZ << ")";
			}
			outFile << "\t\t" << gTheVessels.GetVesselLength(m_iSegmentID) << "\t\t\t"
				<< m_iSumOfAllBranches << "\n";
		}
	}

	for (register int i = 0; i < m_iNumOfChildren; i++)
		m_aChildren[i]->Print(outFile);
}

//////////////////////////////////////////////////
// Method: Draw
// 
// draw the tree in the projection images
//
void CTreeNode::Draw(unsigned char color)
{
	if (m_Type != ROOT)
	{
		CVessel* aVessel = gTheVessels.GetVessel(m_iSegmentID);
		if (aVessel)
		{
			aVessel->DrawCenterline(*TrackImageXY, color); 
			aVessel->DrawCenterlineXZ(*TrackImageXZ, color);
			aVessel->DrawCenterlineYZ(*TrackImageYZ, color);

			aVessel->WriteIDXY(*TrackImageXY, IDColor);
			aVessel->WriteIDXZ(*TrackImageXZ, IDColor);
			aVessel->WriteIDYZ(*TrackImageYZ, IDColor);

			aVessel->m_iDrawFlag = 1;
		}
	}

	for (register int i = 0; i < m_iNumOfChildren; i++)
		m_aChildren[i]->Draw(color);
}

//////////////////////////////////////////////////////////////////////////
// Method: PrintNeurolucidaFormat
//
// Write my subtree to the given file in Neurolucida format using the
// given color
void CTreeNode::PrintNeurolucidaFormat(ofstream& outFile, int iColor)
{
	if (m_Type == ROOT)
		outFile << " (Axon)" << endl;
	else
	{
		CVessel* pVessel = gTheVessels.GetVessel(m_iSegmentID);
		int iFlippedFlag = 0;
		if (pVessel)
		{
			CLNode<CPoint>* head = pVessel->m_Center.head;
			CLNode<CPoint>* tail = pVessel->m_Center.tail;
			CLNode<CPoint>* temp = tail;

			if (! (*(temp->data) == m_Point))
			{
				head = pVessel->m_Center.tail;
				tail = pVessel->m_Center.head;
				iFlippedFlag = 1;
			}
			temp = head;
			while (temp)
			{
				// print the leading spaces
				if (giFlag)
					outFile << gachLeadingSpaces;
				giFlag = 1;

				outFile << "( " << temp->data->m_iX << " " << temp->data->m_iY
					<< " " << temp->data->m_iZ << " )" << endl;

				if (iFlippedFlag)
					temp = temp->before;
				else
					temp = temp->after;
			}
		}
	}

	if (m_iNumOfChildren)
	{
		if (m_Type == ROOT)
		{
			giNumOfLeadingSpaces++;
			memcpy(gachLeadingSpaces,
				gachSpaces,
				sizeof(char) * giNumOfLeadingSpaces);
			gachLeadingSpaces[giNumOfLeadingSpaces] = '\0';
		}
		else
		{
			// becuase I must have more than one child, I must start an object
			if (m_iNumOfChildren >= 2)
			{
				outFile << gachLeadingSpaces << "(";
				giNumOfLeadingSpaces++;
				memcpy(gachLeadingSpaces,
					gachSpaces,
					sizeof(char) * giNumOfLeadingSpaces);
				gachLeadingSpaces[giNumOfLeadingSpaces] = '\0';	
				giFlag = 0;
			}
		}

		for (register int i = 0; i < m_iNumOfChildren; i++)
		{
			m_aChildren[i]->PrintNeurolucidaFormat(outFile, iColor);

			if (m_Type != ROOT && i != m_iNumOfChildren - 1)
				outFile << "|" << endl;
		}

		giNumOfLeadingSpaces--;
		if (giNumOfLeadingSpaces < 0)
			giNumOfLeadingSpaces = 0;

		gachLeadingSpaces[giNumOfLeadingSpaces] = '\0';
		if (m_Type != ROOT && m_iNumOfChildren != 1)
			outFile << gachLeadingSpaces << "); End SubTree" << endl;
	}
}

//////////////////////////////////////////////////////////////////////////////
//									Class CTree    											 //
//////////////////////////////////////////////////////////////////////////////

// default CTOR
CTree::CTree() : root(0), m_iSumOfAllBranches(0), m_iLengthOfLongestPath(0)
{
}

// DTOR
CTree::~CTree()
{
	if (root)
		delete root;
}

///////////////////////////////////
// Method: ConstructTree
//
// Construct a tree from the root intersection point. Such intersection
// point will only have one vessel and one soma data members
void CTree::ConstructTree(int ParentSomaID, CIntersectionPoint* pIntPoint)
{
	// the segment (dendrite) ID
	int segmentID = pIntPoint->m_aVesselIDs[0];
	if (segmentID)
	{
		// to prevent cycling, a flag is raised in each of the intersection
		// points upon processing it. Reset such flags
		gIntersectionPoints.ResetFlags();

		// create the root. The argument 0 represent the ID of the segment
		// connecting to the new node. But since such node is the 
		// root no such segment exists
		root = new CTreeNode(ParentSomaID, 0, pIntPoint);
		if (root->m_aChildren[0] != NULL)
		{
			m_iSumOfAllBranches = root->m_iSumOfAllBranches;

			// find the longest path in this tree
			FindLongestPath();
		}
		else
		{
			delete root;
			root = NULL;
		}
	}
}


//////////////////////////////////
// Method: FindLongestPath
//
// Find the longest path in this tree and return its length
int CTree::FindLongestPath()
{
	m_iLengthOfLongestPath = root->FindLongestPath(m_LongestPath);
	return m_iLengthOfLongestPath;
}

/////////////////////////////////
// Method: Print
//
// Print the contents of the tree
void CTree::Print(ostream& outFile)
{
	// 
	outFile << " (" << root->m_Point.m_iX << ", " << root->m_Point.m_iY << ", "
		<< root->m_Point.m_iZ << ")\n";

	outFile << "\t===================\n";
	outFile << "\tSum Of Branches Lengths: " << m_iSumOfAllBranches << "\n";
	outFile << "\tLongest Path:";
	CLNode<int>* temp = m_LongestPath.head;
	while (temp)
	{
		outFile << * (temp->data) << ", ";
		temp = temp->after;
	}

	outFile << "\n\tLength Of Longest Path: " << m_iLengthOfLongestPath << "\n";
	outFile << "\tLength Of Longest Path Branches: " << m_iSumOfAllBranches - m_iLengthOfLongestPath << "\n\n";
	outFile << "    Segment ID              From-To            Length    Branches Sum\n";
	outFile << "    =================================================================\n";

	root->Print(outFile);

	outFile << "\n\n";
}

///////////////////////////////////////////////////////////
// Method: Draw its vessels in the three projection images
//
// Print the contents of the tree
void CTree::Draw(unsigned char color)
{
	root->Draw(color);
}

//////////////////////////////////////////////////////////////////////////
// Method: PrintNeurolucidaFormat
//
// Write this tree to the given file in Neurolucida format using the
// given color
void CTree::PrintNeurolucidaFormat(ofstream& outFile, int iColor)
{
	if (root)
	{
		// increment the number of spaces to be printed
		giNumOfLeadingSpaces = 0;
		gachLeadingSpaces[giNumOfLeadingSpaces] = '\0';

		// start a tree
		outFile << "((Color RGB(" << gaRGBColors[iColor][0] << " "
			<< gaRGBColors[iColor][1] << " " << gaRGBColors[iColor][2] << "))"
			<< endl;			   

		root->PrintNeurolucidaFormat(outFile, iColor);
		giNumOfLeadingSpaces = 0;
		outFile << "); End Of Tree " << endl;
	}
}



