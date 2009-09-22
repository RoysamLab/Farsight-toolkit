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

#include <iostream>
#include <math.h>
#include <vtksys/hash_map.hxx>
#include "tinyxml/tinyxml.h"

#define VTK_CREATE(type, var) \
  vtkSmartPointer<type> var = vtkSmartPointer<type>::New()

#define MY_ENCODING "ISO-8859-1"

#include "TraceBit.h"
#include "TraceLine.h"

///////////////////////////////////////////////////////////////////////////////
TraceLine::TraceLine()
{
  this->m_parent = NULL;
  this->root= -1;
  this->level = 0;
  this->m_id = -(1<<30);
  this->m_branches.clear();
}

///////////////////////////////////////////////////////////////////////////////
TraceLine::~TraceLine()
{
}

///////////////////////////////////////////////////////////////////////////////
TraceLine* TraceLine::GetParent()
{
  return this->m_parent;
}
int TraceLine::GetParentID()
{
	if (this->m_parent)
	{
		return this->m_parent->GetId();
	}
	else 
	{
		return -1;
	}
}

///////////////////////////////////////////////////////////////////////////////
void TraceLine::SetParent(TraceLine* p)
{
  this->m_parent = p;
}
int TraceLine::GetRootID()
{
	return this->root;
}
int TraceLine::GetLevel()
{
	return this->level;
}
void TraceLine::setRoot(int RootID, int traceLevel)
{
	this->root = RootID;
	this->level = traceLevel;
}
///////////////////////////////////////////////////////////////////////////////
void TraceLine::AddBranch(TraceLine* b)
{ 
  this->m_branches.push_back(b);
}


///////////////////////////////////////////////////////////////////////////////
TraceLine* TraceLine::GetBranch1()
{
  return this->m_branches[0];
}

///////////////////////////////////////////////////////////////////////////////
void TraceLine::SetBranch1(TraceLine* b0)
{
  this->m_branches[0] = b0;
}

///////////////////////////////////////////////////////////////////////////////
TraceLine* TraceLine::GetBranch2()
{
  return this->m_branches[1];
}

///////////////////////////////////////////////////////////////////////////////
void TraceLine::SetBranch2(TraceLine* b1)
{
  this->m_branches[1] = b1;
}

///////////////////////////////////////////////////////////////////////////////
unsigned char TraceLine::GetType()
{
  return this->m_type;
}

///////////////////////////////////////////////////////////////////////////////
void TraceLine::SetType(unsigned char t)
{
  this->m_type = t;
}

///////////////////////////////////////////////////////////////////////////////
void TraceLine::AddTraceBit(TraceBit tbit)
{
  this->m_trace_bits.push_back(tbit);
}

///////////////////////////////////////////////////////////////////////////////
TraceLine::TraceBitsType::iterator TraceLine::GetTraceBitIteratorBegin()
{
  return this->m_trace_bits.begin();
}

///////////////////////////////////////////////////////////////////////////////
TraceLine::TraceBitsType::iterator TraceLine::GetTraceBitIteratorEnd()
{
  return this->m_trace_bits.end();
}

///////////////////////////////////////////////////////////////////////////////
TraceLine::TraceBitsType* TraceLine::GetTraceBitsPointer()
{
  return &m_trace_bits;
}

///////////////////////////////////////////////////////////////////////////////
void TraceLine::SetId(int lid){ m_id = lid;}

///////////////////////////////////////////////////////////////////////////////
int TraceLine::GetId()
{ 
  return m_id;
};

///////////////////////////////////////////////////////////////////////////////
int TraceLine::GetSize()
{
  return m_trace_bits.size();
};

///////////////////////////////////////////////////////////////////////////////
void TraceLine::Print(std::ostream &c,int indent)
{
  for(int counter=0; counter<indent; counter++)
	  c<<" ";
  c<<"TraceLine: "<<std::endl;
  for(int counter=0; counter<indent; counter++)
	  c<<" ";
  c<<"Size: "<<m_trace_bits.size()<<std::endl;
  for(unsigned int counter=0; counter< m_branches.size(); counter++)
  {
	  m_branches[counter]->Print(std::cout,indent+4);
  }
}

///////////////////////////////////////////////////////////////////////////////
std::vector<unsigned int> * TraceLine::GetMarkers()
{
  return &(this->m_markers);
}

///////////////////////////////////////////////////////////////////////////////
std::vector<TraceLine*> * TraceLine::GetBranchPointer()
{
  return &(this->m_branches);
}

///////////////////////////////////////////////////////////////////////////////
void TraceLine::setTraceColor(double newColor)
{
  this->traceColor = newColor;
};

///////////////////////////////////////////////////////////////////////////////
double TraceLine::getTraceColor()
{
  return this->traceColor;
};

///////////////////////////////////////////////////////////////////////////////
 TraceLine::TraceLine(const TraceLine &t)
{
  this->traceColor = t.traceColor;
  this->m_id = t.m_id;
  this->m_markers = t.m_markers;
  this->m_type = t.m_type;
  this->m_parent = NULL;
  for(unsigned int counter=0; counter< t.m_branches.size(); counter++)
  {
    TraceLine *temp = new TraceLine();
    temp->SetParent(this);
    this->AddBranch(temp);
    *temp = *t.m_branches[counter];
  }
  this->m_trace_bits = t.m_trace_bits;
}

///////////////////////////////////////////////////////////////////////////////
void TraceLine::Getstats()
{
  double XF, XB, YF, YB, ZF, ZB;
  XF = m_trace_bits.front().x;
  XB = m_trace_bits.back().x;
  YF= m_trace_bits.front().y;
  YB= m_trace_bits.back().y;
  ZF= m_trace_bits.front().z;
  ZB= m_trace_bits.back().z;
  printf("Trace # %d \t Trace Size: \t %d ", m_id, (int)m_trace_bits.size());
  printf("First bit x: %4.2f y: %4.2f z: %4.2f \t", XF, YF, ZF); 
  printf("Endt bit x: %4.2f y: %4.2f z: %4.2f \n", XB, YB, ZB); 
}


///////////////////////////////////////////////////////////////////////////////
void TraceLine::EndPtDist(TraceLine *Trace2, int &dir1, int &dir2, double &dist,
                          double &maxdist, double &angle) 
{
  //int center1=this->GetSize()/2, center2=Trace2->GetSize()/2;
  double XF, XB, YF, YB, ZF, ZB, XF2,XB2, YF2, YB2, ZF2, ZB2, distances[4];
  double delX1, delX2, delY1, delY2, delZ1, delZ2, min, norm1, norm2;

  //trace 1
  XF = m_trace_bits.front().x;
  YF= m_trace_bits.front().y;
  ZF= m_trace_bits.front().z;
  XB = m_trace_bits.back().x;
  YB= m_trace_bits.back().y;
  ZB= m_trace_bits.back().z;

  //trace 2 
  XF2=Trace2->m_trace_bits.front().x;
  YF2=Trace2->m_trace_bits.front().y;
  ZF2=Trace2->m_trace_bits.front().z;
  XB2=Trace2->m_trace_bits.back().x;
  YB2=Trace2->m_trace_bits.back().y;
  ZB2=Trace2->m_trace_bits.back().z;

  //delta x,y,z 
  delX1=XF-XB;
  delX2=XF2-XB2;
  delY1=YF-YB;
  delY2=YF2-YB2;
  delZ1=ZF-ZB;
  delZ2=ZF2-ZB2;

  //compute the endpt distances
  distances[0]=sqrt(pow((XF-XF2),2)+pow((YF-YF2),2)+pow((ZF-ZF2),2));//0 F-F
  distances[1]=sqrt(pow((XF-XB2),2)+pow((YF-YB2),2)+pow((ZF-ZB2),2));//1 F-B
  distances[2]=sqrt(pow((XB-XF2),2)+pow((YB-YF2),2)+pow((ZB-ZF2),2));//2 B-F
  distances[3]=sqrt(pow((XB-XB2),2)+pow((YB-YB2),2)+pow((ZB-ZB2),2));//3 B-B

  //determine minimum spacing
  min = distances[0];
  int i, mark=0;
  for (i = 1; i<4; i++)
    {
    if (min > distances[i])
      {
      min= distances[i];
      mark=i;
      }
    } 
  norm1=sqrt(pow((delX1),2)+ pow((delY1),2)+ pow((delZ1),2));
  norm2=sqrt(pow((delX2),2)+ pow((delY2),2)+ pow((delZ2),2));
  angle = acos(((delX1 * delX2) + (delY1 *delY2) + (delZ1 * delZ2)) /
          (norm2 * norm1));
  // from min determine orientation and return distance
  if (mark ==0)
    {
    dist = distances[0];
    dir1= m_trace_bits.front().marker;
    dir2= Trace2->m_trace_bits.front().marker;
    angle=PI-angle;
    maxdist= distances[3];
    }
  else if (mark ==1)
    {
    dist = distances[1];
    dir1= m_trace_bits.front().marker;
    dir2= Trace2->m_trace_bits.back().marker;
    maxdist= distances[2];
    }
  else if (mark ==2)
    {
    dist = distances[2];
    dir1= m_trace_bits.back().marker;
    dir2= Trace2->m_trace_bits.front().marker;
    maxdist= distances[1];
    }
  else
    {
    dist = distances[3];
    dir1= m_trace_bits.back().marker;
    dir2= Trace2->m_trace_bits.back().marker;
    angle=PI-angle;
    maxdist= distances[0];
    }
}
bool TraceLine::Orient(TraceLine * Trunk)
{
	double XF, XB, YF, YB, ZF, ZB, XB2, YB2, ZB2, distances[2];
	//trace 1
	XF = m_trace_bits.front().x;
	YF= m_trace_bits.front().y;
	ZF= m_trace_bits.front().z;
	XB = m_trace_bits.back().x;
	YB= m_trace_bits.back().y;
	ZB= m_trace_bits.back().z;

	//Trunk 	
	XB2=Trunk->m_trace_bits.back().x;
	YB2=Trunk->m_trace_bits.back().y;
	ZB2=Trunk->m_trace_bits.back().z;
	//compute the endpt distances
	distances[0]=sqrt(pow((XF-XB2),2)+pow((YF-YB2),2)+pow((ZF-ZB2),2));// F-B
	distances[1]=sqrt(pow((XB-XB2),2)+pow((YB-YB2),2)+pow((ZB-ZB2),2));// B-B
	if(distances[0] < distances[1])
	{
		return true;	//oriented correctly
	}
	return false;		//needs to be flipped
}

///////////////////////////////////////////////////////////////////////////////
std::vector<double> TraceLine::stats()
{
  std::vector<double> thisStats;
  thisStats.push_back(this->m_id);
  thisStats.push_back(this->m_type);
  thisStats.push_back(this->GetSize());
  //thisStats.push_back(this->m_parent->GetId());
  thisStats.push_back(this->m_trace_bits.front().x);
  thisStats.push_back(this->m_trace_bits.front().y);
  thisStats.push_back(this->m_trace_bits.front().z);
  thisStats.push_back(this->m_trace_bits.back().x);
  thisStats.push_back(this->m_trace_bits.back().y);
  thisStats.push_back(this->m_trace_bits.back().z);
  return thisStats;
}

