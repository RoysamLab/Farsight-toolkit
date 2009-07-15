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

#ifndef __TRACE_H
#define __TRACE_H
#include <vector>
#include <list>
#include <iostream>
#include <hash_map>
#include <queue>

#include "vtkSmartPointer.h"
#include "vtkPolyLine.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkAppendPolyData.h"

/*
#include <libxml/xmlreader.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
*/

#define MY_ENCODING "ISO-8859-1"

/* A TraceBit has the x,y,z and id of a tracelet */
struct TraceBit{
	double x,y,z,r;
	int id;
	unsigned int marker;
	void Print(std::ostream &c)
	{
		c<<"\t\tTraceBit:"<<std::endl;
		c<<"\t\tx:"<<x<<" y:"<<y<<" z:"<<z<<" r:"<<r<<" id:"<<id<<" marker:"<<marker;
	}
};

/* A TraceLine has a sequence of TraceBits and pointers to two other TraceLines */
class TraceLine{
public:
	TraceLine(){ m_parent = NULL; m_id = -(1<<30);}
	TraceLine *GetParent(){ return m_parent;}
	void SetParent(TraceLine* p){ m_parent = p;}
	void AddBranch(TraceLine* b)
	{ 
		m_branches.push_back(b);
	}
	
	TraceLine *GetBranch1(){ return m_branches[0];}
	void SetBranch1(TraceLine* b0){ m_branches[0] = b0;}
	TraceLine *GetBranch2(){ return m_branches[1];}
	void SetBranch2(TraceLine* b1){ m_branches[1] = b1;}
	unsigned char GetType(){ return m_type;}
	void SetType(unsigned char t) {m_type = t;}
	void AddTraceBit(TraceBit tbit) { m_trace_bits.push_back(tbit);}
	typedef std::list<TraceBit> TraceBitsType;
	TraceBitsType::iterator GetTraceBitIteratorBegin(){ return m_trace_bits.begin();}
	TraceBitsType::iterator GetTraceBitIteratorEnd(){ return m_trace_bits.end();}
	TraceBitsType * GetTraceBitsPointer(){ return &m_trace_bits;}
	void SetId(int lid){ m_id = lid;}
	int GetId()
	{ 
		return m_id;
	};
	int GetSize()
	{
		return m_trace_bits.size();
	};
	void Getstats();
	
	void Print(std::ostream &c)
	{
		c<<"\tTraceLine: "<<std::endl;
		c<<"\tSize: "<<m_trace_bits.size()<<std::endl;
	}
	std::vector<unsigned int> * GetMarkers() { return &m_markers;}
	std::vector<TraceLine*> * GetBranchPointer(){ return &m_branches;}

private:
	int m_id;
	std::vector<unsigned int> m_markers;
	unsigned char m_type;
	TraceLine *m_parent;
	std::vector<TraceLine* >m_branches;
	TraceBitsType m_trace_bits;
};



/* A Trace has a list of root TraceLines*/
struct TraceObject{

public:
	TraceObject(){}
	~TraceObject()
	{
		for(int counter=0; counter<trace_lines.size(); counter++)
			delete trace_lines[counter];
	}
	bool ReadFromSWCFile(char * filename);
	bool ReadFromRPIXMLFile(char * filename);
	bool ReadFromFeatureTracksFile(char *filename, int type_offset);
	bool ReadFromFeatureTracksFileForKymograph(char *filename,int type_offset);
	bool WriteToSWCFile(char * filename);
	void CreatePolyDataRecursive(TraceLine* , vtkSmartPointer<vtkFloatArray> , vtkSmartPointer<vtkPoints> ,vtkSmartPointer<vtkCellArray>);
	vtkSmartPointer<vtkPolyData> GetVTKPolyData();
	void Print(std::ostream &c)
	{
		c<<"TraceObject:"<<std::endl;
		c<<"Size:"<<trace_lines.size()<<std::endl;
		for(int counter=0; counter<trace_lines.size(); counter++)
		{
			trace_lines[counter]->Print(c);
		}
	}
	std::vector<TraceBit> CollectTraceBits();
	std::vector<TraceLine*>* GetTraceLinesPointer(){ return &trace_lines;}
	stdext::hash_map<unsigned int, unsigned long long int> hashp;
	stdext::hash_map<unsigned int, unsigned long long int> hashc;
private:
	std::vector<TraceLine*> trace_lines;	
};




#endif
