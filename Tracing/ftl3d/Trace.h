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

/**
 \brief Class for the trace objects. 
 \author $ Author: Amit Mukherjee, James Alex Tyrrell $
 \version $Revision 1.0 $
 \date May 05 2009
*/
// Copyright Farsight Toolkit 2009, Rensselaer Polytechnic institute Troy NY 12180.


#ifndef TRACE_H
#define TRACE_H

#include <iostream>
#include "itkImage.h"
#include "TVessel.h"
#include "SegFit.h"
#include <vector>

typedef itk::Image<float,3> ImageType3D;

class Trace {

public:

	Trace();
	Trace(TVessel* , unsigned long );
	~Trace();

	unsigned int numNodes;
    long TraceID;
    double L;

    long NodeAID;
    long NodeBID;
	double dirA[3];
	double dirB[3];

	typedef std::vector<long> TraceIDListType;

    void PrintSelf();
    bool IsInList(TraceIDListType tlist, long id);
    void UpdateTrace(TVessel *&seg, TVessel *&seg1, const int );
    TVessel* Step(TVessel *, ImageType3D::Pointer , unsigned long , char,  double, double, double);
	void Reverse(TVessel *);

private:
	SegFit * fitter;
};


#endif
