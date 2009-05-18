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
    TVessel* Step(TVessel *, ImageType3D::Pointer , unsigned long , char);
	void Reverse(TVessel *);

private:
	SegFit * fitter;
};


#endif
