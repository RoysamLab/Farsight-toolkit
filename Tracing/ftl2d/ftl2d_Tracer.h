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

#ifndef TRACE_H
#define TRACE_H

#include <vector>
#include "ftl2d_BranchPoint.h"
#include "ftl2d_TraceConfig.h"
#include <itkImage.h>
#include "ftl2d_Vessel.h"
#include "ftl2d_SeedContainer2D.h"
#include "ftl2d_Segment2D.h"


#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xmlreader.h>

#define MY_ENCODING "ISO-8859-1"

typedef itk::Image< float, 2 >   ImageType;


class Tracer		{

 public:
	Tracer(TraceConfig*);
	~Tracer();
	void Run(ImageType::Pointer& , SeedContainer2D *, TraceConfig *);
	Segment2D* TraceStep(ImageType::Pointer &, Segment2D* , TraceConfig *);
	void ApplyRules();
	void ReportXML (TraceConfig*);
	
	Segment2D * HitTest( Segment2D *);
	inline Vessel* getVessel(unsigned int i) {return (VesselContainer[i]);}
	inline void AddVessel(Vessel *v) {this->VesselContainer.push_back(v); }
	bool IsPointInsideSE(Vect2 &p, Segment2D *);
	bool getTraces(Vect2&, unsigned int , unsigned int );		
	void AddBranchPoint(Segment2D*, Segment2D*);
	unsigned int getNumberOfVessels(){return(VesselContainer.size());}
	
	inline std::string ToString(double val)	{
	    std::ostringstream strm;
	    strm<< val<<std::endl;
	    return strm.str();
	}

private:	
	std::vector<BranchPoint *>  BranchPointContainer;
	std::vector<Vessel *> VesselContainer;
};

#endif
