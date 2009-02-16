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
