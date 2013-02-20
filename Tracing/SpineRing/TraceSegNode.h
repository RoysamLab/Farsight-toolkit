#ifndef TRACESEGNODE_H
#define TRACESEGNODE_H

//#include <iostream>
//#include <vector>
//#include <algorithm>
//#include "vnl/vnl_math.h"
//
//
//#include <itkImage.h>
//#include "itkObjectFactory.h"
//#include "itkMacro.h"
//#include "itkLightObject.h"
//
//#include "itkPointSet.h"
//
//#include "SpineConsts.h"
//#include "SpineUtils.h"
//#include "CommonTypeDefs.h"
//#include "ImageProc.h"
class ImageDebugger;


class TraceSegNode {
public:
	long ID;
	long TraceID;

	//itk::Vector<double,3> loc;
	//std::vector<long> nbrID;
	double q1[4];
	double a1;
    double a2;
    double a3;
    double f;
    double b;
    vnl_vector_fixed<double, 3> mu;
    double e1;
    double e2;
    double R1[3];
    double R2[3];
    double R3[3];
	int    type;
	std::vector<long int> NbrID;
    void PrintSelf();
   
	void rotation_quat(TRMatrix & R);
	void transpose_matrix(TRMatrix & Result );


	TraceSegNode()	{
		ID = -1;
		TraceID = -1;
		//loc.Fill(0);
	}
};

class TraceContainer {
public:
	TraceContainer (const char *f, ImageDebugger *imdebugger );
	~TraceContainer() {};
	ImageDebugger			*imdbg;
	TraceSegNodeVecType NodeContainer;
	const char* xmlfname;

	TraceIDVecType TraceIDList;
	bool ReadTraceSegNodeXMLFile();
	bool xmlvalid;
	bool xmlconsistent;
	void PrintSelf();
	PointSetPointerMap      DendMuMap;
	void					CheckConsistency(SpineImageType::Pointer image);
private:
	PointSetType::Pointer	currptset;
	PointSetType::PointsContainer::Pointer	pcont;
	PointSetType::PointDataContainer::Pointer pdata;

	PointType currpt;
};

#endif
