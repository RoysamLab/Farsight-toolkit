/**
 \brief A Class for storing the superellipses and also detecting their interactions. 
 \author $ Author: Amit Mukherjee, James Alex Tyrrell $
 \version $Revision 1.0 $
 \date May 05 2009
*/
// Copyright Farsight Toolkit 2009, Rensselaer Polytechnic institute Troy NY 12180.

#ifndef NODE_CONTAINER3D_H
#define NODE_CONTAINER3D_H

#include "itkImage.h"

#include "itkObjectFactory.h"
#include "itkSmartPointer.h"
#include "itkMacro.h"
#include "itkLightObject.h"

#include "TVessel.h"
#include "Trace.h"
#include "tinyxml.h"

#include <iostream>
#include <fstream>

#include <algorithm>
typedef std::vector <TVessel *> StartSegContainerType;
typedef std::vector<TVessel*> NodeContainerType;

class NodeContainer3D: public itk::LightObject {

public:
	typedef NodeContainer3D Self;
	typedef  itk::SmartPointer< Self >  Pointer;
	typedef  itk::SmartPointer< const Self >  ConstPointer;
	itkNewMacro(Self);



	TVessel* getSegment(long);
	void AddNode(TVessel*);
	void AddFront(TVessel*);

	long HitTest (TVessel *);
	long HitTest (TVessel *, Trace*);
	void PrintSelf();

   	void PrintFullTrack(Trace *);
    bool getFullTrack(NodeContainerType&, Trace*, long);
    bool getPropagationDirection(long, Trace *);
    bool getPropagationDirectionR3(long , Trace * );

    bool IsTraceValid(Trace *, TVessel*);
	void WriteSegmentsToTextFile(std::string&);
	void WriteSegmentsToXMLFile(std::string&);

private:
	NodeContainerType NodeList;
	bool IsInList(NodeContainerType&,  long);
	unsigned int NumberChainLinks(Trace *, TVessel* , TVessel* );
	double CurveLength(Trace *, TVessel* , TVessel* );
	bool ProjectedVsActualShiftOK(Trace * , TVessel* );

protected:
	NodeContainer3D(); //copy the config informations
	~NodeContainer3D();
};

#endif
