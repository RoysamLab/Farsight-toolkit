/** @file NodeContainer3D.h
*   @brief Class for Tracing in 3D volume.
*
*   @author Amit Mukherjee,
*/
#ifndef NODE_CONTAINER3D_H
#define NODE_CONTAINER3D_H

#include "itkImage.h"

#include "itkObjectFactory.h"
#include "itkSmartPointer.h"
#include "itkMacro.h"
#include "itkLightObject.h"

#include "TVessel.h"
#include "Trace.h"

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
	void WriteSegmentsToTextFile(std::string SegTxtFname);

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
