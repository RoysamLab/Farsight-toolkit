/** @file TraceContainer3D.h
*   @brief Class for Tracing in 3D volume.
*
*   @author Amit Mukherjee,
*/
#ifndef TRACE_CONTAINER3D_H
#define TRACE_CONTAINER3D_H

#include "itkImage.h"

#include "itkObjectFactory.h"
#include "itkSmartPointer.h"
#include "itkMacro.h"
#include "itkLightObject.h"

#include "itkImage.h"

#include "TVessel.h"
#include "Trace.h"
#include "Seed2Seg.h"
#include "SegFit.h"
#include "NodeContainer3D.h"

#include <algorithm>
#include <vector>

typedef std::vector <TVessel *> StartSegContainerType;
typedef itk::Image<float,3> ImageType3D;

class TraceContainer3D: public itk::LightObject {

public:
	typedef TraceContainer3D Self;
	typedef  itk::SmartPointer< Self >  Pointer;
	typedef  itk::SmartPointer< const Self >  ConstPointer;
	itkNewMacro(Self);

	typedef std::vector<Trace*> TraceContainer;
	void ComputeTrace(ImageType3D::Pointer , Seed2Seg::Pointer );
	void WriteTraceToTxtFile(std::string SegTxtFname);

private:
	NodeContainer3D::Pointer NodeList;
	TraceContainer TraceList;

	//TVessel* Step(TVessel *, double* , ImageType3D::Pointer , unsigned long );
	Trace* InitiazeTracer(TVessel* , unsigned long );
	//TVessel* CopyToNewSegment(TVessel * );
	//bool HitTest (TVessel *);

protected:
	TraceContainer3D(); //copy the config informations
	~TraceContainer3D();
};

#endif
