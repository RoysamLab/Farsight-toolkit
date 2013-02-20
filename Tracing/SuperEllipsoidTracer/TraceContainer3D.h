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
 \brief The Vessel class containing parameters of the superellipse and its spatial neighbors. 
 \author $ Author:  Amit Mukherjee $
 \version $Revision 1.0 $
 \date May 05 2009
*/
// Copyright Farsight Toolkit , Rensselaer Polytechnic institute Troy NY 12180.


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
	void Configure(TraceConfig::Pointer&);
	void ComputeTrace(ImageType3D::Pointer , Seed2Seg::Pointer );
	void WriteTraceToTxtFile(std::string SegTxtFname);
	void WriteTraceToXMLFile(std::string);
	void AddBranchPoint(TVessel *, TVessel *);

private:
	NodeContainer3D::Pointer NodeList;
	TraceContainer TraceList;
	double m_THRESH, m_minL, m_Stepsize, m_AspectRatio;
	itk::FixedArray<double, 3> m_Spacing;

	//TVessel* Step(TVessel *, double* , ImageType3D::Pointer , unsigned long );
	Trace* InitiazeTracer(TVessel* , unsigned long );
	//TVessel* CopyToNewSegment(TVessel * );
	//bool HitTest (TVessel *);

protected:
	TraceContainer3D(); //copy the config informations
	~TraceContainer3D();
};

#endif
