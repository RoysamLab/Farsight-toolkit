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
 \brief The Class for fitting the initial superellipse from Seed Points in 3D volume. 
 \author $ Author: Amit Mukherjee, James Alex Tyrrell $
 \version $Revision 1.0 $
 \date May 05 2009
*/
// Copyright Farsight Toolkit 2009, Rensselaer Polytechnic institute Troy NY 12180.


#ifndef SEED2SEG_H_
#define SEED2SEG_H_

#include "itkObjectFactory.h"
#include "itkSmartPointer.h"
#include "itkMacro.h"
#include "itkLightObject.h"

#include "itkImage.h"


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <vector>
#include <algorithm>
#include <cstdlib>
#include <vnl/vnl_vector_fixed.h>

#include "TVessel.h"
#include "SegInit.h"
#include "SeedContainer3D.h"
#include "TraceConfig.h"

typedef itk::Image<float , 3> ImageType3D;
typedef std::vector<SeedPoint3D *> SeedContainerType;
typedef std::vector <TVessel *> StartSegContainerType;

class Seed2Seg: public itk::LightObject	{
public:
	typedef Seed2Seg Self;
	typedef  itk::SmartPointer< Self >  Pointer;
	typedef  itk::SmartPointer< const Self >  ConstPointer;
	itkNewMacro(Self);
	typedef vnl_vector_fixed<double,3> Vect3;

	void Configure(double FitIterations, double AspectRatio, 
		double minVesselWidth, double startThresh, double minContrast);
	void ComuputeStartSegments(SeedContainer3D::Pointer , ImageType3D::Pointer);
	void SortStartSegments();
	unsigned int getNumberOfStartSegments() {return (unsigned int)SSContainer.size();}
	TVessel* getStartSegment(unsigned int i) {return SSContainer[i];}


private:
	bool HitTestSeedPoints(SeedPoint3D*);
	bool HitTestSeg(TVessel*  );
	bool IsStartSegmentValid(TVessel* );
	//bool Sortfunction(TVessel*s1, TVessel*s2);

	double iterations;
	double AS_RATIO;
	StartSegContainerType SSContainer;
	double min_a, THRESH, minL;



protected:
	Seed2Seg();
	Seed2Seg(TraceConfig::Pointer& );

	~Seed2Seg();
};

#endif /*SEGINIT_H_*/
