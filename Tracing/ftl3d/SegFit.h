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
 \brief Class to update parameters of superellipse from its spatial neighbor. 
 \author $ Author: James Alex Tyrrell, Amit Mukherjee $
 \version $Revision 1.0 $
 \date May 05 2009
*/
// Copyright Farsight Toolkit , Rensselaer Polytechnic institute Troy NY 12180.

#ifndef SEGFIT_H_
#define SEGFIT_H_

#include "itkObjectFactory.h"
#include "itkSmartPointer.h"
#include "itkMacro.h"
#include "itkLightObject.h"

#include "itkImage.h"
#include "itkLinearInterpolateImageFunction.h"
#include "vnl/vnl_math.h"


#include <stdlib.h>
#include <algorithm>
#include <cstdlib>
#include <stdio.h>
#include <string.h>
#include <math.h>


#include "TVessel.h"
#include "SegInit.h"

typedef itk::Image<float , 3> ImageType3D;



#define DIM 3     /* dimension of points, must be < 31 for SIZEcube */
#define SIZEcube (1<<DIM)
#define SIZEdiamond (2*DIM)
#define TOTpoints (SIZEcube + SIZEdiamond)
#define _PI 3.14159265
#define sign(x) (x>0?1.0:-1.0)
#define _INC 19
#define MAX_SIZE 2000000
#define XR(i,j) xr[i+4*j]
#define XI(i,j) xi[i+4*j]



/*class SegFit: public itk::LightObject	{
public:
	typedef SegFit Self;
	typedef  itk::SmartPointer< Self >  Pointer;
	typedef  itk::SmartPointer< const Self >  ConstPointer;
	itkNewMacro(Self);

	*/
class SegFit : public SegInit	{
public:

bool fitSE (ImageType3D::Pointer, TVessel&, double , double , double);

public:
	SegFit()	{
		gNum_facets = 0;
		gNum_points = 0;
	}

	~SegFit()	{
	}

};

#endif /*SEGINIT_H_*/
