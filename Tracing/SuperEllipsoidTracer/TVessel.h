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
 \author $ Author: James Alex Tyrrell, Amit Mukherjee $
 \version $Revision 1.0 $
 \date May 05 2009
*/
// Copyright Farsight Toolkit , Rensselaer Polytechnic institute Troy NY 12180.

#ifndef TVESSEL_H
#define TVESSEL_H

#include <iostream>
#include <vector>
#include "vnl/vnl_math.h"

class TVessel
{

public:

	TVessel();
	~TVessel();

    double q3[4];
    double q2[4];
    double q1[4];
    double a1;
    double a2;
    double a3;
    double f;
    double b;
    double mu[3];
    double e1;
    double e2;
    double last_dir[3];
    double last_pos[3];
    double L;
    double bndy[24];
    double type;
    double index;
    double MAD;
    long int ID;
    long int TraceID;
    double K1;
    double K2;
    double a;
    double k;
    double R1[3];
    double R2[3];
    double R3[3];
	std::vector<long int> NbrID;
    unsigned int numNbrs;

    void PrintSelf();
    bool IsSegmentValid(TVessel*, double, double);
    TVessel* CopyToNewSegment();

};

typedef struct _tagMatrix
{
    double r[3][3];
}TRMatrix;


#endif
