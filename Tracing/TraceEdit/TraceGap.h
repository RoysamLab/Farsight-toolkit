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

#ifndef TRACEGAP_H_
#define TRACEGAP_H_


#include "TraceLine.h"
#include "TraceBit.h"

//class TraceLine;
//class TraceBit;

/**
 * A TraceGap is a simple object that represents the space between two
 * TraceLines
 **/

class VBTTracingCosts{

public:
	double tracingCost;
	double vesselnessCost;
	double scaleVar;
	double vesselnessVar;
};

class TraceGap
{
public:
  TraceGap();
  ~TraceGap();
  std::string stats();
  std::string GapStatHeaders();
	int compID;
	TraceLine *Trace1;
	TraceLine *Trace2;
	int endPT1, endPT2;
	TraceBit startBit, endBit;
	double angle;
	double dist; 
	double maxdist;
	double length;
	double smoothness;
	double cost;
	int gap_label;
	VBTTracingCosts tracingCosts;
};
#endif
