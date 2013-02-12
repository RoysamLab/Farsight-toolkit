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

#include "TraceGap.h"

TraceGap::TraceGap()
{
}

TraceGap::~TraceGap()
{
}
std::string TraceGap::stats()
{
	std::stringstream thisStats;
	thisStats << this->compID << "\t" << this->Trace1->GetId()
		<< "\t" << this->Trace2->GetId() << "\t" << this->angle  
		<< "\t" << this->dist << "\t" << this->maxdist
		<< "\t" << this->length << "\t" << this->smoothness	<<"\t" << this->cost
		<< "\t" << this->tracingCosts.tracingCost;
	return thisStats.str();
}
std::string TraceGap::GapStatHeaders()
{
	return std::string( "Gap:\t Trace:\t to:\t Angle:\t Distance:\t Maximum Distance:\t Length:\t Smoothness:\t cost:\t TracingCost\t");
}
