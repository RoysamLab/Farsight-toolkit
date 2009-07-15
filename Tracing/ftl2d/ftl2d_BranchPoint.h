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

/** @file Branch point.h
*   @Class that storees the results of Hit test
*
*   @author Amit Mukherjee
*/
#ifndef BRANCH_POINT_H
#define BRANCH_POINT_H

#include "ftl2d_Segment2D.h"

class BranchPoint	{
	public:
	BranchPoint() { 
	}
	~BranchPoint()	{
	}
	
	inline void AddBranchPoint(Segment2D* s1, Segment2D* s2)	{
		this->seg1 = s1;
		this->seg2 = s2;
	}
	
	private:
	Segment2D* seg1;
	Segment2D* seg2;
};

#endif
