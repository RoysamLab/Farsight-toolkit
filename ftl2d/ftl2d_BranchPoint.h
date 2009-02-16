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
