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

#ifndef VESSEL_H
#define VESSEL_H

#include <vector>
#include "ftl2d_Segment2D.h"
#include "ftl2d_TraceConfig.h"

class Vessel 	{
	
public:
	Vessel(TraceConfig*, unsigned int);
	~Vessel();

	void AddSegmentToVessel( Segment2D* );
	bool AddSegmentToVesselIfValid(Segment2D*);	
	inline Segment2D* getSegment(unsigned int i) {return (SegmentContainer[i]);}
	unsigned int getNumberOfSegments() {SegmentContainer.size();}

	Segment2D* TraceStep(ImageType::Pointer, Segment2D *);
	void flipDirection();
	bool IsVesselValid();
	bool IsPointInsideSE(Vect2&, Segment2D*);
	inline unsigned int getID() {return(VesselID);}
	inline void SetEndType(unsigned int e) {this->EndType = e;}
	void PrintSelf();

private:
	unsigned int VesselID;
	unsigned int SeedPointID;
	double VesselLength;
	unsigned int SegmentCounter;
	unsigned int self_intersect;
	double MinVesselLength;
	unsigned int EndType;
	double last_predicted_steps[5];  
	double last_traversed_steps[5];
	
	std::vector<Segment2D *> SegmentContainer;

};

#endif
