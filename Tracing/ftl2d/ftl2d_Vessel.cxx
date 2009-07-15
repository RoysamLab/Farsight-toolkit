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

/** @file SEtracing2D.cpp
*   @brief The main program of the tracing algorithm using Superellipsoid.
*
*   @author Amit Mukherjee
*//////////////////////////////////////////////////////////////////////////////////

#include "ftl2d_Vessel.h"

Vessel::Vessel(TraceConfig* config, unsigned int ID)	{
	this->VesselID = ID;
	this->MinVesselLength = config->getMinimumVesselLength();
	this->self_intersect = 0;
	this->VesselLength = 0;
	for(int i=0;i<5;i++)	{
			this->last_predicted_steps[i] = 0;
			this->last_traversed_steps[i] = 0;
	}

	this->SegmentContainer.reserve(100);
	this->SeedPointID = 0;
	this->EndType = 0;
}

Vessel::~Vessel()	{
	//delete SegmentContainer
	std::vector<Segment2D *>::iterator it;
	for (it = SegmentContainer.begin(); it != SegmentContainer.end(); ++it)
	         delete (*it);
	// std::cout << "Vessel DELETED" << std::endl;
}


void Vessel::flipDirection()	{

//	std::vector<Segment2D*>::iterator it = this->SegmentContainer.begin();
	this->SeedPointID = this->SegmentContainer.back()->getID();
	unsigned int m_size = this->SegmentContainer.size()-1;
	Segment2D * temp;
	for (unsigned int i = 0, j = m_size; i < m_size/2; i++, j--)	{
		temp = this->SegmentContainer[i];
		this->SegmentContainer[i] = this->SegmentContainer[j];
		this->SegmentContainer[j] = temp;
	}


	for(int i=0;i<this->SegmentContainer.size();i++)	{
		this->SegmentContainer[i]->setID(i);
	}

	for(int i=0;i<5;i++)	{
		this->last_predicted_steps[i] = 0;
		this->last_traversed_steps[i] = 0;
	}

	this->self_intersect = 0;

}

void Vessel::AddSegmentToVessel( Segment2D* s) {
	s->setVesselID(this->VesselID);
	SegmentContainer.push_back(s);
}


bool Vessel::AddSegmentToVesselIfValid(Segment2D* s)	{
	if (s->getEndType() != 0){
		this->EndType = s->getEndType();
		this->SegmentContainer.back()->setEndType(s->getEndType());
		return 0;			//If segment level rules fail then return 0
	}

	s->setVesselID(this->VesselID);

	//Implement the hittest here
	std::vector<Segment2D*>::iterator it;
	unsigned int minimumIndex = 10000;
	double signR2 = 1.0;
	Vect2 mu, bdry1, bdry2, offset;

/*	offset[0] =  sin(s->getq())*s->geta1();  //Porjection on R2
	offset[1] =  cos(s->getq())*s->geta1();
	mu = s->getmu();
	bdry1 = mu + offset;
	bdry2 = mu - offset;


	for(it = this->SegmentContainer.begin(); it !=this->SegmentContainer.end(); ++it) {
		if (IsPointInsideSE(bdry1, *it) || IsPointInsideSE(bdry2, *it) || IsPointInsideSE(mu, *it))	{
			//obtain the minimum index which is the farthest segment
			if((*it)->getID() < minimumIndex) 	{
				minimumIndex = (*it)->getID();
				signR2 = dot_product(s->getR2(),(*it)->getR2());
			}
		}
	}

	//std::cout <<"Here MIdx " <<minimumIndex << " "<<s->getID()<< " "<<signR2<<" SelfIntersect"<<this->self_intersect<<std::endl;
	if ((abs(minimumIndex - s->getID()) > 3) && (abs(minimumIndex - this->SeedPointID) > 5) && (minimumIndex != 10000)) {

//		std::cout <<"SignR2 " <<signR2 << " Min Index "<< minimumIndex<<" Flipped Seed Segment "<< this->SeedPointID<<std::endl;
		this->EndType = 32;	//setting the vessel endtype
		s->setEndType(32);	//setting the segment endtype
		this->SegmentContainer.push_back(s);
		return 0;
	}

	if ((minimumIndex != 10000)	 && (signR2 < 0) && (abs(minimumIndex - this->SeedPointID) > 5)){
		this->self_intersect++;
//		std::cout <<"SelfIntersect :"<<this->self_intersect<<" " <<signR2 << " "<< minimumIndex<<std::endl;
	}


	if (this->self_intersect >= 3)	{
		this->EndType = 33;	//setting the vessel endtype
		s->setEndType(33);	//setting the segment endtype
		this->SegmentContainer.push_back(s);
		return 0;
	}
*/
	last_predicted_steps[0] = last_predicted_steps[1];
	last_predicted_steps[1] = last_predicted_steps[2];
	last_predicted_steps[2] = last_predicted_steps[3];
	last_predicted_steps[3] = last_predicted_steps[4];
	last_predicted_steps[4] = s->getPredictedStep();

	last_traversed_steps[0] = last_traversed_steps[1];
	last_traversed_steps[1] = last_traversed_steps[2];
	last_traversed_steps[2] = last_traversed_steps[3];
	last_traversed_steps[3] = last_traversed_steps[4];
	last_traversed_steps[4] = s->getLength();

	double local_length_traversed, local_length_predicted;
	local_length_traversed = last_traversed_steps[0] + last_traversed_steps[1] + last_traversed_steps[2] + last_traversed_steps[3] + last_traversed_steps[4];
	local_length_predicted = last_predicted_steps[0] + last_predicted_steps[1] + last_predicted_steps[2] + last_predicted_steps[3] + last_predicted_steps[4];
//	std::cout<<	"local_length : (predicted)" << local_length_predicted <<" (traversed)" << local_length_traversed<< std::endl;
	if (local_length_traversed < 0.2*local_length_predicted)	{
			this->EndType = 5;	//setting the vessel endtype
			s->setEndType(5);	//setting the segment endtype
			this->SegmentContainer.push_back(s);
			return 0;
	}

	else	{
		//Add the segment to vessel
		this->SegmentContainer.push_back(s);
		return 1;
	}

}

bool Vessel::IsVesselValid()	{
	this->VesselLength = 0;
	Vect2 mu1, mu2, dmu;

	std::vector<Segment2D*>::iterator it;
	it = SegmentContainer.begin();

	while(it != SegmentContainer.end())	{
		mu1 = (*it)->getmu();
		++it;

		if(it == SegmentContainer.end())
			break;

		mu2 = (*it)->getmu();
		dmu = mu2 - mu1;
		this->VesselLength += dmu.magnitude();
	//	std::cout << "Lenght: " << dmu.magnitude() << std::endl;
	}




	//Implement MinVesselLength here. Reject a vessel if the length is less than MinVesselLength and
	// the vessel is not a HIT
	if (((this->VesselLength < this->MinVesselLength)||(SegmentContainer.size() < 5))  && (this->EndType != 22) )	{
		//std::cout << "Vessel length: " <<this->VesselLength << " EndType: "<< this->EndType << " ...REJECTED" <<std::endl;
		return 0;
	}
	else	{
		//std::cout << "Vessel length: " <<this->VesselLength << " EndType: "<< this->EndType << " ...ACCEPTED" <<std::endl;
		return 1;
	}

}

bool Vessel::IsPointInsideSE(Vect2 &p, Segment2D *s)	{
	Vect2 dmu, R1, R2,D;
	double a1, a2, d;
	R1[0] = cos(s->getq());
	R1[1] = sin(s->getq());
	R2[0] = -1*R1[1];
	R2[1] = R1[0];
	a1 = s->geta1();
	a2 = s->geta2() / 2;
	dmu = p - (s->getmu() - a2*R2); // this is BK in the code

	D[0] = dot_product(R1,dmu);
	D[1] = dot_product(R2,dmu);
	D[0] = D[0]/a1;
	D[1] = D[1]/a2;

	d = D[0]*D[0] + D[1]*D[1];
	if (d < 1)
		return 1;
	else
		return 0;
}

void Vessel::PrintSelf()	{
	std::cout<<"Vessel ID "<< this->getID() << "  Length "<<this->VesselLength << " ( "<< this->SegmentContainer.size() <<" segs)"<<std::endl;
	std::vector<Segment2D*>::iterator it;
	for (it = SegmentContainer.begin(); it != SegmentContainer.end(); ++it)
		std::cout<< (*it)->getID() <<" (" << (*it)->getEndType() <<") ";

	std::cout<<" Segment EndType "<< SegmentContainer.back()->getEndType()<<" Vessel Endtype "<< this->EndType <<std::endl;
}
