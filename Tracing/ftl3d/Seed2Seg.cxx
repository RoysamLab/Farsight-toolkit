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

# include "Seed2Seg.h"

Seed2Seg::Seed2Seg()	{
	iterations = 100.0;
	AS_RATIO = 2.0;
	min_a = 1.5;
	THRESH = 0.3;
	minL = 3.0;
}

void Seed2Seg::Configure(double FitIterations, double AspectRatio, 
						 double minVesselWidth, double startThresh, double minContrast)	
{
	this->iterations = FitIterations;
	this->AS_RATIO = AspectRatio;
	this->min_a = minVesselWidth;
	this->THRESH = startThresh;
	this->minL = minContrast;
}


Seed2Seg::~Seed2Seg() {
	std::vector<TVessel *>::iterator it;
	for (it = SSContainer.begin(); it != SSContainer.end(); ++it)
         delete *it;
}


void Seed2Seg::ComuputeStartSegments(SeedContainer3D::Pointer seeds, ImageType3D::Pointer im)	
{

	long int ID = 0;
	SeedContainerType SeedContainer = seeds->getContainer();
	SeedContainerType::iterator it;
	for (it = SeedContainer.begin(); it != SeedContainer.end(); ++it)	{
		if (!HitTestSeedPoints(*it))	{
			TVessel* segment = new  TVessel();
			//copy seed information here
			Vect3 mu = (*it)->getPosition() ;
			segment->mu[0] = mu[0];
			segment->mu[1] = mu[1];
			segment->mu[2] = mu[2];
			segment->ID = ID;

			SegInit *fitter = new SegInit();
			bool ret = fitter->fitSE(im, *segment, this->iterations, this->AS_RATIO, this->THRESH);

			if (ret) {
				if (!HitTestSeg(segment) && IsStartSegmentValid(segment))	{
					SSContainer.push_back(segment);
					ID = ID+1;
					std::cout<<"*";  //accepted segments
				}
				else {
					std::cout<<"|";	//good segments but hit 
				}
			}
			else {
				std::cout<<".";
			}
			delete fitter;
		}
	}
	std::cout << "Total number of Start Segments detected " << SSContainer.size() << std::endl;
}



bool Seed2Seg::IsStartSegmentValid(TVessel* seg)	{
	if ((seg->a1 < this->min_a) && (seg->a2 < this->min_a)&& (seg->a3 <= (2 * this->min_a)))	{
		return false;
	}
	if (((seg->L/seg->MAD) < this->THRESH) || (seg->L < this->minL ))	{
		return false;
	}
	return true;
}

bool Seed2Seg::HitTestSeedPoints(SeedPoint3D* sd) {

	Vect3 mu1 = sd->getPosition();
	StartSegContainerType::iterator it;

	double min_D = 1.0;
	double loc = -2.0;

	for (it = SSContainer.begin(); it != SSContainer.end(); ++it)	{

		//(*it)->PrintSelf();
		double* mu2 = (*it)->mu;
		double a1 = (*it)->a1;
		double a2 = (*it)->a2;
		double a3 = (*it)->a3;
		double *R1 = (*it)->R1;
		double *R2 = (*it)->R2;
		double *R3 = (*it)->R3;

		double u = R1[0]*(mu1[0]-mu2[0]) + R1[1]*(mu1[1]-mu2[1]) + R1[2]*(mu1[2]-mu2[2]);
		double v = R2[0]*(mu1[0]-mu2[0]) + R2[1]*(mu1[1]-mu2[1]) + R2[2]*(mu1[2]-mu2[2]);
		double w = R3[0]*(mu1[0]-mu2[0]) + R3[1]*(mu1[1]-mu2[1]) + R3[2]*(mu1[2]-mu2[2]);

		double D = vcl_pow(u/a1,2.0)+vcl_pow(v/a2,2.0)+vcl_pow(w/a3,2.0);


		if(D < min_D) {
		     min_D = D;
		     loc = (*it)->ID;
	    }

	}

	if (loc > -1)	{
		return 1;
	}
	else 	{
		return 0;
	}
}

bool Seed2Seg::HitTestSeg(TVessel*  s) {

	double p[6][3];
	double mu1[3];
	mu1[0] = s->mu[0];
	mu1[1] = s->mu[1];
	mu1[2] = s->mu[2];

	//along R1
	p[0][0] = s->mu[0] + (s->a1*s->R1[0]);
	p[0][1] = s->mu[1] + (s->a1*s->R1[1]);
	p[0][2] = s->mu[2] + (s->a1*s->R1[2]);

	p[1][0] = s->mu[0] - (s->a1*s->R1[0]);
	p[1][1] = s->mu[1] - (s->a1*s->R1[1]);
	p[1][2] = s->mu[2] - (s->a1*s->R1[2]);

	//along R2
	p[2][0] = s->mu[0] + (s->a2*s->R2[0]);
	p[2][1] = s->mu[1] + (s->a2*s->R2[1]);
	p[2][2] = s->mu[2] + (s->a2*s->R2[2]);

	p[3][0] = s->mu[0] - (s->a2*s->R2[0]);
	p[3][1] = s->mu[1] - (s->a2*s->R2[1]);
	p[3][2] = s->mu[2] - (s->a2*s->R2[2]);

	//along R3
	p[4][0] = s->mu[0] + (s->a3*s->R3[0]);
	p[4][1] = s->mu[1] + (s->a3*s->R3[1]);
	p[4][2] = s->mu[2] + (s->a3*s->R3[2]);

	p[5][0] = s->mu[0] - (s->a3*s->R3[0]);
	p[5][1] = s->mu[1] - (s->a3*s->R3[1]);
	p[5][2] = s->mu[2] - (s->a3*s->R3[2]);

	double min_D = 1.0;
	double loc = -2.0;

	StartSegContainerType::iterator it;
	for (it = SSContainer.begin(); it != SSContainer.end(); ++it)	{


		if (vnl_math_abs(s->mu[0] - (*it)->mu[0]) + vnl_math_abs(s->mu[1] - (*it)->mu[1]) + vnl_math_abs(s->mu[2] - (*it)->mu[2]) > 50)	{
			continue;
		}


		double mu2[3] = {0 ,0 ,0};
		double a1 = (*it)->a1;
		double a2 = (*it)->a2;
		double a3 = (*it)->a3;
		double *R1 = (*it)->R1;
		double *R2 = (*it)->R2;
		double *R3 = (*it)->R3;

		for (unsigned int i=0; i<=5; i++)	{
			mu2[0] = p[i][0];
			mu2[1] = p[i][1];
			mu2[2] = p[i][2];

			double u = R1[0]*(mu1[0]-mu2[0]) + R1[1]*(mu1[1]-mu2[1]) + R1[2]*(mu1[2]-mu2[2]);
			double v = R2[0]*(mu1[0]-mu2[0]) + R2[1]*(mu1[1]-mu2[1]) + R2[2]*(mu1[2]-mu2[2]);
			double w = R3[0]*(mu1[0]-mu2[0]) + R3[1]*(mu1[1]-mu2[1]) + R3[2]*(mu1[2]-mu2[2]);

			double D = vcl_pow(u/a1,2.0)+vcl_pow(v/a2,2.0)+vcl_pow(w/a3,2.0);

			if(D < min_D) {
				 min_D = D;
				 loc = (*it)->ID;
			}

		}
	}
	if (loc > -1)	{
		return 1;
	}
	else 	{
		return 0;
	}
}


bool Sortfunction(TVessel*s1, TVessel*s2)	{
	double w1 = s1->L * s1->a1 * s1->a2;
	double w2 = s2->L * s2->a1 * s2->a2;
	return( w1 > w2);
}

void Seed2Seg::SortStartSegments() {

//	P("Start segment container sorting")
//	SSContainer[1]->PrintSelf();
	std::sort(SSContainer.begin(), SSContainer.end(), Sortfunction);
	std::cout << " Start segments Likelihood values, from " << \
	SSContainer[0]->L << " to " << SSContainer[SSContainer.size()-1]->L << "." <<std::endl;
//	SSContainer[1]->PrintSelf();
//	S("Before and after sorting")

}
