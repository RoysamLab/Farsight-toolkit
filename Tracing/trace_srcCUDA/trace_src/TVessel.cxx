///////////////////////////////////////////////////////////////////////////
//Objects of size TVessel encode the actual superellipsoid model fits.
//Class contains all the relevant superellipsoid parameters.  Functionality
//is limited to the detection of valid model fits.
//To estimate parameters see SegInit.cxx and SegFit.cxx.

#include "TVessel.h"


TVessel::TVessel()	{
    for (int i=0; i<4; i++) {
    	q3[i] = 0.0;
    	q2[i] = 0.0;
    	q1[i] = 0.0;
	}

    q1[0] = 1.0;
    q2[0] = 1.0;
    q3[0] = 1.0;
    a1 = 4.0;
    a2 = 4.0;
    a3 = 4.0;
    f = 0;
    b = 255;
    for (int i=0; i<3; i++) {
	  	mu[i] = 0.0;
	  	last_dir[i] = 0.0;
		last_pos[i] = 0.0;
		R1[i] = 0.0;
		R2[i] = 0.0;
		R3[i] = 0.0;
	}
    e1 = 1.0;
    e2 = 1.0;
    L = 0;

    for (int i=0; i<24; i++) {
		bndy[i] = 0.0;
	}
    type = 0;
    index = 0;
    MAD = 0;
    ID = -1;
    TraceID = -1;
    K1 = 0;
    K2 = 0;
    a = 0;
    k = 0;
    R1[0] = 1.0;
    R2[1] = 1.0;
    R3[2] = 1.0;

    numNbrs = 0;
    for (int i=0; i<4; i++) {
		NbrID[i] = -1;
	}

}

TVessel::~TVessel()	{
}

bool TVessel::IsSegmentValid(TVessel* refSeg)	{

	if (this->L < this->MAD*.5 || this->L < 3.0 || this->f < 0 || this->b < 0 )	{
		return 0;
	}

	double maxW = vnl_math_max(refSeg->a1, refSeg->a2);
	double W = vnl_math_max(this->a1, this->a2);
	if (W > 2*maxW)	{
		return 0;
	}

	double Proj = this->R3[0]*refSeg->R3[0] + this->R3[1]*refSeg->R3[1] + this->R3[2]*refSeg->R3[2];
	if (Proj < -0.25)	{
		return 0;
	}

	return 1;
}



void TVessel::PrintSelf()	{
		std::cout << std::endl <<"ID = " << ID << " TraceID = " << TraceID << " @ [" << mu[0] << "," << mu[1] << "," <<  mu[2] << "]" << std::endl;
		std::cout << " A:<" << a1 << ", " << a2 << ", " <<  a3 << ">\n Q:<" << q1[0]<< ", "<< q1[1] << ", "<< q1[2] << ", " << q1[3] << ">" << std::endl;
		std::cout << " R:\t" << R1[0] << " " << R2[0] << " " <<  R3[0] << std::endl;
		std::cout << " \t" << R1[1] << " " << R2[1] << " " <<  R3[1] << std::endl;
		std::cout << " \t" << R1[2] << " " << R2[2] << " " <<  R3[2] << std::endl;
		std::cout << " Foregd:" << f << "  Backgd:" << b << "  Lhood:" <<  L << "  MAD:" <<  MAD << std::endl;
		std::cout << this->numNbrs << " neighbors :" ;
			for (unsigned int i=0; i<this->numNbrs; i++)	{
				std::cout << NbrID[i] <<", ";
		}
		std::cout << std::endl;
}


TVessel* TVessel::CopyToNewSegment()	{
		// also copy the start segment to the Node container
		TVessel *seg1 = new TVessel();

		for (int i=0; i<4; i++)	{
			seg1->q3[i] = this->q3[i];
			seg1->q2[i] = this->q2[i];
			seg1->q1[i] = this->q1[i];
			seg1->NbrID[i] = -1;
		}

		for (int i=0; i<3; i++)	{
			seg1->mu[i] = this->mu[i];
			seg1->last_dir[i] = this->last_dir[i];
			seg1->last_pos[i] = this->last_pos[i];
			seg1->R1[i] = this->R1[i];
			seg1->R2[i] = this->R2[i];
			seg1->R3[i] = this->R3[i];
		}

		for (int i=0; i<24; i++)	{
			seg1->bndy[i] = this->bndy[i];
		}

		seg1->a1 = this->a1;
		seg1->a2 = this->a2;
		seg1->a3 = this->a3;
		seg1->f = this->f;
		seg1->b = this->b;
		seg1->e1 = this->e1;
		seg1->e2 = this->e2;
		seg1->L = this->L;
		seg1->type = this->type;
		seg1->index = this->index;
		seg1->MAD = this->MAD;
		seg1->ID = this->ID;
		seg1->TraceID = this->TraceID;
		seg1->K1 = this->K1;
		seg1->K2 = this->K2;
		seg1->a = this->a;
		seg1->k = this->k;
	    seg1->numNbrs = 0;
	    return seg1;
}
