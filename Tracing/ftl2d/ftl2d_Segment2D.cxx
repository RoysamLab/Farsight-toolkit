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

#include "ftl2d_Segment2D.h"
#include <vnl/vnl_math.h>
#define DEBUG 0

Segment2D::Segment2D(TraceConfig *config) {
	//Loading default values
    this->SegmentID = 0;
    this->VesselID = 1000;	//unassigned
	this->mu[0] = 0.0;
	this->mu[1] = 0.0;
	this->a1 = 1.0;
	this->a2 = 1.0;
	this->e1 = 1.0;
	this->q = 0.0;
	this->ForegdIntensity = 0.0;
	this->BackgdIntensity = 255.0;
	this->MADEstimate = 0;
	this->Likelihood = -1000;
	this->EndType = 0;
	this->last_pos = this->mu;
	this->last_dir[0]=0;this->last_dir[1]=1;
	this->R1[0]=0;		this->R1[1]=0;
	this->R2[0]=1;		this->R2[0]=0;
	this->SelfIntersectCounter = 0;
	this->MinimumWidth = config->getMinimumVesselWidth();
	this->MaximumWidth = config->getMaximumVesselWidth();
	this->PROP = config->getPROP();
	this->PredictedStep = 1.0;
	this->AspectRatio = config->getAspectRatio();
	this->StepSize = config->getStepSize();
	//this->Comment = "Any Comment";

}

Segment2D::Segment2D(Segment2D *s_old)	{
	    this->SegmentID = s_old->SegmentID+1;
	    this->VesselID = s_old->VesselID;	//unassigned
		this->mu[0] = s_old->mu[0] ;
		this->mu[1] = s_old->mu[1];
		this->a1 = s_old->a1;
		this->a2 = s_old->a2;
		this->e1 = s_old->e1;
		this->q = s_old->q;
		this->ForegdIntensity = s_old->ForegdIntensity;
		this->BackgdIntensity = s_old->BackgdIntensity;
		this->MADEstimate = s_old->MADEstimate;
		this->Likelihood = s_old->Likelihood;
		this->EndType = s_old->EndType;
		this->last_pos = s_old->last_pos;
		this->last_dir[0] = s_old->last_dir[0];
		this->last_dir[1] = s_old->last_dir[1];
		this->R1[0] = s_old->R1[0];		
		this->R1[1] = s_old->R1[1];
		this->R2[0] = s_old->R2[0];
		this->R2[0] = s_old->R2[0];
		this->SelfIntersectCounter = s_old->SelfIntersectCounter;
		this->MinimumWidth = s_old->MinimumWidth;
		this->MaximumWidth = s_old->MaximumWidth;
		this->PROP = s_old->PROP;
		this->PredictedStep = s_old->PredictedStep;
		this->AspectRatio = s_old->AspectRatio;
		this->StepSize = s_old->StepSize;
}

Segment2D::~Segment2D() {
	//std::cout<<"Destructor of Segment"<<std::endl;
}

void Segment2D::advance( Segment2D * s_prev)	{
	//this function has to copy all attributes of s_prev to this.
	this->SegmentID = s_prev->getID()+1;
	this->a1 = s_prev->geta1();
	this->a2 = s_prev->geta2();
	this->e1 = s_prev->gete1();
	this->q = s_prev->getq();
	this->ForegdIntensity = s_prev->getForegdIntensity();
	this->BackgdIntensity = s_prev->getBackgdIntensity();
	this->MADEstimate =  s_prev->getMAD();
	this->Likelihood =  s_prev->getLikelihood();
	
	this->last_pos = s_prev->getmu();
	this->last_dir[0] = sin(s_prev->getq());
	this->last_dir[1] = cos(s_prev->getq());
	this->last_min_width = vnl_math_min(s_prev->geta1(),s_prev->geta2()); //Min a1 a2

	double min_width = vnl_math_min(s_prev->geta1(),s_prev->geta2());
	min_width = vnl_math_max(this->StepSize*min_width, 1.0);

	this->mu = this->last_pos + this->last_dir*min_width;
	this->PredictedStep = min_width;

}

void Segment2D::flipDirection(unsigned int ndx)	{
	this->q += vnl_math::pi;
	this->SegmentID = ndx;
}

void Segment2D::CopySeedInformation(SeedPoint2D* seed)	{
	this->mu[0] = static_cast<double> (seed->getx());
	this->mu[1] = static_cast<double> (seed->gety());
	this->a1 = seed->getScale();
	this->a2 = seed->getScale();
	this->ForegdIntensity = seed->getIntensity();
	this->BackgdIntensity = 70;  //fixed
}

bool Segment2D::IsSeedSegmentValid()	{

//	if ((this->Likelihood < this->PROP*this->MADEstimate) || (this->a2 < 1.35*this->a1) || (this->EndType != 0)) {
		//if ((this->Likelihood < 0.3*this->PROP*this->MADEstimate) || (this->EndType != 0)) {
		if ((this->Likelihood < 3*this->PROP) || (this->EndType != 0)) {	
//		std::cout <<"Ratio of a2 to a1 (<1.35)" <<this->a2/this->a1 << std::endl;
		return(0);
	}
	else	{
		return(1);
	}
}

void Segment2D::SegmentValidCheck()	{
	double w = vnl_math_min(this->a1 , this->a2);
	
	if (this->SegmentID >200)
		this->EndType = 10;
		return;
	
	
/*	if (w > 2*this->last_min_width)	{
		this->EndType = 7;
		return;
	}

	if ((this->a1 < this->MinimumWidth) || (this->a1 > this->MaximumWidth))	{
		this->EndType = 8;
		return;
	}
*/
	if (this->Likelihood < this->PROP)	{
//	if (this->Likelihood < this->PROP*this->MADEstimate)	{
		this->EndType = 2;
		return;
	}
	
	this->EndType = 0;
}

bool Segment2D::IsSegmentValid()	{
	if ((this->EndType == 0)||(this->EndType < 50))
		return 1;
	else
		return 0;

}

void Segment2D::Initialize(ImageType::Pointer im)	{

	const unsigned int ITERATION = 60;
	unsigned char IsInit= 1;
	this->e1 = 0.75;
	
	SuperEllipsoid2D * se = new SuperEllipsoid2D(this->AspectRatio);
	
	
	for (int i = 0; i < ITERATION; i++)	{
		
		#if DEBUG == 1
			std::cout<< i <<"B  ";
			this->PrintSelf(2);
			std::cin.get();
		#endif
		
		if ((i%10) == 0)	{
	//		std::cout << i;
			
			UpdateIntensity(im);
		}
				
		se->generateConvexHull(this->a1, this->a2, this-> e1, this->mu, this->q);
		
		
		if (!se->Interpolation(im))	{
			this->a2 = this->a2/1.5;
			se->generateConvexHull(this->a1, this->a2, this-> e1, this->mu, this->q);
			if (!se->Interpolation(im))	{
				this->EndType = 12;
				break;
			}
			this->EndType = 11;
		}
				
		se->UpdateLikelihood(this->ForegdIntensity, this->BackgdIntensity);
		
		se->gradient_superquad2d(this->a1, this->a2, this-> e1, this->mu, this->q);
		
		se->UpdateOrientation(this->q, 0); //0-no debug
		
		
		if (i >= 5)
			se->UpdateShapeParameter (this->e1, this->a1, this->a2, this->q);

		se->UpdateScale(this->q, this->a1, this->a2);
		
			
		se->UpdatePosition(this->mu, IsInit, this->q);
		
		se->UpdateAxis(this->a1, this->a2, this->q, this->last_dir);
		
		//Check for a1, a2 mu, q nan
		if ((std::isnan(a1)) || (std::isnan(a2)) || (std::isnan(mu[0])) || (std::isnan(mu[1])) || (std::isnan(q)))
			return;
	}

	if (this-> EndType != 12)	{
		UpdateIntensity(im);
		this->R1[0] = cos(this->q);
		this->R1[1] = -sin(this->q);
		this->R2[0] = sin(this->q);
		this->R2[1] = cos(this->q);
	}
	else	{
		this->Likelihood = -1000;
		this->MADEstimate = 0;
		this->q = 0;
	}
	
	delete se;

}

////////////////////////////////////////////////////////////////////////////////////////////
void Segment2D::fitSuperEllipsoid(ImageType::Pointer im)	{

	const unsigned int ITERATION = 20;
	unsigned char IsInit= 0;
	this->e1 = 0.75;
	this->a1 = 0.65*this->a1;
	SuperEllipsoid2D * se = new SuperEllipsoid2D(this->AspectRatio);


	for (int i = 0; i < ITERATION; i++)	{

		se->generateConvexHull(this->a1, this->a2, this-> e1, this->mu, this->q);
		
		if (!se->Interpolation(im))	{
			this->a2 = this->a2/1.5;
			se->generateConvexHull(this->a1, this->a2, this-> e1, this->mu, this->q);
			if (!se->Interpolation(im))	{
				this->EndType = 12;
				break;
			}
			this->EndType = 11;
		}

		se->UpdateLikelihood(this->ForegdIntensity, this->BackgdIntensity);
		se->gradient_superquad2d(this->a1, this->a2, this-> e1, this->mu, this->q);
		
		//std::cout<< i <<"B  " << std::isnan(q); this->PrintSelf(2);			std::cin.get();
		
		se->UpdateOrientation(this->q, 0); //1 - debug

		if (i >= 5)
			se->UpdateShapeParameter (this->e1, this->a1, this->a2, this->q);


		se->UpdateScale(this->q, this->a1, this->a2);
		se->UpdatePosition(this->mu, IsInit, this->q);
		se->UpdateAxis(this->a1, this->a2, this->q, this->last_dir);
		
		//Check for a1, a2 mu, q nan
		if ((std::isnan(a1)) || (std::isnan(a2)) || (std::isnan(mu[0])) || (std::isnan(mu[1])) || (std::isnan(q)))
			break;
	}

	if (this-> EndType != 12)	{
		UpdateIntensity(im);
		this->R1[0] = cos(this->q);
		this->R1[1] = -sin(this->q);
		this->R2[0] = sin(this->q);
		this->R2[1] = cos(this->q);

	}
	else	{
		this->Likelihood = -1000;
		this->MADEstimate = 0;
		this->q = 0;
	}
	
	delete se;
}




//This function updates the foreground and background intensities
bool Segment2D::UpdateIntensity(ImageType::Pointer im)	{

	double width = 1.5*vnl_math_max(this->a1,this->a2) ;	//4 * Min a1 a2

	ImageType::RegionType region;
	ImageType::RegionType::SizeType region_size;
	ImageType::RegionType::IndexType region_start;

	region_start[0] = static_cast<unsigned int>(this->mu(0) -  width);
	region_start[1] = static_cast<unsigned int>(this->mu(1) -  width);

	region_size[0] = static_cast<unsigned int>(2*width + 1);
	region_size[1] = region_size[0];

	region.SetSize( region_size );
  	region.SetIndex( region_start );
	
	
	if ( im->GetRequestedRegion().IsInside( region ) )    {

		IteratorType	It( im, region );
		double x, y,xx,yy, x1, y1, dist;
		double cs, sn;

		std::vector <float> FgList, BgList;
		FgList.reserve(2000);
		BgList.reserve(2000);
		cs = cos(this->q);
		sn = sin(this->q);
		
		for ( It.GoToBegin(); !It.IsAtEnd();    ++It)	{
			xx = static_cast<double>( It.GetIndex()[0] );
			yy = static_cast<double>( It.GetIndex()[1] );
			x = xx - this->mu(0);
			y = yy - this->mu(1);
			x1 = x*cs - y*sn;
			y1 = x*sn + y*cs;

			dist = pow((fabs(x1)/this->a1),(2.0/this->e1)) + pow((fabs(y1)/this->a2),(2.0/this->e1));
			
			//if (dist < 1)
			//	std::cout << dist << "\t" <<It.Get() << "\t[" <<xx <<","<<yy <<"]"<<std::endl;

			if (dist < 1)
				FgList.push_back(It.Get());
			else
				BgList.push_back(It.Get());
		
		}
		
		//std::cout << "Intensity TEST" << "\t" <<width << "\tst->" <<region_start <<" sz-> "<<region_size <<std::endl;
		//std::cout << "Intensity TEST" << "\t Fg->" <<FgList.size() << "\tBg->" <<BgList.size() <<std::endl;
		//EstimatorMAD(FgList, BgList);
		EstimatorMAD2(FgList, BgList);
		//std::cout << "Estimator Output {MAD, Likelihood}" << "\t" <<this->MADEstimate<< "\t" <<this->Likelihood <<" !!! "<<std::endl;
		
		FgList.clear();
		BgList.clear();
		//succesfully computed MADEstimate. F. B, Likelihood       
		return(1);
	}

	else
		return (0);
}

void  Segment2D::EstimatorMAD2(std::vector<float>& F, std::vector<float>& B)	{
	//In this implementation of the Intensity estimates, I have emlpoyed pure median as 
	//foregd and backgd estimates. 
	this->BackgdIntensity = ComputeMedian(B);
	this->ForegdIntensity = ComputeMedian(F);
	
	//The MAD has form
/*	std::vector<float> deviation;
	std::vector<float>::iterator it;
	float Bg_MAD = 0, Fg_MAD = 0, L = 0;
	
	//Background Median calculation
	for (it = B.begin(); it < B.end() ; ++it) {
		deviation.push_back(fabs((*it) - this->BackgdIntensity));
	}
	Bg_MAD = ComputeMedian(deviation) + 1.0;
	
	deviation.clear();
	//Foreground Median calculation
	for (it = F.begin(); it < F.end() ; ++it) {
		deviation.push_back(fabs((*it) - this->ForegdIntensity));
	}
	Fg_MAD = ComputeMedian(deviation) + 1.0;
	
	//Likelihood calculation
//	for (it = F.begin(); it < F.end() ; ++it)	{
//		L += (fabs((*it)-this->ForegdIntensity )/Fg_MAD) - (fabs((*it)-this->BackgdIntensity )/Bg_MAD);
//	}
		
//	this->Likelihood = -1*L/F.size();
*/	this->Likelihood = this->BackgdIntensity - this->ForegdIntensity;  
	
	
/*	if (this->a1 < 2.0) 	{
		this->MADEstimate = Bg_MAD; 
	}
	else {
		this->MADEstimate = Bg_MAD + Fg_MAD;   //Max of the two
	}

	if(this->MADEstimate < 1.0) 	{
*/		this->MADEstimate = 1.0;    //min MAD = 1
//	}
	
/*	std::cout << "Estimator Test: " << std::endl;
	std::cout 	<< " Fg: " << this->ForegdIntensity ;
	std::cout	<< " Bg: " << this->BackgdIntensity ;
	std::cout	<< " FgMAD: " << Fg_MAD ;
	std::cout	<< " BgMAD: " << Bg_MAD ; 
	std::cout	<< " L: " << L ;
	std::cout	<< " Likelihood: " << this->Likelihood ;
	std::cout	<<std::endl;
	std::cin.get();				
*/	
}


void  Segment2D::EstimatorMAD(std::vector<float>& F, std::vector<float>& B)	{

	std::vector<float> temp, deviation;
	std::vector<float>::iterator it;
	float T;
	float Bg_MAD = 0, Fg_MAD = 0, m = 0;
	float L = 0;

	T = 0.5*(this->ForegdIntensity + this->BackgdIntensity);
	
	
	//Background Median calculation
	for (it = B.begin(); it < B.end() ; ++it) {
		if ((*it > T) && (*it < 254))		//foreground < background
			temp.push_back(*it);
	}
	this->BackgdIntensity = ComputeMedian(temp);
	//std::sort(temp.begin(), temp.end());
	//this->BackgdIntensity = temp[temp.size()/2+1];
	
	//std::cout << "Estimator TEST 1" << "\t(F+B)/2=" <<T << "\ttemp size= " <<temp.size() <<" B= "<<this->BackgdIntensity <<std::endl;
	
	//recalculate the threshold
	T = 0.5*(this->ForegdIntensity + this->BackgdIntensity);
	
	//Background MAD calculation
	temp.clear();
	for (it = B.begin(); it < B.end() ; ++it)	{
		if ((*it > T) && (*it <=254))
			temp.push_back(*it);
	}

	m = ComputeMedian(temp);
	//std::sort(temp.begin(), temp.end());
	//m = temp[temp.size()/2+1];
	
	for(it = temp.begin(); it<temp.end(); ++it)
		deviation.push_back(fabs((*it) - m));
	
	Bg_MAD = ComputeMedian(deviation);
	//std::sort(deviation.begin(), deviation.end());
	//Bg_MAD = deviation[deviation.size()/2];
	
	//std::cout << "Estimator TEST 2" << "\t(F+B)/2=" <<T << "\ttemp size= " <<temp.size() <<" BgMAD= "<<Bg_MAD <<std::endl;


	//Foreground Median calculation
	temp.clear();
	deviation.clear();
	
	this->ForegdIntensity = ComputeMedian(F);	
//	std::sort(F.begin(), F.end());
//	this->ForegdIntensity = F[F.size()/2];
	
	//std::cout << "Estimator TEST 3: " ;
	//for (it = F.begin(); it < F.end() ; ++it)
	//	std::cout <<*it <<" "; 
	//std::cout<< std::endl;
	
	this->ForegdIntensity = vnl_math_min(this->ForegdIntensity, this->BackgdIntensity);
	
	T = 0.5*(this->ForegdIntensity + this->BackgdIntensity);
	
	for (it = F.begin(); it < F.end() ; ++it)	{
		if (*it < T)		//foreground < background
			temp.push_back(*it);
	}
	//Foreground MAD calculation
	m = ComputeMedian(temp);
	//std::sort(temp.begin(), temp.end());
	//m = temp[temp.size()/2+1];
	
	for(it = temp.begin(); it<temp.end(); ++it)
		deviation.push_back(fabs((*it) - m));

	Fg_MAD = ComputeMedian(temp);
	//std::sort(deviation.begin(), deviation.end());
	//Fg_MAD = deviation[deviation.size()/2];
	
	//std::cout << "Estimator TEST 4" << "\t(F+B)/2=" <<T << "\ttemp size= " <<temp.size() << "  Foregd= " << this->ForegdIntensity <<"  FgMAD= "<<Fg_MAD <<std::endl;

	if (this->a1 < 2.0)
		this->MADEstimate = Bg_MAD;
	else
		this->MADEstimate = vnl_math_max(Bg_MAD,Fg_MAD);	//Max of the two

	if(this->MADEstimate < 2.0)
		this->MADEstimate = 2.0;	//min MAD = 2
	
	//Likelihood calculation
	for (it = F.begin(); it < F.end() ; ++it)
		L += fabs((*it)-this->ForegdIntensity ) - fabs((*it)-this->BackgdIntensity );
	
	this->Likelihood = -1*L/F.size();
}

double Segment2D::getLength()	{
	Vect2 dmu = this->mu - this->last_pos;
	//std::cout<<"GetLengthFuntion ["<< this->mu <<"] ["<< this->last_pos <<"] - " << dmu.magnitude() <<std::endl;
	return (dmu.magnitude());
}

double Segment2D::ComputeMedian(std::vector<float>& v){
	std::sort(v.begin(), v.end());
	int r = v.size();
	if (r>2) {
		if (r%2)//odd
			return(v[(r-1)/2]);
		else
			return(0.5*(v[r/2]+v[r/2-1]));
		}
	else
		return 0;
}

void Segment2D::PrintSelf(unsigned char i)	{
	
	if (i==1)	{
		std::cout <<std::endl;
		std::cout << "Seg ID: "	<<SegmentID << std::endl;
		std::cout << "mu: " <<mu << std::endl;
		std::cout << "q: " 	<<q*180.0/vnl_math::pi <<" degrees" << std::endl;
		std::cout << "a1: " <<a1<< std::endl;
		std::cout << "a2: " <<a2 << std::endl;
		std::cout << "e1: " <<e1 << std::endl;
		std::cout << "R1: " <<R1<< std::endl;
		std::cout << "R2: " <<R2<< std::endl;
		std::cout << "ForegdIntensity: "<<ForegdIntensity	<< std::endl;
		std::cout << "BackgdIntensity: "<<BackgdIntensity	<< std::endl;
		std::cout << "MADEstimate: " <<MADEstimate 	<< std::endl;
		std::cout << "Likelihood: "	<<Likelihood 	<< std::endl;
		std::cout << "EndType: "	<<EndType 	<< std::endl;
	}
	else if (i == 2)	{	
	//	std::cout << "ID\tmu\tq\ta1\ta2\te1\tFg\tBg\tMAD\tL"<< std::endl;
		std::cout << SegmentID<<" mu:"<<mu<<" q:"<<q*180.0/vnl_math::pi<<" a1:"<<a1<<" a2:"<<a2<<" e1:"<<e1<<" F:"<<ForegdIntensity<<" B:"<<BackgdIntensity<<" M:"<<MADEstimate<<" L:"<<Likelihood<<" E:"<<EndType<< std::endl << std::endl;
	}
	
	else	{
		
		std::cout << SegmentID << " mu:"<<mu<< " q:" <<q*180.0/vnl_math::pi<<" step sz "<<this->getLength();
	}
}
