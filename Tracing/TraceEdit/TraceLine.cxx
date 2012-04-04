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

#include <vtksys/hash_map.hxx>
#include "tinyxml/tinyxml.h"

#define VTK_CREATE(type, var) \
  vtkSmartPointer<type> var = vtkSmartPointer<type>::New()

#define MY_ENCODING "ISO-8859-1"

#include "TraceLine.h"

///////////////////////////////////////////////////////////////////////////////
TraceLine::TraceLine()
{
	this->m_parent = NULL;
	this->root= -1;
	this->level = 0;
	this->m_id = -(1<<30);
	this->m_branches.clear();
	this->EuclideanD = -1;
	this->length = 0;
	this->radii = 0; 
	this->sectionArea = 0;
	this->surfaceArea = 0;
	this->volume = 0;
	this->BurkTaper = 0;
	this->HillmanTaper = 0;
	this->HillmanThreshold = 0;
	this->FileName = "default";

	this->daughterRatio = 0;
	this->parentDaughterRatio = 0;
	this->partitionAsymmetry = 0;

	this->rallPower = 0;
	this->Pk = 0;
	this->Pk_2 = 0;
	this->Pk_classic = 0;

	this->BifAmplLocal = -1;
	this->BifAmpRemote = -1;

	this->BifTiltLocal = -1;
	this->BifTiltRemote = -1;
	//there are 4 tilts per segment (need to add two more)
	this->BifTiltLocalBig = -1;
	this->BifTiltLocalSmall = -1;
	this->BifTiltRemoteBig = -1;
	this->BifTiltRemoteSmall = -1;

	this->BifTorqueLocal = -1;
	this->BifTorqueRemote = -1;
	this->BifTorqueLocalBig = -1;
	this->BifTorqueLocalSmall = -1;
	this->BifTorqueRemoteBig = -1;
	this->BifTorqueRemoteSmall = -1;

	this->BifTiltLocalTwoDaughter = -1;
	this->BifTiltRemoteTwoDaughter = -1;
	this->BifTorqueLocalTwoDaughter = -1;
	this->BifTorqueRemoteTwoDaughter = -1;

	this->CellFeatures.clear();
//for Dist to device calculations
	this->DistanceToROI = 0;
	this->ROICoord_X = 0;
	this->ROICoord_Y = 0;
	this->ROICoord_Z = 0;

	this->actualBifurcation = false;

	//// classification results
	//this->prediction = -1;
	//this->confidence = -1;
}

///////////////////////////////////////////////////////////////////////////////
TraceLine::~TraceLine()
{
}

///////////////////////////////////////////////////////////////////////////////
TraceLine* TraceLine::GetParent()
{
  return this->m_parent;
}
unsigned int TraceLine::GetParentID()
{
	if (this->m_parent)
	{
		return this->m_parent->GetId();
	}
	else 
	{
		return -1;
	}
}

///////////////////////////////////////////////////////////////////////////////
void TraceLine::SetParent(TraceLine* p)
{
  this->m_parent = p;
}
int TraceLine::GetRootID()
{
	return this->root;
}
int TraceLine::GetLevel()
{
	return this->level;
}
void TraceLine::setTerminalDegree(int degree)
{
	this->terminalDegree = degree;
}
bool TraceLine::isLeaf()
{
	if (this->m_branches.size() ==0)
	{
		return true;
	}
	else
	{
		return false;
	}
}
bool TraceLine::isRoot()
{
	if (!this->m_parent)
	{
		return true;
	}
	else
	{
		return false;
	}
}
bool TraceLine::isFree()
{
	if (this->isRoot() && this->isLeaf())
	{
		return true;
	}
	else
	{
		return false;
	}
}
bool TraceLine::isBranch()
{
	if (this->m_branches.size() > 1)
	{
		return true;
	}
	else
	{
		return false;
	}
}
void TraceLine::setActualBifurcation(bool bifurcate)
{
	this->actualBifurcation = bifurcate;
}
bool TraceLine::isActualBifurcation()
{
	return this->actualBifurcation;
}
void TraceLine::setRoot(int RootID, int traceLevel, double parentPath)
{
	this->root = RootID;
	this->level = traceLevel;
	this->PathLength = parentPath + this->GetLength();
}
void TraceLine::calculateVol()
{
	if (this->m_trace_bits.size() >1)
	{
		double dist = 0, r = 0, Df = 0, Dl = 0;
		TraceBit pre, cur;
		TraceBitsType::iterator it = this->m_trace_bits.begin();
		pre = *it; 
		it++;
		Df = 2*this->m_trace_bits.front().r;
		Dl = 2*this->m_trace_bits.back().r;
		r = pre.r;
		for (; it != this->m_trace_bits.end(); it++)
		{
			cur = *it;
			r += cur.r;
			dist += Euclidean(pre, cur);
			pre = cur;
		}
		if (Df !=Dl)
		{
			this->BurkTaper = (Df - Dl)/dist;
			this->HillmanTaper = Dl/Df;
		}
		else
		{
			this->BurkTaper = 0;
			this->HillmanTaper = 0;
		}
		this->length = dist;
		this->radii = r / this->m_trace_bits.size(); //ave radii of segment
		this->sectionArea = PI*pow((this->radii),2);
		this->volume = this->sectionArea*this->length;
		this->surfaceArea = 2*this->radii*PI*this->length;
	}//end size >1
	else
	{
		this->length = 0;
		this->radii = this->m_trace_bits.front().r; 
		this->sectionArea = PI*pow((this->radii),2);
		this->surfaceArea = 4*this->sectionArea;
		this->volume = (4/3)*PI*pow((this->radii),3);
		this->BurkTaper = 0;
		this->HillmanTaper = 0;
	}//end else
}//end vol calculation
void TraceLine::calculateBifFeatures()
{
	this->actualBifurcation = true;
	TraceBit BranchBit= this->m_trace_bits.back();
	TraceBit previousBit = this->GetBitXFromEnd(2);
	double BranchBitRadii = BranchBit.r;
	TraceLine* Daughter1 = this->GetBranch1();
	TraceBit  D1F = Daughter1->GetTraceBitsPointer()->front();
	TraceBit  D1B = Daughter1->GetTraceBitsPointer()->back();

	int n1 = Daughter1->GetTerminalDegree();
	double D1Radii = D1F.r;
	TraceLine* Daughter2 = this->GetBranch2();
	TraceBit  D2F = Daughter2->GetTraceBitsPointer()->front();
	TraceBit  D2B = Daughter2->GetTraceBitsPointer()->back();

	int n2 = Daughter2->GetTerminalDegree();
	double D2Radii = D2F.r;
	if (D1Radii > D2Radii)
	{
		this->daughterRatio = D1Radii /D2Radii;
	}
	else
	{
		this->daughterRatio = D2Radii / D1Radii;
	}
	
	Daughter1->SetParentDaughterRatio( D1Radii / BranchBitRadii);
	Daughter2->SetParentDaughterRatio( D2Radii / BranchBitRadii);
	this->HillmanThreshold = .5*BranchBitRadii + .25*D1Radii +.25*D2Radii;

	this->rallPower = this->RallPower(BranchBitRadii, D1Radii, D2Radii);
	//this->rallPower = pow ((1 + pow (ratio, power)), -1/power);	//need to solve ratio and power
	if (this->rallPower != -1)
	{
		this->Pk = this->CalculatePk(BranchBitRadii, D1Radii, D2Radii, this->rallPower);
	}
	this->Pk_classic = this->CalculatePk(BranchBitRadii, D1Radii, D2Radii, 1.5);
	this->Pk_2 = this->CalculatePk(BranchBitRadii, D1Radii, D2Radii, 2);

	if ((n1+n2) > 2)
	{
		this->partitionAsymmetry = abs(n1 - n2) / (n1 + n2 - 2);
	}else
	{
		this->partitionAsymmetry = 1;
	}

	this->BifAmplLocal = this->Angle(D1F, BranchBit, D2F);
	this->BifAmpRemote = this->Angle(D1B, BranchBit, D2B);

	//Identify bigger vs smaller daughter angle
	double angle1 = this->Angle(previousBit, BranchBit, D1F);
	double angle2 = this->Angle(previousBit, BranchBit, D2F);
	if (angle1 >= angle2)
	{
		this->setBifTiltLocalBig(angle1);
		this->setBifTiltLocalSmall(angle2);
	}
	else
	{
		this->setBifTiltLocalBig(angle2);
		this->setBifTiltLocalSmall(angle1);
	}

	this->setBifTiltLocal((angle1+angle2)/2); //average

	angle1 = this->Angle(previousBit, BranchBit, D1B);
	angle2 = this->Angle(previousBit, BranchBit, D2B);
	if (angle1 >= angle2)
	{
		this->setBifTiltRemoteBig(angle1);
		this->setBifTiltRemoteSmall(angle2);
	}
	else
	{
		this->setBifTiltRemoteBig(angle2);
		this->setBifTiltRemoteSmall(angle1);
	}

	this->setBifTiltRemote((angle1+angle2)/2);

	//Daughter1->setBifTiltLocal( this->Angle(previousBit, BranchBit, D1F));
	//Daughter2->setBifTiltLocal( this->Angle(previousBit, BranchBit, D2F));

	//Daughter1->setBifTiltRemote(this->Angle(previousBit, BranchBit, D1B));
	//Daughter2->setBifTiltRemote(this->Angle(previousBit, BranchBit, D2B));

	double* ParentPlaneLocal_ptr  = this->Plane(D1F, BranchBit, D2F);
	double* ParentPlaneRemote_ptr = this->Plane(D1B, BranchBit, D2B);
	double* Daughter1PlaneLocal_ptr;
	double* Daughter1PlaneRemote_ptr;
	double* Daughter2PlaneLocal_ptr;
	double* Daughter2PlaneRemote_ptr;
	double planeAnglelocal1 = -1;
	double planeAnglelocal2 = -1;
	double planeAngleremote1 = -1;
	double planeAngleremote2 = -1;

	if (!Daughter1->isLeaf())
	{
		TraceLine* GrandDaughter1of1 = Daughter1->GetBranch1();
		TraceLine* GrandDaughter2of1 = Daughter1->GetBranch2();
		TraceBit GD1F1 = GrandDaughter1of1->GetTraceBitsPointer()->front();
		TraceBit GD1B1 = GrandDaughter1of1->GetTraceBitsPointer()->back();
		TraceBit GD2F1 = GrandDaughter2of1->GetTraceBitsPointer()->front();
		TraceBit GD2B1 = GrandDaughter2of1->GetTraceBitsPointer()->back();
		Daughter1PlaneLocal_ptr  = this->Plane(GD1F1,D1B,GD2F1);
		Daughter1PlaneRemote_ptr = this->Plane(GD1B1,D1B,GD2B1);
		planeAnglelocal1 = this->PlaneAngle(ParentPlaneLocal_ptr, Daughter1PlaneLocal_ptr);
		planeAngleremote1 = this->PlaneAngle(ParentPlaneRemote_ptr,Daughter1PlaneRemote_ptr);
	}
	if (!Daughter2->isLeaf())
	{
		TraceLine* GrandDaughter1of2 = Daughter2->GetBranch1();
		TraceLine* GrandDaughter2of2 = Daughter2->GetBranch2();
		TraceBit GD1F2 = GrandDaughter1of2->GetTraceBitsPointer()->front();
		TraceBit GD1B2 = GrandDaughter1of2->GetTraceBitsPointer()->back();
		TraceBit GD2F2 = GrandDaughter2of2->GetTraceBitsPointer()->front();
		TraceBit GD2B2 = GrandDaughter2of2->GetTraceBitsPointer()->back();
		Daughter2PlaneLocal_ptr  = this->Plane(GD1F2,D2B,GD2F2);
		Daughter2PlaneRemote_ptr = this->Plane(GD1B2,D2B,GD2B2);
		planeAnglelocal2 = this->PlaneAngle(ParentPlaneLocal_ptr, Daughter2PlaneLocal_ptr);
		planeAngleremote2 = this->PlaneAngle(ParentPlaneRemote_ptr,Daughter2PlaneRemote_ptr);
	}

	if (planeAnglelocal1 >= planeAnglelocal2)
	{
		this->setBifTorqueLocalBig(planeAnglelocal1);
		this->setBifTorqueLocalSmall(planeAnglelocal2);
	}
	else
	{
		this->setBifTorqueLocalBig(planeAnglelocal2);
		this->setBifTorqueLocalSmall(planeAnglelocal1);
	}
	if (!Daughter1->isLeaf() || !Daughter2->isLeaf())
	{
		if (!Daughter1->isLeaf() && !Daughter2->isLeaf())
		{
			this->setBifTorqueLocal((planeAnglelocal1+planeAnglelocal2)/2);
			this->setBifTorqueLocalTwoDaughter(planeAnglelocal1+planeAnglelocal2);
		}
		else if (!Daughter1->isLeaf())
		{
			this->setBifTorqueLocal(planeAnglelocal1);
			this->setBifTorqueLocalTwoDaughter(planeAnglelocal1);
		}
		else
		{
			this->setBifTorqueLocal(planeAnglelocal2);
			this->setBifTorqueLocalTwoDaughter(planeAnglelocal2);
		}
	}

	if (planeAngleremote1 >= planeAngleremote2)
	{
		this->setBifTorqueRemoteBig(planeAngleremote1);
		this->setBifTorqueRemoteSmall(planeAngleremote2);
	}
	else
	{
		this->setBifTorqueRemoteBig(planeAngleremote2);
		this->setBifTorqueRemoteSmall(planeAngleremote1);
	}
	
	if (!Daughter1->isLeaf() || !Daughter2->isLeaf())
	{
		if (!Daughter1->isLeaf() && !Daughter2->isLeaf())
		{
			this->setBifTorqueRemote((planeAngleremote1+planeAngleremote2)/2);
			this->setBifTorqueRemoteTwoDaughter(planeAngleremote1+planeAngleremote2);
		}
		else if (!Daughter1->isLeaf())
		{
			this->setBifTorqueRemote(planeAngleremote1);
			this->setBifTorqueRemoteTwoDaughter(planeAngleremote1);
		}
		else
		{
			this->setBifTorqueRemote(planeAngleremote2);
			this->setBifTorqueRemoteTwoDaughter(planeAngleremote2);
		}
	}
}
void TraceLine::setTraceBitIntensities(vtkSmartPointer<vtkImageData> imageData)
{
	TraceBit curBit;
	TraceBitsType::iterator it = this->m_trace_bits.begin();
	for (; it != this->m_trace_bits.end(); it++)
	{
		int lx = 0, ly = 0, lz = 0; 
		curBit = *it;
		lx = (int) floor(curBit.x + 0.5);
		ly = (int) floor(curBit.y + 0.5);
		lz = (int) floor(curBit.z + 0.5);
		(*it).I = imageData->GetScalarComponentAsDouble(lx,ly,lz,0);
		std::cout<< "\nid:\t"<< (*it).id << "\tI:\t" << (*it).I;
	}
	//should have ended
}
double TraceLine::GetEuclideanLength()
{
	if (this->m_trace_bits.size() <2)
	{
		this->EuclideanD =0;
	}
	else
	{
		TraceBit front = this->m_trace_bits.front();
		TraceBit back  = this->m_trace_bits.back();
		this->EuclideanD = this->Euclidean(front, back);
	}
	return this->EuclideanD;
}
double TraceLine::GetBitDensity()
{
	if (this->GetSize() >1)
	{
		this->BitDensity = this->GetSize() / this->GetLength();
	}
	else
	{
		this->BitDensity = 1;
	}
	return this->BitDensity;
}
double TraceLine::GetDistToParent()
{
	if (this->m_parent)
	{
		this->DistToParent = this->Euclidean(this->m_trace_bits.front(), //
			this->m_parent->m_trace_bits.back());
		if (this->m_trace_bits.size()>1)
		{
			double Leading = 0;//, dist =0;
			TraceBit pre, cur;
			TraceBitsType::iterator it = this->m_trace_bits.begin();
			pre = *it; 
			it++;
			cur = *it;
			Leading = Euclidean(pre, cur);
			if (Leading > 2*this->DistToParent)
			{
				this->DistToParent = Leading;
			}
		}
		return this->DistToParent;
	}
	else
	{
		this->DistToParent = 0;
		return -1;
	}
}
double TraceLine::GetFragmentationSmoothness()
{
	if (!(this->EuclideanD > -1))
	{
		this->GetEuclideanLength();
	}
	double t = -1;
	if ( this->m_trace_bits.size() > 1)
	{
		t = this->length/this->EuclideanD;
	}
	return t;
}
///////////////////////////////////////////////////////////////////////////////
void TraceLine::AddBranch(TraceLine* b)
{ 
  this->m_branches.push_back(b);
}

///////////////////////////////////////////////////////////////////////////////
TraceLine* TraceLine::GetBranch1()
{
  return this->m_branches[0];
}

///////////////////////////////////////////////////////////////////////////////
void TraceLine::SetBranch1(TraceLine* b0)
{
  this->m_branches[0] = b0;
}

///////////////////////////////////////////////////////////////////////////////
TraceLine* TraceLine::GetBranch2()
{
  return this->m_branches[1];
}

///////////////////////////////////////////////////////////////////////////////
void TraceLine::SetBranch2(TraceLine* b1)
{
  this->m_branches[1] = b1;
}

///////////////////////////////////////////////////////////////////////////////
unsigned char TraceLine::GetType()
{
  return this->m_type;
}

///////////////////////////////////////////////////////////////////////////////
void TraceLine::SetType(unsigned char t)
{
  this->m_type = t;
}

///////////////////////////////////////////////////////////////////////////////
void TraceLine::AddTraceBit(TraceBit tbit)
{
  this->m_trace_bits.push_back(tbit);
}
void TraceLine::ExtendTrace(TraceBit tbit)
{
	bool front = false;
	if (this->m_trace_bits.size()>1)
	{
		if (!this->Orient(tbit))
		{
			this->m_trace_bits.push_front(tbit);
			front = true;
		}
	}
	if (!front)
	{
		this->m_trace_bits.push_back(tbit);
	}
}
TraceBit TraceLine::removeLastBit()
{
	TraceBit lastBit = this->m_trace_bits.back();
	if (this->m_trace_bits.size() > 1)
	{
		this->m_trace_bits.pop_back();
	}
	return lastBit;
}
TraceBit TraceLine::removeFirstBit()
{
	TraceBit firstBit = this->m_trace_bits.front();
	if (this->m_trace_bits.size() > 1)
	{
		this->m_trace_bits.pop_front();
	}
	return firstBit;
}

//This is an error fix
bool TraceLine::removeLeadingBit()
{
	if (this->m_trace_bits.size() < 3)
	{
		return false;
	}else
	{		
		//std::vector<double> dist;// = 0, r = 0;
		double Leading = 0, dist =0;
		TraceBit pre, cur;
		TraceBitsType::iterator it = this->m_trace_bits.begin();
		pre = *it; 
		it++;
		cur = *it;
		Leading = Euclidean(pre, cur);
		it++;
		for (; it != this->m_trace_bits.end(); it++)
		{
			cur = *it;
			dist += this->Euclidean(pre, cur);
			pre = cur;
		}//end for
		if(Leading > (dist/ (this->m_trace_bits.size() -1)))
		{
			this->m_trace_bits.pop_front();
			return true;
		}//end dist
		return false;
	}//end else size
}
///////////////////////////////////////////////////////////////////////////////
TraceLine::TraceBitsType::iterator TraceLine::GetTraceBitIteratorBegin()
{
  return this->m_trace_bits.begin();
}
TraceBit TraceLine::GetBitXFromBegin(int x)
{
	TraceBit curBit;
	TraceBitsType::iterator iter = this->m_trace_bits.begin();
	TraceBitsType::iterator iterend = this->m_trace_bits.end();
	if (3 > (int)this->m_trace_bits.size())
	{
		--iterend;
		curBit = * iterend;
		//return  curBit;
	}
	else
	{
		int i = 0;
		while(( i< (int)(this->m_trace_bits.size()-2))&&( i < x ))
		{
			++iter;
			i++;
		}
		curBit = *iter;
	}
	return curBit;
}
///////////////////////////////////////////////////////////////////////////////
TraceLine::TraceBitsType::iterator TraceLine::GetTraceBitIteratorEnd()
{
  return this->m_trace_bits.end();
}

TraceBit TraceLine::GetBitXFromEnd(int x)
{
	TraceBit curBit;
	TraceBitsType::iterator iter = this->m_trace_bits.begin();
	TraceBitsType::iterator iterend = this->m_trace_bits.end();
	if (3 >(int)this->m_trace_bits.size())
	{
		curBit = * iter;
		
	}
	else
	{
		int i = 0;
		while(( i<(int)(this->m_trace_bits.size()-2))&&( i < x ))
		{
			--iterend;
			i++;
		}
		curBit = *iterend;
	}
	return curBit;
}
///////////////////////////////////////////////////////////////////////////////
TraceLine::TraceBitsType* TraceLine::GetTraceBitsPointer()
{
  return &m_trace_bits;
}

///////////////////////////////////////////////////////////////////////////////
void TraceLine::SetId(unsigned int lid){ m_id = lid;}

///////////////////////////////////////////////////////////////////////////////
unsigned int TraceLine::GetId()
{ 
  return m_id;
};

///////////////////////////////////////////////////////////////////////////////
int TraceLine::GetSize()
{
  return (int)m_trace_bits.size();
};

///////////////////////////////////////////////////////////////////////////////
void TraceLine::Print(std::ostream &c,int indent)
{
  for(int counter=0; counter<indent; counter++)
	  c<<" ";
  c<<"TraceLine: "<<std::endl;
  for(int counter=0; counter<indent; counter++)
	  c<<" ";
  c<<"Size: "<<m_trace_bits.size()<<std::endl;
  for(unsigned int counter=0; counter< m_branches.size(); counter++)
  {
	  m_branches[counter]->Print(std::cout,indent+4);
  }
}

///////////////////////////////////////////////////////////////////////////////
std::vector<unsigned int> * TraceLine::GetMarkers()
{
  return &(this->m_markers);
}

///////////////////////////////////////////////////////////////////////////////
std::vector<TraceLine*> * TraceLine::GetBranchPointer()
{
  return &(this->m_branches);
}

///////////////////////////////////////////////////////////////////////////////
void TraceLine::setTraceColor(double newColor)
{
  this->traceColor = newColor;
};

///////////////////////////////////////////////////////////////////////////////
double TraceLine::getTraceColor()
{
  return this->traceColor;
};

///////////////////////////////////////////////////////////////////////////////
 TraceLine::TraceLine(const TraceLine &t)
{
  this->traceColor = t.traceColor;
  this->m_id = t.m_id;
  this->m_markers = t.m_markers;
  this->m_type = t.m_type;
  this->m_parent = NULL;
  for(unsigned int counter=0; counter< t.m_branches.size(); counter++)
  {
    TraceLine *temp = new TraceLine();
    temp->SetParent(this);
    this->AddBranch(temp);
    *temp = *t.m_branches[counter];
  }
  this->m_trace_bits = t.m_trace_bits;
  this->FileName = t.FileName;
}

///////////////////////////////////////////////////////////////////////////////
void TraceLine::Getstats()
{
  double XF, XB, YF, YB, ZF, ZB;
  XF = m_trace_bits.front().x;
  XB = m_trace_bits.back().x;
  YF= m_trace_bits.front().y;
  YB= m_trace_bits.back().y;
  ZF= m_trace_bits.front().z;
  ZB= m_trace_bits.back().z;
  printf("Trace # %d \t Trace Size: \t %d ", m_id, (int)m_trace_bits.size());
  printf("First bit x: %4.2f y: %4.2f z: %4.2f \t", XF, YF, ZF); 
  printf("Endt bit x: %4.2f y: %4.2f z: %4.2f \n", XB, YB, ZB); 
}


///////////////////////////////////////////////////////////////////////////////
double TraceLine::Euclidean(TraceBit bit1, TraceBit bit2)
{
	double distance, x, y, z;
	x = pow((bit1.x -bit2.x),2);
	y = pow((bit1.y -bit2.y),2);
	z = pow((bit1.z -bit2.z),2);
	distance = sqrt(x +y +z);
	return distance;

}
double TraceLine::Angle(TraceBit bit1f, TraceBit bit1b, TraceBit bit2f, TraceBit bit2b)
{
	double delX1, delX2, delY1, delY2, delZ1, delZ2,  norm1, norm2, angle;
	//delta x,y,z 
	angle = 0;
	delX1= fabs (bit1f.x-bit1b.x);	delX2= fabs (bit2f.x-bit2b.x);
	delY1= fabs (bit1f.y-bit1b.y);	delY2= fabs (bit2f.y-bit2b.y);
	delZ1= fabs (bit1f.z-bit1b.z);	delZ2= fabs (bit2f.z-bit2b.z);
	norm1=sqrt(pow((delX1),2)+ pow((delY1),2)+ pow((delZ1),2));
	norm2=sqrt(pow((delX2),2)+ pow((delY2),2)+ pow((delZ2),2));
	angle = acos(((delX1 * delX2) + (delY1 *delY2) + (delZ1 * delZ2)) /
          (norm2 * norm1));
	if (!(angle >= 0))
	{//prevent nan rounding errors
		angle = 0;
	}
	return angle;
}
double TraceLine::Angle(TraceBit bit1, TraceBit vertex, TraceBit bit2)
{
	double delX1, delX2, delY1, delY2, delZ1, delZ2,  norm1, norm2;
	double NewAngle = 0;
	//delta x,y,z 	
	delX1= fabs (vertex.x-bit1.x);	delX2= fabs (vertex.x-bit2.x);
	delY1= fabs (vertex.y-bit1.y);	delY2= fabs (vertex.y-bit2.y);
	delZ1= fabs (vertex.z-bit1.z);	delZ2= fabs (vertex.z-bit2.z);
	norm1=sqrt(pow((delX1),2)+ pow((delY1),2)+ pow((delZ1),2));
	norm2=sqrt(pow((delX2),2)+ pow((delY2),2)+ pow((delZ2),2));
	NewAngle = acos(((delX1 * delX2) + (delY1 *delY2) + (delZ1 * delZ2)) /
          (norm2 * norm1));
	if (!(NewAngle >= 0))
	{//prevent nan rounding errors
		NewAngle = 0;
	}
	return (NewAngle *180 )/PI;
}
double TraceLine::AzimuthAngle(TraceBit vertex, TraceBit bit1)
{
	double delX, delY;
	double NewAngle = 0;
	//delta x,y
	delX = bit1.x - vertex.x;
	delY = bit1.y - vertex.y;

	NewAngle = atan2(delY,delX);
	//gives positive angle if counterclockwise and negative angle if clockwise

	return (NewAngle *180 )/PI;
}
double TraceLine::ElevationAngle(TraceBit vertex, TraceBit bit1)
{
	double delX, delY, delZ;
	double NewAngle = 0;
	//delta x,y,z
	delX = bit1.x - vertex.x;
	delY = bit1.y - vertex.y;
	delZ = bit1.z - vertex.z;
	
	double hypotenuse = sqrt(pow(delX,2) + pow(delY,2));
	NewAngle = atan2(delZ,hypotenuse);
	//gives positive angle if counterclockwise and negative angle if clockwise

	return (NewAngle *180 )/PI;
}
double TraceLine::RallPower(double diamParent, double diamD1, double diamD2)
{
	double m = .001;
	double min = 1000; //min should be less than max
	while (m < 5)
	{
		double a1 = pow(diamParent, m);
		double a2 = pow(diamD1, m) + pow(diamD2, m);
		double a3 = fabs(a1 - a2);
		if (a3 <= .001)
		{
			return m;
		}
		m += .001;
	}
	return -1;
}
bool TraceLine::EndPtDist(TraceLine *Trace2, int &dir1, int &dir2, double &dist,
                          double &maxdist, double &angle) 
{
  //int center1=this->GetSize()/2, center2=Trace2->GetSize()/2;
  double min, distances[4];

  int Xbits = 10; //small step from end to determine 
  //compute the endpt distances
  distances[0]=Euclidean(m_trace_bits.front(), Trace2->m_trace_bits.front());//0 F-F
  distances[1]=Euclidean(m_trace_bits.front(), Trace2->m_trace_bits.back());//1 F-B
  distances[2]=Euclidean(m_trace_bits.back(), Trace2->m_trace_bits.front());//2 B-F
  distances[3]=Euclidean(m_trace_bits.back(), Trace2->m_trace_bits.back());//3 B-B

  //determine minimum spacing
  min = distances[0];
  int i, mark=0;
  for (i = 1; i<4; i++)
    {
    if (min > distances[i])
      {
      min= distances[i];
      mark=i;
      }
    } 
// std::cout << mark << "=mark\n";
  // from min determine orientation and return distance
  if (mark ==0)
    {//F-F
		//std::cout <<"FF\n";
	if (this->m_branches.size() > 0 && Trace2->m_branches.size()> 0)
		{return false;}
	if(this->GetParentID()!=-1 || Trace2->GetParentID()!=-1)
	{return false;}
    dist = distances[0];
    dir1= m_trace_bits.front().marker;
    dir2= Trace2->m_trace_bits.front().marker; 
	angle = Angle(m_trace_bits.front(),this->GetBitXFromBegin(Xbits), 
		Trace2->m_trace_bits.front(),Trace2->GetBitXFromBegin(Xbits));
    //angle=PI-angle;
    maxdist= distances[3];
    }
  else if (mark ==1)
    {//1 F-B
		//std::cout <<"FB\n";
	if ((this->GetParentID() != -1)||!Trace2->isLeaf())
		{return false;}
    dist = distances[1];
    dir1= m_trace_bits.front().marker;
    dir2= Trace2->m_trace_bits.back().marker; 
	angle = Angle(m_trace_bits.front(),this->GetBitXFromBegin(Xbits), 
		Trace2->GetBitXFromEnd(Xbits),Trace2->m_trace_bits.back());
    maxdist= distances[2];
    }
  else if (mark ==2)
    {//2 B-F
		//std::cout <<"BF\n";
		if (!this->isLeaf() || (Trace2->GetParentID()!=-1))
		{return false;}
    dist = distances[2];
    dir1= m_trace_bits.back().marker;
    dir2= Trace2->m_trace_bits.front().marker; 
	angle = Angle(this->GetBitXFromEnd(Xbits),m_trace_bits.back(), 
		Trace2->m_trace_bits.front(),Trace2->GetBitXFromBegin(Xbits));
    maxdist= distances[1];
    }
  else
    {//3 B-B
		//std::cout <<"BB\n";
	if (!this->isLeaf() || !Trace2->isLeaf()) 
		{return false;}
    dist = distances[3];
    dir1= m_trace_bits.back().marker;
    dir2= Trace2->m_trace_bits.back().marker;
	angle = Angle(this->GetBitXFromEnd(Xbits),m_trace_bits.back(), 
		Trace2->GetBitXFromEnd(Xbits),Trace2->m_trace_bits.back());
    //angle=PI-angle;
    maxdist= distances[0];
    }
  return true;
}
double TraceLine::CalculatePk(double Dp, double Da, double Db, double n)
{
	return (pow(Da, n) + pow(Db, n))/pow(Dp, n);
}
bool TraceLine::Orient(TraceLine * Trunk)
{
	double distances[2];
	//compute the endpt distances
	distances[0]= Euclidean(m_trace_bits.front(),	Trunk->m_trace_bits.back());// F-B
	distances[1]= Euclidean(m_trace_bits.back(),	Trunk->m_trace_bits.back());// B-B
	if(distances[0] < distances[1])
	{
		return true;	//oriented correctly
	}
	return false;		//needs to be flipped
}
bool TraceLine::Orient(TraceBit bit)
{
	double distances[2];
	//compute the endpt distances
	distances[0]= Euclidean(m_trace_bits.front(),	bit);// F-Bit
	distances[1]= Euclidean(m_trace_bits.back(),	bit);// B-Bit
	if(distances[0] > distances[1])
	{
		return true;	//oriented correctly
	}
	return false;		//needs to be flipped
}
///////////////////////////////////////////////////////////////////////////////
void TraceLine::SetFileName(char *newFileName)
{
	this->FileName = newFileName;
}
const char * TraceLine::GetFileName()
{
	return this->FileName.c_str();
}
void TraceLine::getEndPtBounds(double bounds[])
{
	//x
	if ( this->m_trace_bits.front().x < this->m_trace_bits.back().x)
	{
		bounds[0] = this->m_trace_bits.front().x;
		bounds[1] = this->m_trace_bits.back().x;
	}
	else
	{
		bounds[1] = this->m_trace_bits.front().x;
		bounds[0] = this->m_trace_bits.back().x;
	}
	//y
	if ( this->m_trace_bits.front().y < this->m_trace_bits.back().y)
	{
		bounds[2] = this->m_trace_bits.front().y;
		bounds[3] = this->m_trace_bits.back().y;
	}
	else
	{
		bounds[3] = this->m_trace_bits.front().y;
		bounds[2] = this->m_trace_bits.back().y;
	}
	//z
	if ( this->m_trace_bits.front().z < this->m_trace_bits.back().z)
	{
		bounds[4] = this->m_trace_bits.front().z;
		bounds[5] = this->m_trace_bits.back().z;
	}
	else
	{
		bounds[5] = this->m_trace_bits.front().z;
		bounds[4] = this->m_trace_bits.back().z;
	}
	//debug
	for (int i = 0; i < 6; i++)
	{
		std::cout<< i <<" Value " << bounds[i];
	}
}
std::string TraceLine::stats()
{
	std::stringstream thisStats;
	thisStats << this->GetId();
	thisStats << "\t";
	thisStats << (int)this->GetType();
	thisStats << "\t" ;
	thisStats << this->GetSize();
	thisStats << "\t";
	thisStats << this->GetLength();
	thisStats << "\t" ;
	thisStats << this->GetEuclideanLength();
	thisStats << "\t" ;
	thisStats << this->GetRadii();
	thisStats << "\t";
	thisStats << this->GetFragmentationSmoothness() ;
	thisStats << "\t";
	thisStats << this->GetParentID();
	return thisStats.str();
}
std::string TraceLine::RootCoord()
{
	std::stringstream thisStats;
	thisStats << this->m_trace_bits.front().x;
	thisStats << "_";
	thisStats << this->m_trace_bits.front().y;
	thisStats << "_";
	thisStats << this->m_trace_bits.front().z;
	return thisStats.str();

}
std::string TraceLine::statHeaders()
{
	std::stringstream thisStatsHeaders; 
	thisStatsHeaders <<"ID"
		<<"\tType"
		<<"\tSize"
		<<"\tLength"
		<<"\tEuclidean L"
		<<"\tRadii"
		<<"\tContraction"
		<<"\tParent ID";
	return thisStatsHeaders.str();
}

double* TraceLine::Plane(TraceBit bit1, TraceBit vertex, TraceBit bit2)
{
	double vector[3];
	vector[0] = bit1.x - vertex.x;
	vector[1] = bit1.y - vertex.y;
	vector[2] = bit1.z - vertex.z;
	return vector;
}

double TraceLine::PlaneAngle(double *plane1, double *plane2)
{
	double top = plane1[0]*plane2[0] + plane1[1]*plane2[1] + plane1[2]*plane2[2];
	double bottom = sqrt(plane1[0]*plane1[0] + plane1[1]*plane1[1]+ plane1[2]*plane1[2])*sqrt(plane2[0]*plane2[0] + plane2[1]*plane2[1]+ plane2[2]*plane2[2]);
	return std::acos(top/bottom)*180/PI;
}

//void TraceLine::calculateAzimuthElevation()
//{
//	TraceBit P1 = this->GetParent()->GetTraceBitsPointer()->front();
//	TraceBit D1 = this->GetTraceBitsPointer()->front();
//	/*TraceBit Ref;
//	
//	Ref.x = D1.x;
//	Ref.y = D1.y;
//	Ref.z = P1.z;*/
//	
//	azimuthangle = this->Azimuth(P1,D1);
//	elevationangle = this->Elevation(P1,D1);
//}

double TraceLine::GetCompartmentCurvature()
{
	double angle = 0;
	if (this->m_parent)
	{
		if (this->m_trace_bits.size()>2)
		{
			TraceBit pre, cur, next;
			TraceBitsType::iterator it = this->m_trace_bits.begin();
			pre = *it;
			it++;
			cur = *it;
			it++;
			next = *it;
			angle = this->Angle(pre,cur,next);
		}
		return angle;
	}
	else
		return -1;
}

double TraceLine::GetAzimuth()
{
	double newAzimuthAngle = 0;
	if (this->m_parent)
	{
		if (this->m_trace_bits.size()>1)
		{
			//TraceBit pre, cur;
			//TraceBitsType::iterator it = this->m_trace_bits.begin();
			//pre = *it;
			//it++;
			//cur = *it;
			//newAzimuthAngle = this->AzimuthAngle(pre,cur);

			TraceBit leadBit = this->m_trace_bits.front();
			TraceBit endBit = this->m_trace_bits.back();
			newAzimuthAngle = this->AzimuthAngle(leadBit,endBit);
		}
		return newAzimuthAngle;
	}
	else
	{
		return -1;
	}
}
double TraceLine::GetElevation()
{
	double newElevationAngle = 0;
	if (this->m_parent)
	{
		if (this->m_trace_bits.size()>1)
		{
			//TraceBit pre, cur;
			//TraceBitsType::iterator it = this->m_trace_bits.begin();
			//pre = *it;
			//it++;
			//cur = *it;
			//newElevationAngle = this->ElevationAngle(pre,cur);
			
			TraceBit leadBit = this->m_trace_bits.front();
			TraceBit endBit = this->m_trace_bits.back();
			newElevationAngle = this->ElevationAngle(leadBit,endBit);
		}
		return newElevationAngle;
	}
	else
	{
		return -1;
	}
}
