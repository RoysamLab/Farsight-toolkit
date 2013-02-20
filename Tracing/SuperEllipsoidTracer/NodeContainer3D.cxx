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

#include "NodeContainer3D.h"
#define S(x) std::cout<< x <<std::endl;
#define SHOWTRACEDIG 0

NodeContainer3D::NodeContainer3D()	{
	NodeList.reserve(10000);
}

NodeContainer3D::~NodeContainer3D()	{
	NodeContainerType::iterator it;
	std::cout << "Deleting nodes" << std::endl;
	for (it = NodeList.begin(); it != NodeList.end(); ++it)	{
			delete (*it);
	}
}

TVessel* NodeContainer3D::getSegment(long i) {
	// change this to match ID if required
	NodeContainerType::iterator it;
	for (it = NodeList.begin(); it != NodeList.end(); ++it)	{
		if((*it)->ID == i)	{
			return(*it);
		}
	}
	return NULL;
}

void NodeContainer3D::AddNode(TVessel* seg)	{
	NodeList.push_back(seg);
}

void NodeContainer3D::AddFront(TVessel* seg)	{
	NodeList.insert( NodeList.begin(), seg );
}

long NodeContainer3D::HitTest (TVessel *s ) {
	Trace* tr = NULL;
	return(HitTest (s, tr));
}

long NodeContainer3D::HitTest (TVessel* s, Trace* tr)	{
	// hit test on NodeContainer
	double p[7][3];

	double fact = 0.7;
	//along R1
	p[0][0] = s->mu[0] + fact*(s->a1*s->R1[0]);
	p[0][1] = s->mu[1] + fact*(s->a1*s->R1[1]);
	p[0][2] = s->mu[2] + fact*(s->a1*s->R1[2]);

	p[1][0] = s->mu[0] - fact*(s->a1*s->R1[0]);
	p[1][1] = s->mu[1] - fact*(s->a1*s->R1[1]);
	p[1][2] = s->mu[2] - fact*(s->a1*s->R1[2]);

	//along R2
	p[2][0] = s->mu[0] + fact*(s->a2*s->R2[0]);
	p[2][1] = s->mu[1] + fact*(s->a2*s->R2[1]);
	p[2][2] = s->mu[2] + fact*(s->a2*s->R2[2]);

	p[3][0] = s->mu[0] - fact*(s->a2*s->R2[0]);
	p[3][1] = s->mu[1] - fact*(s->a2*s->R2[1]);
	p[3][2] = s->mu[2] - fact*(s->a2*s->R2[2]);

	//along R3
	p[4][0] = s->mu[0] + fact*(s->a3*s->R3[0]);
	p[4][1] = s->mu[1] + fact*(s->a3*s->R3[1]);
	p[4][2] = s->mu[2] + fact*(s->a3*s->R3[2]);

	p[5][0] = s->mu[0] - fact*(s->a3*s->R3[0]);
	p[5][1] = s->mu[1] - fact*(s->a3*s->R3[1]);
	p[5][2] = s->mu[2] - fact*(s->a3*s->R3[2]);
	
	//center
	p[6][0] = s->mu[0];
	p[6][1] = s->mu[1];
	p[6][2] = s->mu[2];

	double min_D = 1.0;
	long int loc = -2;


	NodeContainerType::iterator it;
	for (it = NodeList.begin(); it != NodeList.end(); ++it)	{

		//allow hitting the same trace
		// S("Checking "<< (*it)->ID )

		if (tr != NULL) {
			// S("GOT Trace value in Hitest " << tr->TraceID)
			if ((s->TraceID == (*it)->TraceID))	{
				//check the number of chain link < threshold
				long nlink = NumberChainLinks(tr, s, *it);
				// S("Same Trace nlink " << nlink)
				if ( nlink <= 8)	{
					//self intersection but NumChainLink less than threshold - hence ignore
					continue;
				}
			}
		}
		else {
			// S("Using default NULL value in Hitest")
		}

		// S("Consider " << (*it)->ID <<" for Distance test")

		double d = (s->mu[0] - (*it)->mu[0])*(s->mu[0] - (*it)->mu[0]) + \
		(s->mu[1] - (*it)->mu[1])*(s->mu[1] - (*it)->mu[1]) + \
		(s->mu[2] - (*it)->mu[2])*(s->mu[2] - (*it)->mu[2]);
		d = vcl_sqrt(d);

		if( d > 50)	{
			continue;
		}

		// S("Consider " << (*it)->ID <<" for Intersection test")

		double mu2[3] = {0 ,0 ,0};
		double a1 = (*it)->a1;
		double a2 = (*it)->a2;
		double a3 = (*it)->a3/2.0;

		double mu1[3];
		mu1[0] = (*it)->mu[0];
		mu1[1] = (*it)->mu[1];
		mu1[2] = (*it)->mu[2];


		// S("mu1 = " << mu1[0] << "," << mu1[1] << "," << mu1[2] )
		// S("a = " << a1 << "," << a2 << "," << a3 )
		for (unsigned int i=0; i<=6; i++)	{
			mu2[0] = p[i][0];
			mu2[1] = p[i][1];
			mu2[2] = p[i][2];

			// S("mu2 = " << mu2[0] << "," << mu2[1] << "," << mu2[2] )

			double u = (*it)->R1[0]*(mu1[0]-mu2[0]) + (*it)->R1[1]*(mu1[1]-mu2[1]) + (*it)->R1[2]*(mu1[2]-mu2[2]);
			double v = (*it)->R2[0]*(mu1[0]-mu2[0]) + (*it)->R2[1]*(mu1[1]-mu2[1]) + (*it)->R2[2]*(mu1[2]-mu2[2]);
			double w = (*it)->R3[0]*(mu1[0]-mu2[0]) + (*it)->R3[1]*(mu1[1]-mu2[1]) + (*it)->R3[2]*(mu1[2]-mu2[2]);

			// S("uvw = " << u << "," << v << "," << w )

			double D = vcl_pow(u/a1,2.0)+vcl_pow(v/a2,2.0)+vcl_pow(w/a3,2.0);
			// S("D found for " << i << " = " << D )
			if(D < min_D) {
				 min_D = D;
				 loc = (*it)->ID;
			}

		}
	}

	 if ((loc > -1) && (tr != NULL))	{
		S("Hit detected between "<< s->TraceID << " and " << getSegment(loc)->TraceID << " at " << s->ID << "x" << loc)
//	 	getSegment(loc)->PrintSelf();
//	 	s->PrintSelf();
	 }
	return loc;
}


bool NodeContainer3D::getFullTrack( NodeContainerType &tracelist, Trace * tr, long maxNode = 1000)	{

	tracelist.push_back(getSegment(tr->NodeAID));
	long t = tr->NodeAID;
	TVessel *seg, *nseg;
	unsigned char flag;

	if (SHOWTRACEDIG) {std::cout << "Digging "<< tr->TraceID << " [ " << tr->NodeAID << " , " << tr->NodeBID << "]" <<std::endl;}

	long Cnt = 0;
	while (t!= tr->NodeBID)	{
		seg = getSegment(t);
		flag = 0;
		for(unsigned int i=0; i<seg->numNbrs; i++){
			long nt = seg->NbrID[i];
			nseg = getSegment(nt);
			nt = nseg->ID;
			if ((nt == tr->TraceID) && (!IsInList(tracelist, nt)))	{
				tracelist.push_back(nseg);
				flag = 1;
				t = nt;
				if (SHOWTRACEDIG)	{	std::cout <<"{" << nt << "}, ";	}
				break;
			}
			else	{
				if (SHOWTRACEDIG)	{	std::cout << nt << ", ";		}
			}
		}
		if ((flag==0) && (t!=tr->NodeBID))	{
			std::cout << "missing link" << std::endl;
			return 0;
		}
		Cnt++;
		if (Cnt > maxNode)	{
			return 1;
		}

	}
	return 1;
}



bool NodeContainer3D::IsInList(NodeContainerType& tlist,  long id)	{
	NodeContainerType::iterator it;
	//std::cout << "Checking Tracelist for  " << id << " in " << tlist.size() << std::endl;
	for (it = tlist.begin(); it != tlist.end(); ++it)	{
		//std::cout << " " << (*it)->ID ;
		if ((*it)->ID==id)	{
			//std::cout << "----FOUND return 1" << std::endl;
			return 1;
		}
	}
	//std::cout << "----NOT FOUND return 0" << std::endl;
	 return 0;

}


void NodeContainer3D::PrintSelf()	{
	std::cout << std::endl << "Number of Nodes " << NodeList.size() << std::endl;
		NodeContainerType::iterator it;
		for (it = NodeList.begin(); it != NodeList.end(); ++it)	{
			std::cout << "ID:"<< (*it)->ID << " @ [" << (*it)->mu[0] <<"," << (*it)->mu[1] << "," << (*it)->mu[2] << "]   Nbrs: " ;
			for (unsigned int i=0; i<(*it)->numNbrs; i++)	{
				std::cout << (*it)->NbrID[i] <<" ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
}


void NodeContainer3D::PrintFullTrack(Trace * tr)	{

	std::cout << "Digging "<< tr->TraceID << " [ " << tr->NodeAID << " , " << tr->NodeBID << "]" <<std::endl;
	NodeContainerType tracelist;
	bool ret = getFullTrack( tracelist, tr );
	if (ret)	{
		std::cout << "Length of Track" << tracelist.size() << std::endl;
		NodeContainerType::iterator it;
		for(it = tracelist.begin(); it != tracelist.end(); ++it)	{
			std::cout << "ID:"<< (*it)->ID << " @ [" << (*it)->mu[0] <<"," << (*it)->mu[1] << "," << (*it)->mu[2] << "]   Nbrs: " ;
			for (unsigned int i=0; i<(*it)->numNbrs; i++)	{
				std::cout << (*it)->NbrID[i] <<" ";
			}
		std::cout << std::endl;
		}
	}
}


bool NodeContainer3D::getPropagationDirectionR3(long id, long id1, Trace * tr) {
//this is a replica of Alex's Prop direction where R3 was used as prop direction

	if (id == tr->NodeAID)	{
		TVessel *seg = getSegment(id);
		TVessel *seg1 = getSegment(id1);
		tr->dirA[0] = seg->R3[0] + seg1->R3[0] ;
		tr->dirA[1] = seg->R3[1] + seg1->R3[1];
		tr->dirA[2] = seg->R3[2] + seg1->R3[2];
		double nrm = vcl_sqrt(tr->dirA[0]*tr->dirA[0] + tr->dirA[1]*tr->dirA[1] + tr->dirA[2]*tr->dirA[2]) + 0.001f;
		tr->dirA[0] /= nrm;
		tr->dirA[1] /= nrm;
		tr->dirA[2] /= nrm;
	}
	else if (id == tr->NodeBID) {
		TVessel *seg = getSegment(id);
		TVessel *seg1 = getSegment(id1);
		tr->dirB[0] = seg->R3[0] + seg1->R3[0];
		tr->dirB[1] = seg->R3[1] + seg1->R3[1];
		tr->dirB[2] = seg->R3[2] + seg1->R3[2];
		double nrm = vcl_sqrt(tr->dirB[0]*tr->dirB[0] + tr->dirB[1]*tr->dirB[1] + tr->dirB[2]*tr->dirB[2]) + 0.001f;
		tr->dirB[0] /= nrm;
		tr->dirB[1] /= nrm;
		tr->dirB[2] /= nrm;
	}
	else 	{
		std::cerr << "Direction queried NOT from Terminal nodes " << std::endl;
		return false;
	}
	return true;
}


bool NodeContainer3D::IsTraceValid(Trace * tr, TVessel* seg) {

	if (ProjectedVsActualShiftOK(tr, seg) == 1)	{
		return 1;
	}
	else	{
		std::cout  << "Projected vs Actual shift fail" << std::cout;
		return 0;
	}
}


double NodeContainer3D::CurveLength(Trace * tr, TVessel* seg1, TVessel* seg2)	{
	// curve distance from seg1 to seg2
	double maxd = 10000.0;

	//first do a sanity check All IDs match
	if ((seg1->TraceID != tr->TraceID) || (seg2->TraceID != tr->TraceID))	{
		std::cout << "Incorrect Segments IDs donot match" << std::endl;
		return maxd;
	}

	double d = 0;
	unsigned char flag = 0;
	// for each neighbor of the root segment
	for (unsigned int i=0; i<seg1->numNbrs; i++)	{
		long t = seg1->NbrID[i];
		TVessel* s = getSegment(t);
		NodeContainerType  tracelist;
		tracelist.push_back(seg1);
		tracelist.push_back(s);

		d += vcl_sqrt((seg1->mu[0]-s->mu[0])*(seg1->mu[0]-s->mu[0]) +  \
					(seg1->mu[1]-s->mu[1])*(seg1->mu[1]-s->mu[1]) + \
					(seg1->mu[2]-s->mu[2])*(seg1->mu[2]-s->mu[2])) ;

		long Cnt = 0;
		//while the leaf segment is not reached
		while (t!= seg2->ID)	{
			flag = 0;
			for (unsigned int j=0; j<s->numNbrs; j++)	{
				long nt1 = s->NbrID[j];
				TVessel* nseg1 = getSegment(nt1);
				if ((nseg1->TraceID == tr->TraceID) && (!IsInList(tracelist, nt1)))	{

					d += vcl_sqrt(  (nseg1->mu[0]-s->mu[0])*(nseg1->mu[0]-s->mu[0]) +  \
							    	(nseg1->mu[1]-s->mu[1])*(nseg1->mu[1]-s->mu[1]) + \
								    (nseg1->mu[2]-s->mu[2])*(nseg1->mu[2]-s->mu[2]) );
					flag = 1;
					t = nt1;
					s = nseg1;
					tracelist.push_back(s);
					break;
				}
			}

			if (flag==0)	{
				std::cout << "missing link" << std::endl;
				return maxd;
			}
			Cnt++;
			if (Cnt > 1000)	{
				std::cout << "Infinite loop" << std::endl;
				return maxd;
			}

			if ((t==tr->NodeAID) || (t==tr->NodeBID))	{
				break;
			}
		}

		if (t!=seg2->ID)	{
			return maxd;
		}
		else 	{
			return d;
		}
	}
	return maxd;
}

unsigned int NodeContainer3D::NumberChainLinks(Trace * tr, TVessel* seg1, TVessel* seg2)	{
	// link distance from seg1 to seg2
	unsigned int maxd = 10000;

	//first do a sanity check All IDs match
	if (seg1->ID == seg2->ID)	{
		return 0;
	}

	if ((seg1->TraceID != tr->TraceID) || (seg2->TraceID != tr->TraceID))	{
		std::cout << "Incorrect Segments IDs donot match" << std::endl;
		return maxd;
	}

	unsigned int d = 0;
	unsigned char flag = 0;
	long t = 0;
	// for each neighbor of the root segment
	for (unsigned int i=0; i<seg1->numNbrs; i++)	{
		t = seg1->NbrID[i];
		TVessel* s = getSegment(t);
		NodeContainerType  tracelist;
		tracelist.push_back(seg1);
		tracelist.push_back(s);

		d = 1;
		long Cnt = 0;
		//while the leaf segment is not reached
		while (t!= seg2->ID)	{
			flag = 0;
			for (unsigned int j=0; j<s->numNbrs; j++)	{
				long nt1 = s->NbrID[j];
				TVessel* nseg1 = getSegment(nt1);

				if ((nseg1->TraceID == tr->TraceID) && (!IsInList(tracelist, nt1)))	{

					d += 1;
					flag = 1;
					t = nt1;
					s = nseg1;
					tracelist.push_back(s);
					break;
				}
			}

			if (flag==0)	{
				std::cout << "missing link" << std::endl;
				return maxd;
			}
			Cnt++;
			if (Cnt > 1000)	{
				std::cout << "Infinite loop" << std::endl;
				return maxd;
			}

			if ((t==tr->NodeAID) || (t==tr->NodeBID))	{
				break;
			}
		}

		if (t==seg2->ID)	{
			break;
		}
	}

	if (t==seg2->ID) {
		return d;
	}
	else	{
		return maxd;
	}
}



bool NodeContainer3D::ProjectedVsActualShiftOK(Trace * tr, TVessel* seg) {

	// sanity check
	if (seg->TraceID != tr->TraceID)	{
		std::cout << "Incorrect Segments ID in ProjectedVsActualShiftOK" << std::endl;
		return 0;
	}

	long endID = -1;
	if (seg->ID == tr->NodeAID)	{
		endID = tr->NodeBID;
	}
	else if (seg->ID == tr->NodeBID)	{
		endID = tr->NodeAID;
	}
	else {
		std::cout << "End nodes not provided in ProjectedVsActualShiftOK" << std::endl;
		return 0;
	}

	double projected = 0.0;
	double actual = 0.0;
	unsigned char flag = 0;
	long t = seg->ID;
	unsigned int Cnt = 0;

	NodeContainerType  tracelist;
	tracelist.push_back(seg);

	while ((t!= endID) && (Cnt < 5))	{
		flag = 0;
		for (unsigned int i=0; i<seg->numNbrs; i++)	{
			long nt = seg->NbrID[i];
			TVessel* nseg = getSegment(nt);
			if ((nseg->TraceID == tr->TraceID) && (!IsInList(tracelist, nt)))	{

				projected += 0.5*vnl_math_max(nseg->a1, nseg->a2);
				actual += vcl_sqrt( (seg->mu[0]-nseg->mu[0])*(seg->mu[0]-nseg->mu[0]) + \
							        (seg->mu[1]-nseg->mu[1])*(seg->mu[1]-nseg->mu[1]) + \
							        (seg->mu[2]-nseg->mu[2])*(seg->mu[2]-nseg->mu[2]) );

				flag = 1;
				t = nt;
				seg = nseg;
				tracelist.push_back(seg);
				break;
			}
		}

		if (flag==0)	{
			std::cout << "missing link" << std::endl;
			return 0;
		}
		Cnt++;
	}

	if (0.25*projected > actual)	{
		return 0;
	}
	else {
		return 1;
	}
}



void NodeContainer3D::WriteSegmentsToTextFile(std::string& fname) {
	std::ofstream ofile, linkfile;
	ofile.open(fname.c_str(), std::ios::out);

	std::string linkfname = fname+std::string("-link.txt");

	linkfile.open(linkfname.c_str(), std::ios::out);

	std::cout << "Writing output to textfile. " << std::endl << "Filename: " << fname << "  Linkfile: " <<  linkfname << std::endl;
	NodeContainerType::iterator it;
	for (it = NodeList.begin(); it != NodeList.end(); ++it)	{

		TVessel* seg = (*it);
		ofile << "ID:" << seg->ID << " TraceID:" << seg->TraceID << \
		" mu:[" << seg->mu[0] << "," << seg->mu[1] << "," <<  seg->mu[2] << \
		"] A:[" << seg->a1 << "," << seg->a2 << "," <<  seg->a3 << \
		"] Q:[" << seg->q1[0]<< ","<< seg->q1[1] << ","<< seg->q1[2] << "," << seg->q1[3] << \
		"] e1:[" << seg->e1 <<\
		"] R1:[" <<seg->R1[0] << "," << seg->R2[0] << "," <<  seg->R3[0] << \
		"] R2:[" << seg->R1[1] << "," << seg->R2[1] << "," <<  seg->R3[1] << \
		"] R3:[" << seg->R1[2] << "," << seg->R2[2] << "," <<  seg->R3[2] << \
		"] Foregd:" << seg->f <<\
		" Backgd:" << seg->b <<\
		" Lhood:" <<  seg->L <<\
		" MAD:" <<  seg->MAD <<\
		" NumNeighbors:" << seg->numNbrs <<std::endl;

		linkfile << seg->ID << " ";
		for (unsigned int i=0; i<seg->numNbrs; i++)	{
			linkfile << seg->NbrID[i] <<" ";
		}
		linkfile << std::endl;
	}

	ofile.close();
	linkfile.close();
	std::cout << "...Complete."<<std::endl;
}



void NodeContainer3D::WriteSegmentsToXMLFile(std::string& fname) {

	std::cout << "Writing output to XMLfile. " << std::endl << "Filename: " << fname << std::endl;
	TiXmlDocument doc;
	TiXmlDeclaration * decl = new TiXmlDeclaration( "1.0", "", "" );
	doc.LinkEndChild( decl );
	

	NodeContainerType::iterator it;
	for (it = NodeList.begin(); it != NodeList.end(); ++it)	{

		TiXmlElement * element = new TiXmlElement( "Superellipse" );
		element->SetAttribute("ID", (*it)->ID );
		element->SetAttribute("TraceID", (*it)->TraceID );
		element->SetDoubleAttribute("x", (*it)->mu[0] );
		element->SetDoubleAttribute("y", (*it)->mu[1] );
		element->SetDoubleAttribute("z", (*it)->mu[2] );
		element->SetDoubleAttribute("a1", (*it)->a1 );
		element->SetDoubleAttribute("a2", (*it)->a2 );
		element->SetDoubleAttribute("a3", (*it)->a3 );
		element->SetDoubleAttribute("q1", (*it)->q1[0] );
		element->SetDoubleAttribute("q2", (*it)->q1[1]);
		element->SetDoubleAttribute("q3", (*it)->q1[2]);
		element->SetDoubleAttribute("q4", (*it)->q1[3]);
		element->SetDoubleAttribute("R11", (*it)->R1[0]);
		element->SetDoubleAttribute("R12", (*it)->R1[1]);
		element->SetDoubleAttribute("R13", (*it)->R1[2]);
		element->SetDoubleAttribute("R21", (*it)->R2[0]);
		element->SetDoubleAttribute("R22", (*it)->R2[1]);
		element->SetDoubleAttribute("R23", (*it)->R2[2]);
		element->SetDoubleAttribute("R31", (*it)->R3[0]);
		element->SetDoubleAttribute("R32", (*it)->R3[1]);
		element->SetDoubleAttribute("R33", (*it)->R3[2]);
		element->SetDoubleAttribute("Foregd", (*it)->f);
		element->SetDoubleAttribute("Backgd", (*it)->b);
		element->SetDoubleAttribute("Lhood", (*it)->L);
		element->SetDoubleAttribute("MAD", (*it)->MAD);

		for (unsigned int i=0; i<(*it)->numNbrs; i++) {
			TiXmlElement * nbr = new TiXmlElement( "Neighbors" );
			nbr->SetAttribute("ID", (*it)->NbrID[i]);
			element->LinkEndChild(nbr);
		}
	doc.LinkEndChild( element );
	}
	
	
	doc.SaveFile( fname.c_str() );
	std::cout << "...Complete."<<std::endl;
}


void NodeContainer3D::GenerateStatistics(std::string& fname, itk::FixedArray<double,3> spacing) {
	fname.replace(fname.length()-4,11 ,"_Stats.txt");
	std::cout << "Writing Statistics file. " << std::endl << "Filename: " << fname << std::endl;

	NodeContainerType::iterator it;
	double totLen = 0.0;
	unsigned long numBranchPts = 0;

	for (it = NodeList.begin(); it != NodeList.end(); ++it)	{
		double x = (*it)->mu[0];
		double y = (*it)->mu[1];
		double z = (*it)->mu[2];
		
		if ((*it)->numNbrs > 2)	{
			numBranchPts++;
		}

		for (unsigned int i=0; i<(*it)->numNbrs; i++) {
			long nt = (*it)->NbrID[i];
			if((*it)->NbrID[i] < (*it)->ID)	{
				TVessel* nseg = getSegment(nt);
				double x1 = nseg->mu[0];
				double y1 = nseg->mu[1];
				double z1 = nseg->mu[2];

				double dx = (x - x1)*spacing[0];
				double dy = (y - y1)*spacing[1];
				double dz = (z - z1)*spacing[2];
				totLen += vcl_sqrt( dx*dx + dy*dy + dz*dz);
			}
		}
	}
	std::ofstream ofile;
	ofile.open(fname.c_str(), std::ios::out);
	ofile << "Image spacing used: "  << spacing << " microns" << std::endl;
	ofile << "Total length of traces: "  << totLen << " microns" << std::endl;
	ofile << "Total number of branches: "  << numBranchPts << std::endl;
	ofile.close();
	std::cout << "...Complete."<<std::endl;
}


