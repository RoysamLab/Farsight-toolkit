////////////////////////////////////////////////////////////////////////////
//Trace points are stored as nodes in an ordered data structure.
//Entire traces are stored by id, and neighboring trace points are encoded 
//within each node.  Class contains some support functions to determine
//validity of trace by detecting intersection with self or other traces, or
//detecting trace discontinuities.

#include "NodeContainer3D.h"
#include "myDebug.h"

#define SHOWTRACEDIG 0

NodeContainer3D::NodeContainer3D()	{
	NodeList.reserve(10000);
}

NodeContainer3D::~NodeContainer3D()	{
	NodeContainerType::iterator it;
	for (it = NodeList.begin(); it != NodeList.end(); ++it)	{
			delete (*it);
	}
}

//retrieve a trace point by its index
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

//add trace point to end of trace
void NodeContainer3D::AddNode(TVessel* seg)	{
	NodeList.push_back(seg);
}

//add trace point to beginning of trace
void NodeContainer3D::AddFront(TVessel* seg)	{
	NodeList.insert( NodeList.begin(), seg );
}

//detect intersection between candidate fit
//and current set of traces
long NodeContainer3D::HitTest (TVessel *s ) {
	Trace* tr = NULL;
	return(HitTest (s, tr));
}

//do the actual hit testing, trace pointer can be
//used to specify the current trace, which generally
//requires a bit of special handling to detect
//self-intersection.
long NodeContainer3D::HitTest (TVessel* s, Trace* tr)	{
	// hit test on NodeContainer
	double p[6][3];


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
	long int loc = -2;


	NodeContainerType::iterator it;
	for (it = NodeList.begin(); it != NodeList.end(); ++it)	{

		if (tr != NULL) {
			if ((s->TraceID == (*it)->TraceID))	{
				//check the number of chain link < threshold
				long nlink = NumberChainLinks(tr, s, *it);
				// S("Same Trace nlink " << nlink)
				if ( nlink <= 5)	{
					//self intersection but NumChainLink less than threshold - hence ignore
					continue;
				}
			}
		}

		double d = (s->mu[0] - (*it)->mu[0])*(s->mu[0] - (*it)->mu[0]) + \
		(s->mu[1] - (*it)->mu[1])*(s->mu[1] - (*it)->mu[1]) + \
		(s->mu[2] - (*it)->mu[2])*(s->mu[2] - (*it)->mu[2]);
		d = vcl_sqrt(d);

		if( d > 50)	{
			continue;
		}

		double mu2[3] = {0 ,0 ,0};
		double a1 = (*it)->a1;
		double a2 = (*it)->a2;
		double a3 = (*it)->a3/2.0;

		double mu1[3];
		mu1[0] = (*it)->mu[0];
		mu1[1] = (*it)->mu[1];
		mu1[2] = (*it)->mu[2];

		for (unsigned int i=0; i<=5; i++)	{
			mu2[0] = p[i][0];
			mu2[1] = p[i][1];
			mu2[2] = p[i][2];

			double u = (*it)->R1[0]*(mu1[0]-mu2[0]) + (*it)->R1[1]*(mu1[1]-mu2[1]) + (*it)->R1[2]*(mu1[2]-mu2[2]);
			double v = (*it)->R2[0]*(mu1[0]-mu2[0]) + (*it)->R2[1]*(mu1[1]-mu2[1]) + (*it)->R2[2]*(mu1[2]-mu2[2]);
			double w = (*it)->R3[0]*(mu1[0]-mu2[0]) + (*it)->R3[1]*(mu1[1]-mu2[1]) + (*it)->R3[2]*(mu1[2]-mu2[2]);

			double D = vcl_pow(u/a1,2.0)+vcl_pow(v/a2,2.0)+vcl_pow(w/a3,2.0);

			if(D < min_D) {
				 min_D = D;
				 loc = (*it)->ID;
			}

		}
	}

	return loc;
}

//traverse the internal data structure holding the trace
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


//find trace point
bool NodeContainer3D::IsInList(NodeContainerType& tlist,  long id)	{
	NodeContainerType::iterator it;
	
	for (it = tlist.begin(); it != tlist.end(); ++it)	{
		if ((*it)->ID==id)	{
			return 1;
		}
	}
	
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

//propogation direction is given by model fit at trace point
bool NodeContainer3D::getPropagationDirectionR3(long id, Trace * tr) {
//this is a replica of Alex's Prop direction where R3 was used as prop direction

	if (id == tr->NodeAID)	{
		TVessel *seg = getSegment(id);
		tr->dirA[0] = seg->R3[0];
		tr->dirA[1] = seg->R3[1];
		tr->dirA[2] = seg->R3[2];
	}
	else if (id == tr->NodeBID) {
		TVessel *seg = getSegment(id);
		tr->dirB[0] = seg->R3[0];
		tr->dirB[1] = seg->R3[1];
		tr->dirB[2] = seg->R3[2];
	}
	else 	{
		std::cerr << "Direction queried NOT from Terminal nodes " << std::endl;
		return false;
	}
	return true;
}


//interpolate trace direction for neighboring trace points
bool NodeContainer3D::getPropagationDirection(long id, Trace * tr) {
	long stNode, endNode;
	double d[3];
	d[0] = 0.0;  d[1] = 0.0;  d[2] = 0.0;


	if (id == tr->NodeAID)	{
		stNode = tr->NodeAID;
		endNode = tr->NodeBID;
	}
	else if (id == tr->NodeBID) {
		stNode = tr->NodeBID;
		endNode = tr->NodeAID;
	}
	else 	{
		std::cerr << "Direction queried NOT from Terminal nodes " << std::endl;
		return false;
	}

	if (tr->numNodes == 1)	{
		TVessel *seg = getSegment(id);
		d[0] = seg->R3[0];
		d[1] = seg->R3[1];
		d[2] = seg->R3[2];
	}
	else if (tr->numNodes == 2)	{
		TVessel *seg1 = getSegment(stNode);
		TVessel *seg2 = getSegment(endNode);

		d[0] = seg1->mu[0] - seg2->mu[0];
		d[1] = seg1->mu[1] - seg2->mu[1];
		d[2] = seg1->mu[2] - seg2->mu[2];

		d[0] = d[0] / (d[0]+d[1]+d[2]);
		d[1] = d[1] / (d[0]+d[1]+d[2]);
		d[2] = d[2] / (d[0]+d[1]+d[2]);
	}

	else	{

		TVessel *seg1 = getSegment(stNode);
		TVessel *seg2 = NULL, *seg3 = NULL;
		for (unsigned int i=0; i<seg1->numNbrs; i++)	{
			long nt = seg1->NbrID[i];
			if (seg1->TraceID  == tr->TraceID)	{
				seg2 = getSegment(nt);
			}
		}

		for (unsigned int i=0; i<seg2->numNbrs; i++)	{
			long nt = seg2->NbrID[i];
			if ((seg2->TraceID  == tr->TraceID) && (nt != seg1->ID)) 	{
				seg3 = getSegment(nt);
			}
		}

		if((seg1) && (seg2) && (seg3))	{

			double d1[3], d2[3];
			d1[0] = 0.5*(seg1->mu[0] - seg3->mu[0]);
			d1[1] = 0.5*(seg1->mu[1] - seg3->mu[1]);
			d1[2] = 0.5*(seg1->mu[2] - seg3->mu[2]);

			d2[0] = (seg1->mu[0] + seg3->mu[0] - 2.0*seg2->mu[0]);
			d2[1] = (seg1->mu[1] + seg3->mu[1] - 2.0*seg2->mu[1]);
			d2[2] = (seg1->mu[2] + seg3->mu[2] - 2.0*seg2->mu[2]);

			d[0] = d1[0] + 0.5*d2[0];
			d[1] = d1[1] + 0.5*d2[1];
			d[2] = d1[2] + 0.5*d2[2];

			d[0] = d[0] / (d[0]+d[1]+d[2]);
			d[1] = d[1] / (d[0]+d[1]+d[2]);
			d[2] = d[2] / (d[0]+d[1]+d[2]);
		}
	}

	if (id == tr->NodeAID)	{
		tr->dirA[0] = d[0];
		tr->dirA[1] = d[1];
		tr->dirA[2] = d[2];
	}
	else if (id == tr->NodeBID) {
		tr->dirB[0] = d[0];
		tr->dirB[1] = d[1];
		tr->dirB[2] = d[2];
	}
	return true;
}

//check consistency of trace directions
bool NodeContainer3D::IsTraceValid(Trace * tr, TVessel* seg) {

	if (ProjectedVsActualShiftOK(tr, seg) == 1)	{
		return 1;
	}
	else	{
		std::cout  << "Projected vs Actual shift fail" << std::cout;
		return 0;
	}
}

// curve distance from seg1 to seg2
double NodeContainer3D::CurveLength(Trace * tr, TVessel* seg1, TVessel* seg2)	{
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

//count number of edges joining adjacent trace points from the same trace
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
				//S("SegEnd reached still could not find " << seg2->ID)
				break;
			}
		}

		if (t==seg2->ID)	{
			//S("HURRAY found .. breaking " << seg2->ID)
			break;
		}
	}

	if (t==seg2->ID) {
		//S("Returning link-length "<<d)
		return d;
	}
	else	{
		return maxd;
	}
}


//looking for discontinuities in trace direction when traversing nodes
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



void NodeContainer3D::WriteSegmentsToTextFile(std::string fname) {
	std::ofstream ofile, linkfile;
	ofile.open(fname.c_str(), std::ios::out);

	std::string linkfname = fname+std::string("-link.txt");

	linkfile.open(linkfname.c_str(), std::ios::out);


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
}


