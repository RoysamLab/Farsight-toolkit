///////////////////////////////////////////////////////////
//TraceContainer3D performs the actual trace calculation.
//Models are fit from a set of seeds points and are
//stored in a node list.  

#include "TraceContainer3D.h"

#include "myDebug.h"

#define STEP_FWD 1
#define STEP_BWD -1

TraceContainer3D::TraceContainer3D() {
	TraceList.reserve(1000);
}


TraceContainer3D::~TraceContainer3D() {
	TraceContainer::iterator it2;
	for (it2 = TraceList.begin(); it2 != TraceList.end(); ++it2)	{
		delete (*it2);
	}
}


//Function to compute the actual trace from a seed point. Procedure is to
//trace in two directions 180 degrees apart.  Tracing proceeds until
//boundary is detected, another trace is intersected, or likelihood drops
//below threshold.
void TraceContainer3D::ComputeTrace(ImageType3D::Pointer im, Seed2Seg::Pointer sseg) {

	unsigned long tID = 0;
	unsigned long segID = 0;

	NodeList = NodeContainer3D::New();


	S("Tracer start...")
	for (unsigned int i=0; i< sseg->getNumberOfStartSegments(); i++ )	{

		if (NodeList->HitTest(sseg->getStartSegment(i)) > -1)	{
			S("Failed starting HIT  test")
			continue;
		}

		TVessel *seg = sseg->getStartSegment(i)->CopyToNewSegment();
		seg->ID = segID;
		seg->TraceID = tID;
		NodeList->AddNode(seg);
		segID++;

		Trace* tr = InitiazeTracer(seg, tID);

		//travel forward in node container

		long rootSegID = seg->ID;
		int numStep = 0;
		while(numStep < 1000)	{

			numStep++;
			TVessel *seg1 = tr->Step(seg, im, segID, STEP_FWD);

			if (!seg1) {
				S("seg 1 returned empty")
				break;
			}

			if (seg1->IsSegmentValid(seg)==0)	{
				seg1->PrintSelf();
				S("Seg Valid Violation")
				break;
			}

			NodeList->AddNode(seg1);
			segID++;

			tr->UpdateTrace(seg, seg1, STEP_FWD);
			NodeList->getPropagationDirectionR3( seg1->ID, tr);

			if (NodeList->IsTraceValid(tr, seg1)==0)	{
				S("Trace Valid Violation")
				break;
			}

			long hitSeg;
			if (hitSeg = NodeList->HitTest(seg1, tr) > -1)	{
				S("Hit segment" << hitSeg)
				break;
			}

			seg = seg1;
		}

		//travel reverse in node container
		seg = NodeList->getSegment(rootSegID);
		tr->Reverse( seg );

		numStep = 0;
		while(numStep < 1000)	{

			numStep++;
			TVessel *seg1 = tr->Step(seg, im, segID, STEP_BWD);

			if (!seg1) {
				S("seg 1 returned empty")
				break;
			}

			if (seg1->IsSegmentValid(seg)==0)	{
				seg1->PrintSelf();
				S("Seg Valid Violation")
				break;
			}

			NodeList->AddFront(seg1);
			segID++;

			tr->UpdateTrace(seg, seg1, STEP_BWD);
			NodeList->getPropagationDirectionR3( seg1->ID, tr);

			if (NodeList->IsTraceValid(tr, seg1)==0)	{
				S("Trace Valid Violation")
				break;
			}

			long hitSeg;
			if (hitSeg = NodeList->HitTest(seg1, tr) > -1)	{
				S("Hit segment" << hitSeg)
				break;
			}

			seg = seg1;
		}

		tr->PrintSelf();
		TraceList.push_back(tr);
		tID++;
	}

}

//Output the trace to a text file.
void TraceContainer3D::WriteTraceToTxtFile(std::string SegTxtFname)	{
	NodeList->WriteSegmentsToTextFile(SegTxtFname);
	std::string traceFname = SegTxtFname + std::string("-trace.txt");

	std::ofstream ofile;
	ofile.open(traceFname.c_str(), std::ios::out);

	TraceContainer::iterator it;
	for (it = TraceList.begin(); it !=  TraceList.end(); ++it )	{
		ofile << (*it)->TraceID << " " << (*it)->numNodes << " " << \
		(*it)->NodeAID << " " << (*it)->NodeBID << " " << \
		(*it)->L/(*it)->numNodes << std::endl;
	}
	ofile.close();
}

//Set the initial directions from a seed point.  Also, assign an id to the
//trace.
Trace* TraceContainer3D::InitiazeTracer(TVessel* seg, unsigned long tID)	{
	Trace * tr = new Trace();
	tr->TraceID = tID;
	tr->numNodes = 1;
	tr->L = seg->L;
	tr->NodeAID = seg->ID;
	tr->NodeBID = seg->ID;
	for (int i=0; i<3; i++)	{
		tr->dirA[i] = seg->R3[i];
		tr->dirB[i] = -1*seg->R3[i];
	}
	return tr;
}


