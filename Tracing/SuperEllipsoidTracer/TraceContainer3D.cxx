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

#include "TraceContainer3D.h"
#define S(x) std::cout<< x <<std::endl;

#define STEP_FWD 1
#define STEP_BWD -1

TraceContainer3D::TraceContainer3D() {
	TraceList.reserve(1000);
	m_THRESH = 0.2;
	m_minL = 3.0;
	m_Stepsize = 0.5;
	m_Spacing.Fill(1.0);
}

void TraceContainer3D::Configure(TraceConfig::Pointer& config)	{
	m_THRESH = config->getTHRESHOLD();
	m_minL = config->getminContrast();
	m_Stepsize = config->getStepRatio();
	m_AspectRatio = config->getAspectRatio();
	m_Spacing = config->GetSpacing();
}

TraceContainer3D::~TraceContainer3D() {
	TraceContainer::iterator it2;
	std::cout << "Deleting Traces" << std::endl;
	for (it2 = TraceList.begin(); it2 != TraceList.end(); ++it2)	{
		delete (*it2);
	}
}


void TraceContainer3D::ComputeTrace(ImageType3D::Pointer im, Seed2Seg::Pointer sseg) {

	unsigned long tID = 0;
	unsigned long segID = 0;
	double iterations = 25.0;
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
			TVessel *seg1 = tr->Step(seg, im, segID, STEP_FWD, iterations, m_AspectRatio, m_THRESH);
			

			if (!seg1) {
				S("seg 1 returned empty")
				break;
			}

			if (seg1->IsSegmentValid(seg, m_THRESH, m_minL)==0)	{
				seg1->PrintSelf();
				S("Seg Valid Violation")
				break;
			}

			NodeList->AddNode(seg1);
			segID++;

			tr->UpdateTrace(seg, seg1, STEP_FWD);
			NodeList->getPropagationDirectionR3( seg1->ID, seg->ID, tr);

			//seg1->PrintSelf();	tr->PrintSelf();	//NodeList->PrintSelf();	S("Step done ")


			if (NodeList->IsTraceValid(tr, seg1)==0)	{
				S("Trace Valid Violation")
				break;
			}

			long hitSeg = NodeList->HitTest(seg1, tr);
			if (hitSeg > -1)	{
				AddBranchPoint(seg1, NodeList->getSegment(hitSeg));
				S("Hit segment " << hitSeg)
				break;
			}
			seg = seg1;
			std::cout <<">";
		}

		//travel reverse in node container
		seg = NodeList->getSegment(rootSegID);
		tr->Reverse( seg );

		numStep = 0;
		while(numStep < 1000)	{

			numStep++;
			TVessel *seg1 = tr->Step(seg, im, segID, STEP_BWD, iterations, m_AspectRatio, m_THRESH);

			if (!seg1) {
				S("seg 1 returned empty")
				break;
			}

			if (seg1->IsSegmentValid(seg, m_THRESH, m_minL)==0)	{
				seg1->PrintSelf();
				S("Seg Valid Violation")
				break;
			}

			NodeList->AddNode(seg1);
			segID++;

			tr->UpdateTrace(seg, seg1, STEP_BWD);
			NodeList->getPropagationDirectionR3( seg1->ID, seg->ID, tr);

			//seg1->PrintSelf();	tr->PrintSelf();	//NodeList->PrintSelf();	S("Step done ")


			if (NodeList->IsTraceValid(tr, seg1)==0)	{
				S("Trace Valid Violation")
				break;
			}

			long hitSeg = NodeList->HitTest(seg1, tr);
			if (hitSeg > -1)	{
				AddBranchPoint(seg1, NodeList->getSegment(hitSeg));
				S("Hit segment " << hitSeg)
				break;
			}

			seg = seg1;
			std::cout <<"<";
		}

		tr->PrintSelf();
		TraceList.push_back(tr);
		tID++;
	}

}


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


void TraceContainer3D::WriteTraceToXMLFile(std::string SegXMLFname)	{
	NodeList->WriteSegmentsToXMLFile(SegXMLFname);
	NodeList->GenerateStatistics(SegXMLFname , m_Spacing );
}

void TraceContainer3D::AddBranchPoint(TVessel *s1, TVessel *s2)	{
	
	if ((s1->ID) == (s2->ID))	{
		return;
	}
    if ((s1->TraceID) == (s2->TraceID))	{
		return;
	}
	
	s1->numNbrs++;
	s2->numNbrs++;
	s1->NbrID.push_back(s2->ID);
	s2->NbrID.push_back(s1->ID);
	std::cout << "Branch Point added between Trace:" << s1->TraceID 
		<< " (Segment:" << s1->ID <<") and  Trace:" << s2->TraceID 
		<< " (Segment:" << s2->ID <<")" << std::endl;
	s1->PrintSelf();
	s2->PrintSelf();
}