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

#include "ftl2d_Tracer.h"

#define DEBUG 0

Tracer::Tracer(TraceConfig *config)	{
	VesselContainer.reserve(2000);

}

Tracer::~Tracer()	{
	//delete VesselContainer
	std::vector<Vessel *>::iterator it;
	for (it = VesselContainer.begin(); it != VesselContainer.end(); ++it)
	         delete (*it);
//	 std::cout << "Trace DELETED" << std::endl;	
}

void Tracer::Run(ImageType::Pointer &im, SeedContainer2D *sc, TraceConfig *config)	{

	//TraceContainer *TC = new TraceContainer(config);

	for (unsigned int i = 0; i <sc->getNumberOfSeeds(); i++)	{
	//for (unsigned int i = 6; i < 8; i++)	{
		
		Segment2D *s = new Segment2D(config);   //index  = 0
		s->CopySeedInformation(sc->getSeed(i));
		
		//Print the seed and Segment
		#if DEBUG == 2
			std::cout << "*****************************************************************************"<< std::endl;
			std::cout << "Seed Point  "<< std::endl;
			sc->getSeed(i)->PrintSelf();
			s->PrintSelf(2);
		#endif
		
		Segment2D* s_hit = this->HitTest(s);
		//if (!this->HitTest(s, s_hit))	{					
		if (s_hit == NULL)	{	
				
				delete s_hit;
				s->Initialize(im);
				
				#if DEBUG == 2
				std::cout << "Starting Segment [(0)-not valid (1)-valid]   "<<s->IsSeedSegmentValid()<< std::endl;
				s->PrintSelf(2);
				#endif
			
			if (s->IsSeedSegmentValid())	{
					
				Vessel *v = new Vessel(config, this->getNumberOfVessels());
				v->AddSegmentToVessel(s);
				Segment2D *s_prev = s;
				
				while (s_prev->IsSegmentValid() == 1)	{
					Segment2D * s_new = TraceStep(im, s_prev, config);
						#if DEBUG == 2
						std::cout<<"Segment ID " << s_new->getID() << " FWD "<<"Vessel size " << v->getNumberOfSegments()<<std::endl;
						s_new->PrintSelf(2);
						std::cin.get();
						#endif
					if( !v->AddSegmentToVesselIfValid(s_new) )
						break;
					
					Segment2D* s_hit = this->HitTest(s_new);	
					//if (this->HitTest(s_new, s_hit))	{
					if (s_hit != NULL)	{	
							v->AddSegmentToVessel( s_hit);
							v->SetEndType(22);	//The vessel is HIT
							#if DEBUG == 2
								s_hit->PrintSelf(2);
								s_new->PrintSelf(2);
								v->PrintSelf(); std::cin.get();
							#endif
							break;
						}
					else	{
						delete s_hit;
					}
					
					s_prev = s_new;
				}

				s_prev = s;
				s_prev->flipDirection(v->getNumberOfSegments());
				v->flipDirection();

				

				while (s_prev->IsSegmentValid() == 1)	{
					Segment2D *s_new = TraceStep(im, s_prev, config);
						#if DEBUG == 2
						std::cout<<"Segment ID " << s_new->getID() << " REV "<<"Vessel size " << v->getNumberOfSegments()<<std::endl;
						s_new->PrintSelf(2);
						std::cin.get();
						#endif
					if( !v->AddSegmentToVesselIfValid(s_new) )
						break;
					Segment2D* s_hit = this->HitTest(s_new);	
					if (s_hit != NULL)	{
							v->AddSegmentToVessel( s_hit);
							v->SetEndType(22);	//The vessel is HIT
							#if DEBUG == 2
								s_hit->PrintSelf(2);
								s_new->PrintSelf(2);
								v->PrintSelf(); std::cin.get();
							#endif
							break;
						}
					else	{
						delete s_hit;	
					}
				s_prev = s_new;
				}
				
				if (v->IsVesselValid()) {
					#if DEBUG >= 1
						v->PrintSelf(); std::cin.get();
					#endif
					
					this->AddVessel(v);
					
					#if DEBUG == 0
						std::cout << i << "/"<<sc->getNumberOfSeeds()<<"\r";
					#endif
				}
			}
		}
		
	}
}


Segment2D* Tracer::TraceStep(ImageType::Pointer &im, Segment2D* s_prev, TraceConfig *config){
	Segment2D* s = new Segment2D(config);
	s->advance( s_prev );
	s->fitSuperEllipsoid(im);
	s->SegmentValidCheck(); //implements the segment level rules
	return s;
}


void Tracer::ApplyRules()	{

}

void Tracer::AddBranchPoint(Segment2D* s1, Segment2D* s2) {
	BranchPoint *b = new BranchPoint();
	b->AddBranchPoint(s1, s2);
	this->BranchPointContainer.push_back(b);
}

Segment2D* Tracer::HitTest( Segment2D *s) {
		
	std::vector<Vessel *>::iterator it;

	Vect2 mu, bdry1, bdry2, offset;

	offset[0] =  sin(s->getq())*s->geta1();
	offset[1] =  cos(s->getq())*s->geta1();
	mu = s->getmu();
	bdry1 = mu + offset;
	bdry2 = mu - offset;
	
			
	for (it = this->VesselContainer.begin(); it !=this->VesselContainer.end(); ++it)	{
		for (unsigned int j = 0; j < (*it)->getNumberOfSegments(); j++)	{
//			if (IsPointInsideSE(bdry1, (*it)->getSegment(j)) || IsPointInsideSE(bdry2, (*it)->getSegment(j)) || IsPointInsideSE(mu, (*it)->getSegment(j)))	{
			if (IsPointInsideSE(mu, (*it)->getSegment(j)))	{
				//std::cout<< "HIT vessel " << (*it)->getID()<< "Segment" << (*it)->getSegment(j)->getID() << std::endl;
				this->AddBranchPoint((*it)->getSegment(j), s);
				Segment2D *s_hit = new Segment2D((*it)->getSegment(j));
				return s_hit;
			}
		}
	}
	
	return NULL;
}

bool Tracer::IsPointInsideSE(Vect2 &p, Segment2D *s)	{
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

bool Tracer::getTraces(Vect2& p, unsigned int vs, unsigned int seg)	{
	
	if (vs < this->VesselContainer.size())	{
		if (seg < this->VesselContainer[vs]->getNumberOfSegments())	{
			p = this->VesselContainer[vs]->getSegment(seg)->getmu();
			return 1;
		}
	}
	return 0;
}


void Tracer::ReportXML (TraceConfig *config)	{
	
	
	xmlDocPtr doc = NULL;       /* document pointer */
	xmlNodePtr root_node = NULL, vessel_node = NULL, segment_node = NULL, node2 = NULL, node3 = NULL;/* node pointers */
	xmlDtdPtr dtd = NULL;       /* DTD pointer */
	
	std::string str;
	
	LIBXML_TEST_VERSION;
	
	/* 
	 * Creates a new document, a node and set it as a root node
	 */
	doc = xmlNewDoc(BAD_CAST "1.0");
	root_node = xmlNewNode(NULL, BAD_CAST "Tracing2D");
	str = "SuperEllipsoid2D";
	xmlNewProp(root_node, BAD_CAST "Algorihm", BAD_CAST str.c_str());
	str = ToString((double)config->getGridSpacing());
	xmlNewProp(root_node, BAD_CAST "InputFileName", BAD_CAST config->getInputFileName().c_str());
	str = ToString((double)config->getGridSpacing());
	xmlNewProp(root_node, BAD_CAST "GridSpacing", BAD_CAST str.c_str());
	str = ToString((double)config->getInitSeedSize());
	xmlNewProp(root_node, BAD_CAST "InitialSeedSize", BAD_CAST str.c_str());
	str = ToString((double)config->getSeedIntensityThreshold());
	xmlNewProp(root_node, BAD_CAST "SeedIntensityThreshold", BAD_CAST str.c_str());
	str = ToString((double)config->getStepRatio());
	xmlNewProp(root_node, BAD_CAST "StepRatio", BAD_CAST str.c_str());
	str = ToString((double)config->getAspectRatio());
	xmlNewProp(root_node, BAD_CAST "AspectRatio", BAD_CAST str.c_str());
	str = ToString((double)config->getMinimumVesselWidth());
	xmlNewProp(root_node, BAD_CAST "MinimumVesselWidth", BAD_CAST str.c_str());
	str = ToString((double)config->getMaximumVesselWidth());
	xmlNewProp(root_node, BAD_CAST "MaximumVesselWidth", BAD_CAST str.c_str());
	str = ToString((double)config->getPROP());
	xmlNewProp(root_node, BAD_CAST "Sensitivity", BAD_CAST str.c_str());
	xmlDocSetRootElement(doc, root_node);

	std::vector<Vessel *>::iterator it;
	for (it = this->VesselContainer.begin(); it !=this->VesselContainer.end(); ++it)	{
		
		vessel_node = xmlNewChild(root_node, NULL, BAD_CAST "TraceLine", NULL);
		str = ToString((double)(*it)->getID());
		xmlNewProp(vessel_node, BAD_CAST "ID", BAD_CAST  str.c_str());
		
		for (unsigned int j = 0; j < (*it)->getNumberOfSegments(); j++)	{
			segment_node = xmlNewChild(vessel_node, NULL, BAD_CAST "TraceBit", NULL);
			str = ToString((double)(*it)->getSegment(j)->getID());
			xmlNewProp(segment_node, BAD_CAST "ID", BAD_CAST  str.c_str());
			
			//node2 = xmlNewChild(segment_node, NULL, BAD_CAST "MU", NULL);
			str = ToString((double) (*it)->getSegment(j)->getmu()[0]);    
			xmlNewProp(segment_node, BAD_CAST "x", BAD_CAST str.c_str());
			
			str = ToString((double) (*it)->getSegment(j)->getmu()[1]);    
			xmlNewProp(segment_node, BAD_CAST "y", BAD_CAST str.c_str());
			
			str = ToString((double) (*it)->getSegment(j)->geta1());    
			xmlNewProp(segment_node, BAD_CAST "a1", BAD_CAST str.c_str());
			
			str = ToString((double) (*it)->getSegment(j)->geta2());    
			xmlNewProp(segment_node, BAD_CAST "a2", BAD_CAST str.c_str());
			
			str = ToString((double) (*it)->getSegment(j)->gete1());    
			xmlNewProp(segment_node, BAD_CAST "e1", BAD_CAST str.c_str());
			
			str = ToString((double) (*it)->getSegment(j)->getq());    
			xmlNewProp(segment_node, BAD_CAST "q", BAD_CAST str.c_str());
			
			str = ToString((double) (*it)->getSegment(j)->getEndType());    
			xmlNewProp(segment_node, BAD_CAST "EndType", BAD_CAST str.c_str());
			
			str = ToString((double) (*it)->getSegment(j)->getForegdIntensity());    
			xmlNewProp(segment_node, BAD_CAST "ForeGroundIntensity", BAD_CAST str.c_str());
			
			str = ToString((double) (*it)->getSegment(j)->getBackgdIntensity());    
			xmlNewProp(segment_node, BAD_CAST "BackGroundIntensity", BAD_CAST str.c_str());
			
		}
		
	}
	
	xmlSaveFormatFileEnc(config->getOutputFileName().c_str(), doc, MY_ENCODING, 1);
		
		  /*free the document */
	xmlFreeDoc(doc);
		  
	xmlCleanupParser();

}



