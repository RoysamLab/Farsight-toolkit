//#include <iostream>
//#include <vector>
//#include <algorithm>
//#include "vnl/vnl_math.h"

#include <limits.h>

#include <itkImage.h>
//#include "itkObjectFactory.h"
//#include "itkMacro.h"
//#include "itkLightObject.h"

#include "itkPointSet.h"

#include "SpineConsts.h"
#include "SpineUtils.h"
#include "CommonTypeDefs.h"
#include "TraceSegNode.h"
#include "tinyxml.h"
#include "ImageProc.h"

TraceContainer::TraceContainer(const char*f, ImageDebugger *imdebugger ){
	xmlfname = f;
	xmlvalid    = ReadTraceSegNodeXMLFile();
	//TraceIDMap.reserve(1000);
	PrintSelf();
	imdbg = imdebugger;
}

void TraceSegNode::rotation_quat(TRMatrix & R)
{
    double w = q1[0];
    double x = q1[1];
    double y = q1[2];
    double z = q1[3];

    //double w = cos(s/2.0);
    //double x = v1*sin(s/2.0);
    //double y = v2*sin(s/2.0);
    //double z = v3*sin(s/2.0);

    R.r[0][0] = 1-2*(pow(y,2.0)+pow(z,2.0));
    R.r[1][1] = 1-2*(pow(x,2.0)+pow(z,2.0));
    R.r[2][2] = 1-2*(pow(x,2.0)+pow(y,2.0));

    R.r[0][1] = 2*(x*y-w*z);
    R.r[0][2] = 2*(x*z+w*y);
    R.r[1][2] = 2*(y*z-w*x);

    R.r[2][0] = 2*(x*z-w*y);
    R.r[2][1] = 2*(y*z+w*x);
    R.r[1][0] = 2*(x*y+w*z);

}

void TraceSegNode::transpose_matrix(TRMatrix & Result )
{
    double tmp;

    tmp = Result.r[1][0];
    Result.r[1][0] = Result.r[0][1];
    Result.r[0][1] = tmp;

    tmp = Result.r[2][0];
    Result.r[2][0] = Result.r[0][2];
    Result.r[0][2] = tmp;

    tmp = Result.r[2][1];
    Result.r[2][1] = Result.r[1][2];
    Result.r[1][2] = tmp;

}

bool TraceContainer::ReadTraceSegNodeXMLFile() {
	NodeContainer.reserve(1000);
	TiXmlDocument doc(xmlfname);
	if (!doc.LoadFile()) {
		std::cout <<"TraceContainer could not open xml file!!"<<std::endl;
		return false;
	}

	//scan each Superellipse
	TiXmlNode* xmlnode; 
	TraceIDVecType TrIDList; //temp list; made unique at end
	for ( xmlnode = doc.FirstChild(); xmlnode != 0; xmlnode = xmlnode->NextSibling()) 	{

		//verify if the xmlnode is a type element
		if (xmlnode->Type()!=TiXmlNode::ELEMENT)	{
			continue;
		}

		//verify if the xmlnode is a superellipse, if not 
		if (strcmp(xmlnode->Value(),"Superellipse"))	{
			continue;
		}

		TraceSegNode *n = new TraceSegNode();
		TiXmlAttribute* pAttrib = xmlnode->ToElement()->FirstAttribute();
		while (pAttrib)	{
			if (!strcmp(pAttrib->Name(),"ID"))	{
				int temp = -1;
				if (pAttrib->QueryIntValue(&temp)==TIXML_SUCCESS)
					n->ID = temp;
			}
			else if (!strcmp(pAttrib->Name(),"TraceID"))	{
				int temp = -1;
				if (pAttrib->QueryIntValue(&temp)==TIXML_SUCCESS) {
					TrIDList.push_back(temp);
					n->TraceID = temp;
				}
			}
			
			else if (!strcmp(pAttrib->Name(),"x"))	{
				double temp = -1.0;
				if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)
					n->mu[0] = temp;
			}

			else if (!strcmp(pAttrib->Name(),"y"))	{
				double temp = -1.0;
				if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)
					n->mu[1] = temp;
			}
			
			else if (!strcmp(pAttrib->Name(),"z"))	{
				double temp = -1.0;
				if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)
					n->mu[2] = temp;
			}
			
			else if (!strcmp(pAttrib->Name(),"a1"))	{
				double temp = -1.0;
				if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)
					n->a1 = temp;
			}
			else if (!strcmp(pAttrib->Name(),"a2"))	{
				double temp = -1.0;
				if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)
					n->a2 = temp;
			}
			else if (!strcmp(pAttrib->Name(),"a3"))	{
				double temp = -1.0;
				if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)
					n->a3 = temp;
			}
			else if (!strcmp(pAttrib->Name(),"Foregd"))	{
				double temp = -1.0;
				if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)
					n->f = temp;
			}
 			else if (!strcmp(pAttrib->Name(),"Backgd"))	{
				double temp = -1.0;
				if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)
					n->b = temp;
			}
			else if (!strcmp(pAttrib->Name(),"q1"))	{
				double temp = -1.0;
				if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)
					n->q1[0] = temp;
			}
			else if (!strcmp(pAttrib->Name(),"q2"))	{
				double temp = -1.0;
				if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)
					n->q1[1] = temp;
			}
			else if (!strcmp(pAttrib->Name(),"q3"))	{
				double temp = -1.0;
				if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)
					n->q1[2] = temp;
			}
			else if (!strcmp(pAttrib->Name(),"q4"))	{
				double temp = -1.0;
				if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)
					n->q1[3] = temp;
			}
			else if (!strcmp(pAttrib->Name(),"e1"))	{
				double temp = -1.0;
				if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)
					n->e1 = temp;
			}
			else if (!strcmp(pAttrib->Name(),"e2"))	{
				double temp = -1.0;
				if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)
					n->e2 = temp;
			}
			n->type = 0;
			pAttrib=pAttrib->Next();
		}

		TiXmlNode* nbr; 
		for ( nbr = xmlnode->FirstChild(); nbr != 0; nbr = nbr->NextSibling())		{
			TiXmlAttribute* nAttr = nbr->ToElement()->FirstAttribute();
			if (!strcmp(nAttr->Name(),"ID"))	{
				int temp = -1;
				if (nAttr->QueryIntValue(&temp)==TIXML_SUCCESS)
					n->NbrID.push_back(temp);
			}
		}
		// DendMuMap: indexed by dendrite(TraceID) and in each 
		// map element there is a PointSetPointer of mu coords
		PointSetPointerMap::iterator   dmuiter;
		dmuiter = DendMuMap.find(n->TraceID);
		if (dmuiter == DendMuMap.end()) {
			currptset = PointSetType::New();
			DendMuMap[n->TraceID]= currptset;
			pcont = PointSetType::PointsContainer::New();
			pdata = PointSetType::PointDataContainer::New();
			currptset->SetPointData(pdata);
			currptset->SetPoints(pcont);
		}
		else {
			currptset = (*dmuiter).second;
			pcont     = currptset->GetPoints();
			pdata     = currptset->GetPointData();
		}
		for (int ptidx=0; ptidx<spr_SPIMGDIMENSION; ptidx++)
			currpt[ptidx]=n->mu[ptidx];
		pcont->push_back(currpt);
		pdata->push_back(n->ID);

		//store in container
		NodeContainer.push_back(n);
	}
	std::sort(TrIDList.begin(), TrIDList.end());
	std::unique_copy( TrIDList.begin(), TrIDList.end(), std::back_inserter( TraceIDList ) );

	return true;

}
void TraceContainer::PrintSelf() {
	if (spr_DEBUGGING < xml_DBG)
		return;
	std::cout << "Nodes found: " << NodeContainer.size() << std::endl;
	std::cout << "Dendrite-MU Map: "<< DendMuMap.size() << std::endl;
	std::cout << "Trace Lines :" << TraceIDList.size()   << "[" << std::endl;
	for (unsigned int i=0; i<TraceIDList.size(); i++)
		std::cout << TraceIDList[i] << ", ";
	std::cout << ']' << std::endl;
}

void TraceContainer::CheckConsistency(SpineImageType::Pointer image) {
	TraceSegNodeVecType::iterator	niter;
	SpineImageType::IndexType		idx;
	ImagePixelType                  muI;
	SpineImageType::SizeType        maxmu = {{0,0,0}};
	int                             inconsistF = 0, inconsistSize=0;
	ImagePixelType					minf,maxf,minb,maxb;
	ImagePixelType					F, B, muImin,muImax;
	minf=minb=muImin = spr_PXLMAXVAL;
	maxf=maxb=muImax = 0;

	xmlconsistent = false;
	if (!xmlvalid)
		return;
	bool       beforestats = true;
	while (!xmlconsistent) {
		for (niter = NodeContainer.begin(); niter != NodeContainer.end(); ++niter) {
			for (int i=0; i<spr_SPIMGDIMENSION; i++) {
				idx[i] = (*niter)->mu[i];
				if (unsigned(idx[i]) > maxmu[i])
					maxmu[i] = idx[i];
			}
			F = (*niter)->f;
			muI = image->GetPixel(idx);
			if (minf > F)
				minf = F;
			if (muImin > muI)
				muImin = muI;
			if (muImax < muI)
				muImax = muI;
			if (maxf < F)
				maxf = F;
			B = (*niter)->b;
			if (minb > B)
				minb = B;
			if (maxb < B)
				maxb = B;
			if(abs(F- muI) - abs (B - muI) > spr_CONSISTENCYTOL)
				inconsistF++;
		}
		SpineImageType::SizeType imgsize = image->GetBufferedRegion().GetSize();
		for ( int i=0; i<spr_SPIMGDIMENSION; i++) {
			if (imgsize[i] < maxmu[i]) {
				std::cout<< "Image size inconsistency: image size[" <<i<< "]="<<imgsize[i]<<". Maximum XML ["<<i<<"]="<<maxmu[i]<<std::endl;
				inconsistSize = 1;
				image->Print(std::cout);
				return;
			}
		}
		if (inconsistF) {
			if (beforestats) {
				std::cout << "Image inconsistent with XML: Foreground discrepancy:" << 100*inconsistF/NodeContainer.size() << "%" << std::endl;
				std::cout << " XML intensities: F=["<<minf<<","<<maxf<<"]; B=["<<minb<<","<<maxb<<"]"<< std::endl;
				std::cout << " Image MU intensities: =["<<muImin<<","<<muImax<<"]"<< std::endl;
				VERBOSE("Mapping image intensities to [0,255]");
				imdbg->ImageStatistics(image);
				inconsistF = 0;
				beforestats = false;
			}
			else
				return;
		}
		else
			xmlconsistent = true;
	}
}
