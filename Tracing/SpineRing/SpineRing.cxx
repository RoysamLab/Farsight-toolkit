#include "SpineRing.h"
#include "ImageProc.h"
#include <itkRegionOfInterestImageFilter.h>
#include <algorithm>
#include <limits.h>

TraceSegNode* GlobalDetector::getSegment(long i) {
	TraceSegNodeVecType::iterator	niter;
	for (niter = NodeList->begin(); niter != NodeList->end(); ++niter) {
		if((*niter)->ID == i)      {
			return(*niter);
		}
	}
	return NULL;
}


GlobalDetector::SpineRing::SpineRing(GlobalDetector *p, spr_SpineRingType srtype) {
	Parent    = p;
	RingType  = srtype;

}

void GlobalDetector::SpineRing::PrintSelf() {
	// writing out foreground pixels inside ring
	Parent->imdbg->WriteTempObj2D(seg->TraceID, seg->ID, &inring_f_idx);
	// writing out entire ring struct
	DEBUGSTMT(Parent->imdbg->WriteTempObj2D(seg->TraceID+100, seg->ID, &inring_all_idx));
}

void GlobalDetector::SpineRing::RestartRing(TraceSegNode* trbit) {
	seg       = trbit;
	thickness = 0;
	validring = false;		
	inring_f_vals.clear();
	DEBUGSTMT(inring_all_idx.clear());
	inring_f_idx.clear();
	short nbrs = seg->NbrID.size();
	for (int i = 0; i < nbrs; i++) {
		TraceSegNode* s1=Parent->getSegment(seg->NbrID[i]);
		if (!s1)
		{
			VERBOSE("XML file appears inconsistent! Neighbor missing!");
			exit(-1);
		}
		if ((s1->TraceID) != (seg->TraceID)){
			continue;
			// Hussein FIXME: is there anything special to be done around 
			// the branchpoint in terms of distance to previous 
			// for setting cmind ...?
		}
		if (s1->type)
			continue;
		thickness += EDIST(seg->mu, s1->mu);
		s1->type = 1;
	}
	if (thickness==0)
		thickness = seg->a3*2.0;
	else
		if (nbrs==1)
			thickness+= seg->a3;
		else
			thickness+= seg->a3/5.0;
	if (RingType==RingTypeSE)
    {
		if ( (validring = RingSE()) )
      {
			DEBUGSTMT(PrintSelf());
			CCimUpdate();
		  }
	  }
	else {
		////////////////// Possible Future Work: //////////////////
		// in case we needed to create the cylindrical detector:
		// Use these for a cylindrical ring
		// amax  = MAX(seg->a1, seg->a2);
		//radi  = MIN(amax*radiscale, amax+4);
		//rado  = MAX(radoscale*amax, amax+10);
		//dlen  = MAX(amax, thickness); 
		//////////////// end ////////////////////////////
		std::cout<<"Warning: cylindrical torus not implemented yet! "<<std::endl;
		exit(1);
		//thering=GetRingCyl(trbit);
	}
}


void GlobalDetector::SpineRing::CCimUpdate() {
	// Add new forground pixels to conn Comp image
	std::vector<ImagePixelType>::iterator viter;
	IndexVecType::iterator                xiter;
	for (viter  = inring_f_vals.begin(),
		xiter  = inring_f_idx.begin();
		viter != inring_f_vals.end();
	viter++, xiter++) {		
		Parent->CCim->SetPixel(*xiter, *viter);
	}
};


GlobalDetector::SpineRing::~SpineRing() {
	inring_f_vals.clear();
	inring_f_idx.clear();
}



void TranslatePt(double* p1, double* Dir, double a, double *p2) {
	// calculates a translation in the direction of Dir:
	// p2 = p1 + a*Dir
	// same as untranslate.m if a<0
	for (int ii=0; ii< spr_SPIMGDIMENSION; ii++) 
		p2[ii] = p1[ii] + (a*Dir[ii]);
}


GlobalDetector::GlobalDetector(SpineImageType::Pointer im, TraceContainer *TCP, ImageDebugger *imdebugger) {
	//GlobalDetector::GlobalDetector(SpineImageType::Pointer im, NodeContainerType nc) {
	NodeList     = &(TCP->NodeContainer);
	TraceIDList  = &(TCP->TraceIDList);
	DendMuMap    = &(TCP->DendMuMap);
	image        = im;
	CCim         = BinaryImageType::New();
	maxmapsz     = 0;
	spcandtotal  = 0;
	//need to make this into an argument
	up_sampling  = spr_UP_SAMPLING;
	imdbg        = imdebugger;
	imgsize		 = image->GetBufferedRegion().GetSize();
	BinaryImageType::IndexType start;
	for (int i=0; i<spr_SPIMGDIMENSION; i++) 
		start[i]=0;
	//  BinaryImageType::SizeType  size;
	BinaryImageType::RegionType region;
	region.SetSize( im->GetLargestPossibleRegion().GetSize() );
	region.SetIndex( start );
	CCim->SetRegions( region );
	CCim->Allocate();

	//// init node->type =0;
	//for(iter = NodeList.begin(); iter != NodeList.end(); ++iter) 
	//	(*iter)->type=0;
}
void GlobalDetector::ResetCCim() {
	//Hussein FIXME: is this the best way to reset the binary image? 
	CCim->FillBuffer(0);
	CCim->Update();

	bbmax   = image->GetBufferedRegion().GetIndex();
	bbmin   = bbmax + imgsize;

}

void GlobalDetector::Run() {
	TraceSegNodeVecType::iterator	niter;
	TraceIDVecType::iterator		titer;
	SpineRing curring(this, RingTypeSE);
	for (titer = TraceIDList->begin(); titer!=TraceIDList->end(); titer++) {
		VERBOSE2("Rings for dendrite ",(*titer));
		ResetCCim();
		for(niter = NodeList->begin(); niter != NodeList->end();  niter++) {
			if ((*niter)->TraceID!=(*titer))
				continue;
			VERBOSE2NL((*niter)->ID, " ");std::flush(std::cout);
			curring.RestartRing(*niter);
			if (curring.validring)
				curring.CCimUpdate();
			//Q.NewRing(&curring);
			//else {
			//Q.Flush();
		}
		VERBOSE("");
		VERBOSE2("Conn Comp Filter on dendrite ",(*titer));
		DoCC(*titer);
	}
}

void GlobalDetector::DoCC(unsigned short trID) {
	CCFilterType::Pointer   CCfilter = CCFilterType::New();
	RelabelType::Pointer    relabel  = RelabelType::New();

	SpineImageType::RegionType trRegion;
	trRegion.SetIndex(bbmin);
	SpineImageType::SizeType regsize;
	for (int i=0; i<spr_SPIMGDIMENSION; i++)
		regsize[i] = abs(bbmax[i]-bbmin[i]);
	trRegion.SetSize(regsize);

	typedef itk::RegionOfInterestImageFilter<BinaryImageType, BinaryImageType> ROIFilterType;
	ROIFilterType::Pointer ROIfilter = ROIFilterType::New();
	ROIfilter->SetRegionOfInterest(trRegion);
	ROIfilter->SetInput(CCim);

	//	CCim->SetRequestedRegion(RingRegion);
	VERBOSE2("....Region Size=", regsize);
	try
	{
		CCfilter->SetInput (ROIfilter->GetOutput());
		CCfilter->Update();
		DEBUGSTMT(VERBOSE("\n Printing CCfilter"));
		DEBUGSTMT(CCfilter->Print(std::cout));
		//VERBOSE(CCfilter->GetOutput()->GetRequestedRegion());
		//VERBOSE(CCfilter->GetOutput()->GetLargestPossibleRegion());
		relabel->SetInput( CCfilter->GetOutput() );
		relabel->SetMinimumObjectSize( spr_MINSPINESIZE );

		relabel->Update();
		//VERBOSE2("Relabel Req Region:", relabel->GetOutput()->GetRequestedRegion());
		//VERBOSE2("Relabel LP Region:", relabel->GetOutput()->GetLargestPossibleRegion());
		DEBUGSTMT(VERBOSE("\n Printing relabel filter"));
		DEBUGSTMT(relabel->Print(std::cout));
		VERBOSE2("....Relabel: Num of Obj=", relabel->GetNumberOfObjects());
		imdbg->MaxIPGenerator(CCim);
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Relabel: exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}
	SpCandMapType *SpCandMap;
	std::map<unsigned short, SpCandMapType*>::iterator   dciter;
	SpCandMapType::iterator         sciter;
	dciter = DendSpCandMap.find(trID);
	if (dciter == DendSpCandMap.end()) {
		VERBOSE2("	new dend-sp-cand-map for traceID", trID);
		SpCandMap = new SpCandMapType();
		DendSpCandMap[trID]= SpCandMap;
	}
	else
		SpCandMap = (*dciter).second;
	//need to access the CC objects and create spinecandidates for each 	unsigned short numObjects = relabel->GetNumberOfObjects();
	itk::ImageRegionIterator<CCOutputImageType>
		it(relabel->GetOutput(), relabel->GetOutput()->GetRequestedRegion());

	///// GLOBAL: make sure vectors are initialized
	////          make sure use Update and GetBufferedRegion or GetMaxPossible
	SpCandidate *spcand;
	CCOutputImageType::PixelType pxlabel;
	int regfpxls =0 , regpxls= 0;
	while( !it.IsAtEnd() )
	{
		pxlabel  = it.Get();
		if (pxlabel != 0)
		{ //push it to spine candidate number obj
			sciter = SpCandMap->find(pxlabel);
			if (sciter == SpCandMap->end()) {
				DEBUGSTMT(VERBOSE("		New spine Cand"));
				spcand = new SpCandidate(this, trID, pxlabel);
				SpCandMap->insert(std::pair<unsigned short, 
					SpCandidate*>
					(pxlabel, spcand));
			}
			else
				spcand = (*sciter).second;
			spcand->AddPoint(it.GetIndex(), bbmin);
			regfpxls++;
		}
		++it;
		regpxls++;
	}
	DEBUGSTMT(VERBOSE2("	regional foreground pxls =", regfpxls));
	DEBUGSTMT(VERBOSE2("	All regional pxls        =", regpxls));
	PrintDendSpCandMap(trID);
}

void GlobalDetector::PrintDendSpCandMap(unsigned short trID) {
	if (spr_DEBUGGING < cand_DBG)
		return;
	//	SpCandidate *spcand;
	SpCandMapType *SpCandMap;
	std::map<unsigned short, SpCandMapType*>::iterator   dciter;
	SpCandMapType::iterator         sciter;
	dciter = DendSpCandMap.find(trID);
	if (dciter == DendSpCandMap.end()) {
		VERBOSE2("	ERROR: Printing EMPTY DendSpCandMap for traceID ", trID);
		return;
	}
	else
		SpCandMap = (*dciter).second;
	std::cout << "number of spines in traceID " << trID << "is " << SpCandMap->size();
	for (sciter = SpCandMap->begin(); sciter !=SpCandMap->end(); sciter++) {
		(*(*sciter).second).PrintSelf();
	}
}

void SpCandidate::PrintSelf() {
	std::cout<< "spine cand " << candID << " of trace " << traceID << "size=" << pxls->GetNumberOfPoints() <<std::endl;
}

SpCandidate::SpCandidate(GlobalDetector* gd, int trID, unsigned short cid) : valid(true), 
ParentDet(gd), candID(cid), traceID(trID) {
	pxls  = PointSetType::New();
	pcont = PointSetType::PointsContainer::New();
	pdata = PointSetType::PointDataContainer::New();
	pxls->SetPoints(pcont);
	pxls->SetPointData(pdata);
	NbrMus       = PointSetType::New();
	NbrMusPtCont = PointSetType::PointsContainer::New();
	NbrMus->SetPoints(NbrMusPtCont);
	PSDists = 0;
	/*pditer = pdata->Begin();
	pciter = pcont->Begin();*/
	ptID=0;
	sigmaI = 0;
	munode = 0;
	imgsize  = ParentDet->GetImageSize();
	scbbmax  = ParentDet->image->GetBufferedRegion().GetIndex();
	scbbmin   = scbbmax + imgsize;
	SpCandRegion.SetSize(0,0);
	spinesize = 0;
	validPath = false;
}

SpCandidate::~SpCandidate()
{
	delete PSDists;
}

void SpCandidate::AddPoint(SpineImageIndexType pxl, SpineImageIndexType bbmin ) {
	PointType p;
	SpineImageIndexType ptidx;
	for (int j=0; j<3; j++)
	{
		p[j]     = pxl[j] + bbmin[j];
		ptidx[j] = pxl[j] + bbmin[j];
		// update spine pixels bounding box
		if (scbbmin[j] > ptidx[j])
			scbbmin[j] = ptidx[j];
		if (scbbmax[j] < ptidx[j])
			scbbmax[j] = ptidx[j];
	}
	pcont->InsertElement(ptID, p);
	//		.push_back(pxl);
	ImagePixelType ptval = ParentDet->image->GetPixel(ptidx);
	pdata->InsertElement(ptID++, ptval);
	spinesize++;
	//pxlvals.push_back(ParentDet->image->GetPixel(pxl));
	//sigmaI += ptval;

}


void SpCandidate::GetClosestMUs() {

	const PointSetType::Pointer mus = ParentDet->DendMuMap->find(traceID)->second;
	PSDists = new PointSetDistS(pxls, mus);

	ImagePixelType nodeID = 0;
	//TraceSegNodeVecType  *NodeList = ParentDet->NodeList;
	TraceSegNode		 *seg;

	PSDists->set2->GetPointData(PSDists->minidx2, &nodeID);
	munode = ParentDet->getSegment((int)nodeID);

	std::vector<int>::iterator psiter;
	for (psiter = PSDists->minidx21a.begin(); psiter!=PSDists->minidx21a.end();
       psiter++)
    {
		PSDists->set2->GetPointData((*psiter), &nodeID);
		seg = ParentDet->getSegment((int)nodeID);
		GetNbrMUs(seg, MUOFFSET);		
	  }
	//	std::sort(mucands.begin(), mucands.end());
	//	std::unique_copy( mucands.begin(), mucands.end(), 
	//						std::back_inserter( CloseNodeidx ) );
	//sproot = PSDists->pt1m;
	//mu     = PSDists->pt2m;
	//X      = PSDists->pt1x;
	//X2     = PSDists->pt1x2;
	//mui1   = PSDists->minidx2;
	//mucandidx = CloseNodeidx;
}

void SpCandidate::GetNbrMUs(TraceSegNode *seg, int count)
{
	PointType currpt;
	for (int ptidx=0; ptidx<spr_SPIMGDIMENSION; ptidx++)
		currpt[ptidx]=seg->mu[ptidx];
	NbrMusPtCont->push_back(currpt);
	NbrMuIDVec.push_back(seg->ID);
	short nbrs = seg->NbrID.size();
	if (count)
    {
		for (int i = 0; i < nbrs; i++)
      {
			TraceSegNode* s1=ParentDet->getSegment(seg->NbrID[i]);
			if ((s1->TraceID) != (seg->TraceID))
        {
				continue;
        }
			bool found=false;
			std::vector<int>::iterator iter;
			for (iter = NbrMuIDVec.begin(); iter!=NbrMuIDVec.end(); iter++)
        {
        if (*iter == s1->ID) 
          {
          found=true;
          break;
          }
        }
      if (!found)
        {
        GetNbrMUs(s1, count-1);
        }
		  }
	  }
}

bool SpCandidate::Validate() {
	PointSetType::PointDataContainer::Iterator pdataiter;
	PointSetType::PointsContainerIterator      pciter, pciter2;
	if (!munode)
		GetClosestMUs();
	ImagePixelType	f = munode->f;
	ImagePixelType	b = munode->b;
	ImagePixelType	pxlval, minpxlval = spr_PXLMAXVAL;
	int				sumAD=0;

	for (pdataiter=pdata->Begin(), pciter=pcont->Begin(); pdataiter!=pdata->End();	pdataiter++, pciter++)
	{
		pxlval = pdataiter.Value();
		sumAD  += -abs(pxlval-f) + abs(pxlval-b);
		if (pxlval < minpxlval)
		{
			minpxlval = pxlval;
			pciter2   = pciter;
		}
	}
	brtpt = pciter2->Value();
	sumAD /= pxls->GetNumberOfPoints();
	VERBOSE2NL( candID , ": ");
	PathPtsOfInt();
	return (valid = (sumAD >= spr_LIKELIHOODTHRESH));
}


void GlobalDetector::ValidateSpCands() {  
	std::map<unsigned short, SpCandMapType*>::iterator   dciter;
	SpCandMapType::iterator         sciter;
	for (dciter = DendSpCandMap.begin();
		dciter != DendSpCandMap.end(); dciter++) 
	{
		int currmapsz = dciter->second->size();
			for (sciter=dciter->second->begin();
				sciter!=dciter->second->end(); sciter++)
			{
				
				VERBOSE2NL(currmapsz, "-");
				sciter->second->Validate();
			}
		}
}



PointSetDistS::PointSetDistS(PointType p1, PointSetType::Pointer  s2) 
{
	set2 = s2;
	set1 = PointSetType::New();
	set1->SetPoint(0, p1);
	Compute();
}
PointSetDistS::PointSetDistS(PointSetType::Pointer  s1, PointSetType::Pointer  s2)
{
	set1 = s1;
	set2 = s2;
	Compute();
}

void PointSetDistS::Print()
{
		std::cout<<"dmin="<<dmin<<"; dmax= "<<dmax<<"; pt1m="<<pt1m<<"; pt2m="<<pt2m<<std::endl;
		//std::cout<<"mind21="<<mind21<<"; pt1x= "<<pt1x<<"; pt1x2="<< pt1x2       <<std::endl; 
		std::cout<< "pt1x= "<<pt1x<<"; pt1x2="<< pt1x2       <<std::endl; 
		std::cout<<"minidx21a=";
		for (unsigned int i=0; i<minidx21a.size(); i++)
			std::cout<<minidx21a[i]<<" ";
		std::cout<<std::endl;
		std::cout<<"mind21a=";
		for (unsigned int i=0; i<mind21a.size(); i++)
			std::cout<<mind21a[i]<<" ";
		std::cout<<std::endl;
}

void PointSetDistS::Compute(){
	//% takes two sets of 3D points, finds closest 2 points 
	//% then takes the pt from set2 and finds the farthest in set1
	//% set1 is 3xm  set2 is 3xn
	//% usepdist=0 or 1 
	//% dmin = dist between closest 2 pts
	//% pt1m & pt2m are closest 2 pts in set1 and set2 respect.
	//% dmax = dist from pt2m to pt1x      (dist dend center to tip)
	//% pt1x= farthest pt in set1 to set2 (tip)
	//% pt1x2= point in set1 at max integral distance (sum of distances)to set2
	//% minidx2= index of pt2m in set2
	//% minidx2a =vector of indeces in set2 for all min pts in 2
	//% mind21a  = vector of minimum distances from 2 to 1
	// FD = Discrete Frechet dist = min (max distances array from set1 to set2)
	// HD = Discrete Housdorf dist = max(min ................................) 
	//int set1size = set1->GetNumberOfPoints();
	//int set2size = set2->GetNumberOfPoints();
	int mini = 10000000.0;
  int minj = 10000000.0;
  int maxi = -10000000.0;
  int maxj = 10000000.0;
  int minji = 10000000.0;
  int summaxidx = -10000000.0;
	double d, currmax=0, currpt1mind, currmin = 10000000.0;//some big number
	std::vector<int> minidx21;
	std::vector<double> maxd21a;
	double sumd=0, currsummax=0, currpt1maxd;
	PointType p1, p2;
	c1 = set1->GetPoints();
	c2 = set2->GetPoints();
	int i, j;
	for (i1=c1->Begin(), i=0; i1!=c1->End(); i1++, i++) 
	{
		p1 = i1.Value();
		sumd=0;
		currpt1mind=10000000.0;
		currpt1maxd=0.0;
		for (i2=c2->Begin(), j=0; i2!=c2->End(); i2++, j++) 
		{
			p2 = i2.Value();
			d = ED( p1, p2);
			sumd += d;
			if (d<currmin) 
			{
				currmin=d;
				mini=i;minj=j;
			}
			else
				if (d>currmax) 
				{
					currmax=d;
					maxi=i;maxj=j;
				}
			if (d<currpt1mind) 
			{
				currpt1mind=d;
				minji=j;
			}
			else
				if ( d>currpt1maxd)
				{
					currpt1maxd = d;
				}

		}
		minidx21.push_back(minji);
		mind21a.push_back(currpt1mind);
		maxd21a.push_back(currpt1maxd);
		if (sumd>currsummax) {
			currsummax=sumd;
			summaxidx = i;
		}
	}
	dmin  = currmin;
	set1->GetPoint(mini, &pt1m);
	set2->GetPoint(minj, &pt2m);
	dmax  = currmax;
	set1->GetPoint(maxi, &pt1x);
	set1->GetPoint(summaxidx, &pt1x2);
	FD = *(std::min_element(maxd21a.begin(),maxd21a.end()));
	maxd21a.clear();
	HD = *(std::max_element (mind21a.begin(), mind21a.end()));
	minidx2 = minj;
	std::sort(minidx21.begin(), minidx21.end());
	std::unique_copy( minidx21.begin(), minidx21.end(), std::back_inserter( minidx21a ) );
}


int SpCandidate::WritePxlsToRGBImage(RGBImageType::Pointer rgbim2D, RGBPixelType color) 
{
	return ParentDet->imdbg->PlotContOnRGBImage(rgbim2D, color, pcont);
}

void SpCandidate::WritePoIToRGBImage(RGBImageType::Pointer rgbim2D)
{
	RGBPixelType color;
	color.SetRed(255);
	color.SetGreen(0);
	color.SetBlue(0);
	//return ParentDet->imdbg->PlotOnRGBImage(rgbim2D, color, ParentDet->imdbg->MakeBox(brtpt));
	ParentDet->imdbg->PlotPtAsBoxOnRGBImage(rgbim2D, color, brtpt, 1, 1);
	color.SetRed(0);
	color.SetGreen(255);
	color.SetBlue(0);
	ParentDet->imdbg->PlotContAsBoxOnRGBImage(rgbim2D, color, NbrMusPtCont, 1, 1);
}

void GlobalDetector::GetSpCandStats()
{
	//get a count of total spine candidates and 
	// maximum spines per dendrite (AKA map)
	if (spcandtotal !=0 && maxmapsz!=0)
		return;
	unsigned int currmax  = 0;
	unsigned int candtot  = 0;
	std::map<unsigned short, SpCandMapType*>::iterator   dciter;
	SpCandMapType *SpCandMap;
	for (dciter = DendSpCandMap.begin(); dciter != DendSpCandMap.end(); dciter++) 
	{
		SpCandMap = dciter->second;
		if (currmax < SpCandMap->size())
			currmax = SpCandMap->size();
		candtot += SpCandMap->size();
	}
	maxmapsz = currmax;
	spcandtotal = candtot;
}

RGBImageType::Pointer   GlobalDetector::GetRGBPOImage() 
{
	return imdbg->rgbim2D; 
}

int GlobalDetector::WriteDendSpCandMap2RGB(RGBImageType::Pointer rgbim2D) 
{
	std::map<unsigned short, SpCandMapType*>::iterator   dciter;
	SpCandMapType::iterator         sciter;
	int validspc=0;
	for (dciter = DendSpCandMap.begin(); dciter != DendSpCandMap.end(); dciter++) 
	{
		for (sciter=dciter->second->begin();
			sciter!=dciter->second->end(); sciter++)
		{
			validspc += sciter->second->WritePxlsToRGBImage(rgbim2D, GetColorMap(sciter->first));
			
            DEBUGSTMT(std::stringstream tmpstr);
            DEBUGSTMT(tmpstr << imdbg->basefilename << "_spcpxls_" << sciter->second->GetcandID() << "_" << sciter->second->GettraceID() << ".png" );
			DEBUGSTMT(writeImageToFile(rgbim2D, tmpstr.str() ));
		}
	}
	return validspc;
}



////////////////////////////////
////////////////////
/////////////
////
//
//////void PointSetDistS::Compute(){
////////% takes two sets of 3D points, finds closest 2 points 
////////% then takes the pt from set2 and finds the farthest in set1
////////% set1 is 3xm  set2 is 3xn
////////% usepdist=0 or 1 
////////% dmin = dist between closest 2 pts
////////% pt1m & pt2m are closest 2 pts in set1 and set2 respect.
////////% dmax = dist from pt2m to pt1x      (dist dend center to tip)
////////% pt1x= farthest pt in set1 to set2 (tip)
////////% pt1x2= point in set1 at max integral distance (sum of distances)to set2
////////% minidx2= index of pt2m in set2
////////% minidx2a=vector of indeces in set2 for all min pts in 2
////////% mind21 = vector of minimum distances from 2 to 1
////////function [dmin pt1m pt2m dmax pt1x pt1x2 minidx2 minidx2a mind21] =...
//////  //  testClosest2Ptsets(set1, set2)
//////	int set1size = set1->size();
//////	int set2size = set2->size();
//////	int mini, minj,maxi,maxj, minji, summaxidx;
//////	double d, currmax=0, currmini, currmin = 10000000.0;//some big number
//////	std::vector<int> minidx21;
//////	double sumd=0, currsummax=0;
//////	//vnl_matrix_fixed<double, set1size, set2size> dist12;
//////	//vnl_vector_fixed<int, set2size>di;
//////	vnl_vector_fixed <double, 3> p1, p2;
//////	for (int i=0; i<set1size; i++) {
//////		for (int jj=0; jj<3; jj++)
//////				p1[jj] = set1[i][jj];
//////		sumd=0;
//////		currmini=10000000.0;
//////		for (int j=0; j<set2size; j++) {
//////			for (int jj=0; jj<3; jj++)
//////				p2[jj] = set2[j][jj];
//////			d = EDIST(p1,p2);
//////			sumd += d;
//////			if (d<currmin) {
//////				currmin=d;
//////				mini=i;minj=j;
//////				//dist12(i,j) = EDIST(p1,p2);
//////			}
//////			else
//////			if (d>currmax) {
//////				currmax=d;
//////				maxi=i;maxj=j;
//////			}
//////			if (d<currmini) {
//////				currmini=d;
//////				minji=j;
//////			}
//////		}
//////		minidx21.push_back(minji);
//////		mind21a.push_back(currmini);
//////		if (sumd>currsummax) {
//////			currsummax=sumd;
//////			summaxidx = i;
//////		}
//////	}
//////	dmin  = currmin;
//////	pt1m  = set1[mini];
//////	pt2m  = set2[minj];
//////    dmax  = currmax;
//////	//for (int ii=0; ii<3; ii++) {
//////	pt1x  = set1[maxi];
//////	pt1x2 = set1[summaxidx];
//////	//}
//////	minidx2 = minj;
//////	std::sort(minidx21.begin(), minidx21.end());
//////	std::unique_copy( minidx21.begin(), minidx21.end(), std::back_inserter( minidx21a ) );
//////}



//SpineDetector::Run() {
//	// for vv all vessels
//	// for jj all segments
//	num=CurrVess.num;
//	mu = CurrVess.Segments[jj].mu;       
//    q1 = seg.R1;
//	q2 = seg.R2;
//    q3 = seg.R3;
//	a3=seg.a3;
//	nexti = min(num, jj+1);
//	previ = max(1,jj-1);
//    nextseg=CurrVess.Segments(nexti);
//    prevseg=CurrVess.Segments(previ); 
//	Vect3 mu0, nextmu0, prevmu0;
//	TranlatePt(mu, q3, a3, 1, mu0);
//	TranlatePt(nextseg.mu, nextseg.q3, nextseg.a3, nextmu0);
//	TranlatePt(prevseg.mu, prevseg.q3, prevseg.a3, prevmu0);
//
//    if (jj==1)       // cmind = norm(mu0 - nextmu0)+a3;
//		cmind = EDIST(mu0, nextmu0)+a3;
//    else
//		if (jj==num) // cmind = norm(mu0 - prevmu0)+a3;
//			cmind = EDIST(mu0, prevmu0)+a3;
//		else         //  cmind = norm(mu0- nextmu0)+a3/5;
//			cmind = EDIST(mu0, nextmu0)+a3/5;            
//    //R  = q1 q2 q3;
//    segcylmatchcnt  = 0;
//    temp_spines     = [];            
//    vr_detect       = 0;
//    current_minfit  = minfit;
//    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//    %% RING CALCULATION     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//	if (USEDONUT) {
//        //[c3 bio] = donutb(radi, rado, dlen, cylres, 0);
//        //S1  = R*c3; 
//        //bio=R*bio;
//        //ringinfo.ri=radi; ringinfo.ro=rado;ringinfo.dlen=dlen;
//	}
//    else
//        SpineRing(seg, cmind, RingTypeSE);
//           
//
//}

//void DetectorQ::NewRing(SpineRing* curring,
//						IndexVecType mindx,
//						IndexVecType maxdx) 
//{
//	Q.push_back(curring);
//	if (Q.size() > maxq) 
//		Flush();
//}
//
//DetectorQ::DetectorQ(){
//	qminidx = {0,0,0};
//	qmaxidx = {0,0,0};
//}
//
//DetectorQ::~DetectorQ(){
//	Q.clear();
//}

//GlobalDetector::SpineRing::SpineRing(GlobalDetector *p, TraceSegNode* trbit, SpineRingType srtype) {
//	Parent    = p;
//	RingType  = srtype;
//	seg       = trbit;
//	thickness = 0;
//	validring = false;
//	short nbrs = seg->NbrID.size();
//	for (int i = 0; i < nbrs; i++) {
//		TraceSegNode* s1=Parent->getSegment(seg->NbrID[i]);
//		if ((s1->TraceID) != (seg->TraceID)){
//			continue;
//		// Hussein FIXME: is there anything special to be done around 
//		// the branchpoint in terms of distance to previous 
//		// for setting cmind ...?
//        }
//		if (s1->type)
//			continue;
//		thickness += EDIST(seg->mu, s1->mu);
//		s1->type = 1;
//	}
//	if (thickness==0)
//		thickness = seg->a3*2.0;
//	else
//		if (nbrs==1)
//			thickness+= seg->a3;
//		else
//			thickness+= seg->a3/5.0;
//	if (RingType==RingTypeSE) {
//		if (validring = RingSE())
//			CCimUpdate();
//	}
//	else {
//	 //   amax  = MAX(seg->a1, seg->a2);
//		//radi  = MIN(amax*radiscale, amax+4);
//		//rado  = MAX(radoscale*amax, amax+10);
//		//dlen  = MAX(amax, thickness); 
//		std::cout<<"Warning: cylindrical torus not implemented yet! "<<std::endl;
//		exit(1);
//		//thering=GetRingCyl(trbit);
//	}
//}
