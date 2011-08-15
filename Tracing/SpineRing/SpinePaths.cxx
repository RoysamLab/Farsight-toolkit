#include "SpineRing.h"
#include "ImageProc.h"
#include <itkRegionOfInterestImageFilter.h>
//#include "itkRescaleIntensityImageFilter.h"

RGBImageType::Pointer GlobalDetector::GetRGBPathsImage()
{
	if (RGBPathsImage.IsNull())
		RGBPathsImage = imdbg->CreateRGBImFrom2D(imdbg->MIP);
	return RGBPathsImage;
}

void SpCandidate::PathPtsOfInt() 
{
	//%function spgpp = PathPoI(spcell, maxiterinit, MU)
	//% Here we take pxls and extract points of interest such as brightest point
	//% as most likely spine tip, sproot, bmu=mu closest to brtpt, X is farthest
	//% pt from all mus, X2 is integrally farthest (sum of distances to mus)
	//global IMGpe DBGG
	//global CVpe CVIpe
	//% this will get an approximate closest mu and tip
	//%spgpp         = spcell;
	//int muoffset=10; //%before and after current mu to find best path endpt on CV
	//roundpxl      = unique(round(spcell.pxl)','rows')';

	//int maxiter = NUMMAX(spr_MAXPATHITERINIT, 2*pxls->GetNumberOfPoints());
	//[dumd2mu spgpp.sproot spgpp.mu dumu2tip spgpp.X spgpp.X2 mui1 mucandidx] = ...
	//   testClosest2Ptsets(roundpxl, MU);
	//sproot = PSDists.pt1m;
	//mu     = PSDists.pt2m;
	//X      = PSDists.pt1x;
	//X2     = PSDists.pt1x2;
	//mui1   = PSDists.minidx2;
	//mucandidx = CloseNodeidx;
	// get candidate end points
	//mucandidx = max(mucandidx(1)-muoffset,1):min(mucandidx(end)+muoffset,CVpe.num);
	////spgpp.mucand=[CVpe.Segments(mucandidx).mu];
	//% 2nd estimate
	//MUfsbl = spgpp.mucand;
	//[len fpx spgpp.brtpt]= Check_Match4(IMGpe,CVpe, spgpp.mui1, 1, roundpxl, 0, 1);

	if (!valid)
		//warning(['invalid spine found. Size='...
		//  num2str(size(roundpxl,2)) '. Fit=' num2str(len)])
		return;

	PointSetDistS FsblMuDists(brtpt, NbrMus);
	//[dumd2mu dumbrtpt spgpp.bmu] = testClosest2Ptsets(spgpp.brtpt, MUfsbl);
	bmu = FsblMuDists.pt2m;
	if (spr_DEBUGGING >= path_DBG) {
		WritePoIToRGBImage(ParentDet->GetRGBPOImage());
		std::stringstream tmpstr;
		tmpstr << ParentDet->imdbg->basefilename << "_spcand_" << candID << "_" << traceID << ".png" ;
		//writeImageToFile<RGBImageType>(ParentDet->GetRGBImage(), tmpstr.str() );
		writeImageToFile(ParentDet->GetRGBPOImage(), tmpstr.str() );
	}
	RegionToPath();
	//spgpp.Vess = CVIpe;
	///////////////////////////////////////////////////////////////
	////// to implement next
	//////currpath   = get_spine_path35({spgpp}, IMGpe, maxiter);
	//////[vpath  pathL plen gd]  = ValidatePath(currpath, spgpp);
	//////if isempty(vpath)
	//////    warning('empty path')
	//////    spgpp.path = [];
	//////    spgpp.invalid=1;
	//////else
	//////    spgpp.path = vpath;
	//////    if currpath{4,1}
	//////       %this is the final dendrite mu closest to the spine after path
	//////       % extraction to the shortest mu in mucand
	//////        spgpp.mui1 = mucandidx(currpath{4,1});
	//////        spgpp.mu   = spgpp.mucand(:,currpath{4,1});
	//////    end
	//////end
	//////spgpp.pathL = pathL;
	//////spgpp.gd    = gd;
	//////spgpp.plen  = plen;
	////// end to implement next
	//////////////////////////////////////////////////////////////
}
SpineImageType::RegionType SpCandidate::GetSpineRegion() 
{
	if (SpCandRegion.GetSize()[0] == 0)
		ComputeSpineRegion();
	return SpCandRegion;
}

void SpCandidate::ComputeSpineRegion()
{
	//pxls = [spcell{1}.pxl stpt{1} endpt{2}];
	PointSetType::PointsContainerIterator	pciter = NbrMusPtCont->Begin();
	PointType    pt;
	SpineImageIndexType pixelIndex;
	while (pciter != NbrMusPtCont->End()) 
	{
		pt = pciter->Value();
		bool isInside = ParentDet->image->TransformPhysicalPointToIndex( pt, pixelIndex );
		if ( isInside )
		{
			for (int j=0; j<spr_SPIMGDIMENSION; j++) 
			{
				if (scbbmin[j] > pixelIndex[j])
					scbbmin[j] = pixelIndex[j];
				if (scbbmax[j] < pixelIndex[j])
					scbbmax[j] = pixelIndex[j];
			}
		}
		pciter++;
	}

	SpineImageIndexType rmin, rmax;
	rmin[0] = NUMMAX(0,          scbbmin[0]-spr_PATHREGIONOFFSET);
	rmax[0] = NUMMIN(static_cast<int>(imgsize[0]), scbbmax[0]+spr_PATHREGIONOFFSET);

	rmin[1] = NUMMAX(0,          scbbmin[1]-spr_PATHREGIONOFFSET);
	rmax[1] = NUMMIN(static_cast<int>(imgsize[1]), scbbmax[1]+spr_PATHREGIONOFFSET);

	rmin[2] = NUMMAX(0,          scbbmin[2]-spr_PATHREGIONOFFSET);
	rmax[2] = NUMMIN(static_cast<int>(imgsize[2]), scbbmax[2]+spr_PATHREGIONOFFSET);

	SpCandRegion.SetIndex(rmin);
	SpineImageType::SizeType regsize;
	for (int i=0; i<spr_SPIMGDIMENSION; i++)
		regsize[i] = abs(rmax[i]-rmin[i]);
	SpCandRegion.SetSize(regsize);
	return;
}





PointSetContainerType::Pointer SpCandidate::TranslatePts(PointSetContainerType::Pointer cont, short dir)
{
	PointSetContainerType::Pointer regionptcont = PointSetType::PointsContainer::New();
	PointSetType::PointsContainerIterator	pciter = cont->Begin();
	PointType								pt, regionpt;
	//SpineImageIndexType						pixelIndex;
	while (pciter != cont->End()) 
	{
		pt = pciter->Value();
		//bool isInside = ParentDet->image->TransformPhysicalPointToIndex( pt, pixelIndex );
		//if ( isInside )
		//{
		for (int j=0; j<spr_SPIMGDIMENSION; j++) 
		{
			//regionpt[j]=pt[j];
			regionpt[j] = NUMMAX(pt[j] + dir*GetSpineRegion().GetIndex()[j], 0);
		}
		regionptcont->push_back(regionpt);
		//}
		pciter++;
	}
	return regionptcont;
}

PointType SpCandidate::TranslatePts(PointType pt0, short dir)
{
	PointType								regionpt;
	//SpineImageIndexType						pixelIndex;
	//bool isInside = ParentDet->image->TransformPhysicalPointToIndex( pt0, pixelIndex );
	//if ( isInside )
	//{
	for (int j=0; j<spr_SPIMGDIMENSION; j++) 
	{
		regionpt[j] = NUMMAX(pt0[j] + dir*GetSpineRegion().GetIndex()[j], 0);
	}
	//}
	return regionpt;
}

SpeedImageType::Pointer SpCandidate::GetSpeedSubIm()
//SpineImageType::Pointer SpCandidate::GetSpeedFunction(SpineImageType::RegionType spregion)
{
	typedef itk::RegionOfInterestImageFilter<SpineImageType, SpineImageType> ROIFilterType;
	ROIFilterType::Pointer ROIfilter = ROIFilterType::New();
	ROIfilter->SetRegionOfInterest(GetSpineRegion());
	ROIfilter->SetInput(ParentDet->image);
	ROIfilter->Update();
	DEBUGSTMT(VERBOSE2("....SpeedFunc calc: Region Size=", GetSpineRegion().GetSize()));

	SpeedImageType::Pointer SpeedFun = SpeedImageType::New();
	SpeedFun->SetRegions( ROIfilter->GetOutput()->GetRequestedRegion() );
	SpeedFun->CopyInformation( ROIfilter->GetOutput() );
	SpeedImageType::PointType Origin;
	Origin.Fill(0);
	SpeedFun->SetOrigin(Origin);
	SpeedFun->Allocate();
	itk::ImageRegionIterator<SpineImageType>
		it(ROIfilter->GetOutput(), ROIfilter->GetOutput()->GetRequestedRegion());
	itk::ImageRegionIterator<SpeedImageType> spdit(SpeedFun, SpeedFun->GetRequestedRegion());
	SpeedPixelType spdpxl;
	int counter=0;
	for ( it.GoToBegin(), spdit.GoToBegin(); !spdit.IsAtEnd(); ++spdit, ++it)
	{
		spdpxl = (255-it.Get()); //spregion = (spregion-m)/(M-m);  %0 to 1
		//spdpxl = 1.0-spdpxl + 1.0/255.0;         //spregion = 254*(1-spregion) + 1;%reverse
		spdit.Set(spdpxl);
#if 0
		SpeedImageType::IndexType idx;
		idx[0]=6;idx[1]=12;idx[2]=20;
		if (spdit.GetIndex()==idx)
			std::cout<<spdit.GetIndex()<<"==>"<<spdpxl<<std::endl;
		idx[0]=1;idx[1]=29;idx[2]=16;
		if (spdit.GetIndex()==idx)
			std::cout<<spdit.GetIndex()<<"==>"<<spdpxl<<std::endl;
#endif
		//spdit.Set(it.Get());
		counter++;
	}
	return SpeedFun;
}

#if 0
// need to make image go from 0 to 1 with 1 being fastest
// so need to reverse it to become bright fg on dark bg.
SpeedPixelType spmin=USHRT_MAX, spmax=0;

itk::ImageRegionIterator<SpineImageType>
it(ROIfilter->GetOutput(), ROIfilter->GetOutput()->GetRequestedRegion());
ImagePixelType sppxl;

for ( it.GoToBegin(); !it.IsAtEnd(); ++it)
{
	sppxl  = it.Get();
	if (spmin > sppxl)
		spmin = sppxl;
	if (spmax < sppxl)
		spmax = sppxl;
}

SpeedImageType::Pointer SpeedFun = SpeedImageType::New();
SpeedFun->SetRegions( ROIfilter->GetOutput()->GetRequestedRegion() );
SpeedFun->CopyInformation( ROIfilter->GetOutput() );
SpeedImageType::PointType Origin;
Origin.Fill(0);
SpeedFun->SetOrigin(Origin);
SpeedFun->Allocate();
itk::ImageRegionIterator<SpeedImageType> spdit(SpeedFun, SpeedFun->GetRequestedRegion());
SpeedPixelType spdpxl;
int counter=0;
for ( it.GoToBegin(), spdit.GoToBegin(); !spdit.IsAtEnd(); ++spdit, ++it)
{
	spdpxl = (it.Get()-spmin)/(spmax-spmin); //spregion = (spregion-m)/(M-m);  %0 to 1
	spdpxl = 1.0-spdpxl + 1.0/255.0;         //spregion = 254*(1-spregion) + 1;%reverse
	spdit.Set(spdpxl);
	counter++;
}
std::cout << "Spine region: Set "<< counter << "pixels. SpeedFunc Size=" << SpeedFun->GetLargestPossibleRegion().GetSize() << std::endl;

//////////////////////////////////////////////////////////
// DEBUGGING ONLY ////////////////////////////////////////
//////////////////////////////////////////////////////////

typedef itk::RescaleIntensityImageFilter< 
SpeedImageType, 
SpineImageType >   CastFilterType;
typedef  itk::ImageFileWriter<  SpineImageType  > WriterType;
CastFilterType::Pointer caster1 = CastFilterType::New();
WriterType::Pointer writer1 = WriterType::New();
caster1->SetInput( SpeedFun );
writer1->SetInput( caster1->GetOutput() );
writer1->SetFileName("FastMarchingFilterOutput1.tiff");
caster1->SetOutputMinimum(   0 );
caster1->SetOutputMaximum( 255 );
writer1->Update();
//////////////////////////////////////////////////////////
// DEBUGGING END  ////////////////////////////////////////
//////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////
// DEBUGGING ONLY ////////////////////////////////////////
//////////////////////////////////////////////////////////
//SpineImageType::Pointer SpineFun = SpineImageType::New();
//SpineFun->SetRegions( ROIfilter->GetOutput()->GetRequestedRegion() );
//SpineFun->CopyInformation( ROIfilter->GetOutput() );
//SpineFun->Allocate();
//itk::ImageRegionIterator<SpineImageType> spdit2(SpineFun, SpineFun->GetRequestedRegion());
//float spdpxl2;
//for (  it.GoToBegin(), spdit2.GoToBegin(); !spdit2.IsAtEnd(); ++spdit2, ++it)
//{
//	spdpxl2= (it.Get()-spmin)/(spmax-spmin); //spregion = (spregion-m)/(M-m);  %0 to 1
//	spdit2.Set(ImagePixelType(254*(1-spdpxl2) + 1));                //spregion = 254*(1-spregion) + 1;%reverse
//	if (it.Get()<120)
//		std::cout <<"Reg Pxl: Orig, [0,1], streched:"<<it.Get()<<",  "<<spdpxl2<<",  "<<spdit2.Get()<<std::endl;
//}
//SpineImage2DType::Pointer Sp2DReg = MaxIPGenerator(SpineFun);
//std::stringstream tmpstr;
//tmpstr << ParentDet->imdbg->basefilename << "_MaxIPReg_spcand_" << candID << "_" << traceID << ".png" ;
//writeImageToFile(Sp2DReg, tmpstr.str());
//////////////////////////////////////////////////////////
// DEBUGGING END  ////////////////////////////////////////
//////////////////////////////////////////////////////////


//return SpineFun;
return SpeedFun;
}
#endif

//template  int writeImageToFile<SpeedImage2DType>(SpeedImage2DType::Pointer , std::string );
//template
///int writeImageToFile<SpeedImage2DType>(SpeedImage2DType::Pointer , std::string );

//template SpeedImage2DType::Pointer MIPGenerator<SpeedImageType,SpeedImage2DType>(SpeedImageType::Pointer im);



# if 0
/*
function DendSpF = testPathExtract2(DendSp, img, CVi, CV)
% 04/17/09   muoffset changed from 4 to 10 since now we are taking median geodesic distance for endpoints
global IMGpe DBGG DOPATHMERGE
global CVpe CVIpe
CVpe=CV;
IMGpe=img;
CVIpe=CVi;
maxiterinit=50000; %any large-enough number
MU     = [CVpe.Segments(:).mu];

spnum  = length(DendSp);
spgpp  = cell(1, spnum);
disp(['Extracting path and salient pts for ' num2str(spnum) ' spine cndd'])

for i=1:spnum
% this will get an approximate closest mu and tip
spgpp{i}   = PathPoI(DendSp{i}, maxiterinit, MU);
if mod(i, round(spnum/10))==1
disp(['done ' num2str(i) ' of ' num2str(spnum)])
toc
end
end

////disp('Done paths')
////toc
////DBGG=0;
////save('/tmp/testpaths', 'spgpp')

if DOPATHMERGE
disp('Now attempting to merge ...')
for i=1:length(spgpp)
if ~IsMergeable(spgpp{i})///////////////////zzzzzzzzzzzzzzzzzzz
continue
end
mrgin=[];
NS = FindNeighborSpines(spgpp, i);/////////////zzzzzzzzzzzzzzzzzz
while ~isempty(NS)
j=NS(1);
if ~IsMergeable(spgpp{j})
NS = NS(2:end);
continue;
end
pm = PathMatch(spgpp([i,j]));/////////////////zzzzzzzzzzzzzzzzzzzz Frechet?
if pm
disp(['merging spines ', num2str(i), ' and ', num2str(j)])
PLOT_PATHS(spgpp,i,j);
spgpp{i}.pxl   = unique([spgpp{i}.pxl spgpp{j}.pxl]','rows')';
spgpp{i}       = PathPoI(spgpp{i}, maxiterinit, MU);
spgpp{j}.mrg2  = i;
mrgin          = [mrgin j]; %#ok<AGROW>
NS = unique([NS(2:end) FindNeighborSpines(spgpp, i)]);
NS = setdiff(NS,j); %j is already in nbrs and merged!
else
NS = NS(2:end);
end % if pm
end %while
spgpp{i}.mrgin = mrgin;
end   %for
end
DendSpF = spgpp; %(~cellfun(@(x) isinf(x.dendroot(1)), spgpp));







function [vp pathL plen gd] = ValidatePath(cp, spcell)//////////////////zzzzzzzzzzzzzzzzzzzzzzzzzzz
% cp is currentpath, fresh from get_spine_pathxx.m 
% This arg is a 2x15 cell. If the first cell is valid (cp{1,2}==0) then the
% path in the first cell is the desired one cp(1,1). Otherwise, the path
% in the first cell is the last option because it won't reach the dendrite
% as it is extracted within the spine only. But before making that the
% final answer, we had tried extracting paths from a different starting
% point such as sproot to multiple ending points in mucand. These paths are
% in cp(2,:). We pick the most valid in the sense of likelihood and
% shortest length. If they are likely valid, we pick the shortest.
global CVpe IMGpe
vp=[];pathL=-64;plen=0;gd=NaN;
if cp{2}==1
[i1,i2]=size(cp);
pfit=-64*ones(1,i2);
len=inf(1,i2);
for j=1:i2
p=cp{1,j};
if isempty(p)
continue
end
if length(p)>1
error('need to update this! using many enpts?')
end
testp=p{1};
pxl=PathPxls(testp);
pfit(j)=Check_Match4(IMGpe,CVpe, spcell.mui1, 1, pxl, 0, 1);
[dum len(j)]=vessdist(unique(round(pxl'),'rows')');
%         end
end
% get the shorter ones first and test likelihood
[sv si]=sort(len);
for j=si
if pfit(j)>-64 && len<inf
vp=cp{1,j}{1};
pathL=pfit(j);
plen = len(j);
gd=cp{3,j};
break
end
end
else
vp    = cp{1,1}{1};
pxl   = PathPxls(vp);
pathL = Check_Match4(IMGpe,CVpe, spcell.mui1, 1, pxl, 0, 1);
gd    = cp{3,1};
[dum plen] = vessdist(unique(round(pxl'),'rows')');
end


//////void UniquePointList(PointSetContainerType::Pointer  pcont)
//////{
//////	PointSetType::PointType spp;
//////	SpineImageIndexType     pxlidx;
//////	PointSetType::PointsContainerIterator   pciter = pcont->Begin();
//////	while (pciter != pcont->End()) 
//////	{
//////		spp = pciter->Value();
//////		for (int ii=0; ii< spr_Dimension; ii++)
//////			pxlidx[ii] = 
//////}














%% PathMatch check
% inputs: p1:candidate path to be merged with, p2:master path to merge
% outputs: 
%      pm: logical, true if there is a path match
%      spine: find an intersection point between the merged paths...
function pm = PathMatch(sp)///////////////////////////zzzzzzzzzzzzzzzzzzzzzzzzzzzzz
global DBGG
pm = 0;
p = {sp{1}.path, sp{2}.path};
% pxl2=sp{2}.pxl';
p1px2=0;
for i=1:2
if isempty(p{i})
return
end
if iscell(p{i})
[szi,szj]=size(p{i});
if szi>1
error('multiple path 1 ambiguity');
end
if szj~=2
error('multiple path 2 ambiguity');
end
else
error('path fortmat error')
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  TEST B  %%%%%%%%%%%%%%%%
if ~p1px2
p1p22 = TestOvrlp(p{1}{2}, p{2}{2});
if DBGG && p1p22
disp('paths overlap')
end
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pm = p1px2||p1p22;






%% FindNeighborSpines
function j = FindNeighborSpines(spgpp, i)///////////////////zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
dpxl = 20;
maxnbr=50;  %max number of neighbors
j=zeros(1,maxnbr);
jj=0;
for k=[1:i-1 i+1:length(spgpp)]
if isempty(spgpp{k})
continue
end
if isfield(spgpp{k}, 'mrg2')
continue;
end
if inregionsp(spgpp([i,k]), dpxl)
jj=jj+1;
if jj>maxnbr
warning('too many neighbors found')
break
end
j(jj)=k;
end
end
j=j(1:jj);





%% Curve Distance to check for path overlaps
function ok=TestOvrlp(c1, c2)                   //////////////////////zzzzzzzzzzzzzzzz
global PATHOVERLAPRATIO MINPATHSZ
ok=0;
p1 = unique(round(c1),'rows');
p2 = unique(round(c2),'rows');
[intpxl, idx1, idx2] = intersect(p1, p2,'rows');
sz1 = size(p1,1);
sz2 = size(p2,1);
if sz1 < MINPATHSZ
warning('path 1 too short for overlap')
return
end
if sz2 < MINPATHSZ
warning('path 2 too short for overlap')
return
end
if length(idx1)>=PATHOVERLAPRATIO*sz1 ||...
length(idx2)>=PATHOVERLAPRATIO*sz2
ok=1;
else
[mind21 m1 m2 dummaxd21 M2 MS2 m2i m2idxv mind21v] = ...
testClosest2Ptsets(p1', p2');
mind21v=sort(mind21v);
ld21=length(mind21v);
if all(mind21v(1:round(ld21*PATHOVERLAPRATIO))<=2)
ok=1;
end
end






function ok=IsMergeable(sp)////////////////////////zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
global MINPATHSZ
ok = 0;
if isempty(sp)
return
end
if isfield(sp, 'mrg2') || isfield(sp, 'invalid')
return
end
currpath = sp.path;
if isempty(currpath)
return
end
if ~iscell(currpath)
error('weirdness currpath')
end
if size(currpath{1},1)+size(currpath{2},1)<MINPATHSZ
return
end
ok=1;





function pxl=PathPxls(cp)
if iscell(cp)
if length(cp)==2
pxl=cp{2}';
else
pxl=cp{1}';
end
else
pxl=cp';
end
*/
#endif

