//#include "itkPointSet.h"
//// Define Objects
//typedef itk::PointSet <unsigned short, 3> TracePointSet;
//typedef TracePointSet::PointType TracePoint;

#ifndef SPINERING_H
#define SPINERING_H

#include <itkImage.h>
#include "itkObjectFactory.h"
#include "itkMacro.h"
#include "itkLightObject.h"
#include <itkSmartPointer.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageLinearConstIteratorWithIndex.h>
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
//#include "itkImageRegionIterator.h"
//#include "itkSimpleFilterWatcher.h"
#include "itkImageRegionIterator.h"
//#include "itkLinearInterpolateImageFunction.h"
#include "itkPointSet.h"
#include "itkPolyLineParametricPath.h"
//#include "vnl/vnl_sample.h"
#include "SpineConsts.h"
#include "CommonTypeDefs.h"
#include "SpineTypedefs.h"
#include "SpineUtils.h"
#include "TraceSegNode.h"

//class DetectorQ {
//public:
//	DetectorQ();
//	~DetectorQ();
//	NewRing(SpineRing*, IndexVecType, IndexVecType);
//	Flush();
//private:
//	IndexVecType qminidx;
//	IndexVecType qmaxidx;
//	std::vector<SpineRing*> Q;
//	short size;
//};
class ImageDebugger;
class GlobalDetector;

inline
double ED(PointType p1, PointType p2) {
	return sqrt(pow(p1[0]-p2[0],2) + pow(p1[1]-p2[1],2) + pow(p1[2]-p2[2],2));
}

struct PointSetDistS {
	PointSetDistS(PointSetType::Pointer  s1, PointSetType::Pointer  s2);
	PointSetDistS(PointType p1, PointSetType::Pointer  s2);
	PointSetType::Pointer  set1, set2;
	PointSetContainerType::Pointer c1, c2;
	PointSetContIterType        i1, i2;
	PointType  pt1m, pt2m, pt1x, pt1x2;
	int minidx2; //idx of pt2m 
	double FD, HD;
	std::vector<int>minidx21a;//closest pt indxs in 2 to 1; size of set1
	std::vector<double>mind21a;//min distances 2 to 1;   size of set1
	double dmin, dmax; 
	void Compute();
	void Print();
};



class SpCandidate {
public:
	SpCandidate(GlobalDetector* gd, int trID, unsigned short cid);
	~SpCandidate();
	void            AddPoint(SpineImageIndexType pxl, SpineImageIndexType bbmin);
	void			GetClosestMUs();
	void            GetNbrMUs(TraceSegNode *seg, int count);
	bool			Validate();
	void            PrintSelf();
	void            PathPtsOfInt();
	void			RegionToPath();
	unsigned short	GetcandID() {return candID;}
	int				GettraceID() {return traceID;}
	int				WritePxlsToRGBImage(RGBImageType::Pointer rgbim2D, RGBPixelType color);
	void			WritePoIToRGBImage(RGBImageType::Pointer rgbim2D);
private:
	int				spinesize; // pidx;
	PointSetType::Pointer	pxls;
	PointSetContainerType::Pointer	pcont;
	PointSetType::PointDataContainer::Pointer pdata;

	bool			valid;
	bool            validPath;
	GlobalDetector  *ParentDet;
	unsigned short	candID;

	int             traceID;
	ImagePixelType  sigmaI;
	unsigned short	ptID;
	PointSetType::Pointer	NbrMus;
	PointSetContainerType::Pointer	NbrMusPtCont;
	std::vector<int> NbrMuIDVec; //temp holder non unique
	PointSetDistS    *PSDists;
	PointType        brtpt, bmu;
	std::vector<int> CloseNodeidx;// replaces mucandidx
	TraceSegNode    *munode;
	SpineImageIndexType	scbbmin, scbbmax; //bounding box
	SpineImageType::SizeType imgsize;
	//mu;X;X2;path;
	PathType::Pointer				Path;
	SpineImageType::RegionType		GetSpineRegion();
	SpineImageType::RegionType		SpCandRegion;
	void							ComputeSpineRegion() ;
	SpeedImageType::Pointer			GetSpeedFunc(SpeedImageType::Pointer Spregion);
									/*, 	bool UseGradient = false,
											bool EnhanceContrast=true, 
											bool SmoothImage=true,
											const double GaussSigma = 1.0,
											const double SigmoidAlpha = -.5,
											const double SigmoidBeta  = 3.0
									 */
	PathType::Pointer				SpeedToFMPath(SpeedImageType::Pointer,
											PointType, 
											PointSetContainerType::Pointer, 
											std::string , 
											RGBImageType::Pointer rgbim2d);
	//PathType::Pointer				FMPath(SpineImageType::Pointer, PointType, PointSetContainerType::Pointer);
	PointType						TranslatePts(PointType pt0,  short dir);
	PointSetContainerType::Pointer	TranslatePts(PointSetContainerType::Pointer cont, short dir);
	//SpineImageType::Pointer			GetSpeedFunction(SpineImageType::RegionType spregion);
	SpeedImageType::Pointer			GetSpeedSubIm();
};


class GlobalDetector {
// manages spine detection on entire image
public:
	GlobalDetector(SpineImageType::Pointer im, TraceContainer *TCP, ImageDebugger* imdebugger);
	~GlobalDetector(){};
	void Run();
	BinaryImageType::Pointer CCim;
	SpineImageType::Pointer  image; //spine image
	
	ImageDebugger			*imdbg;
	TraceSegNodeVecType     *NodeList;
	const PointSetPointerMap* DendMuMap;
	TraceIDVecType          *TraceIDList;
	TraceSegNode			*getSegment(long i);
	void					DoCC(unsigned short trID);
	void                    ResetCCim();
	void					ValidateSpCands();
	void                    PrintDendSpCandMap(unsigned short trID);
	void                    PrintCCim();
	RGBImageType::Pointer   GetRGBPOImage();
	RGBImageType::Pointer   GetRGBPathsImage();
	SpineImageType::SizeType GetImageSize() {return imgsize;};
	void                    GetSpCandStats();
	RGBPixelType			GetColorMap(int loc);
	unsigned short			GetSpCandTotal() {return spcandtotal;}
	int						WriteDendSpCandMap2RGB(RGBImageType::Pointer rgbim2D);
	const double			*up_sampling;// = UP_SAMPLING;
private:
	
	// this is the hierarchical spine candidate list
	// as a map indexed by the dendrite ID
	// each map element points to the sub-map of spine 
	// candidates for that dendrite ID. 
	// The reason why the spine candidates are in a map 
	// rather than a list is the way they get filled when
	// iterating in the relabeled CCim.
	std::map<unsigned short, SpCandMapType*>    DendSpCandMap;
	SpineImageType::SizeType					imgsize;
	SpineImageIndexType					        bbmin, bbmax;
	RGBImageType::Pointer						RGBPathsImage;
	std::vector<RGBPixelType>					colormap;
	unsigned short								maxmapsz, spcandtotal;
	//RingQ *Q;
	//TraceSegNode *seg;
	unsigned int VID;
	unsigned int segID;
	class SpineRing	{  //: public SegInit {
	//single ring detector class
		public:
			SpineRing(GlobalDetector*, spr_SpineRingType);
			~SpineRing();
			GlobalDetector *Parent;
			void PrintSelf();
			void   CCimUpdate();
			inline float getForegdIntensity();// {return(this->ForegdIntensity);}
			inline float getBackgdIntensity();// {return(this->BackgdIntensity);}
			//inline double getriscale() {return (this->radiscale);}
			//inline double getroscale() {return (this->radoscale);}
			inline void setRingType(spr_SpineRingType i){RingType = i;}
			inline spr_SpineRingType	getRingType ()	{return RingType;}
			inline unsigned int		getTraceID()	{return(seg->TraceID);}
			inline unsigned int		getSegID()		{return(seg->ID);}
			bool					RingSE();
			void                    RestartRing(TraceSegNode*);
			bool					validring;
			// HUSSEIN: Can be implemented later as well
			// bool RingCyl();
		private:
			TraceSegNode*		seg; //TraceBitp; this is the underlying trace segment
			double				thickness;
			// to implement up_sampling properly, may need to 
			// make UP_SAMPLING into an argument
			std::vector<ImagePixelType> inring_f_vals;
			spr_SpineRingType			RingType;
			IndexVecType				inring_f_idx;
			DEBUGSTMT2(IndexVecType	,	inring_all_idx);
	};
};



//struct PointSetDistS {
//	std::vector< InterpPointType >* set1, *set2;
//	InterpPointType  pt1m, pt2m, pt1x, pt1x2;
//	int minidx2; //idx of pt2m 
//	std::vector<int>minidx21a;//closest pt indxs in 2 to 1; size of set1
//	std::vector<double>mind21a;//min distances 2 to 1;   size of set1
//	double dmin, dmax, mind21;
//	void Compute();
//};
//class SERing : public SpineRing {
//public:
//	SERing();
//	~SERing();
//	void SetRscale(double s) {rscale=s;};
//	void SetDist2Center(double d) {dist2center = d;};
//	void SetThickness(double t) {thisckness = t;};
//	void SetUp_Sampling(double* up) {up_sampling = up;};
//	void SetSeg(TVessel* s) {seg = s;};
//private:
//	double rscale;
//	double dist2center;
//    double *up_sampling;
//	double thickness;
//	//bool seout;
//	TVessel* seg;
//}

#endif
