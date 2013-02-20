// ############################################################################################################################################################################
#ifndef FTKVOTING_H
#define FTKVOTING_H 

#include<iostream>
#include<vector>
//#include<stringstream>

#include"ftkVotingGlobal.h"

#include "itkSobelOperator.h"
#include "itkNeighborhoodInnerProduct.h"
#include "itkMinimumMaximumImageCalculator.h"

#include "itkComposeImageFilter.h"
#include "itkVectorImage.h"

#include "itkCannyEdgeDetectionImageFilter.h"


//#include <image.h>
//#include <point.h>
//#include <global.h>

//namespace stdi {

//	namespace voting {

using namespace std;

// Voting Direction X and Y
typedef double VotingDirPixelType;
typedef itk::Image< VotingDirPixelType, 2 > VotingDirType;

// Direction per pixel
typedef int VotingDirPerPixelType;
typedef itk::Image< VotingDirPerPixelType, 2 > VotingDirPerType_scalar;
typedef itk::VectorImage< VotingDirPerPixelType, 2 > VotingDirPerType;
typedef itk::ComposeImageFilter<VotingDirPerType_scalar> ImageToVectorImageFilterType;



// Output Image Type
typedef unsigned char OutputPixelType;
typedef itk::Image< OutputPixelType, 2 > OutputImageType;

// ############################################################################################################################################################################
/** Voting point */
struct VPoint2D {
	short x, y; //(PARA QUE LA DIFERENCIA ENTRE X,Y Y XC Y XC)
	short xc, yc; //< center point (la nueva direccino a la que apunta)
	short angIndex; //< angle index
	int pos; //PARA QUE ES POS
	double mag; //< magnitude
	int scale;
};

// ############################################################################################################################################################################
/** Weighted point */
struct ftkWPoint2D {    
	static int nx, ny; // SE TIENEN QUE DEFINIR EN ALGUN PUNTO EN EL CPP "int WPoint2D::nx = 0;"
	ftkWPoint2D(int xx=0, int yy=0, double ww=0) : x(xx), y(yy), off(0), w(ww) {}
	static void setImageSize(int xsize, int ysize)
	{  
		nx = xsize;
		ny = ysize;
	}
	void setOffset() 
	{  off = x + y*nx; }
	int x, y; //< coordinates
	int off; //< offset position (//PARA QUE ES ESTE PARAMETRO)
	double w; //< weight
};

// ############################################################################################################################################################################
/** An intersection of cone */
//  typedef vector<WPoint2D> ConePlane2D;  
struct ftkBins2D : public vector<ftkWPoint2D> {
	static int nx, ny;
	ftkBins2D(int n, const ftkWPoint2D &v) : vector<ftkWPoint2D>(n, v) { }
	ftkBins2D() : vector<ftkWPoint2D>() { }
	static void setImageSize(int xsize, int ysize)
	{  
		nx = xsize;
		ny = ysize;
	}
	void setOffset()
	{
		for(iterator it = begin(); it!=end(); it++) //NO ENTIENDO SU FUNCION MUY BIEN
			it->setOffset();     
	}
	void nic (){}
	//double center_weight;
};

// ############################################################################################################################################################################
/** \brief A 2D cone */
struct ftkCone2D : public vector<ftkBins2D> {
	static int nx, ny; //< image size
	static int contador_1;
	ftkCone2D() : vector<ftkBins2D>() { }
	void setDirection(double xx, double yy); //No veo que sea util siendo que ya tenamos un vector de estoy de tamano 256
	static void setImageSize(int xsize, int ysize)
	{  
		nx = xsize;
		ny = ysize;
	}
	void setOffset()
	{
		for(iterator it = begin(); it!=end(); it++) {
			it->setOffset();
		}
	}
	void inline vote(VotingDirType::PixelType * p, const VPoint2D& vp);
	void inline vote(VotingDirType::PixelType * p, const VPoint2D& vp, int &dist);
	void inline vote(VotingDirType::PixelType * p, const VPoint2D& vp, int& dist, int &mag);

	//void inline vote_dir(VotingDirType::PixelType * p, const VPoint2D& vp); // Guarda direccion 
	void inline vote_dir(vector<vector<int> >& p_dir, VotingDirType::PixelType * p, const VPoint2D& vp, int& dist, int& offset_1); // Guarda direccion

	//void vote(VotingDirType::PixelType * p, VotingDirType::PixelType * pmask, const VPoint2D& vp);
	double dx, dy; //< direction
};

// ############################################################################################################################################################################
class ftkVoting {

public:

	static const double epsilon;
	static const double pi;

	ftkVoting();
	~ftkVoting();

	void setParams(	int hmin, int hmax,	int radius,	double min_grad, double scale );

	/** main user interface for voting */
	void compute(nftkVot::InputImageType::Pointer I);

	/**
	Set the prefix of the path to store the voting landscope
	*/
	void setPrefix(const string& p);


	// -------------------------------------
	void SetSlice( unsigned int slice){
		_slice = slice;
	}
	unsigned int _slice;

	void setOutBySlice( nftkVot::InputImageType_3D::Pointer input ){
		_votingBySlice = input;
	}
	nftkVot::InputImageType_3D::Pointer _votingBySlice;

	void setOutBySliceLastbigSapan( nftkVot::InputImageType_3D::Pointer input ){
		_votingBySliceLastbigSapan = input;
	}
	nftkVot::InputImageType_3D::Pointer _votingBySliceLastbigSapan;

	void setOutBySliceLastbigSapanProb( nftkVot::InputImageType_3D::Pointer input ){
		_votingBySliceLastbigSapanProb = input;
	}
	nftkVot::InputImageType_3D::Pointer _votingBySliceLastbigSapanProb;
	// -------------------------------------

	double _OriginalScale;


private:

	/**
	Compute voting aperture in the next iteration
	\param h voting range
	\param r current voting aperture
	\return voting aperture in the next iteration
	*/
	int nextConeRadius(int h, int r);


	/** 
	Compute cone for voting (2D cone is a triangle)
	\param hmin lower bound of voting range, measured in pixel
	\param hmax upper bound of voting range, measured in pixel, usually set as 3/4 of the object diameter
	*/
	void computeCones(int hmin, int hmax, int radius);
	void computeCones_prob(int hmin, int hmax, int radius);

	/**
	Update voting direction for a given point: vp; the new direction is from vp to the point with maximum voting response within vp's voting range.
	*/
	//void getMaxDirection(VPoint2D& vp);
	void inline updateDirection(VPoint2D& vp);

	/**
	Update voting direction for a given point: vp; the new direction is from vp to the point with maximum voting response within vp's voting range.
	*/
	//void getMaxDirection(VPoint2D& vp);
	void inline updateDirection_prob(VPoint2D& vp);

//	/**
//	Update voting direction for a given point: vp; the new direction is from vp to the point (x,y) within vp's voting range, in which (x,y) is the average coordinate value weighted by its votes at current iteration.
//	*/
//	void getMeanDirection(VPoint2D& vp);
//
	/**
	\return The angle index in the look-up table
	*/
	int computeAngleIndex(double dx, double dy) const;
//	/** voting */
	void vote();
	void inline votar(VotingDirType::PixelType * p, const VPoint2D& vp, int angl_indx);
	void inline votar(VotingDirType::PixelType * p, const VPoint2D& vp, int angl_indx, int dist);
	void inline votar_prob(VotingDirType::PixelType * p, const VPoint2D& vp, int angl_indx, int dist);
	void updateCones();

//	PointArray _centers; //< detected centers of objects
//	vector<double> _center_weights; //< accumulated votes for each center
//
//	doubleImage _sum, _mask, _vote; //<sum of voting
//
//	
//
//	
//	int bw;//< padding size for border pixels during voting 
//
//


	// PARAMETERS
	//reg_type (have not underestand what is this for)
	int _hmin;
	int _hmax;
	int _radius;
	double _min_grad;
	int _intSpan;
	//_threshold (for picking seed points, not necesary for now)
	double _scale; // Scale for computing the gradient using DoG
	//_zc_only TRUE - voting points should also be the zero-crossing points of the image; FALSE - otherwise

	// Default Parameters for quantizing the direction of voting
	int ntheta;
	double delta_theta;

	//INTERNAL
	string _prefix; //< prefix for path to record voting landscope in each iteration
	int nx; // x-size
	int ny; // y-size
	int npix; // size of image
	int bw;//< padding size for border pixels during voting 

	//// Vectores
	vector<VPoint2D> _voting_points; // voting candidate points
	//vector<Cone2D> _coe::Pointer _votingVotes;

	// Image of votes (with the padding)
	VotingDirType::Pointer _votingSumVotes;

	// Nose que es
	VotingDirType::Pointer _votingMaskVotes;

	//typename VotingDirType::Pointer DerivateImage( int dIrection );

	//vector<vector<int> > _votingMaskVotes_dir;

	void inline votar_dir(vector<vector<int> >& p_dir,VotingDirType::PixelType * p, const VPoint2D& vp, int angl_indx, int dist, int& offset_1);


	//< pre-computed cone structures for voting
	vector<ftkCone2D> _conesPru;
	vector<ftkCone2D> _conesPru_prob;
	vector< pair< int,int > > _voteDirec;
	vector< pair< int,int > > _voteDirec_copy;
	vector< pair< pair< int,int > , vector <int> > > _voteDirec_prob;
	vector< pair< pair< int,int > , vector <int> > > _voteDirec_prob_copy;

	//// Voting Direction X and Y
	//typedef double VotingDirPixelType;
	//typedef itk::Image< VotingDirPixelType, 2 > VotingDirType;

	// Image of votes (without padding)
	//VotingDirTyp



	

	};

#endif
// ############################################################################################################################################################################