// ############################################################################################################################################################################
#ifndef FTKVOTING_3D_H
#define FTKVOTING_3D_H 

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

#include "itkRecursiveGaussianImageFilter.h"
#include "itkImageDuplicator.h"


//#include <image.h>
//#include <point.h>
//#include <global.h>

//namespace stdi {

//	namespace voting {

using namespace std;

// Voting Direction X and Y
typedef double VotingDirPixelType;
typedef itk::Image< VotingDirPixelType, 3 > VotingDirType;

// Direction per pixel ( I think this was tested for storing the direction of each one of the votes to each pixel)
typedef int VotingDirPerPixelType;
typedef itk::Image< VotingDirPerPixelType, 2 > VotingDirPerType_scalar;
typedef itk::VectorImage< VotingDirPerPixelType, 2 > VotingDirPerType;
typedef itk::ComposeImageFilter<VotingDirPerType_scalar> ImageToVectorImageFilterType;



// Output Image Type
typedef unsigned char OutputPixelType;
typedef itk::Image< OutputPixelType, 2 > OutputImageType;





// ############################################################################################################################################################################
/** Voting point */
struct VPoint3D {
	short x, y, z; //(PARA QUE LA DIFERENCIA ENTRE X,Y Y XC Y XC)
	short xc, yc, zc; //< center point (la nueva direccino a la que apunta)
	pair<int,int> angIndex; //< angle index
	int pos; //PARA QUE ES POS
	double mag; //< magnitude
	int scale;
};

// ############################################################################################################################################################################
/** Weighted point */
struct ftkWPoint3D {    
	static int nx, ny, nz; // SE TIENEN QUE DEFINIR EN ALGUN PUNTO EN EL CPP "int WPoint3D::nx = 0;"
	ftkWPoint3D(int xx=0, int yy=0, int zz=0, double ww=0) : x(xx), y(yy), z(zz), off(0), w(ww) {}
	static void setImageSize(int xsize, int ysize, int zsize)
	{  
		nx = xsize;
		ny = ysize;
		nz = zsize;
	}
	void setOffset() 
	{  off = x + y*nx + z*nx*ny; }
	int x, y, z; //< coordinates
	int off; //< offset position (//PARA QUE ES ESTE PARAMETRO)
	double w; //< weight
};

// ############################################################################################################################################################################
/** An intersection of cone */
//  typedef vector<WPoint3D> ConePlane3D;  
struct ftkBins3D : public vector<ftkWPoint3D> {
	static int nx, ny, nz;
	ftkBins3D(int n, const ftkWPoint3D &v) : vector<ftkWPoint3D>(n, v) { }
	ftkBins3D() : vector<ftkWPoint3D>() { }
	static void setImageSize(int xsize, int ysize, int zsize)
	{  
		nx = xsize;
		ny = ysize;
		nz = zsize;
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
/** \brief A 3D cone */
struct ftkCone3D : public vector<ftkBins3D> {
	static int nx, ny, nz; //< image size
	static int contador_1;
	ftkCone3D() : vector<ftkBins3D>() { }
	void setDirection(double xx, double yy, double zz); //No veo que sea util siendo que ya tenamos un vector de estoy de tamano 256
	static void setImageSize(int xsize, int ysize, int zsize)
	{  
		nx = xsize;
		ny = ysize;
		nz = zsize;
	}
	void setOffset()
	{
		for(iterator it = begin(); it!=end(); it++) {
			it->setOffset();
		}
	}
	void inline vote(VotingDirType::PixelType * p, const VPoint3D& vp);
	void inline vote(VotingDirType::PixelType * p, const VPoint3D& vp, int &dist);
	void inline vote(VotingDirType::PixelType * p, const VPoint3D& vp, int& dist, int &mag);

	//void inline vote_dir(VotingDirType::PixelType * p, const VPoint3D& vp); // Guarda direccion 
	void inline vote_dir(vector<vector<int> >& p_dir, VotingDirType::PixelType * p, const VPoint3D& vp, int& dist, int& offset_1); // Guarda direccion

	//void vote(VotingDirType::PixelType * p, VotingDirType::PixelType * pmask, const VPoint3D& vp);
	double dx, dy, dz; //< direction
};

// ############################################################################################################################################################################
class ftkVoting_3D {

public:

	static const double epsilon;
	static const double pi;

	ftkVoting_3D();
	~ftkVoting_3D();

	void setParams(	int hmin, int hmax,	int radius,	double min_grad, double scale );

	/** main user interface for voting */
	void compute(nftkVotingGlobal::InputImageType_3D::Pointer I);

	/**
	Set the prefix of the path to store the voting landscope
	*/
	void setPrefix(const string& p);

private:

	/**
	Compute voting aperture in the next iteration
	\param h voting range
	\param r current voting aperture
	\return voting aperture in the next iteration
	*/
	//int nextConeRadius(int h, int r);


	/** 
	Compute cone for voting (3D cone is a triangle)
	\param hmin lower bound of voting range, measured in pixel
	\param hmax upper bound of voting range, measured in pixel, usually set as 3/4 of the object diameter
	*/
	void computeCones(int hmin, int hmax, int radius);
	void computeCones_prob(int hmin, int hmax, int radius);

	/**
	Update voting direction for a given point: vp; the new direction is from vp to the point with maximum voting response within vp's voting range.
	*/
	//void getMaxDirection(VPoint3D& vp);
	void inline updateDirection(VPoint3D& vp);

	/**
	Update voting direction for a given point: vp; the new direction is from vp to the point with maximum voting response within vp's voting range.
	*/
	//void getMaxDirection(VPoint3D& vp);
	void inline updateDirection_prob(VPoint3D& vp);

//	/**
//	Update voting direction for a given point: vp; the new direction is from vp to the point (x,y) within vp's voting range, in which (x,y) is the average coordinate value weighted by its votes at current iteration.
//	*/
//	void getMeanDirection(VPoint3D& vp);
//
	/**
	\return The angle index in the look-up table
	*/
	int computeAngleIndex(double dx, double dy) const;
	pair<int,int> computeAngleIndex_3D(double dx, double dy, double dz) const;


//	/** voting */
	void vote();
	void inline votar(VotingDirType::PixelType * p, const VPoint3D& vp, int angl_indx);
	void inline votar(VotingDirType::PixelType * p, const VPoint3D& vp, int angl_indx, int dist);
	void inline votar_prob(VotingDirType::PixelType * p, const VPoint3D& vp, int angl_indx, int dist);
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
	int nz; // z-size
	int npix; // size of image
	int bw;//< padding size for border pixels during voting 

	//// Vectores
	vector<VPoint3D> _voting_points; // voting candidate points
	//vector<Cone3D> _coe::Pointer _votingVotes;

	// Image of votes (with the padding)
	VotingDirType::Pointer _votingSumVotes;

	// Nose que es
	VotingDirType::Pointer _votingMaskVotes;

	//typename VotingDirType::Pointer DerivateImage( int dIrection );

	//vector<vector<int> > _votingMaskVotes_dir;

	void inline votar_dir(vector<vector<int> >& p_dir,VotingDirType::PixelType * p, const VPoint3D& vp, int angl_indx, int dist, int& offset_1);


	//< pre-computed cone structures for voting
	vector<ftkCone3D> _conesPru;
	vector<ftkCone3D> _conesPru_prob;
	vector< pair< int,int > > _voteDirec;
	vector< pair< pair< int,int > , vector <int> > > _voteDirec_prob;

	//// Voting Direction X and Y
	//typedef double VotingDirPixelType;
	//typedef itk::Image< VotingDirPixelType, 2 > VotingDirType;

	// Image of votes (without padding)
	//VotingDirTyp
	};

#endif
// ############################################################################################################################################################################