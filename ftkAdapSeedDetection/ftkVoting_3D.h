// ############################################################################################################################################################################
#ifndef FTKVOTING_3D_H
#define FTKVOTING_3D_H 

// ############################################################################################################################################################################
// ??? -> Questions
// !!! -> Testing
// *** -> Someting to be improved
// ############################################################################################################################################################################

#include "ftkVotingGlobal.h"

// Voting Direction X and Y
typedef double VotingDirPixelType;
typedef itk::Image< VotingDirPixelType, 3 > VotingDirType_3D;

// Direction per pixel ( I think this was tested for storing the direction of each one of the votes to each pixel)
typedef int VotingDirPerPixelType;
typedef itk::Image< VotingDirPerPixelType, 2 > VotingDirPerType_scalar;
typedef itk::VectorImage< VotingDirPerPixelType, 2 > VotingDirPerType;
typedef itk::ComposeImageFilter<VotingDirPerType_scalar> ImageToVectorImageFilterType;



// ############################################################################################################################################################################
/** Voting point */
struct VPoint3D {
	// ???
	short x, y, z;		
	// Center point
	short xc, yc, zc;	
	// Actual position on the data structure ("cube")
	int pos;			
	// Magnitude of the vote
	double mag;			
	// Scale of the vote
	int scale;			
	// The direction of the vore
	int direc_vote;		
	// Position on the original data structure ("cube")
	int posOri;			
};

// ############################################################################################################################################################################
/** Weighted point */
struct ftkWPoint3D {    
	// Size of the image, initialized in the cpp file
	static int nx, ny, nz;		
	// Constructor, initialize the variables
	ftkWPoint3D(int xx=0, int yy=0, int zz=0, double ww=0) : x(xx), y(yy), z(zz), off(0), w(ww) {}		
	// Set image size (nx, ny, nz are static)
	static void setImageSize(int xsize, int ysize, int zsize)
	{  
		nx = xsize;
		ny = ysize;
		nz = zsize;
	}
	// To avoid having to calculate x+y*nx+z*nx*ny for every access of the votes
	void setOffset()			
	{  
		off = x + y*nx + z*nx*ny; 
	}
	// ???
	int x, y, z;				
	// Actual position on the data structure ("cube")
	int off;					
	// Weight of the vote
	double w;					
};

// ############################################################################################################################################################################
/** Bins of the cones */
struct ftkBins3D : public std::vector<ftkWPoint3D> {
	// Size of the image, initialized in the cpp file
	static int nx, ny, nz;		
	// Constructor, initialize the variables
	ftkBins3D(int n, const ftkWPoint3D &v) : std::vector<ftkWPoint3D>(n, v) { }		
	ftkBins3D() : std::vector<ftkWPoint3D>() { }
	// Set image size (nx, ny, nz are static)
	static void setImageSize(int xsize, int ysize, int zsize)
	{  
		nx = xsize;
		ny = ysize;
		nz = zsize;
	}
	// To avoid having to calculate x+y*nx+z*nx*ny for every access of the votes
	void setOffset()
	{
		for(iterator it = begin(); it!=end(); it++)
			it->setOffset();     
	}
	// !!!
	void nic (){}
};

// ############################################################################################################################################################################
/** 3D Cone */
struct ftkCone3D : public std::vector<ftkBins3D> {
	// Size of the image, initialized in the cpp file
	static int nx, ny, nz; //< image size
	// !!! Counter of votes
	static int contador_1;
	ftkCone3D() : std::vector<ftkBins3D>() { }
	// Set the direction of the cone
	void setDirection(double xx, double yy, double zz)
	{
		dxx = xx;
		dyy = yy;
		dzz = zz;
	}
	// Set the image size (nx, ny, nz are static)
	static void setImageSize(int xsize, int ysize, int zsize)
	{  
		nx = xsize;
		ny = ysize;
		nz = zsize;
	}
	// To avoid having to calculate x+y*nx+z*nx*ny for every access of the votes
	void setOffset()
	{
		for(iterator it = begin(); it!=end(); it++) {
			it->setOffset();
		}
	}
	// Functions to actually vote vote
	//void inline vote(VotingDirType_3D::PixelType * p, const VPoint3D& vp);
	//void inline vote(VotingDirType_3D::PixelType * p, const VPoint3D& vp, int &dist);
	//void inline vote(VotingDirType_3D::PixelType * p, const VPoint3D& vp, int& dist, int &mag);

	// !!! Store the direction of the votes
	//void inline vote_dir(VotingDirType_3D::PixelType * p, const VPoint3D& vp);
	// !!! Store the direction of the votes
	void inline vote_dir(std::vector<std::vector<int> >& p_dir, VotingDirType_3D::PixelType * p, const VPoint3D& vp, int& dist, int& offset_1);

	// Direction of the cone (
	double dxx;
	double dyy;
	double dzz; 

};

// ############################################################################################################################################################################
/** \brief A 3D cone */
struct ftkCone3D_new : public std::vector<ftkWPoint3D> { // Cones without distance quantized
	static int nx, ny, nz; //< image size
	static int contador_1;
	ftkCone3D_new() : std::vector<ftkWPoint3D>() { }
	void setDirection(double xx, double yy, double zz); //No veo que sea util siendo que ya tenamos un vector de estoy de tamano 256
	static void setImageSize(int xsize, int ysize, int zsize)
	{  
		nx = xsize;
		ny = ysize;
		nz = zsize;
	}
	void setOffset()
	{
		for(iterator it = begin(); it!=end(); it++) { // I think this is not necessary
			it->setOffset();
		}
	}
	//void inline vote(VotingDirType_3D::PixelType * p, const VPoint3D& vp);
	//void inline vote(VotingDirType_3D::PixelType * p, const VPoint3D& vp, int &dist);
	//void inline vote(VotingDirType_3D::PixelType * p, const VPoint3D& vp, int& dist, int &mag);

	//void inline vote_dir(VotingDirType_3D::PixelType * p, const VPoint3D& vp); // Guarda direccion 
	//void inline vote_dir(std::vector<std::vector<int> >& p_dir, VotingDirType_3D::PixelType * p, const VPoint3D& vp, int& dist, int& offset_1); // Guarda direccion

	//void vote(VotingDirType_3D::PixelType * p, VotingDirType_3D::PixelType * pmask, const VPoint3D& vp);
	double dx, dy, dz; //< direction

	double dxx;
	double dyy;
	double dzz; // This is the direction of the center pixel

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
	void setPrefix(const std::string& p);

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
	std::pair<int,int> computeAngleIndex_3D(double dx, double dy, double dz) const;


//	/** voting */
	void vote();
	void inline votar(VotingDirType_3D::PixelType * p, const VPoint3D& vp, int angl_indx);
	void inline votar(VotingDirType_3D::PixelType * p, const VPoint3D& vp, int angl_indx, int dist);
	void inline votar_prob(VotingDirType_3D::PixelType * p, const VPoint3D& vp, int angl_indx, int dist);
	void updateCones();

//	PointArray _centers; //< detected centers of objects
//	std::vector<double> _center_weights; //< accumulated votes for each center
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
	std::string _prefix; //< prefix for path to record voting landscope in each iteration
	int nx; // x-size
	int ny; // y-size
	int nz; // z-size
	int npix; // size of image
	int bw;//< padding size for border pixels during voting 

	//// Vectores
	std::vector< VPoint3D > _voting_points; // voting candidate points
	//std::vector<Cone3D> _coe::Pointer _votingVotes;

	// Image of votes (with the padding)
	VotingDirType_3D::Pointer _votingSumVotes;

	// Nose que es
	VotingDirType_3D::Pointer _votingMaskVotes;

	//typename VotingDirType_3D::Pointer DerivateImage( int dIrection );

	//std::vector<std::vector<int> > _votingMaskVotes_dir;

	void inline votar_dir(std::vector<std::vector<int> >& p_dir,VotingDirType_3D::PixelType * p, const VPoint3D& vp, int angl_indx, int dist, int& offset_1);


	//< pre-computed cone structures for voting
	std::vector<std::vector<ftkCone3D> > _conesPru_3D;
	std::vector < ftkCone3D > _conesPru_3D_new;
	std::vector<ftkCone3D> _conesPru_prob;
	std::vector< std::pair< int,int > > _voteDirec;
	std::vector< std::pair< std::pair< int,int > , std::vector <int> > > _voteDirec_prob;


	std::vector< std::vector< std::vector< int > > > _voteDirec_3D_new;

	//// Voting Direction X and Y
	//typedef double VotingDirPixelType;
	//typedef itk::Image< VotingDirPixelType, 2 > VotingDirType_3D;

	// Image of votes (without padding)
	//VotingDirTyp

	int _NN_dir;
	double _sigmaG;



	// Cones
	std::vector < ftkCone3D > _conesPru_3D_new_para; // A test to do the calculation of the cones in parallel
	std::vector< std::vector< std::vector< int > > > _voteDirec_3D_new_para;

	std::vector < ftkCone3D_new > _conesPru_3D_new_para_notqu; // A test to do the calculation of the cones in parallel and not quantized
	std::vector< std::vector< std::vector< int > > > _voteDirec_3D_new_para_notqu;


	clock_t _t1_begin;
	clock_t _t1_end;


	std::vector< VPoint3D > _voting_points_3D; // List of voting points

	int computeDirecIndex_3D(double dx, double dy, double dz) const;

	int _Z_factor;

	};

#endif
// ############################################################################################################################################################################