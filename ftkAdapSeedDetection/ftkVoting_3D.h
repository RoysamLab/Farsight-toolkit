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
	// The direction of the vore (quantized to the _NN_ bins)
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

	// Direction of the cone ( *** in order to speed up the process I should store of the possible dot products in a matrix)
	double dxx;
	double dyy;
	double dzz; 

};

// ############################################################################################################################################################################
/** 3D Cone without distance quantized */
struct ftkCone3D_new : public std::vector<ftkWPoint3D> { 
	// Size of the image, initialized in the cpp file
	static int nx, ny, nz;
	// !!! Counter of votes
	static int contador_1;
	ftkCone3D_new() : std::vector<ftkWPoint3D>() { }
	// Set the direction of the cone
	void setDirection(double xx, double yy, double zz); //No veo que sea util siendo que ya tenamos un vector de estoy de tamano 256
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
		for(iterator it = begin(); it!=end(); it++) { // I think this is not necessary
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

	// Direction of the cone ( *** in order to speed up the process I should store of the possible dot products in a matrix)
	double dxx;
	double dyy;
	double dzz; 

};

// ############################################################################################################################################################################
/** Main class which performs the voting process */
class ftkVoting_3D {

public:

	ftkVoting_3D();
	~ftkVoting_3D();

	/** Set the prefix of the path to store the voting landscope */
	void setPrefix(const std::string& p);

	/** Set the parameters for the algorithm */
	void setParams(	int hmin, int hmax,	int radius,	double min_grad, double scale );

	/** main user interface for voting */
	void compute(nftkVot::InputImageType_3D::Pointer I);

	/** Useful constants */
	static const double epsilon;
	static const double pi;

private:

	/** Compute cone for voting (3D cone is a triangle) */
	void computeCones_3D(int hmin, int hmax, int radius);

	/** Compute cone for voting (3D cone is a triangle) */
	void computeCones_3D_prob(int hmin, int hmax, int radius);

	/** Function that actually performs the vote operation. */
	void vote();

	/** Use a mapping for the possible directions of voting, Creates the mapping. */
	int computeDirecIndex_3D(double dx, double dy, double dz) const;

	/** This function updates the cones (reduce the span). */
	void updateCones();

	/**	Update voting direction for a given point: vp; the new direction is from vp to the point with maximum voting response within vp's voting range.*/
	void inline updateDirection_3D(VPoint3D& vp);

	/**	Update voting direction for a given point: vp; the new direction is from vp to the point with maximum voting response within vp's voting range. */
	void inline updateDirection_3D_prob(VPoint3D& vp);


	// VECTORS
	/** Voting points (tokens)*/
	std::vector< VPoint3D > _voting_points_3D;
	/** Image of votes (with the padding) */
	VotingDirType_3D::Pointer _votingSumVotes;
	/** Image of votes Mask (with the padding) */
	VotingDirType_3D::Pointer _votingMaskVotes;


	// Cones
	/** A test to do the calculation of the cones in parallel */
	std::vector < ftkCone3D > _conesPru_3D_new_para; 
	std::vector< std::vector< std::vector< int > > > _voteDirec_3D_new_para;

	/** A test to do the calculation of the cones in parallel and not quantized */
	std::vector < ftkCone3D_new > _conesPru_3D_new_notqu; 
	std::vector< std::vector< std::vector< int > > > _voteDirec_3D_new_notqu;



	std::vector< std::vector< std::vector< int > > > _voteDirec_3D_new;


	void inline votar_dir(std::vector<std::vector<int> >& p_dir,VotingDirType_3D::PixelType * p, const VPoint3D& vp, int angl_indx, int dist, int& offset_1);



	// PARAMETERS
	/** Minimumn distance of the votes */
	int _hmin;
	/** Maximum distance of the votes */
	int _hmax;
	/** Opposite catheti distance of the votes */
	int _radius;
	/** Minimum magnitude of the gradient, normlized (0-1) */
	double _min_grad;
	/** Span of the actual vote (for now is from 0-9, where 10 is the maximum span)*/
	int _intSpan;
	/** Range of the span, by defaul is 10*/
	int _maxSpan;
	/** Scale for computing the gradient using DoG */
	double _scale; 
	/** Number of quantized directions on the sphere */
	int _NN_dir;
	/** Factor of the Z direction with respect to X and Y */
	int _Z_factor;
	/** prefix for path to record voting landscope in each iteration */
	std::string _prefix;
	/** Scale of the derivative. */
	double _sigmaG;
	/** Timming variables */
	clock_t _t1_begin;
	clock_t _t2_begin;
	clock_t _t3_begin;
	clock_t _t4_begin;
	clock_t _t1_end;
	clock_t _t2_end;
	clock_t _t3_end;
	clock_t _t4_end;


	// ??? these variables can be private ?
	/** x-size */
	int nx; 
	/** y-size */
	int ny;
	/** z-size */
	int nz;
	/** Size of the image */
	int npix;
	/** padding size for border pixels during voting  */
	int bw;
	/** Original size. */
	int _bw;
	int _nx_org;
	int _ny_org;
	int _nz_org;
	/** Original size */
	int _npix_org;


//	/** voting */

	void inline votar(VotingDirType_3D::PixelType * p, const VPoint3D& vp, int angl_indx);
	void inline votar(VotingDirType_3D::PixelType * p, const VPoint3D& vp, int angl_indx, int dist);
	void inline votar_prob(VotingDirType_3D::PixelType * p, const VPoint3D& vp, int angl_indx, int dist);

//	PointArray _centers; //< detected centers of objects
//	std::vector<double> _center_weights; //< accumulated votes for each center	


	//< pre-computed cone structures for voting
	std::vector<std::vector<ftkCone3D> > _conesPru_3D;
	std::vector < ftkCone3D > _conesPru_3D_new;
	std::vector<ftkCone3D> _conesPru_prob;
	std::vector< std::pair< int,int > > _voteDirec;
	std::vector< std::pair< std::pair< int,int > , std::vector <int> > > _voteDirec_prob;

	// Copy of the pointer of the input image
	nftkVot::InputImageType_3D::Pointer _inputImage;


	};

#endif
// ############################################################################################################################################################################