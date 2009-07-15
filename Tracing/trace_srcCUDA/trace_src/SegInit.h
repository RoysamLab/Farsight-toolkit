#ifndef SEGINIT_H_
#define SEGINIT_H_

#include "gpu.h"

#include "itkObjectFactory.h"
#include "itkSmartPointer.h"
#include "itkMacro.h"
#include "itkLightObject.h"

#include "itkImage.h"
#include "itkLinearInterpolateImageFunction.h"
#include "vnl/vnl_math.h"


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <cstdlib>


#include "TVessel.h"

typedef itk::Image<float , 3> ImageType3D;

#define DIM 3     /* dimension of points, must be < 31 for SIZEcube */
#define SIZEcube (1<<DIM)
#define SIZEdiamond (2*DIM)
#define TOTpoints (SIZEcube + SIZEdiamond)
#define _PI 3.14159265
#define sign(x) (x>0?1.0:-1.0)
#define _INC 19
#define MAX_SIZE 2000000
#define XR(i,j) xr[i+4*j]
#define XI(i,j) xi[i+4*j]

class SegInit	{
public:

bool fitSE (ImageType3D::Pointer, TVessel&, double , double );
void reverse(TVessel & seg);


protected:

	typedef struct _tagMatrix
	{
	    DOUBLE r[3][3];
	}TRMatrix;

	TRMatrix gRT;

	typedef struct _tagVertex
	{
	    DOUBLE p[3];
	}TVertex;

	typedef struct _tagFacet
	{
	    int vertex[3];
	    DOUBLE centroid[3];  // can leave double
	    DOUBLE normal[3];
	    DOUBLE area;         // can leave double
	}TFacet;

	typedef struct _tagTDamp
	{
	    double dt2;
	    DOUBLE dt;
	    double dt_u[3];
	    double dt_a[3];
	    double sign_S[4];
	    double sign_A[3];
	    double sign_U[3];
	}TDamp;


	TFacet  gFacets[2000];
	TVertex gVertexC[6000];
	TVertex gVertexS[6000];
	TVertex gVertexU[6000];
   DOUBLE gF[6000];

	int gNum_facets;
	int gNum_points;

	void rotation_quat(double * q, TRMatrix & R);
	void transpose_matrix(TRMatrix & Result );
	double calc_rho(double x, double y, double z, double e1, double e2 );
	void get_rotation(TVessel & vessel,TRMatrix & Result);
	double get_area(TFacet & facet);
	void generate_convex_hullq(TVessel & vessel,int NN);
	void matrix_point(TVertex & point, TRMatrix & Result );
	void matrix_point(double & x, double & y,double & z, TRMatrix & Result );
	void get_surface2(TVessel & vessel, TFacet * pFacets, int NumFacets, TVertex * pUVertices, TVertex * pSVertices, TVertex * pCVertices, int NumVertices);
    double partial_rho_e1(double x, double y, double z, double e1, double e2 );
	void update_Quaternions2(TVessel & vessel,TFacet * pFacets, int NumFacets, DOUBLE* pF, TRMatrix & R, double total_area, TDamp & damp);
	void update_e1(TVessel & vessel,TFacet * pFacets, int Num_facets,double * pF, TRMatrix & R, TDamp &damp);
	void update_mu(TVessel & vessel,TFacet * pFacets, int Num_Facets, DOUBLE* pF,  TRMatrix & R, TDamp & damp);
	void update_scale( TVessel & vessel, TFacet * pFacets,int Num_Facets, DOUBLE* pF, TRMatrix & R, TDamp & damp, double AS_RATIO );
	void axis_angle(double * q);
	void quat_mult(double * q3, double * q2, double * q1 );
	void update_axes( TVessel & vessel, int AT_END );
	double superquad3d_distq(double tx, double ty, double tz, TVessel & vessel );
	double rho_dot_asym(double x,double A,double B);
	int compare( const void * p1, const void * p2 );
	int compare_int( const void * p1, const void * p2 );

	double getMedian(std::vector<double> arr);
	double getMean(std::vector<double> arr);
	double getRegionMean (ImageType3D::Pointer im, TRMatrix& lim);
	//static void fill_array( double	*xr, double  *xi);

	int interp3(ImageType3D::Pointer im, DOUBLE* pF, TFacet * pFacets, int Num_Facets, TVertex * pVertexC, double  );
	int interp3DefBG(ImageType3D::Pointer  im, DOUBLE* pF, TFacet * pFacets, int Num_Facets, TVertex * pVertexC, double defBG );
	int handle_end( ImageType3D::Pointer im, DOUBLE* pF, TVessel & vessel, TFacet * pFacets, int Num_Facets, TVertex * pUVertices, TVertex * pSVertices, TVertex * pCVertices, int Num_Vertices  );
	double update_FB(ImageType3D::Pointer im, TVessel & vessel );
	double update_FBNew(ImageType3D::Pointer im, TVessel & vessel , int Initialize);
	double update_FBSimple(ImageType3D::Pointer im, TVessel & vessel );

public:
	SegInit()	{
		gNum_facets = 0;
		gNum_points = 0;
	}

	~SegInit()	{
	}
	
};

#endif /*SEGINIT_H_*/
