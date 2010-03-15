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

/**
 \brief Class to initialize the superellipse from seed point. 
 \author $ Author: James Alex Tyrrell, Amit Mukherjee $
 \version $Revision 1.0 $
 \date May 05 2009
*/
// Copyright Farsight Toolkit , Rensselaer Polytechnic institute Troy NY 12180.

#ifndef SEGINIT_H_
#define SEGINIT_H_

#include "itkObjectFactory.h"
#include "itkSmartPointer.h"
#include "itkMacro.h"
#include "itkLightObject.h"

#include "itkImage.h"
#include "itkLinearInterpolateImageFunction.h"
#include "vnl/vnl_math.h"
#include "itkContinuousIndex.h"


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <cstdlib>

#define USEMEAN 1


#include "TVessel.h"

typedef itk::Image<float , 3> ImageType3D;

/*typedef struct _tagVessel
{
    double q3[4];
    double q2[4];
    double q1[4];
    double a1;
    double a2;
    double a3;
    double f;
    double b;
    double mu[3];
    double e1;
    double e2;
    double last_dir[3];
    double last_pos[3];
    double L;
    double bndy[24];
    double type;
    double index;
    double MAD;
    double V_ID;
    double K1;
    double K2;
    double a;
    double k;
    double R1[3];
    double R2[3];
    double R3[3];

}TVessel;
*/


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



/*class SegInit: public itk::LightObject	{
public:
	typedef SegInit Self;
	typedef  itk::SmartPointer< Self >  Pointer;
	typedef  itk::SmartPointer< const Self >  ConstPointer;
	itkNewMacro(Self);

	*/
class SegInit	{
public:

bool fitSE (ImageType3D::Pointer, TVessel&, double , double , double);
void reverse(TVessel & seg);


protected:

	typedef struct _tagMatrix
	{
	    double r[3][3];
	}TRMatrix;

	TRMatrix gRT;


	typedef struct _tagVertex
	{
	    double p[3];
	}TVertex;

	typedef struct _tagFacet
	{
	    int vertex[3];
	    double centroid[3];
	    double normal[3];
	    double area;
	}TFacet;


	typedef struct _tagTDamp
	{
	    double dt2;
	    double dt;
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
	double gF[6000];

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
	void update_Quaternions2(TVessel & vessel,TFacet * pFacets, int NumFacets, double * pF, TRMatrix & R, double total_area, TDamp & damp);
	void update_e1(TVessel & vessel,TFacet * pFacets, int Num_facets,double * pF, TRMatrix & R, TDamp &damp);
	void update_mu(TVessel & vessel,TFacet * pFacets, int Num_Facets, double * pF,  TRMatrix & R, TDamp & damp);
	void update_scale( TVessel & vessel, TFacet * pFacets,int Num_Facets, double * pF, TRMatrix & R, TDamp & damp, double AS_RATIO );
	void axis_angle(double * q);
	void quat_mult(double * q3, double * q2, double * q1 );
	void update_axes( TVessel & vessel, int AT_END );
	double superquad3d_distq(double tx, double ty, double tz, TVessel & vessel );
	double rho_dot_asym(double x,double A,double B);
	int compare( const void * p1, const void * p2 );
	int compare_int( const void * p1, const void * p2 );

	double getMedian(std::vector<double> arr);
	double getMean(std::vector<double> arr);
	double getRegionMean (ImageType3D::Pointer im, double, double, double, double, double, double);
	//static void fill_array( double	*xr, double  *xi);

	int interp3(ImageType3D::Pointer im, double * pF, TFacet * pFacets, int Num_Facets, TVertex * pVertexC );
	bool interp3DefBG(ImageType3D::Pointer  im, double * pF, TFacet * pFacets, int Num_Facets, TVertex * pVertexC, double defBG );
	//double update_FB(ImageType3D::Pointer im, TVessel & vessel );
	double update_FBNew(ImageType3D::Pointer im, TVessel & vessel , int Initialize);
	double update_FBSimple(ImageType3D::Pointer im, TVessel & vessel );
	int handle_end( ImageType3D::Pointer im, double * pF, TVessel & vessel, TFacet * pFacets, int Num_Facets, TVertex * pUVertices, TVertex * pSVertices, TVertex * pCVertices, int Num_Vertices  );

public:
	SegInit()	{
		gNum_facets = 0;
		gNum_points = 0;
	}

	~SegInit()	{
	}

};

#endif /*SEGINIT_H_*/
