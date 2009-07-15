//////////////////////////////////////////////////////////////////////////////
//This module contains the full model estimation required for
//seed points.  This includes multiple intensity estimation steps
//and is typically called ~100 iterations.
//Also serves as the base class for the lightweight SegFit module
//used to refine model when tracing.

#include "SegInit.h"

#define VERBOSEFB 0
#define VERBOSEITER 0

//Flip direction using quaternion multiplication
void SegInit::reverse( TVessel & seg )
{
	TRMatrix R;

	rotation_quat(seg.q1,R);

	double qdir[4];
	qdir[0] = _PI;
	qdir[1] = R.r[0][1];
	qdir[2] = R.r[1][1];
	qdir[3] = R.r[2][1];
	
	axis_angle(qdir);
	double res[4];
	quat_mult(qdir,seg.q1,res);
	memcpy(seg.q1,res,sizeof(double)*4);

	rotation_quat(seg.q1,R);

	seg.R1[0] = R.r[0][0];
    seg.R1[1] = R.r[1][0];
    seg.R1[2] = R.r[2][0];

    seg.R2[0] = R.r[0][1];
    seg.R2[1] = R.r[1][1];
    seg.R2[2] = R.r[2][1];

    seg.R3[0] = R.r[0][2];
    seg.R3[1] = R.r[1][2];
    seg.R3[2] = R.r[2][2];


}

//get 3-D rotation matrix from quaternion
void SegInit::rotation_quat(double * q, TRMatrix & R)
{
    double w = q[0];
    double x = q[1];
    double y = q[2];
    double z = q[3];


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


void SegInit::transpose_matrix(TRMatrix & Result )
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

//rho is the inverse-generalized distance to superellipsoid boundary
double SegInit::calc_rho(double x, double y, double z, double e1, double e2 )
{
    double tmp1 = pow(fabs(x),2.0/e2);
    double tmp2 = pow(fabs(y),2.0/e2);
    double tmp3 = pow(tmp1+tmp2,e2/e1);
    double tmp4 = pow(fabs(z),2.0/e1);

    return pow(tmp3+tmp4,-(e1/2.0));

}

void SegInit::get_rotation(TVessel & vessel,TRMatrix & Result)
{

    rotation_quat(vessel.q1,Result);

}

//area of facet using cross-product
//also get the normal direction
double SegInit::get_area(TFacet & facet)
{
    int i1 = facet.vertex[0];
    int i2 = facet.vertex[1];
    int i3 = facet.vertex[2];

    double a1 = gVertexS[i1].p[0]-gVertexS[i2].p[0];
    double a2 = gVertexS[i1].p[1]-gVertexS[i2].p[1];
    double a3 = gVertexS[i1].p[2]-gVertexS[i2].p[2];


    double b1 = gVertexS[i1].p[0]-gVertexS[i3].p[0];
    double b2 = gVertexS[i1].p[1]-gVertexS[i3].p[1];
    double b3 = gVertexS[i1].p[2]-gVertexS[i3].p[2];

    double c1 = a2*b3-a3*b2;
    double c2 = -a1*b3+a3*b1;
    double c3 = a1*b2-a2*b1;
    double l = sqrt(c1*c1+c2*c2+c3*c3);

    facet.area = l/2.0;
    facet.normal[0] = c1/(l);
    facet.normal[1] = c2/(l);
    facet.normal[2] = c3/(l);


    facet.centroid[0] = gVertexS[i1].p[0]+gVertexS[i2].p[0]+gVertexS[i3].p[0];
    facet.centroid[0] /= 3.0;

    facet.centroid[1] = gVertexS[i1].p[1]+gVertexS[i2].p[1]+gVertexS[i3].p[1];
    facet.centroid[1] /= 3.0;

    facet.centroid[2] = gVertexS[i1].p[2]+gVertexS[i2].p[2]+gVertexS[i3].p[2];
    facet.centroid[2] /= 3.0;

    return facet.area;
}

//create the nodal mesh
void SegInit::generate_convex_hullq(TVessel & vessel,int NN)
{

    double step = 2.0*_PI/(double)(NN+1);
    gNum_facets = 0;
    gNum_points = 0;
    double mu[3];
    double a1,a2,a3;
    double e1,e2;

    a1 = vessel.a1;
    a2 = vessel.a2;
    a3 = vessel.a3;

    e1 = vessel.e1;
    e2 = vessel.e2;

    memcpy(mu,vessel.mu,sizeof(double)*3);

    double x,y,z,x1,y1,z1,xx1,yy1,zz1,rho;

    TRMatrix R;
    get_rotation(vessel,R);


    //int incount = 0;
    //int outcount = 0;

    //add first end point
    x = 0;
    y = 0;
    z = a3;

    gVertexU[gNum_points].p[0] = 0;
    gVertexU[gNum_points].p[1] = 0;
    gVertexU[gNum_points].p[2] = 1;

    gVertexS[gNum_points].p[0] = x;
    gVertexS[gNum_points].p[1] = y;
    gVertexS[gNum_points].p[2] = z;

    xx1 = R.r[0][0]*x + R.r[0][1]*y + R.r[0][2]*z+mu[0];
    yy1 = R.r[1][0]*x + R.r[1][1]*y + R.r[1][2]*z+mu[1];
    zz1 = R.r[2][0]*x + R.r[2][1]*y + R.r[2][2]*z+mu[2];

    gVertexC[gNum_points].p[0] = xx1;
    gVertexC[gNum_points].p[1] = yy1;
    gVertexC[gNum_points].p[2] = zz1;


    for(double n = step; n <= _PI-step; n+=step){
        for(double w = 0; w <= (double)NN*step; w+=step){
            gNum_points++;
            x1 = sin(n)*cos(w);
            y1 = sin(n)*sin(w);
            z1 = cos(n);

            x = a1*x1;
            y = a2*y1;
            z = a3*z1;

            gVertexU[gNum_points].p[0] = x1;
            gVertexU[gNum_points].p[1] = y1;
            gVertexU[gNum_points].p[2] = z1;

            rho = calc_rho(x1,y1,z1,e1,e2);
            x = rho*x;
            y = rho*y;
            z = rho*z;

            gVertexS[gNum_points].p[0] = x;
            gVertexS[gNum_points].p[1] = y;
            gVertexS[gNum_points].p[2] = z;

            xx1 = R.r[0][0]*x + R.r[0][1]*y + R.r[0][2]*z+mu[0];
            yy1 = R.r[1][0]*x + R.r[1][1]*y + R.r[1][2]*z+mu[1];
            zz1 = R.r[2][0]*x + R.r[2][1]*y + R.r[2][2]*z+mu[2];

            gVertexC[gNum_points].p[0] = xx1;
            gVertexC[gNum_points].p[1] = yy1;
            gVertexC[gNum_points].p[2] = zz1;
        }
    }


    //add last end point
    gNum_points++;

    gVertexU[gNum_points].p[0] = 0;
    gVertexU[gNum_points].p[1] = 0;
    gVertexU[gNum_points].p[2] = -1;


    x = 0;
    y = 0;
    z = -a3;

    gVertexS[gNum_points].p[0] = x;
    gVertexS[gNum_points].p[1] = y;
    gVertexS[gNum_points].p[2] = z;

    xx1 = R.r[0][0]*x + R.r[0][1]*y + R.r[0][2]*z+mu[0];
    yy1 = R.r[1][0]*x + R.r[1][1]*y + R.r[1][2]*z+mu[1];
    zz1 = R.r[2][0]*x + R.r[2][1]*y + R.r[2][2]*z+mu[2];

    gVertexC[gNum_points].p[0] = xx1;
    gVertexC[gNum_points].p[1] = yy1;
    gVertexC[gNum_points].p[2] = zz1;

    gNum_points++;
    //all points added

    int ncount = (int)((_PI-2*step)/step)+1;
    int wcount = NN+1;

    for( int u = 1; u < ncount; u++ ){
        for( int w = 1; w <= wcount; w++ ){
            gFacets[gNum_facets].vertex[0] = w + (u-1)*wcount;
            gFacets[gNum_facets].vertex[1] = w + (u)*wcount;
            //	double tmp = floor(((double)w+1)/wcount)*wcount;
            gFacets[gNum_facets].vertex[2] = (int)(((w+1)<=wcount?w+1:(w+1-floor(((double)w+1)/wcount)*wcount)) + (u-1)*wcount);
            get_area(gFacets[gNum_facets]);
            gNum_facets++;
        }
    }


    for( int u = 1; u < ncount; u++ ){
        for( int w = 1; w <= wcount; w++ ){
            gFacets[gNum_facets].vertex[0] = (int)(((w+1)<=wcount?w+1:w+1-floor(((double)w+1)/wcount)*wcount) + (u)*wcount);
            gFacets[gNum_facets].vertex[2] = w + (u)*wcount;
            gFacets[gNum_facets].vertex[1] = (int)(((w+1)<=wcount?w+1:w+1-floor(((double)w+1)/wcount)*wcount) + (u-1)*wcount);
            get_area(gFacets[gNum_facets]);
            gNum_facets++;
        }
    }

    for( int w = 1; w <= wcount; w++ ){
        gFacets[gNum_facets].vertex[0] = w;
        gFacets[gNum_facets].vertex[1] = (int)((w+1)<=wcount?w+1:w+1-floor(((double)w+1)/wcount)*wcount);
        gFacets[gNum_facets].vertex[2] = 0;
        get_area(gFacets[gNum_facets]);
        gNum_facets++;

    }

    for( int w = 1; w <= wcount; w++ ){
        gFacets[gNum_facets].vertex[0] = w+(ncount-1)*wcount;
        gFacets[gNum_facets].vertex[1] = ((w-1)<=0?(w-1)+wcount:(w-1))+(ncount-1)*wcount;
        gFacets[gNum_facets].vertex[2] = gNum_points-1;
        get_area(gFacets[gNum_facets]);
        gNum_facets++;
    }
}

//apply rotation matrix to a point
void SegInit::matrix_point(TVertex & point, TRMatrix & Result )
{

    TVertex pt;

    pt.p[0] = Result.r[0][0]*point.p[0] + Result.r[0][1]*point.p[1] + Result.r[0][2]*point.p[2];
    pt.p[1] = Result.r[1][0]*point.p[0] + Result.r[1][1]*point.p[1] + Result.r[1][2]*point.p[2];
    pt.p[2] = Result.r[2][0]*point.p[0] + Result.r[2][1]*point.p[1] + Result.r[2][2]*point.p[2];

    point.p[0] = pt.p[0];
    point.p[1] = pt.p[1];
    point.p[2] = pt.p[2];


}

//apply to individual components
void c_matrix_point(double x, double y,double z, cTRMatrix *Result )
{

    cTVertex pt,point;
    point.p[0] = x;
    point.p[1] = y;
    point.p[2] = z;

    pt.p[0] = Result->r[0][0]*point.p[0] + Result->r[0][1]*point.p[1] + Result->r[0][2]*point.p[2];
    pt.p[1] = Result->r[1][0]*point.p[0] + Result->r[1][1]*point.p[1] + Result->r[1][2]*point.p[2];
    pt.p[2] = Result->r[2][0]*point.p[0] + Result->r[2][1]*point.p[1] + Result->r[2][2]*point.p[2];

    x = pt.p[0];
    y = pt.p[1];
    z = pt.p[2];
}

//transform the canonical surface based on parameters in model
void SegInit::matrix_point(double & x, double & y,double & z, TRMatrix & Result )
{

    TVertex pt,point;
    point.p[0] = x;
    point.p[1] = y;
    point.p[2] = z;

    pt.p[0] = Result.r[0][0]*point.p[0] + Result.r[0][1]*point.p[1] + Result.r[0][2]*point.p[2];
    pt.p[1] = Result.r[1][0]*point.p[0] + Result.r[1][1]*point.p[1] + Result.r[1][2]*point.p[2];
    pt.p[2] = Result.r[2][0]*point.p[0] + Result.r[2][1]*point.p[1] + Result.r[2][2]*point.p[2];

    x = pt.p[0];
    y = pt.p[1];
    z = pt.p[2];
}

//transform the canonical surface based on parameters in model
void SegInit::get_surface2(TVessel & vessel, TFacet * pFacets, int NumFacets, TVertex * pUVertices, TVertex * pSVertices, TVertex * pCVertices, int NumVertices)
{

    TRMatrix R;
    get_rotation(vessel,R);
    double x1,y1,z1,x,y,z,xx1,yy1,zz1,rho;

    for( int i = 0; i < NumVertices; i++ ){

        x1 = pUVertices[i].p[0];
        y1 = pUVertices[i].p[1];
        z1 = pUVertices[i].p[2];

        x = vessel.a1*x1;
        y = vessel.a2*y1;
        z = vessel.a3*z1;

        rho = calc_rho(x1,y1,z1,vessel.e1,vessel.e2);
        x = rho*x;
        y = rho*y;
        z = rho*z;

        pSVertices[i].p[0] = x;
        pSVertices[i].p[1] = y;
        pSVertices[i].p[2] = z;

        xx1 = R.r[0][0]*x + R.r[0][1]*y + R.r[0][2]*z+vessel.mu[0];
        yy1 = R.r[1][0]*x + R.r[1][1]*y + R.r[1][2]*z+vessel.mu[1];
        zz1 = R.r[2][0]*x + R.r[2][1]*y + R.r[2][2]*z+vessel.mu[2];

        pCVertices[i].p[0] = xx1;
        pCVertices[i].p[1] = yy1;
        pCVertices[i].p[2] = zz1;
    }

    for(int i  = 0; i < NumFacets; i++ ){
        get_area(pFacets[i]);
        x = pFacets[i].normal[0];
        y = pFacets[i].normal[1];
        z = pFacets[i].normal[2];

        pFacets[i].normal[0] = R.r[0][0]*x + R.r[0][1]*y + R.r[0][2]*z;
        pFacets[i].normal[1] = R.r[1][0]*x + R.r[1][1]*y + R.r[1][2]*z;
        pFacets[i].normal[2] = R.r[2][0]*x + R.r[2][1]*y + R.r[2][2]*z;

    }
}

//interpret force from intensity by interpolation
int SegInit::interp3DefBG(ImageType3D::Pointer  im, double * pF, TFacet * pFacets, int Num_Facets, TVertex * pVertexC, double defBG )
{

	//uncomment to use low level implementation
//	return interp3( im, pF, pFacets, Num_Facets, pVertexC, defBG );

	//AMIT MADE THE CHANGES
	typedef itk::LinearInterpolateImageFunction<ImageType3D,double>  InterpolatorType;
	typedef InterpolatorType::PointType PointType3D;
	InterpolatorType::Pointer interp = InterpolatorType::New();
	interp->SetInputImage(im);

	int SUCCESS = 0;
	int i1,i2,i3;
	double tx,ty,tz;


	for( int fi = 0; fi < Num_Facets; fi++ ) {
		i1 = pFacets[fi].vertex[0];
		i2 = pFacets[fi].vertex[1];
		i3 = pFacets[fi].vertex[2];
		tx = (pVertexC[i1].p[0]+pVertexC[i2].p[0]+pVertexC[i3].p[0])/3.0;
		ty = (pVertexC[i1].p[1]+pVertexC[i2].p[1]+pVertexC[i3].p[1])/3.0;
		tz = (pVertexC[i1].p[2]+pVertexC[i2].p[2]+pVertexC[i3].p[2])/3.0;

		PointType3D pt;
		pt[0] = tx;
		pt[1] = ty;
		pt[2] = tz;
		if (interp->IsInsideBuffer( pt ))	{
			pF[fi] = static_cast<double>(interp->Evaluate( pt ));
			SUCCESS = 1;
		}
		else	{
			pF[fi] = defBG ;
		}
	}
	return SUCCESS;
}

//low level interpolation
int SegInit::interp3(ImageType3D::Pointer  im, double * pF, TFacet * pFacets, int Num_Facets, TVertex * pVertexC, double defBG )
{
	
    int f_x,f_y,f_z,i1,i2,i3;
    double tx,ty,tz,s,t,w;
	ImageType3D::SizeType im_dim = im->GetRequestedRegion().GetSize();

    for( int fi = 0; fi < Num_Facets; fi++ ) {

        i1 = pFacets[fi].vertex[0];
        i2 = pFacets[fi].vertex[1];
        i3 = pFacets[fi].vertex[2];

        tx = (pVertexC[i1].p[0]+pVertexC[i2].p[0]+pVertexC[i3].p[0])/3.0;
        ty = (pVertexC[i1].p[1]+pVertexC[i2].p[1]+pVertexC[i3].p[1])/3.0;
        tz = (pVertexC[i1].p[2]+pVertexC[i2].p[2]+pVertexC[i3].p[2])/3.0;

        f_x = (int)floor(tx);
        f_y = (int)floor(ty);
        f_z = (int)floor(tz);

        s = tx-(double)f_x;
        t = ty-(double)f_y;
        w = tz-(double)f_z;

        if( f_x < 1 || f_x >= im_dim[0] || f_y < 1 || f_y >= im_dim[1] || f_z < 1 || f_z >= im_dim[2] )
            pF[fi] = defBG;
		else{
            int ndx = f_x+(f_y-1)*im_dim[0]+(f_z-1)*im_dim[0]*im_dim[1];

			ImageType3D::IndexType ndx1 = {{(long int)f_x, (long int)f_y, (long int)f_z}};
			float pixval1 = im->GetPixel(ndx1);
			ImageType3D::IndexType ndx2 = {{(long int)f_x+1, (long int)f_y, (long int)f_z}};
			float pixval2 = im->GetPixel(ndx2);
			ImageType3D::IndexType ndx3 = {{(long int)f_x, (long int)f_y+1, (long int)f_z}};
			float pixval3 = im->GetPixel(ndx3);
			ImageType3D::IndexType ndx4 = {{(long int)f_x+1, (long int)f_y+1, (long int)f_z}};
			float pixval4 = im->GetPixel(ndx4);
			ImageType3D::IndexType ndx5 = {{(long int)f_x, (long int)f_y, (long int)f_z+1}};
			float pixval5 = im->GetPixel(ndx5);
			ImageType3D::IndexType ndx6 = {{(long int)f_x+1, (long int)f_y, (long int)f_z+1}};
			float pixval6 = im->GetPixel(ndx6);
			ImageType3D::IndexType ndx7 = {{(long int)f_x, (long int)f_y+1, (long int)f_z+1}};
			float pixval7 = im->GetPixel(ndx7);
			ImageType3D::IndexType ndx8 = {{(long int)f_x+1, (long int)f_y+1, (long int)f_z+1}};
			float pixval8 = im->GetPixel(ndx8);

            pF[fi] = ((pixval1*(1-t)+pixval2*t)*(1-s)+
            (pixval3*(1-t) + pixval4*t)*s)*(1-w)+
            ((pixval5*(1-t)+pixval6*t)*(1-s)+
            (pixval7*(1-t) + pixval8*t)*s)*(w);
        }
    }

    return 1;
}

//derivative wrt to shape parameter
double c_partial_rho_e1(double x, double y, double z, double e1, double e2 )
{
    return  pow(pow(pow(fabs(x),2.0/e2)+pow(fabs(y),2.0/e2),1.0*e2/e1)+pow(fabs(z),2.0/e1),-0.5*e1)*(-log(pow(pow(fabs(x),2.0/e2)+pow(fabs(y),2.0/e2),1.0*e2/e1)+pow(fabs(z),2.0/e1)+.000001)/2.0-e1*(-pow(pow(fabs(x),2.0/e2)+pow(fabs(y),2.0/e2),1.0*e2/e1)*e2/(e1*e1)*log(pow(fabs(x),2.0/e2)+pow(fabs(y),2.0/e2)+.0000001)-2.0*pow(fabs(z),2.0/e1)/(e1*e1)*log(fabs(z)+.000001))/(pow(pow(fabs(x),2.0/e2)+pow(fabs(y),2.0/e2),1.0*e2/e1)+pow(fabs(z),2.0/e1))/2.0);
}
double SegInit::partial_rho_e1(double x, double y, double z, double e1, double e2 )
{
    return  pow(pow(pow(fabs(x),2.0/e2)+pow(fabs(y),2.0/e2),1.0*e2/e1)+pow(fabs(z),2.0/e1),-0.5*e1)*(-log(pow(pow(fabs(x),2.0/e2)+pow(fabs(y),2.0/e2),1.0*e2/e1)+pow(fabs(z),2.0/e1)+.000001)/2.0-e1*(-pow(pow(fabs(x),2.0/e2)+pow(fabs(y),2.0/e2),1.0*e2/e1)*e2/(e1*e1)*log(pow(fabs(x),2.0/e2)+pow(fabs(y),2.0/e2)+.0000001)-2.0*pow(fabs(z),2.0/e1)/(e1*e1)*log(fabs(z)+.000001))/(pow(pow(fabs(x),2.0/e2)+pow(fabs(y),2.0/e2),1.0*e2/e1)+pow(fabs(z),2.0/e1))/2.0);
}



//update quaternion parameters using coordinate descent
void SegInit::update_Quaternions2(TVessel & vessel,TFacet * pFacets, int NumFacets, double * pF, TRMatrix & R, double total_area, TDamp & damp)
{

	double s,v1,v2,v3;
    double * sign_S = damp.sign_S;

    s = vessel.q1[0];
    v1 = vessel.q1[1];
    v2 = vessel.q1[2];
    v3 = vessel.q1[3];

    double R1[3],R2[3],R3[3],C1[3],C2[3],C3[3],C4[3];
    double p1,p2,p3;
    double s1=0,s2=0,s3=0,s4=0;

    for( int i = 0; i < NumFacets; i++ ){
        p1 = pFacets[i].centroid[0];
        p2 = pFacets[i].centroid[1];
        p3 = pFacets[i].centroid[2];

        R1[0] = R.r[0][0]*0	+ R.r[0][1]*p3 + R.r[0][2]*-p2;
        R1[1] = R.r[0][0]*-p3	+ R.r[0][1]*0 + R.r[0][2]*p1;
        R1[2] = R.r[0][0]*p2	+ R.r[0][1]*-p1 + R.r[0][2]*0;

        R2[0] = R.r[1][0]*0	+ R.r[1][1]*p3 + R.r[1][2]*-p2;
        R2[1] = R.r[1][0]*-p3	+ R.r[1][1]*0 + R.r[1][2]*p1;
        R2[2] = R.r[1][0]*p2	+ R.r[1][1]*-p1 + R.r[1][2]*0;

        R3[0] = R.r[2][0]*0	+ R.r[2][1]*p3 + R.r[2][2]*-p2;
        R3[1] = R.r[2][0]*-p3	+ R.r[2][1]*0 + R.r[2][2]*p1;
        R3[2] = R.r[2][0]*p2	+ R.r[2][1]*-p1 + R.r[2][2]*0;

        C1[0] = 2*(R1[0]*-v1 + R1[1]*-v2 + R1[2]*-v3);
        C1[1] = 2*(R2[0]*-v1 + R2[1]*-v2 + R2[2]*-v3);
        C1[2] = 2*(R3[0]*-v1 + R3[1]*-v2 + R3[2]*-v3);

        C2[0] = 2*(R1[0]*s + R1[1]*-v3 + R1[2]*v2);
        C2[1] = 2*(R2[0]*s + R2[1]*-v3 + R2[2]*v2);
        C2[2] = 2*(R3[0]*s + R3[1]*-v3 + R3[2]*v2);

        C3[0] = 2*(R1[0]*v3 + R1[1]*s + R1[2]*-v1);
        C3[1] = 2*(R2[0]*v3 + R2[1]*s + R2[2]*-v1);
        C3[2] = 2*(R3[0]*v3 + R3[1]*s + R3[2]*-v1);

        C4[0] = 2*(R1[0]*-v2 + R1[1]*v1 + R1[2]*s);
        C4[1] = 2*(R2[0]*-v2 + R2[1]*v1 + R2[2]*s);
        C4[2] = 2*(R3[0]*-v2 + R3[1]*v1 + R3[2]*s);

        C1[0] = C1[0]*pFacets[i].normal[0] + C1[1]*pFacets[i].normal[1] + C1[2]*pFacets[i].normal[2];
        C2[0] = C2[0]*pFacets[i].normal[0] + C2[1]*pFacets[i].normal[1] + C2[2]*pFacets[i].normal[2];
        C3[0] = C3[0]*pFacets[i].normal[0] + C3[1]*pFacets[i].normal[1] + C3[2]*pFacets[i].normal[2];
        C4[0] = C4[0]*pFacets[i].normal[0] + C4[1]*pFacets[i].normal[1] + C4[2]*pFacets[i].normal[2];

        s1 += C1[0]*pF[i]*damp.dt2/total_area;
        s2 += C2[0]*pF[i]*damp.dt2/total_area;
        s3 += C3[0]*pF[i]*damp.dt2/total_area;
        s4 += C4[0]*pF[i]*damp.dt2/total_area;
    }

    vessel.q1[0] += s1;
    vessel.q1[1] += s2;
    vessel.q1[2] += s3;
    vessel.q1[3] += s4;

    if(s1*damp.sign_S[0] <= 0 || s2*damp.sign_S[1] <= 0 || s3*damp.sign_S[2] <= 0 || s4*damp.sign_S[3] <= 0)
        damp.dt2*=.8;

    sign_S[0] = s1;
    sign_S[1] = s2;
    sign_S[2] = s3;
    sign_S[3] = s4;

    double l = vessel.q1[0]*vessel.q1[0]+vessel.q1[1]*vessel.q1[1]+vessel.q1[2]*vessel.q1[2]+vessel.q1[3]*vessel.q1[3];
    l = sqrt(l);

    vessel.q1[0] /= l;
    vessel.q1[1] /= l;
    vessel.q1[2] /= l;
    vessel.q1[3] /= l;

}

//update shape parameter using coordinate descent
void c_update_e1(void *cpp_vessel,void *cpp_pFacets, int Num_facets,double * pF, void *cpp_R, void *cpp_damp, int gNum_facets, void *cpp_gVertexU)
{

    cTVessel *vessel = (cTVessel*)cpp_vessel;
    cTFacet *pFacets = (cTFacet*)cpp_pFacets;
    cTRMatrix *R = (cTRMatrix*)cpp_R;
    cTDamp *damp = (cTDamp*)cpp_damp;
    cTVertex *gVertexU = (cTVertex*)cpp_gVertexU;
    
    double de1;
    double x,y,z;
    double ss[3];
    double sde1 = 0;

    //TRMatrix R;
    //get_rotation(vessel,R);

    for(int i = 0; i < gNum_facets; i++ ){
        x = gVertexU[pFacets[i].vertex[0]].p[0];
        y = gVertexU[pFacets[i].vertex[0]].p[1];
        z = gVertexU[pFacets[i].vertex[0]].p[2];

        de1 = c_partial_rho_e1(x,y,z,vessel->e1,1.0);
        ss[0] = vessel->a1*x*de1;
        ss[1] = vessel->a2*y*de1;
        ss[2] = vessel->a3*z*de1;


        x = gVertexU[pFacets[i].vertex[1]].p[0];
        y = gVertexU[pFacets[i].vertex[1]].p[1];
        z = gVertexU[pFacets[i].vertex[1]].p[2];

        de1 = c_partial_rho_e1(x,y,z,vessel->e1,1.0);
        ss[0] += vessel->a1*x*de1;
        ss[1] += vessel->a2*y*de1;
        ss[2] += vessel->a3*z*de1;

        x = gVertexU[pFacets[i].vertex[2]].p[0];
        y = gVertexU[pFacets[i].vertex[2]].p[1];
        z = gVertexU[pFacets[i].vertex[2]].p[2];

        de1 = c_partial_rho_e1(x,y,z,vessel->e1,1.0);
        ss[0] += vessel->a1*x*de1;
        ss[1] += vessel->a2*y*de1;
        ss[2] += vessel->a3*z*de1;

        ss[0] /= 3.0;
        ss[1] /= 3.0;
        ss[2] /= 3.0;

        c_matrix_point(ss[0],ss[1],ss[2], R );

        sde1 -= damp->dt/2.0*pF[i]*(ss[0]*pFacets[i].normal[0]+ss[1]*pFacets[i].normal[1]+ss[2]*pFacets[i].normal[2]);

    }

    if(sde1+vessel->e1 > 1.0)
        sde1 = 0;

    vessel->e1 += sde1;
    vessel->e1 = vnl_math_max(.25,vessel->e1);
    vessel->e1 = vnl_math_max(1.0,vessel->e1);
}

void SegInit::update_e1(TVessel & vessel,TFacet * pFacets, int Num_facets,double * pF, TRMatrix & R, TDamp &damp)
{
    double de1;
    double x,y,z;
    double ss[3];
    double sde1 = 0;

    //TRMatrix R;
    //get_rotation(vessel,R);

    for(int i = 0; i < gNum_facets; i++ ){
        x = gVertexU[pFacets[i].vertex[0]].p[0];
        y = gVertexU[pFacets[i].vertex[0]].p[1];
        z = gVertexU[pFacets[i].vertex[0]].p[2];

        de1 = partial_rho_e1(x,y,z,vessel.e1,1.0);
        ss[0] = vessel.a1*x*de1;
        ss[1] = vessel.a2*y*de1;
        ss[2] = vessel.a3*z*de1;


        x = gVertexU[pFacets[i].vertex[1]].p[0];
        y = gVertexU[pFacets[i].vertex[1]].p[1];
        z = gVertexU[pFacets[i].vertex[1]].p[2];

        de1 = partial_rho_e1(x,y,z,vessel.e1,1.0);
        ss[0] += vessel.a1*x*de1;
        ss[1] += vessel.a2*y*de1;
        ss[2] += vessel.a3*z*de1;

        x = gVertexU[pFacets[i].vertex[2]].p[0];
        y = gVertexU[pFacets[i].vertex[2]].p[1];
        z = gVertexU[pFacets[i].vertex[2]].p[2];

        de1 = partial_rho_e1(x,y,z,vessel.e1,1.0);
        ss[0] += vessel.a1*x*de1;
        ss[1] += vessel.a2*y*de1;
        ss[2] += vessel.a3*z*de1;

        ss[0] /= 3.0;
        ss[1] /= 3.0;
        ss[2] /= 3.0;

        matrix_point(ss[0],ss[1],ss[2], R );

        sde1 -= damp.dt/2.0*pF[i]*(ss[0]*pFacets[i].normal[0]+ss[1]*pFacets[i].normal[1]+ss[2]*pFacets[i].normal[2]);

    }

    if(sde1+vessel.e1 > 1.0)
        sde1 = 0;

    vessel.e1 += sde1;
    vessel.e1 = vnl_math_max(.25,vessel.e1);
    vessel.e1 = vnl_math_max(1.0,vessel.e1);
}


//update model center point using coordinate descent
void	SegInit::update_mu(TVessel & vessel,TFacet * pFacets, int Num_Facets, double * pF,  TRMatrix & R, TDamp & damp)
{
    double du1=0,du2=0,du3=0;
    double du1b=0,du2b=0,du3b=0;
    double * dt_u = damp.dt_u;
    //TRMatrix R;

    //get_rotation(vessel,R);
    vessel.R3[0] = R.r[0][2];
    vessel.R3[1] = R.r[1][2];
    vessel.R3[2] = R.r[2][2];

    for(int i = 0; i < Num_Facets; i++){
        du1 -= dt_u[0]*pFacets[i].normal[0]*pF[i];
        du2 -= dt_u[1]*pFacets[i].normal[1]*pF[i];
        du3 -= dt_u[2]*pFacets[i].normal[2]*pF[i];
    }

    du1b = du1-(vessel.R3[0]*du1+vessel.R3[1]*du2+vessel.R3[2]*du3)*vessel.R3[0];
    du2b = du2-(vessel.R3[0]*du1+vessel.R3[1]*du2+vessel.R3[2]*du3)*vessel.R3[1];
    du3b = du3-(vessel.R3[0]*du1+vessel.R3[1]*du2+vessel.R3[2]*du3)*vessel.R3[2];

    vessel.mu[0] += du1b;
    vessel.mu[1] += du2b;
    vessel.mu[2] += du3b;

    if(du1b*damp.sign_U[0] <= 0 )
        damp.dt_u[0]*=.8;

    if(du2b*damp.sign_U[1] <= 0 )
        damp.dt_u[1]*=.8;

    if(du3b*damp.sign_U[2] <= 0 )
        damp.dt_u[2]*=.8;

    damp.sign_U[0] = du1b;
    damp.sign_U[1] = du2b;
    damp.sign_U[2] = du3b;
}


//update scale parameters using coordinate descent
void SegInit::update_scale( TVessel & vessel, TFacet * pFacets,int Num_Facets, double * pF, TRMatrix & R, TDamp & damp, double AS_RATIO )
{
    //TRMatrix R;
    //get_rotation(vessel,R);

    double a1,a2,a3;
    double da1=0,da2=0,da3=0;
    double * dt_a = damp.dt_a;

    a1 = vessel.a1;
    a2 = vessel.a2;
    a3 = vessel.a3;

    for( int i = 0; i < Num_Facets; i++ ){
        da1 += pF[i]*(R.r[0][0]*pFacets[i].centroid[0]/a1*pFacets[i].normal[0]+R.r[1][0]*pFacets[i].centroid[0]/a1*pFacets[i].normal[1]+R.r[2][0]*pFacets[i].centroid[0]/a1*pFacets[i].normal[2]);
        da2 += pF[i]*(R.r[0][1]*pFacets[i].centroid[1]/a2*pFacets[i].normal[0]+R.r[1][1]*pFacets[i].centroid[1]/a2*pFacets[i].normal[1]+R.r[2][1]*pFacets[i].centroid[1]/a2*pFacets[i].normal[2]);
        da3 += pF[i]*(R.r[0][2]*pFacets[i].centroid[2]/a3*pFacets[i].normal[0]+R.r[1][2]*pFacets[i].centroid[2]/a3*pFacets[i].normal[1]+R.r[2][2]*pFacets[i].centroid[2]/a3*pFacets[i].normal[2]);
    }

    da1 *= -dt_a[0];
    da2 *= -dt_a[1];
    da3 *= -dt_a[2];

//    vessel.a1 += da1;
 //   vessel.a2 += da2;

    if ( vessel.a1 < vessel.a2 ){
        vessel.a1 += da1;
        if( a2 + da2 < 4*vessel.a1 && (da2 > 0) )
            vessel.a2 += da2;
        else if( da2 < 0 )
            vessel.a2 += da2;
    }
    else{
        vessel.a2 += da2;
        if( a1 + da1 < 4*vessel.a2 && (da1 > 0) )
            vessel.a1 += da1;
        else if( da1 < 0 )
            vessel.a1 += da1;
    }

    if( (a3 + da3 < AS_RATIO*vnl_math_max(vessel.a2,vessel.a1)) && (da3 > 0) )
        vessel.a3 += da3;
    else if( da3 < 0 )
        vessel.a3 += da3;


    if(da1*damp.sign_A[0] <= 0 )
        damp.dt_a[0]*=.8;

    if(da2*damp.sign_A[1] <= 0 )
        damp.dt_a[1]*=.8;

    if(da3*damp.sign_A[2] <= 0 )
        damp.dt_a[2]*=.8;

    damp.sign_A[0] = da1;
    damp.sign_A[1] = da2;
    damp.sign_A[2] = da3;

    vessel.a1 = fabs(vessel.a1);
    vessel.a2 = fabs(vessel.a2);
    vessel.a3 = fabs(vessel.a3);

    vessel.a1 = vnl_math_max(1.5,vessel.a1);
    vessel.a2 = vnl_math_max(1.5,vessel.a2);
    vessel.a3 = vnl_math_max(3.0,vessel.a3);

    vessel.a3 = vnl_math_min(vessel.a3,AS_RATIO*vnl_math_max(vessel.a2,vessel.a1));

}


//quaternion representation
  void SegInit::axis_angle(double * q)
{
    double a = q[0];

    q[0] = cos(a/2.0);
    q[1] = sin(a/2.0)*q[1];
    q[2] = sin(a/2.0)*q[2];
    q[3] = sin(a/2.0)*q[3];

}

  //quaternion arithmetic
  void SegInit::quat_mult(double * q3, double * q2, double * q1 )
{
    double q[4];


    q[0] = q3[0]*q2[0] - q3[1]*q2[1] - q3[2]*q2[2] - q3[3]*q2[3];
    q[1] = q3[0]*q2[1] + q3[1]*q2[0] + q3[2]*q2[3] - q3[3]*q2[2];
    q[2] = q3[0]*q2[2] - q3[1]*q2[3] + q3[2]*q2[0] + q3[3]*q2[1];
    q[3] = q3[0]*q2[3] + q3[1]*q2[2] - q3[2]*q2[1] + q3[3]*q2[0];

    memcpy(q1,q,sizeof(double)*4);


}

  //if the model's semi major/minor axis change, a change of direction is implied
  void SegInit::update_axes( TVessel & vessel, int AT_END )
{
    TRMatrix R;

    if( vessel.a1 > 1.1*vessel.a3 && vessel.a1 > vessel.a2 && !AT_END ){
        get_rotation(vessel,R);
        double q[4];

        q[0] = _PI/2;
        q[1] = R.r[0][1];
        q[2] = R.r[1][1];
        q[3] = R.r[2][1];
        axis_angle(q);

        quat_mult(q,vessel.q1,vessel.q1);
        double tmp = vessel.a1;
        vessel.a1 = vessel.a3;
        vessel.a3 = tmp;
    }else if ( vessel.a2 > 1.1*vessel.a3 && vessel.a2 > vessel.a1 && !AT_END ){
        get_rotation(vessel,R);
        double q[4];

        q[0] = _PI/2;
        q[1] = R.r[0][0];
        q[2] = R.r[1][0];
        q[3] = R.r[2][0];
        axis_angle(q);

        quat_mult(q,vessel.q1,vessel.q1);
        double tmp = vessel.a2;
        vessel.a2 = vessel.a3;
        vessel.a3 = tmp;
    }
}

  //implicit superellipsoid distance function
double SegInit::superquad3d_distq(double tx, double ty, double tz, TVessel & vessel )
{
    double mu[3];
    double s1,s2,s3;
    double e1,e2;

    s1 = vessel.a1;
    s2 = vessel.a2;
    s3 = vessel.a3;

    e1 = vessel.e1;
    e2 = vessel.e2;

    memcpy(mu,vessel.mu,sizeof(double)*3);

    TRMatrix Result;

    get_rotation(vessel,Result);
    transpose_matrix(Result);

    double x = Result.r[0][0]*(tx-mu[0])+Result.r[0][1]*(ty-mu[1])+Result.r[0][2]*(tz-mu[2]);
    double y = Result.r[1][0]*(tx-mu[0])+Result.r[1][1]*(ty-mu[1])+Result.r[1][2]*(tz-mu[2]);
    double z = Result.r[2][0]*(tx-mu[0])+Result.r[2][1]*(ty-mu[1])+Result.r[2][2]*(tz-mu[2]);

    return pow(pow(fabs(x)/s1,2.0/e2)+pow(fabs(y)/s2,2.0/e2),e2/e1)+pow(fabs(z)/s3,2.0/e1);

}


  double SegInit::rho_dot_asym(double x,double A,double B)
{
    double val = 0;

    if (x < 0 )
        val=x/(1.0+pow(x/A,2.0));
    else
        val=x/(1.0+pow(x/B,2.0));

    return val;
}

  //sort function
int SegInit::compare( const void * p1, const void * p2 )
{
    double * P1 = (double*)p1;
    double * P2 = (double*)p2;

    if(*P1 > * P2)
        return -1;
    if(*P1 < * P2)
        return 1;

    return 0;
}

int SegInit::compare_int( const void * p1, const void * p2 )
{
    int * P1 = (int*)p1;
    int * P2 = (int*)p2;

    if(*P1 > * P2)
        return 1;
    if(*P1 < * P2)
        return -1;

    return 0;
}

//non-robust intensity estimation
double SegInit::update_FBSimple(ImageType3D::Pointer im, TVessel & vessel )
{

	//This is just run once to initialise the super ellipse fitting
	double A = 3.0*vnl_math_max(vessel.a1,vessel.a2);
	A = vnl_math_min(A,50.0);

	std::vector<double> f_tmp; 	f_tmp.reserve(MAX_SIZE);
	std::vector<double> b_tmp; 	b_tmp.reserve(MAX_SIZE);


	ImageType3D::SizeType im_dim = im->GetRequestedRegion().GetSize();

	TRMatrix Result;

	get_rotation(vessel,Result);
	transpose_matrix(Result);

	double mu[3];
	memcpy(mu,vessel.mu,sizeof(double)*3);
	double s1 = vessel.a1;
	double s2 = vessel.a2;
	double s3 = vessel.a3;
	double e1 = vessel.e1;
	double e2 = 1.0;
	double D;


	TRMatrix lim;
	lim.r[0][0] = (int)(vessel.mu[0]+.5)-(int)(A+.5);
	lim.r[0][1] = (int)(vessel.mu[0]+.5)+(int)(A+.5);

	lim.r[1][0] = (int)(vessel.mu[1]+.5)-(int)(A+.5);
	lim.r[1][1] = (int)(vessel.mu[1]+.5)+(int)(A+.5);

	lim.r[2][0] = (int)(vessel.mu[2]+.5)-(int)(A/2.0+.5);
	lim.r[2][1] = (int)(vessel.mu[2]+.5)+(int)(A/2.0+.5);



	double x,y,z;
	double i,j,k;
	double sampling = 1;

	if( A > 20 )
		sampling = 2;


	for( k = lim.r[2][0]; k <= lim.r[2][1]; k+=sampling){
		for( i = lim.r[0][0]; i <= lim.r[0][1]; i+=sampling){
			for( j = lim.r[1][0]; j <= lim.r[1][1]; j+=sampling){
				if( i >= 1 && i < im_dim[0] && j >= 1 && j < im_dim[1] && k >= 1 && k < im_dim[2] ){
					x = Result.r[0][0]*(i-mu[0])+Result.r[0][1]*(j-mu[1])+Result.r[0][2]*(k-mu[2]);
					y = Result.r[1][0]*(i-mu[0])+Result.r[1][1]*(j-mu[1])+Result.r[1][2]*(k-mu[2]);
					z = Result.r[2][0]*(i-mu[0])+Result.r[2][1]*(j-mu[1])+Result.r[2][2]*(k-mu[2]);

					D = pow(pow(fabs(x)/s1,2.0/e2)+pow(fabs(y)/s2,2./e2),e2/e1)+pow(fabs(z)/s3,2.0/e1);

					ImageType3D::IndexType ndx = {{(long int)i, (long int)j, (long int)k}};
					float pixelVal = im->GetPixel(ndx);

					if( D <= 1.0)	{
						f_tmp.push_back(pixelVal);
					}

					else {
						b_tmp.push_back(pixelVal);
					}
				}
			}
		}
	}

	vessel.f = getMedian(f_tmp);
	vessel.b = getMedian(b_tmp);

	//vessel.f = vnl_math_min(vessel.b - 1.0,vessel.f);		//f always < b

	std::vector<double> fdev;
	fdev.reserve(f_tmp.size());
	std::vector<double>::iterator it1;
	for (it1 = f_tmp.begin(); it1<f_tmp.end(); it1++) {
		fdev.push_back(vnl_math_abs(vessel.f - (*it1)));
	}
	double MAD1 = getMedian(fdev);

	std::vector<double> bdev;
	bdev.reserve(b_tmp.size());
	std::vector<double>::iterator it2;
	for (it2 = b_tmp.begin(); it2<b_tmp.end(); it2++) {
		bdev.push_back(vnl_math_abs(vessel.b - (*it2)));
	}
	double MAD2 = getMedian(bdev);

	double MAD = vnl_math_max(MAD1,MAD2);
	vessel.MAD = MAD;


	std::vector<double>::iterator it3;
	double L = 0.0;
	for (it3 = f_tmp.begin(); it3 < f_tmp.end(); it3++) {
		L += -fabs((*it3)-vessel.f)+fabs((*it3)-vessel.b);
	}
	vessel.L = L/(static_cast<double>(f_tmp.size())+.0001f);

	return vessel.L;
}


//Amit's intensity estimation 
double SegInit::update_FBNew(ImageType3D::Pointer im, TVessel & vessel , int Initialize)
{

///AMIT REWROTE THIS MODULE ::
	double A = 3.0*vnl_math_max(vessel.a1,vessel.a2);
	A = vnl_math_min(A,50.0);

	std::vector<double> f_tmp; 	f_tmp.reserve(MAX_SIZE);
	std::vector<double> b_tmp; 	b_tmp.reserve(MAX_SIZE);
	std::vector<double> inside; inside.reserve(MAX_SIZE);


	ImageType3D::SizeType im_dim = im->GetRequestedRegion().GetSize();
	TRMatrix Result;

	get_rotation(vessel,Result);
	transpose_matrix(Result);

	double mu[3];
	memcpy(mu,vessel.mu,sizeof(double)*3);
	double s1 = vessel.a1;
	double s2 = vessel.a2;
	double s3 = vessel.a3;
	double e1 = vessel.e1;
	double e2 = 1.0;
	double D;


	TRMatrix lim;
	lim.r[0][0] = (int)(vessel.mu[0]+.5)-(int)(A+.5);
	lim.r[0][1] = (int)(vessel.mu[0]+.5)+(int)(A+.5);

	lim.r[1][0] = (int)(vessel.mu[1]+.5)-(int)(A+.5);
	lim.r[1][1] = (int)(vessel.mu[1]+.5)+(int)(A+.5);

	lim.r[2][0] = (int)(vessel.mu[2]+.5)-(int)(A/2.0+.5);
	lim.r[2][1] = (int)(vessel.mu[2]+.5)+(int)(A/2.0+.5);

	double x,y,z;
	double i,j,k;
	double sampling = 1;

	if( A > 20 )
		sampling = 2;


	double cutoff;
	if (Initialize==1)	{
		cutoff = getRegionMean(im, lim);
	}
	else	{
		cutoff = 0.5* (vessel.f+vessel.b);
	}

	
	for( k = lim.r[2][0]; k <= lim.r[2][1]; k+=sampling){
		for( i = lim.r[0][0]; i <= lim.r[0][1]; i+=sampling){
			for( j = lim.r[1][0]; j <= lim.r[1][1]; j+=sampling){
				if( i >= 1 && i < im_dim[0] && j >= 1 && j < im_dim[1] && k >= 1 && k < im_dim[2] ){
					x = Result.r[0][0]*(i-mu[0])+Result.r[0][1]*(j-mu[1])+Result.r[0][2]*(k-mu[2]);
					y = Result.r[1][0]*(i-mu[0])+Result.r[1][1]*(j-mu[1])+Result.r[1][2]*(k-mu[2]);
					z = Result.r[2][0]*(i-mu[0])+Result.r[2][1]*(j-mu[1])+Result.r[2][2]*(k-mu[2]);

					D = pow(pow(fabs(x)/s1,2.0/e2)+pow(fabs(y)/s2,2./e2),e2/e1)+pow(fabs(z)/s3,2.0/e1);


					ImageType3D::IndexType ndx = {{(long int)i, (long int)j, (long int)k}};
					float pixelVal = im->GetPixel(ndx);
					

					if( D <= 1.0 && pixelVal <= cutoff ){
						f_tmp.push_back(pixelVal); 
					}
				    else if( D > 1.0 && pixelVal > cutoff   ){
						b_tmp.push_back(pixelVal); 
					}

					if(D <= 1.0){
						inside.push_back(pixelVal);
					}
				}
			}
		}
	}

	if (f_tmp.size()<10)	{
		vessel.f = getMean(inside);
	}
	else {
		vessel.f = getMedian(f_tmp);
	}

	vessel.b = getMedian(b_tmp);

	vessel.f = vnl_math_min(vessel.b,vessel.f-1);

	std::vector<double> fdev;
	fdev.reserve(f_tmp.size());
	std::vector<double>::iterator it1;
	for (it1 = f_tmp.begin(); it1<f_tmp.end(); it1++) {
		fdev.push_back(vnl_math_abs(vessel.f- (*it1)));
	}
	double MAD1 = getMedian(fdev);

	std::vector<double> bdev;
	bdev.reserve(b_tmp.size());
	std::vector<double>::iterator it2;
	for (it2 = b_tmp.begin(); it2<b_tmp.end(); it2++) {
		bdev.push_back(vnl_math_abs(vessel.b- (*it2)));
	}
	double MAD2 = getMedian(bdev);

	double MAD = vnl_math_max(MAD1,MAD2);
	vessel.MAD = MAD;


	std::vector<double>::iterator it3;
	double L = 0.0;
	for (it3 = inside.begin(); it3 < inside.end(); it3++) {
		L += -fabs((*it3)-vessel.f)+fabs((*it3)-vessel.b);
	}
	vessel.L = L/(static_cast<double>(inside.size())+0.0001f);

	if (VERBOSEFB )	{
		std::cout << "Update FB @ [" << vessel.mu[0] << ", " << vessel.mu[1] << ", " << vessel.mu[2] << "] width= "<< A << std::endl;
		std::cout << "\t# inside=" << inside.size() << ",  # Fg=" << f_tmp.size() << ",  # Bg=" << b_tmp.size() << std::endl;
		std::cout << "\tFest=" << vessel.f << ", Best=" << vessel.b << ",  MAD=" << vessel.MAD << ",  L=" << vessel.L << std::endl;
	}

	return vessel.L;
}


//The correct, robust intensity and likelihood estimation
double SegInit::update_FB(ImageType3D::Pointer im, TVessel & vessel )
{

	///AMIT MADE CHANGES IN IMAGE TYPE AND ACCESS PROCEDURES.. ALEX'S CODE COMMENTED
	double A = 3.0*vnl_math_max(vessel.a1,vessel.a2);
	A = vnl_math_min(A,50.0);

	std::vector<double> f_tmp; 	f_tmp.reserve(MAX_SIZE);
	std::vector<double> b_tmp; 	b_tmp.reserve(MAX_SIZE);
	std::vector<double> inside; inside.reserve(MAX_SIZE);


	ImageType3D::SizeType im_dim = im->GetRequestedRegion().GetSize();
	int f_count = 0;
	int b_count = 0;

	TRMatrix Result;

	get_rotation(vessel,Result);
	transpose_matrix(Result);

	double mu[3];
	memcpy(mu,vessel.mu,sizeof(double)*3);
	double s1 = vessel.a1;
	double s2 = vessel.a2;
	double s3 = vessel.a3;
	double e1 = vessel.e1;
	double e2 = 1.0;
	double D;


	TRMatrix lim;
	lim.r[0][0] = (int)(vessel.mu[0]+.5)-(int)(A+.5);
	lim.r[0][1] = (int)(vessel.mu[0]+.5)+(int)(A+.5);

	lim.r[1][0] = (int)(vessel.mu[1]+.5)-(int)(A+.5);
	lim.r[1][1] = (int)(vessel.mu[1]+.5)+(int)(A+.5);

	lim.r[2][0] = (int)(vessel.mu[2]+.5)-(int)(A/2.0+.5);
	lim.r[2][1] = (int)(vessel.mu[2]+.5)+(int)(A/2.0+.5);



	double x,y,z;
	double i,j,k;
	int loop_count = 0;
	int loop_count2  =0;

	double sampling = 1;

	if( A > 20 )
		sampling = 2;

	double cutoff = vessel.f+vessel.b;
	cutoff /= 2;

	//vessel.PrintSelf();

	for( k = lim.r[2][0]; k <= lim.r[2][1]; k+=sampling){
		for( i = lim.r[0][0]; i <= lim.r[0][1]; i+=sampling){
			for( j = lim.r[1][0]; j <= lim.r[1][1]; j+=sampling){
				if( i >= 1 && i < im_dim[0] && j >= 1 && j < im_dim[1] && k >= 1 && k < im_dim[2] ){
					x = Result.r[0][0]*(i-mu[0])+Result.r[0][1]*(j-mu[1])+Result.r[0][2]*(k-mu[2]);
					y = Result.r[1][0]*(i-mu[0])+Result.r[1][1]*(j-mu[1])+Result.r[1][2]*(k-mu[2]);
					z = Result.r[2][0]*(i-mu[0])+Result.r[2][1]*(j-mu[1])+Result.r[2][2]*(k-mu[2]);

					D = pow(pow(fabs(x)/s1,2.0/e2)+pow(fabs(y)/s2,2./e2),e2/e1)+pow(fabs(z)/s3,2.0/e1);

					ImageType3D::IndexType ndx = {{(long int)i, (long int)j, (long int)k}};
					float pixelVal = im->GetPixel(ndx);


					if( D <= 1.0 && pixelVal <= cutoff ){
						f_tmp.push_back(pixelVal); //f_tmp[f_count]=pixelVal;
						f_count++;
					}
				    else if( D > 1.0 && pixelVal > cutoff && pixelVal < 254  ){
						b_tmp.push_back(pixelVal); //b_tmp[b_count]=pixelVal;
						b_count++;
					}

					if(D <= 1.0){
						inside.push_back(pixelVal);  //inside[loop_count]=pixelVal;
						loop_count++;
					}
				}
			}
		}
	}

	vessel.f = getMedian(f_tmp);
	vessel.b = getMedian(b_tmp);

	vessel.f = vnl_math_min(vessel.b,vessel.f);

	std::vector<double> fdev;
	fdev.reserve(f_tmp.size());
	std::vector<double>::iterator it1;
	for (it1 = f_tmp.begin(); it1<f_tmp.end(); it1++) {
		fdev.push_back(vnl_math_abs(vessel.f- (*it1)));
	}
	double MAD1 = getMedian(fdev);

	std::vector<double> bdev;
	bdev.reserve(b_tmp.size());
	std::vector<double>::iterator it2;
	for (it2 = b_tmp.begin(); it2<b_tmp.end(); it2++) {
		bdev.push_back(vnl_math_abs(vessel.b- (*it2)));
	}
	double MAD2 = getMedian(bdev);

	double MAD = vnl_math_max(MAD1,MAD2);
	vessel.MAD = MAD;


	std::vector<double>::iterator it3;
	double L = 0.0;
	for (it3 = inside.begin(); it3 < inside.end(); it3++) {
		L += -fabs((*it3)-vessel.f)+fabs((*it3)-vessel.b);
	}
	vessel.L = L/(static_cast<double>(inside.size())+0.00001f);

	if (VERBOSEFB)	{
		std::cout << "Update FB @ [" << vessel.mu[0] << ", " << vessel.mu[1] << ", " << vessel.mu[2] << "] width= "<< A << std::endl;
		std::cout << "\t# inside=" << inside.size() << ",  # Fg=" << f_tmp.size() << ",  # Bg=" << b_tmp.size() << std::endl;
		std::cout << "\tFest=" << vessel.f << ", Best=" << vessel.b << ",  MAD=" << vessel.MAD << ",  L=" << vessel.L << std::endl;
	}

	return vessel.L;

}


double SegInit::getMedian(std::vector<double> arr)	{
	unsigned int r = arr.size();
	double med = 0;
	if (r > 2)	{
		std::sort(arr.begin(), arr.end());
		if (r%2)	{
			med = arr[(r-1)/2];
		}
		else {
			med = 0.5 * (arr[r/2] + arr[(r+2)/2]);
		}
	}
	return med;
}


//special handling allows model to impinge on the boundary in order to complete trace.
int SegInit::handle_end( ImageType3D::Pointer im, double * pF, TVessel & vessel, TFacet * pFacets, int Num_Facets, TVertex * pUVertices, TVertex * pSVertices, TVertex * pCVertices, int Num_Vertices  )
{
    vessel.a3 /= 1.5;
    get_surface2(vessel,pFacets,Num_Facets,pUVertices,pSVertices,pCVertices,Num_Vertices);

    return interp3( im,  pF, pFacets, Num_Facets, pCVertices, 255.0 );
}


//fit the superellipsoid using coordinate descent, i.e., update parameters one at a time
//AS_RATIO controls how "stretched" the model can become
bool SegInit::fitSE (ImageType3D::Pointer im, TVessel& vessel, double iterations, double AS_RATIO) 
{

    gNum_facets = 0;
    gNum_points = 0;


    bool ret = 1;

    double L;
    double area;
    double total_area = 0.0;
	double total_f = 0.0;
	double total_f2 = 0.0;

    double F_new = vessel.f;
    double B_new = vessel.b;
    double FBden = 0.0;
    double q[4];

    //L = update_FBSimple(im,im_dim,vessel,2);
    vessel.a1 = .75*vessel.a1;
    vessel.a2 = .75*vessel.a2;
    //vessel.e1 = .75;
    TDamp damp;
    damp.dt = .003;
    damp.dt2 = .03;
    damp.dt_u[0] = .003;
    damp.dt_u[1] = .003;
    damp.dt_u[2] = .003;
	damp.sign_A[0] = 0.0;
	damp.sign_A[1] = 0.0;
	damp.sign_A[2] = 0.0;
	damp.sign_U[0] = 0.0;
	damp.sign_U[1] = 0.0;
	damp.sign_U[2] = 0.0;
	damp.sign_S[0] = 0.0;
	damp.sign_S[1] = 0.0;
	damp.sign_S[2] = 0.0;
	damp.sign_S[4] = 0.0;

    memcpy(damp.dt_a, damp.dt_u,sizeof(double)*3);

    generate_convex_hullq(vessel,19);
    L = update_FB(im,vessel);
    
	q[0] = _PI/2;
    q[1] = 0;
    q[2] = 1;
    q[3] = 0;
    axis_angle(q);
    memcpy(vessel.q1,q,sizeof(double)*4);

	TRMatrix R;

	for( int i = 1; i <= iterations; i++){

        if( i%25 == 0 || i == 10 ){
            update_FB(im,vessel);
            
            if( vessel.L < vessel.MAD*.5 && i > 25 && vessel.a3 > 10 )	{
				ret = 0;
			    break;
			}
			if ((vessel.b - vessel.f) < 2.0)	{
				ret = 0;
				break;
			}
        }

        FBden = vessel.b - vessel.f;
        if (vnl_math_abs(FBden) < 1.0)	{
			FBden = 1.0;
		}
		get_rotation(vessel,R);

        get_surface2(vessel,gFacets,gNum_facets,gVertexU,gVertexS,gVertexC,gNum_points);
        int SUCCESS = interp3DefBG( im, gF, gFacets, gNum_facets, gVertexC, 255.0 );
        if (SUCCESS==0)	{
			ret = 0;
			break;
		}

        total_area = 0.0;
        total_f = 0.0;
		total_f2 = 0.0;

        for( int fi = 0; fi < gNum_facets; fi++ ) {

            area = gFacets[fi].area;
            total_area += area;

            total_f += gF[fi];
			gF[fi] = fabs(vessel.f-gF[fi])-fabs(vessel.b-gF[fi]);
            gF[fi] = gF[fi]*area / FBden;
			total_f2 += gF[fi];

        }

        update_Quaternions2(vessel,gFacets,gNum_facets,gF,R,total_area,damp);

        if( i > .2*iterations && i < 3.0/5.0*iterations )
            update_e1(vessel,gFacets,gNum_facets,gF,R,damp);
           // update_e1(vessel,gFacets,gNum_facets,gF,R,damp);

        update_scale(vessel,gFacets,gNum_facets,gF,R,damp,AS_RATIO);
        update_mu(vessel,gFacets,gNum_facets,gF,R,damp);
        update_axes(vessel,0);

		if (VERBOSEITER)	{
			std::cout << "Iteration "<< i << std::endl;
			vessel.PrintSelf();
			std::cin.get();
		}
    }

    L = update_FB(im,vessel);


    TRMatrix RR;
    rotation_quat(vessel.q1, RR);
    vessel.R1[0] = RR.r[0][0];
    vessel.R1[1] = RR.r[1][0];
    vessel.R1[2] = RR.r[2][0];

    vessel.R2[0] = RR.r[0][1];
    vessel.R2[1] = RR.r[1][1];
    vessel.R2[2] = RR.r[2][1];

    vessel.R3[0] = RR.r[0][2];
    vessel.R3[1] = RR.r[1][2];
    vessel.R3[2] = RR.r[2][2];

    return ret;

}

double SegInit::getRegionMean (ImageType3D::Pointer im, TRMatrix& lim) {

	float pixelVal = 0.0;
	unsigned long cnt = 0;
	ImageType3D::SizeType im_dim = im->GetBufferedRegion().GetSize();

	for( double k = lim.r[2][0]; k <= lim.r[2][1]; k++){
		for( double i = lim.r[0][0]; i <= lim.r[0][1]; i++){
			for( double j = lim.r[1][0]; j <= lim.r[1][1]; j++){
				if( i >= 1 && i < im_dim[0] && j >= 1 && j < im_dim[1] && k >= 1 && k < im_dim[2] ){

					ImageType3D::IndexType ndx = {{(long int)i, (long int)j, (long int)k}};
					pixelVal += im->GetPixel(ndx);
					cnt ++;
				}
			}
		}
	}
	return (static_cast<double>(pixelVal)/static_cast<double>(cnt));
}

double SegInit::getMean(std::vector<double> arr)	{
	unsigned int r = arr.size();
	if (r==0)	{
		return 0;
	}
	double mean = 0;
	std::vector<double>::iterator it;
	for (it = arr.begin(); it!= arr.end(); it++ )	{
		mean += (*it);
	}
	return (mean/r);
}
