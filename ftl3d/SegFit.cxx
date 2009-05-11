#include "SegFit.h"
#include "myDebug.h"

#define VERBOSEFB 0
#define VERBOSEITER 0



//void SegInit::MainFunction (int nlhs, mxArray* *plhs, int nrhs, const mxArray* *prhs) {
bool SegFit::fitSE (ImageType3D::Pointer im, TVessel& vessel, double iterations, double AS_RATIO) {


    //if the mu is outside image then return 0 
	ImageType3D::IndexType ndx = {{static_cast<long>(vessel.mu[0]), static_cast<long>(vessel.mu[1]), static_cast<long>(vessel.mu[2]) }};
	if (!im->GetBufferedRegion().IsInside(ndx))	{
		return false;
		}

    gNum_facets = 0;
    gNum_points = 0;

    //double * im = mxGetPr(prhs[0]);
    //ImageType3D::SizeType im_dim = im->GetBufferedRegion().GetSize();

    bool ret = 1;

    double L;
    double area;
    double total_area = 0.0;
	double total_f = 0.0;
	double total_f2 = 0.0;

    double F_new = vessel.f;
    double B_new = vessel.b;
    double FBden = 0.0;
    
    //L = update_FBSimple(im,im_dim,vessel,2);
    //vessel.a1 = .75*vessel.a1;
    //vessel.a2 = .75*vessel.a2;
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

    //if( iterations > 0 )
    generate_convex_hullq(vessel,19);
    L = update_FBNew(im,vessel, 1);

	TRMatrix R;

    //mexPrintf("death 1\n");

	for( int i = 1; i <= iterations; i++){

        if( i%25 == 0 || i == 10 ){
            update_FBNew(im,vessel, 0);
            //update_FBSimple(im,vessel);
//            mexPrintf("%i\n",i);
            if( vessel.L < vessel.MAD*.5 && i > 25 && vessel.a3 > 10 )	{
				//std::cout << "EXIT CRITERION: Likelihood failed" << std::endl;
				ret = 0;
				//vessel.PrintSelf();
                break;
			}
			if ((vessel.b - vessel.f) < 2.0)	{
				//std::cout << "EXIT CRITERION ( "<< i << "): Contrast failed: Bg:" << vessel.b << "  Fg:" << vessel.f << std::endl;
				ret = 0;
				//vessel.PrintSelf();
				break;
			}
        }



        //F_new = vessel.f;
        //B_new = vessel.b;
        FBden = vessel.b - vessel.f;
        if (vnl_math_abs(FBden) < 1.0)	{
			FBden = 1.0;
		}
		get_rotation(vessel,R);

        get_surface2(vessel,gFacets,gNum_facets,gVertexU,gVertexS,gVertexC,gNum_points);
        int SUCCESS = interp3DefBG( im, gF, gFacets, gNum_facets, gVertexC, 255.0 );
        if (SUCCESS==0)	{
			//std::cout << "EXIT CRITERION: End of image reached" << std::endl;
			ret = 0;
			break;
		}

        //END = interp3( im, gF, gFacets, gNum_facets, gVertexC );
        //if(END){
        //    END = handle_end(im,gF,vessel,gFacets,gNum_facets,gVertexU,gVertexS,gVertexC,gNum_points);
        //    if(END){
        //        END=2;
        //        P("END of image")
        //        break;
        //    }else
        //        END=1;
        //}

        total_area = 0.0;
        total_f = 0.0;
		total_f2 = 0.0;

        for( int fi = 0; fi < gNum_facets; fi++ ) {

            area = gFacets[fi].area;
            total_area += area;

            total_f += gF[fi];
			gF[fi] = fabs(F_new-gF[fi])-fabs(B_new-gF[fi]);
            gF[fi] = gF[fi]*area / FBden;
			total_f2 += gF[fi];

        }


        update_Quaternions2(vessel,gFacets,gNum_facets,gF,R,total_area,damp);

        if( i > .2*iterations && i < 3.0/5.0*iterations )
            update_e1(vessel,gFacets,gNum_facets,gF,R,damp);

        update_scale(vessel,gFacets,gNum_facets,gF,R,damp,AS_RATIO);
        update_mu(vessel,gFacets,gNum_facets,gF,R,damp);
        update_axes(vessel,0);

		if (VERBOSEITER)	{
			std::cout << "Iteration "<< i << std::endl;
			vessel.PrintSelf();
			std::cin.get();
		}


    }

    L = update_FBNew(im,vessel, 0);
    //L = update_FBSimple(im,vessel);


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

