#include "SpineRing.h"
bool GlobalDetector::SpineRing::RingSE() {
	 // , bool seout)
	// This is from antise.m and Check_Match4.m
	/* we loop within a region around the dendrite modeled by a superellipsoid.
	So the donut-shaped ring is not circular but rather superellipsoidal
	from the inside. The outside is not too important except for its distance
	from the dendrite center. Usually this is RADOSCALE*radius.
	RADISCALE is a small scale to go above the dendrite surface.
	thickness is calculated from the superellipse centers.
	seg is the current superellipsoid structure in the dendrite modeling.
	*/

	double CVa1 = seg->a1*spr_RADISCALE;
	double CVa2 = seg->a2*spr_RADISCALE;
	double CVa3 = MAX(thickness, seg->a3);
	double CVe1 = seg->e1;
	double CVe2 = seg->e2;
	double A1x  = CVa1*spr_RADOSCALE*1.8;
	double A2x  = CVa2*spr_RADOSCALE*1.8;
	double A3x  = CVa3/2*1.8;
	//double A    = MAX(A1x, A2x)*2;
	
	TRMatrix Result;
	seg->rotation_quat(Result);
	//seg->transpose_matrix(Result);
	double mu[3];
	for (int ii=0; ii<3; ii++)
		mu[ii]=seg->mu[ii];
	double r00, r01,r02, r10, r11, r12, r20, r21, r22;

	r00 = Result.r[0][0];
	r01 = Result.r[0][1];
	r02 = Result.r[0][2];
	
	r10 = Result.r[1][0];
	r11 = Result.r[1][1];
	r12 = Result.r[1][2];
	
	r20 = Result.r[2][0];
	r21 = Result.r[2][1];
	r22 = Result.r[2][2];


	//TRMatrix lim;
	//lim.r[0][0] = (int)(mu[0]+.5)-(int)(A+.5);
	//lim.r[0][1] = (int)(mu[0]+.5)+(int)(A+.5);

	//lim.r[1][0] = (int)(mu[1]+.5)-(int)(A+.5);
	//lim.r[1][1] = (int)(mu[1]+.5)+(int)(A+.5);

	//lim.r[2][0] = (int)(mu[2]+.5)-(int)(A3x*2.0+.5);
	//lim.r[2][1] = (int)(mu[2]+.5)+(int)(A3x*2.0+.5);
	//lim.r[2][0] = (int)(mu[2]+.5)-(int)(A/2.0+.5);
	//lim.r[2][1] = (int)(mu[2]+.5)+(int)(A/2.0+.5);

	double x,y,z;
	double i,j,k;
	double step0 = Parent->up_sampling[0]; 
	double step1 = Parent->up_sampling[1]; 
	double step2 = Parent->up_sampling[2]; 
	//int loop_count = 0;
	//int loop_count2  =0;

	//double sampling = 1;

	inring_f_vals.reserve(spr_MAXCANDSIZE);
	//std::vector<double> b_tmp; 	b_tmp.reserve(MAXCANDSZ);
	inring_f_idx.reserve(spr_MAXCANDSIZE);
	DEBUGSTMT(inring_all_idx.reserve(5000));
	SpineImageType::SizeType im_dim = Parent->image->GetRequestedRegion().GetSize();
	//int f_count = 0;
	//int b_count = 0;
	double forgrnd  = seg->f;
	double bakgrnd  = seg->b;
	//double s1 = seg->a1;
	//double s2 = seg->a2;
	double s3 = seg->a3;
	//double e1 = seg->e1;
	//double e2 = 1.0;
	double D;

    //InterpolatorType::Pointer interp = InterpolatorType::New();
    //interp->SetInputImage(Parent->image);
	
	//InterpPointType ndx;
	SpineImageIndexType pixelIndex;
	PointType point;
	ImagePixelType pixelVal;

	for ( k = -A3x; k <= A3x; k+=step2) 
	{
		for( i = -A1x; i <= A1x; i+=step0)
		{
			for( j = -A2x; j <= A2x; j+=step1)
			{
				D = pow(pow(fabs(i)/CVa1,2.0/CVe2)+pow(fabs(j)/CVa2,2./CVe2),CVe2/CVe1)+pow(fabs(k)/CVa3,2.0/CVe1);
				if (D>spr_RADISCALE)
				{
					x = r00*i + r01*j + r02*k + mu[0] + r02*s3;
					y = r10*i + r11*j + r12*k + mu[1] + r12*s3;
					z = r20*i + r21*j + r22*k + mu[2] + r22*s3;
					if( x >= 0 && x < im_dim[0] && y >= 0 && y < im_dim[1] && z >= 0 && z < im_dim[2] )
					{		
						//ndx[0] = i; ndx[1]= j; ndx[2]= k;
						point[0] = x;
						point[1] = y; 
						point[2] = z; 						
						bool isInside = Parent->image->TransformPhysicalPointToIndex( point, pixelIndex );
						if ( isInside )
						{
							pixelVal = Parent->image->GetPixel( pixelIndex );
							if (abs(pixelVal-forgrnd)-abs(pixelVal-bakgrnd)<0) 
							{
								inring_f_vals.push_back(pixelVal);
								for (int ii=0;ii<3;ii++) 
								{
									//imidx[ii]=(long int)ndx[ii];
									if (Parent->bbmin[ii] > pixelIndex[ii])
										Parent->bbmin[ii] = pixelIndex[ii];
									if (Parent->bbmax[ii] < pixelIndex[ii])
										Parent->bbmax[ii] = pixelIndex[ii];
								}
								// saving pixelIndex rather than physical point!
								// Reason: needed for connected components binary image 
								// to save physical points (FUTURE: Hussein) :
								// need to save them in a points container in the SpCandidate
								// and after Connected Components: query what pt is
								// within which resultant connected comp. object.
								// The physical pts become the spine cand. pixels.
								// Those without a valid Conn comp obj. are discarded.
								// this adds a little accuracy but maybe is not worthwhile.
								inring_f_idx.push_back(pixelIndex);
								//inside[loop_count]=pixelVal;
								//inring_f_ndx_count++;
							}
							DEBUGSTMT(inring_all_idx.push_back(pixelIndex));
						}
					}
					//else {
					//	b_tmp.push_back(pixelVal); //b_tmp[b_count]=pixelVal;
					//	b_count++;
					//}
				}
			}
		}
	}
	return (!inring_f_vals.empty());   
}



	//if (seout) {
	//	//%just to give it a round shape on the outside (which does not affect
	//	//%results):
	//	Do =((abs(S(1,:))/A1x).^(2./CVe2)+...
	//		(abs(S(2,:))./A2x).^(2./CVe2))...
	//		.^(CVe2./CVe1) + ...
	//		(abs(S(3,:))./A3x).^(2./CVe1);
	//   
	//	pxlsin = S(:,(D>dist2center & Do<1));
	//}

