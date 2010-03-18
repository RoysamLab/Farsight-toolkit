


#include "ClebschGordan.h"

#define REAL float


/*************************************************************************************************
*
* void STderivReal[Forward/Backward](REAL *in, int sz[], int L,REAL *out, REAL factor = 1.0)
*
* REAL *in    :  volume of size (L*2 , sz[0], sz[1], sz[2])
* int sz[]    :  size array of length 3
* int L       :  rank of input, rank of spherical tensor is L-1
* REAL *out   :  preallocated output volume of size ((L+1)*2,  sz[0], sz[1], sz[2])
* REAL factor :  result is multiplied by this real factor
*
*
* Remarks: * complex numbers are stored interleaved. 
*	   * rank 0 tensor has also to be complex, i.e. of real dimension (2, sz[0],sz[1],sz[2])
*	   * m-indices of spherical tensor are negative and reversed, i.e. -(L-1),-(L-2)...,0
*      * for central differences the image interpreted periodically (circular)
*      * for backward/forward differences "no" (bad) boundary treatment
*
**************************************************************************************************/



void STderivReal(REAL *in, int sz[], int L,REAL *out, REAL factor = 1.0)
{

	REAL *real_out = out;
	REAL *imag_out = out+1;
	REAL *real_in = in;
	REAL *imag_in = in+1;

	int totsiz = sz[0]*sz[1]*sz[2];
	int sz1 = sz[0];
	int sz2 = sz[0]*sz[1];
	int w = sz[0];
	int h = sz[1];
	int d = sz[2];

	int start = sz[1]*sz[0] + sz[0] + 1;
	int end = totsiz-start;


	int TL = L+1;
	int TLS = L;
	for (int ms = -L; ms <= 0 ;ms++)
	{
		REAL cg_m1 = factor*sqrt( (REAL) (L-ms)*(L-ms-1) /((2*L)*(2*L-1)) ) ;
		REAL cg_0 =  factor*sqrt( (REAL) 2*(L+ms)*(L-ms) /(L*(2*L-1)) );
		REAL cg_p1 = factor*sqrt( (REAL) (L+ms)*(L+ms-1) /((2*L)*(2*L-1)) );

		for (int z = 0; z < sz[2]; z++)
		for (int y = 0; y < sz[1]; y++)
		for (int x = 0; x < sz[0]; x++)
		{
			int k = x+sz[0]*(y+sz[1]*z);

			int dx,dx_,dy,dy_,dz,dz_;

			dx = 2*((x+1)%w + w*(y+h*z))*TLS;
			dx_ = 2*((x-1+w)%w + w*(y+h*z))*TLS;
			dy = 2*(x + w*((y+1)%h+h*z))*TLS;
			dy_ = 2*(x + w*((y-1+h)%h+h*z))*TLS;
			dz = 2*(x + w*(y+h*((z+1)%d)))*TLS;
			dz_ = 2*(x + w*(y+h*((z-1+d)%d)))*TLS;

			int m = 2*(ms+(L-1));
			int iout = 2*(k*TL+ms+L);

			real_out[iout] = 0;	
			imag_out[iout] = 0;

			if (ms-1 >= -(L-1))
			{
				real_out[iout] += cg_p1 * ((real_in[dx+m-2] - real_in[dx_+m-2]) + (imag_in[dy+m-2] - imag_in[dy_+m-2]));
				imag_out[iout] += cg_p1 * ((imag_in[dx+m-2] - imag_in[dx_+m-2]) - (real_in[dy+m-2] - real_in[dy_+m-2]));
			}
			if (ms <= L-1 && ms >= -(L-1))		
			{
				real_out[iout] += cg_0  * (real_in[dz+m] - real_in[dz_+m]);
				imag_out[iout] += cg_0  * (imag_in[dz+m] - imag_in[dz_+m]);
			}
			if (ms+1 <= (L-1))
			{
				if (ms == 0)
				{

					real_out[iout] += cg_m1 * ((real_in[dx+m-2] - real_in[dx_+m-2]) + (imag_in[dy+m-2] - imag_in[dy_+m-2]));
					imag_out[iout] -= cg_m1 * ((imag_in[dx+m-2] - imag_in[dx_+m-2]) - (real_in[dy+m-2] - real_in[dy_+m-2]));
				}
				else
				{
					real_out[iout] -= cg_m1 * ((real_in[dx+m+2] - real_in[dx_+m+2]) - (imag_in[dy+m+2] - imag_in[dy_+m+2]));
					imag_out[iout] -= cg_m1 * ((imag_in[dx+m+2] - imag_in[dx_+m+2]) + (real_in[dy+m+2] - real_in[dy_+m+2]));
				}
			}

		}


	}

}


void STderivRealForward(REAL *in, int sz[], int L,REAL *out, REAL factor = 1.0)
{

	REAL *real_out = out;
	REAL *imag_out = out+1;
	REAL *real_in = in;
	REAL *imag_in = in+1;

	int totsiz = sz[0]*sz[1]*sz[2];
	int sz1 = sz[0];
	int sz2 = sz[0]*sz[1];

	int start = sz[1]*sz[0] + sz[0] + 1;
	int end = totsiz-start;


	int TL = L+1;
	int TLS = L;
	for (int ms = -L; ms <= 0 ;ms++)
	{
		REAL cg_m1 = factor*sqrt( (REAL) (L-ms)*(L-ms-1) /((2*L)*(2*L-1)) ) ;
		REAL cg_0 =  factor*sqrt( (REAL) 2*(L+ms)*(L-ms) /(L*(2*L-1)) );
		REAL cg_p1 = factor*sqrt( (REAL) (L+ms)*(L+ms-1) /((2*L)*(2*L-1)) );

		for (int k = start;k < end;k++)
		{
			int dx = 2*(k+1)*TLS;
			int dx_ = 2*(k)*TLS;
			int dy = 2*(k+sz1)*TLS;
			int dy_ = 2*(k)*TLS;
			int dz = 2*(k+sz2)*TLS;
			int dz_ = 2*(k)*TLS;

			int m = 2*(ms+(L-1));
			int iout = 2*(k*TL+ms+L);

			real_out[iout] = 0;	
			imag_out[iout] = 0;

			if (ms-1 >= -(L-1))
			{
				real_out[iout] += cg_p1 * ((real_in[dx+m-2] - real_in[dx_+m-2]) + (imag_in[dy+m-2] - imag_in[dy_+m-2]));
				imag_out[iout] += cg_p1 * ((imag_in[dx+m-2] - imag_in[dx_+m-2]) - (real_in[dy+m-2] - real_in[dy_+m-2]));
			}
			if (ms <= L-1 && ms >= -(L-1))		
			{
				real_out[iout] += cg_0  * (real_in[dz+m] - real_in[dz_+m]);
				imag_out[iout] += cg_0  * (imag_in[dz+m] - imag_in[dz_+m]);
			}
			if (ms+1 <= (L-1))
			{
				if (ms == 0)
				{

					real_out[iout] += cg_m1 * ((real_in[dx+m-2] - real_in[dx_+m-2]) + (imag_in[dy+m-2] - imag_in[dy_+m-2]));
					imag_out[iout] -= cg_m1 * ((imag_in[dx+m-2] - imag_in[dx_+m-2]) - (real_in[dy+m-2] - real_in[dy_+m-2]));
				}
				else
				{
					real_out[iout] -= cg_m1 * ((real_in[dx+m+2] - real_in[dx_+m+2]) - (imag_in[dy+m+2] - imag_in[dy_+m+2]));
					imag_out[iout] -= cg_m1 * ((imag_in[dx+m+2] - imag_in[dx_+m+2]) + (real_in[dy+m+2] - real_in[dy_+m+2]));
				}
			}

		}


	}

}




void STderivRealBackward(REAL *in, int sz[], int L,REAL *out, REAL factor = 1.0)
{

	REAL *real_out = out;
	REAL *imag_out = out+1;
	REAL *real_in = in;
	REAL *imag_in = in+1;

	int totsiz = sz[0]*sz[1]*sz[2];
	int sz1 = sz[0];
	int sz2 = sz[0]*sz[1];

	int start = sz[1]*sz[0] + sz[0] + 1;
	int end = totsiz-start;


	int TL = L+1;
	int TLS = L;
	for (int ms = -L; ms <= 0 ;ms++)
	{
		REAL cg_m1 = factor*sqrt( (REAL) (L-ms)*(L-ms-1) /((2*L)*(2*L-1)) ) ;
		REAL cg_0 =  factor*sqrt( (REAL) 2*(L+ms)*(L-ms) /(L*(2*L-1)) );
		REAL cg_p1 = factor*sqrt( (REAL) (L+ms)*(L+ms-1) /((2*L)*(2*L-1)) );

		for (int k = start;k < end;k++)
		{
			int dx = 2*(k)*TLS;
			int dx_ = 2*(k-1)*TLS;
			int dy = 2*(k)*TLS;
			int dy_ = 2*(k-sz1)*TLS;
			int dz = 2*(k)*TLS;
			int dz_ = 2*(k-sz2)*TLS;

			int m = 2*(ms+(L-1));
			int iout = 2*(k*TL+ms+L);

			real_out[iout] = 0;	
			imag_out[iout] = 0;

			if (ms-1 >= -(L-1))
			{
				real_out[iout] += cg_p1 * ((real_in[dx+m-2] - real_in[dx_+m-2]) + (imag_in[dy+m-2] - imag_in[dy_+m-2]));
				imag_out[iout] += cg_p1 * ((imag_in[dx+m-2] - imag_in[dx_+m-2]) - (real_in[dy+m-2] - real_in[dy_+m-2]));
			}
			if (ms <= L-1 && ms >= -(L-1))		
			{
				real_out[iout] += cg_0  * (real_in[dz+m] - real_in[dz_+m]);
				imag_out[iout] += cg_0  * (imag_in[dz+m] - imag_in[dz_+m]);
			}
			if (ms+1 <= (L-1))
			{
				if (ms == 0)
				{

					real_out[iout] += cg_m1 * ((real_in[dx+m-2] - real_in[dx_+m-2]) + (imag_in[dy+m-2] - imag_in[dy_+m-2]));
					imag_out[iout] -= cg_m1 * ((imag_in[dx+m-2] - imag_in[dx_+m-2]) - (real_in[dy+m-2] - real_in[dy_+m-2]));
				}
				else
				{
					real_out[iout] -= cg_m1 * ((real_in[dx+m+2] - real_in[dx_+m+2]) - (imag_in[dy+m+2] - imag_in[dy_+m+2]));
					imag_out[iout] -= cg_m1 * ((imag_in[dx+m+2] - imag_in[dx_+m+2]) + (real_in[dy+m+2] - real_in[dy_+m+2]));
				}
			}

		}


	}

}






void STderivRealDown(REAL *in, int sz[], int L,REAL *out, REAL factor = 1.0)
{

	REAL *real_out = out;
	REAL *imag_out = out+1;
	REAL *real_in = in;
	REAL *imag_in = in+1;


	int totsiz = sz[0]*sz[1]*sz[2];
	int sz1 = sz[0];
	int sz2 = sz[0]*sz[1];
	int w = sz[0];
	int h = sz[1];
	int d = sz[2];


	int start = sz[1]*sz[0] + sz[0] + 1;
	int end = totsiz-start;


	int TL = L+1;
	int TLS = L;
	for (int ms = -(L-1); ms <= 0 ;ms++)
	{
		REAL cg_m1 = factor*sqrt( (REAL) (L+ms)*(L+ms+1)  /(2*L*(2*L+1)) ) ;
		REAL cg_0 =  -factor*sqrt( (REAL) 2*(L-ms)*(L+ms) / (L*(2*L+1))  );
		REAL cg_p1 = factor*sqrt( (REAL) (L-ms)*(L-ms+1)  /(2*L*(2*L+1)) );




		for (int z = 0; z < sz[2]; z++)
		for (int y = 0; y < sz[1]; y++)
		for (int x = 0; x < sz[0]; x++)
		{
			int k = x+w*(y+h*z);


			int dx,dx_,dy,dy_,dz,dz_;

			dx = 2*((x+1)%w + w*(y+h*z))*TL;
			dx_ = 2*((x-1+w)%w + w*(y+h*z))*TL;
			dy = 2*(x + w*((y+1)%h+h*z))*TL;
			dy_ = 2*(x + w*((y-1+h)%h+h*z))*TL;
			dz = 2*(x + w*(y+h*((z+1)%d)))*TL;
			dz_ = 2*(x + w*(y+h*((z-1+d)%d)))*TL;

			int m = 2*(ms+L);
			int iout = 2*(k*TLS+ms+(L-1));

			real_out[iout] = cg_p1 * ((real_in[dx+m-2] - real_in[dx_+m-2]) + (imag_in[dy+m-2] - imag_in[dy_+m-2]));
			imag_out[iout] = cg_p1 * ((imag_in[dx+m-2] - imag_in[dx_+m-2]) - (real_in[dy+m-2] - real_in[dy_+m-2]));

			real_out[iout] += cg_0  * (real_in[dz+m] - real_in[dz_+m]);
			imag_out[iout] += cg_0  * (imag_in[dz+m] - imag_in[dz_+m]);
			if (ms == 0)
			{	

				real_out[iout] += cg_m1 * ((real_in[dx+m-2] - real_in[dx_+m-2]) + (imag_in[dy+m-2] - imag_in[dy_+m-2]));
				imag_out[iout] += cg_m1 * ((imag_in[dx+m-2] - imag_in[dx_+m-2]) - (real_in[dy+m-2] - real_in[dy_+m-2]));
			}
			else
			{
	
				real_out[iout] -= cg_m1 * ((real_in[dx+m+2] - real_in[dx_+m+2]) - (imag_in[dy+m+2] - imag_in[dy_+m+2]));
				imag_out[iout] -= cg_m1 * ((imag_in[dx+m+2] - imag_in[dx_+m+2]) + (real_in[dy+m+2] - real_in[dy_+m+2]));
			}

		}


	}

}




void STderivRealDownForward(REAL *in, int sz[], int L,REAL *out, REAL factor = 1.0)
{

	REAL *real_out = out;
	REAL *imag_out = out+1;
	REAL *real_in = in;
	REAL *imag_in = in+1;


	int totsiz = sz[0]*sz[1]*sz[2];
	int sz1 = sz[0];
	int sz2 = sz[0]*sz[1];

	int start = sz[1]*sz[0] + sz[0] + 1;
	int end = totsiz-start;


	int TL = L+1;
	int TLS = L;
	for (int ms = -(L-1); ms <= 0 ;ms++)
	{
		REAL cg_m1 = factor*sqrt( (REAL) (L+ms)*(L+ms+1)  /(2*L*(2*L+1)) ) ;
		REAL cg_0 =  -factor*sqrt( (REAL) 2*(L-ms)*(L+ms) / (L*(2*L+1))  );
		REAL cg_p1 = factor*sqrt( (REAL) (L-ms)*(L-ms+1)  /(2*L*(2*L+1)) );

		for (int k = start;k < end;k++)
		{
			int dx = 2*(k+1)*TL;
			int dx_ = 2*(k)*TL;
			int dy = 2*(k+sz1)*TL;
			int dy_ = 2*(k)*TL;
			int dz = 2*(k+sz2)*TL;
			int dz_ = 2*(k)*TL;

			int m = 2*(ms+L);
			int iout = 2*(k*TLS+ms+(L-1));

			real_out[iout] = cg_p1 * ((real_in[dx+m-2] - real_in[dx_+m-2]) + (imag_in[dy+m-2] - imag_in[dy_+m-2]));
			imag_out[iout] = cg_p1 * ((imag_in[dx+m-2] - imag_in[dx_+m-2]) - (real_in[dy+m-2] - real_in[dy_+m-2]));

			real_out[iout] += cg_0  * (real_in[dz+m] - real_in[dz_+m]);
			imag_out[iout] += cg_0  * (imag_in[dz+m] - imag_in[dz_+m]);
			if (ms == 0)
			{	

				real_out[iout] += cg_m1 * ((real_in[dx+m-2] - real_in[dx_+m-2]) + (imag_in[dy+m-2] - imag_in[dy_+m-2]));
				imag_out[iout] += cg_m1 * ((imag_in[dx+m-2] - imag_in[dx_+m-2]) - (real_in[dy+m-2] - real_in[dy_+m-2]));
			}
			else
			{
	
				real_out[iout] -= cg_m1 * ((real_in[dx+m+2] - real_in[dx_+m+2]) - (imag_in[dy+m+2] - imag_in[dy_+m+2]));
				imag_out[iout] -= cg_m1 * ((imag_in[dx+m+2] - imag_in[dx_+m+2]) + (real_in[dy+m+2] - real_in[dy_+m+2]));
			}

		}


	}

}





void STderivRealDownBackward(REAL *in, int sz[], int L,REAL *out, REAL factor = 1.0)
{

	REAL *real_out = out;
	REAL *imag_out = out+1;
	REAL *real_in = in;
	REAL *imag_in = in+1;


	int totsiz = sz[0]*sz[1]*sz[2];
	int sz1 = sz[0];
	int sz2 = sz[0]*sz[1];

	int start = sz[1]*sz[0] + sz[0] + 1;
	int end = totsiz-start;


	int TL = L+1;
	int TLS = L;
	for (int ms = -(L-1); ms <= 0 ;ms++)
	{
		REAL cg_m1 = factor*sqrt( (REAL) (L+ms)*(L+ms+1)  /(2*L*(2*L+1)) ) ;
		REAL cg_0 =  -factor*sqrt( (REAL) 2*(L-ms)*(L+ms) / (L*(2*L+1))  );
		REAL cg_p1 = factor*sqrt( (REAL) (L-ms)*(L-ms+1)  /(2*L*(2*L+1)) );

		for (int k = start;k < end;k++)
		{
			int dx = 2*(k)*TL;
			int dx_ = 2*(k-1)*TL;
			int dy = 2*(k)*TL;
			int dy_ = 2*(k-sz1)*TL;
			int dz = 2*(k)*TL;
			int dz_ = 2*(k-sz2)*TL;

			int m = 2*(ms+L);
			int iout = 2*(k*TLS+ms+(L-1));

			real_out[iout] = cg_p1 * ((real_in[dx+m-2] - real_in[dx_+m-2]) + (imag_in[dy+m-2] - imag_in[dy_+m-2]));
			imag_out[iout] = cg_p1 * ((imag_in[dx+m-2] - imag_in[dx_+m-2]) - (real_in[dy+m-2] - real_in[dy_+m-2]));

			real_out[iout] += cg_0  * (real_in[dz+m] - real_in[dz_+m]);
			imag_out[iout] += cg_0  * (imag_in[dz+m] - imag_in[dz_+m]);
			if (ms == 0)
			{	

				real_out[iout] += cg_m1 * ((real_in[dx+m-2] - real_in[dx_+m-2]) + (imag_in[dy+m-2] - imag_in[dy_+m-2]));
				imag_out[iout] += cg_m1 * ((imag_in[dx+m-2] - imag_in[dx_+m-2]) - (real_in[dy+m-2] - real_in[dy_+m-2]));
			}
			else
			{
	
				real_out[iout] -= cg_m1 * ((real_in[dx+m+2] - real_in[dx_+m+2]) - (imag_in[dy+m+2] - imag_in[dy_+m+2]));
				imag_out[iout] -= cg_m1 * ((imag_in[dx+m+2] - imag_in[dx_+m+2]) + (real_in[dy+m+2] - real_in[dy_+m+2]));
			}

		}


	}

}



void STmultiply_withCG(REAL *A1, int L1, REAL *A2, int L2, 
		REAL *out, int L3, int imsize, REAL factor = 1.0)
{


	REAL cg[1000];
        int indextups[3*1000];    
        int cnt = 0;
	for (int m = -(L3-1); m <= 0; m++)
		for(int ms = -(L1-1);ms <= L1-1;ms++)
		{
			if (abs(m-ms) <= L2-1)
			{
				int errflag;
				cg[cnt] = factor*ClebschGordan(L1-1,ms,L2-1,m-ms,L3-1,m, errflag);
//                printf("cg[%d] = %f\n",cnt,cg[cnt]);
				indextups[3*cnt] = 2*ms;
				indextups[3*cnt+1] = 2*(m-ms);
				indextups[3*cnt+2] = 2* m;
				cnt++;
			}
		}


	int p1 = -2;
	int p2 = -2;
	int p3 = -2;

	REAL *A1real = A1;
	REAL *A1imag = A1+1;
	REAL *A2real = A2;
	REAL *A2imag = A2+1;
	REAL *Oreal = out;	
	REAL *Oimag = out+1;	

	int L1_2 = 2*L1;
	int L2_2 = 2*L2;
	int L3_2 = 2*L3;
//    printf("cnt = %d L1_2 =%d L2_2 =%d L3_2 = %d\n", cnt, L1_2, L2_2, L3_2);
	for (int i = 0; i < imsize;i++)
	{
		p1 += L1_2;
		p2 += L2_2;
		p3 += L3_2;

		for (int m = 0; m < 2*L3-1; m++)
		{
			out[p3-m] = 0;	
		}

		for (int m = 0; m < cnt; m++)
		{
			int M1 = indextups[3*m];
			int M2 = indextups[3*m+1];
			int M3 = indextups[3*m+2];
//         if(A1real[p1+M1] > 10e10 || A1imag[p2+M2] > 10e10)
//         {
//            printf("A1real[] = %f A1imag[] = %f\n",A1real[p1+M1],A1imag[p2+M2]);
//         }

			if (M1 <= 0)
			{
				if (M2 <= 0)
				{
					Oreal[p3+M3] += cg[m] * (A1real[p1+M1]*A2real[p2+M2] - A1imag[p1+M1]*A2imag[p2+M2]);
					Oimag[p3+M3] += cg[m] * (A1real[p1+M1]*A2imag[p2+M2] + A1imag[p1+M1]*A2real[p2+M2]);
				}
				else
				{			
					if (!(M2&2))
					{		
					Oreal[p3+M3] += cg[m] * (A1real[p1+M1]*A2real[p2-M2] + A1imag[p1+M1]*A2imag[p2-M2]);
					Oimag[p3+M3] -= cg[m] * (A1real[p1+M1]*A2imag[p2-M2] - A1imag[p1+M1]*A2real[p2-M2]);
					}
					else
					{		
					Oreal[p3+M3] -= cg[m] * (A1real[p1+M1]*A2real[p2-M2] + A1imag[p1+M1]*A2imag[p2-M2]);
					Oimag[p3+M3] += cg[m] * (A1real[p1+M1]*A2imag[p2-M2] - A1imag[p1+M1]*A2real[p2-M2]);
					}
				}
			}

			else
			{
				if (M2 <= 0)
				{
					if (!(M1&2))
					{	
					Oreal[p3+M3] += cg[m] * (A1real[p1-M1]*A2real[p2+M2] + A1imag[p1-M1]*A2imag[p2+M2]);
					Oimag[p3+M3] += cg[m] * (A1real[p1-M1]*A2imag[p2+M2] - A1imag[p1-M1]*A2real[p2+M2]);
					}
					else
					{	
					Oreal[p3+M3] -= cg[m] * (A1real[p1-M1]*A2real[p2+M2] + A1imag[p1-M1]*A2imag[p2+M2]);
					Oimag[p3+M3] -= cg[m] * (A1real[p1-M1]*A2imag[p2+M2] - A1imag[p1-M1]*A2real[p2+M2]);
					}
				}
				else
				{
					if (!(M3&2))
					{	
					Oreal[p3+M3] += cg[m] * (A1real[p1-M1]*A2real[p2-M2] - A1imag[p1-M1]*A2imag[p2-M2]);
					Oimag[p3+M3] -= cg[m] * (A1real[p1-M1]*A2imag[p2-M2] + A1imag[p1-M1]*A2real[p2-M2]);
					}
					else
					{	
					Oreal[p3+M3] -= cg[m] * (A1real[p1-M1]*A2real[p2-M2] - A1imag[p1-M1]*A2imag[p2-M2]);
					Oimag[p3+M3] += cg[m] * (A1real[p1-M1]*A2imag[p2-M2] + A1imag[p1-M1]*A2real[p2-M2]);
					}
		
				}

			}


		} 
	}

}






void STmultiplyFourier(REAL *A1, int L1, REAL *A2, int L2, 
		REAL *out, int L3, int imdim[],	
		int *indextups, REAL *cg,int cnt)
{


	REAL *A1real = A1;
	REAL *A1imag = A1+1;
	REAL *A2real = A2;
	REAL *A2imag = A2+1;
	REAL *Oreal = out;	
	REAL *Oimag = out+1;	

	int p1 = -2;
	int p2 = -2;
	int p3 = -2;

	int nx = imdim[0];
	int ny = imdim[1];
	int nz = imdim[2];

	int L1_2 = 2*L1;
	int L2_2 = 2*L2;
	int L3_2 = 2*L3;


	int nxny = nx*ny;


	for (int z = 0; z < nz; z++)
	{
	int sp = 1 + ((nz-z)%nz) * nxny;

	for (int y = 0; y < ny; y++)
	{
	int sp_ = sp + ((ny-y)%ny)*nx;


	for (int x = 0; x < nx; x++)
	{
		p1 += L1_2;
		p2 += L2_2;
		p3 += L3_2;
	
		int sp__ = sp_ + (nx-x)%nx ;
		int sp1 = 2*(sp__ * L1 - 1) ;
		int sp2 = 2*(sp__ * L2 - 1) ;
/*
		int sp = ((nz-z)%nz) *nx*ny + ((ny-y)%ny)*nx + (nx-x)%nx;
		int sp1 = (sp*L1 +L1-1)*2;
		int sp2 = (sp*L2 +L2-1)*2;
*/
		for (int m = 0; m < L3_2-1; m++)
		{
			out[p3-m] = 0;	
		}

		for (int m = 0; m < cnt; m++)
		{
			int M1 = indextups[3*m];
			int M2 = indextups[3*m+1];
			int M3 = indextups[3*m+2];
//            cg[m] = 1;
			if (M1 <= 0)
			{
				if (M2 <= 0)
				{
					Oreal[p3+M3] += cg[m] * (A1real[p1+M1]*A2real[p2+M2] - A1imag[p1+M1]*A2imag[p2+M2]);
					Oimag[p3+M3] += cg[m] * (A1real[p1+M1]*A2imag[p2+M2] + A1imag[p1+M1]*A2real[p2+M2]);
				}
				else
				{			
					if (!(M2&2))
					//if ((M2/2)%2 == 0)
					{
					Oreal[p3+M3] += cg[m] * (A1real[p1+M1]*A2real[sp2-M2] + A1imag[p1+M1]*A2imag[sp2-M2]);
					Oimag[p3+M3] -= cg[m] * (A1real[p1+M1]*A2imag[sp2-M2] - A1imag[p1+M1]*A2real[sp2-M2]);
					}
					else
					{		
					Oreal[p3+M3] -= cg[m] * (A1real[p1+M1]*A2real[sp2-M2] + A1imag[p1+M1]*A2imag[sp2-M2]);
					Oimag[p3+M3] += cg[m] * (A1real[p1+M1]*A2imag[sp2-M2] - A1imag[p1+M1]*A2real[sp2-M2]);
					}
				}
			}

			else
			{
				if (M2 <= 0)
				{
					if (!(M1&2))
					{	
					Oreal[p3+M3] += cg[m] * (A1real[sp1-M1]*A2real[p2+M2] + A1imag[sp1-M1]*A2imag[p2+M2]);
					Oimag[p3+M3] += cg[m] * (A1real[sp1-M1]*A2imag[p2+M2] - A1imag[sp1-M1]*A2real[p2+M2]);
					}
					else
					{	
					Oreal[p3+M3] -= cg[m] * (A1real[sp1-M1]*A2real[p2+M2] + A1imag[sp1-M1]*A2imag[p2+M2]);
					Oimag[p3+M3] -= cg[m] * (A1real[sp1-M1]*A2imag[p2+M2] - A1imag[sp1-M1]*A2real[p2+M2]);
					}
				}
				else
				{
					if (!(M3&2))
					{	
					Oreal[p3+M3] += cg[m] * (A1real[sp1-M1]*A2real[sp2-M2] - A1imag[sp1-M1]*A2imag[sp2-M2]);
					Oimag[p3+M3] -= cg[m] * (A1real[sp1-M1]*A2imag[sp2-M2] + A1imag[sp1-M1]*A2real[sp2-M2]);
					}
					else
					{	
					Oreal[p3+M3] -= cg[m] * (A1real[sp1-M1]*A2real[sp2-M2] - A1imag[sp1-M1]*A2imag[sp2-M2]);
					Oimag[p3+M3] += cg[m] * (A1real[sp1-M1]*A2imag[sp2-M2] + A1imag[sp1-M1]*A2real[sp2-M2]);
					}
		
				}

			}

		} 
	}	}	}
	

}







