// Extract the ridges and valleys feature in 3D vector field
// --- Input: 3D vector field
// --- Output: 3D image-scalar field
// --- Author: Xiaosong Yuan, RPI
// --- Modified Date: 10/6/2005

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <math.h>

using namespace std;

struct  Vector3D
{
	float x;
	float y;
	float z;
};

void PartialDerivative(float *Is, float *Isd, int direc, int L, int M, int N);
void PartialDerivative1(float *Is, float *Isd, int direc, int L, int M, int N);
void ComputeRotMatrix(float RotateMatrix[3][3], Vector3D v);
void RotMatrixFromAngle(float RMatrix[3][3], float phi, float theta, float psi);
void Transpose(float MatTransp[3][3], float Mat[3][3]);
void Matrix3Multiply(float Mat[3][3], float Mat1[3][3], float Mat2[3][3]);

int main (int argc, char *argv[])
{
  ifstream fin;
  FILE *fout;
  char *infilename;
  Vector3D vecin, gradient0;
  float *Iu, *Iv, *Iw;
  float *Iuu, *Iuv, *Iuw, *Ivv, *Ivw, *Iww;
  float *curv;
  int *f;
  long idx, iidx, slsz, sz;
  int i, j, k;
  int ii, jj, kk;
  int x,y,z;
  int cc;
  int L, M, N;
  float k1, k2, gLength;
  float RotMatrix[3][3];
  float RotMatrixTransp[3][3];
  float Hessian[3][3];
  float HessianPrime[3][3];
  float HessianPrime1[3][3];
  int numSeeds;

  if (argc < 6)
  {
    printf("Usage: %s <Input vectors> <xs> <ys> <zs> <Output scalars>.\n", argv[0]);
    exit(1);
  }

  infilename = new char[80];
  infilename = argv[1];

  fin.open(infilename);
  if (!fin)  {
     cerr << "couldn't open " << infilename << " for input" << endl;
     return -1;
  }

  L = atoi(argv[2]);
  M = atoi(argv[3]);
  N = atoi(argv[4]);

  slsz = L*M;        // slice size
  sz = slsz*N;

 /* if ((fout = fopen(argv[5],"w")) == NULL)
  {
    printf("Cannot open %s for writing\n",argv[5]);
    exit(1);
  }
  */

  errno_t err; 
 if((err=fopen_s(&fout,argv[5],"w"))!=NULL)
			{printf("Input file open error!\n");
			 exit(-1);
			}

  Iu = new float[L*M*N];
  Iv = new float[L*M*N];
  Iw = new float[L*M*N];
  Iuu = new float[L*M*N];
  Ivv = new float[L*M*N];
  Iww = new float[L*M*N];
  Iuv = new float[L*M*N];
  Iuw = new float[L*M*N];
  Ivw = new float[L*M*N];
  curv = new float[L*M*N];
  f = new int[L*M*N];

  for(idx=0; idx<sz; idx++)   {  //Initialize to zeros
	  f[idx]=0;
	  Iu[idx]=0;
	  Iv[idx]=0;
	  Iw[idx]=0;
	  curv[idx]=0;
  }

  fin >> x >> y >> z >> vecin.x >> vecin.y >> vecin.z;
  while (!fin.eof() ) 
  {
	    idx = z*slsz + y*L +x;
	    Iu[idx] = vecin.x;
	    Iv[idx] = vecin.y;
	    Iw[idx] = vecin.z;
	    f[idx]=2;
		fin >> x >> y >> z >> vecin.x >> vecin.y >> vecin.z;
	}

  // make surface point f[]=1
  for (k = 1; k < N-1; k++)
     for (j = 1; j < M-1; j++)
        for (i = 1; i < L-1; i++) {
	    idx = k*slsz + j*L +i;
	    if(f[idx]==2) {
	           cc = 0;
		   for (kk=-1; kk<=1; kk++)
		       for (jj=-1; jj<=1; jj++)
		           for (ii=-1; ii<=1; ii++) {
			         iidx = (k+kk)*slsz + (j+jj)*L +(i+ii);
				 if (f[iidx]==0) cc++;
			   }
		   if (cc>=1) f[idx] = 1;
	    }
	}


  PartialDerivative1(Iu, Iuu, 1, L, M, N);
  PartialDerivative1(Iv, Ivv, 2, L, M, N);
  PartialDerivative1(Iw, Iww, 3, L, M, N);
  PartialDerivative1(Iu, Iuv, 2, L, M, N);
  PartialDerivative1(Iu, Iuw, 3, L, M, N);
  PartialDerivative1(Iv, Ivw, 3, L, M, N);


  int DisAway = 2;
  for (k = DisAway; k < N-DisAway; k++)
     for (j = DisAway; j < M-DisAway; j++)
        for (i = DisAway; i < L-DisAway; i++) {
	     idx = k*slsz + j*L + i;
	     if (f[idx] !=2) continue;
	     gradient0.x = Iu[idx];
	     gradient0.y = Iv[idx];
	     gradient0.z = Iw[idx];
	     Hessian[0][0] = Iuu[idx];
	     Hessian[0][1] = Iuv[idx];
	     Hessian[0][2] = Iuw[idx];
	     Hessian[1][0] = Iuv[idx];
	     Hessian[1][1] = Ivv[idx];
	     Hessian[1][2] = Ivw[idx];
	     Hessian[2][0] = Iuw[idx];
	     Hessian[2][1] = Ivw[idx];
	     Hessian[2][2] = Iww[idx];
	     ComputeRotMatrix(RotMatrix, gradient0);
	     Transpose(RotMatrixTransp, RotMatrix);
	     Matrix3Multiply(HessianPrime1, RotMatrixTransp, Hessian);
	     Matrix3Multiply(HessianPrime, HessianPrime1, RotMatrix);
	     k1 = (float)0.5*(HessianPrime[1][1]+HessianPrime[2][2])
	         +(float)0.5*sqrt(pow((Hessian[1][1]-Hessian[2][2]),2)+4*Hessian[1][2]*Hessian[2][1]);
	     k2 = (float)0.5*(HessianPrime[1][1]+HessianPrime[2][2])
	         -(float)0.5*sqrt(pow((Hessian[1][1]-Hessian[2][2]),2)+4*Hessian[1][2]*Hessian[2][1]);

		 ///  add (float) to remove warning 

	     gLength = sqrt(gradient0.x*gradient0.x +gradient0.y*gradient0.y +gradient0.z*gradient0.z);
	     //if (k2>k1) k1=k2;
	     if(gLength !=0) curv[idx] = -(k1)/gLength;   else curv[idx]=9999;//set large to be included in skeleton
	     //if(gLength !=0) curv[idx] = -(k1)/pow(gLength,1.5);   else curv[idx]=10000;//set large to be included in skeleton
		 

	     //case 1: get neg maximum
	     if (curv[idx]>10000) curv[idx]=10000;

	     if (curv[idx]< 5 ) curv[idx]= 0; // first phase: threshold the curvature value (larger->fewer)
		            //< 2.0

		   // thres =   3 -6712   seeds    for "Ecadherin2_T60DkR30.1392x1024x9.Aniso_k800d02t2.grdt.vec"
		   //          20 - 283   seeds
		   //          60 -  52   seeds

		   // thres =  50 - 186   seeds    for "Trach6A_T10DkR1k.512x512x26.grdt.vec"
		   //          20 - 223   seeds
		   //          15 - 234   seeds  
		   //          10 - 251   seeds
		   //           7 - 275   seeds


		   // thres =  0  -  57   seeds    for "phantomSpines.64x40x31.grdt.vec"
		   // thres = -1  -  72   seeds
		   // thres = -3  -  92   seeds 

		   // thres =  0  - 43504 seeds    for "tlapse330_T2DkR1k.512x512x35.Aniso_k800d02t4.grdt.vec" 
		   // thres = -1  - 61375 seeds    for "tlapse330_T2DkR1k.512x512x35.Aniso_k800d02t4.grdt.vec"
		   // thres =  1  - 15169 seeds    for "tlapse330_T2DkR1k.512x512x35.Aniso_k800d02t4.grdt.vec"
		   // thres =0.2  - 29749 seeds    for "tlapse330_T2DkR1k.512x512x35.Aniso_k800d02t4.grdt.vec"
		   // thres =  0  - 12299 seeds    for "tlapse330_T2DkR1k.512x512x35.grdt.vec"
		   // thres = -1  - 12667 seeds    for "tlapse330_T2DkR1k.512x512x35.grdt.vec"
		 

		   // thres =  0  - 14568 seeds
		   // thres =0.2  - 10368 seeds
           // thres =0.4  -  7713 seeds
		   // thres =0.8  -  5835 seeds    for "Cal20m_T2DkR1k.416x256x35.Aniso_k800d02t4.raw.grdt.vec"  

		   // thres = 1   for "jt72D1_T3DkFR1kDk.512x512x45.Aniso_k800d02t8.raw.3-2-gradt.vec" - 2532 seeds
		   // thres =0.8  for "jt72D1_T3DkFR1kDk.512x512x45.Aniso_k800d02t8.grdt.vec" - 10032 seeds
		   // thres = 3   for "jt72D1_T3DkFR1kDk.512x512x45.grdt.vec" - 3659 seeds (good)
	                    // thres = 0.5 for cube19x29x41
					    // thres = 11  for knight Hierarch #1
					    // thres =  8  for knight Hierarch #2
					    // thres = 11 for knightNS Hierarch #1
					    // thres = 8  for knightNS Hierarch #2
					    // thres = 5 for monster
					    // thres = 10 for tetrahedral
					    // thres = 8 for mushroom
					    // thres = 3 for 767
					    // thres = 10 for human
					    // thres = 18 for TightSegColon
					    // thres = 16 /11/ 6 for herachical cow
					    // thres = 16  for shark
					    // thres = 10  for hammerhead
					    // thres = 14  for dogDilate
					    // thres = 14  for woodman
					    // thres =  8  for teapot
					    // thres =     for temple
					    // thres =  9  for toilet
					    // thres = 50  for cannon
					    // thres =7.5  for kD #1
					    // thres =12   for kD #2
					    // thres = 11  for bolt
					    // thres =3.9for obj
					    // thres = 7  for trout
					    // thres =60for mug

//	     fprintf(fout,"%d %d %d %f %f\n", i, j, k, curv[idx], curv[idx]);
	}

  // only select the local maximum voxel of curv[idx] as skeleton

  DisAway = 2; //5;
  numSeeds = 0;
  for (k = DisAway; k < N-DisAway; k++)
     for (j = DisAway; j < M-DisAway; j++)
        for (i = DisAway; i < L-DisAway; i++) {
	     idx = k*slsz + j*L + i;
	     //consider the six face neighbors
	     iidx = k*slsz + j*L + i-1;
	     if (curv[iidx]> curv[idx]) continue;
	     iidx = k*slsz + j*L + i+1;
	     if (curv[iidx]> curv[idx]) continue;
	     iidx = k*slsz + (j-1)*L + i;
	     if (curv[iidx]> curv[idx]) continue;
	     iidx = k*slsz + (j+1)*L + i;
	     if (curv[iidx]> curv[idx]) continue;
	     iidx = (k-1)*slsz + j*L + i;
	     if (curv[iidx]> curv[idx]) continue;
	     iidx = (k+1)*slsz + j*L + i;
	     if (curv[iidx]> curv[idx]) continue;

		 if (curv[idx] != 0)  {
	          //fprintf(fout,"%d %d %d %f %f\n", i, j, k, curv[idx], curv[idx]);
		    fprintf(fout,"%d %d %d\n", i, j, k);
		    numSeeds++;
		 }
  }

  printf("Num of seeds is %d\n", numSeeds);

  fclose(fout);

  printf("done \n");
  return 0;

}


void PartialDerivative(float *Is, float *Isd, int direc, int L, int M, int N)  {
  // direc = 1: x-direction   2: y-direction   3: z-direction
  long idx;
  int i,k,j;
  long slsz = L*M;
  float lessXsum, moreXsum, lessYsum, moreYsum, lessZsum, moreZsum;
  int disaway = 1;
  for (k = disaway; k < N-disaway; k++)
     for (j = disaway; j < M-disaway; j++)
        for (i = disaway; i < L-disaway; i++) {
	    idx = k*slsz + j*L +i;
	    if (direc == 1) {
	         lessXsum = Is[(k-1)*slsz +(j-1)*L +(i-1)] + Is[(k-1)*slsz +j*L +(i-1)] +Is[(k-1)*slsz +(j+1)*L +(i-1)]
		           +Is[  k*slsz + (j-1)*L +(i-1)] + Is[   k*slsz + j*L +(i-1)] + Is[  k*slsz + (j+1)*L +(i-1)]
			   +Is[(k+1)*slsz +(j-1)*L +(i-1)] + Is[(k+1)*slsz +j*L +(i-1)] +Is[(k+1)*slsz +(j+1)*L +(i-1)];
		 moreXsum = Is[(k-1)*slsz +(j-1)*L +(i+1)] + Is[(k-1)*slsz +j*L +(i+1)] +Is[(k-1)*slsz +(j+1)*L +(i+1)]
		           +Is[  k*slsz + (j-1)*L +(i+1)] + Is[   k*slsz + j*L +(i+1)] + Is[  k*slsz + (j+1)*L +(i+1)]
			   +Is[(k+1)*slsz +(j-1)*L +(i+1)] + Is[(k+1)*slsz +j*L +(i+1)] +Is[(k+1)*slsz +(j+1)*L +(i+1)];
		 Isd[idx] = moreXsum - lessXsum;
	    }
	    else if(direc == 2) {
	         lessYsum = Is[(k-1)*slsz +(j-1)*L +(i-1)] + Is[(k-1)*slsz +(j-1)*L + i] +Is[(k-1)*slsz +(j-1)*L +(i+1)]
		           +Is[  k*slsz + (j-1)*L +(i-1)] + Is[   k*slsz + (j-1)*L +i] + Is[  k*slsz + (j-1)*L +(i+1)]
			   +Is[(k+1)*slsz +(j-1)*L +(i-1)] + Is[(k+1)*slsz +(j-1)*L +i] +Is[(k+1)*slsz +(j-1)*L +(i+1)];
	         moreYsum = Is[(k-1)*slsz +(j+1)*L +(i-1)] + Is[(k-1)*slsz +(j+1)*L + i] +Is[(k-1)*slsz +(j+1)*L +(i+1)]
		           +Is[  k*slsz + (j+1)*L +(i-1)] + Is[   k*slsz + (j+1)*L +i] + Is[  k*slsz + (j+1)*L +(i+1)]
			   +Is[(k+1)*slsz +(j+1)*L +(i-1)] + Is[(k+1)*slsz +(j+1)*L +i] +Is[(k+1)*slsz +(j+1)*L +(i+1)];
		 Isd[idx] = moreYsum - lessYsum;
	    }
	    else {
	         lessZsum = Is[(k-1)*slsz +(j-1)*L +(i-1)] + Is[(k-1)*slsz +(j-1)*L + i] +Is[(k-1)*slsz +(j-1)*L +(i+1)]
		           +Is[ (k-1)*slsz + j *L + (i-1)] + Is[(k-1)*slsz + j*L + i ] + Is[(k-1)*slsz + j*L + (i+1)]
			   +Is[(k-1)*slsz +(j+1)*L +(i-1)] + Is[(k-1)*slsz +(j+1)*L +i] +Is[(k-1)*slsz +(j+1)*L +(i+1)];
	         moreZsum = Is[(k+1)*slsz +(j-1)*L +(i-1)] + Is[(k+1)*slsz +(j-1)*L + i] +Is[(k+1)*slsz +(j-1)*L +(i+1)]
		           +Is[ (k+1)*slsz + j *L + (i-1)] + Is[(k+1)*slsz + j*L + i ] + Is[(k+1)*slsz + j*L + (i+1)]
			   +Is[(k+1)*slsz +(j+1)*L +(i-1)] + Is[(k+1)*slsz +(j+1)*L +i] +Is[(k+1)*slsz +(j+1)*L +(i+1)];
		 Isd[idx] = moreZsum - lessZsum;
	    }
        }
}


void PartialDerivative1(float *Is, float *Isd, int direc, int L, int M, int N)  {
  // direc = 1: x-direction   2: y-direction   3: z-direction
  long idx;
  int i,k,j;
  long slsz = L*M;
  float lessXsum, moreXsum, lessYsum, moreYsum, lessZsum, moreZsum;
  int disaway = 1;
  for (k = disaway; k < N-disaway; k++)
     for (j = disaway; j < M-disaway; j++)
        for (i = disaway; i < L-disaway; i++) {
	    idx = k*slsz + j*L +i;
	    if (direc == 1) {
	         lessXsum = Is[k*slsz + j*L +(i-1)];
		 moreXsum = Is[k*slsz + j*L +(i+1)];
		 Isd[idx] = moreXsum - lessXsum;
	    }
	    else if(direc == 2) {
	         lessYsum = Is[k*slsz + (j-1)*L +i];
	         moreYsum = Is[k*slsz + (j+1)*L +i];
		 Isd[idx] = moreYsum - lessYsum;
	    }
	    else {
	         lessZsum = Is[(k-1)*slsz + j*L + i ];
	         moreZsum = Is[(k+1)*slsz + j*L + i ];
		 Isd[idx] = moreZsum - lessZsum;
	    }
        }
}




void ComputeRotMatrix(float RotateMatrix[3][3], Vector3D v) {
   // Find the matrix that can rotate any vector v to x-axis direction
   float phi, theta;
   Vector3D vector1;
   float RotateMatrix1[3][3];

   theta = atan2(v.z, v.y);
   RotMatrixFromAngle(RotateMatrix1, 0, -theta, 0);
   vector1.x = RotateMatrix1[0][0]*v.x +RotateMatrix1[0][1]*v.y +RotateMatrix1[0][2]*v.z;
   vector1.y = RotateMatrix1[1][0]*v.x +RotateMatrix1[1][1]*v.y +RotateMatrix1[1][2]*v.z;
   vector1.z = RotateMatrix1[2][0]*v.x +RotateMatrix1[2][1]*v.y +RotateMatrix1[2][2]*v.z;
   phi = atan2(vector1.y, vector1.x);
   RotMatrixFromAngle(RotateMatrix, -phi, -theta, 0);
}


void RotMatrixFromAngle(float RMatrix[3][3], float phi, float theta, float psi) {
  // rotation matrix is      R(0,0) R(0,1) R(0,2)
  //	        	     R(1,0) R(1,1) R(1,2)
  //   			     R(2,0) R(2,1) R(2,2)
   RMatrix[0][0] = cos(phi)*cos(psi) - cos(theta)*sin(phi)*sin(psi);
   RMatrix[1][0] = cos(psi)*sin(phi) + cos(phi)*cos(theta)*sin(psi);
   RMatrix[2][0] = sin(psi)*sin(theta);
   RMatrix[0][1] = -(cos(psi)*cos(theta)*sin(phi)) - cos(phi)*sin(psi);
   RMatrix[1][1] = cos(phi)*cos(psi)*cos(theta) - sin(phi)*sin(psi);
   RMatrix[2][1] = cos(psi)*sin(theta);
   RMatrix[0][2] = sin(phi)*sin(theta);
   RMatrix[1][2] = -(cos(phi)*sin(theta));
   RMatrix[2][2] = cos(theta);
}


void Transpose(float MatTransp[3][3], float Mat[3][3])  {
   for(int j=0; j<=2; j++)
        for(int i=0; i<=2; i++)  {
             MatTransp[j][i] = Mat[i][j];
	}
}


void Matrix3Multiply(float Mat[3][3], float Mat1[3][3], float Mat2[3][3]) {
   for(int j=0; j<=2; j++)
        for(int i=0; i<=2; i++)  {
             Mat[j][i] = Mat1[j][0]*Mat2[0][i] +Mat1[j][1]*Mat2[1][i] +Mat1[j][2]*Mat2[2][i];
	}
}

