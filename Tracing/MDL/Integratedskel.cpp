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

// Extract the ridges and valleys feature in 3D vector field
// --- Input: 3D vector field
// --- Output: 3D image-scalar field
// --- Author: Xiaosong Yuan, RPI    modified by xiao liang
// --- Modified Date: 10/6/2005

#if defined(_MSC_VER)
#pragma warning(disable : 4996)
#endif
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

#define FLOAT_SKELETON

#define MIN(x,y) (((x) < (y))?(x):(y))
#define MAX(x,y) (((x) > (y))?(x):(y))
#define SURF 100
#define INTERIOR 200
#define EPS 0.001

int sign(float value) {
   if (value > 0) return 1;
   else if(value < 0) return -1;
        else return 0;
}


int sign1(float value) {
   if (value > 1e-5) return 1;
   else if(value < -1e-5) return -1;
        else return 0;
}

double PartialDerivativeLocal(float *Is, int i,int j,int k, int direc, int L, int M, int N);
double interpolation(float x, float y, float z, int sizx, int sizy, int sizz, float *Iu,float *Iv,float *Iw);



int main (int argc, char *argv[])
{
  ifstream fin;
  FILE *CurveSeedfout;
  FILE *fout;
  Vector3D vecin;
  float *Iu, *Iv, *Iw;
 // float *Iuu, *Iuv, *Iuw, *Ivv, *Ivw, *Iww;
  float *curv;
  unsigned char *f;
  long idx, iidx, slsz, sz;
  int i, j, k;
  int ii, jj, kk;
  int x,y,z;
  int cc;
  int L, M, N;
  double k1,k2;
  float gLength;

  double highCurvatureThreshold;
  long numSeeds;


  Vector3D *seeds;
  bool *FlagOnSkeleton;  //  indicator of skeleton or not 
  
  double vectorMagnitude;



  if (argc < 8)
   {
    printf("Usage: %s <vector file> <xs> <ys> <zs> <vector mag> <seeds file> <out skel> [measureTimeFlag].\n",argv[0]);
    exit(1);
   }

  // ------------------------------------  Open Gradient vector file ----------------------//
  
  fin.open(argv[1]);
  if (!fin)  {
     cerr << "couldn't open " << argv[1] << " for input" << endl;
     return -1;
  }
 //----------------------------------------end open -------------------------------------//

  L = atoi(argv[2]);
  M = atoi(argv[3]);
  N = atoi(argv[4]);

  vectorMagnitude = atof(argv[5]);

  slsz = L*M;        // slice size
  sz = slsz*N;

  //highCurvatureThreshold = atof(argv[5]);

  //highCurvatureThreshold = 20;
 
  

  if ((CurveSeedfout= fopen(argv[6],"w")) == NULL)  // arg[6]
  {
    printf("Cannot open %s for writing\n",argv[6]);  
    exit(1);
  }

  if ((fout= fopen(argv[7],"w")) == NULL)  // arg[7]
  {
    printf("Cannot open %s for writing\n",argv[6]);  
    exit(1);
  }

  Iu = new float[L*M*N];
  Iv = new float[L*M*N];
  Iw = new float[L*M*N];
  seeds = new Vector3D[20000000];  


  curv = new float[L*M*N];
  f = new  unsigned char [L*M*N];
  FlagOnSkeleton = new bool [L*M*N];

  for(idx=0; idx<sz; idx++)   {  //Initialize to zeros
	  f[idx] =0;
	  Iu[idx] =0;
	  Iv[idx] =0;
	  Iw[idx] =0;
	  curv[idx] =0;
	  FlagOnSkeleton[idx] = 0;
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
    
  printf("hhhhhh %ld", idx);
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

          
// ------------------------------------------------------------------------------------------------//
// ------------------------------------------------------------------------------------------------// 
// We use mean curvature to coompute the higher curvature points

// Dr. Xheiaosong's method,  to compute the princple curvature; see the original C++ code;
// Dr. xiao think that it is not necessary to compute the principle curvature, just mean curvature is ok!
     
// ------------------------------------------------------------------------------------------------//

 double PIuu, PIvv, PIww, PIuv, PIuw, PIvw;  // temp varibles to save the second derivatives
 double gx,gy,gz;  // temp variables to save the Ix^2,Iy^2,Iz^2;
 int DisAway = 2;

 double PrinC;
  for (k = DisAway; k < N-DisAway; k++)
     for (j = DisAway; j < M-DisAway; j++)
        for (i = DisAway; i < L-DisAway; i++) 
		{
	       idx = k*slsz + j*L + i;

	       if (f[idx] !=2) continue;

		   PIuu = PartialDerivativeLocal(Iu, i,j,k, 1, L, M, N);
           PIvv = PartialDerivativeLocal(Iv, i,j,k, 2, L, M, N);
           PIww = PartialDerivativeLocal(Iw, i,j,k, 3, L, M, N);
           PIuv = PartialDerivativeLocal(Iu, i,j,k, 2, L, M, N);
           PIuw = PartialDerivativeLocal(Iu, i,j,k, 3, L, M, N);
           PIvw = PartialDerivativeLocal(Iv, i,j,k, 3, L, M, N);

		   gx=Iu[idx]*Iu[idx]; gy=Iv[idx]*Iv[idx]; gz=Iw[idx]*Iw[idx];

           gLength = (float) sqrt(gx +gy +gz);// amplitude of gradient 
		   if(gLength<0.0001) continue;
		
		   //k1=PIuu+PIvv+PIww+PIuu*(gy +gz) +PIvv*(gx+gz)+PIww*(gx+gy);
		   k1=gx*(PIvv+PIww)+gy*(PIuu+PIww)+gz*(PIuu+PIvv);
		   k1 -= 2*Iu[idx] * Iv[idx] * PIuv + 2*Iv[idx] * Iw[idx] * PIvw + 2*Iu[idx] * Iw[idx] * PIuw;
		   k1 = -k1/((2*pow(gLength,(float)1.5)));  // this is mean curvature

		   //k2 = PIuu*PIvv + PIuu*PIww + PIvv*PIww-PIuu*PIuu-PIvv*PIvv - PIww*PIww;
		   k2 =Iu[idx]*Iv[idx]*PIuw*PIvw+Iu[idx]*Iw[idx]*PIuv*PIvw+Iv[idx]*Iw[idx]*PIuv*PIuw;
		   k2 -=Iu[idx]*Iw[idx]*PIuw*PIvv+Iv[idx]*Iw[idx]*PIuu*PIvw+Iu[idx]*Iv[idx]*PIuv*PIww;

		   k2=k2*2;
           k2+=gx*PIvv*PIww+gy*PIuu*PIww+gz*PIuu*PIvv;
		   k2-=gx*PIvw*PIvw+gy*PIuw*PIuw+gz*PIuv*PIuv;

		   k2/=(gx+gy+gz)*(gx+gy+gz);//this is Guassian curvature
           if (k1*k1-k2>0)
		   PrinC=k1+sqrt(k1*k1-k2); // this is Principle Curvature;
		   else 
           PrinC=k1;
       	   curv[idx] =(float) (PrinC);    
		   //if (abs(k2)<0.01)curv[idx] =(float) (k1); 
		  // curv[idx]=(float)k2;
	}// end of computing mean curvature for every points.

// ------------------------------------------------------------------------------------------------//
// ----------use the mean value of the mean curvature as the threshold --------------------------- //

   double meanCurvature=0;
 
   for (k = DisAway; k < N-DisAway; k++)
     for (j = DisAway; j < M-DisAway; j++)
        for (i = DisAway; i < L-DisAway; i++) 
		{
	       idx = k*slsz + j*L +i;
		   meanCurvature += curv[idx];
		   
		}

   meanCurvature /=(double)(L*M*N); 

   meanCurvature = (meanCurvature<0 ? 0:meanCurvature);
   highCurvatureThreshold =(meanCurvature>10 ? 10:meanCurvature);  //the threshold is a value between [0,10].

   //highCurvatureThreshold=0;
   //highCurvatureThreshold =10;//;

   printf("mean =%f The estimated  good threshold of highCurvatureThreshold %f \n",meanCurvature, highCurvatureThreshold );

// --------------------------------- end estimation  -------------------------------------------- //


   
// ---------------------------------compute seeds with local maiximal curvature  ---------------------//
  DisAway = 2; //5;
  numSeeds = 0;
  // only select the local maximum voxel of curv[idx] as skeleton
  for (k = DisAway; k < N-DisAway; k++)
     for (j = DisAway; j < M-DisAway; j++)
        for (i = DisAway; i < L-DisAway; i++) {
	     idx = k*slsz + j*L + i;
	     //consider the six face neighbors

         if (curv[idx]< highCurvatureThreshold ) curv[idx]= 0;   // thresholding compared with threshold 

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

		 if (curv[idx] != 0 && numSeeds <L*M*N/1000)  {
	          //fprintf(fout,"%d %d %d %f %f\n", i, j, k, curv[idx], curv[idx]);
		    seeds[numSeeds].x =(float)i;
            seeds[numSeeds].y =(float)j;
			seeds[numSeeds].z =(float)k;
		    fprintf(CurveSeedfout,"%d %d %d\n", i, j, k);
		    numSeeds++;
		 }
		
  }
//------------------------------------------end compute seeds with higher curvature ----------------------// 
  printf("Num of seeds with higher curvature  is %ld\n", numSeeds);

//-------------------------------Normalizing gradient vector---------------------------------------------//
  float Length;

  idx = 0;
  for (k = 1; k < N-1; k++)
     for (j = 1; j < M-1; j++)
        for (i = 1; i < L-1; i++) 
		{
	     idx = k*slsz + j*L +i;
         Length = (float) sqrt(Iu[idx]*Iu[idx]+Iv[idx]*Iv[idx]+Iw[idx]*Iw[idx]);
         Iu[idx] /= (float) (Length+0.0001);
         Iv[idx] /= (float) (Length+0.0001);
         Iw[idx] /= (float) (Length+0.0001);
  }

//------------------------------------------Start compute seeds with critical points  ----------------------//
  //printf(" I am here  FFFF %d",  idx);

  int numBoundSeeds = numSeeds;
  float Ngrid = 8;
  int SemiNgrid= (int)Ngrid/2;

  double totalVecLength;

  
  double  OutForce;
  double divx,divy,divz,div;
  long NumCritPoints =0;

  for (k = 1; k < N-1; k++)
  //for (k = 1; k < -1; k++)  //TEST: skip all critical pts
   for (j = 1; j < M-1; j++)
     for (i = 1; i < L-1; i++) 
	  {
	    idx = k*slsz + j*L +i;
		if (sqrt(Iu[idx]*Iu[idx]+Iv[idx]*Iv[idx]+Iw[idx]*Iw[idx])<0.1)continue;
	    totalVecLength = 0;
		divx = 0;  divy = 0;  divz = 0;
		// to check whether the position is near to boundaries, compute divergence in x,y,z respectively
	    if (f[k*slsz + j*L  +i] == 0) 
			continue;	else  {
				idx = k*slsz + j*L +i;
				totalVecLength +=sqrt(Iu[idx]*Iu[idx]+Iv[idx]*Iv[idx]+Iw[idx]*Iw[idx]);
				divx -= Iu[idx];
				divy -= Iv[idx];
				divz -= Iw[idx];
			}
	    if (f[k*slsz + j*L  +(i+1)] == 0) 
			continue;   else  {
				idx = k*slsz + j*L  +(i+1);
				totalVecLength +=sqrt(Iu[idx]*Iu[idx]+Iv[idx]*Iv[idx]+Iw[idx]*Iw[idx]);
				divx += Iu[idx];
				divy -= Iv[idx];
				divz -= Iw[idx];
			}
	    if (f[k*slsz +(j+1)*L + i] == 0) 
			continue;   else  {
				idx=k*slsz +(j+1)*L + i;
				totalVecLength +=sqrt(Iu[idx]*Iu[idx]+Iv[idx]*Iv[idx]+Iw[idx]*Iw[idx]);
				divx -= Iu[idx];
				divy += Iv[idx];
				divz -= Iw[idx];
			}
	    if (f[k*slsz +(j+1)*L + (i+1)] == 0) 
			continue;   else  {
				idx=k*slsz +(j+1)*L + (i+1);
				totalVecLength += sqrt(Iu[idx]*Iu[idx]+Iv[idx]*Iv[idx]+Iw[idx]*Iw[idx]);
				divx += Iu[idx];
				divy += Iv[idx];
				divz -= Iw[idx];
			}
	    if (f[(k+1)*slsz + j*L  +i] == 0) 
			continue;   else  {
				idx= (k+1)*slsz + j*L  +i;
				totalVecLength += sqrt(Iu[idx]*Iu[idx]+Iv[idx]*Iv[idx]+Iw[idx]*Iw[idx]);
				divx -= Iu[idx];
				divy -= Iv[idx];
				divz += Iw[idx];
			}
	    if (f[(k+1)*slsz + j*L  +(i+1)] == 0) 
			continue;   else  {
				idx=(k+1)*slsz + j*L  +(i+1);
				totalVecLength += sqrt(Iu[idx]*Iu[idx]+Iv[idx]*Iv[idx]+Iw[idx]*Iw[idx]);
				divx += Iu[idx];
				divy -= Iv[idx];
				divz += Iw[idx];
			}
	    if (f[(k+1)*slsz +(j+1)*L + i] == 0) 
			continue;   else  {
				idx=(k+1)*slsz +(j+1)*L + i;
				totalVecLength += sqrt(Iu[idx]*Iu[idx]+Iv[idx]*Iv[idx]+Iw[idx]*Iw[idx]);
				divx -= Iu[idx];
				divy += Iv[idx];
				divz += Iw[idx];
			}
	    if (f[(k+1)*slsz +(j+1)*L + (i+1)] == 0) 
			continue;   else  {
				idx=(k+1)*slsz +(j+1)*L + i;
				totalVecLength += sqrt(Iu[idx]*Iu[idx]+Iv[idx]*Iv[idx]+Iw[idx]*Iw[idx]);
				divx += Iu[idx];
				divy += Iv[idx];
				divz += Iw[idx];
			}

		if (totalVecLength < 4)   {  // skip the zero vector areas
			continue;
		}

		div = divx + divy + divz;
		if (div > -1) {    //skip if the cube possibly contain a repelling point (divergence is too big)
			              // 1;  (good for phantom) 
						  // 0;  (better) 
						  // -1; (not always good)
			//printf("div=%f (%d)  ", div, k);
			continue;
		}
		
     	
	    //for (kk = 0; kk<Ngrid; kk++)
		//  for (jj = 0; jj<Ngrid; jj++)
		 //   for (ii = 0; ii<Ngrid; ii++) 

		for (kk = SemiNgrid; kk<Ngrid; kk++)
		  for (jj = SemiNgrid; jj<Ngrid; jj++)
		    for (ii = SemiNgrid; ii<Ngrid; ii++)
			{
              
		        OutForce=interpolation(i+ii/Ngrid, j+jj/Ngrid, k+kk/Ngrid, L, M, N, Iu,Iv,Iw);
				
			    if( OutForce < vectorMagnitude && NumCritPoints < 2*numBoundSeeds) {   //<0.15
				
				seeds[numSeeds].x= i + ii/Ngrid;
				seeds[numSeeds].y= j + jj/Ngrid;
				seeds[numSeeds].z= k + kk/Ngrid;
				numSeeds++;
				NumCritPoints++;
				//printf("%f ,% f ",OutForce,numSeeds);
			}//end if
	    }// end for
    // printf(" I am here  2222");
  

    // numSeeds++;
    // NumCritPoints++;
}// end outer for 

//-------------------------------------------------------------------------------------------//
 
//printf("Number of critical points is: %d, and number of seeds %d\n", NumCritPoints, numSeeds);


//-------------------------------------------------------------------------------------------//

//------------------------------------Line Path algorithm------------------------------------//

 // seeds[0..numBoundSeeds-1] = boundary seeds, seeds[numBoundSeeds.. numSeeds-1] = critical points
 long idxSeeds = numSeeds-1;
 Vector3D Startpos;

 int Dvoxel = 2;  // 3 is good for many results

 int streamSteps = 0;

 while (idxSeeds >= 0) 
 {

	Startpos.x=seeds[idxSeeds].x; Startpos.y=seeds[idxSeeds].y; Startpos.z=seeds[idxSeeds].z;
	idx = (int)Startpos.z * slsz + (int)Startpos.y *L + (int)Startpos.x;

	//Check whether the high curv points are within D-voxel distance to the existing skeleton
	  
	if (idxSeeds < numBoundSeeds)   Dvoxel= 2;  //4; // 10 <- Found too big on April 12, 2006, which cause less spines
	                                             // 3 is good for tlapse330, 6 is good for Trach
	int FlagWithin = 0;
	if (idxSeeds < numBoundSeeds) 
	{
		for (kk = -Dvoxel+1; kk <= Dvoxel-1; kk++)
		   for (jj = -Dvoxel+1; jj <= Dvoxel-1; jj++)
		      for (ii = -Dvoxel+1; ii <= Dvoxel-1; ii++) 
			  {
				  iidx = idx + kk*slsz + jj*L + ii;
				  if (iidx < 0 || iidx >= sz)  continue;
		          if(FlagOnSkeleton[iidx] == 1)  FlagWithin = 1;
		       }//end for
		if(FlagWithin == 1) {
			idxSeeds--;
			continue;
		}// end if 
	}// end if

	// being able to not show critical point in the streamlines
	if (idxSeeds >= numBoundSeeds)
		{
		   #ifdef FLOAT_SKELETON
		     fprintf(fout,"%f %f %f %d\n", Startpos.x, Startpos.y, Startpos.z, 1);
		    #else
		     fprintf(fout,"%d %d %d %d %d\n",(int)Startpos.x,(int)Startpos.y,(int)Startpos.z, idxSeeds, idxSeeds);
		   #endif
		}
	else {
		   #ifdef FLOAT_SKELETON
	             fprintf(fout,"%f %f %f %d\n", Startpos.x, Startpos.y, Startpos.z, 1);
		    #else
		     fprintf(fout,"%d %d %d %d %d\n",(int)Startpos.x,(int)Startpos.y,(int)Startpos.z, idxSeeds, idxSeeds);
		   #endif
	    }

	FlagOnSkeleton[idx] = 1;

  /*
	while(streamSteps < 1) //4000  
	  {    // < 4000
		rk2(Startpos.x, Startpos.y, Startpos.z, sizeX, sizeY, sizeZ, float(0.8), force, &Nextpos);   //0.2, 0.8, 2     float() added by xiao liang
		streamSteps++;
		Startpos.x = Nextpos.x;
		Startpos.y = Nextpos.y;
		Startpos.z = Nextpos.z;

		idx = (int)Nextpos.z *slsz + (int)Nextpos.y *sizeX + (int)Nextpos.x;
        if (FlagOnSkeleton[idx] != 1) {
			#ifdef FLOAT_SKELETON
	                  fprintf(fout,"%f %f %f %d\n", Nextpos.x, Nextpos.y, Nextpos.z, 1);
			//printf("%f %f %f %d\n", Nextpos.x, Nextpos.y, Nextpos.z, 1);
			 #else
			  fprintf(fout,"%d %d %d %d %d\n",(int)Nextpos.x,(int)Nextpos.y,(int)Nextpos.z, idxSeeds, idxSeeds);
			#endif
			FlagOnSkeleton[idx] = 1;
		} // end if 
	} // end while
   */
	streamSteps = 0;
	idxSeeds--;
 } // end outer while



// ----------------------------------------release memory and file close by xiao liang-------------------//
  
 fclose(CurveSeedfout);
 fclose(fout);
 delete []Iu;
 delete []Iv;
 delete []Iw;
 delete []curv;
 delete []f;
 delete []FlagOnSkeleton;
 return 0;

}


double PartialDerivativeLocal(float *Is, int i,int j,int k, int direc, int L, int M, int N)  {
  // direc = 1: x-direction   2: y-direction   3: z-direction
  long idx;
  long slsz = L*M;
  float lessXsum, moreXsum, lessYsum, moreYsum, lessZsum, moreZsum;
  double Isd;


	    idx = k*slsz + j*L +i; // get the location
        // boundary treatment;
		if (i ==L-1) i= L-2;
		if (j ==M-1) j= M-2;
		if (k ==N-1) k= N-2; 
	    if (direc == 1) {
	         lessXsum = Is[k*slsz + j*L +(i-1)];
		 moreXsum = Is[k*slsz + j*L +(i+1)];
		 Isd = moreXsum - lessXsum;
	    }
	    else if(direc == 2) {  // if it is y-direction
	         lessYsum = Is[k*slsz + (j-1)*L +i];
	         moreYsum = Is[k*slsz + (j+1)*L +i];
		 Isd = moreYsum - lessYsum;
	    }
	    else {  
	         lessZsum = Is[(k-1)*slsz + j*L + i ];
	         moreZsum = Is[(k+1)*slsz + j*L + i ];
		 Isd = moreZsum - lessZsum;
	    }

      return Isd; 
}



double interpolation(float x, float y, float z, int sizx, int sizy, int sizz, float *Iu,float *Iv,float *Iw)
    {
	
	//-------by xiao  new implementation
		
	float alpha, beta, gamma;
	///Vector3D forceInt;
	long slsz;
    int Intx,Inty,Intz;
	Intx=int(x);
	Inty=int(y);
	Intz=int(z);
	alpha = x- Intx;   
	beta = y- Inty;
	gamma = z-Intz;
	slsz=sizy*sizx;

	float tempIu=0, tempIv=0,tempIw=0;
	
	/*
	long NeigberPosition;   
	NeigberPosition=Intz*slsz + Inty*sizx + Intx;
    
	float a[8];
	a[0] = (1-alpha)*(1-beta)*(1-gamma);
	a[1] = (1-alpha)*(1-beta)*gamma;
	a[2] = (1-alpha)*beta*(1-gamma);
	a[3] = alpha*(1-beta)*(1-gamma);
	a[4] = alpha*(1-beta)*gamma;
	a[5] = alpha*beta*(1-gamma);
	a[6] = (1-alpha)*beta*gamma;
	a[7] = alpha*beta*gamma;

	tempIu=Iu[NeigberPosition]*a[0]
			+Iu[NeigberPosition+slsz]*a[1]
			+Iu[NeigberPosition+sizx]*a[2]
			+Iu[NeigberPosition+1]*a[3]
			+Iu[NeigberPosition+slsz+1]*a[4]
			+Iu[NeigberPosition+sizx+1]*a[5]
			+Iu[NeigberPosition+slsz+sizx]*a[6]
			+Iu[NeigberPosition+slsz+sizx+1]*a[7];

	
	tempIv=Iv[NeigberPosition]*a[0]
			+Iv[NeigberPosition+slsz]*a[1]
			+Iv[NeigberPosition+sizx]*a[2]
			+Iv[NeigberPosition+1]*a[3]
			+Iv[NeigberPosition+slsz+1]*a[4]
			+Iv[NeigberPosition+sizx+1]*a[5]
			+Iv[NeigberPosition+slsz+sizx]*a[6]
			+Iv[NeigberPosition+slsz+sizx+1]*a[7];	

    tempIw=Iw[NeigberPosition]*a[0]
			+Iw[NeigberPosition+slsz]*a[1]
			+Iw[NeigberPosition+sizx]*a[2]
			+Iw[NeigberPosition+1]*a[3]
			+Iw[NeigberPosition+slsz+1]*a[4]
			+Iw[NeigberPosition+sizx+1]*a[5]
			+Iw[NeigberPosition+slsz+sizx]*a[6]
			+Iw[NeigberPosition+slsz+sizx+1]*a[7];	
	*/

	//printf("  %f  ", Iu[Intz*slsz + Inty*sizx + Intx]);

	tempIu=Iu[Intz*slsz + Inty*sizx + Intx]*(1-alpha)*(1-beta)*(1-gamma)
			+Iu[(Intz+1)*slsz + Inty*sizx + Intx]*(1-alpha)*(1-beta)*gamma
			+Iu[Intz*slsz + (Inty+1)*sizx + Intx]*(1-alpha)*beta*(1-gamma)
			+Iu[Intz*slsz + Inty*sizx + (Intx+1)]*alpha*(1-beta)*(1-gamma)
			+Iu[(Intz+1)*slsz + Inty*sizx + (Intx+1)]*alpha*(1-beta)*gamma
			+Iu[Intz*slsz + (Inty+1)*sizx + (Intx+1)]*alpha*beta*(1-gamma)
			+Iu[(Intz+1)*slsz + (Inty+1)*sizx + Intx]*(1-alpha)*beta*gamma
			+Iu[(Intz+1)*slsz + (Inty+1)*sizx + (Intx+1)]*(alpha*beta*gamma);

	tempIu=Iv[Intz*slsz + Inty*sizx + Intx]*(1-alpha)*(1-beta)*(1-gamma)
			+Iv[(Intz+1)*slsz + Inty*sizx + Intx]*(1-alpha)*(1-beta)*gamma
			+Iv[Intz*slsz + (Inty+1)*sizx + Intx]*(1-alpha)*beta*(1-gamma)
			+Iv[Intz*slsz + Inty*sizx + (Intx+1)]*alpha*(1-beta)*(1-gamma)
			+Iv[(Intz+1)*slsz + Inty*sizx + (Intx+1)]*alpha*(1-beta)*gamma
			+Iv[Intz*slsz + (Inty+1)*sizx + (Intx+1)]*alpha*beta*(1-gamma)
			+Iv[(Intz+1)*slsz + (Inty+1)*sizx + Intx]*(1-alpha)*beta*gamma
			+Iv[(Intz+1)*slsz + (Inty+1)*sizx + (Intx+1)]*alpha*beta*gamma;

	tempIw=Iw[Intz*slsz + Inty*sizx + Intx]*(1-alpha)*(1-beta)*(1-gamma)
			+Iw[(Intz+1)*slsz + Inty*sizx + Intx]*(1-alpha)*(1-beta)*gamma
			+Iw[Intz*slsz + (Inty+1)*sizx + Intx]*(1-alpha)*beta*(1-gamma)
			+Iw[Intz*slsz + Inty*sizx + (Intx+1)]*alpha*(1-beta)*(1-gamma)
			+Iw[(Intz+1)*slsz + Inty*sizx + (Intx+1)]*alpha*(1-beta)*gamma
			+Iw[Intz*slsz + (Inty+1)*sizx + (Intx+1)]*alpha*beta*(1-gamma)
			+Iw[(Intz+1)*slsz + (Inty+1)*sizx + Intx]*(1-alpha)*beta*gamma
			+Iw[(Intz+1)*slsz + (Inty+1)*sizx + (Intx+1)]*alpha*beta*gamma;
	
   
	    //return  (sqrt(tempIu*tempIu+tempIw*tempIw+tempIv*tempIv));
	      return abs(tempIu)+abs(tempIv)+abs(tempIw);  // because sqrt() is time consuming;
		
		
    }
