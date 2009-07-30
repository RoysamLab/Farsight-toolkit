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


double PartialDerivative3(float *Is, int i,int j,int k, int direc, int L, int M, int N);



int main (int argc, char *argv[])
{
  ifstream fin;
  FILE *fout;
  char *infilename;
  Vector3D vecin;
  float *Iu, *Iv, *Iw;
 // float *Iuu, *Iuv, *Iuw, *Ivv, *Ivw, *Iww;
  float *curv;
  int *f;
  long idx, iidx, slsz, sz;
  int i, j, k;
  int ii, jj, kk;
  int x,y,z;
  int cc;
  int L, M, N;
  double k1;
  float gLength;

  double highCurvatureThreshold;
  int numSeeds;

  if (argc < 6)
  {
    printf("Usage: %s <Input vectors> <xs> <ys> <zs>  <Output scalars>.\n", argv[0]);
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

  //highCurvatureThreshold = atof(argv[5]);

  //highCurvatureThreshold = 20;
 
  

  if ((fout = fopen(argv[5],"w")) == NULL)  // arg[6]
  {
    printf("Cannot open %s for writing\n",argv[5]);  // by xiao  ,,, argc[5]
    exit(1);
  }

  Iu = new float[L*M*N];
  Iv = new float[L*M*N];
  Iw = new float[L*M*N];

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



 double PIuu,PIvv, PIww, PIuv,PIuw,PIvw;
 double gx,gy,gz;
  int DisAway = 2;
  for (k = DisAway; k < N-DisAway; k++)
     for (j = DisAway; j < M-DisAway; j++)
        for (i = DisAway; i < L-DisAway; i++) 
		{
	     idx = k*slsz + j*L + i;
	     if (f[idx] !=2) continue;


		   PIuu = PartialDerivative3(Iu, i,j,k, 1, L, M, N);
           PIvv = PartialDerivative3(Iv, i,j,k, 2, L, M, N);
           PIww = PartialDerivative3(Iw, i,j,k, 3, L, M, N);
           PIuv = PartialDerivative3(Iu, i,j,k, 2, L, M, N);
           PIuw = PartialDerivative3(Iu, i,j,k, 3, L, M, N);
           PIvw = PartialDerivative3(Iv, i,j,k, 3, L, M, N);

          
// /////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////////////////////// 
// Dr. Xheiaosong's method,  to compute the princple curvature; see the original C++ code;
// Dr. xiao think that it is not necessary to compute the principle curvature, just mean curvature is ok!
     
//////////////////////////////////////////////////////////////////////////////////////////////////////


        
		 gx=Iu[idx]*Iu[idx];gy=Iv[idx]*Iv[idx];gz=Iw[idx]*Iw[idx];
         gLength = (float) sqrt(gx +gy +gz);// amplitude of gradient 
		
		 k1=PIuu+PIvv+PIww+PIuu*(gy +gz) +PIvv*(gx+gz)+PIww*(gx+gy);
		 k1 -= 2*Iu[idx] * Iv[idx] * PIuv + 2*Iv[idx] * Iw[idx] * PIvw + 2*Iu[idx] * Iw[idx] * PIuw;
		 k1 = -k1/((2*pow(gLength,(float)1.5))+0.0001);  // this is mean curvature

		 curv[idx] =(float) (k1);  
        
	
	}

//////////////////////////////////////////////////////////////////////////////
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

  highCurvatureThreshold =(meanCurvature>10 ? 10:meanCurvature);
  //highCurvatureThreshold =10;

  printf("mean =%f The estimated  good threshold of highCurvatureThreshold %f \n",meanCurvature, highCurvatureThreshold );


  DisAway = 2; //5;
  numSeeds = 0;
  // only select the local maximum voxel of curv[idx] as skeleton
  for (k = DisAway; k < N-DisAway; k++)
     for (j = DisAway; j < M-DisAway; j++)
        for (i = DisAway; i < L-DisAway; i++) {
	     idx = k*slsz + j*L + i;
	     //consider the six face neighbors

         if (curv[idx]< highCurvatureThreshold ) curv[idx]= 0;   // thresholding to 

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

// release memory and file close by xiao liang
  fclose(fout);
  delete []curv;
  delete []f;
  return 0;

}


double PartialDerivative3(float *Is, int i,int j,int k, int direc, int L, int M, int N)  {
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

