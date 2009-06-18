// ----
// ----  Extract the thinness metric based skeleton. 
// ----  Uses Toriwaki and Saito's DT algorithm. 
// ----  
// ----  Implementation by : Xiaosong, RPI
// ----
// ----  Input : Binary 3D volume with sizes. 
// ----  Output: volume file with distance transform, 
// ----          DT is zero for all object voxels; DT is scaled to 0-255.
// ----
// $Id:  2008/05/02 $

#include "MinSpanTree.h"
#include "distTransform.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <utility>

using namespace std;

#define MIN(x,y) (((x) < (y))?(x):(y))
#define MAX(x,y) (((x) > (y))?(x):(y))

// Arguments of function:
// Input as a volume with non-zero as objects;
// Output to the same volume as distance trasnform 

void distTransform(unsigned char *f, int L, int M, int N) 
{
  
  int i,j,k,n;
  float *buff , df, db, d, w;
  long *fDist;
  long idx, slsz, sz;

  slsz = L*M;		// slice size
  sz = slsz*N;

  fDist = new long[L*M*N];

  for (idx = 0; idx < slsz*N; idx++) {
      if (f[idx] > 0)   fDist[idx] = 0;
	  else fDist[idx] = 5000;
  }

  int maxdim = MAX(L,M);
  maxdim = MAX(maxdim,N);
  buff = new float[maxdim+10];

  // Using Algorithm 3 from Appendix 

  // Step 1  forward scan
  
  for (k = 0; k < N; k++)
    for (j = 0; j < M; j++)
    {
       df = L;
       for (i = 0; i < L; i++)
       {
         idx = k*slsz + j*L + i;
         if (fDist[idx] !=0)
           df = df + 1;
         else
           df = 0;
         fDist[idx] = df*df;
       }
     }
 
  //  Step 1 backward scan
  
  for (k = 0; k < N; k++)
    for (j = 0; j < M; j++)
    {
      db = L;
      for (i = L-1; i >=0; i--)
      {
        idx = k*slsz + j*L + i;
        if (fDist[idx] !=0)
          db = db + 1;
        else
          db = 0; 
        fDist[idx] = MIN(fDist[idx],db*db);
      }
    }

  // Step 2

  for (k = 0; k < N; k++)
    for (i = 0; i < L; i++)
    {
      for (j =0; j < M; j++)
        buff[j] = fDist[k*slsz + j*L +i];
    
      for (j = 0; j < M; j++)
      {
        d = buff[j];
        if (d != 0)
        {
          int rmax, rstart, rend;
          rmax = (int) floor(sqrt(d)) + 1;
          rstart = MIN(rmax, (j-1));
          rend = MIN(rmax, (M-j));
          for (n = -rstart; n < rend; n++)
          {
              if (j+n >= 0 && j+n < M)
              {
                w = buff[j+n] + n*n;
                if (w < d)  d = w;
              }
          }
        }
        idx = k*slsz + j*L +i;
        fDist[idx] = d;
      }
    }

  // Step 3
  for (j = 0; j < M; j++)
    for (i = 0; i < L; i++)
    {
      for (k =0; k < N; k++)
        buff[k] = fDist[k*slsz + j*L +i];
    
      for (k = 0; k < N; k++)
      {
        d = buff[k];
        if (d != 0)
        {
          int rmax, rstart, rend;
          rmax = (int) floor(sqrt(d)) + 1;
          rstart = MIN(rmax, (k-1));
          rend = MIN(rmax, (N-k));
          for (n = -rstart; n < rend; n++)
          {
              if (k+n >= 0 && k+n < N)
              {
                w = buff[k+n] + n*n;
                if (w < d)  d = w;
              }
          }
        }
        idx = k*slsz + j*L +i;
        fDist[idx] = d;
      }
    }

  for (idx = 0; idx < slsz*N; idx++) {
      fDist[idx] = sqrt(float(fDist[idx]));
  }


  double dMax = 0;
  for(idx=0; idx<sz; idx++)   {  // Scale the dist result to 255
	  if (fDist[idx] > dMax)   dMax = fDist[idx];
  }
  for(idx=0; idx<sz; idx++)   {  // Scale the dist result to 255
	  //f[idx] = fDist[idx] * 255/ dMax;
	  f[idx] = fDist[idx];
  }

  return;

}






