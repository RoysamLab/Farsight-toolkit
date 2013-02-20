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
// --- Author: Xiaosong Yuan, improved and optimized by xiao liang, RPI
// --- Modified Date: 03/Sep/2009

#if defined(_MSC_VER)
#pragma warning(disable : 4996)
#endif
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <algorithm>
#include "robustness.h"

using std::cerr;
using std::cout;
using std::endl;

struct  Vector3D
{
  float x;
  float y;
  float z;
};

struct  VoxelPosition
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

int sign(float value)
  {
  if (value > 0)
    {
    return 1;
    }
  else if(value < 0)
    {
    return -1;
    }
  else
    {
    return 0;
    }
  }


int sign1(float value)
  {
  if (value > 1e-5)
    {
    return 1;
    }
  else if(value < -1e-5)
    {
    return -1;
    }
  else
    {
    return 0;
    }
  }

float veclength(Vector3D vecin)
  {
  return sqrt(vecin.x * vecin.x + vecin.y * vecin.y +vecin.z * vecin.z);
  }

//function prototypes
void PartialDerivative(float *Is, float *Isd, int direc, int L, int M, int N);
void PartialDerivative1(float *Is, float *Isd, int direc, int L, int M, int N);
void ComputeRotMatrix(float RotateMatrix[3][3], Vector3D v);
void RotMatrixFromAngle(float RMatrix[3][3], float phi, float theta, float psi);
void Transpose(float MatTransp[3][3], float Mat[3][3]);
void Matrix3Multiply(float Mat[3][3], float Mat1[3][3], float Mat2[3][3]);
void rk2(float x, float y, float z, int sizx, int sizy, int sizz, float steps,
         Vector3D *Force_ini, VoxelPosition *nextPos);
Vector3D interpolation(float x, float y, float z, int sizx, int sizy, int sizz,
                       Vector3D *forcevec);
void RotMatrixFromAngle(float RMatrix[3][3], float cosphi, float sinphi,
                        float costheta, float sintheta, float cospsi,
                        float sinpsi);

int main (int argc, char *argv[])
{
  ifstream fin;
  FILE *CurveSeedfout;
  FILE *fout;
  Vector3D vecin,gradient0;
  float *Iu, *Iv, *Iw;
  float *Iuu, *Iuv, *Iuw, *Ivv, *Ivw, *Iww;
  float *curv;
  unsigned char *f;
  long idx, iidx, slsz, sz;
  int i, j, k;
  int ii, jj, kk;
  int x,y,z;
  int cc;
  int L, M, N;
  //double k1,k2;
  double k1;
  float gLength;
  
  float RotMatrix[3][3];
  float RotMatrixTransp[3][3];
  float Hessian[3][3];
  float HessianPrime[3][3];
  float HessianPrime1[3][3];

  double highCurvatureThreshold;
  long numSeeds;


  Vector3D *seeds;
  bool *FlagOnSkeleton;  //  indicator of skeleton or not 
  
  double vectorMagnitude;



  if (argc < 8)
    {
    cerr << "Usage: " << argv[0] << " <vector file> <xs> <ys> <zs> "
         << "<vector mag> <seeds file> <out skel> [measureTimeFlag]" << endl;
     return 1;
    }

  //----------------------  Open Gradient vector file ----------------------//
  fin.open(argv[1]);
  if (!fin)
    {
    cerr << "couldn't open " << argv[1] << " for input" << endl;
    return -1;
    }
  //----------------------------- end open ---------------------------------//

  L = atoi(argv[2]);
  M = atoi(argv[3]);
  N = atoi(argv[4]);

  vectorMagnitude = atof(argv[5]);

  slsz = L*M;        // slice size
  sz = slsz*N;

  if ((CurveSeedfout= fopen(argv[6],"w")) == NULL)  // arg[6]
    {
    cerr << "Cannot open " << argv[6] << " for writing" << endl;
    return 1;
    }

  if ((fout= fopen(argv[7],"w")) == NULL)  // arg[7]
    {
    cerr << "Cannot open " << argv[7] << " for writing" << endl;
    return 1;
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

  seeds = new Vector3D[20000000];  


  curv = new float[L*M*N];
  f = new  unsigned char [L*M*N];
  FlagOnSkeleton = new bool [L*M*N];

  for(idx=0; idx<sz; idx++)
    {
    //Initialize to zeros
    f[idx] =0;
    Iu[idx] =0;
    Iv[idx] =0;
    Iw[idx] =0;
    curv[idx] =0;
    FlagOnSkeleton[idx] = 0;
    }

  fin >> x >> y >> z >> vecin.x >> vecin.y >> vecin.z;
  cout << "Reading input" << endl;
  while (!fin.eof() ) 
    {
    idx = z*slsz + y*L +x;
    Iu[idx] = vecin.x;
    Iv[idx] = vecin.y;
    Iw[idx] = vecin.z;
    f[idx]=2;
    fin >> x >> y >> z >> vecin.x >> vecin.y >> vecin.z;
    }
  fin.close();
    
  // make surface point f[]=1
  for (k = 1; k < N-1; k++)
    {
    for (j = 1; j < M-1; j++)
      {
      for (i = 1; i < L-1; i++)
        {
        idx = k*slsz + j*L +i;
        if(f[idx]==2)
          {
          cc = 0;
          for (kk=-1; kk<=1; kk++)
            {
            for (jj=-1; jj<=1; jj++)
              {
              for (ii=-1; ii<=1; ii++)
                {
                iidx = (k+kk)*slsz + (j+jj)*L +(i+ii);
                if (f[iidx]==0)
                  {
                  cc++;
                  }
                }
              }
            }
          if (cc>=1)
            {
            f[idx] = 1;
            }
          }
        }
      }
    }
          
  //------------------------------------------------------------------------//
  //Compute iso-gray surface principle surface curvature    
  //------------------------------------------------------------------------//
  cout << "Computing curvature" << endl;
  PartialDerivative1(Iu, Iuu, 1, L, M, N);
  PartialDerivative1(Iv, Ivv, 2, L, M, N);
  PartialDerivative1(Iw, Iww, 3, L, M, N);
  PartialDerivative1(Iu, Iuv, 2, L, M, N);
  PartialDerivative1(Iu, Iuw, 3, L, M, N);
  PartialDerivative1(Iv, Ivw, 3, L, M, N);
/*
  		FILE *fileout;
		if ((fileout = fopen("out.deriv","w")) != NULL)
		{
			for (int k = 0; k < N; k++) {
				for (int j = 0; j < M; j++) {
					for (int i = 0; i < L; i++) 
					{
						long idx = k*L*M + j*L + i;
						//fprintf(fileout, "%4.4f %4.4f %4.4f\n", Iu[idx], Iv[idx], Iw[idx]);
						fprintf(fileout, "%4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n", Iuu[idx], Ivv[idx], Iww[idx], Iuv[idx], Iuw[idx], Ivw[idx]);
					} 
				} 
			} //end triple for loops
			fclose(fileout);
		}
*/
  //double maxCurvature = 0;
  int DisAway = 2;
  for (k = DisAway; k < N-DisAway; k++)
    {
    for (j = DisAway; j < M-DisAway; j++)
      {
      for (i = DisAway; i < L-DisAway; i++)
        {
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
        k1 = 0.5*(HessianPrime[1][1]+HessianPrime[2][2]) +
             0.5*sqrt(pow((Hessian[1][1]-Hessian[2][2]),2) +
             4*Hessian[1][2]*Hessian[2][1]);
		// the second eignvale is not used 
        /*k2 = 0.5*(HessianPrime[1][1]+HessianPrime[2][2]) -
             0.5*sqrt(pow((Hessian[1][1]-Hessian[2][2]),2) +
             4*Hessian[1][2]*Hessian[2][1]);
	    */ 
        gLength =
          sqrt(gradient0.x*gradient0.x + gradient0.y*gradient0.y +
               gradient0.z*gradient0.z);

        if(gLength !=0)
          {
          curv[idx] = (float) -(k1)/gLength;
          }
		else
		{
		   curv[idx]=0;
		}
		/*
        else
          {
          //set large to be included in skeleton
          curv[idx]=9999;
          }

        //case 1: get neg maximum
        if (curv[idx]>10000)
          {
          curv[idx]=10000;
          }
		  
		
        if (curv[idx]< -1 )
          {
          // first phase: threshold the curvature value (larger->fewer)
          curv[idx]= 0; 
          }
		  */  // the above is xiaosong's method

        }
      }
    }

  //-----------------------end of computing curvature-----------------------//
  //----------------------------release memory------------------------------// 
  delete []Iuu;
  delete []Ivv;
  delete []Iww;
  delete []Iuv;
  delete []Iuw;
  delete []Ivw;

  		//DEBUGGING OUTPUT
		FILE *fout2;
		if ((fout2 = fopen("out.curv","w")) != NULL)
		{
			for (int c = 0; c < N; c++) 
			{
				for (int b = 0; b < M; b++) 
				{
					for (int a = 0; a < L; a++) 
					{
						int idx3 = c*L*M + b*L + a;
						fprintf(fout2, "%d %d %d %f\n", a, b, c, curv[idx3]);
					}
				}
			} //end triple for loops
			fclose(fout2);
		}

  //------------------------------------------------------------------------//

  //-------use the mean value of the mean curvature as the threshold--------//

  double meanCurvature=0;
  //double Tmean = 0; 
  long InterVox =1; 
  for (k = DisAway; k < N-DisAway; k++)
    {
    for (j = DisAway; j < M-DisAway; j++)
      {
      for (i = DisAway; i < L-DisAway; i++) 
        {
        idx = k*slsz + j*L +i;
		if (f[idx]!=0)
		{
         meanCurvature += curv[idx];
		 InterVox++;
        }// end if 
      }// end for
	}
  }
   
  //meanCurvature /=(double)(L*M*N); 
  meanCurvature /=(double)(InterVox); 
 
  meanCurvature = (meanCurvature<0 ? 0:meanCurvature);

  //the threshold is a value between [0,10].
  highCurvatureThreshold = (meanCurvature > 10 ? 10 : meanCurvature);
  //cout << "Trimmed mean is " << Tmean << endl;
  cout << "mean = " << meanCurvature << "." << endl;
  cout << "The estimated good threshold of highCurvatureThreshold "
       << highCurvatureThreshold << endl;

  //----------------------------end estimation------------------------------//

  //-------------compute seeds with local maiximal curvature----------------//
  DisAway = 2; //5;
  numSeeds = 0;
  // only select the local maximum voxel of curv[idx] as skeleton
  for (k = DisAway; k < N-DisAway; k++)
    {
    for (j = DisAway; j < M-DisAway; j++)
      {
      for (i = DisAway; i < L-DisAway; i++)
        {
        idx = k*slsz + j*L + i;

        //consider the six face neighbors

        if (curv[idx]< highCurvatureThreshold )
          {
          // thresholding compared with threshold 
          curv[idx]= 0;
          }

        iidx = k*slsz + j*L + i-1;
        if (curv[iidx]> curv[idx])
          {
          continue;
          }
        iidx = k*slsz + j*L + i+1;
        if (curv[iidx]> curv[idx])
          {
          continue;
          }
        iidx = k*slsz + (j-1)*L + i;
        if (curv[iidx]> curv[idx])
          {
          continue;
          }
        iidx = k*slsz + (j+1)*L + i;
        if (curv[iidx]> curv[idx])
          {
          continue;
          }
        iidx = (k-1)*slsz + j*L + i;
        if (curv[iidx]> curv[idx])
          {
          continue;
          }
        iidx = (k+1)*slsz + j*L + i;
        if (curv[iidx]> curv[idx])
          {
          continue;
          }

        if (curv[idx] != 0)
          {
          seeds[numSeeds].x =(float)i;
          seeds[numSeeds].y =(float)j;
          seeds[numSeeds].z =(float)k;
          fprintf(CurveSeedfout,"%d %d %d\n", i, j, k);
          numSeeds++;
          }
        }
      }
    }
  //----------------end compute seeds with higher curvature-----------------//
  cout << "Num of seeds with higher curvature is " << numSeeds << endl;

  //-----------------------Normalizing force vector-------------------------//
  float Length;
  float totalVecLength;
  Vector3D *force;
  int sizeX=L;
  int sizeY=M;
  int sizeZ=N;
  force = new Vector3D[sizeX*sizeY*sizeZ];
 

  for(idx=0; idx<sz; idx++)
    {
    //Initialize to zeros
    //fc[idx]=0;
    force[idx].x = 0;
    force[idx].y = 0;
    force[idx].z = 0;
    FlagOnSkeleton[idx] = 0;
    }

  //-------------------------Generating force vector------------------------//
  idx = 0;
  for (k = 1; k < N-1; k++)
    {
    for (j = 1; j < M-1; j++)
      {
      for (i = 1; i < L-1; i++) 
        {
        idx = k*slsz + j*L +i;
        Length = (float) sqrt(Iu[idx]*Iu[idx]+Iv[idx]*Iv[idx]+Iw[idx]*Iw[idx]);
        Iu[idx] /= (float) (Length+0.0001);
        Iv[idx] /= (float) (Length+0.0001);
        Iw[idx] /= (float) (Length+0.0001);
        force[idx].x = Iu[idx];
        force[idx].y = Iv[idx];
        force[idx].z = Iw[idx]; 
        // this is just a very simple methods, you can use Xiaosong's method,
        //i.e POW()  and also you can use GDF method
        }
      }
    }
  cout << "Generating force vector!\n" << endl;

  //-----------------------Delete Gradient Vector Memeory---------------------//
  delete []Iu;
  delete []Iv;
  delete []Iw;

  // find all critical points: use points with small length of vector
  float Ngrid = 10;
  int semiNgrid=0;
  Vector3D  OutForce;
  double  divx,divy,divz,div;
  long    NumCritPoints =0;

  for (k = 1; k < sizeZ-1; k++)
    {
    //this loop can take a long time to complete.
    //provide some feedback to the user on how we're doing
    if( ( ( (float)k / float(sizeZ) ) >= 0.25 ) &&
        ( ( (float)(k-1) / float(sizeZ) ) < 0.25 ) )
      {
      cout << "25% complete" << endl;
      }
    if( ( ( (float)k / float(sizeZ) ) >= 0.50 ) &&
        ( ( (float)(k-1) / float(sizeZ) ) < 0.50 ) )
      {
      cout << "50% complete" << endl;
      }
    if( ( ( (float)k / float(sizeZ) ) >= 0.75 ) &&
        ( ( (float)(k-1) / float(sizeZ) ) < 0.75 ) )
      {
      cout << "75% complete" << endl;
      }
    for (j = 1; j < sizeY-1; j++)
      {
      for (i = 1; i < sizeX-1; i++)
        {
        idx = k*slsz + j*sizeX +i;
        totalVecLength = 0;
        divx = 0;  divy = 0;  divz = 0;
        // to check whether the position is near to boundaries, compute
        //divergence in x,y,z respectively
        if (f[k*slsz + j*sizeX  +i] == 0) 
          {
          continue;
          }
        else
          {
          totalVecLength += veclength(force[k*slsz + j*sizeX  +i]);
          divx -= force[k*slsz + j*sizeX  +i].x;
          divy -= force[k*slsz + j*sizeX  +i].y;
          divz -= force[k*slsz + j*sizeX  +i].z;
          }
        if (f[k*slsz + j*sizeX  +(i+1)] == 0) 
          {
          continue;
          }
        else
          {
          totalVecLength += veclength(force[k*slsz + j*sizeX  +(i+1)]);
          divx += force[k*slsz + j*sizeX  +(i+1)].x;
          divy -= force[k*slsz + j*sizeX  +(i+1)].y;
          divz -= force[k*slsz + j*sizeX  +(i+1)].z;
          }
        if (f[k*slsz +(j+1)*sizeX + i] == 0) 
          {
          continue;
          }
        else
          {
          totalVecLength += veclength(force[k*slsz + (j+1)*sizeX  +i]);
          divx -= force[k*slsz + (j+1)*sizeX  +i].x;
          divy += force[k*slsz + (j+1)*sizeX  +i].y;
          divz -= force[k*slsz + (j+1)*sizeX  +i].z;
          }
        if (f[k*slsz +(j+1)*sizeX + (i+1)] == 0) 
          {
          continue;
          }
        else
          {
          totalVecLength += veclength(force[k*slsz + (j+1)*sizeX  +(i+1)]);
          divx += force[k*slsz + (j+1)*sizeX  +(i+1)].x;
          divy += force[k*slsz + (j+1)*sizeX  +(i+1)].y;
          divz -= force[k*slsz + (j+1)*sizeX  +(i+1)].z;
          }
        if (f[(k+1)*slsz + j*sizeX  +i] == 0) 
          {
          continue;
          }
        else
          {
          totalVecLength += veclength(force[(k+1)*slsz + j*sizeX  +i]);
          divx -= force[(k+1)*slsz + j*sizeX  +i].x;
          divy -= force[(k+1)*slsz + j*sizeX  +i].y;
          divz += force[(k+1)*slsz + j*sizeX  +i].z;
          }
        if (f[(k+1)*slsz + j*sizeX  +(i+1)] == 0) 
          {
          continue;
          }
        else
          {
          totalVecLength += veclength(force[(k+1)*slsz + j*sizeX  +(i+1)]);
          divx += force[(k+1)*slsz + j*sizeX  +(i+1)].x;
          divy -= force[(k+1)*slsz + j*sizeX  +(i+1)].y;
          divz += force[(k+1)*slsz + j*sizeX  +(i+1)].z;
          }
        if (f[(k+1)*slsz +(j+1)*sizeX + i] == 0) 
          {
          continue;
          }
        else
          {
          totalVecLength += veclength(force[(k+1)*slsz + (j+1)*sizeX  +i]);
          divx -= force[(k+1)*slsz + (j+1)*sizeX  +i].x;
          divy += force[(k+1)*slsz + (j+1)*sizeX  +i].y;
          divz += force[(k+1)*slsz + (j+1)*sizeX  +i].z;
          }
        if (f[(k+1)*slsz +(j+1)*sizeX + (i+1)] == 0) 
          {
          continue;
          }
        else
          {
          totalVecLength += veclength(force[(k+1)*slsz + (j+1)*sizeX  +(i+1)]);
          divx += force[(k+1)*slsz + (j+1)*sizeX  +(i+1)].x;
          divy += force[(k+1)*slsz + (j+1)*sizeX  +(i+1)].y;
          divz += force[(k+1)*slsz + (j+1)*sizeX  +(i+1)].z;
          }

        if (totalVecLength < 4)
          {
          // skip the zero vector areas
          continue;
          }

        div = divx + divy + divz;
        if (div > -1)
          {
          //skip if the cube possibly contain a repelling point
          //(divergence is too big)
          continue;
          }
        for (kk = semiNgrid; kk<Ngrid; kk++)
          {
          for (jj = semiNgrid; jj<Ngrid; jj++)
            {
            for (ii = semiNgrid; ii<Ngrid; ii++)
              {
              OutForce = interpolation(
                i+ii/Ngrid, j+jj/Ngrid, k+kk/Ngrid, sizeX, sizeY, sizeZ, force);
              if(veclength(OutForce) < vectorMagnitude)
                {
                seeds[numSeeds].x= i + ii/Ngrid;
                seeds[numSeeds].y= j + jj/Ngrid;
                seeds[numSeeds].z= k + kk/Ngrid;
                numSeeds++;
                NumCritPoints++;
                }
              }
            }
          }
        }
      }
    }

 cout << "Number of critical points is: " << NumCritPoints << endl;
 cout << "Number of seeds is: " << numSeeds << endl;

 //int numBoundSeeds = numSeeds;  //Changed Isaac 2/04/2010
 // I figured this out by looking at the old skel_streamline.cpp
 int numBoundSeeds = (numSeeds-NumCritPoints);
 //fprintf(fout,"%d %d %d %f %f\n", 1, 1, 1, -1.0, -1.0);

 // seeds[0..numBoundSeeds-1] = boundary seeds, seeds[numBoundSeeds.. numSeeds-1] = critical points
 int idxSeeds = numSeeds-1;
 VoxelPosition Startpos, Nextpos;
 int streamSteps=0;

 while (idxSeeds >= 0) 
 {

  Startpos.x=seeds[idxSeeds].x; Startpos.y=seeds[idxSeeds].y; Startpos.z=seeds[idxSeeds].z;
  idx = (int)Startpos.z * slsz + (int)Startpos.y *sizeX + (int)Startpos.x;

  //Check whether the high curv points are within D-voxel distance to the existing skeleton
  int Dvoxel = 2;    // 3 is good for many results
  if (idxSeeds < numBoundSeeds)   Dvoxel= 2;  //4; // 10 <- Found too big on April 12, 2006, which cause less spines
                                               // 3 is good for tlapse330, 6 is good for Trach
  int FlagWithin = 0;
  if (idxSeeds < numBoundSeeds) 
  {
    for (kk = -Dvoxel+1; kk <= Dvoxel-1; kk++)
       for (jj = -Dvoxel+1; jj <= Dvoxel-1; jj++)
          for (ii = -Dvoxel+1; ii <= Dvoxel-1; ii++) 
        {
          iidx = idx + kk*slsz + jj*sizeX + ii;
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

   //--------------------------------------Line path algorithm-------------------------------------------------------------------//
  while(streamSteps < 4000) //4000  
    {    // < 4000
    rk2(Startpos.x, Startpos.y, Startpos.z, sizeX, sizeY, sizeZ, float(0.8), force, &Nextpos);   //0.2, 0.8, 2     float() added by xiao liang
    streamSteps++;
    Startpos.x = Nextpos.x;
    Startpos.y = Nextpos.y;
    Startpos.z = Nextpos.z;
	//If at edge of image get out of loop
	if(Nextpos.x >= sizeX-1 || Nextpos.y >= sizeY-1 || Nextpos.z >= sizeZ-1)
		break;
	if(Nextpos.x < 0 || Nextpos.y < 0 || Nextpos.z < 0)
		break;

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
   
  streamSteps = 0;
  idxSeeds--;
 } // end outer while 
//----------------------------------------end of line path algorithm --------------------------------------------------------------//

   fclose(fout);

   delete []f;
   delete []force;
   delete []FlagOnSkeleton;
   delete []seeds;
   delete []curv;

   cout << "End" << endl;
   return 0;
}




    
Vector3D interpolation(float x, float y, float z, int sizx, int sizy, int sizz, Vector3D  *forcevec)
    {
  
  //-------by xiao  new implementation
  Vector3D forceInt;  
  float alpha, beta, gamma;
  long slsz;
  int Intx,Inty,Intz;
  x= (x <0 ? 0:x);
  y= (y <0 ? 0:y);
  z= (z <0 ? 0:z);
  x= (x>=sizx-1 ? sizx-2:x);
  y= (y>=sizy-1 ? sizy-2:y);
  z= (z>=sizz-1 ? sizz-2:z);

  Intx=int(x);
  Inty=int(y);
  Intz=int(z);
  alpha = x- Intx;   
  beta = y- Inty;
  gamma = z-Intz;
  slsz=sizy*sizx;
  float a[8];// for interpolation coefficients
  a[0] = (1-alpha)*(1-beta)*(1-gamma);
  a[1] = (1-alpha)*(1-beta)*gamma;
  a[2] = (1-alpha)*beta*(1-gamma);
  a[3] =  alpha*(1-beta)*(1-gamma); 
  a[4] =  alpha*(1-beta)*gamma;
  a[5] =  alpha*beta*(1-gamma);
  a[6] =  (1-alpha)*beta*gamma;
  a[7] =  (alpha*beta*gamma);

  long Nei[8]; // considering N8 neighborhood in 3D image
  
  Nei[0] = Intz*slsz + Inty*sizx + Intx;
  Nei[1] =  Nei[0] +slsz;
  Nei[2] =  Nei[0] +slsz+sizx;
  Nei[3] =  Nei[0] +1;
  Nei[4] =  Nei[0] +slsz+1;
  Nei[5] =  Nei[0] +sizx+1;
  Nei[6] =  Nei[0] +sizx;
  Nei[7] =  Nei[0] +slsz+sizx+1;

   //------------------------------------ compute interpolation ---------------------------------------//
  forceInt.x=forcevec[ Nei[0]].x*a[0]
      +forcevec[ Nei[1]].x*a[1]
      +forcevec[ Nei[2]].x*a[2]
      +forcevec[ Nei[3]].x*a[3]
      +forcevec[ Nei[4]].x*a[4]
      +forcevec[ Nei[5]].x*a[5]
      +forcevec[ Nei[6]].x*a[6]
      +forcevec[ Nei[7]].x*a[7];


  forceInt.y=forcevec[ Nei[0]].y*a[0]
      +forcevec[ Nei[1]].y*a[1]
      +forcevec[ Nei[2]].y*a[2]
      +forcevec[ Nei[3]].y*a[3]
      +forcevec[ Nei[4]].y*a[4]
      +forcevec[ Nei[5]].y*a[5]
      +forcevec[ Nei[6]].y*a[6]
      +forcevec[ Nei[7]].y*a[7];

  forceInt.z=forcevec[ Nei[0]].z*a[0]
      +forcevec[ Nei[1]].z*a[1]
      +forcevec[ Nei[2]].z*a[2]
      +forcevec[ Nei[3]].z*a[3]
      +forcevec[ Nei[4]].z*a[4]
      +forcevec[ Nei[5]].z*a[5]
      +forcevec[ Nei[6]].z*a[6]
      +forcevec[ Nei[7]].z*a[7];

  return(forceInt);
    
  /*float alpha, beta, gamma;
  Vector3D forceInt;
  long slsz;
    int Intx,Inty,Intz;
  Intx=int(x);
  Inty=int(y);
  Intz=int(z);
  alpha = x- Intx;   
  beta = y- Inty;
  gamma = z-Intz;
  slsz=sizy*sizx;

  forceInt.x=forcevec[Intz*slsz + Inty*sizx + Intx].x*(1-alpha)*(1-beta)*(1-gamma)
      +forcevec[(Intz+1)*slsz + Inty*sizx + Intx].x*(1-alpha)*(1-beta)*gamma
      +forcevec[Intz*slsz + (Inty+1)*sizx + Intx].x*(1-alpha)*beta*(1-gamma)
      +forcevec[Intz*slsz + Inty*sizx + (Intx+1)].x*alpha*(1-beta)*(1-gamma)
      +forcevec[(Intz+1)*slsz + Inty*sizx + (Intx+1)].x*alpha*(1-beta)*gamma
      +forcevec[Intz*slsz + (Inty+1)*sizx + (Intx+1)].x*alpha*beta*(1-gamma)
      +forcevec[(Intz+1)*slsz + (Inty+1)*sizx + Intx].x*(1-alpha)*beta*gamma
      +forcevec[(Intz+1)*slsz + (Inty+1)*sizx + (Intx+1)].x*(alpha*beta*gamma);

  forceInt.y=forcevec[Intz*slsz + Inty*sizx + Intx].y*(1-alpha)*(1-beta)*(1-gamma)
      +forcevec[(Intz+1)*slsz + Inty*sizx + Intx].y*(1-alpha)*(1-beta)*gamma
      +forcevec[Intz*slsz + (Inty+1)*sizx + Intx].y*(1-alpha)*beta*(1-gamma)
      +forcevec[Intz*slsz + Inty*sizx + (Intx+1)].y*alpha*(1-beta)*(1-gamma)
      +forcevec[(Intz+1)*slsz + Inty*sizx + (Intx+1)].y*alpha*(1-beta)*gamma
      +forcevec[Intz*slsz + (Inty+1)*sizx + (Intx+1)].y*alpha*beta*(1-gamma)
      +forcevec[(Intz+1)*slsz + (Inty+1)*sizx + Intx].y*(1-alpha)*beta*gamma
      +forcevec[(Intz+1)*slsz + (Inty+1)*sizx + (Intx+1)].y*alpha*beta*gamma;

  forceInt.z=forcevec[Intz*slsz + Inty*sizx + Intx].z*(1-alpha)*(1-beta)*(1-gamma)
      +forcevec[(Intz+1)*slsz + Inty*sizx + Intx].z*(1-alpha)*(1-beta)*gamma
      +forcevec[Intz*slsz + (Inty+1)*sizx + Intx].z*(1-alpha)*beta*(1-gamma)
      +forcevec[Intz*slsz + Inty*sizx + (Intx+1)].z*alpha*(1-beta)*(1-gamma)
      +forcevec[(Intz+1)*slsz + Inty*sizx + (Intx+1)].z*alpha*(1-beta)*gamma
      +forcevec[Intz*slsz + (Inty+1)*sizx + (Intx+1)].z*alpha*beta*(1-gamma)
      +forcevec[(Intz+1)*slsz + (Inty+1)*sizx + Intx].z*(1-alpha)*beta*gamma
      +forcevec[(Intz+1)*slsz + (Inty+1)*sizx + (Intx+1)].z*alpha*beta*gamma;

  return(forceInt);
  */

    }


void rk2(float x, float y, float z, int sizx, int sizy, int sizz, float steps, Vector3D  *Force_ini, VoxelPosition *nextPos)
   {
  
    Vector3D OutForce;
  OutForce=interpolation(x,y,z,sizx,sizy,sizz,Force_ini);
  nextPos->x = x + OutForce.x * steps;
  nextPos->y = y + OutForce.y * steps;
  nextPos->z = z + OutForce.z * steps;

  /*x = x + OutForce.x * steps;
  y = y + OutForce.y * steps;
  z = z + OutForce.z * steps;
  nextPos->x = x;
  nextPos->y = y;
  nextPos->z = z; 
    */

   }



// ----------------------------------- non optimized computing partial Derivativel ---------------------------------//

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



// ------------------------------------ Optimization for partialDerivel ----------------------------------------------//
// ------------------------------------ ---------by Xiao Liang--------- ----------------------------------------------// 
/*
void PartialDerivative1(float *Is, float *Isd, int direc, int L, int M, int N)  {
  // direc = 1: x-direction   2: y-direction   3: z-direction
  long idx;
  int i,k,j;
  long slsz = L*M;
  int disaway = 1;
  for (k = disaway; k < N-disaway; k++)
     for (j = disaway; j < M-disaway; j++)
        for (i = disaway; i < L-disaway; i++) 
    {
    // forward difference;
    idx = k*slsz + j*L +i;
      if (direc == 1) 
     Isd[idx] = Is[idx+1] - Is[idx-1];
      else if(direc == 2) 
       Isd[idx] =  Is[idx+L] - Is[idx-L];
      else   
       Isd[idx] = Is[idx+slsz] - Is[idx-slsz];
        }// end for
}
*/


//--------------------------------compute rotation of a matrix --------------------------------//
/*
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
*/
// -----------------------------------------------------------------------------------------------//
//----------------------------- non optimized Rotation of Matrix --------------------------------//
/*
void RotMatrixFromAngle(float RMatrix[3][3], float phi, float theta, float psi) {
  // rotation matrix is      R(0,0) R(0,1) R(0,2)
  //                 R(1,0) R(1,1) R(1,2)
  //             R(2,0) R(2,1) R(2,2)
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
*/

//--------------------------------optimized Rotation of Matrix --------------------------------------------------------------------//
//-----------------------------------  this code written by Xiao Liang ------------------------------------------------------------//
//--------------------------------------the idea is void computing cos(),sin() so many times --------------------------------------//


void ComputeRotMatrix(float RotateMatrix[3][3], Vector3D v) {
   // Find the matrix that can rotate any vector v to x-axis direction
   float phi, theta;
   Vector3D vector1;
   float RotateMatrix1[3][3];
   float cosphi,sinphi,cospsi,sinpsi,costheta,sintheta;

   theta = atan2(v.z, v.y);
   //RotMatrixFromAngle(RotateMatrix1, 0, -theta, 0);
   cosphi=1;sinphi=0;costheta=cos(theta);sintheta=sin(theta);cospsi=1;sinpsi=0;
   RotMatrixFromAngle(RotateMatrix1,cosphi,sinphi,costheta,sintheta,cospsi,sinpsi);
   vector1.x = RotateMatrix1[0][0]*v.x +RotateMatrix1[0][1]*v.y +RotateMatrix1[0][2]*v.z;
   vector1.y = RotateMatrix1[1][0]*v.x +RotateMatrix1[1][1]*v.y +RotateMatrix1[1][2]*v.z;
   vector1.z = RotateMatrix1[2][0]*v.x +RotateMatrix1[2][1]*v.y +RotateMatrix1[2][2]*v.z;
   phi = atan2(vector1.y, vector1.x);
   
   cosphi=cos(-phi);sinphi=sin(-phi);
   // costheta=costheta; //it is not needed to reseted since them havenot changed
   sintheta=-sintheta;
   //cospsi=1;sinpsi=0;  //it is not needed to reseted since them havenot changed 
   RotMatrixFromAngle(RotateMatrix,cosphi,sinphi,costheta,sintheta,cospsi,sinpsi);
   //RotMatrixFromAngle(RotateMatrix, -phi, -theta, 0);
}

void RotMatrixFromAngle(float RMatrix[3][3], float cosphi, float sinphi, float costheta,float sintheta, float cospsi, float sinpsi) {
  // rotation matrix is      R(0,0) R(0,1) R(0,2)
  //                 R(1,0) R(1,1) R(1,2)
  //             R(2,0) R(2,1) R(2,2)
   RMatrix[0][0] = cosphi*cospsi - costheta*sinphi*sinpsi;
   RMatrix[1][0] = cospsi*sinphi + cosphi*costheta*sinpsi;
   RMatrix[2][0] = sinpsi*sintheta;
   RMatrix[0][1] = -cospsi*costheta*sinphi - cosphi*sinpsi;
   RMatrix[1][1] = cosphi*cospsi*costheta - sinphi*sinpsi;
   RMatrix[2][1] = cospsi*sintheta;
   RMatrix[0][2] = sinphi*sintheta;
   RMatrix[1][2] = -cosphi*sintheta;
   RMatrix[2][2] = costheta;
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

/*
double TrimmedMean(float a[], int elements, double r)
{
   printf("hello!");
   int i =0;
   //int elements = sizeof(a) / sizeof(a[0]); 
   printf("%d",elements);
   float *b = new float[elements];
   //std::sort(a, a + elements);
   for (int i = 0; i < elements; ++i) 
   {// printf("%f ", a[i]);
     b[i] = a[i];
   }

   int g = 0;
   if(r>=0 && r<0.5)
     g = int (r*elements);
   else g = 0;
   
   double Tmean =0;
   for (i=g;i<elements-g;i++)
     Tmean +=(double)b[i];
   Tmean /=(double)(elements-2*g);  
   return Tmean;
  
}

*/
