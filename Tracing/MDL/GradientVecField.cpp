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

// ----
// ----  Compute the gradient vector field of any 3D density map
// ----
// ----  Created by : Xiaosong Yuan, RPI
// ----  Date: Nov. 11, 2005
// ----  Input : 3D volume density map with any sizes.
// ----  Output: ASCII file with vector 3 components for all object voxels
// ----
#if defined(_MSC_VER)
#pragma warning(disable : 4996)
#endif
#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <math.h>
#include <string.h>
using std::cerr;
using std::cout;
using std::endl;

#define MIN(x,y) (((x) < (y))?(x):(y))
#define MAX(x,y) (((x) > (y))?(x):(y))

#define DATATYPEIN unsigned char

#define SURF 100
#define INTERIOR 200


struct  VoxelPosition
{
	short x;
	short y;
	short z;
};

struct  Vector
{
	double xd;   // For large datasets, using double will cause memory shortage
	double yd;
	double zd;
};

int sign(float value) {
   if (value > 1e-5) return 1;
   else if(value < -1e-5) return -1;
        else return 0;
}


int main(int argc, char *argv[])
{
	
  FILE *filein;
  FILE *fileout;
  DATATYPEIN *volin;
  unsigned char *fc;
  Vector *gradVec;

  int sizeX,sizeY,sizeZ;         // Sizes in x,y,z dimensions
  int MidX, MidY, MidZ;
  int i,j,k;
  int d1, d2;
  long idx, iidx, sls, sz;
  long iidx1, iidx2;
  int measureTime = 0;
  int flagBound;
  int numBound=0;
  int border;
  double s;
  double kernelWeight[3][3];

  if (argc < 6)
  {
    cerr << "Usage: " << argv[0] << " <volfile> <xs> <ys> <zs> <outfile> "
         << "[measureTimeFlag]" << endl;
    return 1;
  }


  if ((filein = fopen(argv[1],"rb")) == NULL)
  {
    cerr << "Cannot open " << argv[1] << endl;
    return 1;
  }

  sizeX = atoi(argv[2]);
  sizeY = atoi(argv[3]);
  sizeZ = atoi(argv[4]);
  MidX = (sizeX-1)/2;
  MidY = (sizeY-1)/2;
  MidZ = (sizeZ-1)/2;
  //ThresDiv = atof(argv[6]);
  if (argc > 6)
     measureTime = 1;

  fc = (unsigned char*)malloc(sizeX*sizeY*sizeZ*sizeof(unsigned char));
  volin = (DATATYPEIN*)malloc(sizeX*sizeY*sizeZ*sizeof(DATATYPEIN));
  gradVec = (Vector *)malloc(sizeX*sizeY*sizeZ*sizeof(Vector));
  
  sls = sizeX*sizeY;		// slice size
  sz = sls*sizeZ;

  if ( fread(volin, sizeof(DATATYPEIN), sz, filein) < (unsigned long)sz)
  {
    cerr << "File size is not the same as volume size" << endl;
    return 1;
  }

  if ((fileout = fopen(argv[5],"w")) == NULL)
  {
    cerr << "Cannot open " << argv[5] << " for writing" << endl;
    return 1;
  }


  for (idx = 0; idx < sls*sizeZ; idx++)  {
	  if (volin[idx] > 0) {
	       fc[idx] = INTERIOR;
	  }
	  else
		   fc[idx] = 0;
  }

  fclose(filein);

  // Determine surface voxels

  border = 1;
  for (k = border; k < sizeZ-border; k++)
    for (j = border; j < sizeY-border; j++)
      for (i = border; i < sizeX-border; i++)
      {
	     flagBound = 0;
         idx = k*sls + j*sizeX + i;

		// CASE 1: treat the inner layer
	    if (fc[idx] == 0) continue;
		//consider six face neighbors, if anyone is zero, it is a boundary voxel
		iidx = k*sls + j*sizeX + i-1;
		if (fc[iidx] == 0) flagBound = 1;
		iidx = k*sls + j*sizeX + i+1;
		if (fc[iidx] == 0) flagBound = 1;
		iidx = k*sls + (j-1)*sizeX + i;
		if (fc[iidx] == 0) flagBound = 1;
		iidx = k*sls + (j+1)*sizeX + i;
		if (fc[iidx] == 0) flagBound = 1;
		iidx = (k-1)*sls + j*sizeX + i;
		if (fc[iidx] == 0) flagBound = 1;
		iidx = (k+1)*sls + j*sizeX + i;
		if (fc[iidx] == 0) flagBound = 1;

		if (flagBound == 1) {
			numBound++;
			fc[idx] = SURF;
		 }
	  }

  // define positive half kernel of derivative (Sobel kernel)
  kernelWeight[0][0] = 1.67; kernelWeight[0][1] = 5.80; kernelWeight[0][2] = 1.67;
  kernelWeight[1][0] = 5.80; kernelWeight[1][1] = 20.1; kernelWeight[1][2] = 5.80;
  kernelWeight[2][0] = 1.67; kernelWeight[2][1] = 5.80; kernelWeight[2][2] = 1.67;

  for (k = border; k < sizeZ-border; k++)
    for (j = border; j < sizeY-border; j++)
      for (i = border; i < sizeX-border; i++)
	  {
		  idx = k*sls + j*sizeX + i;
		  if (fc[idx] == 0)  {
 			  gradVec[idx].xd = 0;
			  gradVec[idx].yd = 0;
			  gradVec[idx].zd = 0;
			  continue; //consider interior only
		  }

		  gradVec[idx].xd = 0;
		  gradVec[idx].yd = 0;
		  gradVec[idx].zd = 0;
		  for (d2 = -1; d2 <= 1; d2++)
			  for (d1 = -1; d1 <= 1; d1++)  {
				  iidx1 = (k+d2) *sls + (j+d1) *sizeX + i-1;
				  iidx2 = (k+d2) *sls + (j+d1) *sizeX + i+1;
				  gradVec[idx].xd = gradVec[idx].xd + kernelWeight[d2+1][d1+1]*(volin[iidx2] - volin[iidx1]);

				  iidx1 = (k+d2) *sls + (j-1) *sizeX + (i+d1);
				  iidx2 = (k+d2) *sls + (j+1) *sizeX + (i+d1);
				  gradVec[idx].yd = gradVec[idx].yd + kernelWeight[d2+1][d1+1]*(volin[iidx2] - volin[iidx1]);

				  iidx1 = (k-1) *sls + (j+d2) *sizeX + (i+d1);
				  iidx2 = (k+1) *sls + (j+d2) *sizeX + (i+d1);
				  gradVec[idx].zd = gradVec[idx].zd + kernelWeight[d2+1][d1+1]*(volin[iidx2] - volin[iidx1]);

		  }

  }


	//print force vectors
    s = 0.002; //scale on outputs  0.002
	for (k = 0; k < sizeZ; k++)
	  for (j = 0; j < sizeY; j++)
		for (i = 0; i < sizeX; i++) {
			idx = k*sls + j*sizeX + i;
			
			//if (fc[idx]==0 || fc[idx]==SURF)   continue;  // not output
			if (fc[idx]==0)   continue;  // not output
			fprintf(fileout, "%d %d %d %f %f %f\n", i, j, k, gradVec[idx].xd*s, gradVec[idx].yd*s, gradVec[idx].zd*s);
	}


   fclose(fileout);

   free(fc);// = (unsigned char*)malloc(sizeX*sizeY*sizeZ*sizeof(unsigned char));
   free(volin);// = (DATATYPEIN*)malloc(sizeX*sizeY*sizeZ*sizeof(DATATYPEIN));
   free(gradVec);// = (Vector *)malloc(sizeX*sizeY*sizeZ*sizeof(Vector));
   fc = NULL;
   volin = NULL;
   gradVec = NULL;
   cout << "Computing gradient vector field is over!\n" << endl;

   return 0;
}

