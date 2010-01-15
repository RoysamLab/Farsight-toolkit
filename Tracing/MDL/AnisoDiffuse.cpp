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

/*  Volume Density Iso- Anisotropic Diffusion
 *   Author: Xiaosong Yuan, RPI
 *  Created on Nov. 16, 2005  

 *  Input parameters
 *            
 *                   */
#if defined(_MSC_VER)
#pragma warning(disable : 4996)
#endif
#include <stdlib.h>
#include <iostream>
#include <cstdio>
#include <string.h>
#include "math.h"

using std::cerr;
using std::cout;
using std::endl;

//using namespace std;

#define DATATYPEIN unsigned char
#define DATATYPEOUT unsigned char

#define Aniso 1      // Anisotropic 1; Isotropic 0


DATATYPEIN *volin;
int sizeX, sizeY, sizeZ, sizeTime;
long voxelCount=0;

struct  Position
{
  int x;
  int y;
  int z;
};

struct  Vector
{
  float xd;   // For large datasets, using double will cause memory shortage
  float yd;
  float zd;
};


int main(int argc, char *argv[])
{
   if(argc < 5) {
      cerr << "Usage: <input file> <sizeX> <sizeY> <sizeZ> <output file>"
           << endl;
      return 1;
  }
  //string s;
  FILE *infile;
  FILE *outfile;
    Vector *gradVec;
  DATATYPEOUT *volout;

  int i,j,k;
  //int t;
  
  int d1, d2;
  int sls, sz;
  long idx, iidx1, iidx2;
  //long iidx;
  int timesDiffuse;
  int border;
  //float kernelWeight[3][3];
  double kernelWeight[3][3];
  double diverg;
  double lamda;
  double voxelUpdate;
  double k_factor;
  //float k_factor;
  double Dxy1, Dxy2;
  

  sizeX = atoi(argv[2]);
  sizeY = atoi(argv[3]);
  sizeZ = atoi(argv[4]);
  timesDiffuse =1;

  sls = sizeX*sizeY;    // slice size
  sz = sls*sizeZ;


  volin = (DATATYPEIN*)malloc(sizeX*sizeY*sizeZ*sizeof(DATATYPEIN));
  volout = (DATATYPEOUT*)malloc(sizeX*sizeY*sizeZ*sizeof(DATATYPEOUT));
  gradVec = (Vector *)malloc(sizeX*sizeY*sizeZ*sizeof(Vector));

  if((infile=fopen(argv[1],"rb"))==NULL)
    {
    cerr << "Input file open error!" << endl;
    return -1;
    }

  if((outfile=fopen(argv[5],"wb"))==NULL)
    {
    cerr << "Output file open error!" << endl;
    return -1;
    }

  if (fread(volin, sizeof(DATATYPEIN), sizeX*sizeY*sizeZ, infile) < (unsigned int)(sizeX*sizeY*sizeZ) )
    {
    cerr << "File size is not the same as volume size" << endl;
    return 1;
    }

  // define positive half kernel of derivative 
  kernelWeight[0][0] = 1; kernelWeight[0][1] = 2; kernelWeight[0][2] = 1;
  kernelWeight[1][0] = 2; kernelWeight[1][1] = 3; kernelWeight[1][2] = 2;
  kernelWeight[2][0] = 1; kernelWeight[2][1] = 2; kernelWeight[2][2] = 1;

  double AveGradient=0.0 ;  // by xiao liang

  double tempGradient;

  
  while (timesDiffuse >0 )
    {

    cout << "Anisotropic Diffusion #" << timesDiffuse << " ..." << endl;

    // Pre-Processing 
    for (idx=0; idx<sz; idx++)
      {
      volout[idx] = 0; 
      }  //initial to zeros

    border = 1;
    // Compute Gradient
    for (k = border; k < sizeZ-border; k++)
      for (j = border; j < sizeY-border; j++)
      for (i = border; i < sizeX-border; i++)
      if(volin[k*sls + j*sizeX + i]>0) // by xiao
      {// by xiao
      {
        idx = k*sls + j*sizeX + i;

      
        gradVec[idx].xd = 0;
        gradVec[idx].yd = 0;
        gradVec[idx].zd = 0;
                
        for (d2 = -1; d2 <= 1; d2++)
          for (d1 = -1; d1 <= 1; d1++) 
          {
            iidx1 = (k+d2) *sls + (j+d1) *sizeX + i-1;
            iidx2 = (k+d2) *sls + (j+d1) *sizeX + i+1;
            gradVec[idx].xd = (float) (gradVec[idx].xd + kernelWeight[d2+1][d1+1]*(volin[iidx2] - volin[iidx1]));

            iidx1 = (k+d2) *sls + (j-1) *sizeX + (i+d1);
            iidx2 = (k+d2) *sls + (j+1) *sizeX + (i+d1);
            gradVec[idx].yd =(float) (gradVec[idx].yd + kernelWeight[d2+1][d1+1]*(volin[iidx2] - volin[iidx1]));

            iidx1 = (k-1) *sls + (j+d2) *sizeX + (i+d1);
            iidx2 = (k+1) *sls + (j+d2) *sizeX + (i+d1);
            gradVec[idx].zd = (float)(gradVec[idx].zd + kernelWeight[d2+1][d1+1]*(volin[iidx2] - volin[iidx1]));
            }

          tempGradient = gradVec[idx].xd * gradVec[idx].xd + gradVec[idx].yd *gradVec[idx].yd+ gradVec[idx].zd * gradVec[idx].zd;
          tempGradient = sqrt(tempGradient);

          AveGradient += tempGradient;

           }
             
      } // end if

      if(sz>0)
               AveGradient/=sz;
      else 
               AveGradient = 400;

      k_factor =   AveGradient;
    // Begin the diffusion based on the above gradients
#if Aniso   // Anisotropic diffusion
    border = 2;
    lamda = 0.02;  // 0.02
    for (k=border; k<sizeZ-border; k++) {
      for (j=border; j<sizeY-border; j++) {
        for (i=border; i<sizeX-border; i++)

          if(volin[k*sls + j*sizeX + i]>0) // by xiao
        {// by xiao

        {
          idx = k*sls + j*sizeX + i;
          diverg = 0;
          iidx1 = k *sls + j *sizeX + i-1;
          iidx2 = k *sls + j *sizeX + i+1;
          Dxy1 = exp( -(gradVec[iidx1].xd/k_factor) * (gradVec[iidx1].xd/k_factor) );
          Dxy2 = exp( -(gradVec[iidx2].xd/k_factor) * (gradVec[iidx2].xd/k_factor) );
          diverg = diverg + (Dxy2 * gradVec[iidx2].xd - Dxy1 * gradVec[iidx1].xd);
          
          iidx1 = k *sls + (j-1) *sizeX + i;
          iidx2 = k *sls + (j+1) *sizeX + i;
          Dxy1 = exp( -(gradVec[iidx1].yd/k_factor) * (gradVec[iidx1].yd/k_factor) );
          Dxy2 = exp( -(gradVec[iidx2].yd/k_factor) * (gradVec[iidx2].yd/k_factor) );
          diverg = diverg + (Dxy2 * gradVec[iidx2].yd - Dxy1 * gradVec[iidx1].yd);
          
          iidx1 = (k-1) *sls + j *sizeX + i;
          iidx2 = (k+1) *sls + j *sizeX + i;
          Dxy1 = exp( -(gradVec[iidx1].zd/k_factor) * (gradVec[iidx1].zd/k_factor) );
          Dxy2 = exp( -(gradVec[iidx2].zd/k_factor) * (gradVec[iidx2].zd/k_factor) );
          diverg = diverg + (Dxy2 * gradVec[iidx2].zd - Dxy1 * gradVec[iidx1].zd);

          voxelUpdate = volin[idx] + lamda * diverg;
          if (voxelUpdate<0)   voxelUpdate=0;
          if (voxelUpdate>255)   voxelUpdate=255;

          volout[idx] = (int)(voxelUpdate+0.5);  // rounding to int
        }
        } // end if 
      }
    }

#else  // Isotropic diffusion
    border = 2;
    lamda = 0.02;
    for (k=border; k<sizeZ-border; k++) {
      for (j=border; j<sizeY-border; j++) {
        for (i=border; i<sizeX-border; i++)
        {
          idx = k*sls + j*sizeX + i;
          diverg = 0;
          iidx1 = k *sls + j *sizeX + i-1;
          iidx2 = k *sls + j *sizeX + i+1;
          diverg = diverg + (gradVec[iidx2].xd - gradVec[iidx1].xd);
          iidx1 = k *sls + (j-1) *sizeX + i;
          iidx2 = k *sls + (j+1) *sizeX + i;
          diverg = diverg + (gradVec[iidx2].yd - gradVec[iidx1].yd);
          iidx1 = (k-1) *sls + j *sizeX + i;
          iidx2 = (k+1) *sls + j *sizeX + i;
          diverg = diverg + (gradVec[iidx2].zd - gradVec[iidx1].zd);

          voxelUpdate = volin[idx] + lamda * diverg;
          if (voxelUpdate<0)   voxelUpdate=0;
          if (voxelUpdate>255)   voxelUpdate=255;

          volout[idx] = (int)(voxelUpdate+0.5);  // rounding to int
        }
      }
    }

#endif





    // copy volout back to volin for the next dilation
    for (idx=0; idx<sz; idx++)  {
      volin[idx] = volout[idx];
    }

    timesDiffuse--;
  }  // End of timesDiffuse "while"


  // Post-smoothing

/*
   double blockAve;
   int ii, jj, kk;

  border = 1;
  for (k=border; k<sizeZ-border; k++)
    for (j=border; j<sizeY-border; j++)
      for (i=border; i<sizeX-border; i++)   {
        blockAve = 0;
        for (kk=-1; kk<=1; kk++)
            for (jj=-1; jj<=1; jj++)
            for (ii=-1; ii<=1; ii++) {
              blockAve += volin[(k+kk)*sizeX*sizeY + (j+jj)*sizeX + (i+ii)];
        }
        volout[k *sizeX*sizeY + j *sizeX + i] = blockAve / 27;
  }
*/  // by xiao liang

  fwrite(volout, sizeX*sizeY*sizeZ, sizeof(DATATYPEOUT), outfile);
  

  FILE *mhdfile;
  
  if((mhdfile=fopen("Aniso_Diffused.mhd","w"))==NULL)
    {
    cerr << "output file open error!" << endl;
    return -1;
    }
  fprintf (mhdfile,"ObjectType = Image\n");
  fprintf (mhdfile,"NDims = 3\n");
  fprintf (mhdfile,"BinaryData = True\n");
  fprintf (mhdfile,"BinaryDataByteOrderMSB = False\n");
  fprintf (mhdfile,"CompressedData = False\n");
  fprintf (mhdfile,"TransformMatrix = 1 0 0 0 1 0 0 0 1\n");
  fprintf (mhdfile,"Offset = 0 0 0\n");
  fprintf (mhdfile,"CenterOfRotation = 0 0 0\n");
  fprintf (mhdfile,"AnatomicalOrientation = RAI\n");
  fprintf (mhdfile,"ElementSpacing = 1 1 1\n");
  fprintf (mhdfile,"DimSize = %d %d %d\n",sizeX,sizeY,sizeZ);
  fprintf (mhdfile,"ElementType = MET_UCHAR\n");
  fprintf (mhdfile,"ElementDataFile = Aniso_Diffused.raw\n");
  fclose(mhdfile);

  
    // memory release, added by xiao
  free(volin);// = (DATATYPEIN*)malloc(sizeX*sizeY*sizeZ*sizeof(DATATYPEIN));
  free(volout);//= (DATATYPEOUT*)malloc(sizeX*sizeY*sizeZ*sizeof(DATATYPEOUT));
  free(gradVec);// = (Vector *)malloc(sizeX*sizeY*sizeZ*sizeof(Vector));
  fclose(infile);
  fclose(outfile);
  cout << "Anistropic Diffusion is Done" << endl;
  return 0;
}


