/*  Volume Density Iso- Anisotropic Diffusion
 *   Author: Xiaosong Yuan, RPI
 *  Created on Nov. 16, 2005  

 *  Input parameters
 *            
 *                   */

#include <stdlib.h>
#include <stdio.h>
//#include <fstream.h>
//#include <iostream.h>
#include <string.h>
#include "math.h"

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
	
	if(argc < 6) {
	    printf("\nUsage: <input file> <sizeX> <sizeY> <sizeZ> <output file> <threshold>.\n");
	    exit(1);
	}

	//string s;
	FILE *infile;
	FILE *outfile;
	char *infilename = new char[80];
	char *outfilename = new char[80];
    Vector *gradVec;
	DATATYPEOUT *volout;

	int i,j,k;
	int ii, jj, kk;
	int d1, d2;
	int sls, sz;
	long idx, iidx1, iidx2;
	int timesDiffuse;
	int border;
	//float kernelWeight[3][3];
	double kernelWeight[3][3];
	double diverg;
	double lamda;
	//float threshold;
	double threshold;
	double voxelUpdate;
	float k_factor;
	double Dxy1, Dxy2;
	double blockAve;

	infilename = argv[1];
	sizeX = atoi(argv[2]);
	sizeY = atoi(argv[3]);
	sizeZ = atoi(argv[4]);
	outfilename = argv[5];
	threshold = atof(argv[6]);

	sls = sizeX*sizeY;		// slice size
	sz = sls*sizeZ;


	volin = (DATATYPEIN*)malloc(sizeX*sizeY*sizeZ*sizeof(DATATYPEIN));
	volout = (DATATYPEOUT*)malloc(sizeX*sizeY*sizeZ*sizeof(DATATYPEOUT));
	gradVec = (Vector *)malloc(sizeX*sizeY*sizeZ*sizeof(Vector));
     


	/*
	if((infile=fopen(infilename,"rb"))==NULL)
			{printf("Input file open error!\n");
			 exit(-1);
			}

	if((outfile=fopen(outfilename,"wb"))==NULL)
			{printf("Output file open error!\n");
			 exit(-1);
			}
	*/

	 errno_t err; 
       if((err=fopen_s(&infile,infilename,"rb"))!=NULL)
			{printf("Input file open error!\n");
			 exit(-1);
			}

	if((err=fopen_s(&outfile,outfilename,"wb"))!=NULL)
			{printf("Output file open error!\n");
			 exit(-1);
			}

/*	if 0 	// read another file if necessary
	FILE *infile2;
	unsigned char *volin2;
	volin2 = (unsigned char*)malloc(sizeX*sizeY*sizeZ*sizeof(unsigned char));
	if((infile2=fopen("contrast_thres.64x64x32.img","rb"))==NULL)
			{printf("Input file open error!\n");
			 exit(-1);
			}
	fread(volin2,sizeX*sizeY*sizeZ,sizeof(unsigned char), infile2);
*/
	if (fread(volin, sizeof(DATATYPEIN), sizeX*sizeY*sizeZ, infile) < (unsigned int)(sizeX*sizeY*sizeZ))
	{
		printf("File size is not the same as volume size\n");
		exit(1);
	}

	//fread(volin,sizeX*sizeY*sizeZ,sizeof(DATATYPEIN), infile);

	// define positive half kernel of derivative 
	kernelWeight[0][0] = 1; kernelWeight[0][1] = 2; kernelWeight[0][2] = 1;
	kernelWeight[1][0] = 2; kernelWeight[1][1] = 3; kernelWeight[1][2] = 2;
	kernelWeight[2][0] = 1; kernelWeight[2][1] = 2; kernelWeight[2][2] = 1;

	
    timesDiffuse = 2;
    while (timesDiffuse >0 ) {

		printf("Volume Processing # %d...\n", timesDiffuse);

		// Pre-Processing 
		for (idx=0; idx<sz; idx++)  {
			volout[idx] = 0; 
		}  //initial to zeros

		border = 1;
		// Compute Gradient
		for (k = border; k < sizeZ-border; k++)
		  for (j = border; j < sizeY-border; j++)
			for (i = border; i < sizeX-border; i++)
			{
				idx = k*sls + j*sizeX + i;
				gradVec[idx].xd = 0;
				gradVec[idx].yd = 0;
				gradVec[idx].zd = 0;
				for (d2 = -1; d2 <= 1; d2++)
					for (d1 = -1; d1 <= 1; d1++)  {
						iidx1 = (k+d2) *sls + (j+d1) *sizeX + i-1;
						iidx2 = (k+d2) *sls + (j+d1) *sizeX + i+1;
						gradVec[idx].xd =(float) (gradVec[idx].xd + kernelWeight[d2+1][d1+1]*(volin[iidx2] - volin[iidx1]));

						iidx1 = (k+d2) *sls + (j-1) *sizeX + (i+d1);
						iidx2 = (k+d2) *sls + (j+1) *sizeX + (i+d1);
						gradVec[idx].yd = (float) (gradVec[idx].yd + kernelWeight[d2+1][d1+1]*(volin[iidx2] - volin[iidx1]));

						iidx1 = (k-1) *sls + (j+d2) *sizeX + (i+d1);
						iidx2 = (k+1) *sls + (j+d2) *sizeX + (i+d1);
						gradVec[idx].zd = (float)  (gradVec[idx].zd + kernelWeight[d2+1][d1+1]*(volin[iidx2] - volin[iidx1]));

				}
		}


		// Begin the diffusion based on the above gradients
#if Aniso   // Anisotropic diffusion
		border = 2;
		lamda = 0.02;  // 0.02
		k_factor = 800;
		for (k=border; k<sizeZ-border; k++) {
			for (j=border; j<sizeY-border; j++) {
				for (i=border; i<sizeX-border; i++)
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
	border = 1;
	for (k=border; k<sizeZ-border; k++)
		for (j=border; j<sizeY-border; j++)
			for (i=border; i<sizeX-border; i++) 	{
				blockAve = 0;
				for (kk=-1; kk<=1; kk++)
				    for (jj=-1; jj<=1; jj++)
						for (ii=-1; ii<=1; ii++) {
							blockAve += volin[(k+kk)*sizeX*sizeY + (j+jj)*sizeX + (i+ii)];
				}
				volout[k *sizeX*sizeY + j *sizeX + i] = (unsigned char) (blockAve / 27);  // by xiao liang     blockAve/27
	}


	fwrite(volout, sizeX*sizeY*sizeZ, sizeof(DATATYPEOUT), outfile);

	fclose(infile);
	fclose(outfile);
	printf("Done \n");
    return 0;
}


