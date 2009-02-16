// Label the connected components of the input volume with zero background
// Remove the connected components with small number of voxels
// --- Input: original volume
// --- Output: removed small objects
// --- Author: Xiaosong Yuan, RPI
// --- Date: 10/3/2005

#include <stdlib.h>
#include <stdio.h>
//#include <fstream.h>
//#include <iostream.h>
#include <string.h>
#define OBJ_BEFORE_CONNCOMP -1
#define ThresCompVoxels 2 //200   // Determine size of components removed

//using namespace std;

#define DATATYPEIN unsigned char
#define DATATYPEOUT unsigned char


DATATYPEIN *volin;
int sizeX, sizeY, sizeZ;
long voxelCount=0;
int VertIndex = 0;   // Number of vertex in each component
int CompIndex = 0;   // Number of components in the volume

struct  Position
{
	int x;
	int y;
	int z;
};

struct StackOfPoints
{
	Position ps;
	StackOfPoints *next;
};


void DepthFirstSearch_N6(int *volIndex, int i, int j, int k);


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
	int i,j,k, t;
	int ii, jj, kk;
	DATATYPEOUT *volout;
	int *volIndex;
	long idx, iidx;
	float threshold;
	int vertHistComp[100000];

	infilename = argv[1];
	sizeX = atoi(argv[2]);
	sizeY = atoi(argv[3]);
	sizeZ = atoi(argv[4]);
	outfilename = argv[5];
	threshold = atof(argv[6]);

	volin = (DATATYPEIN*)malloc(sizeX*sizeY*sizeZ*sizeof(DATATYPEIN));
	volout = (DATATYPEOUT*)malloc(sizeX*sizeY*sizeZ*sizeof(DATATYPEOUT));
	volIndex = (int*)malloc(sizeX*sizeY*sizeZ*sizeof(int));

	if((infile=fopen(infilename,"rb"))==NULL)
			{printf("Input file open error!\n");
			 exit(-1);
			}

	if((outfile=fopen(outfilename,"wb"))==NULL)
			{printf("Output file open error!\n");
			 exit(-1);
			}

	fread(volin,sizeX*sizeY*sizeZ, sizeof(DATATYPEIN), infile);



	printf("Label Connected Components ...\n");

	//Initialize volIndex to binary value (either Zeros or OBJ_BEFORE_CONNCOMP as -1)
	for (k=0; k<sizeZ; k++)
		for (j=0; j<sizeY; j++)
			for (i=0; i<sizeX; i++) {
				if (volin[k *sizeX*sizeY + j *sizeX + i] == 0)  {
					volIndex[k *sizeX*sizeY + j *sizeX + i] = 0;
				}
				else
					volIndex[k *sizeX*sizeY + j *sizeX + i] = OBJ_BEFORE_CONNCOMP;
			} 


	// Label the Connected Components with 1,2,..MaxIndex.
	for (k=0; k<sizeZ; k++) {
		for (j=0; j<sizeY; j++) {
			for (i=0; i<sizeX; i++)
			{
				idx = k *sizeX*sizeY + j *sizeX + i;
				if(volIndex[idx] == OBJ_BEFORE_CONNCOMP) {
					CompIndex ++;
					VertIndex = 0;
                    DepthFirstSearch_N6(volIndex, i, j, k);
					vertHistComp[CompIndex] = VertIndex; //save the number of vertex in each component
				}
			}
		}
	}

	// Turn small components to zeros and keep large components
	for (k=0; k<sizeZ; k++) {
		for (j=0; j<sizeY; j++) {
			for (i=0; i<sizeX; i++)
			{
				idx = k *sizeX*sizeY + j *sizeX + i;
				if(volIndex[idx]!= 0 && vertHistComp[volIndex[idx]] >= ThresCompVoxels) {
					volout[idx] = volin[idx];
				}
				else {
					volout[idx] = 0;
				}
			}
		}
	}

	fwrite(volout, sizeX*sizeY*sizeZ, sizeof(DATATYPEOUT), outfile);

	fclose(infile);
	fclose(outfile);
	printf("Done \n");
	return 0;
}


void DepthFirstSearch_N6(int *volIndex, int i, int j, int k) {
  int ii, jj, kk;
  VertIndex++;   //increase the number of vertices traversed
  volIndex[k *sizeX*sizeY + j *sizeX + i] = CompIndex; 
  //for all vertices adj. to (i,j,k)
  for (kk=-1; kk<=1; kk++)
     for (jj=-1; jj<=1; jj++)
        for (ii=-1; ii<=1; ii++) {
	        //not consider the point itself && consider only neighbor-6
	        if (ii==0 && jj==0 && kk==0 || (abs(ii)+abs(jj)+abs(kk))>1)  continue;
	        if ((i+ii)>=0 && (i+ii)<sizeX && (j+jj)>=0 && (j+jj)<sizeY && (k+kk)>=0 && (k+kk)<sizeZ) {
	            if (volIndex[(k+kk)*sizeX*sizeY +(j+jj)*sizeX +(i+ii)] == OBJ_BEFORE_CONNCOMP) {
				    DepthFirstSearch_N6(volIndex, i+ii, j+jj, k+kk);
			}
	    }
	}
}



