/*  Volume dataset processing
/*  accept a sequence of volumes
/*  Windows version, taken in from Linux version
/*   Author: Xiaosong Yuan, RPI
/*  Modified on Sep. 29, 2005  

/*  Input parameters
/*          1. sizeExpand   
/*          2. preproess          */

//#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
//#include <fstream.h>
//#include <iostream.h>
#include <string.h>
#define FillColor 1  //111
#define OBJECT 200
#define COLORseed 0  //usually 0, can be selected as OBJECT

//using namespace std;

#define DATATYPEIN unsigned char
#define DATATYPEOUT unsigned char


DATATYPEIN *volin;
int sizeX, sizeY, sizeZ, sizeTime;
long voxelCount=0;

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

StackOfPoints *stack = new StackOfPoints;
int stacktop = 0;

void lineseed(int, int, int);
void push (Position pos);
Position pop (void);
bool empty(void);
int findstartx(Position ps);
int findendx(Position ps);
void spread(Position pos, int startx, int endx, int direction);
void lineSeed();


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
	int NearObjFlag;
	DATATYPEOUT *volout;
	long idx, iidx;
	float threshold;
	int kmod8, kdiv8;
	int FlagIsolated;
	int NumConnectComp;
	int sizeExpand = 0;  //10;
	DATATYPEOUT blockMax;
	int timesDilate;
	int border;

	infilename = argv[1];
	sizeX = atoi(argv[2]);
	sizeY = atoi(argv[3]);
	sizeZ = atoi(argv[4]);
	//sizeTime = atoi(argv[5]);
	sizeTime = 1;
	outfilename = argv[5];
	threshold = atof(argv[6]);

	volin = (DATATYPEIN*)malloc(sizeX*sizeY*(sizeZ+sizeExpand*2)*sizeof(DATATYPEIN));
	volout = (DATATYPEOUT*)malloc(sizeX*sizeY*(sizeZ+sizeExpand*2)*sizeof(DATATYPEOUT));

	if((infile=fopen(infilename,"rb"))==NULL)
			{printf("Input file open error!\n");
			 exit(-1);
			}

	if((outfile=fopen(outfilename,"wb"))==NULL)
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

	for (t=0; t<sizeTime; t++) {
		
		if (fread(volin, sizeof(DATATYPEIN), sizeX*sizeY*sizeZ, infile) < sizeX*sizeY*sizeZ)
		{
			printf("File size is not the same as volume size\n");
		    exit(1);
		}
		//fread(volin,sizeX*sizeY*sizeZ,sizeof(DATATYPEIN), infile);


		printf("Volume Processing # %d ...\n", t);


		// Pre-Processing 
		for (k=0; k<(sizeZ+sizeExpand*2); k++)
			for (j=0; j<sizeY; j++)
				for (i=0; i<sizeX; i++) {
					volout[k *sizeX*sizeY + j *sizeX + i] = 0; 
				}  //initial to zeros
		for (k=0; k<sizeZ; k++) {  // threshlod first
			for (j=0; j<sizeY; j++) {
				for (i=0; i<sizeX; i++)
				{
					idx = k *sizeX*sizeY + j *sizeX + i;
					if (volin[idx] >= threshold) {
						//volin[idx] = volin[idx];
					}
					else
						volin[idx] = 0;
				}
			}
		}   
/*
		for (k=0; k<sizeZ; k++)
			for (j=0; j<sizeY; j++)
				for (i=0; i<sizeX; i++) {
					volout[(k+sizeExpand)*sizeX*sizeY + j*sizeX + i] = volin[k*sizeX*sizeY + j*sizeX + i];
				}  // expand the size and pad with zeros
		sizeZ = sizeZ + 2* sizeExpand;
		for (k=0; k<sizeZ; k++)
			for (j=0; j<sizeY; j++)
				for (i=0; i<sizeX; i++) {
					volin[k *sizeX*sizeY + j*sizeX + i] = volout[k *sizeX*sizeY + j*sizeX + i];
				}  // save back to volin

*/

		// Method 1: Cutting a part of the object
/*		for (k=0; k<sizeZ; k++) {
			for (j=0; j<sizeY; j++) {
				for (i=0; i<sizeX; i++)
				{
					idx = k *sizeX*sizeY + j *sizeX + i;
					if (volin[idx] == 0)
						volout[idx] = 0;
					    else
						volout[idx] = OBJECT;
					if ((k<=2 || k>=sizeZ-3)||(j<=2 || j>=sizeZ-3)||(i<=2 || i>=sizeZ-3))
						volout[idx] = 0;
				}
			}
		}
*/

		// Method 2: Dilation of the object
	   timesDilate = 1;
	   border = 3;
	   while (timesDilate >0 ) {
		for (k=border; k<sizeZ-border; k++) {
			for (j=border; j<sizeY-border; j++) {
				for (i=border; i<sizeX-border; i++)
				{
				    blockMax = volin[k *sizeX*sizeY + j *sizeX + i];
				    for (kk=-1; kk<=1; kk++)
				      for (jj=-1; jj<=1; jj++)
						for (ii=-1; ii<=1; ii++) {
							if(volin[(k+kk)*sizeX*sizeY + (j+jj)*sizeX + (i+ii)] > blockMax) 
								blockMax = volin[(k+kk)*sizeX*sizeY + (j+jj)*sizeX + (i+ii)];
						}
					// Keep the peak of the original intensity
					if (blockMax == volin[k *sizeX*sizeY + j *sizeX + i] && blockMax != 0)  {
						blockMax = blockMax + 1;
						if (blockMax > 255)   blockMax = 255;
					}
					volout[k *sizeX*sizeY + j *sizeX + i] = blockMax;
				}
			}
		}

		// copy volout back to volin for the next dilation
		for (k=0; k<sizeZ; k++) 
			for (j=0; j<sizeY; j++) 
				for (i=0; i<sizeX; i++)  {
					volin[k *sizeX*sizeY + j *sizeX + i] = volout[k *sizeX*sizeY + j *sizeX + i];
				}
        timesDilate--;
	   }


		// Method 3: Threshold the voxel intensity
/*		for (k=0; k<sizeZ; k++) {
			for (j=0; j<sizeY; j++) {
				for (i=0; i<sizeX; i++)
				{
					idx = k *sizeX*sizeY + j *sizeX + i;
					if (volin[idx] >= threshold) {
						volout[idx] = OBJECT;
					}
					else
						volout[idx] = 0;
				}
			}
		}
*/

		// Method 4: Flood/scanlines fill
	/*	lineSeed();
		for (k=0; k<sizeZ; k++) {
			for (j=0; j<sizeY; j++) {
				for (i=0; i<sizeX; i++)
				{
					idx = k *sizeX*sizeY + j *sizeX + i;
					if(volin[idx] != FillColor) {
						volout[idx] = OBJECT;
					}
					else {
						volout[idx] = 0;
					}
				}
			}
		}


		// Method 5: Put 64x64x32 slice by slice into one whole image
		int NumXImages = 8;
		int NumYImages = 4;
		for (k=0; k<sizeZ; k++)
			for (j=0; j<sizeY; j++)
				for (i=0; i<sizeX; i++) {
					idx = k *sizeX*sizeY + j *sizeX + i;
					kmod8 = k % NumXImages;
					kdiv8 = k / NumXImages;
					ii = kmod8*sizeX + i;
					jj = kdiv8*sizeY + j;
					kk = 0;
					iidx = kk*sizeX*sizeY*NumXImages*NumYImages + jj*sizeX*NumXImages + ii;
					volout[iidx] = volin[idx];

		}
*/

		//Method 6: Take away each isolated dot (noise) and keep clusters
		// i.e., according to #voxels of each connected components must be >= 2
/*		for (k=0; k<sizeZ; k++)
		    for (j=0; j<sizeY; j++)
			for (i=0; i<sizeX; i++)  {
				idx = k *sizeX*sizeY + j *sizeX + i;
				NumConnectComp = 0;
				for (kk=-1; kk<=1; kk++)
				    for (jj=-1; jj<=1; jj++)
					for (ii=-1; ii<=1; ii++) {
						iidx = idx + kk*sizeX*sizeY + jj*sizeX + ii;
						if (iidx == idx || iidx < 0 || iidx >= sizeX*sizeY*sizeZ)  continue;
						if(volin[iidx] != 0)  NumConnectComp++;
				}
				if(NumConnectComp < 3)
					volout[idx] = 0;
				   else
					volout[idx] = volin[idx];
		}


		// Method 7: Subtract one image from another
		for (k=0; k<sizeZ; k++) {
			for (j=0; j<sizeY; j++) {
				for (i=0; i<sizeX; i++)
				{
					idx = k *sizeX*sizeY + j *sizeX + i;
					volout[idx] = volin[idx];
					if (volin[idx] == 0)  continue;
					if (volin2[idx] != 0)  volout[idx] = 0;
				}
			}
		}
*/

		// Method 8: Creat a small test volume
/*		sizeX=5;
		sizeY=5;
		sizeZ=5;
		for (k=0; k<sizeZ; k++) {
			for (j=0; j<sizeY; j++) {
				for (i=0; i<sizeX; i++)
				{
					if (i==j && i==k) {
                        volout[k *sizeX*sizeY + j *sizeX + i] = 100;
					}
					else  {
						volout[k *sizeX*sizeY + j *sizeX + i] = 0;
					}

					if (i==2 && j==2 && k==3)
						volout[k *sizeX*sizeY + j *sizeX + i] = 100;
				}
			}
		}
*/

		fwrite(volout, sizeX*sizeY*sizeZ, sizeof(DATATYPEOUT), outfile);

	}


	fclose(infile);
	fclose(outfile);
	printf("Done \n");
    return 0;
}






void push (Position pos)
{
	StackOfPoints *temp = new StackOfPoints;
	temp->ps = pos;
	temp->next = stack;
	stack = temp;
	stacktop ++;
}

Position pop (void)
{
	StackOfPoints *temp;
	Position outp = stack->ps;
	temp = stack;
	stack = stack->next;
	stacktop --;
	delete temp;
	return outp;
}

bool empty(void)
{
	return stacktop == 0 ;
}

int findstartx(Position ps)
{
	while (volin[ps.z *sizeX*sizeY + ps.y *sizeX + ps.x] == COLORseed && ps.x >=0)
	{
		ps.x --;
	}
	ps.x ++;
	return ps.x;
}

int findendx(Position ps)
{
	while (volin[ps.z *sizeX*sizeY + ps.y *sizeX + ps.x] == COLORseed && ps.x < sizeX)
	{
		ps.x ++;
	}
	ps.x --;
	return ps.x;
}

void spread(Position pos, int startx, int endx, int direction)
{
	Position pos1; // in a new row
	int newy, newz;
	int startx0, endx0;
	int startx1, endx1;
	int laststartx;

	switch (direction)
	{
		case 1:
			newy = pos.y + 1;
			newz = pos.z;
			break;
		case 2:
			newy = pos.y - 1;
			newz = pos.z;
			break;
		case 3:
			newy = pos.y;
			newz = pos.z +1;
			break;
		case 4:
			newy = pos.y;
			newz = pos.z -1;
			break;
	}

	if (newy < sizeY && newy >= 0 && newz < sizeZ && newz >=0) // within boundary
	{
		startx0 = startx;
		endx0 = endx;
		laststartx = -1;
		for (int i = startx0; i<= endx0; i++ )
		{
			if (volin[newz *sizeX*sizeY + newy *sizeX + i] == COLORseed)
			{
				pos1.x = i;
				pos1.y = newy;
				pos1.z = newz;
				startx1 = findstartx(pos1);
				if (startx1 != laststartx)
				{
					pos1.x = startx1;
					pos1.y = newy;
					pos1.z = newz;
					push(pos1);
					laststartx = startx1;
				}
			}


		}

	}

}


void lineSeed()
{
	Position pos;
	int startx, endx;

	Position initialSeed = {0, 0, 0}; // set the seed position
	push(initialSeed);
	while (!empty())
	{
		pos = pop();
		if(volin[pos.z *sizeX*sizeY + pos.y *sizeX + pos.x] == COLORseed)
		{
			startx = pos.x;
			endx = findendx(pos);
			for (int i = startx; i<= endx; i++)
				volin[pos.z *sizeX*sizeY + pos.y *sizeX + i] = FillColor;
			// four directions
			spread(pos, startx, endx, 1);
			spread(pos, startx, endx, 2);
			spread(pos, startx, endx, 3);
			spread(pos, startx, endx, 4);

		}

	}

}

