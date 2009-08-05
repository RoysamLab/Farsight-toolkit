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

/*  Volume dataset processing
 *  accept a sequence of volumes
 *  Windows version, taken in from Linux version
 *   Author: Xiaosong Yuan, RPI
 *  Modified on Sep. 29, 2005  

 *  Input parameters
 *          1. sizeExpand   
 *          2. preproess          */
#pragma warning(disable : 4996)
//#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
//#include <fstream.h>
//#include <iostream.h>
#include <string.h>
#include <math.h>
#define FillColor 1  //111
#define OBJECT 200
#define COLORseed 0  //usually 0, can be selected as OBJECT

//using namespace std;

#define DATATYPEIN unsigned char
#define DATATYPEOUT unsigned char

const unsigned char  m_NumberOfHistogramBins = 128;

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
double  OtsuThreshold (int sizeX,int sizeY,int sizeZ);



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
	//int NearObjFlag;
	DATATYPEOUT *volout;
	long idx;//, iidx;
	double threshold;
	//float threshold;
	//int kmod8, kdiv8;
	//int FlagIsolated;
	//int NumConnectComp;
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
		
		if (fread(volin, sizeof(DATATYPEIN), sizeX*sizeY*sizeZ, infile) < (unsigned int)(sizeX*sizeY*sizeZ))
		{
			printf("File size is not the same as volume size\n");
		    exit(1);
		}
		//fread(volin,sizeX*sizeY*sizeZ,sizeof(DATATYPEIN), infile);


		printf("Volume Processing # %d ...\n", t);


		// Pre-Processing 

		// by xiao liang, using 3 sigma theory to estimate threshold;
        double meanValue =0.0, VarianceValue =  0.0;

		for (k=0; k<sizeZ; k++)
		{  // threshlod first
			for (j=0; j<sizeY; j++) 
			{
				for (i=0; i<sizeX; i++)
				{
					idx = k *sizeX*sizeY + j *sizeX + i;
					meanValue += volin[idx];
				}
			}
		}

       meanValue = meanValue/(double)(sizeX*sizeY*sizeZ);
 
	   for (k=0; k<sizeZ; k++)
		{  // threshlod first
			for (j=0; j<sizeY; j++) 
			{
				for (i=0; i<sizeX; i++)
				{
					idx = k *sizeX*sizeY + j *sizeX + i;
					VarianceValue += (volin[idx]-meanValue)*(volin[idx]-meanValue);
				}
			}
		}

       VarianceValue =  VarianceValue/(double)(sizeX*sizeY*sizeZ);
	   VarianceValue = sqrt(VarianceValue);

	   double m_threshold=OtsuThreshold (sizeX,sizeY,sizeZ);
	   if (m_threshold > 7||m_threshold<0)threshold =(meanValue-VarianceValue/30); 
	   else
		   threshold = m_threshold;

	   printf ("OTSU optimal threshold %f", threshold);
	
	  // cout << meanValue << VarianceValue;
     // printf("hello %f,hello %f, %f, optimalthreshold= %f ",meanValue,VarianceValue,threshold,m_threshold);
	

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

					if (volin[idx] >= threshold) 
					{}
					/*if (volin[idx] >= meanValue-3*VarianceValue) 
					{
						//volin[idx] = volin[idx];
					}
					*/


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
						//if (blockMax > 255)   blockMax = 255;
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

	free(volin);  // by xiao
	free(volout); // by xiao
	volin=NULL;
	volout=NULL;
	delete []infilename;
	delete []outfilename;

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
	int startx1;//, endx1;
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







// This function is designed to compute the opyimal threshold using OTSU method;
// this algoritm is implemented by xiao liang based on ITK's OTSU algorithm 
double  OtsuThreshold (int sizeX,int sizeY,int sizeZ)

{

  int i,j,k;
  double m_Threshold =0;
  double totalPixels = (double)(sizeX*sizeY*sizeZ);
 

  if ( totalPixels == 0 ) { return m_Threshold; }

  unsigned char MinValue = volin[0], MaxValue =  volin[0];
  double meanValue=0.0, varianceValue=0.0;
  int idx; // just for the location0

  for (k=0; k < sizeZ; k++)
	{  // 
	 for (j=0; j < sizeY; j++) 
			{
				for (i=0; i<sizeX; i++)
				{
					idx = k *sizeX*sizeY + j *sizeX + i;
					meanValue += volin[idx];
					if (volin[idx]> MaxValue) MaxValue = volin[idx];
					if (volin[idx]< MinValue) MinValue = volin[idx];

				}
			}
		}

      printf("Max= %d,Minin= %d", MaxValue,MinValue);
       meanValue = meanValue/totalPixels;
 
	   for (k=0; k<sizeZ; k++)
		{  //
			for (j=0; j<sizeY; j++) 
			{
				for (i=0; i<sizeX; i++)
				{
					idx = k *sizeX*sizeY + j *sizeX + i;
					varianceValue += (volin[idx]-meanValue)*(volin[idx]-meanValue);
				}
			}
		} 


    if ( MinValue >= MaxValue)
    {
	   m_Threshold=MinValue;
       return m_Threshold;
    
    }
   
	m_Threshold = (meanValue-varianceValue/30); 
	  // this step is only initialized a good experimental value for m_Threshold, because the 3D image
	  // is sparse, there are lots of zero values; 

  // create a histogram
  double relativeFrequency[m_NumberOfHistogramBins];

  for ( j = 0; j < m_NumberOfHistogramBins; j++ )
    {
       relativeFrequency[j] = 0.0;
    }

  double binMultiplier = (double) m_NumberOfHistogramBins /(double) ( MaxValue - MinValue);
  printf("binmultiple %f \n", binMultiplier);
 
  
  double Voxelvalue; // temp variable
  unsigned int binNumber;

  for (k=0; k<sizeZ; k++) 
	{  // 
	for (j=0; j<sizeY; j++) 
	 {
	for (i=0; i<sizeX; i++)
	 {
		idx = k *sizeX*sizeY + j *sizeX + i;
        Voxelvalue = volin[idx];

        if ( Voxelvalue == MinValue ) 
         {
         binNumber = 0;
         } // end if 

        else
          {
             binNumber = (unsigned int)(((Voxelvalue - MinValue) * binMultiplier ) - 1);

             if ( binNumber == m_NumberOfHistogramBins ) // in case of rounding errors
              {
                binNumber -= 1;
              }// end if 
           }// end else

       //  printf("binNumber???? = %f  MaxValue=%f \n", binNumber+1,MaxValue);
         relativeFrequency[binNumber] += 1.0;
   
	}//
	}
   }
// test 

//	for (i=0;i<m_NumberOfHistogramBins;i++)
//		printf ( "%d bin = %f,  ",i, relativeFrequency[i]);
 
  // normalize the frequencies
  double totalMean = 0.0;
  for ( j = 0; j < m_NumberOfHistogramBins; j++ )
    {
    relativeFrequency[j] /= totalPixels;
    totalMean += (j+1) * relativeFrequency[j];
    }


  // compute Otsu's threshold by maximizing the between-class
  // variance
  double freqLeft = relativeFrequency[0];
  double meanLeft = 1.0;
  double meanRight = ( totalMean - freqLeft ) / ( 1.0 - freqLeft );

  double maxVarBetween = freqLeft * ( 1.0 - freqLeft ) * sqrt( meanLeft - meanRight );
  int maxBinNumber = 0;

  double freqLeftOld = freqLeft;
  double meanLeftOld = meanLeft;

  for ( j = 1; j < m_NumberOfHistogramBins; j++ )
    {
    freqLeft += relativeFrequency[j];
    meanLeft = ( meanLeftOld * freqLeftOld + 
                 (j+1) * relativeFrequency[j] ) / freqLeft;
    if (freqLeft == 1.0)
      {
      meanRight = 0.0;
      }
    else
      {
      meanRight = ( totalMean - meanLeft * freqLeft ) / ( 1.0 - freqLeft );
      }
    double varBetween = freqLeft * ( 1.0 - freqLeft ) * sqrt( meanLeft - meanRight );
   
    if ( varBetween > maxVarBetween )
      {
      maxVarBetween = varBetween;
      maxBinNumber = j;
      }

    // cache old values
    freqLeftOld = freqLeft;
    meanLeftOld = meanLeft; 

    } 

    m_Threshold = double( MinValue + ( maxBinNumber + 1 ) / binMultiplier );
	printf("m_threshold=%f ",m_Threshold );
	return m_Threshold; 
}
