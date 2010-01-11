/* 
 * Copyright 2009 Rensselaer Polytechnic Institute
 * This program is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
 * for more details.
 * 
 * You should have received a copy of the GNU General Public License along 
 * with this program; if not, write to the Free Software Foundation, Inc., 
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

//
#include "local_max_clust_3D.h"
#include<math.h>

using namespace std;

void get_maximum(float* A, int r1, int r2, int c1, int c2, int z1, int z2, int* rx, int* cx, int* zx, int R, int C, int Z)
{
	float mx = A[(z1*R*C)+(r1*C)+c1];//A[r1][c1][z1];
    rx[0] = r1;
    cx[0] = c1;
	zx[0] = z1;
    for(int i=r1; i<=r2; i++)
    {
        for(int j=c1; j<=c2; j++)
        {
			for(int k=z1; k<=z2; k++)
			{
				if(A[(k*R*C)+(i*C)+j]>=mx)
				{
					mx = A[(k*R*C)+(i*C)+j];//A[i][j][k];
					cx[0] = j;
					rx[0] = i;
					zx[0] = k;
				}
            }
        }
    }
}


void local_max_clust_3D(float* im_vals, unsigned short* local_max_vals, unsigned short* bImg, unsigned short* out1, int r, int c, int z, int scale_xy, int scale_z)
{  
	//im_vals is the Laplacian of Gaussian
	//local_max_vals is the seed points (local maximum) with foreground seeds assigned an id > 0 and background seeds id == -1
	// out1 will contain the clustering output

    int*** max_nghbr_im;
    int min_r, min_c, max_r, max_c, min_z, max_z;
    
	//create max_nghbr_im and initialize it with its index (node) value
	max_nghbr_im = (int ***) malloc(r*sizeof(int**)); 
    for(int i=0; i<r; i++)
    {        
        max_nghbr_im[i] = (int **) malloc(c*sizeof(int*));
        for(int j=0; j<c; j++)
        {			
			max_nghbr_im[i][j] = (int *) malloc(z*sizeof(int));
			for(int k=0; k<z; k++)
			{				
				max_nghbr_im[i][j][k] = (k*r*c)+(i*c)+j;//LMX;
            }
        }
    }
    
    cerr << "max_nghbr_im initialized" << endl;

	//In this loop we look in a local region around each point and find the maximum value in the LoG image
	//Set the value to the index of the local maximum, (so if I am a seed point do nothing).
    for(int i=0; i<r; i++)
    {
        for(int j=0; j<c; j++)
        {
			for(int k=0; k<z; k++)
			{        				
				min_r = (int) max((double)(0.0),(double)(i-scale_xy));
				min_c = (int) max((double)(0.0),(double)(j-scale_xy));
				min_z = (int) max((double)(0.0),(double)(k-scale_z));
				max_r = (int) min((double)(r-1),(double)(i+scale_xy));
				max_c = (int) min((double)(c-1),(double)(j+scale_xy));                         
				max_z = (int) min((double)(z-1),(double)(k+scale_z)); 
                            
				//Check for seed point(local maximum)
				if(local_max_vals[(k*r*c)+(i*c)+j] !=0)//local_max_im[i][j][k]!=0)
				{					
                    continue;					
				}
				else
				{            
					int R, C, Z;
					//find maximum value in the LoG withing the search region (min_r->max_r,min_c->max_c,min_z->max_z).
					//R,C,Z will contain the coordinates of this local maximum value. 
					//r,c,z are the image dimensions.

					//Do not comment this line again please... I resoved the warning issue
					get_maximum(im_vals, min_r, max_r, min_c, max_c, min_z, max_z, &R, &C, &Z, r, c, z);                                              
                                                                                 
					/*double ind;
            
					if(max_nghbr_im[R][C][Z]>0)
						ind = max_nghbr_im[R][C][Z];
					else
						ind = (Z*r*c)+(C*r)+R;
            
					max_nghbr_im[i][j][k] = ind;*/	
					
					max_nghbr_im[i][j][k] = max_nghbr_im[R][C][Z];
				}
			}
        }		
    }
    
	
    int change = 1;
    double LM;
	cerr << "Entering main Clustering Loop" << endl;
	//Now continue to update until no more changes occur, eventually will have clusters pointing to seeds	
	int iterr = 0;
    while(change)	
    {			
		//For now, limit it to a maximum of 10 iterations
		iterr++;
		if(iterr == 10)
			break;		
        change=0;
		
        for(int i=0; i<r; i++)
        {
            for(int j=0; j<c; j++)
            {
				for(int k=0; k<z; k++)
				{                   
					LM = max_nghbr_im[i][j][k];
					if(LM==0)
						continue;													    
					
					//Calculate coordinates of local maximum based on its index
					int rem = ((long)LM) % (r*c);
					int Z = ((int)LM-rem) / (r*c); 
				    int C = ((long)rem) % c;
					int R = (rem-C)/c;
					
                
					//Check for seed (already a local max value) or connected to seed (my local maximum is a seed point)					
					if(local_max_vals[(k*r*c)+(i*c)+j] !=0 || local_max_vals[(Z*r*c)+(R*c)+C]!=0 ) 
						continue;
					else
					{
						change=change+1;
						max_nghbr_im[i][j][k]=max_nghbr_im[R][C][Z];
					}
				}
            }
        }
		cerr<<"change="<<change<<endl;
    }
    
	//cerr << "Preparing Output" << endl;
    
    for(int i=0; i<r; i++)
    {        
        for(int j=0; j<c; j++)
        {
			for(int k=0; k<z; k++)
            {
                LM = max_nghbr_im[i][j][k];
								    
                //if(local_max_vals[(int)LM] == -1 || bImg[(k*r*c)+(i*c)+j]==0)
				if(local_max_vals[(int)LM] == 65535 || bImg[(k*r*c)+(i*c)+j]==0)					
                    out1[(k*r*c)+(i*c)+j] = 0;
                else
				{
					//modified by Yousef on 8/21/2009
					//if the distance between me and my seed is more than a threshold.. then ignore me
					/*int rem = ((long)LM) % (r*c);
					int Z = (LM-rem) / (r*c); 
				    int C = ((long)rem) % c;
					int R = (rem-C)/c;
					double d = (i-R)*(i-R) + (j-C)*(j-C) + 3*(k-Z)*(k-Z);
					d = sqrt(d);
					if(d>10)
						out1[(k*r*c)+(i*c)+j] = 0;
					else*/
						out1[(k*r*c)+(i*c)+j] =local_max_vals[(int)LM];
				}
            }
        }
    }    

    
    for(int i=0; i<r; i++)
    {        
        for(int j=0; j<c; j++)
        {			
			free(max_nghbr_im[i][j]);
        }
		free(max_nghbr_im[i]);
    }
	free(max_nghbr_im);
}
 
