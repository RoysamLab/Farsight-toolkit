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

//
#include "local_max_clust_3D.h"

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


void local_max_clust_3D(float* im_vals, int* local_max_vals, int* bImg, int* out1, int r, int c, int z, int scale_xy, int scale_z)
{  
	//im_vals is the Laplacian of Gaussian
	//local_max_vals is the seed points (local maximum) with foreground seeds assigned an id > 0 and background seeds id == -1
	// out1 will contain the clustering output

    double*** max_nghbr_im;
    int min_r, min_c, max_r, max_c, min_z, max_z;
    
	//create max_nghbr_im and initialize it with its index (node) value
	max_nghbr_im = (double ***) malloc(r*sizeof(double**)); 
    for(int i=0; i<r; i++)
    {        
        max_nghbr_im[i] = (double **) malloc(c*sizeof(double*));
        for(int j=0; j<c; j++)
        {			
			max_nghbr_im[i][j] = (double *) malloc(z*sizeof(double));
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
    while(change)
	//for (int g=1; g<4; g++)
    {			
		cerr<<"change="<<change<<endl;
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
					
					/*int rem = ((long)LM)%(r*c);
					int Z = (LM-rem)/(r*c);
					int R = ((long)rem)%r;
					int C = (rem-R)/r;*/
				    
					
					//Calculate coordinates of local maximum based on its index
					int rem = ((long)LM) % (r*c);
					int Z = (LM-rem) / (r*c); 
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
    }
    
	//cerr << "Preparing Output" << endl;
    
    for(int i=0; i<r; i++)
    {        
        for(int j=0; j<c; j++)
        {
			for(int k=0; k<z; k++)
            {
                LM = max_nghbr_im[i][j][k];
								    
                if(local_max_vals[(int)LM] == -1 || bImg[(k*r*c)+(i*c)+j]==0)
                    out1[(k*r*c)+(i*c)+j] = 0;
                else
                    out1[(k*r*c)+(i*c)+j] =local_max_vals[(int)LM];
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
 
