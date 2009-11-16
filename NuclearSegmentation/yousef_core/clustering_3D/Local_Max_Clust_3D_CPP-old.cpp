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

#include "local_max_clust_3D.h"

using namespace std;

double get_maximum(double* A, int r1, int r2, int c1, int c2, int z1, int z2, int* rx, int* cx, int* zx)
{
    //int mx = im_vals[(z1*r*c)+(c1*r)+r1];//A[r1][c1][z1];
	double mx = A[(z1*(*rx)*(*cx))+(c1*(*rx))+r1];//A[r1][c1][z1];
    rx[0] = r1;
    cx[0] = c1;
    for(int i=r1; i<=r2; i++)
    {
        for(int j=c1; j<=c2; j++)
        {
			for(int k=z1; k<=z2; k++)
			{
				if(A[(k*(*rx)*(*cx))+(j*(*rx))+i]/*A[i][j][k]*/>mx)
				{
					//mx = [(k*r*c)+(j*r)+i];//A[i][j][k];
					mx = A[(k*(*rx)*(*cx))+(j*(*rx))+i];//A[i][j][k];
					rx[0] = i;
					cx[0] = j;
					zx[0] = k;
				}
            }
        }
    }
    return mx;
}


void local_max_clust_3D(double* im_vals, int* local_max_vals, int* out1, int r, int c, int z, int scale_xy, int scale_z)
{  
    //double*** im;
    //double*** local_max_im;
    double*** max_nghbr_im;
    int min_r, min_c, max_r, max_c, min_z, max_z;
    
    //im = (double ***) malloc(r, sizeof(double**)); 
    //local_max_im = (double ***) malloc(r, sizeof(double**)); 
    
	//max_nghbr_im = (double ***) malloc(r, sizeof(double**)); 
	max_nghbr_im = (double ***) malloc(r*sizeof(double**)); 
    for(int i=0; i<r; i++)
    {
        //im[i] = (double **) malloc(c, sizeof(double*));
        //local_max_im[i] = (double **) malloc(c, sizeof(double*));
        max_nghbr_im[i] = (double **) malloc(c*sizeof(double*));
        for(int j=0; j<c; j++)
        {
			//im[i][j] = (double *) malloc(z, sizeof(double));
			//local_max_im[i][j] = (double *) malloc(z, sizeof(double));
			max_nghbr_im[i][j] = (double *) malloc(z*sizeof(double));
			for(int k=0; k<z; k++)
			{
				//im[i][j][k]=im_vals[(k*r*c)+(j*r)+i];
				//double LMX = local_max_vals[(k*r*c)+(j*r)+i];			
				//local_max_im[i][j][k] = local_max_vals[(k*r*c)+(j*r)+i];//LMX;
				max_nghbr_im[i][j][k] = (k*r*c)+(j*r)+i;//LMX;
            }
        }
    }
    
    

    for(int i=0; i<r; i++)
    {
        for(int j=0; j<c; j++)
        {
			for(int k=0; k<z; k++)
			{
                //if(im[i][j][k] == 0)
                //    continue;
				min_r = (int) max((double)(0.0),(double)(i-scale_xy));
				min_c = (int) max((double)(0.0),(double)(j-scale_xy));
				min_z = (int) max((double)(0.0),(double)(k-scale_z));
				max_r = (int) min((double)(r-1),(double)(i+scale_xy));
				max_c = (int) min((double)(c-1),(double)(j+scale_xy));                         
				max_z = (int) min((double)(z-1),(double)(k+scale_z)); 
                            
       
				if(local_max_vals[(k*r*c)+(j*r)+i] !=0)//local_max_im[i][j][k]!=0)
				{
                    continue;
					/*for(int m1=min_r; m1<=max_r; m1++)
                    {
						for(int m2=min_c; m2<=max_c; m2++)
                        {
							for(int m3=min_z; m3<=max_z; m3++)
                            {
                                //if(im[m1][m2][m3]==0)
                                //    max_nghbr_im[m1][m2][m3] = 0;
                                //else
                                    max_nghbr_im[m1][m2][m3] =max_nghbr_im[i][j][k];//local_max_im[i][j][k];
                            }
                        }
                    }*/
				}
				else
				{            
					int R, C, Z;
					double mx = get_maximum(im_vals, min_r, max_r, min_c, max_c, min_z, max_z, &R, &C, &Z);                                              
                    
                    //this while loop to solve the problem of a non-local maxima pixel assigned to itself
                    //this happens because of masking the image with the binarization
            
                    /*while(R==i && C==j && Z==k && im[i][j][k]>=0)
                    {
                        min_r = (int) max(0.0,(double)min_r-1);
                        min_c = (int) max(0.0,(double)min_c-1);
                        min_z = (int) max(0.0,(double)min_z-1);
                        max_r = (int)min((double)r-1,(double)max_r+1);
                        max_c = (int)min((double)c-1,(double)max_c+1);                         
                        max_z = (int)min((double)z-1,(double)max_z+1);                               
                        mx = get_maximum(im, min_r, max_r, min_c, max_c, min_z, max_z, &R, &C, &Z);                       
                    }
                    if(mx==0)
                        continue;*/
            
					double ind;
            
					if(max_nghbr_im[R][C][Z]>0)
						ind = max_nghbr_im[R][C][Z];
					else
						ind = (Z*r*c)+(C*r)+R;
            
					max_nghbr_im[i][j][k] = ind;	
				}
			}
        }		
    }
    
	
    //FILE* fid = fopen("LocalmaxStat.txt","w");
    int change = 1;
    double LM;
	cerr << "Entering main Clustering Loop" << endl;
    while(change)
    {
        //fprintf(fid,"%d",change);
		std::cout<<"change="<<change<<std::endl;
        change=0;
		
        for(int i=0; i<r; i++)
        {
            for(int j=0; j<c; j++)
            {
				for(int k=0; k<z; k++)
				{
                    //if(im[i][j][k] == 0)
                    //    continue;
					//if it is a local maxima, or has been assigned to a local maxima ignore it
					LM = max_nghbr_im[i][j][k];
					if(LM==0)
						continue;

					int rem = ((long)LM)%(r*c);
					int Z = (LM-rem)/(r*c);
					int R = ((long)rem)%r;
					int C = (rem-R)/r;
				
                
					if(local_max_vals[(k*r*c)+(j*r)+i] !=0 /*local_max_im[i][j][k] !=0*/ || local_max_vals[(Z*r*c)+(C*r)+R]!=0 /*local_max_im[R][C][Z]!=0*/) 
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
    //fclose(fid);
    //Correct the IDs of the seeds
    /*int obj_ind = 0;
    for(int i=0; i<r; i++)
    {        
        for(int j=0; j<c; j++)
        {
			for(int k=0; k<z; k++)
            {
                if(local_max_im[i][j][k] > 0)
                {
                    obj_ind++;
                    local_max_im[i][j][k] = obj_ind;
                }
            }
				
        }
    }*/    
    
    
    for(int i=0; i<r; i++)
    {        
        for(int j=0; j<c; j++)
        {
			for(int k=0; k<z; k++)
            {
                LM = max_nghbr_im[i][j][k];
					
			    /*int rem = ((long)LM)%(r*c);
			    int Z = (LM-rem)/(r*c);
			    int R = ((long)rem)%r;
				int C = (rem-R)/r;*/
                if(local_max_vals[(int)LM] == -1)
                    out1[(k*r*c)+(j*r)+i] = 0;
                else
                    out1[(k*r*c)+(j*r)+i] =local_max_vals[(int)LM];
            }
        }
    }    
    
}
         
