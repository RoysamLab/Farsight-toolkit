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

#include "Local_Max_Clust_2D.h"

float get_maximum(float** A, int r1, int r2, int c1, int c2, int* rx, int* cx)
{
    float mx = A[r1][c1];
	rx[0] = r1;
	cx[0] = c1;
    for(int i=r1; i<=r2; i++)
    {
        for(int j=c1; j<=c2; j++)
        {
            if(A[i][j]>=mx)
            {
                mx = A[i][j];
                rx[0] = i;
                cx[0] = j;
            }
        }
    }	
    return mx;
}

int local_max_clust_2D(float* im_vals, int r, int c, double scale, unsigned short* out1, int* seeds_x, int* seeds_y, int num_seeds, unsigned short* bin_image)
{		
    float** im;    
    float** local_max_im;
    float** max_nghbr_im;
    int min_r, min_c, max_r, max_c;
    
    im = (float **) malloc(r* sizeof(float*)); 
    local_max_im = (float **) malloc(r* sizeof(float*)); 
    max_nghbr_im = (float **) malloc(r* sizeof(float*)); 
    for(int i=0; i<r; i++)
    {
        im[i] = (float *) malloc(c* sizeof(float));
        local_max_im[i] = (float *) malloc(c* sizeof(float));
        max_nghbr_im[i] = (float *) malloc(c* sizeof(float));
        for(int j=0; j<c; j++)
        {
            im[i][j] = im_vals[(i*c)+j];			
			max_nghbr_im[i][j] = 0.0;			
            local_max_im[i][j] = 0.0;
            
        }
    }
    
	int s_ind = 0;
	for(int i=0; i<num_seeds; i++)
	{
		float LMX = (float)((seeds_y[i]*c)+seeds_x[i]);
		if(bin_image[(int)LMX] == 0)
			local_max_im[seeds_y[i]][seeds_x[i]] = -1;
		else
		{
			s_ind++;
			local_max_im[seeds_y[i]][seeds_x[i]] = float(s_ind);
		}
		max_nghbr_im[seeds_y[i]][seeds_x[i]] = LMX;
	}

    
	        
    for(int i=0; i<r; i++)
    {		
        for(int j=0; j<c; j++)
        {
			//if(im[i][j] == 0)
			//	continue;			

            min_r = (int) std::max(0.0,i-scale);
            min_c = (int) std::max(0.0,j-scale);
            max_r = (int) std::min((double)r-1,i+scale);
            max_c = (int) std::min((double)c-1,j+scale);          
       
            if(local_max_im[i][j]!=0)
            {
				continue;
                //for(int m1=min_r; m1<=max_r; m1++)
                //    for(int m2=min_c; m2<=max_c; m2++)
                //        max_nghbr_im[m1][m2] =(i*c)+j;
            }
            else
            {            
                int R, C;
                float mx = get_maximum(im, min_r, max_r, min_c, max_c, &R, &C);   				
            
                float ind;

				//try this
				if(mx == im[i][j])
				{
					im[i][j] +=1;
					R = i;
					C = j;
				}
            
                if(max_nghbr_im[R][C]>0)
                    ind = max_nghbr_im[R][C];
                else
                    ind = (float)((R*c)+C);
            
                max_nghbr_im[i][j] = ind;	
            }
        }		
    }

	
	//FILE* fid = fopen("clustering_loop.txt","w");
    int change = 1;
    float LM;
	//int count =0;
    //while(change)
	for(int k=0;k<10;k++)
    {		
        change=0;		
        for(int i=0; i<r; i++)
        {
            for(int j=0; j<c; j++)
            {
				if(bin_image[(i*c)+j] == 0)
					continue;
                
                LM = max_nghbr_im[i][j];
                if(LM==0)
                    continue;
                
				int C = ((long)LM)%c;
				int R = ((int)LM-C)/c;
				
                //if it is a local maxima, or has been assigned to a local maxima ignore it
                if(local_max_im[i][j] !=0 || local_max_im[R][C]!=0) 
                    continue;
                else
                {					
                    change=change+1;
                    max_nghbr_im[i][j]=max_nghbr_im[R][C];
                }
            }
        }		
		std::cerr<<"Change= " <<change<<std::endl;
    }
	    
        
    for(int i=0; i<r; i++)
    {        
        for(int j=0; j<c; j++)
        {			            
			int LMX = (int) max_nghbr_im[i][j];
			int C = ((long)LMX)%c;
			int R = (LMX-C)/c;
			if(local_max_im[R][C]>0 && bin_image[(i*c)+j]>0)
				out1[(i*c)+j] = (unsigned short)local_max_im[R][C];
			else
				out1[(i*c)+j] = 0;
        }
    }
 
	free(im);    
    free(local_max_im);
    free(max_nghbr_im);

	return 1;
}
                
                    
                 
