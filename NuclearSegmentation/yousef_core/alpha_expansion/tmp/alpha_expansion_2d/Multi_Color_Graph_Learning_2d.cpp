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

/*This function will be used for 4-color graph learning
 *The output of this function should be used in graph cuts (alpha expansions)
 */

//#include <vector>
//#include "gvc.h"
//#include <iostream>
//#include <math.h>
//#include "mex.h"
#include "Multi_Color_Graph_Learning_2D.h"

using namespace std;

//A 2-D multivariate gaussian
double Multivar_Norm(double X, double Y, double Ux, double Uy, double S00, double S01, double S11)
{
    double det_segma = (S00*S11)-(S01*S01);
    double Sinv00 = S11/det_segma;
    double Sinv01 = -S01/det_segma;
    double Sinv11 = S00/det_segma;
    double Mah_Dist = Sinv11*pow(Y-Uy,2)+2*Sinv01*(Y-Uy)*(X-Ux)+Sinv00*pow(X-Ux,2);
    double Ex = exp(-.5*Mah_Dist);
    double A = 1/((2*3.1416)*sqrt(det_segma));
    
    return (A*Ex);    
}

//void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
float* multiColGraphLearning(float* X_vals, int* labs_vals, int r, int c, int *NC)
{  
    //Initialize variables
    int** labs_im;
    int max_lab,ncolors;
    int** RAG;
    float* out;
    int* ColorOut;
    double** U;
    double*** Segma;
    double* P_I;
    double** Z;
    int* Z_sum;
   	
    //create a 2-D image that will hold the edge points only
    labs_im = (int **) malloc(r*sizeof(int*)); 
    max_lab = 0;    
    for(int i=0; i<r-1; i++)
    {
        labs_im[i] = (int *) malloc(c*sizeof(int));
        for(int j=0; j<c-1; j++)
        {		
            if(labs_vals[(i*c)+j] == 0)
            {
                labs_im[i][j] = 0;
                continue;
            }			
            if(labs_vals[(i*c)+j]!=labs_vals[((i+1)*c)+j] || labs_vals[(i*c)+j]!=labs_vals[(i*c)+(j+1)] || labs_vals[(i*c)+j]!=labs_vals[((i+1)*c)+(j+1)])
            {
                labs_im[i][j] = (int)labs_vals[(i*c)+j];
                if(labs_im[i][j] > max_lab)
                    max_lab = (int)labs_im[i][j];
            }
            else
                labs_im[i][j] = 0;
             
        }
    }
	
	
	labs_im[r-1] = (int *) malloc(c*sizeof(int));
	for(int j=0; j<c; j++)
    {        
		labs_im[r-1][j] = 0;
        if(labs_vals[((r-1)*c)+j]>max_lab)
            max_lab=labs_vals[((r-1)*c)+j];
    }
	
	for(int i=0; i<r; i++)
    {
		labs_im[i][c-1] = 0;
         if(labs_vals[(i*c)+(c-1)]>max_lab)
            max_lab=labs_vals[(i*c)+(c-1)];
    }
		
    //Build the region adjacency graph
	//fprintf(fid,"Build the region adjacency graph..");
    RAG = (int **) malloc(max_lab*sizeof(int*));
    for(int i=0; i<max_lab; i++)
    {
        RAG[i] = (int *) malloc(max_lab*sizeof(int));
        for(int j=0; j<max_lab; j++)
            RAG[i][j] = 0;
    }
    
	int L, L1, L2, L3;
    for(int i=0; i<r-1; i++)
    {        
        for(int j=0; j<c-1; j++)
        {			
            L = labs_im[i][j];
            if( L == 0)
                continue;
            else
            {
				L1 = labs_vals[((i+1)*c)+j];//labs_vals[(j*r)+(i+1)];
				L2 = labs_vals[(i*c)+(j+1)];//labs_vals[((j+1)*r)+i];
				L3 = labs_vals[((i+1)*c)+(j+1)];//labs_vals[((j+1)*r)+(i+1)];
				if((labs_im[i][j]!=L1)&&L1!=0)
                    RAG[L-1][L1-1] = RAG[L1-1][L-1] = 1;
                if((labs_im[i][j]!=L2)&&L2!=0)
                    RAG[L-1][L2-1] = RAG[L2-1][L-1] = 1;
                if((labs_im[i][j]!=L3)&&L3!=0)
                    RAG[L-1][L3-1] = RAG[L3-1][L-1] = 1;
            }
                        
        }
        free(labs_im[i]);        
    }
    free(labs_im);    
    
    
    //copy the RAG into an std vector of vectors
	std::vector< std::vector<int> > MAP;
	MAP.resize(max_lab);
    std::vector< std::vector<int> > MAP2;
	MAP2.resize(max_lab);
	ColorOut = (int *) malloc(max_lab*sizeof(int));
	for(int i=0; i<max_lab; i++)
	{	
		ColorOut[i] = 0;
		for(int j=0; j<max_lab; j++)
		{
			if(RAG[i][j]==1)
            {
				MAP[i].push_back(j+1);            
                MAP2[i].push_back(j+1);            
            }
		}
        free(RAG[i]);
	}
    free(RAG);
	
    //Added by Yousef on 4-20-2008: Add the second layer of neighbors, 
    //i.e. neighbors of my neighbors are also my neighbors.
    //I added that for situations when two alphas are seperated by just one cell
    //and expanding both alphas will result in merging them if the cell in between
    //is a small one.
    int NG1, NG2;
    for(int i=0; i<max_lab; i++)
	{	        
		for(int j=0; j<MAP2[i].size() ; j++)
		{
            NG1 = MAP2[i][j];
			//for each jth neighbor of the ith cell
            for(int k=0; k<MAP2[NG1-1].size() ; k++)
            {
                //add each kth neighbor of the jth cell to the ith cell
                NG2 = MAP2[NG1-1][k];
                if(NG2!=(i+1))
                    MAP[i].push_back(NG2);
            }				            
		}
	}
    
    //start the graph coloring using Sumit's sequential coloring code
    GVC* Gcol = new GVC();
 	ncolors = (int) NC[0];	
 	Gcol->sequential_coloring(max_lab,  ncolors, ColorOut, MAP );
     //Now, the resulting number of colors could be less than ncolors, so update ncolors
     int mx_col = 1;
     for(int i=0; i<max_lab; i++)
     {
         int cc = ColorOut[i]+1;
         if(cc>mx_col)
             mx_col = cc;
     }
     ncolors = mx_col;
	 NC[0] = ncolors;        
    
    //Now let's do the learning step
    //here, each pixel will have a probability of belonging to each one of the 
    //colors 
    //To do that, we assume that each object (cell) can be represented by a gaussian model
	Segma = (double ***) malloc(max_lab*sizeof(double**));
    U = (double **) malloc(max_lab*sizeof(double*));
    P_I = (double *) malloc(max_lab*sizeof(double));
    Z_sum = (int *) malloc(max_lab* sizeof(int));    
    for(int i=0; i<max_lab; i++)
    {
        Z_sum[i] = 0;
        P_I[i] = 0;
        U[i] = (double *) malloc(2*sizeof(double));
        Segma[i] = (double **) malloc(2*sizeof(double*));
        for(int j=0; j<2; j++)
        {
            U[i][j] = 0;
            Segma[i][j] = (double *) malloc(2*sizeof(double));
            for(int k=0; k<2; k++)
                Segma[i][j][k] = 0;
        }
    }
           
    //Now, compute the estimates of the parameters
    //Start with the means
    double val;
    int label;
    double ZS = 0;
    for(int i=0; i<r; i++)
    {
        for(int j=0; j<c; j++)
        {
            label = labs_vals[(i*c)+j]; 
            if(label == 0 || label>max_lab)
                continue;
            val = 1;//X_vals[(i*c)+j]; //Note,, for this function to work, LPG_im should be normalized between 0 and the maximum repetition		
            U[label-1][0] += (j*val);
            U[label-1][1] += (i*val);
            Z_sum[label-1]+= (int)val;
            ZS+=val;
        }
    }
           
    for(int i=0; i<max_lab; i++)
    {                    
        U[i][0] = (U[i][0] / (double)Z_sum[i]);
        U[i][1] = (U[i][1] / (double)Z_sum[i]);
        P_I[i] = (double)Z_sum[i]/ZS;
    }
    
    //Then the convariance matrices
    for(int i=0; i<r; i++)
    {
        for(int j=0; j<c; j++)
        {
            label = labs_vals[(i*c)+j]; 
            if(label == 0 || label>max_lab)
                continue;
			val = 1;//X_vals[(i*c)+j];
            Segma[label-1][0][0] += ((j-U[label-1][0])*(j-U[label-1][0])*val);
            Segma[label-1][0][1] += ((j-U[label-1][0])*(i-U[label-1][1])*val);
            Segma[label-1][1][0] += ((i-U[label-1][1])*(j-U[label-1][0])*val);
            Segma[label-1][1][1] += ((i-U[label-1][1])*(i-U[label-1][1])*val);            
        }
    }
    
    for(int i=0; i<max_lab; i++)
    {
        Segma[i][0][0] /= (double)Z_sum[i];
        Segma[i][0][1] /= (double)Z_sum[i];
        Segma[i][1][0] = Segma[i][0][1];
        Segma[i][1][1] /= (double)Z_sum[i];
    }
	
    
    //Now, for each pixel in the image, compute its probability to belong to each on of the classes (ideally 4 classes)
    //And set the data term as -ln(Pr)    
	out = (float *) malloc(r*c*(ncolors+1)*sizeof(float));
    double *Pr = (double *) malloc((ncolors+1)*sizeof(double));
    int C, cl, V;
    double P;
    int intst;
    for(int i=0; i<r; i++)
    {
        //Dterm[i] = (double **) malloc(c*sizeof(double*));
        for(int j=0; j<c; j++)
        {
            //Dterm[i][j] = (double *) malloc((ncolors+1)* sizeof(double));
            val = labs_vals[(i*c)+j];              
            for(cl=/*1*/0; cl<ncolors+1;cl++)
                Pr[cl] = /*F_H[intst];*/0;
            if(val == 0 || val>max_lab)
                Pr[0] = 1;				
            else
            {
                val = val-1;
                //check the prob for this point to be in its current object
                //and all the adjacent objects
                C = ColorOut[(int) val];				
                Pr[C+1] /*+*/= /*P_I[(int)val]*/Multivar_Norm(j, i, U[(int)val][0], U[(int)val][1], Segma[(int)val][0][0], Segma[(int)val][0][1], Segma[(int)val][1][1]);
                for (int k = 0; k < MAP[(int)val].size(); k++)
                {
                    V = MAP[(int)val][k]-1;
                    C = ColorOut[V];
                    P = /*P_I[V]*/Multivar_Norm(j, i, U[V][0], U[V][1], Segma[V][0][0], Segma[V][0][1], Segma[V][1][1]);//+F_H[intst];
                    Pr[C+1] = max(Pr[C+1],P);
                }               
			}
			for(int cc=0; cc<ncolors+1; cc++)//[(x+y*width)*#labels + l]				
				out[(j+i*c)*(ncolors+1) + cc] = min(-log(Pr[cc]),100.0);                
        }
    }
   
	return out;
}
