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
#include "cell_binarization.h"
#include <limits.h>

#ifndef SIZE_MAX
#define SIZE_MAX ((size_t)-1)
#endif

using namespace std;

//Main function for 2-D binarization
int Cell_Binarization_2D(unsigned char* imgIn, unsigned short *imgOut, int R, int C, int shd)
{			
	//Now, to do the binarization, follow these steps:
	//1- Assuming that the histogram of the image is modeled by a mixture of two 
	//poisson distributions, estimate the parameters of the mixture model 
	float alpha_B, alpha_F, P_I;
	alpha_B = alpha_F = P_I = 0;	
	//MinErrorThresholding(imgIn, &alpha_B, &alpha_F, &P_I, R, C, 1, shd,imgOut); 		
	float alpha_C, P_I2;
	alpha_C = P_I2 = 0.0;
	threeLevelMinErrorThresh(imgIn, &alpha_B, &alpha_F, &alpha_C, &P_I, &P_I2, R, C, 1);

	//2- Apply binarization refinement using graph-cuts	
	int n_nodes, n_edges;
	//alpha_C = -1;
	if(alpha_C == -1)//this is a two levels case
	{
		//graph-cuts
		Seg_GC_Full_2D(imgIn, R, C, alpha_F, alpha_B, P_I, &n_nodes, &n_edges, imgOut);
	}
	else  //three levels case
	{
		//try this
		alpha_C = (alpha_C+alpha_F)/2;
		alpha_B = (alpha_B+alpha_F)/2;
		//graph-cuts
		Seg_GC_Full_2D_Three_Level(imgIn, R, C, alpha_C,alpha_B, alpha_F, P_I, P_I2,&n_nodes, &n_edges, imgOut);
	}		
	return 1;	
}

//Main function for 3-D binarization
int Cell_Binarization_3D(unsigned char *imgIn, unsigned short* imgOut, int R, int C, int Z, int shd, int div)
{			
	//Now, to do the binarization, follow these steps:
	//1- Assuming that the histogram of the image is modeled by a mixture of two 
	//poisson distributions, estimate the parameters of the mixture model 
	float alpha_B, alpha_F, P_I;
	alpha_B = alpha_F = P_I = 0;	
	//CompMixPoss3D(imgIn, &alpha_B, &alpha_F, &P_I, R, C, Z, shd); 	
	MinErrorThresholding(imgIn, &alpha_B, &alpha_F, &P_I, R, C, Z, shd,imgOut); 	
	
	//Modified by Isaac on 1-22-2010: Check memory before choosing divisor, and make smarter block divisions!!
	//120 bytes are needed for each pixel to create a graph!!!!
	int block_divisor = 1;
	if(div)	//means I'm going to try to find best number of blocks
	{
		bool done = false;
		while(!done)		//Make sure my requested size is less than maximum for machine:
		{
			double totalPix = (double)(R*C*Z) / ((double)block_divisor*(double)block_divisor);
			//120 bytes are needed for each pixel to create a graph!!!!
			unsigned long max = (unsigned long)SIZE_MAX;
			if( 120.0*totalPix > (double)max )
			{
				std::cerr << "Increasing divisor" << std::endl;
				block_divisor = block_divisor*2;
			}
			else
			{
				done = true;
			}
		}
		done = false;
		while(!done)	//Now make sure I can allocate memory
		{
			unsigned long totalPix = (unsigned long)(R*C*Z) / ((unsigned long)block_divisor*(unsigned long)block_divisor);
			//120 bytes are needed for each pixel to create a graph!!!!
			unsigned char *tmpp = (unsigned char*)malloc(120*totalPix);
			if(!tmpp)			//unAble to allocate
			{
				std::cerr << "Increasing divisor" << std::endl;
				block_divisor = block_divisor*2;
			}
			else
			{
				done = true;
			}

			free(tmpp); //delete it
			tmpp = NULL;
		}
	}
	
	std::cerr << "block_divisor = " << block_divisor << std::endl;

	int *subImgBlock = new int[6];//[x1,y1,z1,x2,y2,z1]	
	subImgBlock[4] = 0;
	subImgBlock[5] = Z;
	int blk = 1;
	int cntr = 0;
	for(int i=0; i<R; i+=R/block_divisor)
		for(int j=0; j<C; j+=C/block_divisor)
			cntr++;

	for(int i=0; i<R; i+=R/block_divisor)
	{
		for(int j=0; j<C; j+=C/block_divisor)
		{			
			std::cout<<"    Binarizing block "<<blk<<" of "<<cntr<<std::endl;
			subImgBlock[0] = j;
			subImgBlock[1] = (int)j+C/block_divisor+1;
			subImgBlock[2] = i;
			subImgBlock[3] = (int)i+R/block_divisor+1;
			if(subImgBlock[1] > C)
				subImgBlock[1] = C;
			if(subImgBlock[3] > R)
				subImgBlock[3] = R;

			Seg_GC_Full_3D_Blocks(imgIn, R, C, Z, alpha_F, alpha_B, P_I, imgOut,subImgBlock);
			blk++;
		}
	}
	delete [] subImgBlock;

	return 1;//num_objects;	
}

int Neuron_Binarization_3D(unsigned char *imgIn, unsigned short* imgOut, int R, int C, int Z, int shd, int div) //modifed by Yousef on 5-20-2008.. The first input change from uchar* to int*
{			
	float alpha_B, alpha_F, P_I;
	alpha_B = alpha_F = P_I = 0;		
	MinErrorThresholding(imgIn, &alpha_B, &alpha_F, &P_I, R, C, Z, shd,imgOut); 	
	
	int block_divisor = 1;
	int tmp1;
	
	if (R>128||C>128)
	{
       tmp1 = max(R,C);
	   block_divisor = tmp1/96;
	}
         
	std::cout << "Blocked-Graph cuts segmentation " << std::endl;

	int *subImgBlock = new int[6];//[x1,y1,z1,x2,y2,z1]	
	subImgBlock[4] = 0;
	subImgBlock[5] = Z;
	int blk = 1;
	int cntr = 0;
	for(int i=0; i<R; i+=R/block_divisor)
		for(int j=0; j<C; j+=C/block_divisor)
			cntr++;

	for(int i=0; i<R; i+=R/block_divisor)
	{
		for(int j=0; j<C; j+=C/block_divisor)
		{			
			// std::cout<<"    Binarizing block "<<blk<<" of "<<cntr<<std::endl;
			subImgBlock[0] = j;
			subImgBlock[1] = (int)j+C/block_divisor+1;
			subImgBlock[2] = i;
			subImgBlock[3] = (int)i+R/block_divisor+1;
			if(subImgBlock[1] > C)
				subImgBlock[1] = C;
			if(subImgBlock[3] > R)
				subImgBlock[3] = R;

			Seg_GC_Full_3D_Blocks(imgIn, R, C, Z, alpha_F, alpha_B, P_I, imgOut,subImgBlock);
			blk++;
		}
	}
			
	delete [] subImgBlock;

	return 1;//num_objects;	
}

double compute_poisson_prob(double intensity, double alpha)
{
    /*this is the equation
      P = (alpha^intensity)*exp(-alpha)/factorial(intensity);
      however, since alpha and the intensity could be large, computing P in that
      way will result in infinity values from (alpha^intensity) and
      factorial(intensity) as a result of matlab's limitations of the data types*/

    //here is the solution
    double A, P;
    A = exp(-alpha);
    P = 1;
    for (int i=1; i<= intensity; i++)
        P = P * (alpha/i);
    
    P = P*A;

	if(P < numeric_limits<long double>::epsilon())
		P = numeric_limits<long double>::epsilon();
    
    return P;
}

void Seg_GC_Full_2D(unsigned char* IM,
                    int r, 
                    int c, 
                    double alpha_F, 
                    double alpha_B, 
                    double P_I, 
                    int* num_nodes, 
                    int* num_edges, 
                    unsigned short* Seg_out)
{    
    int curr_node;
    int rght_node;
    int down_node;
    int diag_node;
    double Df;
    double Db;
    double Dr;
    double Dd; 
    double Dg; 
    double sig;
    double F_H[256];
    double B_H[256];
    typedef Graph_B<int,int,int> GraphType;
       
    //Set the number of edges and the number of nodes and open the files that
    //will be used to save the weights
    num_nodes[0] = r*c;
    num_edges[0] = 3*r*c-2*r-2*c+1;
        	
    
    //Before entering the loop, compute the poisson probs 
    for(int i=0; i<256; i++)
    {
		if(i>=alpha_F)
			F_H[i] = (1-P_I)*compute_poisson_prob((int)alpha_F,alpha_F);
		else
			F_H[i] = (1-P_I)*compute_poisson_prob(i,alpha_F);
		if(i<=alpha_B)
			B_H[i] = P_I*compute_poisson_prob(int(alpha_B),alpha_B);
		else
        B_H[i] = P_I*compute_poisson_prob(i,alpha_B);
    }
	
    //Here is the main loop.. 
    //For each point, compute the terminal and neighbor edge weights
	GraphType *g = new GraphType(/*estimated # of nodes*/ num_nodes[0], /*estimated # of edges*/ num_edges[0]); 
    for(int i=0; i<r; i++)
    {
        for(int j=0; j<c; j++)
        {			
			/*Get the terminal edges capacities and write them to a file
			These capacities represent penalties of assigning each point to
			either the fg or the bg. Now, I am using the distance of the point
			to the fg and bg. Then a short distance between a point and a class (fg or bg)
			means the penalty of the assignement is low and vice versa*/ 
								
		   //Isaac change on 5/9/08
		   //int intst = (int) IM[i][j];
			int intst = (int) IM[i*c + j];

            //Added by Yousef on jan 17, 2008
            //check if this is a seed point
            if(intst == 255)
            {
                Df = 0;
                Db = 1000;
            }
            else if(intst == 0)
            {
                Df = 1000;
                Db = 0;
            }
            else
            {                
                Df = -log(F_H[intst]);  //it was multiplied by .5          
                if(Df>1000.0)
                    Df = 1000;
                Db = -log(B_H[intst]);
                if(Db>1000.0)
                    Db=1000;
                
            }            
            
            curr_node = (i*c)+j; 
            
			g -> add_node();
			g -> add_tweights( curr_node,   /* capacities */ Df,Db);                            			
		}
	}

    sig = 15.0;
	double w = 2.0;
    for(int i=0; i<r-1; i++)
      {
      for(int j=0; j<c-1; j++)
        {			
        // get the neighbor edges capacities and write them to a file.
	// Now, each edge capacity between two neighbors p and q represent
	// the penalty for discontinuety. Since I am using the difference in
	// the intensities, I should take the inverse so that very similar
	// objects should have a large discontinuety penalty between them*/

		curr_node = (i*c)+j; 
        rght_node = curr_node+1;
        down_node = curr_node+c;
        diag_node = curr_node+c+1;
            
            
        //from Boykov's paper instead
		Dr = w*exp(-pow((double)IM[i*c + j]-(double)IM[i*c + j+1],2)/(2*pow(sig,2)));
		g -> add_edge( curr_node, rght_node,    /* capacities */  Dr, Dr );	
            
                   
        //from Boykov's paper instead
		Dd = w*exp(-pow((double)IM[i*c + j]-(double)IM[(i+1)*c + j],2)/(2*pow(sig,2)));
		g->add_edge( curr_node, down_node,    /* capacities */  Dd, Dd );
                 
		Dg = w*exp(-pow((double)IM[i*c + j]-(double)IM[(i+1)*c + (j+1)],2)/(2*pow(sig,2)));
		g->add_edge( curr_node, diag_node,    /* capacities */  Dg, Dg );            
               
        }
      }    
	
    //Compute the maximum flow:
	g->maxflow();		//Alex DO NOT REMOVE


	int RR,CC;
	for(int i=0; i<num_nodes[0]; i++)
	{
		CC = ((long)i)%c;
        RR = (i-CC)/c;
		if(g->what_segment(i) == GraphType::SOURCE)
			Seg_out[RR*c + CC]=0;
		else
			Seg_out[RR*c + CC]=255;
	}

	delete g;
}


void Seg_GC_Full_3D(unsigned char* IM, int r, int c, int z, double alpha_F, double alpha_B, double P_I, int* Seg_out)
{   
    int curr_node, nbr_node;
    double Df, Db, Dn, sig, w;
    double F_H[256], B_H[256];
	double intst, intst2;
    unsigned long int num_nodes, num_edges;
	typedef Graph_B<short,short,short> GraphType;  
      
    //Set the number of edges and the number of nodes and open the files that
    //will be used to save the weights
    num_nodes = r*c*z;
    num_edges = (r-1)*(c-1)*(z-1)*3;   
        
    //Before entering the loop, compute the poisson probs    
	for(int i=0; i<256; i++)
    {
		if(i>=alpha_F)
			F_H[i] = (1-P_I)*compute_poisson_prob((int)alpha_F,alpha_F);
		else
			F_H[i] = (1-P_I)*compute_poisson_prob(i,alpha_F);
		if(i<=alpha_B)
			B_H[i] = P_I*compute_poisson_prob(int(alpha_B),alpha_B);
		else
        B_H[i] = P_I*compute_poisson_prob(i,alpha_B);
    }
	//std::cerr << "Poisson Probabilities Computed" << std::endl;
    
	//Construct the graph
	GraphType *g = new GraphType(/*estimated # of nodes*/ num_nodes, /*estimated # of edges*/ num_edges); 

	//std::cerr << "Graph Memory Allocated" << std::endl;

	//Here is the main loop.. 
    //For each point, compute the terminal and neighbor edge weights
	for(int k=0; k<z; k++)
    {		
        for(int j=0; j<r; j++)
        {			
			for(int i=0; i<c; i++)
			{
				
				/*Get the terminal edges capacities and write them to a file
				These capacities represent penalties of assigning each point to
				either the fg or the bg. Now, I am using the distance of the point
				to the fg and bg. Then a short distance between a point and a class (fg or bg)
				means the penalty of the assignement is low and vice versa*/ 
								

				//unsigned long int index = mxCalcSingleSubscript(prhs[0], nsubs, subs);
				//intst = (mxGetPr(prhs[0]))[index];
				// CHANGE BY ISAAC 4/14/08
				
				//CHANGE BY ISAAC 4/17/08
				//curr_node = (k*r*c)+(i*c)+j; 
				curr_node = (k*r*c)+(j*c)+i;

				int intst = (int) IM[curr_node];	//Changed 5/9/08


				//Added by Yousef on jan 17, 2008
				//check if this is a seed point
				        
					Df = -log(F_H[(int)intst]);  //it was multiplied by .5                              
					if(Df>1000.0)
						Df = 1000;
					Db = -log(B_H[(int)intst]);                    
					if(Db>1000.0)
						Db=1000;     			
			        				

				g -> add_node();
				
				g -> add_tweights( curr_node,   /* capacities */ Df,Db);                         
			}
		}
	}

	//std::cout << "First Loop Complete" << std::endl;
	
	sig = 30.0;
	w=10.0;
	for(int k=0; k<z; k++)
    {		
        for(int j=0; j<r; j++)
        {			
			for(int i=0; i<c; i++)
			{
				/*get the neighbor edges capacities and write them to a file.
				Now, each edge capacity between two neighbors p and q represent
				the penalty for discontinuety. Since I am using the difference in
				the intensities, I should take the inverse so that very similar
				objects should have a large discontinuety penalty between them*/

				if(i>=c-1 || j>=r-1 || k>=z-1)
					continue;

				//Try this: Just define 3 edges to the neigbors				
				curr_node = (k*r*c)+(j*c)+i;
				intst = (int) IM[curr_node];
				
				//1.
				nbr_node = (k*r*c)+(j*c)+(i+1);				
				//if(Seg_out[nbr_node] == Seg_out[nbr_node])
				//	Dn = 0;
				//else
				//{
					intst2 = (int) IM[nbr_node];
					Dn = w*exp(-pow((double)intst-intst2,2)/(2*pow(sig,2)));				
				//}
				g -> add_edge( curr_node, nbr_node,    /* capacities */  Dn, Dn );
				//2.				
				nbr_node = (k*r*c)+((j+1)*c)+i;
				//if(Seg_out[nbr_node] == Seg_out[nbr_node])
				//	Dn = 0;
				//else
				//{
					intst2 = (int) IM[nbr_node];
					Dn = w*exp(-pow((double)intst-intst2,2)/(2*pow(sig,2)));				
				//}												
				g -> add_edge( curr_node, nbr_node,    /* capacities */  Dn, Dn );
				//3.				
				nbr_node = ((k+1)*r*c)+(j*c)+i;
				//if(Seg_out[nbr_node] == Seg_out[nbr_node])
				//	Dn = 0;
				//else
				//{
					intst2 = (int) IM[nbr_node];
					Dn = w*exp(-pow((double)intst-intst2,2)/(2*pow(sig,2)));				
				//}												
				g -> add_edge( curr_node, nbr_node,    /* capacities */  Dn, Dn );
			}
		}
	}    
	
	//std::cout << "Second Loop Complete" << std::endl;
	
	//Compute the maximum flow:
	g->maxflow();		//Alex DO NOT REMOVE
	
	for(int k=0; k<z; k++)
    {		
        for(int j=0; j<r; j++)
        {			
			for(int i=0; i<c; i++)
			{
				if(g->what_segment((k*r*c)+(j*c)+i) == GraphType::SOURCE)
					Seg_out[(k*r*c)+(j*c)+i]=0;
				else
					Seg_out[(k*r*c)+(j*c)+i]=255;
			}
		}
	}

	delete g;
}

void CompMixPoss(unsigned char* img, float* alpha_B, float* alpha_A, float* P_I, int R, int C, int shiftDown)
{
	//Copmute the normalized histogram of the image
	double Hist[256];
	long int K[256];
	double R1, R2, R3, a, b, c,r;
	for(int i=0; i<256; i++)
	{
		Hist[i] = 0;		
	}

	for(int i=0; i<R; i++)
	{
		for(int j=0; j<C; j++)
		{
			//Hist[img[i][j]]=Hist[img[i][j]]+1.0;	//BY ISAAC 4/9/08
			Hist[img[i*C + j]]=Hist[img[i*C + j]]+1.0;
		}
	
	}

	//Normalize the histogram and compute the first 3 moments
	
	R1=R2=R3=0;
	for(int i=0; i<256; i++)
	{
		Hist[i] = (Hist[i]/(R*C));// + eps;
		//fprintf(fid,"%f\n",Hist[i]);
		K[i] = i;
		R1 += K[i]*Hist[i];
		R2 += K[i]*K[i]*Hist[i];
		R3 += K[i]*K[i]*K[i]*Hist[i];
	}
	
	//compute the coefficients (a,b,c) or the quadratic eqn: a*Th^2+b*Th+c=0
	
	a = (R1*R1)+R1-R2;
	b = (R1*R1)-(R1*R2)+(2*R1)-(3*R2)+R3;
	c = (R2*R2)-(R1*R1)+(R1*R2)-(R1*R3);
	
	//Now, alpha_B and alpha_F are the roots of the above eqn:
	r = sqrt((b*b)-(4*a*c));
	alpha_B[0] = (-b+r)/(2*a);
	alpha_A[0] = (-b-r)/(2*a); 	

	//finally, compute the weighting coef. P_I:
	P_I[0] = (R1-alpha_A[0])/(alpha_B[0]-alpha_A[0]);

	
	//Some times you need to shift the means down. The next two lines are optional
	if(shiftDown == 1)
	{
		alpha_A[0] = max(alpha_A[0]/2,(alpha_B[0]+alpha_A[0])/2);
		alpha_B[0] = alpha_B[0]/1.5;
	}
	
}

void CompMixPoss3D(unsigned char* img, float* alpha_B, float* alpha_A, float* P_I, int R, int C, int Z, int shiftDown)
{
	//Copmute the normalized histogram of the image
	double Hist[256];
	long int K[256];
	double R1, R2, R3, a, b, c,r;
	for(int i=0; i<256; i++)
	{
		Hist[i] = 0;		
	}

	for(int i=0; i<R; i++)
	{
		for(int j=0; j<C; j++)
		{
			for(int k=0; k<Z; k++)
			{
				//Edit by Isaac 4/14/08
				//Hist[img[i][j][k]]=Hist[img[i][j][k]]+1.0;
				//Hist[img[k][i][j]]=Hist[img[k][i][j]]+1.0;	again 5/9/08
				Hist[img[(k*R*C)+(i*C)+j]]=Hist[img[(k*R*C)+(i*C)+j]]+1.0;
			}
		}
	
	}

	//Normalize the histogram and compute the first 3 moments
	
	R1=R2=R3=0;
	for(int i=0; i<256; i++)
	{
		Hist[i] = (Hist[i]/(R*C*Z));// + eps;
		//fprintf(fid,"%f\n",Hist[i]);
		K[i] = i;
		R1 += K[i]*Hist[i];
		R2 += K[i]*K[i]*Hist[i];
		R3 += K[i]*K[i]*K[i]*Hist[i];
	}
	
	//compute the coefficients (a,b,c) or the quadratic eqn: a*Th^2+b*Th+c=0
	
	a = (R1*R1)+R1-R2;
	b = (R1*R1)-(R1*R2)+(2*R1)-(3*R2)+R3;
	c = (R2*R2)-(R1*R1)+(R1*R2)-(R1*R3);
	
	//Now, alpha_B and alpha_F are the roots of the above eqn:
	r = sqrt((b*b)-(4*a*c));
	alpha_B[0] = (-b+r)/(2*a);
	alpha_A[0] = (-b-r)/(2*a); 	

	//finally, compute the weighting coef. P_I:
	P_I[0] = (R1-alpha_A[0])/(alpha_B[0]-alpha_A[0]);

	
	//Some times you need to shift the means down. The next two lines are optional
	if(shiftDown == 1)
	{
		alpha_A[0] = max(alpha_A[0]/2,(alpha_B[0]+alpha_A[0])/2);
		alpha_B[0] = alpha_B[0]/1.5;
	}
	
}

void MinErrorThresholding(unsigned char* img, float* alpha_B, float* alpha_A, float* P_I, int R, int C, int Z, int shiftDown, unsigned short *imgOut)
{
	typedef  short  InputPixelType;
	typedef  short  OutputPixelType;

	typedef itk::Image< InputPixelType,  3 >   InputImageType;
	typedef itk::Image< OutputPixelType, 3 >   OutputImageType;

	typedef itk::MinErrorThresholdImageFilter<
               InputImageType, OutputImageType >  FilterType;
  	//Create an itk image
	InputImageType::Pointer im;
	im = InputImageType::New();
	InputImageType::PointType origin;
    origin[0] = 0.; 
    origin[1] = 0.;    
	//if(Z>1)
		origin[2] = 0.;    
    im->SetOrigin( origin );

    InputImageType::IndexType start;
    start[0] =   0;  // first index on X
    start[1] =   0;  // first index on Y 
	//if(Z > 1)
		start[2] =   0;  // first index on Z 
    InputImageType::SizeType  size;
    size[0]  = C;  // size along X
    size[1]  = R;  // size along Y
	//if(Z > 1)
		size[2]  = Z;  // size along Z

    InputImageType::RegionType region;
    region.SetSize( size );
    region.SetIndex( start );
    
    im->SetRegions( region );
    im->Allocate();
    im->FillBuffer(0);
	im->Update();
	
	short minv = 255;
	short maxv = 0;
	//copy the input image into the ITK image
	typedef itk::ImageRegionIteratorWithIndex< InputImageType > IteratorType;
	IteratorType iterator1(im,im->GetRequestedRegion());
	int lng = R*C*Z;
	for(int i=0; i<lng; i++)
	{
		short val = (short)img[i];
		if (val > maxv)
			maxv = val;
		if (val < minv)
			minv = val;

		iterator1.Set((short)img[i]);
		++iterator1;	
	}

	//binarize
	FilterType::Pointer filter = FilterType::New();
	filter->SetInput( im );
	filter->SetNumberOfHistogramBins (128);
	filter->Update();

	//get output
	alpha_B[0] = ((float) filter->GetAlphaLeft());
	alpha_A[0] =  (float)filter->GetAlphaRight();
	P_I[0] = (float)filter->GetPriorLeft();
	
	//Some times you need to shift the means down. The next two lines are optional
	if(shiftDown == 1)
	{
		alpha_A[0] = max(alpha_A[0]/2,(alpha_B[0]+alpha_A[0])/2);
		alpha_B[0] = alpha_B[0]/1.5;
	}
	IteratorType iterate2(filter->GetOutput(),filter->GetOutput()->GetRequestedRegion());	
	for(int i=0; i<lng; i++)
	{				
		imgOut[i] = (int)iterate2.Get();
		++iterate2;			
	}
}


void Seg_GC_Full_3D_Blocks(unsigned char* IM, int r, int c, int z, double alpha_F, double alpha_B, double P_I, unsigned short* Seg_out, int* imBlock)
{   
    int curr_node, nbr_node;
    double Df, Db, Dn, sig, w;
    double F_H[256], B_H[256];
	double intst, intst2;
    unsigned long int num_nodes, num_edges;
	typedef Graph_B<short,short,short> GraphType;  
      
    //Set the number of edges and the number of nodes and open the files that
    //will be used to save the weights	
    num_nodes = (imBlock[1]-imBlock[0])*(imBlock[3]-imBlock[2])*(imBlock[5]-imBlock[4]);
    num_edges = (imBlock[1]-imBlock[0]-1)*(imBlock[3]-imBlock[2]-1)*(imBlock[5]-imBlock[4]-1)*3;   
        
    //Before entering the loop, compute the poisson probs    
	for(int i=0; i<256; i++)
    {
		if(i>=alpha_F)
			F_H[i] = (1-P_I)*compute_poisson_prob((int)alpha_F,alpha_F);
		else
			F_H[i] = (1-P_I)*compute_poisson_prob(i,alpha_F);
		if(i<=alpha_B)
			B_H[i] = P_I*compute_poisson_prob(int(alpha_B),alpha_B);
		else
			B_H[i] = P_I*compute_poisson_prob(i,alpha_B);
    }
	//std::cerr << "Poisson Probabilities Computed" << std::endl;
    
	//Construct the graph
	GraphType *g = new GraphType(/*estimated # of nodes*/ num_nodes, /*estimated # of edges*/ num_edges); 

	//std::cerr << "Graph Memory Allocated" << std::endl;

	//Here is the main loop.. 
    //For each point, compute the terminal and neighbor edge weights
	int IND = -1;
	for(int k=imBlock[4]; k<imBlock[5]; k++)
    {		
        for(int j=imBlock[2]; j<imBlock[3]; j++)
        {			
			for(int i=imBlock[0]; i<imBlock[1]; i++)
			{
				IND++;
				/*Get the terminal edges capacities and write them to a file
				These capacities represent penalties of assigning each point to
				either the fg or the bg. Now, I am using the distance of the point
				to the fg and bg. Then a short distance between a point and a class (fg or bg)
				means the penalty of the assignement is low and vice versa*/ 
								

				//unsigned long int index = mxCalcSingleSubscript(prhs[0], nsubs, subs);
				//intst = (mxGetPr(prhs[0]))[index];
				// CHANGE BY ISAAC 4/14/08
				
				//CHANGE BY ISAAC 4/17/08
				//curr_node = (k*r*c)+(i*c)+j; 
				curr_node = (k*r*c)+(j*c)+i;

				int intst = (int) IM[curr_node];	//Changed 5/9/08


				//Added by Yousef on jan 17, 2008
				//check if this is a seed point
				        
					Df = -log(F_H[(int)intst]);  //it was multiplied by .5                              
					if(Df>1000.0)
						Df = 1000;
					Db = -log(B_H[(int)intst]);                    
					if(Db>1000.0)
						Db=1000;     			
			        				

				g -> add_node();
				
				g -> add_tweights( IND,   /* capacities */ Df,Db);                         
			}
		}
	}

	//std::cout << "First Loop Complete" << std::endl;
	
	sig = 50.0;
	w=10.0;
	IND = -1;
	for(int k=imBlock[4]; k<imBlock[5]; k++)
    {		
        for(int j=imBlock[2]; j<imBlock[3]; j++)
        {			
			for(int i=imBlock[0]; i<imBlock[1]; i++)
			{
				IND++;
				/*get the neighbor edges capacities and write them to a file.
				Now, each edge capacity between two neighbors p and q represent
				the penalty for discontinuety. Since I am using the difference in
				the intensities, I should take the inverse so that very similar
				objects should have a large discontinuety penalty between them*/

				if(i>=imBlock[1]-1 || j>=imBlock[3]-1 || k>=imBlock[5]-1)
					continue;

				//Try this: Just define 3 edges to the neigbors				
				curr_node = (k*r*c)+(j*c)+i;
				intst = (int) IM[curr_node];
				
				//1.
				nbr_node = (k*r*c)+(j*c)+(i+1);				
				//if(Seg_out[nbr_node] == Seg_out[nbr_node])
				//	Dn = 0;
				//else
				//{
					intst2 = (int) IM[nbr_node];
					Dn = w*exp(-pow((double)intst-intst2,2)/(2*pow(sig,2)));				
				//}
				g -> add_edge( IND, IND+1,    /* capacities */  Dn, Dn );
				//2.				
				nbr_node = (k*r*c)+((j+1)*c)+i;
				//if(Seg_out[nbr_node] == Seg_out[nbr_node])
				//	Dn = 0;
				//else
				//{
					intst2 = (int) IM[nbr_node];
					Dn = w*exp(-pow((double)intst-intst2,2)/(2*pow(sig,2)));				
				//}												
				g -> add_edge( IND, IND+(imBlock[1]-imBlock[0]),    /* capacities */  Dn, Dn );
				//3.				
				nbr_node = ((k+1)*r*c)+(j*c)+i;
				//if(Seg_out[nbr_node] == Seg_out[nbr_node])
				//	Dn = 0;
				//else
				//{
					intst2 = (int) IM[nbr_node];
					Dn = w*exp(-pow((double)intst-intst2,2)/(2*pow(sig,2)));				
				//}												
				g -> add_edge( IND, IND+(imBlock[1]-imBlock[0])*(imBlock[3]-imBlock[2]),    /* capacities */  Dn, Dn );
			}
		}
	}    
	
	//std::cout << "Second Loop Complete" << std::endl;
	
	//Compute the maximum flow:
	g->maxflow();		//Alex DO NOT REMOVE
	
	IND = -1;
	for(int k=imBlock[4]; k<imBlock[5]; k++)
    {		
        for(int j=imBlock[2]; j<imBlock[3]; j++)
        {			
			for(int i=imBlock[0]; i<imBlock[1]; i++)
			{
				IND++;
				if(g->what_segment(IND) == GraphType::SOURCE)
					Seg_out[(k*r*c)+(j*c)+i]=0;
				else
					Seg_out[(k*r*c)+(j*c)+i]=255;
			}
		}
	}

	delete g;
}

void threeLevelMinErrorThresh(unsigned char* im, float* Alpha1, float* Alpha2, float* Alpha3, float* P_I1, float* P_I2, int r, int c, int z)
{
	//create a normalized image histogram
	float Hst[256];
	for(int i=0; i<256; i++)
		Hst[i] = 0.0;
	
	for(int i=0; i<r*c*z; i++)
	{
		int v = (int) im[i];
		Hst[v]++;
	}

	for(int i=0; i<256; i++)
		Hst[i] /= (r*c*z);


	//The three-level min error thresholding algorithm
	float P0, U0, P1, U1, P2, U2, U, J, min_J;
	min_J = 1000000.0;
	// Try this: we need to define a penalty term that depends on the number of parameters
	//The penalty term is given as 0.5*k*ln(n)
	//where k is the number of parameters of the model and n is the number of samples
	//In this case, k=6 and n=256
	double PenTerm3 =  sqrt(6.0)*log(256.0);
	for(int i=0; i<254; i++)//to set the first threshold
	{
		//compute the current parameters of the first component
		P0 = U0 = 0.0;		
		for(int l=0; l<=i; l++)
		{
			P0+=Hst[l];
			U0+=(l+1)*Hst[l];
		}
		U0 /= P0;

		for(int j=i+1; j<255; j++)//to set the second threshold
		{
			//compute the current parameters of the second component
			P1 = U1 = 0.0;		
		    for(int l=i+1; l<=j; l++)
			{
				P1+=Hst[l];
				U1+=(l+1)*Hst[l];
			}
			U1 /= P1;

			//compute the current parameters of the third component
			P2 = U2 = 0.0;		
		    for(int l=j+1; l<=255; l++)
			{
				P2+=Hst[l];
				U2+=(l+1)*Hst[l];
			}
			U2 /= P2;
			
			//compute the overall mean
			U = P0*U0 + P1*U1 + P2*U2;

			//Compute the current value of the error criterion function
			J =  U - (P0*(log(P0)+U0*log(U0))+ P1*(log(P1)+U1*log(U1)) + P2*(log(P2)+U2*log(U2)));
			//Add the penalty term
			J +=PenTerm3;
			
			if(J<min_J)
			{
				min_J = J;
				Alpha1[0] = U0;
				P_I1[0] = P0;
				Alpha2[0] = U1;
				P_I2[0] = P1;
				Alpha3[0] = U2;				
			}
		}
	}

	//try this: see if using two components is better
	//The penalty term is given as sqrt(k)*ln(n)	
	//In this case, k=4 and n=256
	double PenTerm2 =  2*log(256.0);
	for(int i=0; i<254; i++)//to set the first threshold
	{
		//compute the current parameters of the first component
		P0 = U0 = 0.0;		
		for(int l=0; l<=i; l++)
		{
			P0+=Hst[l];
			U0+=(l+1)*Hst[l];
		}
		U0 /= P0;

		for(int j=i+1; j<255; j++)//to set the second threshold
		{
			//compute the current parameters of the second component
			P1 = U1 = 0.0;		
		    for(int l=j; l<=255; l++)
			{
				P1+=Hst[l];
				U1+=(l+1)*Hst[l];
			}
			U1 /= P1;

			//compute the overall mean
			U = P0*U0 + P1*U1;

			//Compute the current value of the error criterion function
			J =  U - (P0*(log(P0)+U0*log(U0))+ P1*(log(P1)+U1*log(U1)));
			//Add the penalty term
			J +=PenTerm2;
			if(J<min_J)
			{
				min_J = J;
				Alpha1[0] = U0;
				P_I1[0] = P0;
				Alpha2[0] = U1;
				P_I2[0] = P1;
				Alpha3[0] = -1; //Just a negative number to let the program knows that two levels will be used		
			}
		}
	}
	
}

void Seg_GC_Full_2D_Three_Level(unsigned char* IM,
                    int r, 
                    int c, 
                    double alpha_F, 
                    double alpha_B1, 
					double alpha_B2,
                    double P_I1, 
					double P_I2,
                    int* num_nodes, 
                    int* num_edges, 
                    unsigned short* Seg_out)
{    
    int curr_node;
    int rght_node;
    int down_node;
    int diag_node;
    double Df;
    double Db;
    double Dr;
    double Dd; 
    double Dg; 
    double sig;
    double F_H[256];
    double B_H[256];
    typedef Graph_B<int,int,int> GraphType;
       
    //Set the number of edges and the number of nodes and open the files that
    //will be used to save the weights
    num_nodes[0] = r*c;
    num_edges[0] = 3*r*c-2*r-2*c+1;
        	
    
    //Before entering the loop, compute the poisson probs 
    for(int i=0; i<256; i++)
    {
		if(i>=alpha_F)
			F_H[i] = (1-P_I1-P_I2)*compute_poisson_prob((int)alpha_F,alpha_F);
		else
			F_H[i] = (1-P_I1-P_I2)*compute_poisson_prob(i,alpha_F);
		if(i<=alpha_B1)
			B_H[i] = (P_I1)*compute_poisson_prob(int(alpha_B1),alpha_B1)+(P_I2)*compute_poisson_prob(int(alpha_B2),alpha_B2);
		else
        B_H[i] = (P_I1)*compute_poisson_prob(i,alpha_B1)+(P_I2)*compute_poisson_prob(i,alpha_B2);
    }
	
    //Here is the main loop.. 
    //For each point, compute the terminal and neighbor edge weights
	GraphType *g = new GraphType(/*estimated # of nodes*/ num_nodes[0], /*estimated # of edges*/ num_edges[0]); 
    for(int i=0; i<r; i++)
    {
        for(int j=0; j<c; j++)
        {			
			/*Get the terminal edges capacities and write them to a file
			These capacities represent penalties of assigning each point to
			either the fg or the bg. Now, I am using the distance of the point
			to the fg and bg. Then a short distance between a point and a class (fg or bg)
			means the penalty of the assignement is low and vice versa*/ 
								
		   //Isaac change on 5/9/08
		   //int intst = (int) IM[i][j];
			int intst = (int) IM[i*c + j];

            //Added by Yousef on jan 17, 2008
            //check if this is a seed point
            if(intst == 255)
            {
                Df = 0;
                Db = 1000;
            }
            else if(intst == 0)
            {
                Df = 1000;
                Db = 0;
            }
            else
            {                
                Df = -log(F_H[intst]);  //it was multiplied by .5          
                if(Df>1000.0)
                    Df = 1000;
                Db = -log(B_H[intst]);
                if(Db>1000.0)
                    Db=1000;
                
            }            
            
            curr_node = (i*c)+j; 
            
			g -> add_node();
			g -> add_tweights( curr_node,   /* capacities */ Df,Db);                            			
		}
	}

    sig = 25.0;
    for(int i=0; i<r-1; i++)
      {
      for(int j=0; j<c-1; j++)
        {			
        // get the neighbor edges capacities and write them to a file.
	// Now, each edge capacity between two neighbors p and q represent
	// the penalty for discontinuety. Since I am using the difference in
	// the intensities, I should take the inverse so that very similar
	// objects should have a large discontinuety penalty between them*/

	curr_node = (i*c)+j; 
        rght_node = curr_node+1;
        down_node = curr_node+c;
        diag_node = curr_node+c+1;
            
            
        //from Boykov's paper instead
	Dr = 10*exp(-pow((double)IM[i*c + j]-(double)IM[i*c + j+1],2)/(2*pow(sig,2)));
	g -> add_edge( curr_node, rght_node,    /* capacities */  Dr, Dr );	
            
                   
        //from Boykov's paper instead
Dd = 20*exp(-pow((double)IM[i*c + j]-(double)IM[(i+1)*c + j],2)/(2*pow(sig,2)));
g->add_edge( curr_node, down_node,    /* capacities */  Dd, Dd );
                 
Dg = 20*exp(-pow((double)IM[i*c + j]-(double)IM[(i+1)*c + (j+1)],2)/(2*pow(sig,2)));
g->add_edge( curr_node, diag_node,    /* capacities */  Dg, Dg );            
               
        }
      }    
	
    //Compute the maximum flow:
	g->maxflow();		//Alex DO NOT REMOVE


	int RR,CC;
	for(int i=0; i<num_nodes[0]; i++)
	{
		CC = ((long)i)%c;
        RR = (i-CC)/c;
		if(g->what_segment(i) == GraphType::SOURCE)
			Seg_out[RR*c + CC]=0;
		else
			Seg_out[RR*c + CC]=255;
	}

	delete g;
}