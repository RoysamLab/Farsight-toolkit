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

#include "cell_binarization.h"
//#include "maxflow.cpp"
//#include "graph.cpp"

using namespace std;

//Main function for 2-D binarization
int Cell_Binarization_2D(unsigned char* imgIn, int *imgOut, int R, int C, int shd)
{			
	//Now, to do the binarization, follow these steps:
	//1- Assuming that the histogram of the image is modeled by a mixture of two 
	//poisson distributions, estimate the parameters of the mixture model 
	float alpha_B, alpha_F, P_I;
	alpha_B = alpha_F = P_I = 0;
	//CompMixPoss(imgIn, &alpha_B, &alpha_F, &P_I, R, C, shd); 	
	MinErrorThresholding(imgIn, &alpha_B, &alpha_F, &P_I, R, C, 1, shd,imgOut); 	

	//2- Use graph cuts to binarize the image
	//this function will start by graph building (learning step) and then it will 
	//do the inference step (max-flow)
	int n_nodes, n_edges;
	Seg_GC_Full_2D(imgIn, R, C, alpha_F, alpha_B, P_I, &n_nodes, &n_edges, imgOut);
	//cout<<"Finalizing Binarization..";
	//post_binarization(imgOut, 2, 2, R, C, 1);
	//cout<<"done"<<endl;
	return 1;	
}

//Main function for 3-D binarization
int Cell_Binarization_3D(unsigned char *imgIn, int* imgOut, int R, int C, int Z, int shd) //modifed by Yousef on 5-20-2008.. The first input change from uchar* to int*
{			
	//Now, to do the binarization, follow these steps:
	//1- Assuming that the histogram of the image is modeled by a mixture of two 
	//poisson distributions, estimate the parameters of the mixture model 
	float alpha_B, alpha_F, P_I;
	alpha_B = alpha_F = P_I = 0;	
	//CompMixPoss3D(imgIn, &alpha_B, &alpha_F, &P_I, R, C, Z, shd); 	
	MinErrorThresholding(imgIn, &alpha_B, &alpha_F, &P_I, R, C, Z, shd,imgOut); 	
	
	//2- Use graph cuts to binarize the image
	//this function will start by graph building (learning step) and then it will 
	//do the inference step (max-flow)	
	//Seg_GC_Full_3D(imgIn, R, C, Z, alpha_F, alpha_B, P_I, imgOut);
	//cout<<"Finalizing Binarization..";
	//post_binarization(imgOut, 2, 3, R, C, Z);
	//cout<<"done"<<endl;

	//Added by Yousef on 11-18-2008: To save memory, divide the image into Blocks and 
	//apply GC on each independently	
	int *subImgBlock = new int[6];//[x1,y1,z1,x2,y2,z1]
	int block_divisor = 4;
	subImgBlock[4] = 0;
	subImgBlock[5] = Z;
	for(int i=0; i<R; i+=R/block_divisor)
	{
		for(int j=0; j<C; j+=C/block_divisor)
		{
			subImgBlock[0] = j;
			subImgBlock[1] = (int)j+C/block_divisor;
			subImgBlock[2] = i;
			subImgBlock[3] = (int)i+R/block_divisor;
			if(subImgBlock[1] > C)
				subImgBlock[1] = C;
			if(subImgBlock[3] > R)
				subImgBlock[3] = R;

			Seg_GC_Full_3D_Blocks(imgIn, R, C, Z, alpha_F, alpha_B, P_I, imgOut,subImgBlock);
		}
	}
	//Added By Yousef: May 20, 2008
	//3-Convert the binary image into a connected component image
	//int num_objects = getConnCompImage(imgOut, 26, 25, R, C, Z);
	//num_objects = getConnCompImage(imgOut, 26, 25, R, C, Z);
		
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
                    int* Seg_out)
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
Dd = 10*exp(-pow((double)IM[i*c + j]-(double)IM[(i+1)*c + j],2)/(2*pow(sig,2)));
g->add_edge( curr_node, down_node,    /* capacities */  Dd, Dd );
                 
Dg = 10*exp(-pow((double)IM[i*c + j]-(double)IM[(i+1)*c + (j+1)],2)/(2*pow(sig,2)));
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
	std::cerr << "Poisson Probabilities Computed" << std::endl;
    
	//Construct the graph
	GraphType *g = new GraphType(/*estimated # of nodes*/ num_nodes, /*estimated # of edges*/ num_edges); 

	std::cerr << "Graph Memory Allocated" << std::endl;

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

	std::cout << "First Loop Complete" << std::endl;
	
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
	
	std::cout << "Second Loop Complete" << std::endl;
	
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

void MinErrorThresholding(unsigned char* img, float* alpha_B, float* alpha_A, float* P_I, int R, int C, int Z, int shiftDown, int *imgOut)
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


void Seg_GC_Full_3D_Blocks(unsigned char* IM, int r, int c, int z, double alpha_F, double alpha_B, double P_I, int* Seg_out, int* imBlock)
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
	std::cerr << "Poisson Probabilities Computed" << std::endl;
    
	//Construct the graph
	GraphType *g = new GraphType(/*estimated # of nodes*/ num_nodes, /*estimated # of edges*/ num_edges); 

	std::cerr << "Graph Memory Allocated" << std::endl;

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

	std::cout << "First Loop Complete" << std::endl;
	
	sig = 30.0;
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
	
	std::cout << "Second Loop Complete" << std::endl;
	
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

/*int post_binarization(int* bin_im, int kernel_size, int dim, int r, int c, int z)
{	  	

	int const dim2 = 3;
	typedef unsigned char                PixelType;
	//if(dim == 2)
	//	typedef itk::Image< PixelType, 2 > ImageType;
	//else
		typedef itk::Image< PixelType, dim2 > ImageType;

	//Create an itk image
	ImageType::Pointer im;
	im = ImageType::New();
	ImageType::PointType origin;
    origin[0] = 0.; 
    origin[1] = 0.;    
	if(dim == 3)
		origin[2] = 0.;    
    im->SetOrigin( origin );

    ImageType::IndexType start;
    start[0] =   0;  // first index on X
    start[1] =   0;  // first index on Y 
	if(dim == 3)
		start[2] =   0;  // first index on Z 
    ImageType::SizeType  size;
    size[0]  = c;  // size along X
    size[1]  = r;  // size along Y
	if(dim == 3)
		size[2]  = z;  // size along Z

    ImageType::RegionType region;
    region.SetSize( size );
    region.SetIndex( start );
    
    im->SetRegions( region );
    im->Allocate();
    im->FillBuffer(0);
	im->Update();
	
	//copy the input image into the ITK image
	typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
	IteratorType iterator1(im,im->GetRequestedRegion());
	int lng = r*c*z;
	for(int i=0; i<lng; i++)
	{		
		iterator1.Set((unsigned char)bin_im[i]);
		++iterator1;	
	}
	
	//create the structuring elemet
	typedef itk::BinaryBallStructuringElement< PixelType, dim2 > KernelType;
	KernelType ball;
	KernelType ball2;
	KernelType::SizeType ballSize;
	KernelType::SizeType ballSize2;
	ballSize.Fill( kernel_size );
	ballSize2.Fill( kernel_size-1 );
	ball.SetRadius( ballSize );
	ball2.SetRadius( ballSize2 );
	ball.CreateStructuringElement();
	ball2.CreateStructuringElement();

	//start by opening
	typedef itk::BinaryMorphologicalOpeningImageFilter< ImageType, ImageType, KernelType > FilterType;
	FilterType::Pointer filter1 = FilterType::New();
	filter1->SetInput( im );
	filter1->SetKernel( ball );
	if( filter1->GetBackgroundValue() != 0 )
	{
		std::cerr << "Wrong Background default value" << std::endl;
		return EXIT_FAILURE;
    }
	filter1->SetBackgroundValue( 0 );

	if( filter1->GetForegroundValue() != 255 )
    {
		std::cerr << "Wrong Foreground default value" << std::endl;
		return EXIT_FAILURE;
	}
	filter1->SetForegroundValue( 255 ) ;
	
	//try
	//{
	//	filter1->Update() ;
	//}
	//catch( itk::ExceptionObject & err )
	//{
	//	std::cerr << "Error in performing binary morphological opening: " << err << endl ;
	//	return -1;
	//}

	//then closing
	typedef itk::BinaryMorphologicalClosingImageFilter< ImageType, ImageType, KernelType > FilterType2;
	FilterType2::Pointer filter2 = FilterType2::New();
	filter2->SetInput( filter1->GetOutput() );
	filter2->SetKernel( ball2 );
	 // test the default attribute values, and exercise the accesors
	if( !filter2->GetSafeBorder() )
	{
		std::cerr << "Wrong SafeBorder default value" << std::endl;
		return EXIT_FAILURE;
    }
	filter2->SafeBorderOff();
	filter2->SafeBorderOn();
	filter2->SetSafeBorder( 0 );
	if( filter2->GetForegroundValue() != 255 )
    {
		std::cerr << "Wrong Foreground default value" << std::endl;
		return EXIT_FAILURE;
    }
	filter2->SetForegroundValue( 255 );
	try
	{
		filter2->Update() ;
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "Error in performing binary morphological opening: " << err << endl ;
		return -1;
	}

	IteratorType iterate2(filter2->GetOutput(),filter2->GetOutput()->GetRequestedRegion());	
	for(int i=0; i<lng; i++)
	{		
		int ii = bin_im[i];
		bin_im[i] = (int)iterate2.Get();
		++iterate2;			
	}
	return 1;
}*/

