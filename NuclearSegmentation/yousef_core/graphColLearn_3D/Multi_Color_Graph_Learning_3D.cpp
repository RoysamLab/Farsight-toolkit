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

#include "Multi_Color_Graph_Learning_3D.h"
#include "gvc.h"

//added by yousef on 11/3/2008
#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
//#include "itkApproximateSignedDistanceMapImageFilter.h"
#include "itkImageFileWriter.h"
//#include "itkSignedDanielssonDistanceMapImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
//////////////////////////////

//A 3-D multivariate gaussian
double Multivar_Norm(double X, double Y, double Z, double Ux, double Uy, double Uz, double S00, double S01, double S11, double S02, double S12, double S22)
{
	double det_segma = S00*(S11*S22-S12*S12) - S01*(S01*S22-S12*S02) + S02*(S01*S12-S11*S02);
	double Sinv00 = (S22*S11-S12*S12)/det_segma;
	double Sinv01 = -(S22*S01-S12*S02)/det_segma;
	double Sinv02 = (S12*S01-S11*S02)/det_segma;
	double Sinv11 = (S22*S00-S02*S02)/det_segma;
	double Sinv12 = -(S12*S00-S01*S02)/det_segma;
	double Sinv22 = (S11*S00-S01*S01)/det_segma;
	X = X-Ux;
	Y = Y-Uy;
	Z = Z-Uz;
	double Mah_Dist = X*(X*Sinv00+Y*Sinv01+Z*Sinv02) + Y*(X*Sinv01+Y*Sinv11+Z*Sinv12) + Z*(X*Sinv02+Y*Sinv12+Z*Sinv22);
	double Ex = exp(-.5*Mah_Dist);
	double A = 1/(pow((2*3.1416),1.5)*sqrt(det_segma));

	return (A*Ex);

}


float* multiColGraphLearning(float* X_vals, unsigned short* labs_vals, unsigned short* color_im, size_t r, size_t c, size_t z, int *NC, int refinemetRange)
{  
	//inputs are:
	//1- The LoG response image
	//2- The labeled image (1 to numofobjs) 
	//3- Image dimensions
	//4- The number of reqiured colors (ideally >= 4)  

	int*** labs_im;
	int max_lab,ncolors;
	int** RAG;
	float* out;
	int* ColorOut;
	double** U;
	double*** Segma;
	double* P_I;
	//double** Z;
	int* Z_sum;    
	int L, L1, L2, L3, L4, L5, L6, L7;



	//create a 3-D image that will hold the edge points only    
	labs_im = (int ***) malloc(r*sizeof(int**)); 
	max_lab = 0;
	for(int i=0; i<r-1; i++)
	{        
		labs_im[i] = (int **) malloc(c*sizeof(int*));
		for(int j=0; j<c-1; j++)
		{			
			labs_im[i][j] = (int *) malloc(z*sizeof(int));
			for(int k=0; k<z-1; k++)
			{
				L = /*(int)*/labs_vals[(k*r*c)+(j*r)+i]; 
				if( L == 0)
				{
					labs_im[i][j][k] = 0;
					continue;
				}	
				L1 = /*(int)*/labs_vals[(k*r*c)+(j*r)+(i+1)]; 
				L2 = /*(int)*/labs_vals[(k*r*c)+((j+1)*r)+i];
				L3 = /*(int)*/labs_vals[(k*r*c)+((j+1)*r)+(i+1)];
				L4 = /*(int)*/labs_vals[((k+1)*r*c)+(j*r)+i];
				L5 = /*(int)*/labs_vals[((k+1)*r*c)+((j+1)*r)+i];
				L6 = /*(int)*/labs_vals[((k+1)*r*c)+(j*r)+(i+1)];
				L7 = /*(int)*/labs_vals[((k+1)*r*c)+((j+1)*r)+(i+1)];

				if((L!=L1 && L1!=0) || (L!=L2 && L2!=0) ||(L!=L3 && L3!=0) || (L!=L4 && L4!=0) || (L!=L5 && L5!=0) || (L!=L6 && L6!=0) || (L!=L7 && L7!=0))
				{					
					labs_im[i][j][k] = L;
					if(labs_im[i][j][k] > max_lab)
						max_lab = /*(int)*/labs_im[i][j][k];
				}
				else
					labs_im[i][j][k] = 0;
			}

		}
	}

	labs_im[r-1] = (int **) malloc(c*sizeof(int*));
	for(int j=0; j<c; j++)
	{        
		labs_im[r-1][j] = (int *) malloc(z*sizeof(int));
		for(int k=0; k<z; k++)
		{		
			labs_im[r-1][j][k] = 0;
			if(labs_vals[(k*r*c)+(j*r)+(r-1)]>max_lab)
				max_lab=/*(int)*/labs_vals[(k*r*c)+(j*r)+(r-1)];
		}
	}

	for(int i=0; i<r; i++)
	{
		labs_im[i][c-1] = (int *) malloc(z*sizeof(int));
		for(int k=0; k<z; k++)
		{			
			labs_im[i][c-1][k] = 0;			
			if(labs_vals[(k*r*c)+((c-1)*r)+i]>max_lab)
				max_lab=/*(int)*/labs_vals[(k*r*c)+((c-1)*r)+i];
		}
	}


	for(int i=0; i<r; i++)
	{
		for(int j=0; j<c; j++)
		{
			labs_im[i][j][z-1] = 0;				
			if(labs_vals[((z-1)*r*c)+(j*r)+i]>max_lab)
				max_lab=/*(int)*/labs_vals[((z-1)*r*c)+(j*r)+i];
		}
	}


	//Build the region adjacency graph    
	std::cout<<"Building region adjacency graph...";
	RAG = (int **) malloc(max_lab*sizeof(int*));
	for(int i=0; i<max_lab; i++)
	{        
		RAG[i] = (int *) malloc(max_lab*sizeof(int));
		for(int j=0; j<max_lab; j++)
			RAG[i][j] = 0;
	}


	for(int i=0; i<r-1; i++)
	{        
		for(int j=0; j<c-1; j++)
		{	
			for(int k=0; k<z-1; k++)
			{
				L = labs_im[i][j][k];
				if( L == 0)
					continue;
				else
				{			
					L1 = /*(int)*/labs_vals[(k*r*c)+(j*r)+(i+1)]; 
					L2 = /*(int)*/labs_vals[(k*r*c)+((j+1)*r)+i];
					L3 = /*(int)*/labs_vals[(k*r*c)+((j+1)*r)+(i+1)];
					L4 = /*(int)*/labs_vals[((k+1)*r*c)+(j*r)+i];
					L5 = /*(int)*/labs_vals[((k+1)*r*c)+((j+1)*r)+i];
					L6 = /*(int)*/labs_vals[((k+1)*r*c)+(j*r)+(i+1)];
					L7 = /*(int)*/labs_vals[((k+1)*r*c)+((j+1)*r)+(i+1)];

					if(L!=L1 && L1!=0 && L1<=max_lab)
						RAG[L-1][L1-1] = RAG[L1-1][L-1] = 1;
					if(L!=L2 && L2!=0 && L2<=max_lab)
						RAG[L-1][L2-1] = RAG[L2-1][L-1] = 1;
					if(L!=L3 && L3!=0 && L3<=max_lab)
						RAG[L-1][L3-1] = RAG[L3-1][L-1] = 1;
					if(L!=L4 && L4!=0 && L4<=max_lab)
						RAG[L-1][L4-1] = RAG[L4-1][L-1] = 1;
					if(L!=L5 && L5!=0 && L5<=max_lab)
						RAG[L-1][L5-1] = RAG[L5-1][L-1] = 1;
					if(L!=L6 && L6!=0 && L6<=max_lab)
						RAG[L-1][L6-1] = RAG[L6-1][L-1] = 1;
					if(L!=L7 && L7!=0 && L7<=max_lab)
						RAG[L-1][L7-1] = RAG[L7-1][L-1] = 1;
				}
			}                
			//free(labs_im[i][j]);         
		}
		//free(labs_im[i]);        
	}    
	//free(labs_im);
	std::cout<<"done"<<std::endl;

	//copy the RAG into an std vector of vectors	
	std::vector<std::vector<int> > MAP;
	MAP.resize(max_lab);
	std::vector<std::vector<int> > MAP2;
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
		for(unsigned int j=0; j<MAP2[i].size() ; j++)
		{
			NG1 = MAP2[i][j];
			//for each jth neighbor of the ith cell
			for(unsigned int k=0; k<MAP2[NG1-1].size() ; k++)
			{
				//add each kth neighbor of the jth cell to the ith cell
				NG2 = MAP2[NG1-1][k];
				if(NG2!=(i+1))
					MAP[i].push_back(NG2);
			}				            
		}
	}

	//start the graph coloring using Sumit's sequential coloring code
	std::cout<<"Starting graph coloring...";
	GVC* Gcol = new GVC();
	ncolors = NC[0];	
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
	std::cout<<"done with "<<ncolors<<" colors"<<std::endl;

	//Now let's do the learning step
	//here, each pixel will have a probability of belonging to each one of the 
	//colors (4 colors ideally)
	//To do that, we assume that each object (cell) can be represented by a gaussian model
	Segma = (double ***) malloc(max_lab*sizeof(double**));    
	U = (double **) malloc(max_lab*sizeof(double*));    
	P_I = (double *) malloc(max_lab*sizeof(double));
	Z_sum = (int *) malloc(max_lab*sizeof(int));    
	for(int i=0; i<max_lab; i++)
	{
		Z_sum[i] = 0;
		P_I[i] = 0; 
		U[i] = (double *) malloc(3*sizeof(double));        
		Segma[i] = (double **) malloc(3*sizeof(double*));
		for(int j=0; j<3; j++)
		{
			U[i][j] = 0;            
			Segma[i][j] = (double *) malloc(3*sizeof(double));
			for(int k=0; k<3; k++)
				Segma[i][j][k] = 0;
		}
	}

	//Now, compute the estimates of the parameters
	//Start with the means
	float val;
	int label;
	double ZS = 0;
	for(int i=0; i<r; i++)
	{
		for(int j=0; j<c; j++)
		{
			for(int k=0; k<z; k++)
			{
				label = /*(int)*/labs_vals[(k*r*c)+(j*r)+i]; 
				if(label == 0 || label>max_lab)
					continue;
				val = X_vals[(k*r*c)+(j*r)+i]; //Note,, for this function to work, LPG_im should be normalized between 0 and the maximum repetition		
				//try this:
				//val = 1;
				//////////
				U[label-1][0] += (j*val);
				U[label-1][1] += (i*val);
				U[label-1][2] += (k*val);
				Z_sum[label-1]+= (int)val;
				ZS+=val;
			}
		}
	}

	for(int i=0; i<max_lab; i++)
	{                    
		U[i][0] = (U[i][0] / (double)Z_sum[i]);
		U[i][1] = (U[i][1] / (double)Z_sum[i]);
		U[i][2] = (U[i][2] / (double)Z_sum[i]);
		P_I[i] = (double)Z_sum[i]/ZS;
	}

	//Then the convariance matrices
	for(int i=0; i<r; i++)
	{
		for(int j=0; j<c; j++)
		{
			for(int k=0; k<z; k++)
			{
				label = /*(int)*/labs_vals[(k*r*c)+(j*r)+i];            
				if(label == 0 || label>max_lab)
					continue;
				val = X_vals[(k*r*c)+(i*c)+j];
				//try this:
				//val = 1;
				//////////
				Segma[label-1][0][0] += ((j-U[label-1][0])*(j-U[label-1][0])*val);
				Segma[label-1][0][1] += ((j-U[label-1][0])*(i-U[label-1][1])*val);
				Segma[label-1][0][2] += ((j-U[label-1][0])*(k-U[label-1][2])*val);
				Segma[label-1][1][0] += ((i-U[label-1][1])*(j-U[label-1][0])*val);
				Segma[label-1][1][1] += ((i-U[label-1][1])*(i-U[label-1][1])*val);
				Segma[label-1][1][2] += ((i-U[label-1][1])*(k-U[label-1][2])*val);
				Segma[label-1][2][0] += ((k-U[label-1][2])*(j-U[label-1][0])*val);
				Segma[label-1][2][1] += ((k-U[label-1][2])*(i-U[label-1][1])*val);
				Segma[label-1][2][2] += ((k-U[label-1][2])*(k-U[label-1][2])*val);
			}
		}
	}

	for(int i=0; i<max_lab; i++)
	{
		Segma[i][0][0] /= (double)Z_sum[i];
		Segma[i][0][1] /= (double)Z_sum[i];
		Segma[i][0][2] /= (double)Z_sum[i];
		Segma[i][1][0] = Segma[i][0][1];
		Segma[i][1][1] /= (double)Z_sum[i];
		Segma[i][1][2] /= (double)Z_sum[i];
		Segma[i][2][0] = Segma[i][0][2];
		Segma[i][2][1] /= Segma[i][1][2];
		Segma[i][2][2] /= (double)Z_sum[i];
	}



	//Added by yousef on 11/3/2008
	//Compute the distance map from the edges image
	distToEdge(labs_im, r, c, z);
	///////////////////////////////////////////////////////



	//Now, for each pixel in the image, compute its probability to belong to each on of the classes (colors)
	//And set the data term as -ln(Pr)    
	out = (float *) malloc(r*c*z*(ncolors+1)*sizeof(float));
	if(out == NULL)
	{
		double mem_required = r*c*z*(ncolors+1)*sizeof(float) / (1024 * 1024 * 1024.0F);
		std::cerr<<"failed\ncannot allocate " << mem_required << " GB of memory to the data term matrix"<<std::endl;
		exit(0);
	}

	double *Pr = (double *) malloc((ncolors+1)*sizeof(double));
	int C, cl, V;
	double P;
	//int intst;
	for(int i=0; i<r; i++)
	{        
		for(int j=0; j<c; j++)
		{			     
			for(int h=0; h<z; h++)
			{				
				val = labs_vals[(h*r*c)+(j*r)+i];				          
				for(cl=0; cl<ncolors+1;cl++)
					Pr[cl] =0;
				if(val == 0 || val>max_lab)
				{
					Pr[0] = 1;
					//Added By Yousef on 10-27-2008
					//make sure no label is larger than max_lab
					color_im[(h*r*c)+(i*c)+j] = 0 ;
				}
				else
				{
					val = val-1;
					//check the prob for this point to be in its current object
					//and all the adjacent objects
					C = ColorOut[(int) val];				
					//Modified by Yousef on 11/3/2008
					//If the distance of the point is larger than the refinement range
					//then, just set the probability to to belong to the current object to one
					if(labs_im[i][j][h]>refinemetRange)
						Pr[C+1] = 1;
					else
						Pr[C+1] = Multivar_Norm(j,i,h,U[(int)val][0],U[(int)val][1],U[(int)val][2],Segma[(int)val][0][0],Segma[(int)val][1][0], Segma[(int)val][1][1],Segma[(int)val][0][2],Segma[(int)val][1][2],Segma[(int)val][2][2]);
					//Added By Yousef on 10-27-2008
					//Relable the initially segmented image using the colors
					color_im[(h*r*c)+(i*c)+j] = C+1;
					for (unsigned int k = 0; k < MAP[(int)val].size(); k++)
					{
						V = MAP[(int)val][k]-1;
						C = ColorOut[V];
						//Modified by Yousef on 11/3/2008
						//If the distance of the point is larger than the refinement range
						//then, just set the probability to belong to other objects to zero
						//if(labs_im[i][j][h]>refinemetRange)
						//	P = 0;
						//else
						P = Multivar_Norm(j,i,h,U[V][0],U[V][1],U[V][2],Segma[V][0][0],Segma[V][1][0], Segma[V][1][1],Segma[V][0][2],Segma[V][1][2],Segma[V][2][2]);
						Pr[C+1] = max(Pr[C+1],P);
					}               
				}			
				for(int cc=0; cc<ncolors+1; cc++)					
					/*out[(cc*r*c*z)+(h*r*c)+(j*r)+i]*/out[(j+i*c+h*c*r)*(ncolors+1) + cc]=.25*min(-log(Pr[cc]),100.0); 
			}
		}
	}

	//added by yousef on 11/3/2008
	for(int i=0; i<r; i++)
	{
		for(int j=0; j<c; j++)		
			free(labs_im[i][j]);                 
		free(labs_im[i]);           
	}
	free(labs_im);

	for(int i=0; i<max_lab; i++)
	{        				
		for(int j=0; j<3; j++)
		{            
			free(Segma[i][j]);
		}
		Segma[i] = (double **) malloc(3*sizeof(double*));
		free(U[i]);
	}
	free(Z_sum);
	free(P_I);

	return out;
}

//added by Yousef on 11/3/2008
void distToEdge(int *** edge_im, int r, int c, int z)
{
	//  Types should be selected on the desired input and output pixel types.
	typedef    double     InputPixelType;
	typedef    double     OutputPixelType;


	//  The input and output image types are instantiated using the pixel types.
	typedef itk::Image< /*unsigned char, 3 */InputPixelType,  3 >   InputImageType;
	typedef itk::Image< OutputPixelType, 3 >   OutputImageType;

	//Create an itk image
	InputImageType::Pointer im;
	im = InputImageType::New();
	InputImageType::PointType origin;
	origin[0] = 0.; 
	origin[1] = 0.;    
	origin[2] = 0.;    
	im->SetOrigin( origin );

	InputImageType::IndexType start;
	start[0] =   0;  // first index on X
	start[1] =   0;  // first index on Y    
	start[2] =   0;  // first index on Z    
	InputImageType::SizeType  size;
	size[0]  = c;  // size along X
	size[1]  = r;  // size along Y
	size[2]  = z;  // size along Z

	InputImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );

	double spacing[3];
	spacing[0] = 1; //spacing along x
	spacing[1] = 1; //spacing along y
	spacing[2] = 1; //spacing along z

	im->SetRegions( region );
	im->SetSpacing(spacing);
	im->Allocate();
	im->FillBuffer(0);
	im->Update();	
	//copy the input image into the ITK image
	typedef itk::ImageRegionIteratorWithIndex< InputImageType > IteratorType;
	IteratorType iterator1(im,im->GetRequestedRegion());
	for(int k=0; k<z; k++)	
	{				
		for(int i=0; i<r; i++)		
		{
			for(int j=0; j<c; j++)
			{				
				if(edge_im[i][j][k]>0)
					iterator1.Set(255.0);//IM[i]);
				else
					iterator1.Set(0.0);
				++iterator1;
			}
		}
	}


	/*typedef itk::ImageFileWriter< InputImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	std::string hh = "lookatme.tif";
	writer->SetFileName(hh.c_str());
	writer->SetInput( im );
	writer->Update();*/


	/*typedef itk::ApproximateSignedDistanceMapImageFilter<InputImageType, OutputImageType > DTFilter ;
	DTFilter::Pointer dt_obj= DTFilter::New() ;
	dt_obj->SetInput(im) ;
	dt_obj->SetInsideValue(0.0);
	dt_obj->SetOutsideValue(255.0);*/

	//typedef itk::SignedDanielssonDistanceMapImageFilter<InputImageType, OutputImageType > DTFilter ;
	//DTFilter::Pointer dt_obj= DTFilter::New() ;
	typedef itk::SignedMaurerDistanceMapImageFilter<InputImageType, OutputImageType>  DTFilter;
	DTFilter::Pointer dt_obj= DTFilter::New() ;
	dt_obj->SetInput(im) ;	
	try{
		dt_obj->Update() ;
	}
	catch( itk::ExceptionObject & err ){
		std::cerr << "Error calculating distance transform: " << err << endl ;
		//exit(0);
	}

	//   Copy the resulting image into the input array
	typedef itk::ImageRegionIteratorWithIndex< InputImageType > IteratorType;
	IteratorType iterate(dt_obj->GetOutput(),dt_obj->GetOutput()->GetRequestedRegion());

	for(int k=0; k<z; k++)	
	{				
		for(int i=0; i<r; i++)		
		{
			for(int j=0; j<c; j++)
			{				
				edge_im[i][j][k] = fabs(iterate.Get());	       	
				++iterate;
			}
		}
	}	
}
