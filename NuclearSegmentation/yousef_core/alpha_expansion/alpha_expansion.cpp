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

#include "alpha_expansion.h"
#include "GraphCutConstr.cpp"
#include "GraphCutMex.cpp"

void start_alpha_expansion(float* im, unsigned short* seg_im, float* Dterms, int R, int C, int Z, int K)
{
	float w = 10.0;
	//initialize the smoothness constat part
	float* SC = new float[K*K];
	float val;
	for(int i=0; i<K; i++)
	{
		for(int j=0; j<K; j++)
		{
			if(i==j)
				val = 0.0;
			else
				val = 1.0;

			SC[(j*K)+i] = w*val;
		}
	}

	//start by calling the function that builds the graph (OPEN)
	GCoptimization *MyGraph = GraphCut3dConstr(im, Dterms, SC, C, R, Z, K);

	//pass the initial labels
	//GraphCutMex(MyGraph, 's', seg_im, 1);

	//then, call the function the applys the alpha expansion (EXPAND)
	int iter = 1; //if zero, then epand till convergence, else it gives the number of iterations
	GraphCutMex(MyGraph, 'e', seg_im, iter); //the name of this function will be changed to something without MEX
	
	//finally, call the function the destroys (deallocate) the graph (CLOSE)
	GraphCutMex(MyGraph, 'c',seg_im, iter); //the name of this function will be changed to something without MEX
	free(SC);
}
 
//These functions can be used if you want more separation between the main alpha expansion function and the other functions
//void graph_cut_open(float* im, float* Dterm, float* SC, int R, int C, int Z, int K, GraphHandle* gh)
//{
//}

//void graph_cut_expand(GraphHandle* gh, char mode, int* seg_im)
//{
//}

void alpha_expansion_2d( float *im, float *sublogImg, unsigned short *subclustImg, int R, int C )
{
	float *hCue, *vCue;
	int K;	
	//starting graph-roloring based learning
	K = 10000;
	float* Dterms  = multiColGraphLearning(sublogImg, subclustImg, R, C, &K);
	K = K+1;
	std::cerr<<"    Graph Coloring done with "<<K-1<<" colors"<<std::endl;
	
	std::cerr<<"    Starting alpha-expansion..";
	float w = 10.0;
	//initialize the smoothness constat part
	float* SC = new float[K*K];
	float val;
	for(int i=0; i<K; i++)
	{
		for(int j=0; j<K; j++)
		{
			if(i==j)
				val = 0.0;
			else
				val = 1.0;

			SC[(j*K)+i] = w*val;
		}
	}
	float sig = 30.0;
	hCue = new float[R*C];
	vCue = new float[R*C];
	for(int i=0; i<R-1; i++)
	{
		for(int j=0; j<C-1; j++)
		{
			hCue[(j*R)+i] = exp(-fabs(im[(i*C)+j] - im[(i*C)+j+1])/sig);
			vCue[(j*R)+i] = exp(-fabs(im[(i*C)+j] - im[((i+1)*C)+j])/sig);
		}
	}	
	int i=R-1;
	for(int j=0; j<C-1; j++)
	{
		hCue[(j*R)+i] = exp(-fabs(im[(i*C)+j] - im[(i*C)+j+1])/sig);
		vCue[(j*R)+i] = vCue[(j*R)+(i-1)];
	}
	int j=C-1;
	for(int i=0; i<R-1; i++)
	{
		vCue[(j*R)+i] = exp(-fabs(im[(i*C)+j] - im[((i+1)*C)+j])/sig);
		hCue[(j*R)+i] = hCue[((j-1)*R)+i];
	}	
	vCue[((C-1)*R)+(R-1)] = vCue[((R-1)*C)+(C-2)];
	hCue[((C-1)*R)+(R-1)] = hCue[((R-1)*C)+(C-2)];


	//start by calling the function the builds the graph (OPEN)
	GCoptimization *MyGraph = GraphCutConstr(Dterms, SC, hCue, vCue, R, C, K);

	//then, call the function the applys the alpha expansion (EXPAND)	
	int iter = 0; //if zero, then epand till convergence, else it gives the number of iterations
	GraphCutMex(MyGraph, 'e', subclustImg, iter); //the name of this function will be changed to something without MEX
	//finally, call the function the destroys (deallocate) the graph (CLOSE)	
	GraphCutMex(MyGraph, 'c',subclustImg, iter); //the name of this function will be changed to something without MEX
	
	delete[] hCue;
	delete[] vCue;
}
