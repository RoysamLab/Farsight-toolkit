#include "alpha_expansion_3d.h"
#include "GraphCut3dConstr.cpp"
#include "GraphCutMex.cpp"

void start_alpha_expansion(float* im, int* seg_im, float* Dterms, int R, int C, int Z, int K)
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

//void graph_cut_close(GraphHandle* gh, char mode);
//{
//}

