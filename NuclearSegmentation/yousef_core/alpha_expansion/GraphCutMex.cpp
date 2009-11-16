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

////////////////////////////////////////////////////////////////////
//This file was modified by Yousef Al-Kofahi (RPI) on May 22nd 2008
//The function was intended to be used in Matlab as a Mex file 
//Details from the original Authors are below
////////////////////////////////////////////////////////////////////

//#include "mex.h"
#include "GCoptimization.h"
#include "GraphCut.h"
#include <stdlib.h>
#include <iostream>

/* Declarations */
//void SetLabels(GCoptimization *MyGraph, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void SetLables(GCoptimization *MyGraph, unsigned short* seg_im);
//void GetLabels(GCoptimization *MyGraph, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
//void ABswaps(GCoptimization *MyGraph, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
//void Expand(GCoptimization *MyGraph, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void Expand(GCoptimization *MyGraph, unsigned short* seg_im, int nrhs, int max_iterations);
//void Energy(GCoptimization *MyGraph, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
//GCoptimization* GetGCHandle(const mxArray *x);    /* extract ahndle from mxArry */


/*
 * Matlab wrapper for Weksler graph cut implementation
 * 
 * GCoptimization cde by Olga Veksler.
 * Wrapper code by Shai Bagon.
 *
 *
 *   Performing Graph Cut operations on a 2D grid.
 *   
 *   Usage:
 *       [gch ...] = GraphCutMex(gch, mode ...);
 *   
 *   Notes:
 *   1. Data types are crucial!
 *   2. The embedded implementation treat matrices as a row stack, while
 *       matlab treats them as column stack; thus there is a need to
 *       transpose the label matrices and the indices passed to expand
 *       algorithm.
 *
 *   Inputs:
 *   - gch: a valid Graph Cut handle (use GraphCutConstr to create a handle).
 *   - mode: a char specifying mode of operation. See details below.
 *
 *   Output:
 *   - gch: A handle to the constructed graph. Handle this handle with care
 *              and don't forget to close it in the end!
 *
 *   Possible modes:
 *
 *   - 's': Set labels
 *           [gch] = GraphCutMex(gch, 's', labels)
 *
 *       Inputs:
 *           - labels: a width*height array of type int32, containing a
 *              label per pixel. note that labels values must be is range
 *              [0..num_labels-1]. 
 *
 *   - 'g': Get current labeling
 *           [gch labels] = GraphCutMex(gch, 'g')
 *
 *       Outputs:
 *           - labels: a width*height array of type int32, containing a
 *              label per pixel. note that labels values must be is range
 *              [0..num_labels-1].
 *
 *   - 'n': Get current values of energy terms
 *           [gch se de] = GraphCutMex(gch, 'n')
 *
 *       Outputs:
 *           - se: Smoothness energy term.
 *           - de: Data energy term.
 *
 *   - 'e': Perform labels expansion
 *           [gch labels] = GraphCutMex(gch, 'e')
 *           [gch labels] = GraphCutMex(gch, 'e', iter)
 *           [gch labels] = GraphCutMex(gch, 'e', [], label)
 *           [gch labels] = GraphCutMex(gch, 'e', [], label, indices)
 *
 *       When no inputs are provided, GraphCut performs expansion steps
 *       until it converges.
 *
 *       Inputs:
 *           - iter: a double scalar, the maximum number of expand
 *                      iterations to perform.
 *           - label: int32 scalar denoting the label for which to perfom
 *                        expand step.
 *           - indices: int32 array of linear indices of pixels for which
 *                            expand step is computed. indices _MUST_ be zero offset and not
 *                            one offset like matlab usuall indexing system, i.e., to include
 *                            the first pixel in the expand indices array must contain 0 and
 *                            NOT 1!
 *
 *       Outputs:
 *           - labels: a width*height array of type int32, containing a
 *              label per pixel. note that labels values must be is range
 *              [0..num_labels-1].
 *
 *   - 'a': Perform alpha - beta swappings
 *           [gch labels] = GraphCutMex(gch, 'a')
 *           [gch labels] = GraphCutMex(gch, 'a', iter)
 *           [gch labels] = GraphCutMex(gch, 'a', label1, label2)
 *
 *       When no inputs are provided, GraphCut performs alpha - beta swaps steps
 *       until it converges.
 *
 *       Inputs:
 *           - iter: a double scalar, the maximum number of swap
 *                      iterations to perform.
 *           - label1, label2: int32 scalars denoting two labels for swap
 *                                       step.
 *
 *       Outputs:
 *           - labels: a width*height array of type int32, containing a
 *              label per pixel. note that labels values must be is range
 *              [0..num_labels-1].
 *
 *   - 'c': Close the graph and release allocated resources.
 *       [gch] = GraphCutMex(gch,'c');
 *
 *
 *   See Also:
 *       GraphCutConstr
 *
 *   This wrapper for Matlab was written by Shai Bagon (shai.bagon@weizmann.ac.il).
 *   Department of Computer Science and Applied Mathmatics
 *   Wiezmann Institute of Science
 *   http://www.wisdom.weizmann.ac.il/
 *
 *   The core cpp application was written by Veksler Olga:
 * 	  
 *   [1] Efficient Approximate Energy Minimization via Graph Cuts
 *       Yuri Boykov, Olga Veksler, Ramin Zabih,
 *       IEEE transactions on PAMI, vol. 20, no. 12, p. 1222-1239, November 2001.
 * 
 *   [2] What Energy Functions can be Minimized via Graph Cuts?
 *       Vladimir Kolmogorov and Ramin Zabih.
 *       To appear in IEEE Transactions on Pattern Analysis and Machine Intelligence (PAMI).
 *       Earlier version appeared in European Conference on Computer Vision (ECCV), May 2002.
 * 
 *   [3] An Experimental Comparison of Min-Cut/Max-Flow Algorithms
 *       for Energy Minimization in Vision.
 *       Yuri Boykov and Vladimir Kolmogorov.
 *       In IEEE Transactions on Pattern Analysis and Machine Intelligence (PAMI),
 *       September 2004
 *
 *   [4] Matlab Wrapper for Graph Cut.
 *       Shai Bagon.
 *       in www.wisdom.weizmann.ac.il/~bagon, December 2006.
 *     
 * 	This software can be used only for research purposes, you should cite
 * 	the aforementioned papers in any resulting publication.
 * 	If you wish to use this software (or the algorithms described in the aforementioned paper)
 * 	for commercial purposes, you should be aware that there is a US patent: 
 * 
 * 		R. Zabih, Y. Boykov, O. Veksler, 
 * 		"System and method for fast approximate energy minimization via graph cuts ", 
 * 		United Stated Patent 6,744,923, June 1, 2004
 *
 *
 *   The Software is provided "as is", without warranty of any kind.
 *
 *
 */

void GraphCutMex(GCoptimization *MyGraph, char mode, unsigned short* seg_im, int iter)
{    
    if ( ! MyGraph->IsClassValid() ) {
		std::cout<<"GraphCut:handle GC handle is not valid\n";
    }
 
    switch (mode) {
        case 'c': /* close */
            delete MyGraph; /* ->~GCoptimization(); explicitly call the destructor */
            MyGraph = NULL;
            break;
        case 's': /* set labels */
            /* setting the labels: we have an extra parameter - int32 array of size w*l
             * containing the new labels */
            //SetLabels(MyGraph, nlhs, plhs, nrhs, prhs);
			SetLables(MyGraph, seg_im);
            break;
        //case 'g': /* get labels */
        //    GetLabels(MyGraph, nlhs, plhs, nrhs, prhs);
        //    break;
        //case 'a': /* a-b swaps */
        //    ABswaps(MyGraph, nlhs, plhs, nrhs, prhs);
        //    break;
        case 'e': /* expand steps */
            //Expand(MyGraph, nlhs, plhs, nrhs, prhs);
			if(iter==0)				
				Expand(MyGraph, seg_im, 2, iter);
			else
				Expand(MyGraph, seg_im, 3, iter);
            break;
        //case 'n': /* get the current labeling energy */
        //    Energy(MyGraph, nlhs, plhs, nrhs, prhs);
        //    break;
        default:
            //mexErrMsgIdAndTxt("GraphCut:mode","unrecognized mode");
			std::cout<<"unrecognized mode\n";
    }
}



/**************************************************************************************/
/* support several kinds of expansion 
 * 1. default - expand untill convergence no extra parameters
 * 2. #iterations - performs #iterations:   GraphCut(gch, mode, #iteration)
 * 3. expand label - expand specific label  GraphCut(gch, mode, [], label)
 * 4. expand label at specific indices      GraphCut(gch, mode, [], label, indices) indices start with zero and not 1 like usuall matlab inices!!
 */
void Expand(GCoptimization *MyGraph, unsigned short* seg_im, int nrhs, int max_iterations)
{   
    //int num(0), max_iterations(0);
    //GCoptimization::LabelType label(0);
    //GCoptimization::PixelType *indices = NULL;
    
    /* check inputs */
    switch (nrhs) {
        case 2:
            /* default - expand till convergence */
            MyGraph->expansion();
            break;
        case 3:   
            if ( max_iterations == 1 ) 
                MyGraph->oneExpansionIteration();
            else
                MyGraph->expansion(max_iterations);
            break;                
        default:        
			std::cout<<"GraphCut:Expand incorrect number of input arguments\n";			
    }
    
	MyGraph->ExportLabels( (GCoptimization::LabelType*) seg_im );
}

/**************************************************************************************/
/* set user defined labels to graph */
void SetLables(GCoptimization *MyGraph, unsigned short* seg_im)
{               
    MyGraph->SetAllLabels( (GCoptimization::LabelType*) seg_im );
}
