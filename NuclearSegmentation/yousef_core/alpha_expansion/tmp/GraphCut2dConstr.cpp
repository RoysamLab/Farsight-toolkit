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

///////////////////////////////////////////////////////////////////
//This file was modified by Yousef Al-Kofahi (RPI) on May 28th 2008
//The function was intended to be used in Matlab as a Mex file 
//Details from the original Authors are below
////////////////////////////////////////////////////////////////////

//#include "mex.h"
#include "GCoptimization.h"
#include "GraphCut.h"
#include <stdlib.h>

/* Defines */


/*
 * Matlab wrapper for Weksler graph cut implementation
 *
 * usage:
 * [gch] = GraphCutConstr(width, height, num_labels, DataCost, SmoothnessCost,[vCost,hCost])
 *
 * Note that data types are crucials!
 * 
 * Inputs:
 *  width, hieght - 2D grid dimensions. (do not cast to int)
 *  num_labels - number of labels.  (do not cast to int)
 *  DataCost - of type float, array size [width*height*#labels], the data term for pixel
 *             x,y recieving label l is stroed at [(x+y*width)*#labels + l]
 *  SmoothnessCost - of type float, array size [#labels.^2] the cost between l1 and l2 is
 *                   stored at [l1+l2*#labels] = Vpq(lp,lq)
 *  vCost, hCost - of type float, array size[width*height] each, 
 *                 horizontal and vertical cues for smoothness term
 *
 * Outputs:
 *  gch - of type int32, graph cut handle - do NOT mess with it!
 */
//void mexFunction(
//    int		  nlhs, 	/* number of expected outputs */
//    mxArray	  *plhs[],	/* mxArray output pointer array */
//    int		  nrhs, 	/* number of inputs */
//    const mxArray	  *prhs[]	/* mxArray input pointer array */
//    )

GCoptimization * GraphCutConstr(float* DataCostIn, float* SmoothnessCostIn, float* hCueIn, float* vCueIn, int RIn, int CIn, int num_labels)
{
    
    GCoptimization::PixelType width, height;
    //int num_labels;
    Graph::captype *DataCost;
    Graph::captype *SmoothnessCost;
    Graph::captype *hCue = NULL;
    Graph::captype *vCue = NULL;
    GCoptimization::LabelType *Labels;
    GCoptimization *MyGraph = NULL;
            
	width = CIn;        
	height = RIn;
        
    DataCost = (Graph::captype*)DataCostIn;   
	SmoothnessCost = (Graph::captype*)SmoothnessCostIn;

	vCue = (Graph::captype*)vCueIn;
    hCue = (Graph::captype*)hCueIn;
    
    MyGraph = new GCoptimization(width, height, num_labels, SET_ALL_AT_ONCE, SET_ALL_AT_ONCE);
    MyGraph->setData(DataCost);
    if ( vCue != NULL && vCue != NULL ) 
        MyGraph->setSmoothness(SmoothnessCost, hCue, vCue);
    else
        MyGraph->setSmoothness(SmoothnessCost);
        
	return MyGraph;
}
    

