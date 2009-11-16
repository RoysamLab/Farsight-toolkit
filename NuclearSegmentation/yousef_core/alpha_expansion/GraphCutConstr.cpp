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

#include "GCoptimization.h"
#include "GraphCut.h"
#include <stdlib.h>
#include <math.h>

/* Defines */
// define A64BITS for compiling on a 64 bits machine...

/*
 * Matlab wrapper for Weksler graph cut implementation
 *
 * usage:
 * [gch] = GraphCut3dConstr(R, C, Z, num_labels, DataCost, SmoothnessCost, [Contrast])
 *
 * Note that data types are crucials!
 * 
 * Inputs:
 *  R, C, Z - size of 3D grid. 
 *  num_labels - number of labels.  
 *  DataCost - of type float, array size [width*height*depth*#labels], the data term for pixel
 *             x,y,z recieving label l is stroed at [(x+y*width+z*width*depth)*#labels + l] //by yousef: I think depth should be replaced by height
 *  SmoothnessCost - of type float, array size [#labels.^2] the cost between l1 and l2 is
 *                   stored at [l1+l2*#labels] = Vpq(lp,lq)
 *  Contrast - of type float, array size[width*height*depth], the weight Wpq will be determined by the contrast:
 *                  Wpq = exp(-|C(p)-C(q)|)
 *
 * Outputs:
 *  gch - of type int32, graph cut handle - do NOT mess with it!
 */

GCoptimization * GraphCut3dConstr(float* ContrastIn, float* DataCostIn, float* SmoothnessCostIn, int RIn, int CIn, int ZIn, int num_labels)
{       
    Graph::captype *DataCost;
    Graph::captype *SmoothnessCost;
    Graph::captype *Contrast = NULL;
    //GCoptimization::LabelType *Labels;
    GCoptimization::PixelType R, C, Z; 
    
    GCoptimization *MyGraph = NULL;
           
	R = RIn;
	C = CIn;
	Z = ZIn;
       
	DataCost = (Graph::captype*)DataCostIn;
    
	SmoothnessCost = (Graph::captype*)SmoothnessCostIn;

	Contrast = (Graph::captype*)ContrastIn;
	//int uu1 = sizeof(ContrastIn[1]);
	//int uu2 = sizeof(Contrast[1]);
	
    
    //By Yousef: Estimate the number of needed edges
	long num_ed = 0;
	for ( int r = 0 ; r <= R - 2;  r++ )
	{
		for ( int c = 0 ; c <= C - 2; c++ )
		{
			for ( int z = 0 ; z <= Z - 2; z++ ) 
			{
				num_ed+=3;
			}
			num_ed+=2;
		}
		for ( int z = 0 ; z <= Z -2 ; z++ )
		{
			num_ed+=2;
		}
		num_ed+=1;
	}
	for ( int c = 0 ; c <= C - 2; c++ )
	{
		for ( int z = 0 ; z <= Z - 2; z++ )
		{
			num_ed+=2;
		}
		num_ed+=1;
	}
	for ( int z = 0 ; z <= Z - 2; z++ )
	{
		num_ed+=1;
	}

	MyGraph = new GCoptimization(R*C*Z, num_labels, SET_ALL_AT_ONCE, SET_ALL_AT_ONCE, num_ed);		

    /* neighborhod setup */
    GCoptimization::PixelType c(0), r(0), z(0), p(0), q(0);
    if (Contrast) {
		double sig = 30;//added by Yousef on 10-27-2008
        for ( r = 0 ; r <= R - 2;  r++ ) {
            for ( c = 0 ; c <= C - 2; c++ ) {
                for ( z = 0 ; z <= Z - 2; z++ ) {
                    /* all nodes with 3 nieghbors */					
                    p = r+c*R+z*R*C;/*c+r*C+z*R*C;*/
                    q = r+1+c*R+z*R*C;/*c+(r+1)*C+z*R*C;*/
                    MyGraph->setNeighbors(p, q, exp(-fabs(Contrast[p]-Contrast[q])/sig));
                    q = r+(c+1)*R+z*R*C;/*c+1+r*C+z*R*C;*/
                    MyGraph->setNeighbors(p, q, exp(-fabs(Contrast[p]-Contrast[q])/sig));
                    q = r+c*R+(z+1)*R*C;/*c+r*C+(z+1)*R*C;*/
                    MyGraph->setNeighbors(p, q, exp(-fabs(Contrast[p]-Contrast[q])/sig));
                }
               /* top of Z has 2 n in c and r */				
                p = r+c*R+z*R*C;/*c+r*C+z*R*C;*/				
                q = r+1+c*R+z*R*C;/*c+(r+1)*C+z*R*C;*/
                MyGraph->setNeighbors(p, q, exp(-fabs(Contrast[p]-Contrast[q])/sig));
                q = r+(c+1)*R+z*R*C;/*c+1+r*C+z*R*C;*/
                MyGraph->setNeighbors(p, q, exp(-fabs(Contrast[p]-Contrast[q])/sig));
            }
            /* end of c has 2 n in z and r */
            for ( z = 0 ; z <= Z -2 ; z++ ) {				
                p = r+c*R+z*R*C;/*c+r*C+z*R*C;*/
                q = r+1+c*R+z*R*C;/*c+(r+1)*C+z*R*C;*/
                MyGraph->setNeighbors(p, q, exp(-fabs(Contrast[p]-Contrast[q])/sig));
                q = r+c*R+(z+1)*R*C;/*c+r*C+(z+1)*R*C;*/
                MyGraph->setNeighbors(p, q, exp(-fabs(Contrast[p]-Contrast[q])/sig));
            }
            /* end of c abd z has n in r */
            p = r+c*R+z*R*C;/*c+r*C+z*R*C;*/
            q = r+1+c*R+z*R*C;/*c+(r+1)*C+z*R*C;*/
            MyGraph->setNeighbors(p, q, exp(-fabs(Contrast[p]-Contrast[q])/sig));
        }
        /* end of r has n in z and c */
        for ( c = 0 ; c <= C - 2; c++ ) {
            for ( z = 0 ; z <= Z - 2; z++ ) {
                p = r+c*R+z*R*C;/*c+r*C+z*R*C;*/
                q = r+c*R+(z+1)*R*C;/*c+r*C+(z+1)*R*C;*/
                MyGraph->setNeighbors(p, q, exp(-fabs(Contrast[p]-Contrast[q])));
                q = r+(c+1)*R+z*R*C;/*c+1+r*C+z*R*C;*/
                MyGraph->setNeighbors(p, q, exp(-fabs(Contrast[p]-Contrast[q])));
            }
            p = r+c*R+z*R*C;/*c+r*C+z*R*C;*/
            q = r+(c+1)*R+z*R*C;/*c+1+r*C+z*R*C;*/
            MyGraph->setNeighbors(p, q, exp(-fabs(Contrast[p]-Contrast[q])));
        }
        for ( z = 0 ; z <= Z - 2; z++ ) {
            p = r+c*R+z*R*C;/*c+r*C+z*R*C;*/
            q = r+c*R+(z+1)*R*C;/*c+r*C+(z+1)*R*C;*/
            MyGraph->setNeighbors(p, q, exp(-fabs(Contrast[p]-Contrast[q])));
        }
        /* end of graph construction with contrast */
    } else {
        for ( r = 0 ; r <= R - 2;  r++ ) {
            for ( c = 0 ; c <= C - 2; c++ ) {
                for ( z = 0 ; z <= Z - 2; z++ ) {
                /* all nodes with 3 nieghbors */
                    MyGraph->setNeighbors(r+c*R+z*R*C,r+1+c*R+z*R*C);
                    MyGraph->setNeighbors(r+c*R+z*R*C,r+(c+1)*R+z*R*C);
                    MyGraph->setNeighbors(r+c*R+z*R*C,r+c*R+(z+1)*R*C);
                }
            /* top of Z has 2 n in c and r */
                MyGraph->setNeighbors(r+c*R+z*R*C,r+1+c*R+z*R*C);
                MyGraph->setNeighbors(r+c*R+z*R*C,r+(c+1)*R+z*R*C);
            }
        /* end of c has 2 n in z and r */
            for ( z = 0 ; z <= Z -2 ; z++ ) {
                MyGraph->setNeighbors(r+c*R+z*R*C,r+1+c*R+z*R*C);
                MyGraph->setNeighbors(r+c*R+z*R*C,r+c*R+(z+1)*R*C);
            }
        /* end of c abd z has n in r */
            MyGraph->setNeighbors(r+c*R+z*R*C,r+1+c*R+z*R*C);
        }
    /* end of r has n in z and c */
        for ( c = 0 ; c <= C - 2; c++ ) {
            for ( z = 0 ; z <= Z - 2; z++ ) {
                MyGraph->setNeighbors(r+c*R+z*R*C,r+c*R+(z+1)*R*C);
                MyGraph->setNeighbors(r+c*R+z*R*C,r+(c+1)*R+z*R*C);
            }
            MyGraph->setNeighbors(r+c*R+z*R*C,r+(c+1)*R+z*R*C);
        }
        for ( z = 0 ; z <= Z - 2; z++ ) {
            MyGraph->setNeighbors(r+c*R+z*R*C,r+c*R+(z+1)*R*C);
        }
    }
    MyGraph->setData(DataCost);
    MyGraph->setSmoothness(SmoothnessCost);
           
	return MyGraph;
}
    
GCoptimization * GraphCutConstr(float* DataCostIn, float* SmoothnessCostIn, float* hCueIn, float* vCueIn, int RIn, int CIn, int num_labels)
{    
    GCoptimization::PixelType width, height;
    //int num_labels;
    Graph::captype *DataCost;
    Graph::captype *SmoothnessCost;
    Graph::captype *hCue = NULL;
    Graph::captype *vCue = NULL;
    //GCoptimization::LabelType *Labels;
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
