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

    
    MyGraph = new GCoptimization(R*C*Z, num_labels, SET_ALL_AT_ONCE, SET_ALL_AT_ONCE);
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
    

