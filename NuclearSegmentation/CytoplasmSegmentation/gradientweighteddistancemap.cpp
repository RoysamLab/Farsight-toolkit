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

/*=========================================================================

  Program:   Farsight Biological Image Segmentation and Visualization Toolkit
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision: 0.00 $

=========================================================================*/
#ifndef __GradientEnhancedDistanceMapImageFilter_cxx
#define __GradientEnhancedDistanceMapImageFilter_cxx

#define Inp_Im(a,b) INP_IM_2D[(a)+(b)*size1]
#define grad_imw(a,b) GRAD_IMW[(a+1)+(b+1)*(size1+2)]

#include <iostream>
#include <queue>
#include <algorithm>
#include <float.h>
#include <math.h>

//This function takes in an image(array GRAD_IMW) with the background pixels at -inf,
// the pixels from which the gradient enhanced distance map is to be computed are initialized to 0,
// and the pixels where they have to be computed are at +inf or some really really large number
// and the gradient weighted map is returned in the same array. The other input is a pre-smoothed
// grayscale image(array Inp_Im).
// Comment: Works well for 2D images. If extended for 3D Images, please a) use a mapped queue and
// b) change the floating point computations into integers(long integers??) a)->pixels
// are not iterated over multiple times and b)->the processing is faster
// If you have the time, formulate it as a graph building problem and use an open src implementation
// Dijkstra's algorithm
// Original:Kedar 06/10/09
// Last edited:Kedar 10/30/09

struct im_ind{
	int i;
	int j;
};

int gradient_enhanced_distance_map( float *INP_IM_2D, float *GRAD_IMW, int size1, int size2 ){

//Copy data type floats' lower limit into a and get root2 into a local variable
	float flot_mini,root2;

	flot_mini = FLT_MAX*-1;
	root2 = (float)sqrt(2.0);

	std::queue<im_ind> icy_needs_a_change;

//Initialize queue with pixels that are neighbours of the pixels at the initializing regions
	for( int j=1; j<size2-1; j++ )
		for( int i=1; i<size1-1; i++ )
			if(	( grad_imw(i-1,j-1) >= 0  && grad_imw(i,j) > grad_imw(i-1,j-1) ) || \
				( grad_imw(i-1,j)   >= 0 && grad_imw(i,j) > grad_imw(i-1,j)   ) || \
				( grad_imw(i-1,j+1) >= 0 && grad_imw(i,j) > grad_imw(i-1,j+1) ) || \
				( grad_imw(i,j-1)   >= 0 && grad_imw(i,j) > grad_imw(i,j-1)   ) || \
				( grad_imw(i,j+1)   >= 0 && grad_imw(i,j) > grad_imw(i,j+1)   ) || \
				( grad_imw(i+1,j-1) >= 0 && grad_imw(i,j) > grad_imw(i+1,j-1) ) || \
				( grad_imw(i+1,j)   >= 0 && grad_imw(i,j) > grad_imw(i+1,j)   ) || \
				( grad_imw(i+1,j+1) >= 0 && grad_imw(i,j) > grad_imw(i+1,j+1) ) )
			{
					im_ind temp_inds;
					temp_inds.i = i;
					temp_inds.j = j;
					icy_needs_a_change.push( temp_inds );
			}

	im_ind temp_inds;
	int i,j;
	float neigh_vals[8],neigh_vals_cpy[8];
	int best_val_not_found;

//Pop and push pixels indices into a queue till all the pixels are set, i.e, queue is empty
	while( !icy_needs_a_change.empty() ){
		best_val_not_found = 1;
		temp_inds = icy_needs_a_change.front();
		icy_needs_a_change.pop();

		i = temp_inds.i; j = temp_inds.j;

		if( i == 0 || j == 0 || i == (size1-1) || j == (size2-1) ){
			continue;
		}
		
		//Calculate the weighted distances from the neighbors
		neigh_vals[0] = (fabs(Inp_Im(i-1,j-1)-Inp_Im(i,j))+1)*root2+grad_imw(i-1,j-1);
		neigh_vals[1] = (fabs(Inp_Im(i-1,j)-Inp_Im(i,j))+1)+grad_imw(i-1,j);
		neigh_vals[2] = (fabs(Inp_Im(i-1,j+1)-Inp_Im(i,j))+1)*root2+grad_imw(i-1,j+1);
		neigh_vals[3] = (fabs(Inp_Im(i,j-1)-Inp_Im(i,j))+1)+grad_imw(i,j-1);
		neigh_vals[4] = (fabs(Inp_Im(i,j+1)-Inp_Im(i,j))+1)+grad_imw(i,j+1);
		neigh_vals[5] = (fabs(Inp_Im(i+1,j-1)-Inp_Im(i,j))+1)*root2+grad_imw(i+1,j-1);
		neigh_vals[6] = (fabs(Inp_Im(i+1,j)-Inp_Im(i,j))+1)+grad_imw(i+1,j);
		neigh_vals[7] = (fabs(Inp_Im(i+1,j+1)-Inp_Im(i,j))+1)*root2+grad_imw(i+1,j+1);

		//Duplicate and order the list
		for( int k=0; k<8; k++ ) neigh_vals_cpy[k] = neigh_vals[k];
		int num_elements = sizeof(neigh_vals_cpy) / sizeof(neigh_vals_cpy[0]); 
		std::sort(neigh_vals_cpy, neigh_vals_cpy + num_elements);

		for( int k=0; k<8; k++ ){
			//Skip if the computed gradient enhanced distance is either from the background or it is a repeated value or the change is less than a threshold
			if( neigh_vals_cpy[k] <= 0 ) continue;
			//At this point the negative vals have been skipped and the lowest +ve value is the current value on the sorted list
			if( best_val_not_found ){
				if( ( grad_imw(i,j) > neigh_vals_cpy[k] ) && ( (fabs(grad_imw(i,j)-neigh_vals_cpy[k])) > 0.01 )  ){
					grad_imw(i,j) = neigh_vals_cpy[k];
					best_val_not_found = 0;
					continue;
				}
				else break;
			}

			if( (neigh_vals_cpy[k] ==  neigh_vals[0]) ){
				temp_inds.i = i-1; temp_inds.j = j-1;
				icy_needs_a_change.push( temp_inds );
			}
			if( neigh_vals_cpy[k] ==  neigh_vals[1] ){
				temp_inds.i = i-1; temp_inds.j = j;
				icy_needs_a_change.push( temp_inds );
			}
			if( neigh_vals_cpy[k] ==  neigh_vals[2] ){
				temp_inds.i = i-1; temp_inds.j = j+1;
				icy_needs_a_change.push( temp_inds );
			}
			if( neigh_vals_cpy[k] ==  neigh_vals[3] ){
				temp_inds.i = i; temp_inds.j = j-1;
				icy_needs_a_change.push( temp_inds );
			}
			if( neigh_vals_cpy[k] ==  neigh_vals[4] ){
				temp_inds.i = i; temp_inds.j = j+1;
				icy_needs_a_change.push( temp_inds );
			}
			if( neigh_vals_cpy[k] ==  neigh_vals[5] ){
				temp_inds.i = i+1; temp_inds.j = j-1;
				icy_needs_a_change.push( temp_inds );
			}
			if( neigh_vals_cpy[k] ==  neigh_vals[6] ){
				temp_inds.i = i+1; temp_inds.j = j;
				icy_needs_a_change.push( temp_inds );
			}
			if( neigh_vals_cpy[k] ==  neigh_vals[7] ){
				temp_inds.i = i+1; temp_inds.j = j+1;
				icy_needs_a_change.push( temp_inds );
			}
		}
	}
	return 1;
}


#endif
