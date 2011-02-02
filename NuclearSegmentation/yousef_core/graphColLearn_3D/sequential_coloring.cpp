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

/*----------------------------------------------------------------------------*/
/* function : sequential_coloring                                             */
/*                                                                            */
/* Description : This function sequentially colors the vertices of a planar   */
/* graph, whose vertices are passed through the "vector, vec"                 */
/*                                                                            */
/*  Author:                                                                   */
/*                                                                            */
/*    Sumit K Nath,                                                           */
/*    Department of Computer Science,                                         */
/*    University of Missouri-Columbia,                                        */
/*    Columbia, Missouri, USA 65211                                           */
/*                                                                            */
/*                                                                            */
/* Input: bufferOut(the previous mask as an unsigned char)                    */
/*        rows (#of rows of the image)                                        */
/*        vec (a vector, containing statistics of the binary mask)            */ 
/* Output: img (The new re-labeled image as an unsigned char)                 */ 
/*----------------------------------------------------------------------------*/
#include "gvc.h"

//void GVC::sequential_coloring( unsigned char *bufferOut, unsigned char *img,
//                               int rows, int cols, int ncolors,
//                               vector< vector<int> > &vec,
//                               vector< vector<int> > &map )
void GVC::sequential_coloring(int nobjects, int ncolors, int* ColorOut, vector< vector<int> > &map )                                                                                             
{
    //int nobjects = vec.size();
    
    /*allocate memory for colors*/
    ColorClass = new int[nobjects];
    
    for (int i = 0; i < nobjects; ++i) ColorClass[i] = 0;
    
    /*group all objects into "ncolors" level sets*/
    solve( map, 0, ncolors );

	for (int i = 0; i < nobjects; ++i) ColorOut[i] = ColorClass[i];
    
    ///*clear out contents of img*/
    //for (int i = 0; i < rows*cols; ++i) img[i] = 0;
    //
    ///*now re-color img, with new colors stored in "colors"*/
    //for (int i = 0; i < nobjects; ++i)
    //{
    //    for (int y = vec[i][0]; y <= vec[i][2]; ++y)
    //    {
    //        for (int x = vec[i][1]; x <= vec[i][3]; ++x)
    //        {
    //            /*check the label of the connected component.If 
    //              equal then assign all pixels to a color from  
    //              the 'colors' array*/
    //            int index = y + x *rows;
    //            
    //            if ( bufferOut[index] == (i+1))
    //            {
    //                img[index] = (ColorClass[i]+1);
    //            }
    //        }
    //    }
    //}
    
    delete ColorClass;
}


