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

/*void GVC::sequential_coloring( unsigned char *bufferOut, unsigned char *img,\
                               int rows, int cols, int ncolors,\
                               vector< vector<int> > &vec,
                               vector< vector<int> > &map )
*/
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


/*----------------------------------------------------------------------------*/
/* function : solve, valid, goal                                              */
/*                                                                            */
/* Description : These are the three functions that are used to iteratively   */
/* color vertices such that the four-color criterion is satisfied. This is a  */
/* recursive backtracking algorithm. The original code was presented in       */
/* a set of class notes, in JAVA, by R.W.Topor at                             */
/* http://www.cit.gu.edu.au/~rwt/p2.02.1/                                     */
/* in the file "recursion.pdf"                                                */
/*                                                                            */
/*  Author:                                                                   */
/*                                                                            */
/*    Sumit K Nath,                                                           */
/*    Department of Computer Science,                                         */
/*    University of Missouri-Columbia,                                        */
/*    Columbia, Missouri, USA 65211                                           */
/*                                                                            */
/*----------------------------------------------------------------------------*/
bool GVC::solve(vector< vector<int> > &map, int v, int ncolors)
{
    int nobjects = map.size();
    
    
    
    if (goal(v, nobjects))
    {
        return true;
    }
    else
    {
        /*try to color vertex v with each color
          c in turn*/
        for (int c = 0; c < ncolors; c++)
        {
            if (valid(map, v, c))
            {
                ColorClass[v] = c;
                if (solve(map, (v+1), ncolors)) return true;
            }
        }
        return false;
    }
}

bool GVC::valid(vector< vector<int> > &map, int v, int c)
{
    for (int i = 0; i < (int)map[v].size(); i++)
    {
        int u = map[v][i]-1;
        if (u != v && u < v && ColorClass[u] == c) return false;
    }
    return true;
}


bool GVC::goal(int v, int nobjects)
{
    return (v==nobjects);
}

