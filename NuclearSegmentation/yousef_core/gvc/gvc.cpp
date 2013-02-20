#include "gvc.h"

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
bool GVC::solve(std::vector< std::vector<int> > &map, int v, int ncolors)
{
    int nobjects = (int)map.size();
    
    
    
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

bool GVC::valid(std::vector< std::vector<int> > &map, int v, int c)
{
    for (unsigned int i = 0; i < map[v].size(); i++)
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

