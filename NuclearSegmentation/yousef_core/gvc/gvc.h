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

#ifndef _GVC_H_
#define _GVC_H_

#include <iostream>
#include <cstdlib>
#include <vector>

#define MAX_RAND (2.0*(1 << 30))

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef INF
#define INF 100000.0
#endif

class GVC
{
    public:
        GVC(){}
       ~GVC(){}
                      
        //Sequential algorithm
/*         void  sequential_coloring( unsigned char *bufferOut,\
                                    unsigned char *img,\
                                    int rows, int cols, int ncolors,\
                                    vector< vector<int> > &vec,
                                    vector< vector<int> > &map );
*/
          void sequential_coloring(int nobjects, int ncolors, int* ColorOut, std::vector< std::vector<int> > &map );
    private:
        int  *ColorClass;
        
        bool  solve(std::vector< std::vector<int> > &map, int v, int ncolors);
        bool  valid(std::vector< std::vector<int> > &map, int v, int c);
        bool  goal(int v, int nobjects);
};

#endif
