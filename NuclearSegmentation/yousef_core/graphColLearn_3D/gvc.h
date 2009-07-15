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

using namespace std;

class GVC
{
    public:
        GVC(){}
       ~GVC(){}
                      
        //Sequential algorithm
//         void  sequential_coloring( unsigned char *bufferOut,
//                                    unsigned char *img,
//                                    int rows, int cols, int ncolors,
//                                    vector< vector<int> > &vec,
//                                    vector< vector<int> > &map );
          void sequential_coloring(int nobjects, int ncolors, int* ColorOut, vector< vector<int> > &map );
    private:
        int  *ColorClass;
        
        bool  solve(vector< vector<int> > &map, int v, int ncolors);
        bool  valid(vector< vector<int> > &map, int v, int c);
        bool  goal(int v, int nobjects);
};

#endif
