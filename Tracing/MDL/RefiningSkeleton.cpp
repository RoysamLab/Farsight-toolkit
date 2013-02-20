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

// New Skeleton 
// Input format: .backbone spine and ExtraSpine  (3D points)
// Input format: .vtk  (3D skeleton file)
// Author: Liang XIAO, RPI
// Date: 12/Nov/2009


#if defined(_MSC_VER)
#pragma warning(disable : 4996)
#endif

#include <cstring>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/algo/vnl_matrix_inverse.h>

using std::cerr;
using std::cout;
using std::endl;

#define Dimention  3;

struct  VoxelPosition
{
  float x;
  float y;
  float z;
}; 

#define DATATYPEIN unsigned char

int main(int argc, char **argv)
{
  //ifstream fin;
  FILE *smoothbackbone = 0;
  FILE *spine = 0;
  FILE *ExtraSpine = 0;
  FILE *outrefineskel = 0;
  FILE *outseed = 0;
  
  std::string filedir;
  std::string infilename;
  //file names
  std::string spineSkeletonFile;
  std::string extraSpineFile;
  std::string skeletonFile;
  std::string seedFile;

  VoxelPosition *AllBackbonepoints;
  VoxelPosition *AllSpinepoints = 0;
  VoxelPosition *AllExtraSpinepoints = 0;

  bool *RealPointsID = 0;
  bool *RealSpineID = 0;
  bool *ExtraSpineID = 0;
  
  int NumBackbonePoints;
  int NumExtraSpinePoints = 0;
  int NumSpinePoints = 0;
  
  int BackboneOnly = 0; 

  int i;
  int tmp1, tmp2,tmp3;
  int NumLines, SkipNmuber;
  
  if (argc < 7)
    {
    cerr << "Usage: " << argv[0] << " <dir> <smooth backbone Skeleton file>"
         << "  <Spine skeleton file> <ExtraSpine file> <out refine seed file>"
         << "<out refine skeleton file>  <OnlyBackbone (0 or 1)>" << endl;
    return 1;
    }

  filedir = argv[1];

  // open smoothbackbone file
  infilename = filedir + argv[2];
  if(( smoothbackbone=fopen(infilename.c_str(), "rb")) == NULL)  
    {
    cerr << "couldn't open smoothbackbone file " << infilename << " for input" << endl;
    return -1;
    }
  
  BackboneOnly =  atoi(argv[7]);
  if(!BackboneOnly)
    {

    // open spineskeleton file
    spineSkeletonFile = filedir + argv[3];
    if((spine=fopen(spineSkeletonFile.c_str(), "rb")) == NULL)
      {
      cerr << "couldn't open spine skeleton file " << spineSkeletonFile
           << " for input" << endl;
      return -1;
      }
 
    // open extra spine file
    extraSpineFile = filedir + argv[4];
    if((ExtraSpine=fopen(extraSpineFile.c_str(), "rb")) == NULL)  
      {
      cerr << "couldn't open Extra spine skeleton file " << extraSpineFile
           << " for input" << endl;
      return -1;
      }

    }//end if (!BackboneOnly)

  // open seed output file
  seedFile = filedir + argv[5];
  if((outseed=fopen(seedFile.c_str(), "w")) == NULL) 
    {
    cerr << "couldn't open seed file " << seedFile << " for output" << endl;
    return -1;
    }

  // open skeleton output file
  skeletonFile = filedir + argv[6];
  if((outrefineskel=fopen(skeletonFile.c_str(), "wb")) == NULL)  
    {
    cerr << "couldn't open skeleton file " << skeletonFile << " for output"
         << endl;
    return -1;
    }
   
  //------------------Read From smooth backbone skeleton VTK file ---------------------//

  // Skip first 4 lines of skeleton VTK file
  char str[200]; 
  for (i=0;i<12;i++)
    {
    if( fscanf(smoothbackbone,"%s",str) == EOF )
      {
      cerr << "End of file encountered in fscanf!" << endl;
      }
    }
 
  if( fscanf(smoothbackbone,"%s",str) == EOF )
    {
    cerr << "End of file encountered in fscanf!" << endl;
    }
  NumBackbonePoints = atoi(str);
  cout << "There are " << NumBackbonePoints << " skeleton points" << endl;
  if( fscanf(smoothbackbone,"%s",str) == EOF )
    {
    cerr << "End of file encountered in fscanf!" << endl;
    }

  AllBackbonepoints = new VoxelPosition[ NumBackbonePoints];
  float temp;  
  for (i=0;i<NumBackbonePoints;i++)
    {
     if( fscanf (smoothbackbone,"%f",&temp) == EOF ) 
       {
       cerr << "End of file encountered in fscanf!" << endl;
       }
     AllBackbonepoints[i].x = temp;

     if( fscanf (smoothbackbone,"%f",&temp) == EOF ) 
       {
       cerr << "End of file encountered in fscanf!" << endl;
       }
     AllBackbonepoints[i].y = temp; // change the sign,

     if( fscanf (smoothbackbone,"%f",&temp) == EOF ) 
       {
       cerr << "End of file encountered in fscanf!" << endl;
       }
     AllBackbonepoints[i].z =temp;
    }
 
  //- choose the real backpoints from the Backbone.VTK file 

  RealPointsID = new bool[ NumBackbonePoints];
  for (i=0;i<NumBackbonePoints;i++)
  {
   RealPointsID[i]=false;// init flag
  }
 
  if( fscanf(smoothbackbone,"%s",str) == EOF ) // read LINES
    {
    cerr << "End of file encountered in fscanf!" << endl;
    }
  if( fscanf(smoothbackbone,"%d",&NumLines) == EOF )
    {
    cerr << "End of file encountered in fscanf!" << endl;
    }
  // skip the last number in this line
  if( fscanf(smoothbackbone,"%d",&SkipNmuber) == EOF )
    {
    cerr << "End of file encountered in fscanf!" << endl;
    }
  for (i=0;i<NumLines;i++)
    {
    // Read one line into an array;
    // skip first number in the Line
    if( fscanf (smoothbackbone,"%d",&tmp1) == EOF )
      {
      cerr << "End of file encountered in fscanf!" << endl;
      }
    // skip second number in the Line
    if( fscanf (smoothbackbone,"%d",&tmp2) == EOF )
      {
      cerr << "End of file encountered in fscanf!" << endl;
      }
    // skip third number in the Line
    if( fscanf (smoothbackbone,"%d",&tmp3) == EOF )
      {
      cerr << "End of file encountered in fscanf!" << endl;
      }
    RealPointsID[tmp2]= true;
	  RealPointsID[tmp3]= true;
    }
  fclose (smoothbackbone);

//------------------Read From spine skeleton VTK file ---------------------//


   // Skip first 4 lines of skeleton VTK file
 if(!BackboneOnly)
  {
   for (i=0;i<12;i++)
    {
    if( fscanf(spine,"%s",str) == EOF )
      {
      cerr << "End of file encountered in fscanf!" << endl;
      }
    }
   if( fscanf(spine,"%s",str) == EOF )
     {
     cerr << "End of file encountered in fscanf!" << endl;
     }
   NumSpinePoints = atoi(str);
   cout << "There are " << NumSpinePoints << " spine skeleton points" << endl;
   if( fscanf(spine,"%s",str) == EOF )
     {
     cerr << "End of file encountered in fscanf!" << endl;
     }
   AllSpinepoints = new VoxelPosition[NumSpinePoints];
   float temp3;  
   for (i=0;i<NumSpinePoints;i++)
     {
     if( fscanf(spine,"%f",&temp3) == EOF )
       {
       cerr << "End of file encountered in fscanf!" << endl;
       }
     AllSpinepoints[i].x = temp3;
     if( fscanf(spine,"%f",&temp3) == EOF )
       {
       cerr << "End of file encountered in fscanf!" << endl;
       }
     AllSpinepoints[i].y = temp3; // change the sign,
     if( fscanf(spine,"%f",&temp3) == EOF )
       {
       cerr << "End of file encountered in fscanf!" << endl;
       }
     AllSpinepoints[i].z =temp3;
   }

   RealSpineID =  new bool [NumSpinePoints];
   for (i=0;i<NumSpinePoints;i++)
  {
   RealSpineID[i]=false;// init flag
  }
 
  if( fscanf(spine,"%s",str) == EOF ) // read LINES
    {
    cerr << "End of file encountered in fscanf!" << endl;
    }
  if( fscanf(spine,"%d",&NumLines) == EOF )
    {
    cerr << "End of file encountered in fscanf!" << endl;
    }
  // skip the last number in this line
  if( fscanf(spine,"%d",&SkipNmuber) == EOF ) 
    {
    cerr << "End of file encountered in fscanf!" << endl;
    }
  for (i=0;i<NumLines;i++)
    {
    // Read one line into an array;
    if( fscanf (spine,"%d",&tmp1) == EOF ) // skip first number in the Line
      {
      cerr << "End of file encountered in fscanf!" << endl;
      }
    if( fscanf (spine,"%d",&tmp2) == EOF ) // read second number in the line 
      {
      cerr << "End of file encountered in fscanf!" << endl;
      }
    if( fscanf (spine,"%d",&tmp3) == EOF ) // read third number in the line 
      {
      cerr << "End of file encountered in fscanf!" << endl;
      }
    RealSpineID[tmp2]= true;
	  RealSpineID[tmp3]= true;
    }
   fclose (spine);

//------------------Read From Extra-spine skeleton VTK file ---------------------//

 // Skip first 4 lines of skeleton VTK file
  //rewind (ExtraSpine);

   for (i=0;i<12;i++)
    {
    if( fscanf(ExtraSpine,"%s",str) == EOF )
      {
      cerr << "End of file encountered in fscanf!" << endl;
      }
    }
 
   if( fscanf(ExtraSpine,"%s",str) == EOF )
     {
     cerr << "End of file encountered in fscanf!" << endl;
     }
   NumExtraSpinePoints = atoi(str);
   if( fscanf(ExtraSpine,"%s",str) == EOF )
     {
     cerr << "End of file encountered in fscanf!" << endl;
     }
   AllExtraSpinepoints = new VoxelPosition[NumExtraSpinePoints];
   for (i=0;i<NumExtraSpinePoints;i++)
     {
     if( fscanf (ExtraSpine,"%f",&temp) == EOF )
       {
       cerr << "End of file encountered in fscanf!" << endl;
       }
     AllExtraSpinepoints[i].x = temp;
     if( fscanf (ExtraSpine,"%f",&temp) == EOF )
       {
       cerr << "End of file encountered in fscanf!" << endl;
       }
     AllExtraSpinepoints[i].y = temp; // change the sign,
     if( fscanf (ExtraSpine,"%f",&temp) == EOF )
       {
       cerr << "End of file encountered in fscanf!" << endl;
       }
     AllExtraSpinepoints[i].z =temp;
     }

   ExtraSpineID =  new bool [NumExtraSpinePoints];
   for (i=0;i<NumExtraSpinePoints;i++)
     {
     ExtraSpineID[i]=false;// init flag
     }
 
  if( fscanf(ExtraSpine,"%s",str) == EOF ) // read LINES
    {
    }
  if( fscanf(ExtraSpine,"%d",&NumLines) == EOF )
    {
    }
  // skip the last number in this line
  if( fscanf(ExtraSpine,"%d",&SkipNmuber) == EOF ) 
    {
    }
  for (i=0;i<NumLines;i++)
    {
    // Read one line into an array;
    if( fscanf (ExtraSpine,"%d",&tmp1) == EOF ) // skip first number in the Line
      {
      cerr << "End of file encountered in fscanf!" << endl;
      }
    if( fscanf (ExtraSpine,"%d",&tmp2) == EOF ) // read second number in the line
      {
      cerr << "End of file encountered in fscanf!" << endl;
      }
    if( fscanf (ExtraSpine,"%d",&tmp3) == EOF ) // read third number in the line
      {
      cerr << "End of file encountered in fscanf!" << endl;
      }
    ExtraSpineID[tmp2]= true;
	ExtraSpineID[tmp3]= true;

    }
   fclose (ExtraSpine);
  } // end if 
 //---------------------------Output to Skel File
  for (i=0;i<NumBackbonePoints;i++)
  { if (RealPointsID[i]) // if it is a real backbone points
    fprintf(outrefineskel,"%f %f %f %d\n", AllBackbonepoints[i].x, AllBackbonepoints[i].y , AllBackbonepoints[i].z,1);
    fprintf(outseed,"%d %d %d\n", (int) AllBackbonepoints[i].x, (int) AllBackbonepoints[i].y , (int) AllBackbonepoints[i].z);
  }
  if(!BackboneOnly)
  { 
    for (i=0;i<NumSpinePoints;i++)
	{if(RealSpineID[i])
     fprintf(outrefineskel,"%f %f %f %d\n", AllSpinepoints[i].x, AllSpinepoints[i].y , AllSpinepoints[i].z,1);
	 //fprintf(outseed,"%d %d %d\n", (int) AllBackbonepoints[i].x, (int) AllBackbonepoints[i].y , (int) AllBackbonepoints[i].z);
	}
    for (i=0;i<NumExtraSpinePoints;i++)
	{if(ExtraSpineID[i])
     fprintf(outrefineskel,"%f %f %f %d\n", AllExtraSpinepoints[i].x, AllExtraSpinepoints[i].y , AllExtraSpinepoints[i].z,1);
	 //fprintf(outseed,"%d %d %d\n", (int) AllBackbonepoints[i].x, (int) AllBackbonepoints[i].y , (int) AllBackbonepoints[i].z);
	}
	
  } // end if 
  fclose (outrefineskel);
  fclose (outseed);
  //release memnory
  delete []AllBackbonepoints;
  delete []RealPointsID;
  if(!BackboneOnly)
  { delete []AllSpinepoints;
    delete []AllExtraSpinepoints;
	delete []RealSpineID;
	delete []ExtraSpineID;
  }// end if
 
}
