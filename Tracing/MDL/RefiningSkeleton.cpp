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
#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/algo/vnl_matrix_inverse.h>


using namespace std;
//using namespace boost;

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
  FILE  *smoothbackbone,*spine, *ExtraSpine;
  FILE  *outrefineskel;
  
  std::string filedir;
  std::string infilename;
  std::string tempfile1; //  temp file, 
  std::string tempfile2; //  temp file, 
  std::string tempfile3; //  temp file, 

  VoxelPosition *AllBackbonepoints;
  VoxelPosition *AllSpinepoints;
  VoxelPosition *AllExtraSpinepoints;
  
  
  int i;
  
  if (argc < 5)
    {
    if(argc > 1 && strcmp(argv[1], "debug") == 0)
      {
      //---- for debug test----//
      argv[0]="RefiningSkeleton.exe";
      argv[1]="C:\\Farsight\\Bin\\Tracing\\MDL\\release\\";   
      argv[2]="smoothBB.vtk";
      argv[3]="-spines.vtk";
	  argv[4]="smoothSpine.vtk";
	  argv[5]="RefineSkel.skel";
      //argv[9]="flag.txt";
      }
    else
      {
      cerr << "Usage: " << argv[0] << " <dir> <smooth backbone Skeleton file>"
           << "  <Spine skeleton file> <ExtraSpine file> <out refine skeleton file >"
           << " [measureTimeFlag]" << endl;
      exit(1);
      }
    }

  filedir = argv[1];
  infilename = filedir + argv[2];
  if(( smoothbackbone=fopen(infilename.c_str(), "rb")) == NULL)  // open smoothbackbone file
    {
    cerr << "couldn't open smoothbackbone file " << infilename << " for input" << endl;
    exit(-1);
    }

  tempfile1 = filedir + argv[3];
  cout << "second file name" << tempfile1 << endl;

  if((spine=fopen(tempfile1.c_str(), "rb")) == NULL)  // open spineskeleton file
    {
    cerr << "couldn't open spine skeleton file " << filedir << " for input" << endl;
    exit(-1);
    }

  tempfile2 = filedir + argv[4];
  cout << "second file name" << tempfile2 << endl;

  if((ExtraSpine=fopen(tempfile2.c_str(), "rb")) == NULL)  // open spineskeleton file
    {
    cerr << "couldn't open Extra spine skeleton file " << filedir << " for input" << endl;
    exit(-1);
    }
  
  tempfile3 = filedir + argv[5];
  cout << "second file name" << tempfile3 << endl;

  if((outrefineskel=fopen(tempfile3.c_str(), "wb")) == NULL)  // open spineskeleton file
    {
    cerr << "couldn't open  skeleton file " << filedir << " for output" << endl;
    exit(-1);
    }
   
 
   
  //------------------Read From smooth backbone skeleton VTK file ---------------------//

  // Skip first 4 lines of skeleton VTK file
  char str[200]; 
  for (i=0;i<12;i++)
    {
    fscanf(smoothbackbone,"%s",str);
    printf("%s\n", str); 
    }
  int NumBackbonePoints;
  fscanf(smoothbackbone,"%s",str);
  NumBackbonePoints = atoi(str);
  printf("There are %d skeleton' points \n", NumBackbonePoints);
  fscanf(smoothbackbone,"%s",str);

  AllBackbonepoints = new VoxelPosition[ NumBackbonePoints];
  float temp;  
  for (i=0;i<NumBackbonePoints;i++)
  {
   fscanf (smoothbackbone,"%f",&temp);
   AllBackbonepoints[i].x = temp;
   fscanf (smoothbackbone,"%f",&temp);
   AllBackbonepoints[i].y = temp; // change the sign,
   fscanf (smoothbackbone,"%f",&temp);
   AllBackbonepoints[i].z =temp;
  }
  fclose (smoothbackbone);

//------------------Read From spine skeleton VTK file ---------------------//

   // Skip first 4 lines of skeleton VTK file

  for (i=0;i<12;i++)
    {
    fscanf(spine,"%s",str);
    printf("%s\n", str); 
    }
  int NumSpinePoints;
  fscanf(spine,"%s",str);
  NumSpinePoints = atoi(str);
  printf("There are %d spine skeleton' points \n", NumSpinePoints);
  fscanf(spine,"%s",str);

  AllSpinepoints = new VoxelPosition[NumSpinePoints];
  float temp3;  
  for (i=0;i<NumSpinePoints;i++)
  {
   fscanf (spine,"%f",&temp3);
   AllSpinepoints[i].x = temp3;
   fscanf (spine,"%f",&temp3);
   AllSpinepoints[i].y = temp3; // change the sign,
   fscanf (spine,"%f",&temp3);
   AllSpinepoints[i].z =temp3;
  }
  fclose (spine);

//------------------Read From Extra-spine skeleton VTK file ---------------------//

 // Skip first 4 lines of skeleton VTK file
  //rewind (ExtraSpine);

  for (i=0;i<12;i++)
    {
    fscanf(ExtraSpine,"%s",str);
    printf("%s\n", str); 
    }
  int NumExtraSpinePoints;
  fscanf(ExtraSpine,"%s",str);
  NumExtraSpinePoints = atoi(str);
  
  fscanf(ExtraSpine,"%s",str);

  AllExtraSpinepoints = new VoxelPosition[NumExtraSpinePoints];
  //float temp;  
  for (i=0;i<NumExtraSpinePoints;i++)
  {
   fscanf (ExtraSpine,"%f",&temp);
   AllExtraSpinepoints[i].x = temp;
   fscanf (ExtraSpine,"%f",&temp);
   AllExtraSpinepoints[i].y = temp; // change the sign,
   fscanf (ExtraSpine,"%f",&temp);
   AllExtraSpinepoints[i].z =temp;
  }
  fclose (ExtraSpine);
 
 //---------------------------Output to Skel File
  for (i=0;i<NumBackbonePoints;i++)
    fprintf(outrefineskel,"%f %f %f %d\n", AllBackbonepoints[i].x, AllBackbonepoints[i].y , AllBackbonepoints[i].z,1);
  for (i=0;i<NumSpinePoints;i++)
    fprintf(outrefineskel,"%f %f %f %d\n", AllSpinepoints[i].x, AllSpinepoints[i].y , AllSpinepoints[i].z,1);
  for (i=0;i<NumExtraSpinePoints;i++)
    fprintf(outrefineskel,"%f %f %f %d\n", AllExtraSpinepoints[i].x, AllExtraSpinepoints[i].y , AllExtraSpinepoints[i].z,1);
  
  fclose (outrefineskel);
  //release memnory
  delete []AllBackbonepoints;
  delete []AllSpinepoints;
  delete []AllExtraSpinepoints;

 
}