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

// B-spline fitting code  
// Input format: .raw  (3D points)
// Input format: .vtk  (3D skeleton file)
// Author: Liang XIAO, RPI, According to the Xiaosong's Matlab Code
// Date: 12/Aug/2009


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
//#include "BSpline.h"

using namespace std;
//using namespace boost;

#define Dimention  3;

struct  VoxelPosition
{
  float x;
  float y;
  float z;
}; 

struct Graphprop{
   int deg;
   int outVert[8];
};

vnl_vector<float> Q[3];

//-------------------------some B-Spline function-----------------------------//
float CubicNPSpline(float t, int n)
{
  float value;
  float s = t - n;
  if (s >= -2 && s < -1)
    {
    value = (pow((2+s),3))/6;
    }
  else if (s >= -1 && s < 0)
    {
    value = (4 - 6*pow(s,2) - 3*pow(s,3))/6;
    }
  else  if (s >= 0 && s < 1)
    {
    value = (4 - 6*pow(s,2) + 3*pow(s,3))/6;
    }
  else  if (s >= 1 && s <= 2)  
    {
    value = (pow((2-s),3))/6;
    }
  else
    {
    value = 0;
    }
  return value;
}

void FitNPSpline(int Nb,int datalength, VoxelPosition *datapoints)
{
  int Ns = datalength;
  vnl_vector<float> knots(Ns);
  knots(0) = 0;
  float step = static_cast<float>(Nb)/static_cast<float>(Ns-1);
  for(int i=1; i<Ns; i++)
    {
    knots(i) = knots(i-1) + step;
    }

  // Construct basis matrix
  vnl_matrix<float> Basis(Ns,Nb+5);
  for(int i=0; i<Ns; i++)
    {
    for(int j = -3; j< Nb+2; j++)
      {
      Basis(i,j+3) = CubicNPSpline(knots(i),j+1);
      }
    }
  vnl_matrix<float> Binv = vnl_matrix_inverse<float>(Basis);
  vnl_vector<float> x(Ns); 
  vnl_vector<float> y(Ns);
  vnl_vector<float> z(Ns);
  for(int ij=0; ij<Ns; ij++)
    {
    x(ij) = datapoints[ij].x;
    y(ij) = datapoints[ij].y;
    z(ij) = datapoints[ij].z;
    }
  vnl_vector<float> Qx = Binv * x;
  vnl_vector<float> Qy = Binv * y;
  vnl_vector<float> Qz = Binv * z;
  Q[0] = Qx;
  Q[1] = Qy;
  Q[2] = Qz;
}

void SampleNPSpline(int NPointsSample,VoxelPosition *PSample, float P1, float P2)
{
  vnl_vector<float> Qx = Q[0];
  vnl_vector<float> Qy = Q[1];
  vnl_vector<float> Qz = Q[2];
  int Nb = Qx.size() - 5;
  
  // sample points with Qx and Qy
  vnl_vector<float> t(NPointsSample);
  t(0) = P1 * Nb;
  float steps = static_cast<float>((P2-P1)*Nb)/static_cast<float>(NPointsSample-1);
  for(int i=1; i<NPointsSample; i++)
    {
    t(i) = t(i-1) + steps;
    }
  vnl_vector<float> xs(NPointsSample); 
  vnl_vector<float> ys(NPointsSample); 
  vnl_vector<float> zs(NPointsSample); 
  xs.fill(0);
  ys.fill(0);
  zs.fill(0);
  for(int j=-3; j<Nb+2; j++)
    {
    for(int i=0; i<NPointsSample; i++)
      {
      xs(i) = xs(i) + Qx(j+3) * CubicNPSpline(t(i), j+1);
      ys(i) = ys(i) + Qy(j+3) * CubicNPSpline(t(i), j+1);
      zs(i) = zs(i) + Qz(j+3) * CubicNPSpline(t(i), j+1);
      }
    }

  for(int k =0; k<NPointsSample; k++)
    {
    PSample[k].x=xs(k);
    PSample[k].y=ys(k);
    PSample[k].z=zs(k);
    }
  //return PSample;
}

#define DATATYPEIN unsigned char
//------------------------------some sub function-----------------------------//

int round (float number);
int selffloor(float number);

void spap2(int PolyNum, int order, int *t_val, VoxelPosition *points,
           VoxelPosition * coef, double *knots);
double dist2pts(double x1,double y1,double z1,double x2,double y2, double z2);

int main(int argc, char **argv)
{
  //ifstream fin;
  FILE *infile, *inskeleton;
  FILE *outbackbone, *outExspine;
  long sz;  // the variable: sz= sizex*sizey
  
  //DATATYPEIN *volin, *volvessel, *somaDist;
  std::string filedir;
  std::string infilename;
  std::string tempfile; //  temp file, 

  VoxelPosition *Allpoints;
  //VoxelPosition *posExtraSpines;
  //VoxelPosition *pointsVTK;

  DATATYPEIN  * volin; 
    
  int sizeX,sizeY,sizeZ;         // Sizes in x,y,z dimensions
  int i,j,k;
  
  if (argc < 7)
    {
    if(argc > 1 && strcmp(argv[1], "debug") == 0)
      {
      //---- for debug test----//
      argv[0] = (char *)"BSplineFitting.exe";
      argv[1] = (char *)"Trach11A.512x512x18.raw";
      argv[2] = (char *)"Backbone.vtk";
      argv[3] = (char *)"512";
      argv[4] = (char *)"512";
      argv[5] = (char *)"18";
      argv[6] = (char *)"smoothBB.vtk";
      argv[7] = (char *)"Exspine.vtk";
      //argv[9]="flag.txt";
      }
    else
      {
      cerr << "Usage: " << argv[0] << " <3Ddata raw file> <Skeleton file>"
           << " <xs> <ys> <zs> <out backbonevtk graph> <out spinevtk graph>"
           << " [measureTimeFlag]" << endl;
      exit(1);
      }
    }

  infilename = argv[1];
  if(( infile=fopen(infilename.c_str(), "rb")) == NULL)  // open volume file
    {
    cerr << "couldn't open volume file " << infilename << " for input" << endl;
    exit(-1);
    }

  tempfile = argv[2];
  cout << "second file name" << tempfile << endl;

  if((inskeleton=fopen(tempfile.c_str(), "rb")) == NULL)  // open skeleton file
    {
    cerr << "couldn't open skeleton file " << filedir << " for input" << endl;
    exit(-1);
    }
  
  sizeX = atoi(argv[3]); // read the size parameters of input image 
  sizeY = atoi(argv[4]);
  sizeZ = atoi(argv[5]);
  sz= sizeX*sizeY;
  
  //----------------------------read volume data-----------------------------//

  volin = (DATATYPEIN*)malloc(sizeX*sizeY*sizeZ*sizeof(DATATYPEIN)); //  the memory application for volume data

  if ( fread(volin, sizeof(DATATYPEIN), sizeX*sizeY*sizeZ, infile) < (unsigned int)(sizeX*sizeY*sizeZ))  // read in vol file
  {
    printf("File size is not the same as volume size\n");
    exit(1);
  }

  //------------------Processing with skeleton VTK file ---------------------//

  // Skip first 4 lines of skeleton VTK file
  char str[200]; 
  for (i=0;i<12;i++)
    {
    fscanf(inskeleton,"%s",str);
    printf("%s\n", str); 
    }
  int NumAllPoints;
  fscanf(inskeleton,"%s",str);
  NumAllPoints = atoi(str);
  printf("There are %d skeleton' points \n", NumAllPoints);
  fscanf(inskeleton,"%s",str);
  //--------------------------------------------------------------------------//
  Allpoints = new VoxelPosition[NumAllPoints];

  //ifstream fin;
  //fin = inskeleton;
  float temp;
  //------------------------read 3D skeleton points---------------------------//
  for (i=0;i<NumAllPoints;i++)
  {
   fscanf (inskeleton,"%f",&temp);
   Allpoints[i].x = temp;
   fscanf (inskeleton,"%f",&temp);
     Allpoints[i].y = -temp; // change the sign,
   fscanf (inskeleton,"%f",&temp);
   Allpoints[i].z =temp;
  }
  //--------------------------------------------------------------------------//
  float *flagOnBackbone;

  Graphprop * graphInfPoints;
  graphInfPoints = new Graphprop[NumAllPoints];

  flagOnBackbone = new float[NumAllPoints];

  //------Initialize the two arrays: flagOnBackbone (for backbone index)------//

  //----------------graphInfPoints(for graph connection)----------------------//

  for (i=0;i<NumAllPoints;i++)
    {
    flagOnBackbone[i]=0;
    graphInfPoints[i].deg=0;
    for (j=0;j<8;j++)
      {
      graphInfPoints[i].outVert[j]=0;
      }
    } 
  //--------------------------------------------------------------------------//
  //---------------------- read in the 'LINES' line --------------------------//
  int NumLines, SkipNmuber;
  fscanf(inskeleton,"%s",str); // read LINES
  fscanf(inskeleton,"%d",&NumLines);
  fscanf(inskeleton,"%d",&SkipNmuber); // skip the last number in this line
   
  //------------------------------ sub - end ---------------------------------//

  //Method 2: Assign index to backbone according to graph branch connection//
  int tmpdeg;
  int tmp1, tmp2,tmp3;
  int mod;

  for (i=0;i<NumLines;i++)
    {
    // Read one line into an array;
    fscanf (inskeleton,"%d",&tmp1); // skip first number in the Line
    fscanf (inskeleton,"%d",&tmp2); // read second number in the line 
    fscanf (inskeleton,"%d",&tmp3); // read third number in the line 
  
    tmpdeg = graphInfPoints[tmp2].deg; // get the corresponding deg of the node;
    graphInfPoints[tmp2].deg = tmpdeg+1;
    mod = tmpdeg % 8;
    graphInfPoints[tmp2].outVert[mod] = tmp3;

    tmpdeg = graphInfPoints[tmp3].deg;
    graphInfPoints[tmp3].deg = tmpdeg+1;
    mod = tmpdeg % 8;
    graphInfPoints[tmp3].outVert[mod] = tmp2;
    } 

  //------------------------------- sub end ----------------------------------//
  fclose(inskeleton); 

  //-----Get rid of junction point so that breaking up backbone segments------//
  int outpt;
  for (i=0;i<NumAllPoints;i++)
    {
    if (graphInfPoints[i].deg >=3)
      {
      for (j=0;j<graphInfPoints[i].deg;j++)
        {
        outpt = graphInfPoints[i].outVert[j];
        graphInfPoints[outpt].deg = graphInfPoints[outpt].deg - 1;
        if ( graphInfPoints[outpt].deg ==1)
          {
          if (graphInfPoints[outpt].outVert[0] == i)
            {
            graphInfPoints[outpt].outVert[0] = graphInfPoints[outpt].outVert[1]; 
            // Remove junction point from neighbor points
            }
          }
        } // end for
      }
    }// end for

  //---------------------------------sub end ---------------------------------//

  //--------------In Method 2: Doing assigning index to backones--------------//
  // Note: flagOnBackbone(i) is used to save in format xx.ooo where xx is index
  //of backbone, ooo is order of pts.
  int NumBranches = 0;
  float indexBBpts =1;
  //%NumPoly = [0, 0, 0, 0, 0];
  int NumPoly[200];// = zeros(1, 200);

  int lastpoint, nextpoint; 

  for (k=0;k<200;k++)
    {
    NumPoly[k] = 0;  
    }

  for (i=0;i<NumAllPoints;i++)
    {
    if (graphInfPoints[i].deg == 1)
      {
      NumBranches = NumBranches + 1;
      indexBBpts = 1;
      flagOnBackbone[i] = NumBranches + indexBBpts/1000;
      lastpoint = i;
      nextpoint = graphInfPoints[i].outVert[0]; // the first is [0]
      // if not end of branch, keep going.
      while (graphInfPoints[nextpoint].deg == 2) 
        {
        indexBBpts = indexBBpts + 1;
        flagOnBackbone[nextpoint] = NumBranches + indexBBpts/1000;
        NumPoly[NumBranches] = NumPoly[NumBranches] + 1;
        // outVert[0]  --- the first
        if (graphInfPoints[nextpoint].outVert[0] == lastpoint)
          {
          lastpoint = nextpoint;
          // outVert[1] ---- the second
          nextpoint = graphInfPoints[nextpoint].outVert[1];
          }// end if 
        else 
          {
          lastpoint = nextpoint;
          // outVert[0]  --- the first
          nextpoint = graphInfPoints[nextpoint].outVert[0]; 
          } // end else
        }// end while 
      if (graphInfPoints[nextpoint].deg == 1)
        {
        graphInfPoints[nextpoint].deg = 0;
        indexBBpts = indexBBpts + 1;
        flagOnBackbone[nextpoint] = NumBranches+ indexBBpts/1000;
        } // end if
      }// end if
    }// enf for 

  //--------------------NumPoly = round(NumPoly * 0.05)-----------------------//
  // for (i = 0; i< NumAllPoints; i++) printf("%f  ",flagOnBackbone[i]);

 /* FILE * flagfile;
  tempfile = filedir + argv[9];
  cout << "flagbone file name" << tempfile << endl;

  if(( flagfile=fopen(tempfile.c_str(), "w"))==NULL)  // open skeleton file
    {
    cerr << "couldn't open flag file " << filedir << " for input" << endl;
    exit(-1);
    }
  */
  //for (i = 0; i< NumAllPoints; i++) printf("%f  ",flagOnBackbone[i]);
  //for (i = 0; i< NumAllPoints; i++)
    //{
    //fprintf (flagfile,"%f \n",flagOnBackbone[i]);
    //}
  //flagOnBackbone[i]=temp1;}

  //---------------------NumPoly = round(NumPoly * 0.05)----------------------//

  for (k=0;k<200;k++)
    {
    NumPoly[k] =  round(float ((NumPoly[k])*0.05)); 
    //0.04 , how to set this parameter ?    float 
    }

  //---------------------------------sub end----------------------------------//

  int numVTKpts = 0;
  int numExtraSpines = 0;

  //-----------Begin to curve fitting for each branch of backbones------------//
 
  VoxelPosition *points;
  points = NULL;
  int tmp=0;
  int NumPoints; // belong to one branch;
  int *numBranchpts; // need to further
  numBranchpts = new int[NumBranches+1]; // this is to record the points in each branches;
  float discretePt=1; // sampling rate
  int NewNumAllPoints = (int) (((double) NumAllPoints)/discretePt)+1; // the number points which to be write in VTK file ;
  cout << "the number of total new points  " <<  NewNumAllPoints <<  endl;
  VoxelPosition *pointsVTK; 
  VoxelPosition *posExtraSpines;
  pointsVTK = new VoxelPosition[NewNumAllPoints];
  posExtraSpines = new VoxelPosition[NewNumAllPoints]; 

  for (i=0;i<NumBranches+1;i++)
      numBranchpts[i]=0;

  for (i=0;i<NewNumAllPoints-1;i++) //initialization 
    {
    pointsVTK[i].x =0;
    pointsVTK[i].y =0;
    pointsVTK[i].z =0;
	posExtraSpines[i].x=-9999;
	posExtraSpines[i].y=-9999;
	posExtraSpines[i].z=-9999;
    }
  pointsVTK[NewNumAllPoints-1].x =-9999;
  pointsVTK[NewNumAllPoints-1].y =-9999;
  pointsVTK[NewNumAllPoints-1].z =-9999; //end flag

  cout << "There are " <<  NumBranches  << " NumBranches" << endl; 

  int tmpindex; // 
  int x_val_index;
  float minDist;
  float dis;
  int minIndex;
  float mintemp;
  int *Index_Min_dist2spline;
  float  *Min_dist2spline;
  int IsLocalMax;
  int x1,y1,z1;
  int aveInt;
  for (k=1;k<=NumBranches;k++)
    {
    //----------------------Note there will a huge loop-----------------------//
    if (points != NULL)
      {
      delete [] points;
      }
    points = new VoxelPosition[NumAllPoints];
    tmp=0;
    for (i = 0; i< NumAllPoints; i++) // i=0 means the first number 
      {
      if (selffloor(flagOnBackbone[i])== k)
        {
        //indexBBpts =
          //(int) ((flagOnBackbone[i] - selffloor(flagOnBackbone[i])) * 1000 ); 
        tmpindex =
          round((flagOnBackbone[i] - selffloor(flagOnBackbone[i])) * 1000)-1;
            
        //cout << "tmpindex" << tmpindex << endl;
        //tmpindex = (int) indexBBpts;
        // In order to exchange x-coord and y-coord
        points[tmpindex].x = Allpoints[i].x;  
        points[tmpindex].y = Allpoints[i].y; 
        points[tmpindex].z = Allpoints[i].z;
        //printf("%f %f %f\n", points[tmpindex].x, points[tmpindex].y,
                 //points[tmpindex].z);
        tmp++;
        }   
      }// end for 

    numBranchpts[k] = tmp;
    NumPoints = tmp; 
    //------------------------------------------------------------------------//
    //for (i=0;i<NumPoints; i++)
      //{
      // for test 
      //printf("%f %f %f\n", points[i].x , points[i].y , points[i].z); 
      //}
    // cout << "come to Fit here \n" << endl;
    FitNPSpline(6,NumPoints, points);// B-Spline coeff;  // NB =4
   
    //--------------the very important function: Spline sampling--------------//
    
    NumPoints = (int) (NumPoints / discretePt); // new sampling points;

    VoxelPosition *xyz_vals;
    xyz_vals = new VoxelPosition [NumPoints];
     
    //-------------------------Call the resampling----------------------------//
    //cout << "come to Sample here \n" <<endl; 
    SampleNPSpline(NumPoints,xyz_vals, 0, 1);

    //for (i=0;i<NumPoints;i++)
    //  {
    //  printf("%f %f %f\n",xyz_vals[i].x,xyz_vals[i].y,xyz_vals[i].z);
    //  }
    // copy these points which locate at this branch to pointsVTK;
    // currently there are numVTKpts VTK pts to be write 

    // copy these points which locate at this branch to pointsVTK;
    for (i=numVTKpts;i<numVTKpts+NumPoints;i++)
      {
      pointsVTK[i].x =xyz_vals[i-numVTKpts].x;
      pointsVTK[i].y =xyz_vals[i-numVTKpts].y;
      pointsVTK[i].z =xyz_vals[i-numVTKpts].z;
      
      //printf("%f %f %f\n",pointsVTK[i].x,pointsVTK[i].y,pointsVTK[i].z);
      }
  

  
    //------------------------------- Compute missing spines from spline function and original 3D points-------------------------------//
       mintemp=points[0].x;
	   for (i = 0;i<NumPoints-1;i++)
	   {	 
		   if (points[i].x<=mintemp)
			   mintemp = points[i].x;
	   }
        

        Min_dist2spline = new float [NumPoints];
        Index_Min_dist2spline = new int [NumPoints];

	    for (i = 0;i<NumPoints-1;i++)
		{
		  x_val_index = selffloor((points[i].x-mintemp)/discretePt);  //
          if(x_val_index <= 0)       
			  x_val_index=0;
          if (x_val_index >= NumPoints-1)
			  x_val_index=NumPoints-1; 

        //Initialize the min dist and corresponding index with a parallel point on the spline function
		  minDist = (float) dist2pts(points[i].y, points[i].x, points[i].z,  xyz_vals[x_val_index].y,  xyz_vals[x_val_index].x, xyz_vals[x_val_index].z);
          minIndex = x_val_index;
          for( j = 0;j<NumPoints-1;j++)
		  {  if(abs(points[i].x - xyz_vals[j].x)< minDist)  // If the point on spline falls within the local range
			  {//evaluate the current distance of two points
                dis = (float) dist2pts(points[i].y, points[i].x, points[i].z, xyz_vals[j].y,  xyz_vals[j].x, xyz_vals[j].z);
                if (dis < minDist)
				{
                    minDist = dis;
                    minIndex = j;
				}
			  }
		  }
        Min_dist2spline[i] = minDist;
        Index_Min_dist2spline[i] = minIndex;
        // plot3([points(1,i), x_val(minIndex)]+DispBias,  [points(2,i), yz_vals(1,minIndex)]+DispBias,   [points(3,i), yz_vals(2,minIndex)]+DispBias, 'r:'); 
        
		}// end for 

	 //----------------------------------- Find all the spines with local max distance---------------------------------------------------//  
	    for (i = 0;i<NumPoints-1;i++)
		{
         IsLocalMax = 1;
         for(j = 0;j<NumPoints-1;j++)
		 {
			 if (dist2pts(points[j].y, points[j].x, points[j].z,  points[i].y, points[i].x, points[i].z ) < 5)
                //Consider only neighbor points
                if (Min_dist2spline[j] > Min_dist2spline[i])  
				{IsLocalMax=0;  break; }
		 }// end subfor

		 x1 = round(points[i].x);  
		 y1 = round(points[i].y);
		 z1 = round(points[i].z);
         long idx = z1*sz + y1*sizeX +x1;
         aveInt = volin[idx];//[z1][y1][x1];//D(z, y, x); 
		 idx=Index_Min_dist2spline[i];
		 x1 = round(xyz_vals[idx].x);
         y1 = round(xyz_vals[idx].y);
		 z1 = round(xyz_vals[idx].z);

		 idx = z1*sz + y1*sizeX +x1;

		 aveInt = (aveInt + volin[idx])/2; 

        // printf("%d \n", aveInt);
		
      
		 if ( IsLocalMax == 1 && Min_dist2spline[i]>= 1.3 && aveInt > 20)   // 2, 1.5 - threshold for the detect of missing spines  
		 {
			
            numExtraSpines = numExtraSpines + 1;
			
			posExtraSpines[numExtraSpines*2-2].x = points[i].x;
            posExtraSpines[numExtraSpines*2-2].y = points[i].y;
		    posExtraSpines[numExtraSpines*2-2].z = points[i].z;
        
			idx = Index_Min_dist2spline[i];
			posExtraSpines[numExtraSpines*2-1].x =  xyz_vals[idx].x;
            posExtraSpines[numExtraSpines*2-1].y =  xyz_vals[idx].y;
			posExtraSpines[numExtraSpines*2-1].z =  xyz_vals[idx].z;
			 
		 }
        
		} 

	 numVTKpts += NumPoints;
     delete []xyz_vals; // release memory
	 delete []Index_Min_dist2spline;
     delete []Min_dist2spline;
     cout << "there are " << numVTKpts << "points on branch " << k << endl;
    } // End of NumBranches of curve fitting //   end of loop

  if (points != NULL)
    {
    delete [] points;
    }

  cout << "the loop is over" << endl; 
 
  //-----------------------construct two files to write-----------------------//
  if ((outbackbone = fopen(argv[6], "w")) == NULL)
    {
    printf("Cannot open  backbone file %s for writing\n",argv[6]);
    exit(1);
    }

 if ((outExspine = fopen(argv[7], "w")) == NULL)
    {
    printf("Cannot open  backbone file %s for writing\n",argv[7]);
    exit(1);
    }
  //--------------------------------------------------------------------------//
  //------------------Output the extra spines to a vtk file-------------------//
   fprintf(outExspine, "# vtk DataFile Version 3.0\n");
  fprintf(outExspine,"MST of skel\n");
  fprintf(outExspine,"ASCII\n");
  fprintf(outExspine,"DATASET POLYDATA\n");
  fprintf(outExspine,"POINTS %d float\n", 2 * numExtraSpines);

  for (i = 0; i<numExtraSpines*2;i++)
	  fprintf(outExspine,"%f %f %f\n",posExtraSpines[i].x,-posExtraSpines[i].y,posExtraSpines[i].z);
   
  fprintf(outExspine, "LINES %d %d\n", numExtraSpines, numExtraSpines*3);

  for (i = 0;i<numExtraSpines;i++)
    fprintf(outExspine, "2 %d %d\n", 2*i, 2*i+1); // 
  fclose(outExspine);

  //-------------------------end write to outExspine--------------------------//

  //-----------------Output the smooth Backbone to a vtk file-----------------//
  fprintf(outbackbone, "# vtk DataFile Version 3.0\n");
  fprintf(outbackbone,"MST of skel\n");
  fprintf(outbackbone,"ASCII\n");
  fprintf(outbackbone,"DATASET POLYDATA\n");
  fprintf(outbackbone,"POINTS %d float\n", numVTKpts);

  cout << "come to write VTK file" << endl;

  for (i = 0; i<numVTKpts;i++)
    {
    // - pointsVTK[i].y  is just for view
    fprintf(outbackbone,"%f %f %f\n",pointsVTK[i].x, -pointsVTK[i].y,
            pointsVTK[i].z);
    //printf("%f %f %f\n",pointsVTK[i].x,pointsVTK[i].y,pointsVTK[i].z);
    }
   
  fprintf(outbackbone, "LINES %d %d\n", numVTKpts-NumBranches, (numVTKpts-NumBranches)*3);

  int indexPts = 0;
  cout << "I am come to here, one step I will be succeed!";
  for (i = 1; i<= NumBranches; i++)
  //for (k=1;k<=1;k++)
  //i=2;
    {
    for(j =0; j< numBranchpts[i]-1;j++)
      {
      fprintf(outbackbone, "2 %d %d\n", j+indexPts, j+indexPts+1);
      }
    indexPts = indexPts + (int)numBranchpts[i];
    }

  cout << "I am come to here ";
  fclose(outbackbone);

  //--------------------------release global memory---------------------------//
  delete [] pointsVTK;
  delete [] Allpoints;
  delete [] graphInfPoints;
  delete [] numBranchpts;
  return 1;  
}  // end 

//---------------------Round to nearest integer function----------------------//
int round (float number)
{  
  int temp;
  temp = int (number);
  if ((number-temp)>0.5) 
    {
    return (int(number+float(0.5)));
    }
  else 
    {
    return (temp);
    }
}

//-------------------------------selffloor end--------------------------------//
// y = floor(a) rounds fi object a to the nearest integer in the direction of
//negative infinity and returns the result in fi object y.
int selffloor(float number)
{
  if (number >0 || number ==0) 
    {
    return (int)number;
    }
  else
    {
    return (int)(number-0.5); 
    }
}

//--------------------h = dist2pts(x1, y1, z1, x2, y2, z2)--------------------//
double dist2pts(double x1,double y1,double z1,double x2,double y2, double z2)
{
  double h;
  h=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
  return h;
}

