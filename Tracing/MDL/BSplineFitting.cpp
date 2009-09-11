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
// Input format: .vtk  (3D sketon file)
// Author: Liang XIAO, RPI, According to the Xiaosong's Matlab Code
// Date: 12/Aug/2009


#if defined(_MSC_VER)
#pragma warning(disable : 4996)
#endif

#include <iostream>
#include <fstream>
//#include <io.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;
//using namespace boost;

struct  VoxelPosition
{
	float x;
	float y;
	float z;
};

struct Graphprop{
	 float deg;
	 float outVert[8];
};

struct BSPLINE
{
     char  B;
     double *knots; //: [1x17 double]
     VoxelPosition *coefs; //: [3x13 double]
     int number; //: 13
     int order;  // 4
     int dim;    // 3
};

#define DATATYPEIN unsigned char

//int main(int argc, char *argv[])

//----------------------------------------sub function ----------------------------//

int round (float number);
int selffloor(float number);


int main(void)
{
   
  char *argv[9];
  int argc=9; // test
  //ifstream fin;
  FILE *infile, *insketon;
  FILE *outbackbone, *outExspine;
  long sz;  // the variable: sz= sizx*sizey
  

  //DATATYPEIN *volin, *volvessel, *somaDist;
  char *filedir;
  char *infilename;
  char *tempfile; //  temp file, 

  VoxelPosition *Allpoints;


  DATATYPEIN  * volin; 

	  
  int sizeX,sizeY,sizeZ;         // Sizes in x,y,z dimensions
  int i,j,k;
  
  //---- for debug test----//
  argv[0]="BSplineFitting.exe";
  argv[1]="D:\\MDL0804\\B-Spline\\Debug\\";
  argv[2]="Trach11A.512x512x18.raw";
  argv[3]="Backbone.vtk";
  argv[4]="512";
  argv[5]="512";
  argv[6]="18";
  argv[7]="smoothBB.vtk";
  argv[8]="Exspine.vtk";
  if (argc < 9)
  {
    printf("Usage: %s <dir> < 3Ddata raw file> <Sketon file> <xs> <ys> <zs>  <out backbonevtk graph> <out spinevtk graph>[measureTimeFlag].\n", argv[0]);
    exit(1);
  }
   

  filedir = new char[400];
  infilename = new char[300];
  tempfile = new char [300];
  filedir = argv[1];//
  
  strcpy(infilename, filedir);// read path
  
  strcat(infilename, argv[2]);// path + filename
  
   if(( infile=fopen(infilename,"rb"))==NULL)  // open volume file
  { cerr << "couldn't open volume file " << infilename << " for input" << endl;
	exit(-1);
  }

  strcpy(tempfile, filedir);
  strcat(tempfile, argv[3]);

  cout << "second file name" << tempfile << endl;

  if(( insketon=fopen(tempfile,"rb"))==NULL)  // open sketon file
  { cerr << "couldn't open sketon file " << filedir << " for input" << endl;
	exit(-1);
  }
  
  sizeX = atoi(argv[4]); // read the size parameters of input image 
  sizeY = atoi(argv[5]);
  sizeZ = atoi(argv[6]);
  sz= sizeX*sizeY;


  
 //------------------------------------ read volume data ----------------------------------//

  
  volin = (DATATYPEIN*)malloc(sizeX*sizeY*sizeZ*sizeof(DATATYPEIN)); //  the memory application for volume data

  if ( fread(volin, sizeof(DATATYPEIN), sz, infile) < (unsigned long)sz)  // read in vol file
  {
    printf("File size is not the same as volume size\n");
    exit(1);
  }

  // ------------------------------------------------------------------------------------------//

  
  //---------------------------  construct two files to write -----------------------------------//
  
  if ((outbackbone = fopen(argv[7], "w")) == NULL)
  {
    printf("Cannot open  backbone file %s for writing\n",argv[7]);
    exit(1);
  }

 if ((outExspine = fopen(argv[8], "w")) == NULL)
  {
    printf("Cannot open  backbone file %s for writing\n",argv[8]);
    exit(1);
  }


 //--------------------------------Processing with sketon VTK file -------------------------------------//
 // Skip first 4 lines of sketon VTK file
 
 char str[200]; 
 for (i=0;i<12;i++)
 {fscanf(insketon,"%s",str);
  printf("%s\n", str); 
 }
 int NumAllPoints;
 fscanf(insketon,"%s",str);
 NumAllPoints = atoi(str);
 printf("There are %d sketon' points \n", NumAllPoints);
 fscanf(insketon,"%s",str);
 //---------------------------------------------------------------------------------------------------//
 Allpoints = new VoxelPosition[NumAllPoints];

 //ifstream fin;
 //fin = insketon;
  float temp;
 //------------------------------- read 3D sketon points ---------------------------------------------//
 for (i=0;i<NumAllPoints;i++)
 {
	 fscanf (insketon,"%f",&temp);
	 Allpoints[i].x = temp;
	 fscanf (insketon,"%f",&temp);
     Allpoints[i].y = -temp; // change the sign,
	 fscanf (insketon,"%f",&temp);
	 Allpoints[i].z =temp;
 }
 // --------------------------------------------------------------------------------------------------//
 float *flagOnBackbone;


 Graphprop * graphInfPoints;
 graphInfPoints = new Graphprop[NumAllPoints];

 flagOnBackbone = new float[NumAllPoints];
 
 //--------------Initialize the two arrays: flagOnBackbone (for backbone index); --------------------//

 //----------------  graphInfPoints(for graph connection) -------------------------------------------//
 
 for (i=0;i<NumAllPoints;i++)
 {
	 flagOnBackbone[i]=0;
     graphInfPoints[i].deg=0;
	 for (j=0;j<8;j++)
	  graphInfPoints[i].outVert[j]=0;
 } 
//---------------------------------------------------------------------------------------------------//

//---------------------------------- read in the 'LINES' line ---------------------------------------//
int NumLines, SkipNmuber;
fscanf(insketon,"%s",str); // read LINES
fscanf(insketon,"%d",&NumLines);
fscanf(insketon,"%d",&SkipNmuber); // skip the last number in this line
 
//---------------------------------- sub - end -----------------------------------------------------//

// --------------Method 2: Assign index to backbone according to graph branch connection-----------//
int tmpdeg;
int tmp1, tmp2,tmp3;
int mod;

 for (i=0;i<NumLines;i++)
 {
	 // Read one line into an array;

	 fscanf (insketon,"%d",&tmp1); // skip first number in the Line
	 fscanf (insketon,"%d",&tmp2); // read second number in the line 
	 fscanf (insketon,"%d",&tmp3); // read third number in the line 
	 
	
	 tmpdeg = graphInfPoints[tmp2].deg; // get the corresponding deg of the node;
     graphInfPoints[tmp2].deg = tmpdeg+1;
     mod = tmpdeg % 8;
     graphInfPoints[tmp2].outVert[mod] = tmp3;

	 tmpdeg = graphInfPoints[tmp3].deg;
     graphInfPoints[tmp3].deg = tmpdeg+1;
     mod = tmpdeg % 8;
     graphInfPoints[tmp3].outVert[mod] = tmp2;
 } 

 //-------------------------------------- sub end ---------------------------------------//
 fclose(insketon); 


// --------- Get rid of junction point so that breaking up backbone segments -----------//

int outpt;
for (i=0;i<NumAllPoints;i++)
{
  if (graphInfPoints[i].deg >=3)
	  for (j=0;j<graphInfPoints[i].deg;j++)
	  {
            outpt = graphInfPoints[i].outVert[j];
            graphInfPoints[outpt].deg = graphInfPoints[outpt].deg - 1;
			if ( graphInfPoints[outpt].deg ==1)
				if (graphInfPoints[outpt].outVert[0] == i)
					graphInfPoints[outpt].outVert[0] = graphInfPoints[outpt].outVert[1]; 
				  // Remove junction point from neighbor points
      } // end for
}// end for

// ---------- -------------------------------sub end ---------------------------------//


// ---------- ------------In Method 2: Doing assigning index to backones--------------//
// Note: flagOnBackbone(i) is used to save in format xx.ooo where xx is index of backbone, ooo is order of pts.
int NumBranches = 0;
float indexBBpts =1;
//%NumPoly = [0, 0, 0, 0, 0];
int NumPoly[200];// = zeros(1, 200);

int lastpoint, nextpoint; 

for (k=0;k<200;k++)
     NumPoly[k] = 0;  

for (i=0;i<NumAllPoints;i++)
{
  if (graphInfPoints[i].deg == 1)
  {
    NumBranches = NumBranches + 1;
    indexBBpts = 1;
    flagOnBackbone[i] = NumBranches + indexBBpts/1000;
	lastpoint = i;
    nextpoint = graphInfPoints[i].outVert[0]; // the first is [0]
	while (graphInfPoints[nextpoint].deg == 2) // if not end of branch, keep going.
	{
           indexBBpts = indexBBpts + 1;
           flagOnBackbone[nextpoint] = NumBranches + indexBBpts/1000;
           NumPoly[NumBranches] = NumPoly[NumBranches] + 1;
		   if (graphInfPoints[nextpoint].outVert[0] == lastpoint) // outVert[0]  --- the first
		   {
                lastpoint = nextpoint;
                nextpoint = graphInfPoints[nextpoint].outVert[1]; //  outVert[1] ---- the second
		   }// end if 
		   else 
		   {
                lastpoint = nextpoint;
                nextpoint = graphInfPoints[nextpoint].outVert[0]; // outVert[0]  --- the first
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

//-- -------------------  NumPoly = round(NumPoly * 0.05) ------------------------------ -------//

for (k=0;k<200;k++)
     NumPoly[k] = round(NumPoly[k]*0.05); // 0.04 , how to set this parameter ? 

// ---------------- ---------------sub end ----------------------------------------------------//

int numVTKpts = 0;
int numExtraSpines = 0;

// -------------------------------- Begin to curve fitting for each branch of backbones------//
VoxelPosition *points;
points = NULL;
// points = new VoxelPosition[NumAllPoints];
//

int tmp=0;
int NumPoints;

float *t_val;

float discretePt; // sampling rate

for (k=1;k<NumBranches;k++)
{

	//----------------------------Note there will a huge loop --------------------------------//
    if (points != NULL) delete []points;
	for (i = 0; i< NumAllPoints; i++) // i=0 means the first number 
	{
     
		if (selffloor(flagOnBackbone[i])== k)
		{
			indexBBpts = ((flagOnBackbone[i] - selffloor(flagOnBackbone[i])) * 1000 );
            points = new VoxelPosition[NumAllPoints];
		    points[i].x = Allpoints[i].x;  // In order to exchange x-coord and y-coord
			points[i].y = Allpoints[i].y; 
            points[i].z = Allpoints[i].z; 
			tmp++;
		}   
            // %plot3(points(1,:), points(2,:), points(3,:), 'b:', 'LineWidth', 2); hold on;
            //%In order to exchange x-coord and y-coord
            // %plot3(points(2,:)+DispBias, points(1,:)+DispBias,  points(3,:)+DispBias, 'b*', 'MarkerSize', 2); hold on;
            // plot3(Allpoints(1,i)+DispBias, Allpoints(2,i)+DispBias, Allpoints(3,i)+DispBias, 'b*', 'MarkerSize', 1); hold on;
	}// end for 

    // NumPoints = size(points, 2);
     NumPoints = tmp;
	 t_val = new float[NumPoints];
	 for (j=0;j<NumPoints;j++)
		 t_val[j]=j+1;

     
	 // sp = spap2(NumPoly[k], 4, t_val, points)
	 
	 delete []t_val;
	 //    ----------------the very important function --------Spline fitting -----------------------//
     discretePt = 0.5;
	 t_val = new float[NumPoints *2];
	 for (j=0;j<NumPoints;j++)
	 {t_val[2*j]   = j*2;
	  t_val[2*j+1] = 2*j+0.5;
	 } //  this must be modified 



     

} // End of NumBranches of curve fitting 

// ints = size(points, 2);



  
 return 1;  
}


// ---------------------------- Round to nearest integer function-------------------------------//

int round (float number)
{  
   int temp;
   temp = int (number);
   if ((number-temp)>0.5) 
     return (int(number+0.5));
   else 
     return (temp);
}

//---------------------------------------------selffloor  end -----------------------------------------//
// y = floor(a) rounds fi object a to the nearest integer in the direction of negative infinity and returns the result in fi object y.//

int selffloor(float number)
{
   if (number >0 || number ==0) 
	   return (int)number;
   else
	   return (int)(number-0.5); 
}

// -------------------------------------------B-Spline fitting function ------------------------------//
 //sp = spap2(NumPoly[k], 4, t_val, points)

BSPLINE spap2(int PolyNum, int order, int *t_val, VoxelPosition *points)
{
	// ---------------------
	
	int coefslength;
	int knotslength;
	BSPLINE sp;
	
	sp.B='B';
	sp.coefs = new VoxelPosition [coefslength];
	sp.dim=  3;
	sp.knots = new double [knotslength];
	sp.number =13; // .....
	sp.order =4;

	// ----------------------------------------------//
	{
      cout << " do B-Spline fitting ----/n ";
		//..........................
	}
    return sp;
}
