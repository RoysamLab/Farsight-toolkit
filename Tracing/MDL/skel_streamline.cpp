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

// Form streamlines from saddle points and seed points
// --- Input: 1. 3D vector field
//            2. Seed points
// --- Output: skelelton file
// --- Author: Xiaosong Yuan
// --- Modify Date: 7/28/2007

//#include <io.h>
#if defined(_MSC_VER)
#pragma warning(disable : 4996)
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>

using namespace std;

#define FLOAT_SKELETON

#define MIN(x,y) (((x) < (y))?(x):(y))
#define MAX(x,y) (((x) > (y))?(x):(y))
#define SURF 100
#define INTERIOR 200
#define EPS 0.001

struct  VoxelPosition
{
	float x;
	float y;
	float z;
};

struct  Vector
{
	float xd;   // For large datasets, using double will cause memory shortage
	float yd;
	float zd;
};

int sign(float value) {
   if (value > 0) return 1;
   else if(value < 0) return -1;
        else return 0;
}


int sign1(float value) {
   if (value > 1e-5) return 1;
   else if(value < -1e-5) return -1;
        else return 0;
}

float veclength(Vector vecin) {
    return sqrt(vecin.xd * vecin.xd + vecin.yd * vecin.yd +vecin.zd * vecin.zd);
}

double dotProduct(Vector v1, Vector v2) {
    double dotProd;
    dotProd = v1.xd * v2.xd;
    dotProd += v1.yd * v2.yd;
    dotProd += v1.zd * v2.zd;
    return dotProd;
 //   return sign1(dotProd);  // test
}


/* Normalize a vector */
Vector normalize(Vector r)
{
    double mag;
    mag = r.xd * r.xd + r.yd * r.yd + r.zd * r.zd;
    if (mag != 0.0) {
        mag = 1.0 / sqrt(mag);
        r.xd *= mag;
        r.yd *= mag;
        r.zd *= mag;
        }
    return r;
}


Vector interpolation(float x, float y, float z, int sizx, int sizy, int sizz, Vector *forcevec);
void rk2(float x, float y, float z, int sizx, int sizy, int sizz, double steps, Vector *Force_ini, VoxelPosition *nextPos);



int main(int argc, char *argv[])
{
  ifstream fin, fseeds;
  FILE *fout;
  char *infilename, *seedsfilename;
  int sizeX,sizeY,sizeZ;         // Sizes in x,y,z dimensions
  int MidX, MidY, MidZ;
  int i,j,k;
  int *fc;
  Vector *force;
  
  long idx, iidx, slsz, sz;
  int measureTime = 0;
  //int numBound = 0;
  //int flagBound;
  //Vector pointForce, totalForce;
  VoxelPosition Startpos, Nextpos;
  //VoxelPosition seeds[80000];    // 90000 crash
  VoxelPosition *seeds;
  int idxSeeds;
  Vector vecin;
  //long iidx1, iidx2;
  int cc;
  int ii,jj,kk;
  int *FlagOnSkeleton;
  double vecLength;
  int NumCritPoints=0;
  float Ngrid;
  Vector OutForce;
  int streamSteps=0;
  //int FlagCloseToSkeleton;
  int x, y, z;
  double totalVecLength;
  double div, divx, divy, divz;
  float vectorMagnitude;

  if (argc < 8)
  {
    printf("Usage: %s <vector file> <xs> <ys> <zs> <vector mag> <seeds file> <out skel> [measureTimeFlag].\n",argv[0]);
    exit(1);
  }

  infilename = new char[80];
  infilename = argv[1];
  fin.open(infilename);
  if (!fin)  {
     cerr << "couldn't open " << infilename << " for input" << endl;
     return -1;
  }

  sizeX = atoi(argv[2]);
  sizeY = atoi(argv[3]);
  sizeZ = atoi(argv[4]);
  vectorMagnitude = atof(argv[5]);
  MidX = (sizeX-1)/2;
  MidY = (sizeY-1)/2;
  MidZ = (sizeZ-1)/2;

  seedsfilename = new char[80];
  seedsfilename = argv[6];
  fseeds.open(seedsfilename);
  if (!fseeds)  {
     cerr << "couldn't open " << seedsfilename << " for input" << endl;
     return -1;
  }

  if ((fout = fopen(argv[7],"w")) == NULL)
  {
    printf("Cannot open %s for writing\n",argv[6]);
    exit(1);
  }


  if (argc > 8)
    measureTime = 1;

  fc = new int[sizeX*sizeY*sizeZ];
  force = new Vector[sizeX*sizeY*sizeZ];
  FlagOnSkeleton = new int[sizeX*sizeY*sizeZ];

  slsz = sizeX*sizeY;		// slice size
  sz = slsz*sizeZ;
  seeds = new VoxelPosition[2000000];   //  VoxelPosition[200000]   
                                         // 15800000  when crit thre = 0.72
                                         // 35800000  when crit thre = 0.88 

  for(idx=0; idx<sz; idx++)   {  //Initialize to zeros
	  fc[idx]=0;
	  force[idx].xd = 0;
	  force[idx].yd = 0;
	  force[idx].zd = 0;

	  FlagOnSkeleton[idx] = 0;
  }

  // read in force vectors
  fin >> x >> y >> z >> vecin.xd >> vecin.yd >> vecin.zd;
  while (!fin.eof() ) 
  {
	    idx = z*slsz + y*sizeX +x;
	    // normalize the vectors
	    vecLength = sqrt(vecin.xd * vecin.xd + vecin.yd * vecin.yd +vecin.zd * vecin.zd);

		/* this implemention is time consuming, comments by xiao liang
	    if (vecLength ==0)  vecLength=0.000001;
	    force[idx].xd = vecin.xd / pow(vecLength, 5.0/6);
	    force[idx].yd = vecin.yd / pow(vecLength, 5.0/6);
	    force[idx].zd = vecin.zd / pow(vecLength, 5.0/6);
        */

		// by xiao liang
	    force[idx].xd = vecin.xd / (vecLength+0.0001);
	    force[idx].yd = vecin.yd / (vecLength+0.0001);
	    force[idx].zd = vecin.zd / (vecLength+0.0001);

	    //force[idx].xd = vecin.xd / vecLength;
	    //force[idx].yd = vecin.yd / vecLength;
	    //force[idx].zd = vecin.zd / vecLength;

	    fc[idx]=2;

		fin >> x >> y >> z >> vecin.xd >> vecin.yd >> vecin.zd;
	} //end while

  // make fc[]=2 of critical points

  //int seedsNumber=0;

  // make surface point fc[]=1
  for (k = 1; k < sizeZ-1; k++)
     for (j = 1; j < sizeY-1; j++)
        for (i = 1; i < sizeX-1; i++) 
		{
			idx = k*slsz + j*sizeX +i;
			if(fc[idx]==2) 
			{
				cc = 0;
				for (kk=-1; kk<=1; kk++)
					for (jj=-1; jj<=1; jj++)
						for (ii=-1; ii<=1; ii++) 
						{
							iidx = (k+kk)*slsz + (j+jj)*sizeX +(i+ii);
							if (fc[iidx]==0) cc++;
					    }//end for
				if (cc>=1) 
				{fc[idx] = 1;
				 //seedsNumber++;  // this is added by xiaoliang
				}
			}//end outer if 
	}// end for

   int numSeeds = 9058;
   //int numSeeds =seedsNumber/48;  // 6x4x2;
   //printf("Hello, Xiao liang,Seeds number is %d", numSeeds);
  //int numSeeds = 6690;
  //int numSeeds = 2521;

      // 
      // 9688   for "ThyCell1_T20DkR2.643x596x50.Aniso_k800d02t2.grdt.vec_curv0.2.seed"
      //   52   for "Ecadherin2_T60DkR30.1392x1024x9.Aniso_k800d02t2.grdt.vec_curv60.seed"
      //  283   for "Ecadherin2_T60DkR30.1392x1024x9.Aniso_k800d02t2.grdt.vec_curv20.seed"
      // 6712   for "Ecadherin2_T60DkR30.1392x1024x9.Aniso_k800d02t2.grdt.vec_curv3.seed"
      //11840   for "angioB1_T70DkR10.288x320x199.Aniso_k800d02t2.grdt.vec_curv0.2.seed"
      // 3398   for "Ecadherin_T60DkR30.1392x1024x9.Aniso_k800d02t2.grdt.vec_curv3.seed"
      //12028   for "Ecadherin_T60DkR30.1392x1024x9.Aniso_k800d02t2.grdt.vec_curv0.2.seed"
      //15194   for "Cal100m_T2DkR100.416x256x35.Aniso_k800d02t2.grdt.vec_curv0.seed"
      // 7303   for "Cal100m_T2DkR100.416x256x35.Aniso_k800d02t2.grdt.vec_curv4.seed"
      // 1557   for "Ania3sp1_T7DkR100.800x600x71.Aniso_k800d02t2.grdt.vec_curv100.seed"
      // 3162   for "MBFsp6_T7DkR100.512x512x68.Aniso_k800d02t2.grdt.vec_curv1.seed"
      //18337   for "MBFsp6_T7DkR100.512x512x68.Aniso_k800d02t2.grdt.vec_curv0.2.seed"]
      // 9506   for "MBFsp5_T2DkR100.308x512x49.Aniso_k800d02t2.grdt.vec_curv0.seed"
      // 5989   for "MBFsp5_T2DkR100.308x512x49.Aniso_k800d02t2.grdt.vec_curv0.2.seed"
      // 1347   for "MBFsp1_T2DkR100.1032x628x28.Aniso_k800d02t2.grdt.vec_curv0.9.seed"    

      // 783    for"Trach11D_T2DkR100.512x512x12.Aniso_k800d02t2.grdt.vec_curv0.2.seed"
      //  67    for"Trach11D_T10DkR1k.512x512x12.grdt.vec_curv10.seed" 
      //1036    for "Trach11C_T2DkR100.512x512x14.Aniso_k800d02t2.grdt.vec_curv0.2.seed"
      // 106    for"Trach11C_T10DkR1k.512x512x14.grdt.vec_curv10.seed" 
      // 1363   for"Trach11B_T2DkR100.512x512x14.Aniso_k800d02t2.grdt.vec_curv0.2.seed"
      //  149   for"Trach11B_T10DkR1k.512x512x14.grdt.vec_curv10.seed" 
      // 2521   for"Trach11A_T2DkR100.512x512x18.Aniso_k800d02t2.grdt.vec_curv0.seed"
      // 2239   for"Trach11A_T2DkR100.512x512x18.Aniso_k800d02t2.grdt.vec_curv0.06.seed"
      // 2077   for"Trach11A_T2DkR100.512x512x18.Aniso_k800d02t2.grdt.vec_curv0.1.seed"
      // 1312   for"Trach11A_T2DkR100.512x512x18.Aniso_k800d02t2.grdt.vec_curv0.2.seed"
      // 343 (212)    for"Trach11A_T10DkR1k.512x512x18.grdt.vec_curv2(10).seed" 
      // 167 (233)    for"Trach14E_T10DkR1k.512x512x18.grdt.vec_curv10.seed"   
      // 201 (228)    for"Trach14D_T10DkR1k.512x512x18.grdt.vec_curv10.seed"  
      // 149 (329)    for"Trach14C_T10DkR1k.512x512x17.grdt.vec_curv10.seed"  
      // 216 (445)    for"Trach14B_T10DkR1k.512x512x19.grdt.vec_curv10.seed" 
      // 229 (502)    for"Trach14A_T10DkR1k.512x512x25.grdt.vec_curv10.seed"   (502)is for un-deconvolved data
      //1468    for "Trach8E_T2DkR100.512x512x19.Aniso_k800d02t2.grdt.vec_curv0.2.seed"
      // 113    for "Trach8E_T10DkR1k.512x512x19.grdt.vec_curv10.seed"
      //1585    for "Trach8D_T2DkR100.512x512x19.Aniso_k800d02t2.grdt.vec_curv0.2.seed"
      // 122    for "Trach8D_T10DkR1k.512x512x19.grdt.vec_curv10.seed"
      // 342    for "Trach8C_T10DkR1k.512x512x22.grdt.vec_curv10.seed"
      //  811   for "                                                     curv0.4.see"
      // 1131   for "Trach8B_T2DkR100.512x512x18.Aniso_k800d02t2.grdt.vec_curv0.2.seed"
      // 108    for "Trach8B_T10DkR1k.512x512x18.grdt.vec_curv10.seed"
      // 421 (322)    for "Trach8A_T10DkR1k.512x512x23.grdt.vec_curv2(10).seed"
      // 468    for "Trach8A_T2DkR100.512x512x23.Aniso_k800d02t2.grdt.vec_curv6.seed"
      // 448    for "Trach8A_T2DkR100.512x512x23.Aniso_k800d02t2.grdt.vec_curv12.seed"
      // 643    for "Trach8A_T2DkR100.512x512x23.Aniso_k800d02t2.grdt.vec_curv2.seed"
      // 966    for "Trach8A_T2DkR100.512x512x23.Aniso_k800d02t2.grdt.vec_curv1.seed"
      // 457    for "Trach7E_T2DkR100.512x512x13.Aniso_k800d02t2.grdt.vec_curv0.2.seed"
      //  65    for "Trach7E_T10DkR1k.512x512x13.grdt.vec_curv10.seed"
	  // 160    for "Trach7D_T10DkR1k.512x512x17.grdt.vec_curv10.seed"
      // 410    for "Trach7C_T2DkR100.512x512x13.Aniso_k800d02t2.grdt.vec_curv0.2.seed"
      //  43    for "Trach7C_T10DkR1k.512x512x13.grdt.vec_curv10.seed"
      //  836   for "Trach7B_T2DkR100.512x512x13.Aniso_k800d02t2.grdt.vec_curv0.2.seed"
      //  44    for "Trach7B_T10DkR1k.512x512x13.grdt.vec_curv10.seed"
      //  91 (76)    for "Trach7A_T10DkR1k.512x512x14.grdt.vec_curv2(10).seed"
      // 967    for "Trach7A_T2DkR100.512x512x14.Aniso_k800d02t2.grdt.vec_curv-1.seed"
      // 764    for "Trach7A_T2DkR100.512x512x14.Aniso_k800d02t2.grdt.vec_curv0.seed"
      // 675    for "Trach7A_T2DkR100.512x512x14.Aniso_k800d02t2.grdt.vec_curv0.1.seed"
      // 233    for "Trach6E_T10DkR1k.512x512x27.grdt.vec_curv10.seed"
      // 408    for "Trach6D_T10DkR1k.512x512x28.grdt.vec_curv10.seed"
      // 404    for "Trach6C_T10DkR1k.512x512x36.grdt.vec_curv10.seed"
      // 317    for "Trach6B_T10DkR1k.512x512x30.grdt.vec_curv10.seed"
      // 275    for "Trach6A_T10DkR1k.512x512x26.grdt.vec_curv7.seed"
      // 528    for "Trach6A_T10DkR1k.512x512x26.grdt.vec_curv1.seed"
      // 250*   for "Trach6A_T10DkR1k.512x512x26.grdt.vec_curv3.seed"   * based on deconvolved data
      // 189 (251)    for "Trach6A_T10DkR1k.512x512x26.grdt.vec_curv10.seed"
      // 234    for "Trach6A_T10DkR1k.512x512x26.grdt.vec_curv15.seed"
      // 223    for "Trach6A_T10DkR1k.512x512x26.grdt.vec_curv20.seed"
      // 186    for "Trach6A_T10DkR1k.512x512x26.grdt.vec_curv50.seed"
      //  57    for "phantomSpines.64x40x31.grdt.vec_curv.0.seed"

      //19107   for "tlapse330_T2DkR10.512x512x35.Aniso_k800d02t2.grdt.vec_curv0.2.seed"     -- 2/13/2008
      //11992   for "tlapse330_T2DkR10.512x512x35.Aniso_k800d02t2.grdt.vec_curv1.seed" 
      //43504   for "tlapse330_T2DkR1k.512x512x35.Aniso_k800d02t4.grdt.vec_curv.0.seed"
      //61375   for "tlapse330_T2DkR1k.512x512x35.Aniso_k800d02t4.grdt.vec_curv-1.seed"
      //15169   for "tlapse330_T2DkR1k.512x512x35.Aniso_k800d02t4.grdt.vec_curv1.seed"
      //29749   for "tlapse330_T2DkR1k.512x512x35.Aniso_k800d02t4.grdt.vec_curv0.2.seed"

      //12299   for "tlapse330_T2DkR1k.512x512x35.grdt.vec_curv.0.seed"
      // 5835   for "Cal20m_T2DkR1k.416x256x35.Aniso_k800d02t4.raw.grdt.vec_curv.8.seed"
      // 2532   for "jt72D1_T3DkFR1kDk.512x512x45.Aniso_k800d02t8.raw.3-2-gradt_curv1.seed"
      // 3659   for "jt72D1_T3DkFR1kDk.512x512x45.grdt.vec_curv3.0.seed"
      // 1771   for "jt72D1_T3DkFR1kDk.512x512x45.grdt.vec_curv7.0.seed"

            // 28 for knight Hierarch #1
			// 59 for knight Hierarch #2
			// 49 for knightNS Hierarch #1
			// 64 for knightNS Hierarch #2
			// 39 for cube19x29x41
			// 129 for monster
			// 52 for tetrahedral
			// 40 for mushroom
			// 39 for 767
			// 20 for human
			// 885 for TightSegColon
			// 141 for TightSegColon_Cut
			//
			// 18 for cow with curvature thresh 16
			// 34 for cow with curvature thresh 11
			// 84 for cow with curvature thresh 6
			//
			// 25 for shark
			//  9 for hammerhead
			// 35 for dogDilate
			// 37 for woodman
			// 16 for teapot
			// 31 for toilet
			// 45 for cannon
			//
			// 0  for kD Hierarchical without seeds
			// 80 for kD Hierarchical #1
			// 58 for kD Hierarchical #2
			// 41 for trout
			// 18 for bolt
			// 24 for obj
			// 443 for mug
			//


  int numBoundSeeds = numSeeds;
  for (i = 0; i< numSeeds; i++)  {
	fseeds >> seeds[i].x >> seeds[i].y >> seeds[i].z;
  }



  // find all critical points -- method 1: consider two sides
/*  for (k = 1; k < sizeZ-1; k++)
     for (j = 1; j < sizeY-1; j++)
        for (i = 1; i < sizeX-1; i++) {
	    idx = k*slsz + j*sizeX +i;
	    if (fc[idx] == 0 || fc[idx]==1) continue;
	    iidx1 = k*slsz + j*sizeX +(i-1);
	    iidx2 = k*slsz + j*sizeX +(i+1);
	    if(sign(force[iidx1].xd) == sign(force[iidx2].xd)) continue;
	    iidx1 = k*slsz + (j-1)*sizeX + i;
	    iidx2 = k*slsz + (j+1)*sizeX + i;
	    if(sign(force[iidx1].yd) == sign(force[iidx2].yd)) continue;
	    iidx1 = (k-1)*slsz + j*sizeX + i;
	    iidx2 = (k+1)*slsz + j*sizeX + i;
	    if(sign(force[iidx1].zd) == sign(force[iidx2].zd)) continue;
	    seeds[numSeeds].x= i;
	    seeds[numSeeds].y= j;
	    seeds[numSeeds].z= k;
	    numSeeds++;
	}
*/

  // find all critical points -- method 2: consider one side
/*  for (k = 1; k < sizeZ-1; k++)
     for (j = 1; j < sizeY-1; j++)
        for (i = 1; i < sizeX-1; i++) {
	    idx = k*slsz + j*sizeX +i;
	    if (fc[idx] == 0) continue;

	    iidx1 = k*slsz + j*sizeX +(i+1);
	    if( fc[iidx1]==0 || sign(force[idx].xd) == sign(force[iidx1].xd)) continue;

	    iidx1 = k*slsz + (j+1)*sizeX + i;
	    if( fc[iidx1]==0 || sign(force[idx].yd) == sign(force[iidx1].yd)) continue;

	    iidx1 = (k+1)*slsz + j*sizeX + i;
	    if( fc[iidx1]==0 || sign(force[idx].zd) == sign(force[iidx1].zd)) continue;

	    NumCritPoints++;

            for (kk=0; kk<=1; kk++)
	        for (jj=0; jj<=1; jj++)
		     for (ii=0; ii<=1; ii++) {
	    		seeds[numSeeds].x= i+ii;
	    		seeds[numSeeds].y= j+jj;
	    		seeds[numSeeds].z= k+kk;
	    		numSeeds++;
	    }
	}
*/

  // find all critical points -- method 3: consider diagonal point
/*  for (k = 1; k < sizeZ-1; k++)
     for (j = 1; j < sizeY-1; j++)
        for (i = 1; i < sizeX-1; i++) {
	    idx = k*slsz + j*sizeX +i;
	    if (fc[idx] == 0) continue;

	    iidx1 =  k*slsz + j*sizeX + i;
	    iidx2 = (k+1)*slsz + (j+1)*sizeX +(i+1);
	    if( fc[iidx1]==0 || fc[iidx2]==0) continue;
	    if(!(sign(force[iidx1].xd)!=sign(force[iidx2].xd) &&
	         sign(force[iidx1].yd)!=sign(force[iidx2].yd) &&
		 sign(force[iidx1].zd)!=sign(force[iidx2].zd))) continue;

	    iidx1 =  k*slsz + j*sizeX + (i+1);
	    iidx2 = (k+1)*slsz + (j+1)*sizeX + i;
	    if( fc[iidx1]==0 || fc[iidx2]==0) continue;
	    if(!(sign(force[iidx1].xd)!=sign(force[iidx2].xd) &&
	         sign(force[iidx1].yd)!=sign(force[iidx2].yd) &&
		 sign(force[iidx1].zd)!=sign(force[iidx2].zd))) continue;

	    iidx1 =  k*slsz + (j+1)*sizeX + i;
	    iidx2 = (k+1)*slsz + j*sizeX + (i+1);
	    if( fc[iidx1]==0 || fc[iidx2]==0) continue;
	    if(!(sign(force[iidx1].xd)!=sign(force[iidx2].xd) &&
	         sign(force[iidx1].yd)!=sign(force[iidx2].yd) &&
		 sign(force[iidx1].zd)!=sign(force[iidx2].zd))) continue;

	    iidx1 =  k*slsz + (j+1)*sizeX + (i+1);
	    iidx2 = (k+1)*slsz + j*sizeX + i;
	    if( fc[iidx1]==0 || fc[iidx2]==0) continue;
	    if(!(sign(force[iidx1].xd)!=sign(force[iidx2].xd) &&
	         sign(force[iidx1].yd)!=sign(force[iidx2].yd) &&
		 sign(force[iidx1].zd)!=sign(force[iidx2].zd))) continue;


	    NumCritPoints++;

            for (int kk=0; kk<=1; kk++)
	        for (int jj=0; jj<=1; jj++)
		     for (int ii=0; ii<=1; ii++) {
	    		seeds[numSeeds].x= i+ii;
	    		seeds[numSeeds].y= j+jj;
	    		seeds[numSeeds].z= k+kk;
	    		numSeeds++;
	    }
	}
*/

  // find all critical points -- method 4: use points with small length of vector
  Ngrid = 10;
  for (k = 1; k < sizeZ-1; k++)
  //for (k = 1; k < -1; k++)  //TEST: skip all critical pts
     for (j = 1; j < sizeY-1; j++)
        for (i = 1; i < sizeX-1; i++) {
	    idx = k*slsz + j*sizeX +i;
	    totalVecLength = 0;
		divx = 0;  divy = 0;  divz = 0;
		// to check whether the position is near to boundaries, compute divergence in x,y,z respectively
	    if (fc[k*slsz + j*sizeX  +i] == 0) 
			continue;	else  {
				totalVecLength += veclength(force[k*slsz + j*sizeX  +i]);
				divx -= force[k*slsz + j*sizeX  +i].xd;
				divy -= force[k*slsz + j*sizeX  +i].yd;
				divz -= force[k*slsz + j*sizeX  +i].zd;
			}
	    if (fc[k*slsz + j*sizeX  +(i+1)] == 0) 
			continue;   else  {
				totalVecLength += veclength(force[k*slsz + j*sizeX  +(i+1)]);
				divx += force[k*slsz + j*sizeX  +(i+1)].xd;
				divy -= force[k*slsz + j*sizeX  +(i+1)].yd;
				divz -= force[k*slsz + j*sizeX  +(i+1)].zd;
			}
	    if (fc[k*slsz +(j+1)*sizeX + i] == 0) 
			continue;   else  {
				totalVecLength += veclength(force[k*slsz + (j+1)*sizeX  +i]);
				divx -= force[k*slsz + (j+1)*sizeX  +i].xd;
				divy += force[k*slsz + (j+1)*sizeX  +i].yd;
				divz -= force[k*slsz + (j+1)*sizeX  +i].zd;
			}
	    if (fc[k*slsz +(j+1)*sizeX + (i+1)] == 0) 
			continue;   else  {
				totalVecLength += veclength(force[k*slsz + (j+1)*sizeX  +(i+1)]);
				divx += force[k*slsz + (j+1)*sizeX  +(i+1)].xd;
				divy += force[k*slsz + (j+1)*sizeX  +(i+1)].yd;
				divz -= force[k*slsz + (j+1)*sizeX  +(i+1)].zd;
			}
	    if (fc[(k+1)*slsz + j*sizeX  +i] == 0) 
			continue;   else  {
				totalVecLength += veclength(force[(k+1)*slsz + j*sizeX  +i]);
				divx -= force[(k+1)*slsz + j*sizeX  +i].xd;
				divy -= force[(k+1)*slsz + j*sizeX  +i].yd;
				divz += force[(k+1)*slsz + j*sizeX  +i].zd;
			}
	    if (fc[(k+1)*slsz + j*sizeX  +(i+1)] == 0) 
			continue;   else  {
				totalVecLength += veclength(force[(k+1)*slsz + j*sizeX  +(i+1)]);
				divx += force[(k+1)*slsz + j*sizeX  +(i+1)].xd;
				divy -= force[(k+1)*slsz + j*sizeX  +(i+1)].yd;
				divz += force[(k+1)*slsz + j*sizeX  +(i+1)].zd;
			}
	    if (fc[(k+1)*slsz +(j+1)*sizeX + i] == 0) 
			continue;   else  {
				totalVecLength += veclength(force[(k+1)*slsz + (j+1)*sizeX  +i]);
				divx -= force[(k+1)*slsz + (j+1)*sizeX  +i].xd;
				divy += force[(k+1)*slsz + (j+1)*sizeX  +i].yd;
				divz += force[(k+1)*slsz + (j+1)*sizeX  +i].zd;
			}
	    if (fc[(k+1)*slsz +(j+1)*sizeX + (i+1)] == 0) 
			continue;   else  {
				totalVecLength += veclength(force[(k+1)*slsz + (j+1)*sizeX  +(i+1)]);
				divx += force[(k+1)*slsz + (j+1)*sizeX  +(i+1)].xd;
				divy += force[(k+1)*slsz + (j+1)*sizeX  +(i+1)].yd;
				divz += force[(k+1)*slsz + (j+1)*sizeX  +(i+1)].zd;
			}


		if (totalVecLength < 4)   {  // skip the zero vector areas
			continue;
		}

		div = divx + divy + divz;
		if (div > -1) {    //skip if the cube possibly contain a repelling point (divergence is too big)
			              // 1;  (good for phantom) 
						  // 0;  (better) 
						  // -1; (not always good)
			//printf("div=%f (%d)  ", div, k);
			continue;
		}
		
	    for (kk = 0; kk<Ngrid; kk++)
		  for (jj = 0; jj<Ngrid; jj++)
		    for (ii = 0; ii<Ngrid; ii++) {
		        OutForce=interpolation(i+ii/Ngrid, j+jj/Ngrid, k+kk/Ngrid, sizeX, sizeY, sizeZ, force);
			    if(veclength(OutForce) < vectorMagnitude) {   //<0.15
			    //if(veclength(OutForce) < 0.08) {   //<0.15

					//Mont10EunCh5_T30DkR10.1024x1024x51.Aniso_k800d02t2.grdt.vec 0.1-750142 crt pts

					//Mont10EunCh3_T50DkR10.1024x1024x51.Aniso_k800d02t2.grdt.vec 0.1-185377 crt pts

					//ThyCell1_T15DkR2.643x596x50.Aniso_k800d02t2.grdt.vec    0.05 -  8731 crt pts
					//                                                        0.08 - 34327 crt pts

					//ThyCell1_T20DkR2.643x596x50.Aniso_k800d02t2.grdt.vec    0.12 - 77582 crt pts
					//                                                        0.06 - 73336 crt pts

					//Ecadherin2_T60DkR30.1392x1024x9.Aniso_k800d02t2.grdt.vec 0.15-102912 crt pts
					//                                                         0.12- 46146 crt pts
					//                                                         0.10- 30475 crt pts
					//                                                         0.08- 19308 crt pts
					 
					//                                                         0.15-52079 crt pts
					//                                                         0.1 -12384 crt pts
					// Ecadherin_T60DkR30.1392x1024x9.Aniso_k800d02t2.grdt.vec 0.05 -3774 crt pts

					//                                                       0.09 -  4905 crt pts
					// Cal100m_T2DkR100.416x256x35.Aniso_k800d02t2.grdt.vec  0.12 - 13734 crt pts

					// Ania3sp1_T7DkR100.800x600x71.Aniso_k800d02t2.grdt.vec 0.06 -170880 crt pts
					//                                                       0.006-  828 crt pts
					// MBFsp6_T2DkR100.512x512x68.Aniso_k800d02t2.grdt.vec   0.06 - 3083 crt pts
					// MBFsp5_T2DkR100.308x512x49.Aniso_k800d02t2.grdt.vec   0.04 - 1565 crt pts
					//                                                       0.05 - 2746 crt pts
					//                                                       0.06 - 4271 crt pts
					// MBFsp1_T2DkR100.1032x628x28.Aniso_k800d02t2.grdt.vec  0.06 -36604 crt pts
					//                                                       0.04 -17446 crt pts
					//                                                       0.02 - 5092 crt pts

					// Trach14A_T2DkR100.512x512x25.Aniso_k800d02t2.grdt.vec 0.14 -14479 crt pts 
					//                                                       0.1  - 6725 crt pts
					// Trach6B_T2DkR100.512x512x30.Aniso_k800d02t2.grdt.vec  0.14 - 3896 crt pts
					// Trach6A_T2DkR100.512x512x26.Aniso_k800d02t2.grdt.vec  0.14 - 5693 crt pts
					//                                                       0.13 - 4804 crt pts
					//                                                       0.12 - 4070 crt pts 
					// Trach6A_T4DkR100.512x512x26.Aniso_k800d02t2.grdt.vec  0.15 - 4737 crt pts
					// Trach6A_T4DkR200.512x512x26.Aniso_k800d02t2.grdt.vec  0.15 - 4312 crt pts
					// Trach6A_T4DkR1k.512x512x26.Aniso_k800d02t2.grdt.vec   0.15 - 3941 crt pts
					//                                                       0.18 - 6646 crt pts

					//"Trach14D_T10DkR1k.512x512x12.Aniso_k800d02t2.grdt.vec" 0.15 - 5061 critical pts 
					//"Trach14C_T10DkR1k.512x512x14.Aniso_k800d02t2.grdt.vec" 0.15 - 3903 critical pts 
					//"Trach11B_T10DkR1k.512x512x14.Aniso_k800d02t2.grdt.vec" 0.15 - 4721 critical pts 
					//"Trach11A_T10DkR1k.512x512x18.Aniso_k800d02t2.grdt.vec" 0.15 - 5869 critical pts 
					//"Trach14E_T10DkR1k.512x512x18.Aniso_k800d02t2.grdt.vec" 0.15 - 5483 critical pts 
					//"Trach14D_T10DkR1k.512x512x18.Aniso_k800d02t2.grdt.vec" 0.15 - 4121 critical pts 
					//"Trach14C_T10DkR1k.512x512x17.Aniso_k800d02t2.grdt.vec" 0.15 -10917 critical pts 
					//"Trach14B_T10DkR1k.512x512x19.Aniso_k800d02t2.grdt.vec" 0.15 -14057 critical pts
					//"Trach14A_T10DkR1k.512x512x25.Aniso_k800d02t2.grdt.vec" 0.15 -14891 critical pts
					//"Trach8E_T10DkR1k.512x512x19.Aniso_k800d02t2.grdt.vec"  0.15 - 2209 critical pts
					//"Trach8D_T10DkR1k.512x512x19.Aniso_k800d02t2.grdt.vec"  0.15 - 2644 critical pts
					//"Trach8C_T10DkR1k.512x512x22.Aniso_k800d02t2.grdt.vec"  0.15 - 5045 critical pts
					//"Trach8B_T10DkR1k.512x512x18.Aniso_k800d02t2.grdt.vec"  0.15 - 3570 critical pts
					//"Trach8A_T10DkR1k.512x512x23.Aniso_k800d02t2.grdt.vec"  0.15 - 3826 critical pts
					//"Trach7E_T10DkR1k.512x512x13.Aniso_k800d02t2.grdt.vec"  0.15 - 2031 critical pts
					//"Trach7D_T10DkR1k.512x512x17.Aniso_k800d02t2.grdt.vec"  0.15 - 2636 critical pts
					//"Trach7C_T10DkR1k.512x512x13.Aniso_k800d02t2.grdt.vec"  0.15 - 2835 critical pts
					//"Trach7B_T10DkR1k.512x512x13.Aniso_k800d02t2.grdt.vec"  0.15 - 2981 critical pts
					//"Trach7A_T10DkR1k.512x512x14.Aniso_k800d02t2.grdt.vec"  0.15 - 4751 critical pts
					//"Trach6E_T10DkR1k.512x512x27.Aniso_k800d02t2.grdt.vec"  0.15 - 4751 critical pts  
					//"Trach6D_T10DkR1k.512x512x28.Aniso_k800d02t2.grdt.vec"  0.15 - 5514 critical pts  
					//"Trach6B_T10DkR1k.512x512x30.Aniso_k800d02t2.grdt.vec"  0.15 - 5310 critical pts
					//"Trach6A_T2DkR1k.512x512x26.Aniso_k800d02t4.grdt.vec" 0.06 - 4619 critical pts (missing many cri pts)
					                                                     // 0.1  -13678 (too many branches)

					//"Trach6A_T10DkR1k.512x512x26.Aniso_k800d02t2.grdt.vec" 0.1 - 1789 critical pts (OK, miss some)
					                                                     // 0.15 - 5574 (good)
					                                                     // 0.18 - 9341 (same as 0.15)
					                                                     // 0.22 -17587 (same as 0.15) 

					//"phantomSpines.64x40x31.grdt.vec"   0.12 -  195 critical pts (missing two spines)
					                                   // 0.16 -  420 critical pts (still missing two spines)
													   // 0.2  -  797 critical pts (missing one spine)
					                                   // 0.25 - 1877 critical pts (missing one spine)
													   // 0.3  - 3277 critical pts (missing one spine)
													   // 0.4  - 7802 critical pts (missing one spine)
													   // 0.5  -17015 critical pts (missing one spine)
													   // 0.6  -34296 critical pts (some extra pts)
													   // 0.55 -23911 critical pts (relatively good)

					//"tlapse330_T2DkR1k.512x512x35.Aniso_k800d02t4.grdt.vec"   0.12 - 23177 critical pts
					                                                         // 0.2  - 85048 critical pts
												      //(After change div to 1) 0.12 - 25264 critical pts
																			 // 0.2  - 88704 critical pts

					// "tlapse330_T2DkR10.512x512x35.Aniso_k800d02t2.grdt.vec"  0.12 - 75474 critical pts
					                                                          //0.08 - 25295 critical pts


					//"Cal20m_T2DkR1k.416x256x35.Aniso_k800d02t4.raw.grdt.vec"  0.12 -  7160 critical pts (not connected)
																		     // 0.16 - 13593 critical pts (missing spines)
												

					// 0.02 for "jt72D1_thrs3_dilate_flf_rmCC1k_dilat.512x512x45.raw.vec"-  5634 critical pts
                    // 0.04 for "jt72D1_thrs3_dilate_flf_rmCC1k_dilat.512x512x45.raw.vec"- 12960 critical pts
					// 0.06 for "jt72D1_thrs3_dilate_flf_rmCC1k_dilat.512x512x45.raw.vec"- 23226 critical pts
					// 0.08 for "jt72D1_thrs3_dilate_flf_rmCC1k_dilat.512x512x45.raw.vec"- 36100 critical pts
					// 0.10 for "jt72D1_thrs3_dilate_flf_rmCC1k_dilat.512x512x45.raw.vec"- 51504 critical pts (good)

					// 0.10 for "jt72D1_thrs3_dilate_flf_rmCC1k_dilat.512x512x45.raw.Rg20_1-7Samp.vec"- 33317 critical pts
					// 0.04 for "jt72D1_thrs3_dilate_flf_rmCC1k_dilat.512x512x45.raw.Rg20_1-7Samp.vec"-  3827 critical pts
					// 0.06 for "jt72D1_thrs3_dilate_flf_rmCC1k_dilat.512x512x45.raw.Rg20_1-7Samp.vec"-  9888 critical pts
					// 0.08 for "jt72D1_thrs3_dilate_flf_rmCC1k_dilat.512x512x45.raw.Rg20_1-7Samp.vec"- 19515 critical pts (good)

					// 0.09 for "jt72D1_thrs3_dilate_flf_rmCC1k_dilat.512x512x45.raw.Rg20_1-3Samp.vec"- 40668 critical pts (good)

					// 0.09 for "jt72D1_thrs3_dilate_flf_rmCC1k_dilat.512x512x45.raw.Rg20_1-1Samp.vec"- 43326 critical pts


					// 0.08 for "jt72D1_thrs3_dilaK_flf_rmCC1k_dilaK.512x512x45.raw.all-1-gradt.vec" - 24180 critical pts
					// 0.04 for "jt72D1_thrs3_dilaK_flf_rmCC1k_dilaK.512x512x45.raw.all-1-gradt.vec" -  9779 critical pts

					// 0.08 for "jt72D1_thrs3_dilaK_flf_rmCC1k_dilaK.512x512x45.raw.3-2-gradt.vec"  -  26756 critical pts
					// 0.07 for "jt72D1_thrs3_dilaK_flf_rmCC1k_dilaK.512x512x45.raw.3-2-gradt.vec"  -  21546 critical pts

					// 0.08 for "jt72D1_T3DkFR1kDk.512x512x45.IsoD_d02t8.raw.3-2-gradt.vec"   -  4704 critical pts
					// 0.12 for "jt72D1_T3DkFR1kDk.512x512x45.IsoD_d02t8.raw.3-2-gradt.vec"   - 15668 critical pts (good)

					// 0.12 for "jt72D1_T3DkFR1kDk.512x512x45.Anis_k800d02t8.raw.3-2-gradt.vec" - 17088 critical pts
													   // remove zero area (totalVecLength < 4) - 16751 critical pts (similar)
													   // remove zero area (totalVecLength <4.1)- 16701 critical pts
													   // remove zero area (totalVecLength <4.3)- 16366 critical pts (similar)

													   // remove cubes that div > 5.0   - 15926 critical pts
													   // remove cubes that div > 4.0   - 15678 critical pts
													   // remove cubes that div > 3.0   - 15385 critical pts
													   // remove cubes that div > 2.0   - 15063 critical pts
													   // remove cubes that div > 1.0   - 14732 critical pts
													   // remove cubes that div > 0     - 14197 critical pts
													   // remove cubes that div > -1    - 13747 critical pts (good)
													   // remove cubes that div > -2    - 12959 critical pts (miss one spine)
													   // remove cubes that div > -4    - 11524 critical pts (miss many)

					// 0.12 for "jt72D1_T3DkFR1kDk.512x512x45.Anis_k800d02t16.raw.3-2-gradt.vec" 
													   // remove cubes that div > -1    - 22070 critical pts (too many crit pts)
					// 0.08 for "jt72D1_T3DkFR1kDk.512x512x45.Anis_k800d02t16.raw.3-2-gradt.vec" 
													   // remove cubes that div > -1    -  6669 critical pts (not good)

					// 0.12 for "jt72D1_T3DkFR1kDk.512x512x45.Anis_k800d02t8sm.raw.3-2-gradt.vec" 
													   // remove cubes that div > -1    - 45816 critical pts (too many crit pts)
					// 0.08 for "jt72D1_T3DkFR1kDk.512x512x45.Anis_k800d02t8sm.raw.3-2-gradt.vec" 
													   // remove cubes that div > -0.5   - 9492 critical pts 
					// 0.10 for "jt72D1_T3DkFR1kDk.512x512x45.Anis_k800d02t8sm.raw.3-2-gradt.vec" 
													   // remove cubes that div > 0    -  29566 critical pts (some skel miss)
					// 0.09 for "jt72D1_T3DkFR1kDk.512x512x45.Anis_k800d02t8sm.raw.3-2-gradt.vec" 
													   // remove cubes that div > 1    -  10510 critical pts (some skel miss)
													   // remove cubes that div > 10   -  11528 critical pts (some skel miss)


					// 0.12 for "jt72D1_T3DkFR1kDk.512x512x45.Anis_k80d02t8.raw.3-2-gradt.vec" - not good, crash
					// 0.04 for "jt72D1_T3DkFR1kDk.512x512x45.Anis_k80d02t8.raw.3-2-gradt.vec" - 5662 critical pts (not good)
					// 0.06 for "jt72D1_T3DkFR1kDk.512x512x45.Anis_k80d02t8.raw.3-2-gradt.vec" -14899 critical pts 
					    





					        //vector length threshold for critical points
							// 0.05   for knight
							// 0.05   for knightNS
							// 0.1    for TightSegColon_Cut
							// 0.1    for shark
							// 0.05   for hammerhead
							// 0.2    for dogDilate
							// 0.06   for woodman
							// 0.05   for teapot
							// 0.1    for toilet
							// 0.07   for cannon
							// 0.065  for kD
							//  0.1   for enforcer
							// 0.08   for plant
							// 0.002  for 767
							// 0.05   for torus
							// 0.1   for bolt
							// 0.1   for obj
							// 0.008 for mug
							// 0.08  for trout
							// 0.1  for SukBenchCut
							// 0.02 for achilleasCut



				seeds[numSeeds].x= i + ii/Ngrid;
				seeds[numSeeds].y= j + jj/Ngrid;
				seeds[numSeeds].z= k + kk/Ngrid;
				numSeeds++;
				NumCritPoints++;
			}
	    }
  }


 printf("Number of critical points is: %d, and number of seeds %d\n", NumCritPoints, numSeeds);


 //fprintf(fout,"%d %d %d %f %f\n", 1, 1, 1, -1.0, -1.0);

 // seeds[0..numBoundSeeds-1] = boundary seeds, seeds[numBoundSeeds.. numSeeds-1] = critical points
 idxSeeds = numSeeds-1;
 while (idxSeeds >= 0) {

	Startpos.x=seeds[idxSeeds].x; Startpos.y=seeds[idxSeeds].y; Startpos.z=seeds[idxSeeds].z;
	idx = (int)Startpos.z * slsz + (int)Startpos.y *sizeX + (int)Startpos.x;

	//Check whether the high curv points are within D-voxel distance to the existing skeleton
	int Dvoxel = 2;    // 3 is good for many results
	if (idxSeeds < numBoundSeeds)   Dvoxel= 2;  //4; // 10 <- Found too big on April 12, 2006, which cause less spines
	                                             // 3 is good for tlapse330, 6 is good for Trach
	int FlagWithin = 0;
	if (idxSeeds < numBoundSeeds) {
		for (kk = -Dvoxel+1; kk <= Dvoxel-1; kk++)
		   for (jj = -Dvoxel+1; jj <= Dvoxel-1; jj++)
		      for (ii = -Dvoxel+1; ii <= Dvoxel-1; ii++) {
				  iidx = idx + kk*slsz + jj*sizeX + ii;
				  if (iidx < 0 || iidx >= sz)  continue;
		          if(FlagOnSkeleton[iidx] == 1)  FlagWithin = 1;
		}
		if(FlagWithin == 1) {
			idxSeeds--;
			continue;
		}
	}

	// being able to not show critical point in the streamlines
	if (idxSeeds >= numBoundSeeds)
		{
		   #ifdef FLOAT_SKELETON
		     fprintf(fout,"%f %f %f %d\n", Startpos.x, Startpos.y, Startpos.z, 1);
		    #else
		     fprintf(fout,"%d %d %d %d %d\n",(int)Startpos.x,(int)Startpos.y,(int)Startpos.z, idxSeeds, idxSeeds);
		   #endif
		}
	    else {
		   #ifdef FLOAT_SKELETON
	             fprintf(fout,"%f %f %f %d\n", Startpos.x, Startpos.y, Startpos.z, 1);
		    #else
		     fprintf(fout,"%d %d %d %d %d\n",(int)Startpos.x,(int)Startpos.y,(int)Startpos.z, idxSeeds, idxSeeds);
		   #endif
	    }

	FlagOnSkeleton[idx] = 1;

	while(streamSteps < 4000)   {    // < 4000
		rk2(Startpos.x, Startpos.y, Startpos.z, sizeX, sizeY, sizeZ, 0.8, force, &Nextpos);   //0.2, 0.8, 2
		streamSteps++;
		Startpos.x = Nextpos.x;
		Startpos.y = Nextpos.y;
		Startpos.z = Nextpos.z;

		idx = (int)Nextpos.z *slsz + (int)Nextpos.y *sizeX + (int)Nextpos.x;
        if (FlagOnSkeleton[idx] != 1) {
			#ifdef FLOAT_SKELETON
	                  fprintf(fout,"%f %f %f %d\n", Nextpos.x, Nextpos.y, Nextpos.z, 1);
			//printf("%f %f %f %d\n", Nextpos.x, Nextpos.y, Nextpos.z, 1);
			 #else
			  fprintf(fout,"%d %d %d %d %d\n",(int)Nextpos.x,(int)Nextpos.y,(int)Nextpos.z, idxSeeds, idxSeeds);
			#endif
			FlagOnSkeleton[idx] = 1;
		}
	}

	streamSteps = 0;
	idxSeeds--;
 }


   fclose(fout);
   // release memory by xiao liang
   delete []fc;
   delete []force;
   delete []FlagOnSkeleton;
   delete []seeds;

   printf("End \n");
   return 0;

}




Vector interpolation(float x, float y, float z, int sizx, int sizy, int sizz, Vector *forcevec)
    {
	float alpha, beta, gamma;
	Vector forceInt;
	long slsz;

	alpha=x-int(x);
	beta=y-int(y);
	gamma=z-int(z);
	slsz=sizy*sizx;

	forceInt.xd=forcevec[int(z)*slsz + int(y)*sizx + int(x)].xd*(1-alpha)*(1-beta)*(1-gamma)
			+forcevec[(int(z)+1)*slsz + int(y)*sizx + int(x)].xd*(1-alpha)*(1-beta)*gamma
			+forcevec[int(z)*slsz + (int(y)+1)*sizx + int(x)].xd*(1-alpha)*beta*(1-gamma)
			+forcevec[int(z)*slsz + int(y)*sizx + (int(x)+1)].xd*alpha*(1-beta)*(1-gamma)
			+forcevec[(int(z)+1)*slsz + int(y)*sizx + (int(x)+1)].xd*alpha*(1-beta)*gamma
			+forcevec[int(z)*slsz + (int(y)+1)*sizx + (int(x)+1)].xd*alpha*beta*(1-gamma)
			+forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + int(x)].xd*(1-alpha)*beta*gamma
			+forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + (int(x)+1)].xd*(alpha*beta*gamma);

	forceInt.yd=forcevec[int(z)*slsz + int(y)*sizx + int(x)].yd*(1-alpha)*(1-beta)*(1-gamma)
			+forcevec[(int(z)+1)*slsz + int(y)*sizx + int(x)].yd*(1-alpha)*(1-beta)*gamma
			+forcevec[int(z)*slsz + (int(y)+1)*sizx + int(x)].yd*(1-alpha)*beta*(1-gamma)
			+forcevec[int(z)*slsz + int(y)*sizx + (int(x)+1)].yd*alpha*(1-beta)*(1-gamma)
			+forcevec[(int(z)+1)*slsz + int(y)*sizx + (int(x)+1)].yd*alpha*(1-beta)*gamma
			+forcevec[int(z)*slsz + (int(y)+1)*sizx + (int(x)+1)].yd*alpha*beta*(1-gamma)
			+forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + int(x)].yd*(1-alpha)*beta*gamma
			+forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + (int(x)+1)].yd*alpha*beta*gamma;

	forceInt.zd=forcevec[int(z)*slsz + int(y)*sizx + int(x)].zd*(1-alpha)*(1-beta)*(1-gamma)
			+forcevec[(int(z)+1)*slsz + int(y)*sizx + int(x)].zd*(1-alpha)*(1-beta)*gamma
			+forcevec[int(z)*slsz + (int(y)+1)*sizx + int(x)].zd*(1-alpha)*beta*(1-gamma)
			+forcevec[int(z)*slsz + int(y)*sizx + (int(x)+1)].zd*alpha*(1-beta)*(1-gamma)
			+forcevec[(int(z)+1)*slsz + int(y)*sizx + (int(x)+1)].zd*alpha*(1-beta)*gamma
			+forcevec[int(z)*slsz + (int(y)+1)*sizx + (int(x)+1)].zd*alpha*beta*(1-gamma)
			+forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + int(x)].zd*(1-alpha)*beta*gamma
			+forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + (int(x)+1)].zd*alpha*beta*gamma;

	return(forceInt);
    }


void rk2(float x, float y, float z, int sizx, int sizy, int sizz, double steps, Vector *Force_ini, VoxelPosition *nextPos)
   {
	long slsz;
	Vector OutForce;
	//float x1, y1, z1;
	slsz=sizy*sizx;

	OutForce=interpolation(x,y,z,sizx,sizy,sizz,Force_ini);

	x = x + OutForce.xd * steps;
	y = y + OutForce.yd * steps;
	z = z + OutForce.zd * steps;

	nextPos->x = x;
	nextPos->y = y;
	nextPos->z = z; 
   }




