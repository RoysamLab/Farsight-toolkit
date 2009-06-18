// ----
// ----  Compute the vector field of any volume objects
// ----
// ----  Modification by : Xiaosong Yuan, RPI
// ----  Date: Oct. 5, 2005
// ----  Input : Binary 3D volume with sizes.
// ----  Output: ASCII file with vector 3 components for all object voxels
// ----

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MIN(x,y) (((x) < (y))?(x):(y))
#define MAX(x,y) (((x) > (y))?(x):(y))

#define DATATYPEIN unsigned char

#define SURF 100
#define INTERIOR 200

#define SurfSampling 1
#define OUTSIDE_FIELD 0   // compute also the vector field outside the objects

#define ForceRange 60   //20

struct  VoxelPosition
{
	short x;
	short y;
	short z;
};

struct  Vector
{
	double xd;   // For large datasets, using double will cause memory shortage
	double yd;
	double zd;
};

Vector ptForce(float x1, float y1, float z1, float x2, float y2, float z2) {
    double r;
    Vector pointForce;
    r = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));
    pointForce.xd = (x2-x1)/pow(r, 7);
//printf("r=%f, x2-x1=%f,y2-y1=%f,z2-z1=%f,xd=%f ", r,x2-x1,y2-y1,z2-z1,pointForce.xd);
    pointForce.yd = (y2-y1)/pow(r, 7);
    pointForce.zd = (z2-z1)/pow(r, 7);
    return pointForce;
}


void createLUTable (Vector *vecLUTable) {
	int x, y, z;
	int idx;
	double r;
	for (z = 0; z< ForceRange; z++)
		for (y = 0; y< ForceRange; y++)
			for (x = 0; x< ForceRange; x++) {
				if (x==0 && y==0 && z==0) {
					vecLUTable[0].xd = 0;
					vecLUTable[0].yd = 0;
					vecLUTable[0].zd = 0;
					continue;
				}
                r = sqrt(double(x*x + y*y +z*z));
				idx = z*ForceRange*ForceRange + y*ForceRange + x;
				vecLUTable[idx].xd = x/pow(r,7);
				vecLUTable[idx].yd = y/pow(r,7);
				vecLUTable[idx].zd = z/pow(r,7);
			}
}

Vector ptForceLUTable (Vector *vecLUTable, int x1, int y1, int z1, int x2, int y2, int z2) {
	Vector pointForce;
	int idx;
	idx = abs(z1-z2)*ForceRange*ForceRange + abs(y1-y2)*ForceRange + abs(x1-x2);
	pointForce = vecLUTable[idx];
	if (x2 < x1)  pointForce.xd *= -1;
	if (y2 < y1)  pointForce.yd *= -1;
	if (z2 < z1)  pointForce.zd *= -1;
	return pointForce;  
}


int sign(float value) {
   if (value > 1e-5) return 1;
   else if(value < -1e-5) return -1;
        else return 0;
}

double dotProduct(Vector v1, Vector v2) {
    double dotProd;
    dotProd = v1.xd * v2.xd;
    dotProd += v1.yd * v2.yd;
    dotProd += v1.zd * v2.zd;
//    return dotProd;
    return sign(dotProd);  // test
}


int main(int argc, char *argv[])
{
	
  FILE *filein;
  FILE *fileout;
  DATATYPEIN *volin;
  VoxelPosition *boundVoxs;
  Vector *force;
  Vector *vecLUTable;

  int sizeX,sizeY,sizeZ;         // Sizes in x,y,z dimensions
  int MidX, MidY, MidZ;
  int i,j,k,n,s,ss;
  long idx, iidx, sls, sz;
  int measureTime = 0;
  int numBound = 0;
  int flagBound;
  Vector pointForce, totalForce;
  int border;

  if (argc < 6)
  {
    printf("Usage: %s <volfile> <xs> <ys> <zs> <outfile> [measureTimeFlag].\n",argv[0]);
    exit(1);
  }


  if ((filein = fopen(argv[1],"rb")) == NULL)
  {
    printf("Cannot open %s\n",argv[1]);
    exit(1);
  }

  sizeX = atoi(argv[2]);
  sizeY = atoi(argv[3]);
  sizeZ = atoi(argv[4]);
  MidX = (sizeX-1)/2;
  MidY = (sizeY-1)/2;
  MidZ = (sizeZ-1)/2;
  //ThresDiv = atof(argv[6]);
  if (argc > 6)
     measureTime = 1;

  volin = (DATATYPEIN*)malloc(sizeX*sizeY*sizeZ*sizeof(DATATYPEIN));
  force = (Vector *)malloc(sizeX*sizeY*sizeZ*sizeof(Vector));
  vecLUTable = (Vector *)malloc(ForceRange*ForceRange*ForceRange*sizeof(Vector));

  sls = sizeX*sizeY;		// slice size
  sz = sls*sizeZ;

  if ( fread(volin, sizeof(DATATYPEIN), sz, filein) < (unsigned int)sz)
  {
    printf("File size is not the same as volume size\n");
    exit(1);
  }

  if ((fileout = fopen(argv[5],"w")) == NULL)
  {
    printf("Cannot open %s for writing\n",argv[5]);
    exit(1);
  }


  for (idx = 0; idx < sls*sizeZ; idx++)  {
	  if (volin[idx] > 0) {
	       volin[idx] = INTERIOR;
	  }
  }

  fclose(filein);

  // Compute force vectors

 // get all the boundary voxels, saved in array boundVoxs[]
  border = 1;
  for (k = border; k < sizeZ-border; k++)
    for (j = border; j < sizeY-border; j++)
      for (i = border; i < sizeX-border; i++)
      {
	     flagBound = 0;
         idx = k*sls + j*sizeX + i;

		// CASE 1: treat the inner layer
	    if (volin[idx] == 0) continue;
		//consider six face neighbors, if anyone is zero, it is a boundary voxel
		iidx = k*sls + j*sizeX + i-1;
		if (volin[iidx] == 0) flagBound = 1;
		iidx = k*sls + j*sizeX + i+1;
		if (volin[iidx] == 0) flagBound = 1;
		iidx = k*sls + (j-1)*sizeX + i;
		if (volin[iidx] == 0) flagBound = 1;
		iidx = k*sls + (j+1)*sizeX + i;
		if (volin[iidx] == 0) flagBound = 1;
		iidx = (k-1)*sls + j*sizeX + i;
		if (volin[iidx] == 0) flagBound = 1;
		iidx = (k+1)*sls + j*sizeX + i;
		if (volin[iidx] == 0) flagBound = 1;
	

		// CASE 2: treat the outer layer
	/*	if (volin[idx] != 0) continue;
		//consider six face neighbors, if anyone is zero, it is a boundary voxel
		iidx = k*sls + j*sizeX + i-1;
		if (volin[iidx] != 0) flagBound = 1;
		iidx = k*sls + j*sizeX + i+1;
		if (volin[iidx] != 0) flagBound = 1;
		iidx = k*sls + (j-1)*sizeX + i;
		if (volin[iidx] != 0) flagBound = 1;
		iidx = k*sls + (j+1)*sizeX + i;
		if (volin[iidx] != 0) flagBound = 1;
		iidx = (k-1)*sls + j*sizeX + i;
		if (volin[iidx] != 0) flagBound = 1;
		iidx = (k+1)*sls + j*sizeX + i;
		if (volin[iidx] != 0) flagBound = 1;
     */

		if (flagBound == 1) {
	//	 if (flagBound == 1 && i==8) {  // Testing: only see a slice in 2D
			numBound++;
			volin[idx] = SURF;
		 }
	  }

    boundVoxs = (VoxelPosition *)malloc( numBound * sizeof(VoxelPosition));       

	n = 0;
	for (k = 0; k < sizeZ; k++)
	   for (j = 0; j < sizeY; j++)
		  for (i = 0; i < sizeX; i++)  {
			  if (volin[k*sls +j*sizeX +i] == SURF) {
  					boundVoxs[n].x = i;
					boundVoxs[n].y = j;
					boundVoxs[n].z = k;
					n++;
			  }
		  }


	//printf("numBound = %d \n", numBound);
	
	createLUTable(vecLUTable);
	printf("Finished creating LUTable. \n");

	for (k = 0; k < sizeZ; k++)  {
	  for (j = 0; j < sizeY; j++)
		for (i = 0; i < sizeX; i++)
		{
		  idx = k*sls + j*sizeX + i;
		  #if OUTSIDE_FIELD
		  if (volin[idx] == SURF) {
		  #else
		  if (volin[idx] == 0 || volin[idx] == SURF) {
		  #endif
			force[idx].xd = 0;
			force[idx].yd = 0;
			force[idx].zd = 0;
			continue; //consider interior only
		  }

		  totalForce.xd = totalForce.yd = totalForce.zd = 0;
		  s = 0;
		  while(s < numBound) {
			    ss = s;
				#if SurfSampling
				s+= 7;      //3; only consider 1/3 surf points
				#else
				s++;
				#endif
				if(abs(boundVoxs[ss].x-i)>=ForceRange)
						 continue;
					else if(abs(boundVoxs[ss].y-j)>=ForceRange)
							  continue;
						 else if(abs(boundVoxs[ss].z-k)>=ForceRange)
								  continue;             // Not consider surface point too far
				pointForce = ptForceLUTable(vecLUTable, boundVoxs[ss].x, boundVoxs[ss].y, boundVoxs[ss].z, i, j, k);
			//if (i==4) printf("%d/%d/  %f \n", j, k, pointForce.xd);
				totalForce.xd += pointForce.xd;
				totalForce.yd += pointForce.yd;
				totalForce.zd += pointForce.zd;
		  }
		  force[idx].xd = totalForce.xd;
		  force[idx].yd = totalForce.yd;
		  force[idx].zd = totalForce.zd;

		//print out critical points
		  //float thresh_critpoint = 0.0001;
		//printf("print out critical points \n");
	//	 if (fabs(force[idx].xd)<thresh_critpoint && fabs(force[idx].yd)<thresh_critpoint
	//	                   && fabs(force[idx].zd)<thresh_critpoint)
	//		printf("%d %d %d - %f %f %f\n",i, j, k, force[idx].xd, force[idx].yd, force[idx].zd);
	//if (i==5 && j==7 && k==9)
	//      printf("force[5,7,9]= %f %f %f\n", force[idx].xd, force[idx].yd, force[idx].zd);
	//if(i==4) printf("%d/%d/ = %f %f %f\n",j,k,force[idx].xd,force[idx].yd,force[idx].zd);

		}
		printf(" processing slice # %d ..\n", k);
	}

	//print force vectors

	for (k = 0; k < sizeZ; k++)
	  for (j = 0; j < sizeY; j++)
		for (i = 0; i < sizeX; i++) {
			idx = k*sls + j*sizeX + i;
			if (volin[idx]==0 || volin[idx]==SURF) continue;  // not output
			// Test: generate 2D vector field -- add zero slice above and below for AVS display reason
		//if (i == 7)  fprintf(fileout,"%f %f %f\n", 0.0, 0.0, 0.0);
			//if (i == 8)  fprintf(fileout,"%f %f %f\n", 0.0, force[k*sls+j*sizeX+i].yd*5, force[k*sls+j*sizeX+i].zd*5);
		//if (i == 9)  fprintf(fileout,"%f %f %f\n", 0.0, 0.0, 0.0);
			fprintf(fileout, "%d %d %d %f %f %f\n", i, j, k, force[idx].xd*5, force[idx].yd*5, force[idx].zd*5);
	}


	//for (int w = 0; w < sizeZ; w++)  //print out one component of force vectors
	//   printf("%f==",force[w*sls + 7*sizeX + 5].zd);


   fclose(fileout);

   printf("End \n");

   return 0;
}

