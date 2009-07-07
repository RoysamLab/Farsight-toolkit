#ifndef __GPU_H__
#define __GPU_H__

#include <pthread.h>

#ifndef GPU
typedef double DOUBLE;
#else
typedef float DOUBLE;
#endif


/*
typedef struct {
   float Intensity;
   float Scale;
   DOUBLE Position[3];
} cSeedPoint3D;

typedef struct {
   long GridSpacing;
   int IntensityThreshold;
   int numSeeds;
   cSeedPoint3D **SeedContainer;
} cSeedContainer3D;
*/

typedef struct {
    double q3[4];
    double q2[4];
    double q1[4];
    DOUBLE a1;
    DOUBLE a2;
    DOUBLE a3;
    double f;
    double b;
    double mu[3];
    DOUBLE e1;
    DOUBLE e2;
    double last_dir[3];
    double last_pos[3];
    double L;
    double bndy[24];
    double type;
    double index;
    double MAD;
    long int ID;
    long int TraceID;
    double K1;
    double K2;
    double a;
    double k;
    double R1[3];
    double R2[3];
    double R3[3];
    long int NbrID[4];
    unsigned int numNbrs;
} cTVessel;

typedef struct {
    DOUBLE r[3][3];
} cTRMatrix;

typedef struct {
    DOUBLE p[3];
} cTVertex;

typedef struct {
    int vertex[3];
    DOUBLE centroid[3];
    DOUBLE normal[3];
    DOUBLE area;
} cTFacet;

typedef struct {
    DOUBLE dt2;
    DOUBLE dt;
    DOUBLE dt_u[3];
    DOUBLE dt_a[3];
    double sign_S[4];
    double sign_A[3];
    double sign_U[3];
} cTDamp;

/*
typedef struct {
   int id;
   int start;
   int end;
   long int dims[3];
   float *image;
   pthread_mutex_t *lock;
   cSeedContainer3D *seeds;
   cTVessel **SSContainer;
   int *SSContainerSize;
//   Seed2Seg::Pointer ss;
//   SeedContainer3D::Pointer seeds;
//   ImageType3D::Pointer image;
//   TraceConfig::Pointer config;
//   cTraceConfig *config;
//   long int *ID;
} ThreadArgs;
*/

#define MAX_SIZE 2000000
#define _PI 3.14159265

void cuda_update_e1(void*,void*,int,DOUBLE*,void*,void*,void*);

#endif 
