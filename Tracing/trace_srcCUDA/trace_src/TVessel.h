#ifndef TVESSEL_H
#define TVESSEL_H

#include <iostream>
#include "vnl/vnl_math.h"
#include "gpu.h"

class TVessel
{

public:

	TVessel();
	~TVessel();

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
    double e2;
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

    void PrintSelf();
    bool IsSegmentValid(TVessel*);
    TVessel* CopyToNewSegment();

};

typedef struct _tagMatrix
{
    double r[3][3];
}TRMatrix;


#endif
