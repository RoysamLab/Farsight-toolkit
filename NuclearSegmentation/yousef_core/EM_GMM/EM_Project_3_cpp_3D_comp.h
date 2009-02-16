#ifndef _EM_PROJECT_CPP_3D_COMP_H_
#define _EM_PROJECT_CPP_3D_COMP_H_

#include<math.h>
#include<iostream>
#include<limits>
#include <vector>

using namespace std;

#define MAX_SZ 50

void Initialize_Parameters(std::vector<std::vector<double> > *X,int** SEEDS,double** U,double*** Segma,double* PI,double** Z,int num_points,int num_components);

double Multivar_Norm3D(double X, double Y, double Z, double Ux, double Uy, double Uz, double S00, double S01, double S11, double S02, double S12, double S22);

void Post_Expect(double **Gamma, std::vector<std::vector<double> > *X,double *PI,double **U,double ***Segma,int num_points,int num_components);

void Param_Maximization(std::vector<std::vector<double> > *X,double **Gamma,double** U,double*** Segma,double* PI, int num_points, int num_components);

void Compute_Probabilities(std::vector<std::vector<double> > *X, int num_points, double** U,double*** Segma,double* PI);

void EM_Gmm(std::vector<std::vector<double> > *X,int** SEEDS,int num_points, int num_components);


#endif

