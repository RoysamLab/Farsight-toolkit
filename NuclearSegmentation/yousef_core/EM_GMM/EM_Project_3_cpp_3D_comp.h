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

//added by Yousef on 09/06/2009
void Initialize_Parameters_V2(std::vector<std::vector<double> > *X,int** SEEDS,double** U,double*** Segma,double* PI,double** Z,int num_points,int num_components);
void Param_Maximization_V2(std::vector<std::vector<double> > *X,double **Gamma,double** U,double*** Segma,double* PI, int num_points, int num_components);
void AssignToComponent(std::vector<std::vector<double> > *X, int num_points, int num_components, double** U,double*** Segma,double* PI);

#endif

