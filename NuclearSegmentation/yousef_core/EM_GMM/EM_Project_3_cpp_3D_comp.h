/* 
 * Copyright 2009 Rensselaer Polytechnic Institute
 * This program is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
 * for more details.
 * 
 * You should have received a copy of the GNU General Public License along 
 * with this program; if not, write to the Free Software Foundation, Inc., 
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#ifndef _EM_PROJECT_CPP_3D_COMP_H_
#define _EM_PROJECT_CPP_3D_COMP_H_

#include<math.h>
#include<iostream>
#include<limits>
#include <vector>

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

