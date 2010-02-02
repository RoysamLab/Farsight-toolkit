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
/**************************************************************************  
// Extract the ridges and valleys feature in 3D vector field
// --- Input: 3D vector field
// --- Output: 3D image-scalar field
// --- Author: Xiaosong Yuan, improved and optimized by xiao liang, RPI
// --- Modified Date: 03/Sep/2009
 *  Adapted Jan. 2010 by Isaac Abbott 
 *        
 *************************************************************************/
#ifndef __mdlIntegratedSkeleton_h
#define __mdlIntegratedSkeleton_h

#include "mdlTypes.h"

#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include <algorithm>

#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"

namespace mdl
{

#define MIN(x,y) (((x) < (y))?(x):(y))
#define MAX(x,y) (((x) > (y))?(x):(y))

class IntegratedSkeleton
{
public:
	struct  Vector3D{float x; float y; float z;};

	IntegratedSkeleton();
	~IntegratedSkeleton();
	//Setup:
	void SetVectorMagnitude(double vM){ vectorMagnitude = vM; };
	void SetDebug(bool inp = true){ debug = inp; };
	void SetInput(ImageType::Pointer inImage);
	//Methods:
	bool Update();

	//Get Result:

private:
	//Parameters
	bool debug;				//If debug is true, process in steps and print stuff
	double vectorMagnitude;

	//Images
	ImageType::Pointer m_inputImage;

	//Gradient Vectors:
	float *Iu;
	float *Iv;
	float *Iw;
	unsigned char *fc;	//100 for surface voxel, 200 for interior, 0 else

	//Surface Curvature
	float *curv;		//Iso-gray surface curvature

	//Key functions:
	void clearGradientVectorField();
	bool createGradientVectorField();
	bool createGradientVectorFieldWithITK(){return false;};
	bool computeIsoGraySurfaceCurvature();

	//Internal methods
	void partialDerivative1(float *Is, float *Isd, int direc, int L, int M, int N);
	void ComputeRotMatrix(float RotateMatrix[3][3], Vector3D v);
	void RotMatrixFromAngle(float RMatrix[3][3], float cosphi, float sinphi,float costheta, float sintheta, float cospsi,float sinpsi);
	void Transpose(float MatTransp[3][3], float Mat[3][3]);
	void Matrix3Multiply(float Mat[3][3], float Mat1[3][3], float Mat2[3][3]);
	
	//useful functions:
	int sign(float value);
	int sign1(float value);
	float veclength(Vector3D vecin);

	static const unsigned char SURFACE = 100;
	static const unsigned char INTERIOR = 200;
};

}  // end namespace mdl

#endif
