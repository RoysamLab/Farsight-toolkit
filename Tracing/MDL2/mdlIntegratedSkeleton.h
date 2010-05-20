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
#include "mdlUtils.h"

#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <algorithm>

#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkRecursiveGaussianImageFilter.h" 
#include "itkImageFileWriter.h" 


#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>

namespace mdl
{

#define MIN(x,y) (((x) < (y))?(x):(y))
#define MAX(x,y) (((x) > (y))?(x):(y))

class IntegratedSkeleton
{
public:
	IntegratedSkeleton(ImageType::Pointer inImage);
	~IntegratedSkeleton();
	//Setup:
	void SetDebug(bool inp = true){ debug = inp; };
	void SetVectorMagnitude(double vM){ vectorMagnitude = vM; };
	void SetLinePathStepSize(float v){ linePathStepSize = v; };

		//Methods:
	bool Update();

	//Optionally use Xiao Liang's method to compute surface curvature:
	void SetUseXiaoLiangMethod(bool inp = true){ useXiaoLiangMethod = inp; };
	//A TEST METHOD BY XIAO LIANG:
	//bool RunXiaoLSkeletonPoints(float sigma=0,int MoveStep =10);

	//Get Result:
	std::vector<fPoint3D> GetOutput(){ return skeletonPoints; };

private:
	typedef fPoint3D Vector3D;
	typedef fPoint3D VoxelPosition;

	//Parameters
	bool debug;				//If debug is true, process in steps and print stuff
	double vectorMagnitude;
	float linePathStepSize;
	bool useXiaoLiangMethod;	//For computing iso-gray surface curvature

	//Images & size
	ImageType::Pointer m_inputImage;
	ImageType::RegionType region;
	int sizeX;
	int sizeY;
	int sizeZ;
	long numPix;

	//Gradient Vectors:
	float *Iu;
	float *Iv;
	float *Iw;
	unsigned char *fc;	//100 for surface voxel, 200 for interior, 0 else
	Vector3D *force;

	//Surface Curvature
	float *curv;		//Iso-gray surface curvature

	//Seeds
	std::vector<fPoint3D> curvSeeds; //Maximum Curvature Seeds
	std::vector<fPoint3D> critSeeds; //Critical Point Seeds

	//Skeleton points (ouput)
	std::vector<fPoint3D> skeletonPoints;

	// New Skeleton points Computation,  add in March,2010 
	//bool XiaoLComputeSkeletonPoints(void);
	//bool MovingSkeletonPointsAlongGVF(int Step);
    
	//Key functions:
	bool createGradientVectorField(void);
	bool createGradientVectorField(float sigma);
	//bool createGradientVectorFieldWithITK(){return false;};
	bool XiaoLiangComputeIsoGraySurfaceCurvature();
	bool XiaosongComputeIsoGraySurfaceCurvature();
	bool computeSeedsWithMaxCurvature();
	bool convertGradVectorToForceVector();
	bool computeCriticalPointSeeds();
	bool computeSkeleton();
	void cleanUpMemory();

	//Internal methods
	void clearGradientVectorField();
	void partialDerivative1(float *Is, float *Isd, int direc, int L, int M, int N);
	void ComputeRotMatrix(float RotateMatrix[3][3], Vector3D v);
	void RotMatrixFromAngle(float RMatrix[3][3], float cosphi, float sinphi,float costheta, float sintheta, float cospsi,float sinpsi);
	void Transpose(float MatTransp[3][3], float Mat[3][3]);
	void Matrix3Multiply(float Mat[3][3], float Mat1[3][3], float Mat2[3][3]);
	Vector3D interpolation(float x, float y, float z, int sizx, int sizy, int sizz, Vector3D *forcevec);
	void rk2(float x, float y, float z, int sizx, int sizy, int sizz, float steps, Vector3D  *Force_ini, VoxelPosition *nextPos);

	double GetMean(float *buf, int nx, int ny, int nz);
	void printProgress(int pos, int tot);

	//useful functions:
	int sign(float value);
	int sign1(float value);
	float veclength(Vector3D vecin);

	static const unsigned char SURFACE = 100;
	static const unsigned char INTERIOR = 200;
};

}  // end namespace mdl

#endif
