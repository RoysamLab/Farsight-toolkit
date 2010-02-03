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
#if defined(_MSC_VER)
#pragma warning(disable : 4996)
#endif

#include "mdlIntegratedSkeleton.h"

namespace mdl
{

//Constructor
IntegratedSkeleton::IntegratedSkeleton()
{
	m_inputImage = NULL;
	debug = false;

	Iu = NULL;
	Iv = NULL;
	Iw = NULL;
	fc = NULL;
}

IntegratedSkeleton::~IntegratedSkeleton()
{
	this->clearGradientVectorField();
	if(curv)
	{
		delete[] curv;
		curv = NULL;
	}
	m_inputImage=NULL;
}

void IntegratedSkeleton::SetInput(ImageType::Pointer inImage)
{
	m_inputImage = inImage;
}

bool IntegratedSkeleton::Update()
{
	if(!this->createGradientVectorField())
		return false;
	if(!this->computeIsoGraySurfaceCurvature())
		return false;
	return true;
}

int IntegratedSkeleton::sign(float value)
{
	if (value > 0)
		return 1;
	else if(value < 0)
		return -1;
	else
		return 0;
}


int IntegratedSkeleton::sign1(float value)
{
	if (value > 1e-5)
		return 1;
	else if(value < -1e-5)
		return -1;
	else
		return 0;
}

float IntegratedSkeleton::veclength(Vector3D vecin)
{
	return sqrt(vecin.x*vecin.x + vecin.y*vecin.y + vecin.z*vecin.z);
}


bool IntegratedSkeleton::createGradientVectorField()
{
	ImageType::RegionType region = m_inputImage->GetBufferedRegion();
	int sizeX = region.GetSize(0);
	int sizeY = region.GetSize(1);
	int sizeZ = region.GetSize(2);
	long numPix = sizeX*sizeY*sizeZ;
	
	//Allocate memory for gradient field:
	this->clearGradientVectorField();
	Iu = new float[numPix];
	Iv = new float[numPix];
	Iw = new float[numPix];
	fc = new unsigned char[numPix];	//This keeps track of interior/exterior and surface voxels
	if(!Iu || !Iv || !Iw || !fc)
	{
		if(debug) std::cerr << "Unable to allocate memory for gradients" << std::endl;
		return false;
	}

	//Init fc
	int numBound = 0;
	itk::ImageRegionIteratorWithIndex< ImageType > itr( m_inputImage, region );
	for(itr.GoToBegin(); !itr.IsAtEnd(); ++itr)
	{
		ImageType::IndexType index = itr.GetIndex();
		long idx = index[2]*sizeY*sizeX + index[1]*sizeX + index[0];

		if(itr.Get() == 0)
			fc[idx] = 0;
		else
		{
			//if on edge of image continue
			if(index[0] <= 0 || index[1] <= 0 || index[2] <= 0)
			{
				fc[idx]=0;
				continue;
			}
			if(index[0]>=(sizeX-1) || index[1]>=(sizeY-1) || index[2]>=(sizeZ-1)) 
			{
				fc[idx]=0;
				continue;
			}

			//consider six face neighbors, 
			//if anyone is zero, it is a boundary voxel
			bool flagBound = false;
			for(int d=0; d<3; ++d)
			{
				ImageType::IndexType m_ind = index; //
				ImageType::IndexType p_ind = index;
				m_ind[d] -= 1;
				p_ind[d] += 1;
				ImageType::PixelType m_pix = m_inputImage->GetPixel(m_ind);
				ImageType::PixelType p_pix = m_inputImage->GetPixel(p_ind);
				if(m_pix == 0 || p_pix == 0) 
				{
					flagBound = true;
					break;
				}
			}

			if(flagBound)
			{
				numBound++;
				fc[idx] = SURFACE;
			}
			else
			{
				fc[idx] = INTERIOR;
			}
		} // end else (not exterior)
	} // end for(itr)

	if(debug)
		std::cerr << "Num Bounds = " << numBound << std::endl;


	// define positive half kernel of derivative (Sobel kernel)
	double kernelWeight[3][3];
	kernelWeight[0][0] = 1.67; kernelWeight[0][1] = 5.80; kernelWeight[0][2] = 1.67;
	kernelWeight[1][0] = 5.80; kernelWeight[1][1] = 20.1; kernelWeight[1][2] = 5.80;
	kernelWeight[2][0] = 1.67; kernelWeight[2][1] = 5.80; kernelWeight[2][2] = 1.67;
	float s = 0.002; //scale on outputs  0.002

	//Compute gradient of interior/surface voxels:
	ImageType::IndexType m_ind;
	ImageType::IndexType p_ind;
	ImageType::PixelType m_pix;
	ImageType::PixelType p_pix;
	for(itr.GoToBegin(); !itr.IsAtEnd(); ++itr)
	{
		ImageType::IndexType index = itr.GetIndex();
		long idx = index[2]*sizeY*sizeX + index[1]*sizeX + index[0];

		Iu[idx] = 0;
		Iv[idx] = 0;
		Iw[idx] = 0;

		if(fc[idx] == 0)
			continue; //consider interior only

		//if on edge of image continue
		if(index[0] <= 0 || index[1] <= 0 || index[2] <= 0) 
			continue;
		if(index[0]>=(sizeX-1) || index[1]>=(sizeY-1) || index[2]>=(sizeZ-1)) 
			continue;

		for (int d2 = -1; d2 <= 1; d2++)
		{
			for (int d1 = -1; d1 <= 1; d1++)  
			{
				//x gradient:
				m_ind = index;
				m_ind[2]+=d2;
				m_ind[1]+=d1;
				p_ind = m_ind;
				m_ind[0]-=1;
				p_ind[0]+=1;
				m_pix = m_inputImage->GetPixel(m_ind);
				p_pix = m_inputImage->GetPixel(p_ind);
				Iu[idx] += kernelWeight[d2+1][d1+1]*(p_pix - m_pix);

				//y gradient:
				m_ind = index;
				m_ind[2]+=d2;
				m_ind[0]+=d1;
				p_ind = m_ind;
				m_ind[1]-=1;
				p_ind[1]+=1;
				m_pix = m_inputImage->GetPixel(m_ind);
				p_pix = m_inputImage->GetPixel(p_ind);
				Iv[idx] += kernelWeight[d2+1][d1+1]*(p_pix - m_pix);

				//z gradient:
				m_ind = index;
				m_ind[1]+=d2;
				m_ind[0]+=d1;
				p_ind = m_ind;
				m_ind[2]-=1;
				p_ind[2]+=1;
				m_pix = m_inputImage->GetPixel(m_ind);
				p_pix = m_inputImage->GetPixel(p_ind);
				Iw[idx] += kernelWeight[d2+1][d1+1]*(p_pix - m_pix);
			} //end for d1
		} //end for d2
		
		Iu[idx]*= s;
		Iv[idx]*= s;
		Iw[idx]*= s;
	} // end for(itr)

	if(debug)	//write out the vector field
	{
		FILE *fileout;
		if ((fileout = fopen("out.vec","w")) != NULL)
		{
			for (int k = 0; k < sizeZ; k++) {
				for (int j = 0; j < sizeY; j++) {
					for (int i = 0; i < sizeX; i++) 
					{
						long idx = k*sizeX*sizeY + j*sizeX + i;
						if (fc[idx]==0)   continue;  // not output
						fprintf(fileout, "%d %d %d %f %f %f\n", i, j, k, Iu[idx], Iv[idx], Iw[idx]);
					} 
				} 
			} //end triple for loops
			fclose(fileout);
		}
	}

	return true;
}

void IntegratedSkeleton::clearGradientVectorField()
{
	if(Iu)
	{
		delete[] Iu;
		Iu = NULL;
	}
	if(Iv)
	{
		delete[] Iv;
		Iv = NULL;
	}
	if(Iw)
	{
		delete[] Iw;
		Iw = NULL;
	}
	if(fc)
	{
		delete[] fc;
		fc = NULL;
	}	
}

//------------------------------------------------------------------------//
//Compute iso-gray surface principle surface curvature    
//------------------------------------------------------------------------//
bool IntegratedSkeleton::computeIsoGraySurfaceCurvature()
{
	if(!Iu || !Iv || !Iw || !fc)
		return false;

	ImageType::RegionType region = m_inputImage->GetBufferedRegion();
	int L = region.GetSize(0); //x
	int M = region.GetSize(1); //y
	int N = region.GetSize(2); //z
	long numPix = L*M*N;

	//Vectors for partial derivatives:
	float *Iuu = new float[numPix];
	float *Ivv = new float[numPix];
	float *Iww = new float[numPix];
	float *Iuv = new float[numPix];
	float *Iuw = new float[numPix];
	float *Ivw = new float[numPix];

	if(!Iuu || !Ivv || !Iww || !Iuv || !Iuw || !Ivw)
	{
		if(debug)
			std::cerr << "Could not allocate memory for iso gray surface" << std::endl;
		return false;
	}

	if(debug)
		std::cerr << "Computing curvature" << std::endl;

	this->partialDerivative1(Iu, Iuu, 1, L, M, N);
	this->partialDerivative1(Iv, Ivv, 2, L, M, N);
	this->partialDerivative1(Iw, Iww, 3, L, M, N);
	this->partialDerivative1(Iu, Iuv, 2, L, M, N);
	this->partialDerivative1(Iu, Iuw, 3, L, M, N);
	this->partialDerivative1(Iv, Ivw, 3, L, M, N);

	//Temp containers:
	float RotMatrix[3][3];
	float RotMatrixTransp[3][3];
	float Hessian[3][3];
	float HessianPrime[3][3];
	float HessianPrime1[3][3];
	Vector3D gradient0;

	//Allocate memory for curvature:
	if(curv)
	{
		delete[] curv;
		curv = NULL;
	}
	curv = new float[numPix];


	int DisAway = 2;
	for (int k = DisAway; k < N-DisAway; k++)
    {
		for (int j = DisAway; j < M-DisAway; j++)
		{
			for (int i = DisAway; i < L-DisAway; i++)
			{
				long idx = k*L*M + j*L + i;
				if (fc[idx] != INTERIOR) continue;

				gradient0.x = Iu[idx];
				gradient0.y = Iv[idx];
				gradient0.z = Iw[idx];
				Hessian[0][0] = Iuu[idx];
				Hessian[0][1] = Iuv[idx];
				Hessian[0][2] = Iuw[idx];
				Hessian[1][0] = Iuv[idx];
				Hessian[1][1] = Ivv[idx];
				Hessian[1][2] = Ivw[idx];
				Hessian[2][0] = Iuw[idx];
				Hessian[2][1] = Ivw[idx];
				Hessian[2][2] = Iww[idx];
    
				ComputeRotMatrix(RotMatrix, gradient0);
				Transpose(RotMatrixTransp, RotMatrix);
				Matrix3Multiply(HessianPrime1, RotMatrixTransp, Hessian);
				Matrix3Multiply(HessianPrime, HessianPrime1, RotMatrix);

				double k1 = 0.5*(HessianPrime[1][1]+HessianPrime[2][2]) +
					0.5*sqrt(pow((Hessian[1][1]-Hessian[2][2]),2) +
					4*Hessian[1][2]*Hessian[2][1]);
                /* // the second eignvalue is not used, but please remain this code
				double k2 = 0.5*(HessianPrime[1][1]+HessianPrime[2][2]) -
					0.5*sqrt(pow((Hessian[1][1]-Hessian[2][2]),2) +
					4*Hessian[1][2]*Hessian[2][1]);
					*/
				float gLength = sqrt(gradient0.x*gradient0.x + 
					gradient0.y*gradient0.y +
					gradient0.z*gradient0.z);

				if(gLength !=0)
				{
					curv[idx] = (float) -(k1)/gLength;
				}
				else
				{
					//set large to be included in skeleton
					curv[idx]=9999;
				}

				//case 1: get neg maximum
				if (curv[idx]>10000)
				{
					curv[idx]=10000;
				}
				if (curv[idx]< -1 )
				{
					// first phase: threshold the curvature value (larger->fewer)
					curv[idx]= 0; 
				}
			}
		}
    }

	delete[] Iuu;
	delete[] Ivv;
	delete[] Iww;
	delete[] Iuv;
	delete[] Iuw;
	delete[] Ivw;

	return true;
}

// --- Computing partial Derivativel -------------------------//
void IntegratedSkeleton::partialDerivative1(float *Is, float *Isd, int direc, int L, int M, int N)  
{
	// direc = 1: x-direction   2: y-direction   3: z-direction
	long idx;
	int i,k,j;
	long slsz = L*M;
	float lessXsum, moreXsum, lessYsum, moreYsum, lessZsum, moreZsum;
	int disaway = 1;
	for (k = disaway; k < N-disaway; k++)
	{
		for (j = disaway; j < M-disaway; j++)
		{
			for (i = disaway; i < L-disaway; i++) 
			{
				idx = k*slsz + j*L +i;
				if (direc == 1) 
				{
					lessXsum = Is[k*slsz + j*L +(i-1)];
					moreXsum = Is[k*slsz + j*L +(i+1)];
					Isd[idx] = moreXsum - lessXsum;
				}
				else if(direc == 2) 
				{
					lessYsum = Is[k*slsz + (j-1)*L +i];
					moreYsum = Is[k*slsz + (j+1)*L +i];
					Isd[idx] = moreYsum - lessYsum;
				}
				else 
				{
					lessZsum = Is[(k-1)*slsz + j*L + i ];
					moreZsum = Is[(k+1)*slsz + j*L + i ];
					Isd[idx] = moreZsum - lessZsum;
				}
			} //end for i
		} // end for j
	} // end for k
}

//----------------optimized Rotation of Matrix -------------------------------------------------//
//------------  this code written by Xiao Liang -----------------------------------------//
//---------the idea is void computing cos(),sin() so many times ------------------//
void IntegratedSkeleton::ComputeRotMatrix(float RotateMatrix[3][3], Vector3D v)
{
	// Find the matrix that can rotate any vector v to x-axis direction
	float phi, theta;
	Vector3D vector1;
	float RotateMatrix1[3][3];
	float cosphi,sinphi,cospsi,sinpsi,costheta,sintheta;

	theta = atan2(v.z, v.y);
	//RotMatrixFromAngle(RotateMatrix1, 0, -theta, 0);
	cosphi=1;sinphi=0;costheta=cos(theta);sintheta=sin(theta);cospsi=1;sinpsi=0;
	this->RotMatrixFromAngle(RotateMatrix1,cosphi,sinphi,costheta,sintheta,cospsi,sinpsi);
	vector1.x = RotateMatrix1[0][0]*v.x +RotateMatrix1[0][1]*v.y +RotateMatrix1[0][2]*v.z;
	vector1.y = RotateMatrix1[1][0]*v.x +RotateMatrix1[1][1]*v.y +RotateMatrix1[1][2]*v.z;
	vector1.z = RotateMatrix1[2][0]*v.x +RotateMatrix1[2][1]*v.y +RotateMatrix1[2][2]*v.z;
	phi = atan2(vector1.y, vector1.x);
   
	cosphi=cos(-phi);sinphi=sin(-phi);
	// costheta=costheta; //it is not needed to reseted since them havenot changed
	sintheta=-sintheta;
	//cospsi=1;sinpsi=0;  //it is not needed to reseted since them havenot changed 
	RotMatrixFromAngle(RotateMatrix1,cosphi,sinphi,costheta,sintheta,cospsi,sinpsi);
	//RotMatrixFromAngle(RotateMatrix, -phi, -theta, 0);
}

void IntegratedSkeleton::RotMatrixFromAngle(float RMatrix[3][3], float cosphi, float sinphi, float costheta,float sintheta, float cospsi, float sinpsi) 
{
  // rotation matrix is      R(0,0) R(0,1) R(0,2)
  //                 R(1,0) R(1,1) R(1,2)
  //             R(2,0) R(2,1) R(2,2)
   RMatrix[0][0] = cosphi*cospsi - costheta*sinphi*sinpsi;
   RMatrix[1][0] = cospsi*sinphi + cosphi*costheta*sinpsi;
   RMatrix[2][0] = sinpsi*sintheta;
   RMatrix[0][1] = -cospsi*costheta*sinphi - cosphi*sinpsi;
   RMatrix[1][1] = cosphi*cospsi*costheta - sinphi*sinpsi;
   RMatrix[2][1] = cospsi*sintheta;
   RMatrix[0][2] = sinphi*sintheta;
   RMatrix[1][2] = -cosphi*sintheta;
   RMatrix[2][2] = costheta;
}

void IntegratedSkeleton::Transpose(float MatTransp[3][3], float Mat[3][3])  
{
   for(int j=0; j<=2; j++)
        for(int i=0; i<=2; i++)  {
             MatTransp[j][i] = Mat[i][j];
  }
}


void IntegratedSkeleton::Matrix3Multiply(float Mat[3][3], float Mat1[3][3], float Mat2[3][3]) 
{
   for(int j=0; j<=2; j++)
        for(int i=0; i<=2; i++)  {
             Mat[j][i] = Mat1[j][0]*Mat2[0][i] +Mat1[j][1]*Mat2[1][i] +Mat1[j][2]*Mat2[2][i];
  }
}


}



