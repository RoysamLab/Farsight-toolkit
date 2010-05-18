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
IntegratedSkeleton::IntegratedSkeleton(ImageType::Pointer inImage)
{
	m_inputImage = inImage;

	region = m_inputImage->GetBufferedRegion();
	sizeX = region.GetSize(0);
	sizeY = region.GetSize(1);
	sizeZ = region.GetSize(2);
	numPix = sizeX*sizeY*sizeZ;

	debug = false;
	useXiaoLiangMethod = false;
	linePathStepSize = 0.5;

	Iu = NULL;
	Iv = NULL;
	Iw = NULL;
	fc = NULL;
	curv = NULL;
	force = NULL;

	curvSeeds.clear();
	critSeeds.clear();
	skeletonPoints.clear();
}

IntegratedSkeleton::~IntegratedSkeleton()
{
	this->cleanUpMemory();
	m_inputImage=NULL;
	skeletonPoints.clear();
}

void IntegratedSkeleton::cleanUpMemory()
{
	this->clearGradientVectorField();
	if(curv)
	{
		delete[] curv;
		curv = NULL;
	}
	if(fc)
	{
		delete[] fc;
		fc = NULL;
	}
	if(force)
	{
		delete[] force;
		force = NULL;
	}
	curvSeeds.clear();
	critSeeds.clear();
}

bool IntegratedSkeleton::Update()
{
	if(!this->createGradientVectorField())
		return false;
    //if(!createGradientVectorField(3))
	//	return false;

	if (useXiaoLiangMethod)
	{
		if(!this->XiaoLiangComputeIsoGraySurfaceCurvature())
			return false;
	}
	else 
	{   
		if(!this->XiaosongComputeIsoGraySurfaceCurvature())
			return false;
	}

	if(curvSeeds.size() == 0)
		this->computeSeedsWithMaxCurvature();
	if(curvSeeds.size() == 0)
		return false;

	if(!this->convertGradVectorToForceVector())
		return false;

	this->computeCriticalPointSeeds(); //If fails we are not dead in the water

	return this->computeSkeleton();
	//return this->MovingSkeletonPointsAlongGVF(10);
}


bool IntegratedSkeleton::RunXiaoLSkeletonPoints(float sigma, int MoveStep)
{
	if (sigma==0)
	{
       if(!this->createGradientVectorField()) // have two choice
		return false;
	}
	else 
	{
	   if(!this->createGradientVectorField(sigma)) // have two choice
		return false;
	}

	if(!this->XiaoLComputeSkeletonPoints()) // Hessian based Critical points Detection
		return false;
	
	if(!this->convertGradVectorToForceVector()) // normarize the GVF
		return false;
    
	if(!this->MovingSkeletonPointsAlongGVF(MoveStep)) // Move to the centle-line
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


bool IntegratedSkeleton::createGradientVectorField(void)
{
	if(!m_inputImage)
		return false;
	
	//Allocate memory for gradient field:
	this->clearGradientVectorField();
	Iu = new float[numPix];
	Iv = new float[numPix];
	Iw = new float[numPix];
	if(fc)
		delete[] fc;
	fc = new unsigned char[numPix];	//This keeps track of interior/exterior and surface voxels
	if(!Iu || !Iv || !Iw || !fc)
	{
		if(debug) std::cerr << "Unable to allocate memory for gradients" << std::endl;
		return false;
	}

	//Init variables to zeros:
	for(int i=0; i<numPix; ++i)
	{
		Iu[i] = 0;
		Iv[i] = 0;
		Iw[i] = 0;
		fc[i] = 0;
	}

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
			/* 6 connectivity
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
			*/
			// 18 connectivity
			for (int kk=-1; kk<=1; kk++)
            {
				for (int jj=-1; jj<=1; jj++)
				{
					for (int ii=-1; ii<=1; ii++)
					{
						ImageType::IndexType ind = index;
						ind[0] += ii;
						ind[1] += jj;
						ind[2] += kk;
						ImageType::PixelType pix = m_inputImage->GetPixel(ind);
						if(pix == 0)
						{
							flagBound = true;
							break;
						}
					}
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
	double s = 0.; 
	for (int i=0;i<3;i++)
		for (int j=0;j<3;j++)
			s += kernelWeight[i][j];

	s= 1/s; // normarized to 1;

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
						fprintf(fileout, "%d %d %d %4.4f %4.4f %4.4f\n", i, j, k, Iu[idx], Iv[idx], Iw[idx]);
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
}

//------------------------------------------------------------------------//
//Compute iso-gray surface principle surface curvature    
//------------------------------------------------------------------------//
bool IntegratedSkeleton::XiaosongComputeIsoGraySurfaceCurvature()
{
	
	if(debug)
		std::cout << "This is xiaosong' method!" << std::endl;
	if(!Iu || !Iv || !Iw || !fc || !m_inputImage)
		return false;

	int L = sizeX; //x
	int M = sizeY; //y
	int N = sizeZ; //z

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

	//Init partial derivs to zeros:
	for(int i=0; i<numPix; ++i)
	{
		Iuu[i]=0;
		Ivv[i]=0;
		Iww[i]=0;
		Iuv[i]=0;
		Iuw[i]=0;
		Ivw[i]=0;
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

	//Init temp containers to zeros:
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			RotMatrix[i][j] = 0;
			RotMatrixTransp[i][j] = 0;
			Hessian[i][j] = 0;
			HessianPrime[i][j] = 0;
			HessianPrime1[i][j] = 0;
		}
	}

	//Allocate memory for curvature:
	if(curv)
	{
		delete[] curv;
		curv = NULL;
	}
	curv = new float[numPix];
	if(!curv)
		return false;

	//Init to zeros
	for (long k = 0; k < numPix; ++k)
		curv[k] = 0;

	int border = 2;	//Why does this need to by 2?
	for (int k = border; k < N-border; k++)
    {
		for (int j = border; j < M-border; j++)
		{
			for (int i = border; i < L-border; i++)
			{
				long idx = k*L*M + j*L + i;
				if (fc[idx] != INTERIOR)
					continue;

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

				double term1 = 0.5*(HessianPrime[1][1]+HessianPrime[2][2]);
				double pwr = pow( (double)(Hessian[1][1] - Hessian[2][2]), 2.0 );
				double temp = pwr + 4.0*Hessian[1][2]*Hessian[2][1];
				double term2 = 0.5*sqrt(temp);
				double k1 = term1 + term2;

                /* // the second eignvalue is not used, but please remain this code
				double k2 = 0.5*(HessianPrime[1][1]+HessianPrime[2][2]) -
					0.5*sqrt(pow((Hessian[1][1]-Hessian[2][2]),2) +
					4*Hessian[1][2]*Hessian[2][1]);
					*/

				float gLength = sqrt(gradient0.x*gradient0.x + 
					gradient0.y*gradient0.y +
					gradient0.z*gradient0.z);

				if(gLength !=0 )
				{
					curv[idx] = (float) -(k1)/gLength;
				}
				/*else
				{
					//set large to be included in skeleton
					curv[idx]=9999;
				}

				//case 1: get neg maximum
				if (curv[idx]>10000)
				{
					curv[idx]=10000;
				}*/
				else
				{
					curv[idx]=0;
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

	if(debug)
	{
		FILE *fileout;
		if ((fileout = fopen("out.curv","w")) != NULL)
		{
			for (int k = 0; k < N; k++) {
				for (int j = 0; j < M; j++) {
					for (int i = 0; i < L; i++) 
					{
						long idx = k*L*M + j*L + i;
						fprintf(fileout, "%d %d %d %f\n", i, j, k, curv[idx]);
					} 
				} 
			} //end triple for loops
			fclose(fileout);
		}
	}

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
	cosphi=1;
	sinphi=0;
	costheta=cos(theta);
	sintheta=sin(theta);
	cospsi=1;
	sinpsi=0;
	this->RotMatrixFromAngle(RotateMatrix1,cosphi,sinphi,costheta,sintheta,cospsi,sinpsi);
	vector1.x = RotateMatrix1[0][0]*v.x +RotateMatrix1[0][1]*v.y +RotateMatrix1[0][2]*v.z;
	vector1.y = RotateMatrix1[1][0]*v.x +RotateMatrix1[1][1]*v.y +RotateMatrix1[1][2]*v.z;
	vector1.z = RotateMatrix1[2][0]*v.x +RotateMatrix1[2][1]*v.y +RotateMatrix1[2][2]*v.z;
	phi = atan2(vector1.y, vector1.x);
   
	cosphi=cos(-phi);
	sinphi=sin(-phi);
	// costheta=costheta; //it is not needed to reseted since them havenot changed
	sintheta=-sintheta;
	//cospsi=1;sinpsi=0;  //it is not needed to reseted since them havenot changed 
	RotMatrixFromAngle(RotateMatrix,cosphi,sinphi,costheta,sintheta,cospsi,sinpsi);
}

void IntegratedSkeleton::RotMatrixFromAngle(float RMatrix[3][3], float cosphi, float sinphi, float costheta,float sintheta, float cospsi, float sinpsi) 
{
  // rotation matrix is      R(0,0) R(0,1) R(0,2)
  //						R(1,0) R(1,1) R(1,2)
  //						R(2,0) R(2,1) R(2,2)
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

bool IntegratedSkeleton::computeSeedsWithMaxCurvature()
{
	curvSeeds.clear();;

	if(!curv || !m_inputImage)
		return false;

	int L = sizeX; //x
	int M = sizeY; //y
	int N = sizeZ; //z
	long slsz = L*M;

	double meanCurvature = this->GetMean(curv, L, M, N);

	//the threshold is a value between [0,10].
	double highCurvatureThreshold = (meanCurvature>10 ? 10 : meanCurvature);
	if(debug)
	{
		std::cerr << "mean = " << meanCurvature << std::endl;
		std::cerr << "estimated highCurvatureThreshold = " << highCurvatureThreshold << std::endl;
	}

	//-------------compute seeds with local maiximal curvature----------------//
	int DisAway = 1;
	// only select the local maximum voxel of curv[idx] as skeleton
	for (int k = DisAway; k < N-DisAway; k++)
    {
		for (int j = DisAway; j < M-DisAway; j++)
		{
			for (int i = DisAway; i < L-DisAway; i++)
			{
				long idx = k*slsz + j*L + i;
				long iidx;

				if (curv[idx] < highCurvatureThreshold )
				{
					// thresholding compared with threshold 
					curv[idx]= 0;
				}

				//consider the six face neighbors 
				//(if any are greater break out of loop)
				iidx = k*slsz + j*L + i-1;
				if (curv[iidx] > curv[idx])
				{
					continue;
				}
				iidx = k*slsz + j*L + i+1;
				if (curv[iidx] > curv[idx])
				{
					continue;
				}
				iidx = k*slsz + (j-1)*L + i;
				if (curv[iidx] > curv[idx])
				{
					continue;
				}
				iidx = k*slsz + (j+1)*L + i;
				if (curv[iidx] > curv[idx])
				{
					continue;
				}
				iidx = (k-1)*slsz + j*L + i;
				if (curv[iidx] > curv[idx])
				{
					continue;
				}
				iidx = (k+1)*slsz + j*L + i;
				if (curv[iidx] > curv[idx])
				{
					continue;
				}

				//curv[idx] is bigger than all 6 neighbors:
				if (curv[idx] != 0)
				{
					fPoint3D newSeed;
					newSeed.x =(float)i;
					newSeed.y =(float)j;
					newSeed.z =(float)k;
					curvSeeds.push_back(newSeed);
				} //end if (add seed)
			}// end for i
		}// end for j
    }// end for k

	if(debug)
	{
		std::cerr << "# Seeds with high curvature = " << curvSeeds.size() << std::endl;
		std::cerr << "Writing curvSeeds out to 'out.seed'" << std::endl;
		FILE *fileout;
		if ((fileout = fopen("out.seed","w")) != NULL)
		{
			for (int k = 0; k < (int)curvSeeds.size(); k++) 
			{
				fPoint3D s = curvSeeds.at(k);
				fprintf(fileout,"%d %d %d\n", (int)(s.x), (int)(s.y), (int)(s.z));
			}
			fclose(fileout);
		}
	}
	return true;
}

double IntegratedSkeleton::GetMean(float *buf, int nx, int ny, int nz)
{
	if(!buf)
		return 0.0;

	double mean=0.;
	if(nz * ny * nx == 1463*1033*75) // for the neocortical layer 6 
		return 0.02;
	
	// only considering the pixel which locate on the segmented object
	long InterVox =1; 
	for (int k = 0; k < nz; k++)
    {
		for (int j = 0; j < ny; j++)
		{
			for (int i = 0; i < nx; i++) 
			{
				long idx = k*nx*ny + j*nx + i;
				if(this->fc[idx]!=0)
				{  mean += buf[idx];
				   InterVox++;
				}
			}
		}
    }
	//mean /= (double)(nx*ny*nz); 
	mean /= (double)InterVox;
	mean = (mean < 0) ? 0 : mean;
	return mean;
}


bool IntegratedSkeleton::convertGradVectorToForceVector()
{
	if(!Iu || !Iv || !Iw || !m_inputImage)
		return false;

	int L = sizeX; //x
	int M = sizeY; //y
	int N = sizeZ; //z
	long slsz = L*M;

	if(force)
		delete[] force;
	force = new Vector3D[L*M*N];
	if(!force)
		return false;


	for (int k = 0; k < N; k++)
    {
		for (int j = 0; j < M; j++)
		{
			for (int i = 0; i < L; i++)
			{
				long idx = k*slsz + j*L +i;
				if (this->fc[idx] == 0)
				{
                 force[idx].x = 0;
				 force[idx].y = 0;
				 force[idx].z = 0;
				}
				else
				{
				 force[idx].x = Iu[idx];
				 force[idx].y = Iv[idx];
				 force[idx].z = Iw[idx];
				 float length = this->veclength(force[idx]);
				 force[idx].x = force[idx].x / (length+0.0001);
				 force[idx].y = force[idx].y / (length+0.0001);
				 force[idx].z = force[idx].z / (length+0.0001);
				}
				// this is just a very simple method, 
				// you can use Xiaosong's method i.e POW(),
				// and also you can use GDF method
			}
		}
	}
	//Delete grad vector field (free memory)!!
	this->clearGradientVectorField();
	return true;
}

void IntegratedSkeleton::printProgress(int pos, int tot)
{
	if( ( ( (float)pos / float(tot) ) >= 0.25 ) &&
        ( ( (float)(pos-1) / float(tot) ) < 0.25 ) )
	{
		std::cerr << "25% complete" << std::endl;
	}
    if( ( ( (float)pos / float(tot) ) >= 0.50 ) &&
        ( ( (float)(pos-1) / float(tot) ) < 0.50 ) )
	{
		std::cerr << "50% complete" << std::endl;
	}
    if( ( ( (float)pos / float(tot) ) >= 0.75 ) &&
        ( ( (float)(pos-1) / float(tot) ) < 0.75 ) )
	{
		std::cerr << "75% complete" << std::endl;
	}
	if( ( ( (float)pos / float(tot) ) >= 1.00 ) &&
        ( ( (float)(pos-1) / float(tot) ) < 1.00 ) )
	{
		std::cerr << "100% complete" << std::endl;
	}
}

// find all critical points: use points with small length of vector
bool IntegratedSkeleton::computeCriticalPointSeeds()
{
	int grid = 10;
	int maxDiv = -1; //max divergence: 1(good for phantom), 0(better), -1(not always good)

	critSeeds.clear();

	if(!force || !m_inputImage)
		return false;

	long slsz = sizeX*sizeY;

	float totalVecLength = 0;
	double divx,divy,divz,div;
	long idx;

	for (int k = 1; k < sizeZ-1; k++)
    {
		//this loop can take a long time to complete.
		//provide some feedback to the user on how we're doing
		if(debug)
			this->printProgress(k, sizeZ);
		
		for (int j = 1; j < sizeY-1; j++)
		{
			for (int i = 1; i < sizeX-1; i++)
			{
				totalVecLength = 0;
				divx = 0;  
				divy = 0;  
				divz = 0;

				//To check whether the position is near to boundaries, compute
				//divergence in x,y,z respectively
				idx = k*slsz + j*sizeX +i;
				if (fc[idx] == 0) 
					continue;
				totalVecLength += veclength(force[idx]);
				divx -= force[idx].x;
				divy -= force[idx].y;
				divz -= force[idx].z;

				idx = k*slsz + j*sizeX  +(i+1);
				if (fc[idx] == 0) 
					continue;
				totalVecLength += veclength(force[idx]);
				divx += force[idx].x;
				divy -= force[idx].y;
				divz -= force[idx].z;
	
				idx = k*slsz +(j+1)*sizeX + i;
				if (fc[idx] == 0)
					continue;
				totalVecLength += veclength(force[idx]);
				divx -= force[idx].x;
				divy += force[idx].y;
				divz -= force[idx].z;

				idx = k*slsz +(j+1)*sizeX + (i+1);
				if (fc[idx] == 0) 
					continue;
				totalVecLength += veclength(force[idx]);
				divx += force[idx].x;
				divy += force[idx].y;
				divz -= force[idx].z;

				idx = (k+1)*slsz + j*sizeX  +i;
				if (fc[idx] == 0) 
					continue;
				totalVecLength += veclength(force[idx]);
				divx -= force[idx].x;
				divy -= force[idx].y;
				divz += force[idx].z;
        
				idx = (k+1)*slsz + j*sizeX  +(i+1);
				if (fc[idx] == 0) 
					continue;
				totalVecLength += veclength(force[idx]);
				divx += force[idx].x;
				divy -= force[idx].y;
				divz += force[idx].z;

				idx = (k+1)*slsz +(j+1)*sizeX + i;
				if (fc[idx] == 0) 
					continue;
				totalVecLength += veclength(force[idx]);
				divx -= force[idx].x;
				divy += force[idx].y;
				divz += force[idx].z;

				idx = (k+1)*slsz +(j+1)*sizeX + (i+1);
				if (fc[idx] == 0)
					continue;
				totalVecLength += veclength(force[idx]);
				divx += force[idx].x;
				divy += force[idx].y;
				divz += force[idx].z;

				if (totalVecLength < 4)// skip the zero vector areas
					continue;

				div = divx + divy + divz;
				if (div > maxDiv)
				{
					//skip if the cube possibly contain a repelling point
					//(divergence is too big)
					continue;
				}

				Vector3D OutForce;
				for (int kk = 0; kk<grid; kk++)
				{
					for (int jj = 0; jj<grid; jj++)
					{
						for (int ii = 0; ii<grid; ii++)
						{
							float fx = i+ii/(float)grid;
							float fy = j+jj/(float)grid;
							float fz = k+kk/(float)grid;
							
						//Don't add seeds that are outside of image:
							if( (fx>=sizeX-1) || (fy>=sizeY-1) || (fz>=sizeX-1) )
								continue;
							if( (fx<0) || (fy<0) || (fz<0) )
								continue;
                        
					    // don't add newseeds that are within in a voxel;
						   fPoint3D dif;
						   dif.x = fx-i;
						   dif.y = fy-j;
						   dif.z = fz-k;
						   if (veclength(dif)<1)
							   continue;

							OutForce = interpolation(fx, fy, fz, sizeX, sizeY, sizeZ, force);
							if(veclength(OutForce) < vectorMagnitude)
							{
								fPoint3D newSeed;
								newSeed.x = fx;
								newSeed.y = fy;
								newSeed.z = fz;
								critSeeds.push_back(newSeed);
							} //end if less than max vector  magnitude
						}//end for ii
					} //end for jj
				}//end for kk
			}//end for i
		}//end for j
    }//end for k

	if(debug)
	{
		std::cerr << "Number of critical points is: " << critSeeds.size() << std::endl;
	}
	return true;
}

IntegratedSkeleton::Vector3D IntegratedSkeleton::interpolation(float x, float y, float z, int sizx, int sizy, int sizz, Vector3D *forcevec)
{
	//-------by xiao liang:  new implementation  
	

	//ADDED BY ISAAC 2/04/2010
	//Need to make sure that none of the Nei locations are outside the image:
	if(x >= sizx-1)
	{
		x = sizx-2;
	}
	if(y >= sizy-1)
	{
		y = sizy-2;
	}
	if(z >= sizz-1)
	{
		z = sizz-2;
	}
	if(x<0)
		x=0;
	if(y<0)
		y=0;
	if(z<0)
		z=0;

	int Intx=int(x);
	int Inty=int(y);
	int Intz=int(z);
	
	float alpha = x - Intx;   
	float beta = y - Inty;
	float gamma = z -Intz;
	
	long slsz = sizy*sizx;

	float a[8];// for interpolation coefficients
	a[0] = (1-alpha)*(1-beta)*(1-gamma);
	a[1] = (1-alpha)*(1-beta)*gamma;
	a[2] = (1-alpha)*beta*(1-gamma);
	a[3] =  alpha*(1-beta)*(1-gamma); 
	a[4] =  alpha*(1-beta)*gamma;
	a[5] =  alpha*beta*(1-gamma);
	a[6] =  (1-alpha)*beta*gamma;
	a[7] =  (alpha*beta*gamma);

	long Nei[8]; // considering N8 neighborhood in 3D image
    Nei[0] = Intz*slsz + Inty*sizx + Intx;
	Nei[1] =  Nei[0] + slsz;
	Nei[2] =  Nei[0] + slsz + sizx;
	Nei[3] =  Nei[0] + 1;
	Nei[4] =  Nei[0] + slsz + 1;
	Nei[5] =  Nei[0] + sizx + 1;
	Nei[6] =  Nei[0] + sizx;
	Nei[7] =  Nei[0] + slsz + sizx + 1;

	Vector3D forceInt;

	//------ compute interpolation ---------------------------------------//
	forceInt.x=forcevec[Nei[0]].x*a[0]
      +forcevec[Nei[1]].x*a[1]
      +forcevec[Nei[2]].x*a[2]
      +forcevec[Nei[3]].x*a[3]
      +forcevec[Nei[4]].x*a[4]
      +forcevec[Nei[5]].x*a[5]
      +forcevec[Nei[6]].x*a[6]
      +forcevec[Nei[7]].x*a[7];


	forceInt.y=forcevec[Nei[0]].y*a[0]
      +forcevec[Nei[1]].y*a[1]
      +forcevec[Nei[2]].y*a[2]
      +forcevec[Nei[3]].y*a[3]
      +forcevec[Nei[4]].y*a[4]
      +forcevec[Nei[5]].y*a[5]
      +forcevec[Nei[6]].y*a[6]
      +forcevec[Nei[7]].y*a[7];

	forceInt.z=forcevec[Nei[0]].z*a[0]
      +forcevec[Nei[1]].z*a[1]
      +forcevec[Nei[2]].z*a[2]
      +forcevec[Nei[3]].z*a[3]
      +forcevec[Nei[4]].z*a[4]
      +forcevec[Nei[5]].z*a[5]
      +forcevec[Nei[6]].z*a[6]
      +forcevec[Nei[7]].z*a[7];

	return(forceInt);
}

bool IntegratedSkeleton::computeSkeleton()
{
	skeletonPoints.clear();

	if(!m_inputImage)
		return false;

	long slsz = sizeX*sizeY;
	long sz = numPix;
	
	//Combine both types of seeds:
	std::vector<fPoint3D> seeds = curvSeeds;
	seeds.insert(seeds.end(),critSeeds.begin(),critSeeds.end());

	int numSeeds = (int)seeds.size();
	int numBoundSeeds = (int)curvSeeds.size();

	int idxSeeds = numSeeds-1;			//the current seed index

	VoxelPosition Startpos, Nextpos;
	int streamSteps=0;

	bool *FlagOnSkeleton = new bool[sz];
	if(!FlagOnSkeleton)		//couldn't allocate
		return false;
	for(int i=0; i<sz; i++)	//init to false
		FlagOnSkeleton[i] = false;

	while(idxSeeds >= 0)
	{
		Startpos.x = seeds[idxSeeds].x; 
		Startpos.y = seeds[idxSeeds].y; 
		Startpos.z = seeds[idxSeeds].z;
		long idx = (int)Startpos.z * slsz + (int)Startpos.y *sizeX + (int)Startpos.x;

		//Check whether the high curv points are within 
		// D-voxel distance to the existing skeleton
		int Dvoxel = 2;    // 3 is good for many results
		if (idxSeeds < numBoundSeeds)	//These seeds always come first in seeds list
		{
			Dvoxel= 2;  //4; 
			// 10 <- Found too big on April 12, 2006, which cause less spines
            // 3 is good for tlapse330, 
			// 6 is good for Trach
		}

		int FlagWithin = 0;
		if (idxSeeds < numBoundSeeds) 
		{
			for (int kk = -Dvoxel+1; kk <= Dvoxel-1; kk++)
			{
				for (int jj = -Dvoxel+1; jj <= Dvoxel-1; jj++)
				{
					for (int ii = -Dvoxel+1; ii <= Dvoxel-1; ii++) 
					{
						long iidx = idx + kk*slsz + jj*sizeX + ii;
						if (iidx < 0 || iidx >= sz)  
							continue;
						if(FlagOnSkeleton[iidx] == true)  
							FlagWithin = 1;
					}//end for ii
				}//end for jj
			} //end for kk
	
			if(FlagWithin == 1) 
			{
				idxSeeds--;
				continue;
			}// end if  FlagWithin
		}// end if (idxSeeds < numSeeds)

		VoxelPosition newPos;	//To be added to skeleton points
		newPos.x = Startpos.x;
		newPos.y = Startpos.y;
		newPos.z = Startpos.z;

		// being able to not show critical point in the streamlines
		if (idxSeeds >= numBoundSeeds)
		{
			//Add the start pos to the skeleton points
			skeletonPoints.push_back(newPos);
		}
		else
		{
			skeletonPoints.push_back(newPos);
		}

		//-------Line path algorithm------------------------------------------------//
		FlagOnSkeleton[idx] = true;

		while(streamSteps < 4000) //4000  
		{
			rk2(Startpos.x, Startpos.y, Startpos.z, sizeX, sizeY, sizeZ, linePathStepSize, force, &Nextpos);
			
			//If at edge of image get out of loop
			if(Nextpos.x >= sizeX-1 || Nextpos.y >= sizeY-1 || Nextpos.z >= sizeZ-1)
				break;
			if(Nextpos.x < 0 || Nextpos.y < 0 || Nextpos.z < 0)
				break;

			streamSteps++;
			Startpos.x = Nextpos.x;
			Startpos.y = Nextpos.y;
			Startpos.z = Nextpos.z;

			idx = (int)Nextpos.z *slsz + (int)Nextpos.y *sizeX + (int)Nextpos.x;
			if (FlagOnSkeleton[idx] != true) 
			{
				fPoint3D newPos;
				newPos.x = Nextpos.x;
				newPos.y = Nextpos.y;
				newPos.z = Nextpos.z;
				skeletonPoints.push_back(newPos);
	      
				FlagOnSkeleton[idx] = true;
			} // end if !FlagOnSkeleton 
		} // end while steamSteps
	   
		streamSteps = 0;
		idxSeeds--;
	}	//end while idxSeeds (main while loop)

	delete[] FlagOnSkeleton;
	seeds.clear();

	this->cleanUpMemory();

	if(debug)	//write out the vector field
	{
		FILE *fileout;
		if ((fileout = fopen("out.skel","w")) != NULL)
		{
			for (int k = 0; k < (int)skeletonPoints.size(); k++) 
			{
				fprintf(fileout,"%1.1f %1.1f %1.1f %d\n", skeletonPoints.at(k).x, skeletonPoints.at(k).y, skeletonPoints.at(k).z, 1);
			}
			fclose(fileout);
		}
		std::cerr << "Number of skeleton points = " << (int)skeletonPoints.size() << std::endl;
	}

	//Write out skeleton file:
	if(debug)
	{
		vtkFileHandler * fhdl = new vtkFileHandler();
		fhdl->SetNodes(&skeletonPoints);
		fhdl->Write("SkeletonPoints.vtk");
		delete fhdl;
	}//end if debug

	//getchar();

	return true;
}

void IntegratedSkeleton::rk2(float x, float y, float z, int sizx, int sizy, int sizz, float steps, Vector3D  *Force_ini, VoxelPosition *nextPos)
{
	Vector3D OutForce;
	OutForce=interpolation(x,y,z,sizx,sizy,sizz,Force_ini);
	nextPos->x = x + OutForce.x * steps;
	nextPos->y = y + OutForce.y * steps;
	nextPos->z = z + OutForce.z * steps;
}


//------------------------------------------------------------------------//
// Compute iso-gray surface principle curvature using Xiao Liang's 
// Direct Method without using Matrix Rotation, Transpose,Multiplication,etc.
// The Method: the principle curvature of the isosurface can be computed  
// directly from the first and second derivatives of the images.
// According to the Xiaoliang's Direct formula derivation
//------------------------------------------------------------------------//
bool IntegratedSkeleton::XiaoLiangComputeIsoGraySurfaceCurvature()
{ 
  if (debug)
	 std::cout << "This is xiaoliang's iso-curvature computation!" << std::endl;

  if(!Iu || !Iv || !Iw || !fc || !m_inputImage)
		return false;

	int L = sizeX; //x
	int M = sizeY; //y
	int N = sizeZ; //z

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

	//Init partial derivs to zeros:
	for(int i=0; i<numPix; ++i)
	{
		Iuu[i]=0;
		Ivv[i]=0;
		Iww[i]=0;
		Iuv[i]=0;
		Iuw[i]=0;
		Ivw[i]=0;
	}

	if(debug)
		std::cout << "Computing curvature" << std::endl;

	this->partialDerivative1(Iu, Iuu, 1, L, M, N);
	this->partialDerivative1(Iv, Ivv, 2, L, M, N);
	this->partialDerivative1(Iw, Iww, 3, L, M, N);
	this->partialDerivative1(Iu, Iuv, 2, L, M, N);
	this->partialDerivative1(Iu, Iuw, 3, L, M, N);
	this->partialDerivative1(Iv, Ivw, 3, L, M, N);

	//Allocate memory for curvature:
	if(curv)
	{
		delete[] curv;
		curv = NULL;
	}
	curv = new float[numPix];
	if(!curv)
		return false;

	//Init to zeros
	for (long k = 0; k < numPix; ++k)
		curv[k] = 0;

	int border = 2;	//Why does this need to by 2?
	for (int k = border; k < N-border; k++)
    {
		for (int j = border; j < M-border; j++)
		{
			for (int i = border; i < L-border; i++)
			{
				long idx = k*L*M + j*L + i;
				if (fc[idx] != INTERIOR)
					continue;
                
				double term1 = Iu[idx]*Iu[idx]*(Ivv[idx]*Iww[idx]-Ivw[idx]*Ivw[idx])
					+ 2*Iv[idx]*Iw[idx]*(Iuw[idx]*Iuv[idx]-Iuu[idx]*Ivw[idx]);

				double term2 = Iv[idx]*Iv[idx]*(Iuu[idx]*Iww[idx]-Iuw[idx]*Iuw[idx])
					+ 2*Iu[idx]*Iw[idx]*(Ivw[idx]*Iuv[idx]-Ivv[idx]*Iuw[idx]); 

				double term3 = Iw[idx]*Iw[idx]*(Iuu[idx]*Ivv[idx]-Iuv[idx]*Iuv[idx])
					+ 2*Iu[idx]*Iv[idx]*(Iuw[idx]*Ivw[idx]-Iww[idx]*Iuv[idx]); 
				
				double GaussianCurvature = term1 + term2 + term3;

				double MeanCurevature = Iu[idx]*Iu[idx]*(Ivv[idx]+Iww[idx]) -2*Iv[idx]*Iw[idx]*Ivw[idx]
				    + Iv[idx]*Iv[idx]*(Iuu[idx]+Iww[idx]) -2*Iu[idx]*Iw[idx]*Iuw[idx]
					+ Iw[idx]*Iw[idx]*(Iuu[idx]+Ivv[idx]) -2*Iu[idx]*Iv[idx]*Iuv[idx];

				double gLength = (Iu[idx]*Iu[idx] + Iv[idx]*Iv[idx] + Iw[idx]*Iw[idx]);

				if(gLength !=0 )
				{   GaussianCurvature /= (gLength*gLength);
				    MeanCurevature /= 2*pow(gLength,1.5);
			    // Gaussian curveture = k1*k2, is the multiplication of the first and second principle curvature
				// Mean curvature = (k1+k2)/2, is the mean of the  first and second principle curvature
				// thus, we can compute the first curvature from the Gaussian curveture and  Mean curvature

					double delta = MeanCurevature*MeanCurevature - GaussianCurvature;
                    if (delta >= 0)// the first principle curvature
					{    //curv[idx] = MeanCurevature + sqrt(delta);
						 double k1 = abs(MeanCurevature + sqrt(delta)); 
						 double k2 = abs(MeanCurevature - sqrt(delta)); 
                         if (k1 > k2)// the first principle curvature
						   curv[idx] = k1;
						 else 
                           curv[idx] = k2;
                          
					}
					else 
					{
						curv[idx] = MeanCurevature;
						if(debug)
						  std::cout << " Some Points with Mean curvature are used!" << std::endl;
					}
				}
				else
				{
					curv[idx]=0;
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

	if(debug)
	{
		FILE *fileout;
		if ((fileout = fopen("out.curv","w")) != NULL)
		{
			for (int k = 0; k < N; k++) {
				for (int j = 0; j < M; j++) {
					for (int i = 0; i < L; i++) 
					{
						long idx = k*L*M + j*L + i;
						fprintf(fileout, "%d %d %d %f\n", i, j, k, curv[idx]);
					} 
				} 
			} //end triple for loops
			fclose(fileout);
		}
	}

	return true;
}

bool IntegratedSkeleton:: createGradientVectorField(float sigma)
{
	if(!m_inputImage)
		return false;
	
	//Allocate memory for gradient field:
	this->clearGradientVectorField();
	Iu = new float[numPix];
	Iv = new float[numPix];
	Iw = new float[numPix];
	if(fc)
		delete[] fc;
	fc = new unsigned char[numPix];	//This keeps track of interior/exterior and surface voxels
	if(!Iu || !Iv || !Iw || !fc)
	{
		if(debug) std::cerr << "Unable to allocate memory for gradients" << std::endl;
		return false;
	}

	//Init variables to zeros:
	for(int i=0; i<numPix; ++i)
	{
		Iu[i] = 0;
		Iv[i] = 0;
		Iw[i] = 0;
		fc[i] = 0;
	}
    
	typedef itk::Image< float, 3 > FImageType;
	typedef itk::RecursiveGaussianImageFilter<ImageType, FImageType>  GaussianFilterType;
	GaussianFilterType::Pointer ga = GaussianFilterType::New();
	GaussianFilterType::Pointer gb = GaussianFilterType::New();
	GaussianFilterType::Pointer gc = GaussianFilterType::New();
	ga->SetDirection(0);
	gb->SetDirection(1);
	gc->SetDirection(2);
	ga->SetSigma(sigma);
	gb->SetSigma(sigma);
	gc->SetSigma(sigma);
	ga->SetFirstOrder();
	gb->SetFirstOrder();
	gc->SetFirstOrder();
	ga->SetInput(m_inputImage);
	gb->SetInput(m_inputImage);
	gc->SetInput(m_inputImage);

	FImageType::Pointer IuImage;
	FImageType::Pointer IvImage;
	FImageType::Pointer IwImage;

	ga->Update();
	IuImage = ga->GetOutput();
	gb->Update();
    IvImage = gb->GetOutput();
	gc->Update();
    IwImage = gc->GetOutput();
    /*
	if (debug)
	{	typedef itk::ImageFileWriter<FImageType> WriterType;
	    WriterType::Pointer writer = WriterType::New();
		writer->SetInput(IuImage);
		writer->SetFileName("Iu.mhd");
		writer->Update();
        writer->SetInput(IvImage);
		writer->SetFileName("Iv.mhd");
		writer->Update();
		writer->SetInput(IwImage);
		writer->SetFileName("Iw.mhd");
		writer->Update();
	}
    */

	itk::ImageRegionIterator< FImageType > itr0( IuImage, IuImage->GetLargestPossibleRegion() );

    long idx = 0;
	for(itr0.GoToBegin(); !itr0.IsAtEnd(); ++itr0)
	{
		Iu[idx++]=itr0.Get();
	}
    
	itk::ImageRegionIterator< FImageType > itr1( IvImage, IvImage->GetLargestPossibleRegion() );
	idx = 0;
	for(itr1.GoToBegin(); !itr1.IsAtEnd(); ++itr1)
	{
		Iv[idx++]=itr1.Get();
	}

	itk::ImageRegionIterator< FImageType > itr2( IwImage, IwImage->GetLargestPossibleRegion() );
	idx = 0;
	for(itr2.GoToBegin(); !itr2.IsAtEnd(); ++itr2)
	{
		Iw[idx++]=itr2.Get();
	}

    return true;

}

bool IntegratedSkeleton::XiaoLComputeSkeletonPoints(void)
{
    skeletonPoints.clear();
	if(!Iu || !Iv || !Iw || !fc || !m_inputImage)
		return false;

	int L = sizeX; //x
	int M = sizeY; //y
	int N = sizeZ; //z

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

	//Init partial derivs to zeros:
	for(int i=0; i<numPix; ++i)
	{
		Iuu[i]=0;
		Ivv[i]=0;
		Iww[i]=0;
		Iuv[i]=0;
		Iuw[i]=0;
		Ivw[i]=0;
	}

	if(debug)
		std::cout << "Computing curvature" << std::endl;

	this->partialDerivative1(Iu, Iuu, 1, L, M, N);
	this->partialDerivative1(Iv, Ivv, 2, L, M, N);
	this->partialDerivative1(Iw, Iww, 3, L, M, N);
	this->partialDerivative1(Iu, Iuv, 2, L, M, N);
	this->partialDerivative1(Iu, Iuw, 3, L, M, N);
	this->partialDerivative1(Iv, Ivw, 3, L, M, N);

	int border = 2;	//Why does this need to by 2?
    //typedef float  Precision;
	for (int k = border; k < N-border; k++)
    {
		for (int j = border; j < M-border; j++)
		{
			for (int i = border; i < L-border; i++)
			{
				long idx = k*L*M + j*L + i;
				if (fc[idx] != INTERIOR)
					continue;
                
			    vnl_matrix<float> H(3,3);
                H(0,0) = Iuu[idx];
                H(0,1) = H(1,0) = Iuv[idx];
                H(0,2) = H(2,0) = Iuw[idx];
                H(1,1) = Ivv[idx];
                H(1,2) = H(2,1) = Ivw[idx];
                H(2,2) = Iww[idx];
                vnl_symmetric_eigensystem <float> ES(H);
                vnl_vector<float> ev(3);
                ev[0] = ES.get_eigenvalue(0); 
                ev[1] = ES.get_eigenvalue(1); 
                ev[2] = ES.get_eigenvalue(2);
                int First=0;
				int Second=1;
				int Third=2;
				if ( ev[0] > ev[1]  ) 
				{
					std::swap(ev[0], ev[1]);
					//First =1;Second=0;
					std::swap(First,Second);

				}
                if ( ev[1] > ev[2]  ) 
				{
					std::swap(ev[1], ev[2]);
                    //Second=2;Third=1;
					std::swap(Second,Third);
				}
                if ( ev[0] > ev[1]  )
				{
	                std::swap(ev[0], ev[1]);
					//First =1;Second=0;
					std::swap(First,Second);
				}		
				vnl_vector<float> EVecFirst(3);
                vnl_vector<float> EVecSecond(3);
				//vnl_vector<float> EVecThird(3);

				EVecFirst = ES.get_eigenvector(First);
				EVecSecond = ES.get_eigenvector(Second);
                
				// using the eigenvalue and eigenvector to select the ridge points
				if (ev[1] < 0 )
				{
				  float SumofInnerProduct1;
                  float SumofInnerProduct2;
				  SumofInnerProduct1 = fabs(EVecFirst.get(0)*Iu[idx]+EVecFirst.get(1)*Iv[idx]+EVecFirst.get(2)*Iw[idx]);
                  SumofInnerProduct2 = fabs(EVecSecond.get(0)*Iu[idx]+EVecSecond.get(1)*Iv[idx]+EVecSecond.get(2)*Iw[idx]);
				  if (SumofInnerProduct1 < 0.01 && SumofInnerProduct1 < 0.01)
				  {
				  fPoint3D newPos;
				  newPos.x = i;
				  newPos.y = j;
				  newPos.z = k;
				  skeletonPoints.push_back(newPos);
				 }
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

	if(debug)	//write out the vector field
	{
		FILE *fileout;
		if ((fileout = fopen("out.skel","w")) != NULL)
		{
			for (int k = 0; k < (int)skeletonPoints.size(); k++) 
			{
				fprintf(fileout,"%1.1f %1.1f %1.1f %d\n", skeletonPoints.at(k).x, skeletonPoints.at(k).y, skeletonPoints.at(k).z, 1);
			}
			fclose(fileout);
		}
		std::cerr << "Number of skeleton points = " << (int)skeletonPoints.size() << std::endl;
	}

	//Write out skeleton file:
	if(debug)
	{
		vtkFileHandler * fhdl = new vtkFileHandler();
		fhdl->SetNodes(&skeletonPoints);
		fhdl->Write("CriticalPoints.vtk");
		delete fhdl;
	}//end if debug

	//getchar();

	return true;

}

bool IntegratedSkeleton::MovingSkeletonPointsAlongGVF(int Step)
{
  // first we normalize the GVF

  std::vector<fPoint3D> seeds = skeletonPoints;

  skeletonPoints.clear();

  int SeedsNumber = (int )seeds.size();
  

  fPoint3D Startpos;
  long slsz = (this->sizeX)*(this->sizeY);

  int idxSeeds = SeedsNumber-1;

  while(idxSeeds >= 0)
	{
		Startpos.x = seeds[idxSeeds].x; 
		Startpos.y = seeds[idxSeeds].y; 
		Startpos.z = seeds[idxSeeds].z;
		
        for (int i=0;i<Step;i++)
		{   
		long idx = (int)Startpos.z * slsz + (int)Startpos.y *sizeX + (int)Startpos.x;
		Startpos.x = Startpos.x + this->force[idx].x;
		Startpos.y = Startpos.y + this->force[idx].y;
        Startpos.z = Startpos.z + this->force[idx].z;
		}

		skeletonPoints.push_back(Startpos);

		idxSeeds--;

	}

  	//Write out skeleton file:
	if(debug)
	{
		vtkFileHandler * fhdl = new vtkFileHandler();
		fhdl->SetNodes(&skeletonPoints);
		fhdl->Write("SkeletonPointsMovedByGVF.vtk");
		delete fhdl;
	}//end if debug

	return true;
}

}