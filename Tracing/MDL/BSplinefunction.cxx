//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
//                                                                   
//   File Name:      test.cxx                       
//   Author:         Xiao liang, Yu Wang,  Sep.15,2009
//   Description:    Bspline fitting Fucntion
//  
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

#include <cstring>
#include <stdio.h>
#include <vector>
#include <fstream>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/algo/vnl_matrix_inverse.h>
//#include "BSpline.h"
#define Dimention  3;

struct  VoxelPosition
{
	float x;
	float y;
	float z;
}; 

vnl_vector<float> Q[3];



float CubicNPSpline(float t, int n)
{
  float value;
  float s = t - n;
  if (s >= -2 && s < -1)
    value = (pow((2+s),3))/6;
  else if (s >= -1 && s < 0)
    value = (4 - 6*pow(s,2) - 3*pow(s,3))/6;
  else  if (s >= 0 && s < 1)
    value = (4 - 6*pow(s,2) + 3*pow(s,3))/6;
  else  if (s >= 1 && s <= 2)  
    value = (pow((2-s),3))/6;
  else
    value = 0;
  return value;
}


 void FitNPSpline(int Nb,int datalength, VoxelPosition *datapoints)
{
  int Ns = datalength;
  vnl_vector<float> knots(Ns);
  knots(0) = 0;
  float step = static_cast<float>(Nb)/static_cast<float>(Ns-1);
  for(int i=1; i<Ns; i++)
  {
   knots(i) = knots(i-1) + step;
  }

  // Construct basis matrix
  vnl_matrix<float> Basis(Ns,Nb+5);
  for(int i=0; i<Ns; i++)
  {
    for(int j = -3; j< Nb+2; j++)
	{
	  Basis(i,j+3) = CubicNPSpline(knots(i),j+1);
	}
  }
  vnl_matrix<float> Binv = vnl_matrix_inverse<float>(Basis);
  vnl_vector<float> x(Ns); 
  vnl_vector<float> y(Ns);
  vnl_vector<float> z(Ns);
  for(int ij=0; ij<Ns; ij++)
  {
      x(ij) = datapoints[ij].x;
	  y(ij) = datapoints[ij].y;
	  z(ij) = datapoints[ij].z;
  }
  vnl_vector<float> Qx = Binv * x;
  vnl_vector<float> Qy = Binv * y;
  vnl_vector<float> Qz = Binv * z;
  Q[0] = Qx;
  Q[1] = Qy;
  Q[2] = Qz;

}


void SampleNPSpline(int NPointsSample,VoxelPosition PSample[1024], float P1, float P2)
{
  vnl_vector<float> Qx = Q[0];
  vnl_vector<float> Qy = Q[1];
  vnl_vector<float> Qz = Q[2];
  int Nb = Qx.size() - 5;
  
 // sample points with Qx and Qy
  vnl_vector<float> t(NPointsSample);
  t(0) = P1 * Nb;
  float steps = static_cast<float>((P2-P1)*Nb)/static_cast<float>(NPointsSample-1);
  for(int i=1; i<NPointsSample; i++)
  {
   t(i) = t(i-1) + steps;
  }
  vnl_vector<float> xs(NPointsSample); 
  vnl_vector<float> ys(NPointsSample); 
  vnl_vector<float> zs(NPointsSample); 
  xs.fill(0);
  ys.fill(0);
  zs.fill(0);
  for(int j=-3; j<Nb+2; j++)
  {
    for(int i=0; i<NPointsSample; i++)
	{
	 xs(i) = xs(i) + Qx(j+3) * CubicNPSpline(t(i), j+1);
     ys(i) = ys(i) + Qy(j+3) * CubicNPSpline(t(i), j+1);
	 zs(i) = zs(i) + Qz(j+3) * CubicNPSpline(t(i), j+1);
	}
  }

  for(int k =0; k<NPointsSample; k++)
  {
      PSample[k].x=xs(k);
	  PSample[k].y=ys(k);
	  PSample[k].z=zs(k);
  }
  //return PSample;
}



int main( void )
{  
   float Xdata[14]={1.,2.,3.,4.,6.,7.,9.,7.,6.,4.,5.,3.,2.,1.};
   float Ydata[14]={1.,2.,3.,5.,6.,7.,9.,7.,6.,5.,5.,3.,2.,1.};
   float Zdata[14]={1.,2.,3.,5.,6.,7.,9.,7.,6.,5.,5.,3.,2.,1.};
   int i;

   VoxelPosition *PSample;
   PSample = new VoxelPosition[14];
   for (i=0;i<14;i++)
   {PSample[i].x=Xdata[i];
    PSample[i].y=Ydata[i];
	PSample[i].z=Zdata[i];
   }
   int Nb = 14;
   //XFitPSpline(Nb,14, PSample);
   FitNPSpline(Nb,14, PSample);
   VoxelPosition *ResultSample;
   int length = 28;
   ResultSample = new VoxelPosition[length];
   // XSamplePDDSpline(length,ResultSample);
   SampleNPSpline(length,ResultSample,0,1);
   for (i=0;i<length;i++)
	   printf("%f ", ResultSample[i].x);

   printf ("\n");

   for (i=0;i<length;i++)
	   printf("%f ", ResultSample[i].y);

    printf ("\n");

   for (i=0;i<length;i++)
	   printf("%f ", ResultSample[i].z);
   delete []ResultSample;
   delete []PSample;

  
   return 1;
}
