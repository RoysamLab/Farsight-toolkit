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

//#pragma warning(disable:4786)  // disable STL-related warnings
//#pragma warning(disable:4710)  // compiler didn't inline function that got selected for inline expansion
//#pragma warning(disable:4514)  // The optimizer removed an inline function that is not called
//#pragma warning(disable:4702)  // unreachable STLport code

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <list>
#include <queue>
#include <sstream>


#include "CONSTANTS.h"
#include "Mytypes.h"
#include "Config.h"
#include "Cpoints.h"
#include "Dllist.h"
#include "Cimage.h"
#include "C3dimage.h"
#include "Cvessel.h"
#include "Cvector.h"
#include "Template.h"
#include "Extern.h"
#include "StrEle.h"
#include "Branches.h"



//using namespace std;

CBranches::CBranches(void)
{
	int i,j,k;
	NumberOfBranches = 0;	
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			Rot_3D[i][j] = 0.0;


	//

	//Create the traces image. This image is black everywhere except for the centerline points (white)	
	traces_image = new int**[The3DImage->m_iSlices];
	for(i=0;i<The3DImage->m_iSlices;i++)
		traces_image[i] = new int*[The3DImage->m_iRows];
	for(i=0;i<The3DImage->m_iSlices;i++)
	{
		for(j=0;j<The3DImage->m_iRows;j++)
		{
			traces_image[i][j] = new int[The3DImage->m_iCols];
			for(k=0; k<The3DImage->m_iCols; k++)
				traces_image[i][j][k] = 0;
		}
	}
	

	for(i=0;i<gTheVessels.m_iNumOfElements;i++)
	{
		if (gTheVessels.m_apData[i])
			{
				CLNode<CPoint>* temp = gTheVessels.m_apData[i]->m_Center.head;
				while (temp)
				{
					traces_image[temp->data->m_iZ][temp->data->m_iY][temp->data->m_iX] = i+1;					
					temp = temp->after;
				}
			}
	}		

	int index = 0;
	int d1 = 21;
	int d2 = 10;
	int d3 = 10;
	for(i=0;i<d1;i++)
	{
		for(j=-d2;j<=d2;j++)
		{
			for(k=-d3;k<=d3;k++)
			{
				VOL1[0][index] = i;
				VOL1[1][index] = j;
				VOL1[2][index] = k;

				VOL2[0][index] = 0;
				VOL2[1][index] = 0;
				VOL2[2][index] = 0;

				tmp1[k+d3][j+d2][i] = 0;
				tmp1[k+d3][j+d2][i] = 0;
				tmp1[k+d3][j+d2][i] = 0;
				
				index++;
			}
		}
	}	
}

CBranches::~CBranches(void)
{
}

void CBranches::Search()
{
	cout << "\tDetecting Branches ... ";
	for (int i = 0; i < gTheVessels.m_iNumOfElements; i++)
	{		
		if (gTheVessels.m_apData[i])
		{
			//First, ignore the vessel if it is too short (less than 10 voxels)
			if(gTheVessels.m_apData[i]->GetLength() < 10)
				continue;
			//Then:
			//1- Process the Top of the neurite
			bool procTip = ProcessTip(gTheVessels.m_apData[i],Top,i+1);
			//By Yousef on (3-24-2009)
			//A segment should not be connected from both sides
			//so don't process the end unless you have not found a point at the top
			//2- Process the End of the neurite
			if(!procTip)
				ProcessTip(gTheVessels.m_apData[i],End,i+1);
		}
	}
	cout<<NumberOfBranches<<" branches detected"<<endl;
}

bool CBranches::ProcessTip(CVessel* CurrentVessel, TopEnd Flag, int VessID)
{	
	if(CurrentVessel->GetLength() < 15)
		return false;
	int i,j,k;
	//First of all,,,, Reset the vars
	//ResetVars();
	//1. Select the Tip point and the anchor point
	if(Flag == Top) //In case we work on the first point (TOP)
	{
		//Create a temp pointer to the linked list of the center points of the current vessle. The pointer points to the head of the vessel
		CLNode<CPoint>* temp = CurrentVessel->m_Center.head;
		//Set the Tip point to be the first point
		TipPoint    = new CPoint(temp->data->m_iX,temp->data->m_iY,temp->data->m_iZ,temp->data->m_iHDir,temp->data->m_iVDir,temp->data->m_iValue);
		//move the pointer 4 points ahead
		for(i=1; i<=4; i++)
			temp = temp->after;
		//set the anchor point to be the 5th point from top
		AnchorPoint = new CPoint(temp->data->m_iX,temp->data->m_iY,temp->data->m_iZ,temp->data->m_iHDir,temp->data->m_iVDir,temp->data->m_iValue);
	}
	else //In case we work in the last point (END)
	{
		/////////////////////////////////////////////////////////////
		//Here we do the same as shown above but for the last point// 
		/////////////////////////////////////////////////////////////

		CLNode<CPoint>* temp = CurrentVessel->m_Center.tail;
		TipPoint    = new CPoint(temp->data->m_iX,temp->data->m_iY,temp->data->m_iZ,temp->data->m_iHDir,temp->data->m_iVDir,temp->data->m_iValue);
		for(i=1; i<=4; i++)
			temp = temp->before;
		AnchorPoint = new CPoint(temp->data->m_iX,temp->data->m_iY,temp->data->m_iZ,temp->data->m_iHDir,temp->data->m_iVDir,temp->data->m_iValue);
	}

	//2. Compute the entries of the rotation matrix
	//Find the differences vectors between the tip and the anchor points
	int d_x = TipPoint->m_iX - AnchorPoint->m_iX;
	if(d_x == 0)
		d_x++;
	int d_y = TipPoint->m_iY - AnchorPoint->m_iY;
	int d_z = TipPoint->m_iZ - AnchorPoint->m_iZ;	

	//Call The function that computes the 3D rotation matrix (Rot_3D)
	Calc3DRotMatrix((double)d_x,(double)d_y,(double)d_z);

	//3. Transform VOL1, which is oriented along the x-axis
	//3.1 Find the inverse of Rot_3D, which is also its transpose
	double Rot_3D_tr[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
			Rot_3D_tr[i][j] = Rot_3D[j][i];
	}

	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
			Rot_3D[i][j] = Rot_3D_tr[i][j];
	}

	//3.2 Rotate VOL1 (Multiply Rot_3D_tr by VOL1 and save the results in VOL2)
	double v = 0;
	for(j=0;j<9261;j++)
	{
		for(i=0;i<3;i++)
		{
			for(k=0;k<3;k++)
			{
				v += Rot_3D[i][k] * VOL1[k][j];
			}
			VOL2[i][j] = (int)(v);
			v = 0;
		}		
	}
	//3.3 Translate  VOL2 from the origen to the original position of the tip point
	for(i=0; i<9261; i++)
	{
		VOL2[0][i] += TipPoint->m_iX;
		VOL2[0][i] = min(VOL2[0][i],The3DImage->m_iCols-1);
		VOL2[0][i] = max(VOL2[0][i],0);
		VOL2[1][i] += TipPoint->m_iY;
		VOL2[1][i] = min(VOL2[1][i],The3DImage->m_iRows-1);
		VOL2[1][i] = max(VOL2[1][i],0);
		VOL2[2][i] += TipPoint->m_iZ;
		VOL2[2][i] = min(VOL2[2][i],The3DImage->m_iSlices-1);
		VOL2[2][i] = max(VOL2[2][i],0);
	}

	
	//4. Fill the two temporary 3D templates tmp1 and tmp2 using corresponding voxels in The3DImage and tmpImage respectively
	//Also, see if tmp2 (the template of traces) contains any other trace segment
	int r1,r2,c1,c2,s1,s2;
	int hit = 0;
	for(i=0;i<9261;i++)
	{
		r1 = VOL1[1][i]+10;
		r2 = VOL2[1][i];
		c1 = VOL1[0][i];
		c2 = VOL2[0][i];
		s1 = VOL1[2][i]+10;
		s2 = VOL2[2][i];
		tmp1[s1][r1][c1] = (int) (The3DImage->data[s2][r2][c2]);
		tmp2[s1][r1][c1] = traces_image[s2][r2][c2];
		if(tmp2[s1][r1][c1] !=0 && tmp2[s1][r1][c1]!=VessID)
			hit = 1;
	}

	//If no other trace segment was found in tmp2, then there is no need to go further
	if(hit == 0)
		return false;

	int lng = 20;
	///////////////////////////////////////////////////////////////////////////////////
	//By Yousef Nov. 23 2006
	//Try this
	int diff = 0;
	//bool br = false;
	for(i=0;i<21;i++)
	{
		if(i%2 == 0)
			diff = i/2;
		for(j=10-diff; j<=10+diff; j++)
			for(k=10-diff; k<=10+diff;k++)
			{
				if(tmp2[k][j][i]!=0 && tmp2[k][j][i]!=VessID)
				{
					if(i<15)
						lng = i+2;
					else
						lng = 20;

						k = 10+diff;
						j = 10+diff;
						i = 20;
					break;
				}
			}
	}

	//try this
	//lng = 10;
	///////////////////////////////////////////////////////////////////////////////////
	//5. Detect the branch
		
	double threshold = 1.;
	int target = 0;
	bool brs;
	lng = 20;
	if(Flag == Top)
	{
		brs = DetectBranch(lng,threshold,VessID, OnTop, target);
	}
	else
	{
		brs = DetectBranch(lng,threshold,VessID, OnEnd, target);
	}

	if(!brs)
	{
		//there is no branch detected close to the current tip point
		return false;
	}
	
	//6. Add Intersection point
	//In (DetectBranch), the points extended from the tip point to the detected branch point were all added to the current vessel
	//Now, add the detected branch point to the set of intersection points.
	IntrPoint = new CIntersectionPoint(*BrPnt);	

	gIntersectionPoints.Add(IntrPoint);

	return true;
}


void CBranches::Calc3DRotMatrix(double d_x, double d_y, double d_z)
{
	int i, j, k;
	//1. Fill the entries of the three rotation matrices
	double Theta_X = 0;
	double Rot_X[3][3] = 
	{
		{1, 0, 0}, {0, cos(Theta_X), sin(Theta_X)}, {0,-sin(Theta_X), cos(Theta_X)}
	};
	double x_y = d_y/d_x; 
	double Theta_Z = atan(x_y);
	double Rot_Z[3][3] = 
	{
		{cos(Theta_Z),sin(Theta_Z),0},{-sin(Theta_Z),cos(Theta_Z),0},{0,0,1}
	};
	double Theta_Y = 0;
	//double pxz[3] = {{Rot_Z[0][0]*d_x + Rot_Z[0][1]*d_y + Rot_Z[0][2]*d_z}, {Rot_Z[1][0]*d_x + Rot_Z[1][1]*d_y + Rot_Z[1][2]*d_z}, {Rot_Z[2][0]*d_x + Rot_Z[2][1]*d_y + Rot_Z[2][2]*d_z} };
	double pxz[3] = {(Rot_Z[0][0]*d_x + Rot_Z[0][1]*d_y + Rot_Z[0][2]*d_z), (Rot_Z[1][0]*d_x + Rot_Z[1][1]*d_y + Rot_Z[1][2]*d_z), (Rot_Z[2][0]*d_x + Rot_Z[2][1]*d_y + Rot_Z[2][2]*d_z) };
	if(pxz[0]<0)
	{
		if(pxz[2]>=0)
			Theta_Y = 180 - (abs(atan(pxz[2]/pxz[0]))*180/3.142);
		else
			Theta_Y = 180 + (abs(atan(pxz[2]/pxz[0]))*180/3.142);
	}
	else
		Theta_Y = (atan(pxz[2]/pxz[0])*180/3.142);

	Theta_Y = -Theta_Y;
	Theta_Y *= (3.142/180);

	double Rot_Y[3][3] = 
	{
		{cos(Theta_Y),0,-sin(Theta_Y)}, {0,1,0}, {sin(Theta_Y),0,cos(Theta_Y)}
	};

	//2. Create the 3D rotation matrix
	double Rot_XY[3][3] = {{0,0,0},{0,0,0},{0,0,0}};	
	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			for(k=0;k<3;k++)
			{
				Rot_XY[i][j] += (Rot_X[i][k] * Rot_Y[k][j]);
			}
		}
	}
	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			Rot_3D[i][j] = 0;
			for(k=0;k<3;k++)
			{
				Rot_3D[i][j] += (Rot_XY[i][k] * Rot_Z[k][j]);
			}
		}
	}
}


bool CBranches::DetectBranch(int lng, double threshold, int SegID, TopEndMiddleNeither Flag, int& target)
{	
	double F_D_L[256] = {0.000017,0.000034,0.000035,0.000036,0.000038,0.000079,0.000200,0.000201,0.000442,0.000444,0.000685,0.000487,0.000369,0.000490,0.000452,0.000813,0.000415,0.000417,0.000499,0.000700,0.000662,0.000704,0.000946,0.000588,0.000870,0.000552,0.000754,0.000756,0.000519,0.000601,0.000684,0.000486,0.000489,0.000851,0.000733,0.000936,0.000819,0.000861,0.000904,0.000867,0.001349,0.001153,0.001036,0.001358,0.001402,0.001045,0.001168,0.001172,0.001295,0.001259,0.001742,0.001266,0.001549,0.001713,0.001757,0.001880,0.002204,0.001529,0.002212,0.001937,0.002261,0.002185,0.001910,0.002554,0.002279,0.002244,0.002568,0.002653,0.002419,0.002823,0.002389,0.002474,0.002640,0.002086,0.002571,0.002856,0.002782,0.002628,0.002954,0.002641,0.002567,0.002773,0.002460,0.002387,0.005189,0.002440,0.002248,0.002175,0.003340,0.002549,0.002517,0.002125,0.002971,0.002660,0.002827,0.002676,0.002644,0.003132,0.002981,0.002710,0.003318,0.002528,0.002497,0.002746,0.002556,0.003204,0.003494,0.002865,0.002356,0.003085,0.002376,0.003226,0.002398,0.002449,0.002939,0.002631,0.002802,0.003014,0.002786,0.002758,0.002691,0.003302,0.003195,0.002849,0.003341,0.002396,0.003248,0.003222,0.002997,0.003171,0.003145,0.002840,0.003095,0.003230,0.002646,0.002982,0.003717,0.003214,0.003390,0.003207,0.002945,0.003361,0.003099,0.003237,0.003615,0.003993,0.004091,0.004190,0.004010,0.002952,0.003451,0.003751,0.003532,0.003194,0.003375,0.003756,0.003379,0.003522,0.003625,0.003528,0.003792,0.003856,0.004160,0.003825,0.003452,0.003637,0.004223,0.004769,0.004037,0.004385,0.004652,0.004122,0.004231,0.004141,0.004691,0.004243,0.004234,0.004546,0.004379,0.004173,0.004766,0.004361,0.004596,0.004911,0.004149,0.004266,0.005023,0.004423,0.004503,0.004743,0.004744,0.004826,0.004390,0.005272,0.005077,0.004963,0.004970,0.004978,0.005107,0.004997,0.005048,0.005260,0.005034,0.005368,0.005744,0.004843,0.005541,0.005681,0.005544,0.005168,0.005553,0.005220,0.006008,0.005680,0.005593,0.005827,0.005745,0.005745,0.005867,0.006192,0.005480,0.006929,0.006304,0.006042,0.006143,0.006886,0.006716,0.007069,0.006788,0.006791,0.006800,0.007575,0.007437,0.007586,0.007704,0.007671,0.008567,0.008195,0.008836,0.008412,0.008684,0.009414,0.010167,0.010268,0.010321,0.010493,0.010794,0.011277,0.012401,0.013319,0.013639,0.015129,0.016891,0.024860,0.021752,0.022125};
	double B_D_L[256] = {0.169375,0.079679,0.060477,0.099617,0.214613,0.095008,0.052177,0.036614,0.027409,0.021641,0.017268,0.013791,0.011469,0.009418,0.007916,0.006722,0.005876,0.004923,0.004444,0.003885,0.003505,0.003071,0.002911,0.002601,0.002267,0.002154,0.001964,0.001902,0.001798,0.001611,0.001576,0.001437,0.001459,0.001293,0.001274,0.001197,0.000978,0.001041,0.001040,0.000995,0.000931,0.000882,0.000001,0.000827,0.000717,0.000732,0.000719,0.000655,0.000637,0.000593,0.000549,0.000554,0.000532,0.000525,0.000475,0.000490,0.000431,0.000393,0.000402,0.000444,0.000365,0.000358,0.000378,0.000345,0.000332,0.000292,0.000268,0.000294,0.000277,0.000228,0.000270,0.000261,0.000224,0.000220,0.000222,0.000217,0.000200,0.000257,0.000217,0.000222,0.000182,0.000176,0.000202,0.000174,0.000000,0.000143,0.000130,0.000154,0.000141,0.000149,0.000125,0.000114,0.000134,0.000114,0.000112,0.000092,0.000083,0.000079,0.000083,0.000079,0.000099,0.000068,0.000097,0.000068,0.000079,0.000075,0.000057,0.000072,0.000059,0.000075,0.000068,0.000066,0.000064,0.000035,0.000042,0.000040,0.000048,0.000024,0.000046,0.000024,0.000042,0.000057,0.000024,0.000031,0.000044,0.000029,0.000020,0.000022,0.000018,0.000031,0.000026,0.000018,0.000033,0.000044,0.000018,0.000013,0.000018,0.000024,0.000020,0.000020,0.000018,0.000020,0.000013,0.000007,0.000020,0.000022,0.000009,0.000007,0.000015,0.000011,0.000015,0.000011,0.000007,0.000013,0.000020,0.000015,0.000009,0.000015,0.000018,0.000004,0.000004,0.000002,0.000002,0.000007,0.000007,0.000004,0.000007,0.000007,0.000002,0.000000,0.000013,0.000004,0.000007,0.000000,0.000004,0.000007,0.000004,0.000002,0.000009,0.000000,0.000007,0.000004,0.000000,0.000009,0.000004,0.000004,0.000002,0.000004,0.000000,0.000002,0.000000,0.000000,0.000002,0.000000,0.000002,0.000004,0.000004,0.000000,0.000000,0.000000,0.000000,0.000002,0.000000,0.000000,0.000000,0.000000,0.000002,0.000004,0.000000,0.000000,0.000000,0.000002,0.000002,0.000000,0.000000,0.000000,0.000000,0.000002,0.000000,0.000000,0.000000,0.000002,0.000002,0.000002,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000002,0.000002,0.000002,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000002,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000002,0.000000,0.000000};
	double Resp[17] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
	double edpnts2[25/*9*/][3],step,stp;
	double theta1,theta2;
	step = 1.0/lng;	
	
	theta1 = 3.142 * 7.5/180;
	theta2 = 3.142 * 15/180;
	int i,j,k,h,hit,Gr_Mean,d1, d2, d3, d4, index;
	//threshold = 0;
	int pnt[3];
	//int spnt[3];
	//double test = tan(theta1);
	d1 = (int) (lng*tan(theta1));
	d2 = (int) (lng*tan(theta2));
	d3 = (int) (d1/sqrt(2.0));
	d4 = (int) (d2/sqrt(2.0));
	//int stpnts[9][3] = {{0,9,9},{0,9,10},{0,9,11},{0,10,9},{0,10,10},{0,10,11},{0,11,9},{0,11,10},{0,11,11}};
	int stpnts[25][3] = {{0,8,8},{0,8,9},{0,8,10},{0,8,11},{0,8,12},{0,9,8},{0,9,9},{0,9,10},{0,9,11},{0,9,12},{0,10,8},{0,10,9},{0,10,10},{0,10,11},{0,10,12},{0,11,8},{0,11,9},{0,11,10},{0,11,11},{0,11,12},{0,12,8},{0,12,9},{0,12,10},{0,12,11},{0,12,12}};
	//int edpnts[9][3] = {{lng,9,9},{lng,9,10},{lng,9,11},{lng,10,9},{lng,10,10},{lng,10,11},{lng,11,9},{lng,11,10},{lng,11,11}};
	int edpnts[25][3] = {{lng,8,8},{lng,8,9},{lng,8,10},{lng,8,11},{lng,8,12},{lng,9,8},{lng,9,9},{lng,9,10},{lng,9,11},{lng,9,12},{lng,10,8},{lng,10,9},{lng,10,10},{lng,10,11},{lng,10,12},{lng,11,8},{lng,11,9},{lng,11,10},{lng,11,11},{lng,11,12},{lng,12,8},{lng,12,9},{lng,12,10},{lng,12,11},{lng,12,12}};
	int changes[17][2] = {{0,0},{d1,0},{d2,0},{-d1,0},{-d2,0},{0,d1},{0,d2},{0,-d1},{0,-d2},{d3,d3},{d4,d4},{-d3,-d3},{-d4,-d4},{d3,-d3},{d4,-d4},{-d3,d3},{-d4,d4}};	
		

	//Calculate the 17 different responses (we have 17 different orientations of the tested volume
	for(i=0;i<17;i++)
	{
		Gr_Mean = 0;
		index = 0;
		for(j=0;j<25/*9*/;j++)
		{
			for(k=0;k<3;k++)
			{
				edpnts2[j][k] = edpnts[j][k];
			}
			edpnts2[j][1] += changes[i][0];
			edpnts2[j][2] += changes[i][1];
		}
		
		for(k=0;k<lng;k++)
		{
			stp = k*step;
			for(j=0;j<25/*9*/;j++)
			{
				for(h=0;h<3;h++)
				{
					pnt[h] = (int)((1-stp)*stpnts[j][h] + stp*edpnts2[j][h]); 
				}
				Gr_Mean += tmp1[pnt[2]][pnt[1]][pnt[0]];
				index++;
			}			
		}	
		Gr_Mean /= index;
		Resp[i] = F_D_L[Gr_Mean] / B_D_L[Gr_Mean];		
	}

	//Get the maximum
	double MX_RSP = Resp[0];
	int MX_IND = 0;
	for(i=1;i<17;i++)
	{
		if(Resp[i] > MX_RSP)
		{
			MX_RSP = Resp[i];
			MX_IND = i;
		}
	}

	//The maximum must also be greater than the given threshold
	threshold = 0.6;
	if(MX_RSP < threshold)
		return false;	
	
	
	for(j=0;j<25/*9*/;j++)
	{
		for(k=0;k<3;k++)
		{
			edpnts2[j][k] = edpnts[j][k];
		}
		edpnts2[j][1] += changes[MX_IND][0];
		edpnts2[j][2] += changes[MX_IND][1];
	}

	hit = 0;
	index = 0;
	lng = min(lng+5,20);
	for(k=0;k<lng;k++)
	{		
		index++;
		stp = k*step;
		for(j=0;j<25/*9*/;j++)
		{		
			for(h=0;h<3;h++)
			{
				pnt[h] = (int)((1-stp)*stpnts[j][h] + stp*edpnts2[j][h]); 				
			}			
			//if(tmp2[pnt[2]][pnt[1]][pnt[0]]!=0 && tmp2[pnt[2]][pnt[1]][pnt[0]]!=SegID)
			//{				
				pnt[1] -= 10;
				pnt[2] -= 10;
				int p0 = (int)(Rot_3D[0][0]*pnt[0] + Rot_3D[0][1]*pnt[1] + Rot_3D[0][2]*pnt[2]);
				int p1 = (int)(Rot_3D[1][0]*pnt[0] + Rot_3D[1][1]*pnt[1] + Rot_3D[1][2]*pnt[2]);
				int p2 = (int)(Rot_3D[2][0]*pnt[0] + Rot_3D[2][1]*pnt[1] + Rot_3D[2][2]*pnt[2]);
				BrPnt = new CPoint(p0+TipPoint->m_iX, p1+TipPoint->m_iY, p2+TipPoint->m_iZ
, TipPoint->m_iHDir, TipPoint->m_iVDir, TipPoint->m_iValue);
				//get the target vessel
				target = traces_image[BrPnt->m_iZ][BrPnt->m_iY][BrPnt->m_iX];
				//make sure that there is a target vessel
				if(target == 0 || target == SegID)
					continue;
				CVessel* targetVess = gTheVessels.m_apData[target-1];
				//Make sure the branch is not among the first 5 points in the target segment
				CLNode<CPoint>* trgt_strt = targetVess->m_Center.head;
				int fnd = 0;
				for(i=0; i<5; i++)
				{
					if(BrPnt->m_iX == trgt_strt->data->m_iX && BrPnt->m_iY == trgt_strt->data->m_iY && BrPnt->m_iZ == trgt_strt->data->m_iZ)
					{
						fnd = 1;
						break;
					}
					trgt_strt = trgt_strt->after;
				}
				if(fnd == 1)
					return false;
				//Make sure the branch is not among the last 5 points in the target segment
				CLNode<CPoint>* trgt_end = targetVess->m_Center.tail;
				fnd = 0;
				for(i=0; i<5; i++)
				{
					if(BrPnt->m_iX == trgt_end->data->m_iX && BrPnt->m_iY == trgt_end->data->m_iY && BrPnt->m_iZ == trgt_end->data->m_iZ)
					{
						fnd = 1;
						break;
					}
					trgt_end = trgt_end->before;
				}
				if(fnd == 1)
					return false;
				
				//if you reach here, it means that we have found the branch
				NumberOfBranches++;
				hit = 1;				
				break;
			//}			
		}
		if(hit == 1)
			break;
	}		

	if(hit == 0)
		//there is no branch
		return false;
	
	lng = (int) sqrt((double)((BrPnt->m_iX - TipPoint->m_iX)*(BrPnt->m_iX - TipPoint->m_iX) + (BrPnt->m_iY - TipPoint->m_iY)*(BrPnt->m_iY - TipPoint->m_iY) + (BrPnt->m_iZ - TipPoint->m_iZ)*(BrPnt->m_iZ - TipPoint->m_iZ)));
	CPoint** tmpPoints = new CPoint*[lng];	
	step = 1.0/lng;
	for(i=0;i<lng;i++)
	{
		stp = i*step;
		pnt[0] = (int)((1-stp)* TipPoint->m_iX + stp*BrPnt->m_iX); 
		pnt[1] = (int)((1-stp)* TipPoint->m_iY + stp*BrPnt->m_iY);
		pnt[2] = (int)((1-stp)* TipPoint->m_iZ + stp*BrPnt->m_iZ);

		tmpPoints[i] = new CPoint(pnt[0],pnt[1],pnt[2],TipPoint->m_iHDir,TipPoint->m_iVDir,TipPoint->m_iValue);
		gTheVessels.m_apData[SegID-1]->ExtendVesselCenter(tmpPoints[i],Flag);
		int tID = traces_image[BrPnt->m_iZ][BrPnt->m_iY][BrPnt->m_iX];
		gTheVessels.m_apData[SegID-1]->SetParentID(tID);
		if(Flag == OnTop)
			gTheVessels.m_apData[SegID-1]->SetParentLocation(0);
		else
			gTheVessels.m_apData[SegID-1]->SetParentLocation(1);
		gTheVessels.m_apData[SegID-1]->ParentBranchPoint = new CPoint();
		gTheVessels.m_apData[SegID-1]->ParentBranchPoint->m_iX = BrPnt->m_iX;
		gTheVessels.m_apData[SegID-1]->ParentBranchPoint->m_iY = BrPnt->m_iY;
		gTheVessels.m_apData[SegID-1]->ParentBranchPoint->m_iZ = BrPnt->m_iZ;
	}	
	//tmpPoints[i] = new CPoint(BrPnt->m_iX,BrPnt->m_iY,BrPnt->m_iZ,TipPoint->m_iHDir,TipPoint->m_iVDir,TipPoint->m_iValue);
	//gTheVessels.m_apData[SegID-1]->ExtendVesselCenter(tmpPoints[i],Flag);
	return true;
}
