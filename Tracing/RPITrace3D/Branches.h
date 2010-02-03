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

#ifndef Branches_H
#define Branches_H


class C3DImage;
class CVessel;
class CPoint;
class CIntersectionPoint;

enum TopEnd{Top = 0, End = 1};



class CBranches
{
public:
	CBranches(void);
	~CBranches(void);	
	void Search();		
private:
	bool ProcessTip(CVessel* CurrentVessel, TopEnd Flag, int VessID);
	void Calc3DRotMatrix(double d_x, double d_y, double d_z);
	bool DetectBranch(int lng, double threshold, int VessID, TopEndMiddleNeither Flag, int& targetID);
	CPoint* TipPoint;
	CPoint* AnchorPoint;
	CPoint* BrPnt;
	CIntersectionPoint* IntrPoint;
	int NumberOfBranches;
	double Rot_3D[3][3];
	int VOL1[3][9261]; //The coordinates of a 21x21x21 volume that is oriented along the x-axis
	int VOL2[3][9261]; //The coordinates of the volume resulting from rotating VOL1
	int tmp1[21][21][21];
	int tmp2[21][21][21];
	int*** traces_image;		
};
#endif
