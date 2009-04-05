#ifndef Branches_H
#define Branches_H


class C3DImage;
class CVessel;
class CPoint;
class CIntersectionPoint;

typedef enum TopEnd{Top = 0, End = 1};



class CBranches
{
public:
	CBranches(void);
	~CBranches(void);	
	void Search();		
private:
	bool ProcessTip(CVessel* CurrentVessel, TopEnd Flag, int VessID);
	void Calc3DRotMatrix(double d_x, double d_y, double d_z);
	bool DetectBranch(int lng, int threshold, int VessID, TopEndMiddleNeither Flag, int& targetID);
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
