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

// StrEle.h: interface for the StrEle class.
//
//////////////////////////////////////////////////////////////////////

#ifndef StrEle_h
#define StrEle_h

#include <vector>
//#include <cmatrix>

//typedef techsoft::matrix<double> Matrix;
//typedef std::valarray<double> Vector;

class StrEle  
{
public:
	
	StrEle();
	StrEle(C3DImage * data3D, const char * ele_shape, int diameter_xy, int diameter_z);
	virtual ~StrEle();

	void NormalizeGaussianWeight();
	//unsigned char LineFilterNorm(float sigma_f, Matrix * Hessian, Vector * d);
	unsigned char LineFilterNorm(float sigma_f);
	unsigned char LineFilter();
	void SetWeight(int s, int r, int c, float weight);
	void Fill(unsigned char color);
	
	int GetXYBound();
	int GetZBound();
	unsigned char Median();
	unsigned char Min();
	unsigned char Max();
	unsigned char Convolve();
	void AddGaussianWeight(float sigma_x, float sigma_y=0.0, float sigma_z=0.0);
	void SetCenter(int s, int r, int c);
	void SetXYWidth(int new_width) { *this = StrEle(data, shape.c_str(), new_width, d_z); }
	void SetZWidth(int new_width) {	*this = StrEle(data, shape.c_str(), d_xy, new_width); }
	void Set_XY_And_Z_Width(int xy, int z) { *this = StrEle(data, shape.c_str(), xy, z); }
	int ** GetOffsets() { return offsets2; }
	float * GetWeights() { return weights; }
	int GetSize() { return size; }
	bool WithinImagePadding();

private:
	void CopyVectorToPointers();

	void CreateSphere();
	void CreateCube();
	void DisplayVector(std::vector< std::vector<int> > a_vector, std::ostream & out);

	C3DImage * data;
	int d_xy;
	int d_z;
	std::string shape;
	int center_s;	// center slice position
	int center_r;	// center row position
	int center_c;	// center column position
	vector< vector<int> > offsets;
	int **offsets2;
	float *weights;
	int size;
};

#endif

 
