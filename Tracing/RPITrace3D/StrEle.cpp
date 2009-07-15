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

// StrEle.cpp: implementation of the StrEle class.
// Author: Muhammad-Amri Abdul-Karim (abdulm@rpi.edu)
//
//////////////////////////////////////////////////////////////////////
//#pragma warning(disable:4786)  // disable STL-related warnings
//#pragma warning(disable:4710)  // compiler didn't inline function that got selected for inline expansion
//#pragma warning(disable:4514)  // The optimizer removed an inline function that is not called
//#pragma warning(disable:4702)  // unreachable STLport code

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <list>
#include <deque>
#include <queue>
#include <ctime>

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


using namespace std;

#define PI 3.142

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

// default constructor
StrEle::StrEle()
{
	data = NULL;
	weights = NULL;
	center_c = center_r = center_s = 0;
	d_xy = d_z = 0;
}

// default destructor
StrEle::~StrEle()
{
}

// this gaussian assumes u = 0
float Gaussian(float x, float sigma)
{
	return static_cast<float>((1.0 / (sqrt(2.0 * PI) * sigma) * exp(-(pow(static_cast<double>(x), 2.0) / (2 * pow(static_cast<double>(sigma), 2.0))))));
}

// this function creates offsets that defines a sphere
void StrEle::CreateSphere()
{

	vector<int> xyz_offset;
	xyz_offset.push_back(0);xyz_offset.push_back(0);xyz_offset.push_back(0);

	if (d_z > d_xy)
	{
		d_z = d_xy;
	}

	// calculate the z_scale (with round-off correction)
	float fz_scale = (float) d_xy / (float) d_z;
	int z_scale = Round(fz_scale);

	int radius = d_xy / 2;

	register int i, j, k;
	int rSquared = radius* radius;

	float dx, dy, dz;  
	float cx, cy,cz;
	cx = cy = cz = 0.0;

	// from middle, go z-up according to z-scale
	k = static_cast<int>(cz);
	int a = 0;
	while (k <= cz + radius)
	{
		for (j = static_cast<int>(cy) - radius; j <= static_cast<int>(cy) + radius; j++)
		{
			for (i = static_cast<int>(cx) - radius; i <= static_cast<int>(cx) + radius; i++)
			{
				dx = ((float) i - cx) * ((float) i - cx);
				dy = ((float) j - cy) * ((float) j - cy);
				dz = ((float) k - cz) * ((float) k - cz);
				if ((dx + dy + dz) <= rSquared)
				{
					xyz_offset[0] = a;
					xyz_offset[1] = j;
					xyz_offset[2] = i;
					offsets.push_back(xyz_offset);
				}
			}
		}
		a++;
		k = Round((float) a * fz_scale);
	}

	// from middle, go z-down according to z-scale
	k = static_cast<int>(cz) - z_scale;
	a = -1;
	while (k >= cz - radius)
	{
		for (j = static_cast<int>(cy) - radius; j <= cy + radius; j++)
		{
			for (i = static_cast<int>(cx) - radius; i <= cx + radius; i++)
			{
				dx = ((float) i - cx) * ((float) i - cx);
				dy = ((float) j - cy) * ((float) j - cy);
				dz = ((float) k - cz) * ((float) k - cz);
				if ((dx + dy + dz) <= rSquared)
				{
					xyz_offset[0] = a;
					xyz_offset[1] = j;
					xyz_offset[2] = i;
					offsets.push_back(xyz_offset);
				}
			}
		}
		a--;
		k = Round((float) a * fz_scale);
	}
	this->CopyVectorToPointers();
}

void StrEle::Fill(unsigned char color)
{
	for (int i = 0; i < this->GetSize(); i++)
	{
		if(data->WithinImageMargin(center_s + offsets[i][0],center_r + offsets[i][1],center_c + offsets[i][2],giMARGIN))
			this->data->data[center_s + offsets[i][0]][center_r + offsets[i][1]][center_c + offsets[i][2]] = color;
	}
}

// this function creates offsets that define a cube
void StrEle::CreateCube()
{
	vector<int> xyz_offset;

	xyz_offset.push_back(0);xyz_offset.push_back(0);xyz_offset.push_back(0);

	//	cout << "Creating a Cube now" << endl;
	for (int z = -d_z / 2; z < d_z / 2 + 1; z++)
	{
		for (int y = -d_xy / 2; y < d_xy / 2 + 1; y++)
		{
			for (int x = -d_xy / 2; x < d_xy / 2 + 1; x++)
			{
				xyz_offset[0] = z;
				xyz_offset[1] = y;
				xyz_offset[2] = x;
				offsets.push_back(xyz_offset);
			}
		}
	}
	this->CopyVectorToPointers();
}

// constructor with 
StrEle::StrEle(C3DImage* data3D, const char* ele_shape, int diameter_xy,
	int diameter_z)
{
	data = data3D;
	shape = ele_shape;
	d_xy = diameter_xy;
	d_z = diameter_z;
	weights = NULL;

	// do some error checking
	if (!(diameter_xy % 2) || diameter_xy < 0)
	{
		diameter_xy += 1;
	}
	if (!(diameter_z % 2) || diameter_z < 0)
	{
		diameter_z += 1;
	}

	if (!strcmp(ele_shape, "sphere"))
		CreateSphere();

	if (!strcmp(ele_shape, "cube"))
		CreateCube();


	if (!(!strcmp(ele_shape, "sphere") || !strcmp(ele_shape, "cube")))
	{
		cerr << "Invalid structuring element shape " << shape << endl;
		exit(0);
	}
}

void StrEle::SetCenter(int s, int r, int c)
{
	center_s = s;
	center_r = r;
	center_c = c;
}

void StrEle::DisplayVector(vector<vector<int> > a_vector, ostream& out)
{
	std::vector<int> values;
	// check list size and display the content of the list
	if (a_vector.size() == 0)
	{
		out << "[ (empty vector) ]" << endl;
	}
	else
	{
		out << "{ " << endl;
		for (vector<vector<int> >::iterator i = a_vector.begin();
			i != a_vector.end();
			++i)
		{
			values = *i;
			out << "[\t";
			for (int j = 0; j < static_cast<int>(values.size()); j ++)
			{
				out << values[j] << "\t";
			}
			out << "] " << endl;
		}
		out << "}" << endl;
	}
}


// return the maximum value within the structuring element
unsigned char StrEle::Max()
{
	vector<int> values;
	unsigned char curr_max;
	unsigned char temp;
	curr_max = 0;
	for (int i = 0; i < size; i++)
	{
		temp = data->data[center_s + offsets2[i][0]][center_r + offsets2[i][1]][center_c + offsets2[i][2]];
		if (temp > curr_max)
		{
			curr_max = temp;
		}
	}
	return curr_max;
}

// return the minimum value within the structuring element
unsigned char StrEle::Min()
{
	vector<int> values;
	unsigned char curr_min;
	unsigned char temp;
	curr_min = 255;
	for (int i = 0; i < size; i++)
	{
		temp = data->data[center_s + offsets2[i][0]][center_r + offsets2[i][1]][center_c + offsets2[i][2]];
		if (temp < curr_min)
		{
			curr_min = temp;
		}
	}
	return curr_min;
}

// this returns the median value within the structuring element
unsigned char StrEle::Median()
{
	vector<int> values;
	list<unsigned char> intensities;
	unsigned char median = 0;
	unsigned char temp;

	/*	for(vector< vector<int> >::iterator i = offsets.begin(); i != offsets.end(); ++i)  {
	values = *i;
	intensities.push_back(data->data[center_s + values[0]][center_r + values[1]][center_c + values[2]]);
	}*/
	for (int i = 0; i < size; i++)
	{
		temp = data->data[center_s + offsets2[i][0]][center_r + offsets2[i][1]][center_c + offsets2[i][2]];
		intensities.push_back(temp);
	}

	list<unsigned char>::iterator a = intensities.begin();
	intensities.sort();
	int size = intensities.size();
	if (size % 2 == 0)
	{
		for (int j = 0; j < size / 2 + 1 ; j++)
		{
			++a;
			if (j == size / 2)
				median = static_cast<unsigned char>(median + *a);
			if (j == size / 2 + 1)
				median = static_cast<unsigned char>(median + *a);
		}
		median /= 2;
		cout << "median, list size is an even number" << endl;
	}
	else
	{
		for (int j = 0; j < size / 2; j++)
			++a;
		median = *a;
	}

	return median;
}

// get the z bound of the structuring element
int StrEle::GetZBound()
{
	return d_z / 2;
}

// returns the x and/or y bound of the structuring element
int StrEle::GetXYBound()
{
	return d_xy / 2;
}

bool StrEle::WithinImagePadding() 
{
	// check all sides coordinates to be within the image padding 
	return (data->WithinImagePadding(center_s + GetZBound(), center_r + GetXYBound(), center_c + GetXYBound(), giMARGIN) &&
		data->WithinImagePadding(center_s - GetZBound(), center_r + GetXYBound(), center_c + GetXYBound(), giMARGIN) &&
		data->WithinImagePadding(center_s + GetZBound(), center_r + GetXYBound(), center_c + GetXYBound(), giMARGIN) &&
		data->WithinImagePadding(center_s + GetZBound(), center_r - GetXYBound(), center_c + GetXYBound(), giMARGIN) &&
		data->WithinImagePadding(center_s + GetZBound(), center_r + GetXYBound(), center_c + GetXYBound(), giMARGIN) &&
		data->WithinImagePadding(center_s + GetZBound(), center_r + GetXYBound(), center_c - GetXYBound(), giMARGIN));
}

// a function to copy a vector to pointers
void StrEle::CopyVectorToPointers()
{
	size = offsets.size();
	vector<int> values;
	// copy this vector of vectors to a simple 2D pointer array
	offsets2 = new int * [size];

	for (int i = 0; i < size; i++)
	{
		values = offsets[i];
		offsets2[i] = new int[3];
		offsets2[i][0] = values[0];
		offsets2[i][1] = values[1];
		offsets2[i][2] = values[2];
	}
}

// this function adds gaussian weights to each offset of the StrEle
void StrEle::AddGaussianWeight(float sigma_x, float sigma_y, float sigma_z)
{
	if (weights != NULL)
	{
		cerr << "StrEle:AddGaussianWeight: Weights already defined" << endl;
		cout << weights[0] << endl;
		exit(0);
	}
	int i;
	ofstream out("gaussweights.txt");
	weights = new float[size];
	float sum = 0.0;
	for (i = 0; i < size; i++)
	{
		weights[i] = Gaussian(static_cast<float>(offsets2[i][2]), sigma_x) 
			* Gaussian(static_cast<float>(offsets2[i][1]), sigma_y) 
			* Gaussian(static_cast<float>(offsets2[i][0]), sigma_z);
		
		out << "[\t";
		for (int j = 0; j < 3; j ++)
		{
			out << offsets[i][j] << "\t";
		}
		out << "] ";
		out << "[\t";
		out << (Gaussian(static_cast<float>(offsets2[i][2]), sigma_x)) << "\t"
			<< (Gaussian(static_cast<float>(offsets2[i][1]), sigma_y)) << "\t"
			<< (Gaussian(static_cast<float>(offsets2[i][0]), sigma_z)) << "\t";
		out << "] " << "\t";
		out << " Weight = " << weights[i] << endl;

		sum += weights[i];
	}

	// normalize
	float inv_sum = 1.0f / sum;
	for (i = 0; i < size; i++)
	{
		//weights[i] = weights[i] / sum;
		weights[i] = weights[i] * inv_sum;
	}


	if (sum > 1.0 && sum < 0.0)
		cout << "kernel sum = " << sum << endl;
}

// this function will normalize the gaussian weights so that 
// the value of the weights at (0,0,0) will be 1.0
// i.e. the maximum weight will be 1.0
void StrEle::NormalizeGaussianWeight()
{
	int i;
	float curr_max = 0.0;
	for (i = 0; i < size; i++)
	{
		if (weights[i] > curr_max)
		{
			curr_max = weights[i];
		}
	}

	float inv_curr_max = 1.0f / curr_max;
	for (i = 0; i < size; i++)
	{
		weights[i] = weights[i] * inv_curr_max;
	}
}

// this function convolves the weights with the 3D image data
unsigned char StrEle::Convolve()
{
	//ofstream out("Convolve.txt");
	if (weights == NULL)
	{
		cerr << "StrEle:Convolve: No weights associated with kernel" << endl;
		exit(0);
	}

	float sum = 0.0;
	for (int i = 0; i < size; i++)
	{
		sum += (float)
			data->data[center_s + offsets2[i][0]][center_r + offsets2[i][1]][center_c + offsets2[i][2]] * weights[i];
		//	outC << "i = " << i << "\tsum = " << sum << "\tdata = " <<  (float)data->data[center_s + offsets2[i][0]][center_r + offsets2[i][1]][center_c + offsets2[i][2]] << "\tweight[i] = " << weights[i] << endl;
	}

	//cout << sum << endl;
	sum = static_cast<float>(Round(sum));
	//outC << "SUM = " << (int)((unsigned char) sum) << endl;
	return (unsigned char) sum;
}

// function retired: 3-2-04: amri
// Returns the line filter response using Kikinis vessel model
//unsigned char StrEle::LineFilter()
//{
//	using namespace std;
//	using namespace techsoft;
//
//	typedef double Type;
//	typedef matrix<Type> Matrix;
//	typedef valarray<Type> Vector;
//
//	// set up the hessian matrix (3x3 initialized to 0.0)
//	Matrix Hessian = Matrix(3, 3, 0.0);
//	Matrix ev(3, 3);
//	Vector d(3);
//
//	// Ixx = dI^2 / dx^2
//	Hessian[0][0] = (float) data->data[center_s][center_r][center_c + 1] -
//		2.0 * (float) data->data[center_s][center_r][center_c] +
//		(float) data->data[center_s][center_r][center_c - 1];
//
//	// Iyy = dI^2 / dy^2
//	Hessian[1][1] = (float) data->data[center_s][center_r + 1][center_c] -
//		2.0 * (float) data->data[center_s][center_r][center_c] +
//		(float) data->data[center_s][center_r - 1][center_c];
//
//	// Izz = dI^2 / dz^2
//	Hessian[2][2] = (float) data->data[center_s + 1][center_r][center_c] -
//		2.0 * (float) data->data[center_s][center_r][center_c] +
//		(float) data->data[center_s - 1][center_r][center_c];
//
//	// Ixy = Iyx = dI^2 / dxdy = dI^2 / dydx
//	Hessian[0][1] = Hessian[1][0] = (float)
//		data->data[center_s][center_r + 1][center_c + 1] -
//		(float)
//		data->data[center_s][center_r + 1][center_c] -
//		(float)
//		data->data[center_s][center_r][center_c + 1] +
//		(float)
//		data->data[center_s][center_r][center_c];
//
//	// Ixz = Izx = dI^2 / dxdz = dI^2 / dzdx
//	Hessian[0][2] = Hessian[2][0] = (float)
//		data->data[center_s + 1][center_r][center_c + 1] -
//		(float)
//		data->data[center_s + 1][center_r][center_c] -
//		(float)
//		data->data[center_s][center_r][center_c + 1] +
//		(float)
//		data->data[center_s][center_r][center_c];
//
//	// Iyz = Izx = dI^2 / dxdz = dI^2 / dzdx
//	Hessian[1][2] = Hessian[2][1] = (float)
//		data->data[center_s + 1][center_r + 1][center_c] -
//		(float)
//		data->data[center_s + 1][center_r][center_c] -
//		(float)
//		data->data[center_s][center_r + 1][center_c] +
//		(float)
//		data->data[center_s][center_r][center_c];
//
//
//	bool ret;
//	ret = Hessian.eigen(d, ev); 		 // Finds both eigen values and eigen vectors
//
//	int i;
//	int temp1 = 0, temp2 = 0;
//	float l1, l2, l3;
//	float curr_min = 999.9f;
//	for (i = 0; i < 3; i++)
//	{
//		if (d[i] < curr_min)
//		{
//			curr_min = static_cast<float>(d[i]);
//			temp1 = i;
//		}
//	}
//	l3 = curr_min;
//
//	float curr_max = - 999.9f;
//	for (i = 0; i < 3; i++)
//	{
//		if (d[i] > curr_max)
//		{
//			curr_max = static_cast<float>(d[i]);
//			temp2 = i;
//		}
//	}
//	l1 = curr_max;
//
//	l2 = static_cast<float>(d[3 - temp1 - temp2]);
//
//	/*	cout << l1 << endl;
//	cout << l2 << endl;
//	cout << l3 << endl;
//	*/
//	//	Hessian.~matrix();
//	//	ev.~matrix();
//	//	d.~valarray();
//
//	float l3_abs = - l3;
//	float l2_abs = - l2;
//	float sigma_23 = 0.5;
//	float sigma_12 = 1.0;
//	float alpha = 0.25;
//
//	// the line measure
//	float lambda_123;
//	if (l1 <= 0.0 && l2 < l1 && l3 < l2)
//	{
//		lambda_123 = static_cast<float>(l3_abs 
//			* pow(static_cast<double>(l2 / l3),static_cast<double>(sigma_23)) 
//			* pow(static_cast<double>(1.0 + l1 / l2_abs),static_cast<double>(sigma_12)));
//	}
//	else
//	{
//		if (l3 < l2 && l2 <0.0 && l1> 0 && l1 < l2_abs / alpha)
//		{
//			lambda_123 = static_cast<float>(l3_abs * 
//				pow(static_cast<double>(l2 / l3), static_cast<double>(sigma_23)) * 
//				pow(static_cast<double>(1.0 - alpha * (l1 / l2_abs)), static_cast<double>(sigma_12)));
//		}
//		else
//			lambda_123 = 0.0;
//	}
//
//	lambda_123 = static_cast<float>(Round(lambda_123));
//
//	return static_cast<unsigned char>(lambda_123);
//}


// function retired: 3-2-04: amri
/*unsigned char 
StrEle::LineFilterNorm(float sigma_f, Matrix * Hessian, Vector * d)*/
//unsigned char StrEle::LineFilterNorm(float sigma_f)
//{
//	using namespace std;
//	using namespace techsoft;
//
//	typedef double Type;
//	typedef matrix<Type> Matrix;
//	typedef valarray<Type> Vector;
//
//	// set up the hessian matrix (3x3 initialized to 0.0)
//	Matrix* Hessian = new Matrix(3, 3, 0.0);
//	//	Matrix * ev = new Matrix(3,3);
//	Vector* d = new Vector(3);
//	//	Vector d(3);
//
//	float delta_z = 6.0;
//
//	// Ixx = dI^2 / dx^2
//	(*Hessian)[0][0] = (float) data->data[center_s][center_r][center_c + 1] -
//		2.0 * (float) data->data[center_s][center_r][center_c] +
//		(float) data->data[center_s][center_r][center_c - 1];
//
//	// Iyy = dI^2 / dy^2
//	(*Hessian)[1][1] = (float) data->data[center_s][center_r + 1][center_c] -
//		2.0 * (float) data->data[center_s][center_r][center_c] +
//		(float) data->data[center_s][center_r - 1][center_c];
//
//	// Izz = dI^2 / dz^2
//	(*Hessian)[2][2] = ((float) data->data[center_s + 1][center_r][center_c] -
//		2.0 * (float) data->data[center_s][center_r][center_c] +
//		(float) data->data[center_s - 1][center_r][center_c]) /
//		(delta_z * delta_z);
//
//	// Ixy = Iyx = dI^2 / dxdy = dI^2 / dydx
//	(*Hessian)[0][1] = (*Hessian)[1][0] = (float)
//		data->data[center_s][center_r + 1][center_c + 1] -
//		(float)
//		data->data[center_s][center_r + 1][center_c] -
//		(float)
//		data->data[center_s][center_r][center_c + 1] +
//		(float)
//		data->data[center_s][center_r][center_c];
//
//	// Ixz = Izx = dI^2 / dxdz = dI^2 / dzdx
//	(*Hessian)[0][2] = (*Hessian)[2][0] = ((float)
//		data->data[center_s + 1][center_r][center_c + 1] -
//		(float)
//		data->data[center_s + 1][center_r][center_c] -
//		(float)
//		data->data[center_s][center_r][center_c + 1] +
//		(float)
//		data->data[center_s][center_r][center_c]) /
//		delta_z;
//
//	// Iyz = Izx = dI^2 / dxdz = dI^2 / dzdx
//	(*Hessian)[1][2] = (*Hessian)[2][1] = ((float)
//		data->data[center_s + 1][center_r + 1][center_c] -
//		(float)
//		data->data[center_s + 1][center_r][center_c] -
//		(float)
//		data->data[center_s][center_r + 1][center_c] +
//		(float)
//		data->data[center_s][center_r][center_c]) /
//		delta_z;
//
//	// since we use delta_x = delta_y and delta_z = 6*delta_x
//	// normalize whenever we got delta_z
//
//	bool ret;
//	ret = Hessian->eigen(*d);   	   // Finds both eigen values and eigen vectors
//	delete Hessian;
//
//	int i;
//	int temp1 = 0, temp2 = 0;
//	float l1, l2, l3;
//	float curr_min = 999.9f;
//	for (i = 0; i < 3; i++)
//	{
//		if ((*d)[i] < curr_min)
//		{
//			curr_min = static_cast<float>((*d)[i]);
//			temp1 = i;
//		}
//	}
//	l3 = curr_min;
//
//	float curr_max = - 999.9f;
//	for (i = 0; i < 3; i++)
//	{
//		if ((*d)[i] > curr_max)
//		{
//			curr_max = static_cast<float>((*d)[i]);
//			temp2 = i;
//		}
//	}
//	l1 = curr_max;
//
//	l2 = static_cast<float>((*d)[3 - temp1 - temp2]);
//	delete d;
//
//	/*	cout << l1 << endl;
//	cout << l2 << endl;
//	cout << l3 << endl;
//	*/
//	//	Hessian.~matrix();
//	//	ev.~matrix();
//	//	d.~valarray();
//
//	float l3_abs = - l3;
//	float l2_abs = - l2;
//	float sigma_23 = 0.5;
//	float sigma_12 = 1.0;
//	float alpha = 0.25;
//
//	// the line measure
//	float lambda_123;
//
//	if (l1 <= 0.0 && l2 < l1 && l3 < l2)
//	{
//		lambda_123 = static_cast<float>(l3_abs 
//			* pow(static_cast<double>(l2 / l3), static_cast<double>(sigma_23)) 
//			* pow(static_cast<double>(1.0f + l1 / l2_abs), static_cast<double>(sigma_12)));
//	}
//	else
//	{
//		if (l3 < l2 && l2 <0.0 && l1> 0 && l1 < l2_abs / alpha)
//		{
//			lambda_123 = static_cast<float>(l3_abs 
//				* pow(static_cast<double>(l2 / l3), static_cast<double>(sigma_23)) 
//				* pow(static_cast<double>(1.0f - alpha * (l1 / l2_abs)), static_cast<double>(sigma_12)));
//		}
//			
//		else
//			lambda_123 = 0.0;
//	}
//
//	// normalize with respect to sigma_f
//	lambda_123 = lambda_123 * sigma_f * sigma_f;
//
//	// saturate at 0 and 255
//	if (lambda_123 > 255.0)
//		lambda_123 = 255.0;
//	if (lambda_123 < 0.0)
//		lambda_123 = 0.0;
//
//	lambda_123 = static_cast<float>(Round(lambda_123));
//	return static_cast<unsigned char>(lambda_123);
//}

void StrEle::SetWeight(int s, int r, int c, float weight)
{
	int i;
	//if no weights defined, define one and initialize
	if (weights == NULL)
	{
		weights = new float[size];
		for (i = 0; i < size; i++)
		{
			weights[i] = 0.0;
		}
	}
	for (i = 0; i < size; i++)
	{
		if (offsets[i][0] == s && offsets[i][1] == r && offsets[i][2] == c)
		{
			weights[i] = weight;
		}
	}
	ofstream out("weights.txt");
	out << "start" << endl;
	for (i = 0; i < size; i++)
	{
		out << "[\t";
		for (int j = 0; j < 3; j ++)
		{
			out << offsets[i][j] << "\t";
		}
		out << "] " << "\t";
		out << "weight = " << weights[i] << endl;
	}
}


//////////////////////////////////////
//     Morphological Operators  	//
//////////////////////////////////////

C3DImage* mtk_Dilate(C3DImage* source, StrEle* S)
{
	//C3DImage * result = new C3DImage(*source);
	C3DImage* result = new C3DImage(source->m_iSlices,
						   	source->m_iRows,
						   	source->m_iCols);
	for (int s = S->GetZBound(); s < source->m_iSlices - S->GetZBound(); s++)
	{
		cout << "\r\tMax Operator: slice " << s << " of "
			<< (source->m_iSlices - S->GetZBound() - 1) << flush;
		for (int r = S->GetXYBound();
			r < source->m_iRows - S->GetXYBound();
			r++)
		{
			for (int c = S->GetXYBound();
				c < source->m_iCols - S->GetXYBound();
				c++)
			{
				S->SetCenter(s, r, c);
				result->data[s][r][c] = S->Max();
			}
		}
	}
	cout << endl;
	return result;
}

C3DImage* mtk_Erode(C3DImage* source, StrEle* S)
{
	//C3DImage * result = new C3DImage(*source);
	C3DImage* result = new C3DImage(source->m_iSlices,
						   	source->m_iRows,
						   	source->m_iCols);
	for (int s = S->GetZBound(); s < source->m_iSlices - S->GetZBound(); s++)
	{
		cout << "\r\tMin Operator: slice " << s << " of "
			<< (source->m_iSlices - S->GetZBound() - 1) << flush;
		for (int r = S->GetXYBound();
			r < source->m_iRows - S->GetXYBound();
			r++)
		{
			for (int c = S->GetXYBound();
				c < source->m_iCols - S->GetXYBound();
				c++)
			{
				S->SetCenter(s, r, c);
				result->data[s][r][c] = S->Min();
			}
		}
	}
	cout << endl;
	return result;
}

C3DImage* mtk_Median(C3DImage* source, StrEle* S)
{
	//C3DImage * result = new C3DImage(*source);
	C3DImage* result = new C3DImage(source->m_iSlices,
						   	source->m_iRows,
						   	source->m_iCols);
	for (int s = S->GetZBound(); s < source->m_iSlices - S->GetZBound(); s++)
	{
		cout << "Median Operator: slice " << s << " of "
			<< (source->m_iSlices - S->GetZBound() - 1) << endl;
		for (int r = S->GetXYBound();
			r < source->m_iRows - S->GetXYBound();
			r++)
		{
			for (int c = S->GetXYBound();
				c < source->m_iCols - S->GetXYBound();
				c++)
			{
				S->SetCenter(s, r, c);
				result->data[s][r][c] = S->Median();
			}
		}
	}
	cout << endl;
	return result;
}

// Open:
// Erode then Dilate
C3DImage* mtk_Open(C3DImage* source, StrEle* sErode, StrEle* sDilate = NULL)
{
	if (sDilate == NULL)
		sDilate = sErode;

	C3DImage* temp = mtk_Erode(source, sErode);
	C3DImage* result = mtk_Dilate(temp, sDilate);

	delete temp;
	return result;
}

// Close:
// Dilate then Erode
C3DImage* mtk_Close(C3DImage* source, StrEle* sDilate, StrEle* sErode = NULL)
{
	if (sErode == NULL)
		sErode = sDilate;

	C3DImage* temp = mtk_Dilate(source, sDilate);
	C3DImage* result = mtk_Erode(temp, sErode);

	delete temp;
	return result;
}

// OpenClose:
// Open then Close: suppress peaks
C3DImage* OpenClose(C3DImage* source, StrEle* sOpen, StrEle* sClose = NULL)
{
	if (sClose == NULL)
		sClose = sOpen;

	C3DImage* temp = mtk_Open(source, sOpen);
	C3DImage* result = mtk_Close(temp, sClose);

	delete temp;
	return result;
}

// CloseOpen:
// Close then Open: suppress peaks
C3DImage* CloseOpen(C3DImage* source, StrEle* sClose, StrEle* sOpen = NULL)
{
	if (sOpen == NULL)
		sOpen = sClose;

	C3DImage* temp = mtk_Close(source, sClose);
	C3DImage* result = mtk_Open(temp, sOpen);

	delete temp;
	return result;
}

C3DImage* GaussianFilter(C3DImage* source, StrEle* sGaussian)
{
	C3DImage* result = new C3DImage(source->m_iSlices,
						   	source->m_iRows,
						   	source->m_iCols);

	for (int s = sGaussian->GetZBound();
		s < source->m_iSlices - sGaussian->GetZBound();
		s++)
	{
		cout << "Gaussian Operator: slice " << s << " of "
			<< (source->m_iSlices - sGaussian->GetZBound() - 1) << endl;
		for (int r = sGaussian->GetXYBound();
			r < source->m_iRows - sGaussian->GetXYBound();
			r++)
		{
			for (int c = sGaussian->GetXYBound();
				c < source->m_iCols - sGaussian->GetXYBound();
				c++)
			{
				//cout << "s: " << s << "r: " << r << "c: " << c << endl;
				sGaussian->SetCenter(s, r, c);
				result->data[s][r][c] = sGaussian->Convolve();
			}
		}
	}
	return result;
}

enum DiffMode
{
	dx,
	dy,
	dz,
	dxdy,
	dxdz,
	dydz,
	dx2,
	dy2,
	dz2
};

C3DImage* Differentiation(C3DImage* source, int mode)
{
	// set-up
	StrEle* sDiff = new StrEle(source, "cube", 3, 3);
	switch (mode)
	{
	case dx:
		// equivalent to dc (column)
		cout << "dx" << endl;
		sDiff->SetWeight(0, 0, -1, -1.0); 
		sDiff->SetWeight(0, 0, 0, 1.0);
		break;
	case dy:
		// equivalent to dr (row)
		sDiff->SetWeight(0, -1, 0, -1);
		sDiff->SetWeight(0, 0, 0, 1);
		break;
	case dz:
		// equivalent to ds (slice)
		sDiff->SetWeight(-1, 0, 0, -1);
		sDiff->SetWeight(0, 0, 0, 1); 
	default:
		cerr << "Differentiation: Unrecognizable differentiation mode "
			<< mode << endl;
		exit(0);
	}
	// action
	C3DImage* result = new C3DImage(source->m_iSlices,
						   	source->m_iRows,
						   	source->m_iCols);

	for (int s = sDiff->GetZBound();
		s < source->m_iSlices - sDiff->GetZBound();
		s++)
	{
		cout << "Differential Operator: slice " << s << " of "
			<< (source->m_iSlices - sDiff->GetZBound() - 1) << endl;
		for (int r = sDiff->GetXYBound();
			r < source->m_iRows - sDiff->GetXYBound();
			r++)
		{
			for (int c = sDiff->GetXYBound();
				c < source->m_iCols - sDiff->GetXYBound();
				c++)
			{
				//cout << "s: " << s << "r: " << r << "c: " << c << endl;
				sDiff->SetCenter(s, r, c);
				result->data[s][r][c] = sDiff->Convolve();
			}
		}
	}
	return result;
}

// external interface to differentiation operation

C3DImage* Differentiation(C3DImage* source, char* mode)
{
	/*	source->Write("diffsource1.pic");
	CImage * ProjectionImage;
	ProjectionImage = new CImage(iRows, iCols);
	source->MaxProjectXY(*ProjectionImage);
		ProjectionImage->WriteTIFF("diffsource1.tif");*/

	if (strcmp(mode, "dx") == 0 || strcmp(mode, "dc") == 0)
		return Differentiation(source, dx);
	if (strcmp(mode, "dy") == 0 || strcmp(mode, "dr") == 0)
		return Differentiation(source, dy);
	if (strcmp(mode, "dz") == 0 || strcmp(mode, "ds") == 0)
		return Differentiation(source, dz);
	// symmetry: i.e. dxdy = dydx
	if (strcmp(mode, "dxdy") == 0 || strcmp(mode, "dydx") == 0)
		return Differentiation(source, dxdy);
	if (strcmp(mode, "dxdz") == 0 || strcmp(mode, "dzdx") == 0)
		return Differentiation(source, dxdz);
	if (strcmp(mode, "dydz") == 0 || strcmp(mode, "dzdy") == 0)
		return Differentiation(source, dydz);
	if (strcmp(mode, "dx2") == 0)
		return Differentiation(source, dx2);
	if (strcmp(mode, "dy2") == 0)
		return Differentiation(source, dy2);
	if (strcmp(mode, "dz2") == 0)
		return Differentiation(source, dz2);
	else
	{
		cerr << "Differentiation: Unrecognizable differentiation mode "
			<< mode << endl;
		exit(0);
	}
}

C3DImage* Profile1D(C3DImage* source, int slice, int row, int column,
	char direction)
{
	char fName[200];
	int s,r,c; 
	C3DImage* location = new C3DImage(*source);

	std::string output_path = gConfig.GetOutputPath();
	std::string image_name = gConfig.GetImageName();

	sprintf(fName, "%s%s1DProfile_S%d_R%d_C%d_%c.txt", output_path.c_str(),
		image_name.c_str(), slice, row, column, direction);

	ofstream out(fName);

	switch (direction)
	{
	case 's':
		for (s = 0; s < source->m_iSlices; s++)
		{
			out << (int) source->data[s][row][column] << endl;
			location->data[s][row][column] = 255;
		}
		break;
	case 'r':
		for (r = 0; r < source->m_iRows; r++)
		{
			out << (int) source->data[slice][r][column] << endl;
			location->data[slice][r][column] = 255;
		}
		break;
	case 'c':
		for (c = 0; c < source->m_iCols; c++)
		{
			out << (int) source->data[slice][row][c] << endl;
			location->data[slice][row][c] = 255;
		}
		break;
	default:
		cerr << "Profile1D: Unrecognized direction " << direction << endl;
		exit(0);
	}

	return location;
}

// this function creates a sine tube with a gaussian cross sectional profile
C3DImage* SineImage(char* plane, int diameter, float sigma, char* option)
{
	int** offsets;
	float* weights;
	int s,r,c, i;

	diameter = static_cast<int>(2.0 * Round(3.0 * sigma) + 1.0);

	C3DImage* result = new C3DImage(256, 256, 512);

	StrEle* sGaussian;
	sGaussian = new StrEle(result, "sphere", diameter, diameter);
	sGaussian->AddGaussianWeight(sigma, sigma, sigma);
	sGaussian->NormalizeGaussianWeight();

	offsets = sGaussian->GetOffsets();
	weights = sGaussian->GetWeights();

	if (strcmp(plane, "xy") == 0)
	{
		// create the centerline sine
		s = 256 / 2;
		for (c = 0 + diameter; c < 512 - diameter; c++)
		{
			r = static_cast<int>(sin(0.02 * static_cast<int>(c) ) * 64 + 128);

			for (i = 0; i < sGaussian->GetSize(); i++)
			{
				if (strcmp(option, "gaussian") == 0)
				{
					if ((int)
						result->data[s + offsets[i][0]][r + offsets[i][1]][c + offsets[i][2]] <
						Round(weights[i] * 255.0))
					{
						result->data[s + offsets[i][0]][r + offsets[i][1]][c + offsets[i][2]] = (unsigned char)
							Round(weights[i] * 255.0);
					}
				}
				else
				{
					result->data[s + offsets[i][0]][r + offsets[i][1]][c + offsets[i][2]] = 200;
				}
			}
		}
	}

	if (strcmp(plane, "xz") == 0)
	{
		// create the centerline sine
		r = 256 / 2;
		for (c = 0 + diameter; c < 512 - diameter; c++)
		{
			s = static_cast<int>(sin(0.02 * static_cast<float>(c)) * 64 + 128);

			for (i = 0; i < sGaussian->GetSize(); i++)
			{
				if (strcmp(option, "gaussian") == 0)
				{
					if ((int)
						result->data[s + offsets[i][0]][r + offsets[i][1]][c + offsets[i][2]] <
						Round(weights[i] * 255.0))
					{
						result->data[s + offsets[i][0]][r + offsets[i][1]][c + offsets[i][2]] = (unsigned char)
							Round(weights[i] * 255.0);
					}
				}
				else
				{
					result->data[s + offsets[i][0]][r + offsets[i][1]][c + offsets[i][2]] = 200;
				}
			}
		}
	}

	return result;
}

C3DImage* LineImage(int diameter)
{
	int** offsets;
	StrEle* Ball;
	int s,r,c, i;
	C3DImage* result = new C3DImage(256, 256, 512);
	Ball = new StrEle(result, "sphere", diameter, diameter);

	offsets = Ball->GetOffsets();

	s = 256 / 2;
	r = 256 / 2;
	for (c = 0 + diameter * 2; c < 512 - diameter * 2; c++)
	{
		for (i = 0; i < Ball->GetSize(); i++)
		{
			result->data[s + offsets[i][0]][r + offsets[i][1]][c + offsets[i][2]] = 200;
		}
	}
	return result;
}

C3DImage* WobblyLineImage(int diameter)
{
	int s,r,c;
	C3DImage* result = new C3DImage(64, 256, 512);

	int d;
	s = 64 / 2;
	r = 256 / 2;
	double d_sum = 0.0;
	int d_count = 0;
	for (c = 0 + diameter * 2; c < 512 - diameter * 2; c++)
	{
		cout << c << "\t" << flush;
		d = max(Round(sin(0.05 * (float) c) * 16), diameter);
		cout << d << endl;
		
		StrEle Ball(result, "sphere", d, d);
		Ball.SetCenter(s,r,c);
		Ball.Fill(128);
		d_sum += d;
		d_count++;
	}
	cout << d_sum / (double)d_count << endl;
	return result;
}


