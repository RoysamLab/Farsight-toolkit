#ifndef _SPHERICAL_HARMONIC_TRABSFORM_H_
#define _SPHERICAL_HARMONIC_TRABSFORM_H_

#include <boost/math/special_functions/legendre.hpp>
#include <vnl/vnl_matrix.h>
#include <vnl_svd.h>
#include <vnl_vector.h>

#include <algorithm>
#include <vector>
#include <math.h>
#include <stdio.h>

#define PI 3.14159265

class SPH_Trasform
{
public:
	SPH_Trasform();
	~SPH_Trasform();
	void Initialize( vnl_vector<double> x, vnl_vector<double> y, vnl_vector<double> z);
	void Initialize( std::vector<std::vector<double> > data, int L);
	void DataFilter();
	void Legendre_P();
	void Transform();
	void WriteFile();

	int L;

private:
	void xyz2rll();
	std::vector<int > SortwithIndex( vnl_vector<double> elearray, int start, int end);
	void costheta();
	void GetCoefficientSize();
	double legendre(int l, int m, double x);
	void CombiningLon();
	void DataFitting();
	void lmCalculating();

	int numpoints;
	int numcoord;
	int size_coefficient;

	vnl_vector<double> x;
	vnl_vector<double> y;
	vnl_vector<double> z;
	vnl_vector<double> R;
	vnl_vector<double> theta;
	vnl_vector<double> phi;
	vnl_vector<double> fphth;
	vnl_vector<double> lat;
	vnl_vector<double> lon;
	vnl_vector<double> coslat;
	vnl_vector<double> C;

	vnl_matrix<double> lmcosi;
	vnl_matrix<double> plm;

};
#endif