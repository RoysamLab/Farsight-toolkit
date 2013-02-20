#ifndef _SPHERICAL_HARMONIC_TRABSFORM_H_
#define _SPHERICAL_HARMONIC_TRABSFORM_H_

#include <boost/math/special_functions/legendre.hpp>
#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/vnl_vector.h>

#include <algorithm>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#define PI 3.14159265
#define EPS 1e-6 
#define SQROOT 1.41421356

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
	void SortwithIndex( vnl_vector<double> &elearray , int start, int end, std::vector<int > &index);
	void costheta();
	void GetCoefficientSize();
	double legendre(int l, int m, double x);
	void CombiningLon();
	void DataFitting();
	void lmCalculating();
	double Fractial(int l, int m);

	unsigned int numpoints;
	unsigned int numcoord;
	unsigned int size_coefficient;

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