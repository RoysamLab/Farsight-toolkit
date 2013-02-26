#ifndef _LOCALGEOMETRYREF_H_
#define _LOCALGEOMETRYREF_H_

//std_inclide
#include <vector>
#include <string>
#include <math.h>
#include <stdio.h>
#include <sstream>
#include <iostream>
#include <algorithm>

//vxl_include
#include "vnl/vnl_real.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"
#include "vnl/algo/vnl_real_eigensystem.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"
#include "mbl/mbl_stats_nd.h"

//vtk_include
#include <vtkSmartPointer.h>
#include <vtkTable.h>
#include <vtkVariant.h>

//boost_include
#include <boost/math/special_functions.hpp>

//ftk_include
#include "ftkGraphs/kNearestObjects.h"


class LocalGeometryRef
{
public:
	LocalGeometryRef();
	~LocalGeometryRef();
	void Initialize(vnl_matrix<double> data);
	void Normalize();
	void Initialize(vtkSmartPointer< vtkTable > table);
	void Initialize(std::vector<std::vector<double > > data);
	void ComputeSimilarityMatrix(char dis = 'G', int k = 10);
	void ComputeProbabilityMatrix();
	void SVD(int num_eigen);
	vnl_matrix<double> GetSimilarityMatrix();
	vnl_matrix<double> GetProbabilityMatrix();
	vnl_vector<double> GetEigenValues();
	vnl_matrix<double> GetEigenVectors();
	std::vector<std::vector<double > > GetEigenVectors(bool type);

private:
	int num_row;
	int num_col;

	vnl_vector<double> eigVals;
	vnl_matrix<double> eigVecs;
	vnl_matrix<double> data_matrix;
	vnl_matrix<double> nor_data_matrix;
	vnl_matrix<double> similarity_matrix;
	vnl_matrix<double> probability_matrix;
	//vnl_matrix<vcl_complex<double> > D12;
	//vnl_matrix<vcl_complex<double> > Dn12;
	vnl_matrix< double > D12;
	vnl_matrix< double > Dn12;

	double ComputeSimilarity(vnl_vector<double> dt1, vnl_vector<double> dt2, char dis);
	std::vector<std::vector< std::pair<unsigned int, double> > > KnnSearch(int k);

};
#endif
