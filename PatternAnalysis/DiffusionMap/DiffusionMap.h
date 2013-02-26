
#ifndef DIFFUSION_MAP_H
#define DIFFUSION_MAP_H


#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_real.h"
#include "vnl/algo/vnl_real_eigensystem.h"
#include "mbl/mbl_stats_nd.h"
#include "ftkGraphs/kNearestObjects.h"

#include <vtkSmartPointer.h>
#include <vtkTable.h>
#include <vtkVariant.h>

#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <stdio.h>
#include <math.h>

#ifdef _OPENMP
    #include "omp.h"
#endif

class DiffusionMap
{
public:
	DiffusionMap();
	~DiffusionMap();

	void Initialize(vtkSmartPointer< vtkTable > tbl, bool val = false);
	void ComputeDiffusionMap();
	vtkSmartPointer< vtkTable > GetKDiffusionNeighbors(std::vector<unsigned int> IDs, unsigned int k);

private:
	vtkSmartPointer< vtkTable > featureTable;
	vnl_matrix<double> data_matrix;
	vnl_matrix<double> similarity_matrix;
	vnl_matrix<double> random_walk_matrix;
	vnl_vector<double> EigVals;
	vnl_matrix<double> EigVecs;
	bool multi_cell_struct;
	std::map<int, int> rowIdMap, idRowMap; 
	std::map< unsigned int, std::vector<double> > centroidMap;

};

#endif
