#ifndef TREE_H
#define TREE_H

#include "Kmeans.h"
#include "clusclus.h"
#include "LocalGeometryRef.h"

typedef struct 
{
	int numm_levels;
	std::vector<std::vector<std::vector<int > > > level_id;

} Level_id;

class BiTree
{
public:
	BiTree();
	~BiTree();
	void constructTree();
	void constructSpecTree();
	void WriteFile(const char *filename1);
	void setFirstIterationClusterNumber(int firstK);
	void setReorderFlag(bool value);
	void setDataToTree(std::vector<std::vector<double > > & originalmatrix,std::vector<std::vector<double > > & matrix);

	bool reorderFlag;
	Level_id getLevels();
	std::vector<std::vector<double > > getPartitionTree();

private:
	int firstK;
	int num_rows;
	int num_cols;
	int num_rows_original;
	int num_cols_original;

	Level_id levels;
	std::vector<std::vector<double > > matrix;
	std::vector<std::vector<double > > originalmatrix;

	void InitialLevelId();
	void append(std::vector<std::vector<double > > & datatobekmeans);

	void reOrder(std::vector<std::vector<double > > submatrix, std::vector<int > & ctp);
	std::vector<double > computerMeanSubmatrix(std::vector<std::vector<double > > & submatrix);
	std::vector<std::set<int > > KmeansClustering(std::vector<std::vector<double > > & datatobekmeans, int num_cluster);
	void normalize(std::vector<std::vector<double > > datatobenormalized, std::vector<std::vector<double > > & normalized);
};

#endif