#ifndef BICLUSTER_H
#define BICLUSTER_H

#include <queue>
#include "Tree.h"

class Bicluster
{
public:
	Bicluster();
	~Bicluster();
	void setDataToBicluster(std::vector<std::vector<double > > & data);
	void biclustering();
	void WriteFile(const char *filename1, const char *filename2);

	std::vector<int > order1;
	std::vector<int > order2;
	Level_id levels1;
	Level_id levels2;


private:
	int num_rows;
	int num_cols;
	int num_iteration;

	std::queue<int> qorder1;
	std::queue<int> qorder2;

	std::vector<std::vector<double > > data;
	std::vector<std::vector<double > > Data1;
	std::vector<std::vector<double > > Data2;

	std::vector<std::vector<double > > transpose(std::vector<std::vector<double > > & M);
	void reOrder(Level_id levels, int flag);
};
#endif