#ifndef LCSANALYSIS_H
#define LCSANALYSIS_H
#include <cstring>
#include <vector>
class Cell
{
public:
	Cell()
	{
	};
	void set(int r, int c, int s, Cell *pt = NULL)
	{
		row = r;
		col = c;
		score = s;
		prevCell = pt;
	};
	Cell *prevCell;
	int score;
	int row;
	int col;
};

class LCSTable
{
public:
   LCSTable();
   virtual ~LCSTable();
   void Initialize( std::vector<int> &first, std::vector<int> &second);
   void FillInTable();
   int GetTraceback(std::vector<int> &trace);
   
private:
   int row_num;
   int col_num;
   Cell **scoreTable;
   int *row_data;
   int *col_data;
   int length;
   Cell *traceCell;
};
#endif