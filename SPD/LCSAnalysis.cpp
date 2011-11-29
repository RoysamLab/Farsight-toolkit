#include "LCSAnalysis.h"

#define DELETE_POINTER(p) if( p != NULL) {delete p; p = NULL;}

LCSTable::LCSTable()
{
    scoreTable = NULL;
    row_data = NULL;
    col_data = NULL;
    length = 0;
    traceCell = NULL;
}

LCSTable::~LCSTable()
{
    DELETE_POINTER(row_data);
    DELETE_POINTER(col_data);
	if( scoreTable != NULL)
	{
		for( int i = 0; i < row_num; i++)
		{
			delete scoreTable[i];
		}
	}
    DELETE_POINTER(scoreTable);
}

void LCSTable::Initialize(std::vector<int> &first, std::vector<int> &second)
{
	DELETE_POINTER(row_data);
    DELETE_POINTER(col_data);
	if( scoreTable != NULL)
	{
		for( int i = 0; i < row_num; i++)
		{
			delete scoreTable[i];
		}
	}
    DELETE_POINTER(scoreTable);

	row_num = first.size() + 1;
	col_num = second.size() + 1;
    row_data = new int[row_num];
    col_data = new int[col_num];
    row_data[0] = -1;
    col_data[0] = -1;
    for( int i = 1; i < row_num; i++)
    {
        row_data[i] = first[i-1];
    }
    for( int i = 1; i < col_num; i++)
    {
        col_data[i] = second[i - 1];
    }

    scoreTable = new Cell*[ row_num];
    for( int i = 0; i < row_num; i++)
    {
		scoreTable[i] = new Cell[col_num];
        for( int j = 0; j < col_num; j++)
        {
            scoreTable[i][j].set(i,j,0,NULL);
        }
    }
}

void LCSTable::FillInTable()
{
    for( int nr = 1; nr < row_num; nr++)
    {
        for( int nc = 1; nc < col_num; nc++)
        {
            Cell cell = scoreTable[nr][nc];
            Cell aboveCell = scoreTable[nr -1][nc];
            Cell leftCell = scoreTable[nr][nc - 1];
            Cell leftCornerCell = scoreTable[nr -1][nc-1];
            int leftCornerScore = leftCornerCell.score;
            if( row_data[nr] == col_data[nc])
            {
                leftCornerScore += 1;
            }
            int max = aboveCell.score >= leftCell.score ? aboveCell.score : leftCell.score;
            cell.prevCell = aboveCell.score >= leftCell.score ? &aboveCell : &leftCell;
            max = leftCornerCell.score >= max ? leftCornerCell.score : max;
            cell.prevCell = leftCornerCell.score >= max ? &leftCornerCell : cell.prevCell;
            cell.score = max;
            
            if( max > length)
            {
                length = max;
                traceCell = &cell;
            }
        }
    }
}

int LCSTable::GetTraceback(std::vector<int> &trace)
{
    trace.clear();
    if( length == 0 || traceCell == NULL)
    {
        return length;
    }
    Cell *preCell = traceCell->prevCell;
    int nr = traceCell->row;
    int nc = traceCell->col;
    trace.push_back( row_data[nr]);
    while( preCell != NULL)
    {
        if( preCell->score == traceCell->score - 1)
        {
            nr = preCell->row;
            trace.push_back( row_data[nr]);
        }
        traceCell = preCell;
        preCell = preCell->prevCell;
    }
    return length;
}