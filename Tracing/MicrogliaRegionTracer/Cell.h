#ifndef Cell_H
#define Cell_H

#include "itkIntTypes.h"

class Cell
{
public:
	Cell(itk::uint64_t cell_x, itk::uint64_t cell_y, itk::uint64_t cell_z);
	
	itk::uint64_t getX() const;
	itk::uint64_t getY() const;
	itk::uint64_t getZ() const;

private:
	itk::uint64_t cell_x;
	itk::uint64_t cell_y;
	itk::uint64_t cell_z;
};

#endif