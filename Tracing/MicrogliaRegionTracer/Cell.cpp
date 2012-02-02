#include "Cell.h"

Cell::Cell(itk::uint64_t x, itk::uint64_t y, itk::uint64_t z)
{
	this->x = x;
	this->y = y;
	this->z = z;
}

size_t Cell::getX() const
{
	return cell_x;
}

size_t Cell::getY() const
{
	return cell_y;
}

size_t Cell::getZ() const
{
	return cell_z;
}