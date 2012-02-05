#include "Cell.h"

Cell::Cell(itk::uint64_t cell_x, itk::uint64_t cell_y, itk::uint64_t cell_z)
{
	this->cell_x = cell_x;
	this->cell_y = cell_y;
	this->cell_z = cell_z;
}

itk::uint64_t Cell::getX() const
{
	return cell_x;
}

itk::uint64_t Cell::getY() const
{
	return cell_y;
}

itk::uint64_t Cell::getZ() const
{
	return cell_z;
}

void Cell::setSize(Cell::ImageType::SizeType cell_size)
{
	this->cell_size = cell_size;
}

Cell::ImageType::SizeType Cell::getSize()
{
	return cell_size;
}

void Cell::setRequestedSize(Cell::ImageType::SizeType cell_requested_size)
{
	this->cell_requested_size = cell_requested_size;
}

Cell::ImageType::SizeType Cell::getRequestedSize()
{
	return cell_requested_size;
}

void Cell::setShiftIndex(Cell::ImageType::IndexType shift_index)
{
	this->shift_index = shift_index;
}

Cell::ImageType::IndexType Cell::getShiftIndex()
{
	return shift_index;
}


