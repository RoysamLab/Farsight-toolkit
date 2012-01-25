#include "Seed.h"

Seed::Seed(size_t x, size_t y, size_t z)
{
	this->x = x;
	this->y = y;
	this->z = z;
}

size_t Seed::getX() const
{
	return x;
}

size_t Seed::getY() const
{
	return y;
}

size_t Seed::getZ() const
{
	return z;
}