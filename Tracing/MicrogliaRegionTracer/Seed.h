#ifndef SEED_H
#define SEED_H

#include "itkIntTypes.h"

class Seed
{
public:
	Seed(itk::uint64_t x, itk::uint64_t y, itk::uint64_t z);
	
	itk::uint64_t getX() const;
	itk::uint64_t getY() const;
	itk::uint64_t getZ() const;

private:
	itk::uint64_t x;
	itk::uint64_t y;
	itk::uint64_t z;


};

#endif