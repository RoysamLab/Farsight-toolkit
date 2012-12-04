#include "BoundingBox.h"

BoundingBox::BoundingBox(itk::uint64_t minX, itk::uint64_t minY, itk::uint64_t minZ, itk::uint64_t maxX, itk::uint64_t maxY, itk::uint64_t maxZ) :
minX(minX), minY(minY), minZ(minZ),
maxX(maxX), maxY(maxY), maxZ(maxZ)
{
}

BoundingBox::~BoundingBox()
{
}

itk::uint64_t BoundingBox::GetMinX() const
{
    return this->minX;
}

itk::uint64_t BoundingBox::GetMinY() const
{
    return this->minY;
}

itk::uint64_t BoundingBox::GetMinZ() const
{
    return this->minZ;
}

itk::uint64_t BoundingBox::GetMaxX() const
{
    return this->maxX;
}

itk::uint64_t BoundingBox::GetMaxY() const
{
    return this->maxY;
}

itk::uint64_t BoundingBox::GetMaxZ() const
{
    return this->maxZ;
}

std::ostream & operator<< (std::ostream & lhs, const BoundingBox & rhs)
{
    // Print the data members of rightOp using leftOp like you would using cout
    
    lhs << "X: [" << rhs.GetMinX() << ", " << rhs.GetMaxX() << "] Y: [" << rhs.GetMinY() << ", " << rhs.GetMaxY() << "] Z: [" << rhs.GetMinZ() << ", " << rhs.GetMaxZ() << "]";
    
    return lhs;
}