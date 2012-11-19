#include "BoundingBox.h"

BoundingBox::BoundingBox(itk::uint64_t minX, itk::uint64_t minY, itk::uint64_t minZ, itk::uint64_t maxX, itk::uint64_t maxY, itk::uint64_t maxZ) :
minX(minX), minY(minY), minZ(minZ),
maxX(maxX), maxY(maxY), maxZ(maxZ)
{
}

itk::uint64_t BoundingBox::GetMinX()
{
    return this->minX;
}

itk::uint64_t BoundingBox::GetMinY()
{
    return this->minY;
}

itk::uint64_t BoundingBox::GetMinZ()
{
    return this->minZ;
}

itk::uint64_t BoundingBox::GetMaxX()
{
    return this->maxX;
}

itk::uint64_t BoundingBox::GetMaxY()
{
    return this->maxY;
}

itk::uint64_t BoundingBox::GetMaxZ()
{
    return this->maxZ;
}
