#include "BoundingBox.h"

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
