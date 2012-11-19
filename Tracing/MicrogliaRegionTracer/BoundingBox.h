#ifndef BoundingBox_H
#define BoundingBox_H

#include "itkIntTypes.h"

class BoundingBox
{
private:
    itk::uint64_t minX;
    itk::uint64_t minY;
    itk::uint64_t minZ;
    itk::uint64_t maxX;
    itk::uint64_t maxY;
    itk::uint64_t maxZ;
public:
    explicit BoundingBox(itk::uint64_t minX, itk::uint64_t minY, itk::uint64_t minZ, itk::uint64_t maxX, itk::uint64_t maxY, itk::uint64_t maxZ);
    ~BoundingBox();
    
public:
    itk::uint64_t GetMinX();
    itk::uint64_t GetMinY();
    itk::uint64_t GetMinZ();
    itk::uint64_t GetMaxX();
    itk::uint64_t GetMaxY();
    itk::uint64_t GetMaxZ();
};

#endif //BoundingBox_H