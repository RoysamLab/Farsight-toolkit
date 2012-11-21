#ifndef BoundingBox_H
#define BoundingBox_H

#include "itkIntTypes.h"
#include <ostream>

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
    itk::uint64_t GetMinX() const;
    itk::uint64_t GetMinY() const;
    itk::uint64_t GetMinZ() const;
    itk::uint64_t GetMaxX() const;
    itk::uint64_t GetMaxY() const;
    itk::uint64_t GetMaxZ() const;
};

std::ostream & operator<< (std::ostream & lhs, const BoundingBox & rhs); //Support printing from std streams

#endif //BoundingBox_H