#ifndef BoundingBoxFromTraceCalculator_h
#define BoundingBoxFromTraceCalculator_h

#include "BoundingBox.h"
#include "Tree.h"

class BoundingBoxFromTraceCalculator
{
private:
    typedef Tree::NodeVectorType NodeVectorType;
    
private:
    explicit BoundingBoxFromTraceCalculator(); 
    ~BoundingBoxFromTraceCalculator();
    
public:    
    static const BoundingBox GetBoundingBoxFromTree(const Tree & tree);
};

#endif //#ifndef BoundingBoxFromTraceCalculator_h