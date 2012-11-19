#ifndef BoundingVolumeHierarchyNode_h
#define BoundingVolumeHierarchyNode_h

#include "BoundingBox.h"

class BoundingVolumeHierarchyNode
{
private:
    BoundingBox * box;
    
    BoundingVolumeHierarchyNode * smaller_child;
    BoundingVolumeHierarchyNode * larger_child;

public:
    explicit BoundingVolumeHierarchyNode();
    ~BoundingVolumeHierarchyNode();
    
    
};

#endif //#ifndef BoundingVolumeHierarchyNode_h