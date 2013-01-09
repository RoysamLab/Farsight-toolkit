#include "BoundingBoxFromTraceCalculator.h"
#include <limits>
#include <algorithm>
#include <iostream>

const BoundingBox BoundingBoxFromTraceCalculator::GetBoundingBoxFromTree(const Tree & tree)
{
    itk::uint16_t minX = std::numeric_limits< itk::uint16_t >::max();
    itk::uint16_t maxX = std::numeric_limits< itk::uint16_t >::min();
    itk::uint16_t minY = std::numeric_limits< itk::uint16_t >::max();
    itk::uint16_t maxY = std::numeric_limits< itk::uint16_t >::min();
    itk::uint16_t minZ = std::numeric_limits< itk::uint16_t >::max();
    itk::uint16_t maxZ = std::numeric_limits< itk::uint16_t >::min();
    
    NodeVectorType tree_nodes = tree.GetMemberNodes();
    
	std::cout << "Tree has " << tree_nodes.size() << " nodes" << std::endl;

    NodeVectorType::iterator tree_nodes_iterator;
    
    for (tree_nodes_iterator = tree_nodes.begin(); tree_nodes_iterator != tree_nodes.end(); ++tree_nodes_iterator)
    {
        Node * node = *tree_nodes_iterator;
        minX = std::min< itk::uint64_t >(minX, node->GetX());
        maxX = std::max< itk::uint64_t >(maxX, node->GetX());
        minY = std::min< itk::uint64_t >(minY, node->GetY());
        maxY = std::max< itk::uint64_t >(maxY, node->GetY());
        minZ = std::min< itk::uint64_t >(minZ, node->GetZ());
        maxZ = std::max< itk::uint64_t >(maxZ, node->GetZ());
    }
    
    return BoundingBox(minX, minY, minZ, maxX, maxY, maxZ); //Note the order...
}