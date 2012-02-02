#ifndef NODE_H
#define NODE_H

#include "itkIntTypes.h"
#include <vector>

class Node
{
private:
	Node *parent;
	std::vector<Node*> children;
	
	itk::uint64_t id;

public:
	itk::uint64_t x;
	itk::uint64_t y;
	itk::uint64_t z;

	

public:
	Node(itk::uint64_t x, itk::uint64_t y, itk::uint64_t z, itk::uint64_t id);

	void AddChild(Node *child);
	void SetParent(Node *parent);
	std::vector<Node*> GetChildren();

	itk::uint64_t getID();
	Node* GetParent();
};

#endif