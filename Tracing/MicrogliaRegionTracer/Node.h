#ifndef NODE_H
#define NODE_H

#include "itkIntTypes.h"
#include <vector>
#include <cstddef>

class Node
{
private:
	Node *parent;
	std::vector<Node*> children;
	
	itk::uint64_t id;

public:
	double x;
	double y;
	double z;

	

public:
	Node(double x, double y, double z, itk::uint64_t id);

	//Adds the child
	//Note: This does not make the child update its reference to the parent
	void AddChild(Node *child);
	
	//Removes the child from this node, also removing the child's link to the parent
	//Returns true if the child was successfully removed, else returns false if the child to be removed was not found in the list of children
	bool RemoveChild(Node* child_to_be_removed);

	void SetParent(Node *parent);
	
	std::vector<Node*> GetChildren();

	itk::uint64_t getID();
	
	Node* GetParent();
};

#endif