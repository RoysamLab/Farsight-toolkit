#include "Node.h"

#include <cstddef> // for NULL

//Constructor
Node::Node(double x, double y, double z, itk::uint64_t id)
{
	this->x = x;
	this->y = y;
	this->z = z;
	this->id = id;

	this->parent = NULL;
}

//Copy Constructor
//Note: 
//parent is set to old_node parent, the tree copy constructor must take care to regenerate the links to any new parents
//children is set to the old_node vector of children, the tree copy constructor must take care to regenerate the links to the children as well
Node::Node(const Node& old_node)
{
	this->x = old_node.x;
	this->y = old_node.y;
	this->z = old_node.z;
	this->id = old_node.id;

	this->parent = old_node.parent;
	this->children = old_node.children;
}	

void Node::AddChild(Node *child)
{
	children.push_back(child);
}

void Node::SetParent(Node *parent)
{
	this->parent = parent;	//Set the new parent
}

std::vector<Node*> Node::GetChildren()
{
	return children;
}

itk::uint64_t Node::getID()
{
	return id;
}

Node* Node::GetParent()
{
	return parent;
}

bool Node::RemoveChild(Node* child_to_be_removed)
{
	std::vector< Node* >::iterator children_iterator;
	for (children_iterator = children.begin(); children_iterator != children.end(); ++children_iterator)
	{
		Node* child = *children_iterator;

		if (child_to_be_removed == child)
		{
			child->parent = NULL;
			children.erase(children_iterator);
			return true;
		}
	}

	return false;
}
