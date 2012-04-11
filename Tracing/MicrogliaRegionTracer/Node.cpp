#include "Node.h"

#include <cstddef> // for NULL

Node::Node(double x, double y, double z, itk::uint64_t id)
{
	this->x = x;
	this->y = y;
	this->z = z;
	this->id = id;

	this->parent = NULL;
}

void Node::AddChild(Node *child)
{
	children.push_back(child);
}

void Node::SetParent(Node *parent)
{
	////Search for the reference in the old parent to this child and delete it
	//Node* old_parent = this->parent;

	//parent->AddChild(this);	//Add the child to the new parent
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
