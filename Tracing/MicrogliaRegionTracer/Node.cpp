#include "Node.h"

Node::Node(itk::uint64_t x, itk::uint64_t y, itk::uint64_t z, itk::uint64_t id)
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
	//Search for the reference in the old parent to this child and delete it
	Node* old_parent = this->parent;
	if (old_parent != NULL)
	{
		std::vector<Node*>::iterator children_iter;
		for (children_iter = old_parent->children.begin(); children_iter != old_parent->children.end(); ++children_iter)
		{
			Node* old_parent_child = *children_iter;

			if (this == old_parent_child)
				old_parent_child->children.erase(children_iter);
		}
	}

	parent->AddChild(this);	//Add the child to the new parent
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