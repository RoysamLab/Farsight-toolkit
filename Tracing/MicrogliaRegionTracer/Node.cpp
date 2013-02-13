#include "Node.h"

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <cstddef> // for NULL

//Constructor
Node::Node(double x, double y, double z, itk::uint64_t id) :
parent(NULL),
id(id),
x(x),
y(y),
z(z)
{
}

//Copy Constructor
Node::Node(const Node & old_node) :
parent(old_node.parent), //parent is set to old_node parent, the tree copy constructor must take care to regenerate the links to any new parents
id(old_node.id),
children(old_node.children), //children is set to the old_node vector of children, the tree copy constructor must take care to regenerate the links to the children as well
x(old_node.x),
y(old_node.y),
z(old_node.z)
{
}	

//Destructor
Node::~Node()
{
    //Go to the children (if they exist) and removing the link (set to NULL) to their parent (this node)
    std::vector< Node* >::iterator child_iter;
    
    for (child_iter = this->children.begin(); child_iter != this->children.end(); ++child_iter)
    {
        Node* child = *child_iter;
        child->SetParent(NULL);
    }
    
    //Go to the parent (if it exist) and remove ourselves from their children
    Node * parent = this->GetParent();
    const NodeVectorType & parents_children = parent->GetChildren();   //the vector of all our parent's children, CAREFUL: Note the type, if you do not declare the type a reference, the operator= will just copy the pointers instead of referencing the original vector
    
    if (parent != NULL)
    {
        NodeVectorType::const_iterator parents_children_iter;
        
        bool erased = false;
        for (parents_children_iter = parents_children.begin(); parents_children_iter != parents_children.end(); ++parents_children_iter)
        {
            const Node * parents_child = *parents_children_iter;

            if (parents_child->getID() == this->id)
            {
                parent->RemoveChild(parents_children_iter);
                erased = true;
                break;
            }
        }
        if (!erased)
            throw std::runtime_error("Node destructor has problems since since it couldn't find itself inside its parent's list of children");
    }
    
    //The std::vector named children's destructor is implicitly called here and the destructors on its elements are also called.
    //However, since the vector holds pointers, the destructors for pointers do nothing and the children nodes are not destroyed.
}

void Node::AddChild(Node * const child)
{
	this->children.push_back(child);
}

void Node::SetParent(Node * const parent)
{
	this->parent = parent;	//Set the new parent
}

const Node::NodeVectorType & Node::GetChildren() const //remember to change this function not to return member variables
{
	return this->children;
}

itk::uint64_t Node::getID() const
{
	return this->id;
}

Node * Node::GetParent() const
{
	return this->parent;
}

void Node::RemoveChild(Node * child_to_be_removed)
{
	std::vector< Node* >::iterator children_iterator;
	for (children_iterator = children.begin(); children_iterator != children.end(); ++children_iterator)
	{
		Node* child = *children_iterator;

		if (child_to_be_removed == child)
		{
			child->parent = NULL;
			children.erase(children_iterator);
            return;
		}
	}

    throw std::runtime_error("Error, trying to remove a child that does not exist");
}

void Node::RemoveChild(const NodeVectorType::const_iterator & child_const_iter)
{
    NodeVectorType::iterator child_iter(this->children.begin());
    std::advance(child_iter, std::distance< NodeVectorType::const_iterator >(child_iter, child_const_iter));
    
    this->children.erase(child_iter);
}

void Node::ClearChildren()
{
    this->children.clear();
}

const double Node::GetX()
{
    return this->x;
}

const double Node::GetY()
{
    return this->y;
}

const double Node::GetZ()
{
    return this->z;
}

