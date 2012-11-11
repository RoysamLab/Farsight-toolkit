#include "Node.h"

#include <algorithm>
#include <stdexcept>
#include <iostream>
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
Node::Node(const Node & old_node)
{
	this->x = old_node.x;
	this->y = old_node.y;
	this->z = old_node.z;
	this->id = old_node.id;

	this->parent = old_node.parent;     //parent is set to old_node parent, the tree copy constructor must take care to regenerate the links to any new parents
	this->children = old_node.children; //children is set to the old_node vector of children, the tree copy constructor must take care to regenerate the links to the children as well
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
    std::vector< Node * > & parents_children = parent->GetChildren();   //the vector of all our parent's children, CAREFUL: Note the type, if you do not declare the type a reference, the operator= will just copy the pointers instead of referencing the original vector
    
//    std::cerr << "Our ID: " << this->id << std::endl;
    if (parent != NULL)
    {
//        std::cerr << "Parent's ID: " << parent->getID() << std::endl;
        std::vector< Node * >::iterator parents_children_iter;
        
//        for (parents_children_iter = parents_children.begin(); parents_children_iter != parents_children.end(); ++parents_children_iter)
//        {
//            Node* parents_child = *parents_children_iter;
//            std::cerr << "Parent's Child's ID: " << parents_child->getID() << std::endl;
//        }
        
        bool erased = false;
        for (parents_children_iter = parents_children.begin(); parents_children_iter != parents_children.end(); ++parents_children_iter)
        {
            Node * parents_child = *parents_children_iter;
//            std::cerr << "Testing Parent's Child's ID: " << parents_child->getID() << std::endl;
            if (parents_child->getID() == this->id)
            {
                parents_children.erase(parents_children_iter);
//                std::cerr << "Erasing ourselves from the parent's list of children" << std::endl;
                std::vector< Node * >::iterator remaining_children_iter;
                
//                std::vector<Node*>& remaining_children = parent->GetChildren();
//                for (remaining_children_iter = remaining_children.begin(); remaining_children_iter != remaining_children.end(); ++remaining_children_iter)
//                {
//                    std::cerr << "Remaining children: " << (*remaining_children_iter)->getID() << std::endl;
//                }
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
	children.push_back(child);
}

void Node::SetParent(Node * const parent)
{
	this->parent = parent;	//Set the new parent
}

std::vector< Node * > & Node::GetChildren() //remember to change this function not to return member variables
{
	return children;
}

itk::uint64_t Node::getID() const
{
	return id;
}

Node * Node::GetParent() const
{
	return parent;
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
