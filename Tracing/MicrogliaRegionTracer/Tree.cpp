#include "Tree.h"

#include <iostream>
#include <stdexcept>

Tree::Tree()
{
	this->root = NULL;
}

//Copy constructors, for making deep copies
Tree::Tree(const Tree& old_tree)
{
	Node* new_root = new Node(*(old_tree.root));	//Deep copy of node "root", see Node.h
	this->root = new_root;
	
	this->member_nodes.clear();
	CopyConstructorHelper(new_root, old_tree.root);
}

void Tree::CopyConstructorHelper(Node* const new_node, const Node* const old_node)
{
	new_node->children.clear();
	this->member_nodes.push_back(new_node);

	std::vector<Node *>::const_iterator old_node_children_iter;
	for (old_node_children_iter = old_node->children.begin(); old_node_children_iter != old_node->children.end(); ++old_node_children_iter)
	{
		Node* child = *old_node_children_iter;
		Node* new_child = new Node(*child);
        
        new_child->SetParent(new_node); //Update the reference to point to the new parent
        
        CopyConstructorHelper(new_child, child);

		new_node->AddChild(new_child); //Add each new child
	}
}

Tree::~Tree()
{
    TreeDestructorHelper(root);
    
    //vector and other things implicitly destroyed here
}

void Tree::TreeDestructorHelper(Node* node)
{
    std::vector<Node*> children = node->GetChildren();
    std::vector<Node*>::iterator child_iter;
    
    for (child_iter = children.begin(); child_iter != children.end(); ++child_iter)
    {
        Node* child = *child_iter;
        TreeDestructorHelper(child);
    }
    
    //std::cerr << "Deleting node: " << node->getID() << std::endl;
    delete node;
}

void Tree::SetRoot(Node *root)
{
	member_nodes.push_back(root);
	this->root = root;
}

void Tree::AddNode(Node* node, Node* parent)
{
	//The node itself is responsible for keeping track of its parent, as well as the parent keeping track of the node
	
	if (parent == NULL)
		throw std::runtime_error("AddNode needs a parent");

	member_nodes.push_back(node);
}

std::vector<Node*> Tree::GetMemberNodes()
{
	return member_nodes;
}

Node* Tree::GetRoot()
{
	return root;
}

bool Tree::RemoveNode(Node* node)
{
	std::vector< Node* >::iterator member_nodes_iter;
	for (member_nodes_iter = member_nodes.begin(); member_nodes_iter != member_nodes.end(); ++member_nodes_iter)
	{
		Node* member_node = *member_nodes_iter;

		if (member_node == node)
		{
			member_nodes.erase(member_nodes_iter);
			return true;
		}
	}

	return false;
}

void Tree::GetLeafNodes(std::vector<Node *> &leaf_nodes)
{
	Node* root_node = GetRoot();

	VisitChildrenForLeafNodes(root_node, leaf_nodes);
}

void Tree::VisitChildrenForLeafNodes(Node* node, std::vector<Node*> &leaf_nodes)
{
	std::vector<Node*> children = node->GetChildren();

	if (children.size() == 0)
		leaf_nodes.push_back(node);
	else
	{
		std::vector<Node*>::iterator children_iter;
		for (children_iter = children.begin(); children_iter != children.end(); ++children_iter)
			VisitChildrenForLeafNodes(*children_iter, leaf_nodes);
	}
}
