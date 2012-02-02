#ifndef TREE_H
#define TREE_H

#include "Node.h"
#include <stdexcept>

class Tree
{
private:
	Node *root;
	std::vector<Node*> member_nodes;

public:
	Tree();
	~Tree();

	void SetRoot(Node *root);

	void AddNode(Node* node, Node* parent);

	std::vector<Node*> GetMemberNodes();

	Node* getRoot();
};

#endif
