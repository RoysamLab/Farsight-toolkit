#ifndef TREE_H
#define TREE_H

#include "Node.h"

class Tree
{
private:
	typedef Node::NodeVectorType    NodeVectorType;
    
    Node *                  root;
	std::vector<Node*>      member_nodes;

public:
                            explicit Tree();
                            Tree(const Tree & old_tree);	//Copy constructor, for making deep copies!
                            ~Tree();
    void                    TreeDestructorHelper(Node * const node);

	void                    SetRoot(Node * const root);
	void                    AddNode(Node * const node, const Node* parent);
	NodeVectorType          GetMemberNodes() const;
	Node *                  GetRoot() const;
	void                    RemoveNode(const Node * const node);
	void                    GetLeafNodes(NodeVectorType & leaf_nodes) const;
	void                    VisitChildrenForLeafNodes(Node * node, NodeVectorType & leaf_nodes) const;

private:
	void                    CopyConstructorHelper(Node * const new_node, const Node * const old_node);
};

#endif
