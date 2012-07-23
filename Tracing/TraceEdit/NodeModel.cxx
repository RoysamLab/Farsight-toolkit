#include "NodeModel.h"

NodeModel::NodeModel()
{
	this->DataTable = vtkSmartPointer<vtkTable>::New();
	this->Selection = new ObjectSelection();
}

NodeModel::~NodeModel()
{
}

void NodeModel::nodeHeaders()
{
	this->headers.push_back("Bit ID");
}