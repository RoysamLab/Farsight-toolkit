#include "TraceBit.h"
#include "TraceLine.h"
#include "NodeModel.h"

NodeModel::NodeModel()
{
	this->DataTable = vtkSmartPointer<vtkTable>::New();
	this->Selection = new ObjectSelection();

	this->additionalHeaders.clear();
}

NodeModel::NodeModel(std::vector<TraceLine*> trace_lines)
{
	this->DataTable = vtkSmartPointer<vtkTable>::New();	
	this->Selection = new ObjectSelection();
	this->SetLines(trace_lines);
}

NodeModel::~NodeModel()
{
	delete this->Selection;
	this->Selection = NULL;
}

void NodeModel::SetLines(std::vector<TraceLine*> trace_lines)
{
	this->TraceLines.clear();
	this->TraceLines = trace_lines;
	this->SyncModel();
}

void NodeModel::SetupHeaders()
{
	this->headers.clear();
	this->headers.push_back("Bit ID");
	this->headers.push_back("X");
	this->headers.push_back("Y");
	this->headers.push_back("Z");
	this->headers.push_back("Radius");
	this->headers.push_back("Segment ID");
	this->headers.push_back("Root ID");
	
	int size = this->additionalHeaders.size();
	for (int k = 0; k < size; k++)
	{	
		this->headers.push_back(QString(this->additionalHeaders[k].c_str()));
	}

	int numHeaders = (int)this->headers.size();
	
	vtkSmartPointer<vtkVariantArray> column = vtkSmartPointer<vtkVariantArray>::New();
	for(int i=0; i < numHeaders; ++i)
    {		
		column = vtkSmartPointer<vtkVariantArray>::New();
		column->SetName( this->headers.at(i).toStdString().c_str() );
		this->DataTable->AddColumn(column);
    }
	
}

void NodeModel::AddNodeHeader(std::string NewFeatureHeader)
{
	/*!
	*
	*/
	this->additionalHeaders.push_back(NewFeatureHeader);
	this->headers.push_back(QString(NewFeatureHeader.c_str()));
}

void NodeModel::SyncModel()
{
	int size = this->additionalHeaders.size();
	this->DataTable->Initialize();
	this->Selection->clear();
	this->SetupHeaders();
	for (int i = 0; i < (int)this->TraceLines.size(); ++i)
	{
		TraceLine::TraceBitsType::iterator iter = this->TraceLines.at(i)->GetTraceBitIteratorBegin();
		TraceLine::TraceBitsType::iterator iterend = this->TraceLines.at(i)->GetTraceBitIteratorEnd();
		while(iter != iterend)
		{
			if (size == 0)
			{
				vtkSmartPointer<vtkVariantArray> DataRow = vtkSmartPointer<vtkVariantArray>::New();
				DataRow = (*iter).DataRow();
				DataRow->InsertNextValue(this->TraceLines.at(i)->GetId());
				DataRow->InsertNextValue(this->TraceLines.at(i)->GetRootID());
				this->DataTable->InsertNextRow(DataRow);
			}
			iter++;
		}
	}
}

void NodeModel::setDataTable(vtkSmartPointer<vtkTable> table)
{
	this->DataTable = table;
	this->Selection->clear();
}

vtkSmartPointer<vtkTable> NodeModel::getDataTable()
{
	return this->DataTable;
}

void NodeModel::SelectByIDs(std::vector<int> IDs)
{
	//this->Selection->clear();
	std::set<long int> ID;
	std::vector<int>::iterator IDs_iterator;
	
	//std::cout << "IDs.size(): " << IDs.size() << std::endl;
	for (IDs_iterator = IDs.begin(); IDs_iterator != IDs.end(); IDs_iterator++)
	{
		ID.insert(*IDs_iterator);
	}
	Selection->select(ID);
}

std::vector<long int> NodeModel::GetSelectedIDs()
{
	std::vector<long int> SelectedIDs;
	std::set<long> selected = this->Selection->getSelections();
	std::set<long>::iterator it;
	for (it = selected.begin(); it != selected.end(); ++it)
	{		
		SelectedIDs.push_back(*it);
	}
	return SelectedIDs;
}

//std::vector<TraceBit> NodeModel::GetSelectedNodes()
//{
//	std::vector<TraceBit> selectedNodes;
//	std::set<long> selected = this->Selection->getSelections();
//	std::set<long>::iterator it;
//	for (it = selected.begin(); it != selected.end(); ++it)
//	{
//		this->NodeIDLookupIter = this->TraceBits.find((int) *it);
//		if (this->NodeIDLookupIter != this->TraceBits.end())
//		{
//			selectedNodes.push_back((*this->NodeIDLookupIter));
//		}
//	}//end for selected
//	
//	//std::vector<long int> IDList = this->GetSelectedIDs();
//
//	////Search for nodes
//	//for ( unsigned int i = 0; i< IDList.size(); i++)
//	//{
//	//	this->NodeIDLookupIter = this->NodeIDLookupMAP.find(IDList[i]);
//	//	if (this->NodeIDLookupIter != this->NodeIDLookupMAP.end())
//	//	{
//	//		selectedNodes.push_back((*this->NodeIDLookupIter).second);
//	//	}
//	//}//finished with id search
//	return selectedNodes;
//}

unsigned int NodeModel::getNodeCount()
{
	return (unsigned int) this->TraceBits.size();
}

//std::map< int ,TraceBit>::iterator NodeModel::GetNodeiterator()
//{
//	return this->TraceBits.begin();
//}
//
//std::map< int ,TraceBit>::iterator NodeModel::GetNodeiteratorEnd()
//{
//	return this->TraceBits.end();
//}

ObjectSelection * NodeModel::GetObjectSelection()
{
	return this->Selection;
}