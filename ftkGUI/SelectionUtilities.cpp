#include "SelectionUtilities.h"
vtkSelection * SelectionUtilities::ConvertIDsToVTKSelection(vtkIdTypeArray * vtkIDs)
{
	/*!
	* create a vtk selection from a vtk id array
	*/
	vtkSmartPointer<vtkSelectionNode> selectNodeList = vtkSmartPointer<vtkSelectionNode>::New();
	selectNodeList->SetSelectionList( vtkIDs );
	selectNodeList->SetFieldType( vtkSelectionNode::VERTEX );
	selectNodeList->SetContentType( vtkSelectionNode::INDICES );

	vtkSelection * TableRowSelection = vtkSelection::New();
	TableRowSelection->RemoveAllNodes();
	TableRowSelection->AddNode(selectNodeList);
	return TableRowSelection;
}
vtkIdTypeArray * SelectionUtilities::ConvertVTKSelectionToIDArray(vtkSelection * inputSelection)
{
	/*!
	* Convert a vtk selection into vtk id array
	*/
	vtkSelectionNode* vertices = NULL;
	vtkIdTypeArray* vertexList;
	if( inputSelection->GetNode(0))
	{
		if( inputSelection->GetNode(0)->GetFieldType() == vtkSelectionNode::VERTEX)
		{
			vertices = inputSelection->GetNode(0);
		}
	}

	if( inputSelection->GetNode(1))
	{
		if( inputSelection->GetNode(1)->GetFieldType() == vtkSelectionNode::VERTEX)
		{
			vertices = inputSelection->GetNode(1);
		}
	}

	if( vertices != NULL)
	{
		vertexList = vtkIdTypeArray::SafeDownCast(vertices->GetSelectionList());
	}
	return vertexList;
}

std::set< vtkIdType > SelectionUtilities::ObjectSelectionToIDSet(std::set<long int> curSel)
{	
	/*! 
	* convert objectSelection into form selective
	* clustering can use for operations
	*/
	std::set<long int>::iterator iter = curSel.begin();
	std::set< vtkIdType > Selection;
	for (; iter != curSel.end(); iter++)
	{
		vtkIdType id = (vtkIdType) (*iter);
		Selection.insert(id);
	}
	return Selection;
}

void SelectionUtilities::RemoveRowAndReMapTable(vtkIdType key, vtkSmartPointer<vtkTable> modTable, std::map< vtkIdType, vtkIdType> TableIDMap)
{
	/*!
	* remove a row from table
	* rebuild map
	*/
	std::map< vtkIdType, vtkIdType>::iterator TableIDIter = TableIDMap.find(key);
	if (TableIDIter != TableIDMap.end())
	{
		modTable->RemoveRow((*TableIDIter).second);
	}
	for (vtkIdType row = 0; row != modTable->GetNumberOfRows(); row++)
	{
		vtkIdType value = modTable->GetValue(row,0).ToTypeInt64();
		TableIDMap[value] = row;
	}
}