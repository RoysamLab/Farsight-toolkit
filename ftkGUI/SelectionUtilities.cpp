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