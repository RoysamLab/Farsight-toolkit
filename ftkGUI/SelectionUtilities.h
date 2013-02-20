/*=========================================================================

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 
=========================================================================*/
#ifndef SELECTIONUTILITIES_H
#define SELECTIONUTILITIES_H

#include <set>
#include <vector>
#include <string>
#include <sstream>
#include <map>

#include "vtkTable.h"
#include "vtkVariant.h"
#include "vtkIdTypeArray.h"
#include "vtkIntArray.h"
#include "vtkStringArray.h"

#include "vtkAnnotationLink.h"
#include "vtkSelection.h"
#include "vtkSelectionNode.h"

#include "vtkSmartPointer.h"

#include "ftkGUI/ObjectSelection.h"

namespace SelectionUtilities
{
	vtkSelection * ConvertIDsToVTKSelection(vtkIdTypeArray * vtkIDs);

	vtkIdTypeArray * ConvertVTKSelectionToIDArray(vtkSelection * inputSelection);

	std::set< vtkIdType > ObjectSelectionToIDSet(std::set<long int> curSel);

	void RemoveRowAndReMapTable(vtkIdType key, vtkSmartPointer<vtkTable> modTable, 
		std::map< vtkIdType, vtkIdType> TableIDMap);
	
}
#endif