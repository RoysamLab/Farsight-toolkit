#include <QVTKWidget.h>
#include "vtkDelimitedTextReader.h"
#include "vtkGraphLayoutStrategy.h"
#include "vtkGraphLayoutView.h"
#include "vtkTable.h"
#include "vtkTree.h"
#include "vtkTableToGraph.h"
#include "vtkTableToTreeFilter.h"
#include "vtkRenderWindow.h"
#include "vtkSmartPointer.h"


#include "vtkGroupLeafVertices.h"
#include "vtkTreeLayoutStrategy.h"

class FeatureRelation
{
public:
	FeatureRelation(void);
	~FeatureRelation(void);

	void FeatureGraph();
	QVTKWidget * mainQTRenderWidget;
};
