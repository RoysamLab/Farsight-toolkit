
#ifndef DENDROGRAMSELECTION_H_
#define DENDROGRAMSELECTION_H_

#include <string>
#include <set>
#include <fstream>
//QT INCLUDES
#include <QtCore>
#include <QtGui>
#include "vtkTable.h"
#include "vtkSmartPointer.h"
#include "vtkDoubleArray.h"
#include "vtkAbstractArray.h"
#include "vtkVariantArray.h"
#include <ftkGUI/ObjectSelection.h>
#include "Dendrogram.h"
#include <vtkInteractorStyleRubberBand2D.h>
//#include <vtkIdTypeArray.h>

class Dendrogram_Selection : public QMainWindow
{
	Q_OBJECT;

public:
	Dendrogram_Selection(QWidget * parent = 0, Qt::WindowFlags flags = 0);
	~Dendrogram_Selection();
	
	ObjectSelection * GetObjectSelection();
	Dendrogram *dendro;
	/*void SelectByIDs(std::vector<int> IDs);*/
	std::set<long int> GetSelecectedIDs();
	vtkSmartPointer<vtkInteractorStyleRubberBand2D > Rubber_Band ;
	//static void SelectionCallbackFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData); ///Interactor callback
	//vtkSmartPointer<vtkCallbackCommand> selectionCallback;
	
	
protected slots:
	void updateSelection();
	
signals:
	void selection_Changed(void);

private:
	ObjectSelection * Selection;
	QVTKWidget *QVTK;
	
};
#endif
