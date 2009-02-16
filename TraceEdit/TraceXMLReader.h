#ifndef TRACEXMLREADER_H_
#define TRACEXMLREADER_H_

#include "itkObjectFactory.h"
#include "itkSmartPointer.h"
#include "itkMacro.h"
#include "itkLightObject.h"

#include <libxml/xmlreader.h>
#include <libxml/parser.h>
#include <libxml/tree.h>

#include "vtkPolyLine.h"
#include <vector>
#include <string>
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"

#include "vtkFloatArray.h"
#include "vtkAppendPolyData.h"

#include <vector>
#include <string>

class TraceXMLReader:public itk::LightObject
{
public:
	typedef TraceXMLReader Self;
	typedef  itk::SmartPointer< Self >  Pointer;
	typedef  itk::SmartPointer< const Self >  ConstPointer;
	itkNewMacro(Self);

	vtkPolyData* GetTraces();
	void SetFileName(std::string);
	void Read();

private:
	std::string m_XMLFileName;
	vtkPolyData *m_Traces;
	vtkFloatArray* m_lineScalars;
	vtkPoints *linePts ;
	vtkCellArray *lineCell ;

protected:
	TraceXMLReader();
	~TraceXMLReader();
};

#endif /*TRACEXMLREADER_H_*/
