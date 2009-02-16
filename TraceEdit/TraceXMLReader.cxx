#include "TraceXMLReader.h"
#define MY_ENCODING "ISO-8859-1"
/*
TraceXMLReader v2
Trace XML Reader parses an XML file for coordinate positions of points along a line trace
the return type is in polydata that can be viewed in 3d
*/

TraceXMLReader::TraceXMLReader()	{
}

TraceXMLReader::~TraceXMLReader()	{

	//not initalized yet will clean up memory after xml read
}

vtkPolyData* TraceXMLReader::GetTraces()	{
	//std::cout<< "traces data /n" << *m_Traces <<std::endl;
	m_Traces->SetPoints(linePts);
	m_Traces->SetLines(lineCell);
	m_Traces->GetPointData()->SetScalars(m_lineScalars);
	std::cout<< "m_traces set"<< std::endl;
	return m_Traces;
	//polyData container m_Traces filled with the points and cell lines
	//return polydata when GetTraces called
}
void TraceXMLReader::SetFileName(std::string s)	{
	m_XMLFileName = s;
}
void TraceXMLReader::Read() 	{
	//reads the XML file 
	int count = 0;
/*test if libxml works and if xml is parsed*/
	LIBXML_TEST_VERSION;
	xmlDocPtr doc;
	doc = xmlReadFile(m_XMLFileName.c_str(), NULL, 0);
	if (doc == NULL)
	{
		fprintf(stderr, "Failed to parse %s\n", m_XMLFileName.c_str());
		return;
	}
//finds the root of the xml tree
	xmlNodePtr root_node = NULL, trace_node = NULL, line_node = NULL, bit_node = NULL;
	root_node = xmlDocGetRootElement(doc);
	std::cout << "Root node: " << root_node->name <<std::endl;
	
	/*	mtraces is the container linepts and cell added each time*/
	m_Traces = vtkPolyData::New();			
	linePts = vtkPoints::New();
	lineCell = vtkCellArray::New();
	m_lineScalars = vtkFloatArray::New();

	line_node =  root_node->children;
	//reading lines
	for ( ; line_node; line_node = line_node->next)
	{		
		if (line_node->type == XML_ELEMENT_NODE && !xmlStrcmp(line_node->name, BAD_CAST "TraceLine"))		
		{			
			float lid = atof((char*)xmlGetProp(line_node, BAD_CAST "ID"));
			
			//std::cout << "Line node: " << line_node->name << " " <<lid << std::endl;
			int counter= 0, pid =0, ret; double cur[3] , x, y, z;
			bit_node = line_node->children;
			//reading points
			for ( ;bit_node; bit_node = bit_node->next)
			{
				if ( !xmlStrcmp(bit_node->name, BAD_CAST "TraceBit") )	
				{
					pid = atof((char*)xmlGetProp(bit_node, BAD_CAST"ID"));
				    x = atof((char*)xmlGetProp(bit_node, BAD_CAST"X"))/2;
				    y = atof((char*)xmlGetProp(bit_node, BAD_CAST"Y"))/2;
				    z = atof((char*)xmlGetProp(bit_node, BAD_CAST"Z"));
 
					//printf("pid %d\n",pid);
					cur[0] = x; cur[1] = y; cur[2] = z;

					counter++;			
					/*	points and line color set; ret is point id*/
					ret = linePts->InsertNextPoint(cur);	
					m_lineScalars->SetNumberOfComponents(1);	
					m_lineScalars->InsertNextTuple1(1/(lid/10));
					if (counter >1)
					{		//inserting new line between points
						lineCell->InsertNextCell(2);
						lineCell->InsertCellPoint(ret -1);
						lineCell->InsertCellPoint(ret);
						//std::cout<< " add cell to line: " << lid <<std::endl;
					}
				}//end if tracebit 				
			}//end	for bit_node			
		}//end if line_node
	}	
}
