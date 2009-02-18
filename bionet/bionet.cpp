#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkDataSetAttributes.h"
#include <vtkSmartPointer.h>
#include "vtkActor.h"
#include "vtkGraphToPolyData.h"
#include "vtkLabeledDataMapper.h"
#include "vtkActor2D.h"
#include "vtkTextProperty.h"
#include "vtkPolyDataMapper.h"
//#include "vtkIntArray.h"
#include "vtkDoubleArray.h"
#include "vtkProperty.h"
#include "vtkVertexGlyphFilter.h"
#include "vtkVertexListIterator.h"
#include "vtkSphereSource.h"
// include libxml2
#include "libxml/parser.h"
#include <stdio.h>
#include <map>
#include <vector>
//Required for Graph operations
#include "vtkGraph.h"
#include "vtkMutableUndirectedGraph.h"
#include "vtkGraphToPolyData.h"
#include "vtkGraphWriter.h"
#include <vtkBoostKruskalMinimumSpanningTree.h>
#include <vtkBoostPrimMinimumSpanningTree.h>
#include "vtkExtractSelectedGraph.h"
#include <vtkAdjacentVertexIterator.h>
#include "vtkInEdgeIterator.h"
#include "vtkOutEdgeIterator.h"
#include <set> //for detecting cycles in the graph
// Required for implementing Remove vertex/edge
#include "vtkSelectionSource.h"
#include "vtkSelectionNode.h"
#include "vtkVertexDegree.h"

#include "vtkBMPReader.h"
#include "vtkLookupTable.h" // Required for defining colors
#include <vtkWindowToImageFilter.h>
#include <vtkCamera.h>
#include <vtkTIFFWriter.h>
#include "vtkGlyph3D.h"
#include "vtkGlyphSource2D.h"

#include <vtkRenderView.h> //Mouse drag and selections
#include <vtkCommand.h>

#include <iostream>
#include <fstream>


#define VTK_CREATE(type, name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

#define myVTKInt vtkSmartPointer<vtkIntArray>


using namespace std;
using namespace boost;

class TestRenderViewUpdater : public vtkCommand
{
public:
  static TestRenderViewUpdater* New()
  { return new TestRenderViewUpdater; }
  
  void AddView(vtkView* view)
  {
    this->Views.push_back(view);
    view->AddObserver(vtkCommand::SelectionChangedEvent, this);
  }
  
  virtual void Execute(vtkObject*, unsigned long, void*)
  {
    for (unsigned int i = 0; i < this->Views.size(); i++)
      {
      this->Views[i]->Update();
      }
  }
private:
  TestRenderViewUpdater() { }  
  ~TestRenderViewUpdater() { }
  vector<vtkView*> Views;
};


/**************************************************************************************/
// Sets color codes for biological objects. A color array for objects in a graph has   /
// to be constructed and attached to the graph by using the color codes here           /
/**************************************************************************************/


map<char,int> setColorCode() {
	map<char,int> channelColorsTable;
	// Returned Scalars will be mapped to a color during the display
	channelColorsTable['M']=1; //yellow for Microglia
	channelColorsTable['E']=2; //cyan for Endothelials
	channelColorsTable['A']=3; //Red for Astrocytes
	channelColorsTable['N']=4; //blue Astrocytes
							//5 if the vertex is a member of out-pyramidial region
	return channelColorsTable;
};

/**************************************************************************************/
//We need an associative sorted container to be able to relate vertex ids to           /
//node ids given in the xml file. vtk forces us to do it this way. We cannot use       /
//node ids as vertex ids in the graphs. Vertex ids have to start from 0                /
/**************************************************************************************/
map<int,int> vertexNodeLabelToID;	

/************************************************************************************/
// Store Edge Ids with source and target of the edge                                   /
// Will be used for edge deletion                                                      /
/**************************************************************************************/

map<pair<int,int>, vtkEdgeType> edgeSourgeTargetToID;
double avg=0; //average edge length

/**************************************************************************************/
// A helper-function that is used for detecting cycles. If the vertex v is already in  /
// vertices, it means there is a back-edge from it causing a cycle.
/**************************************************************************************/
bool cycleDetected(set<int>* vertices, vtkIdType v) {
	set<int>::iterator pos;
	pos=vertices->find(v);
	if (pos==vertices->end()) {
		vertices->insert(v);
		return false;
	} else
		return true;
}	

/**************************************************************************************/
// A helper-function that is used for computing averages for pyramidial region         /
// The result is the average length of edges outgoing from each node in the graph      /
// For undirected graphs (u,v) and (v,u) are the same edges. But in VTK they are not   /
// Note also that this function can be used to produce averages too for pyramidials        /
/**************************************************************************************/
void pyramidial(vtkGraph* g, float cutoff) {

  //numEdges = 0;
	int inNodes=0, outNodes=0;
  vtkIdType neighborCount =0;
  float subtotal=0, avg=0;

  ofstream myfile ("pyramidial.txt");

  /*
  if (myfile.is_open())
  {
    myfile << "This is NOT a line.\n";
    myfile << "This is another line.\n";
    myfile.close();
  }
  else cout << "Unable to open file";
*/


  VTK_CREATE(vtkVertexListIterator, vertices);
  VTK_CREATE(vtkOutEdgeIterator, outEdges);
  g->GetVertices(vertices);
  vtkOutEdgeType e;
  ///////////////////////////////////////////////////

  if (myfile.is_open())
  {
  //cout<<"Average Edge Lengths for pyramidial Cell Layer:"<<endl;
  while (vertices->HasNext())
    {
    vtkIdType v = vertices->Next();
    g->GetOutEdges(v, outEdges);
    //vtkIdType index = 0;
	neighborCount =0;
	subtotal = 0;
	//cout<<"source: "<<v<<endl;
	set<int>* targetvertices=new set<int>(); //will be used for taking care of cycles. 
							 //Edge weights will be the same in a cycle
	///////////////////////////////////////////////////
    while (outEdges->HasNext()) {
		e = outEdges->Next();
		if (!cycleDetected(targetvertices, e.Target)) {
			//cout<<"\t Target: "<<e.Target<<endl;
			//cout<<"Weight "<<g->GetEdgeData()->GetArray("EdgeWeights")->GetTuple1(e.Id)<<endl;
			subtotal+=g->GetEdgeData()->GetArray("EdgeWeights")->GetTuple1(e.Id);
			neighborCount++;
		}
	}//////////////////////////////////////////////////

	if (neighborCount > 0) {
		avg = subtotal/neighborCount;
		// Out pyramidial - 5 defines the color
		if (avg > cutoff) {g->GetVertexData()->GetArray("ChannelColors")->SetTuple1(v,5); 
		outNodes++;}
		else inNodes++;
		// Use the following line if the averages are required
		myfile <<avg<<endl;
		//cout<<avg<<endl<<endl;

	}
  }
  myfile.close();
  cout<<"Averages for Pyramidial Cell Layer has been written to Pyramidial.txt"<<endl<<endl;
  cout<<"Number of nodes in Pyramidial Cell Layer: "<<inNodes<<endl;  
  cout<<"Number of nodes located out of Pyramidial Cell Layer: "<<outNodes<<endl;

  } else cout << "Unable to open file to write averages for Pyramidial Cell Layer";

}

/*
void pyramidial(vtkGraph* g, double avg_edge_length) {
	
	//VTK_CREATE(vtkOutEdgeIterator, edges);
	VTK_CREATE(vtkEdgeListIterator, edges);
	g->GetEdges(edges);
	g->GetE
	
	while (edges->HasNext()) {
		
		vtkEdgeType e = edges->Next();
		if ((g->GetEdgeData()->GetArray("EdgeWeights")->GetTuple1(e.Id)) > avg) {
			//cout<<"iste avea : "<<g->GetEdgeData()->GetArray("EdgeWeights")->GetTuple1(e.Id)<<endl;
			//cout<<g->GetVertexData()->GetArray("ChannelColors")->GetTuple1(e.Source)<<endl;
			//cout<<"Source "<<e.Source<<endl;
			g->GetVertexData()->GetArray("ChannelColors")->SetTuple1(e.Source,5); 
				//cout<<g->GetVertexData()->GetArray("ChannelColors")->GetTuple1(e.Source)<<endl;
	
		}
				//5 means it is member of out pyramidial
	};
};
*/
// XML Parsing:Begin

bool readXML(char* graphFileName, float distance, vtkMutableUndirectedGraph* g) {
VTK_CREATE(vtkPoints,pts);
	xmlDocPtr doc;
	xmlNodePtr nodeLevel1;
	xmlNodePtr nodeLevel2;
	xmlChar* attr=NULL;
	char* nodeName=NULL;
	//Node attributes
		xmlChar *id=(xmlChar*)"id";
		xmlChar* label=(xmlChar*)"label";
		xmlChar *weight=(xmlChar*)"weight";
		//Attributes for graphics element
				xmlChar *type=(xmlChar*)"type"; //ex:circle
				//coordinates of nodes
				xmlChar* X=(xmlChar*)"x"; 
				xmlChar* Y=(xmlChar*)"y";
				xmlChar* Z=(xmlChar*)"z";

	//Attributes for Edge elements
		xmlChar *source=(xmlChar*)"source";
		xmlChar *target=(xmlChar*)"target";
	char* cnodeValue, *csource, *ctarget, *cweight;
	vtkIdType iID,iX,iY,iZ; // used for graph coordinates
	vtkIdType isource, itarget; //used for defining graph edges
	double iweight; 
	int e=0;    
	vtkIdType vertexID; 
	vtkEdgeType edgeID;	
	vtkIntArray* labels = vtkIntArray::New();
	labels->SetName("Label");
	
	map<char,int> channelColors; 
	//map<char,vtkDoubleArray*> channelColors; //used for finding color codes for channels
	vtkIntArray* colorCodeArr = vtkIntArray::New();
	//vtkDoubleArray* colorCodeArr = vtkDoubleArray::New();
	colorCodeArr->SetName("ChannelColors");

	vtkIntArray* edgeColors = vtkIntArray::New();
	edgeColors->SetName("EdgeColors");

	// This is used for finding minimum-spanning trees
	vtkDoubleArray* edgeWeightsArr = vtkDoubleArray::New();
	edgeWeightsArr->SetName("EdgeWeights");

	vtkIdType i=0; //index for vertex colors
	//Set the colors: 0 if Neurons, 1 if Microglia, 
	vector<char> firstLetterOfLabels;
	pair<int,int> key;
	//double avg=0;
		int mic=0;
	channelColors = setColorCode();
	doc = xmlParseFile(graphFileName);
	for(	nodeLevel1 = doc->children;
		nodeLevel1 != NULL;
		nodeLevel1 = nodeLevel1->next)
	{

		nodeName = (char*)nodeLevel1->name;
		//Check if this is a graph
		//begin:outer-if

		if (strcmp(nodeName,"graph") == 0) {
			for(nodeLevel2 = nodeLevel1->children;
				nodeLevel2 != NULL;
				nodeLevel2 = nodeLevel2->next)
			{
				nodeName = (char*)nodeLevel2->name;
				//Check if this is a node
				//begin:inner-if
				if (strcmp(nodeName,"node")== 0) {
					cnodeValue = (char *)(xmlGetProp(nodeLevel2,id));
					iID = atoi(cnodeValue);		
					cnodeValue = (char *)(xmlGetProp(nodeLevel2->children->next,X));
					iX = atoi(cnodeValue);
					cnodeValue = (char *)(xmlGetProp(nodeLevel2->children->next,Y));
					iY = atoi(cnodeValue);
					cnodeValue = (char *)(xmlGetProp(nodeLevel2->children->next,Z));
					iZ = atoi(cnodeValue);
					
					//Add one vertex for this point
					vertexID=g->AddVertex();
					vertexNodeLabelToID[iID] = vertexID;
					
					//Show each node id in the graph
					labels->InsertValue(vertexID, iID);
					//Set the colors of vertices
					colorCodeArr->InsertValue(i,channelColors[((char *)(xmlGetProp(nodeLevel2,label)))[0]]);
					i++;
					pts->InsertPoint(vertexID,iX,iY,iZ);
					if (channelColors[((char *)(xmlGetProp(nodeLevel2,label)))[0]] == 1) {
						//cout<<"burdaaaaa "<<endl;
						mic++;
					}

					/*
					if (vertexID==1) cout<<vertexID<<" Label: "<<iID<<endl;
					if (vertexID==26) cout<<vertexID<<" Label: "<<iID<<endl;
					if (vertexID==38) cout<<vertexID<<" Label: "<<iID<<endl;
					if (vertexID==18) cout<<vertexID<<" Label: "<<iID<<endl;
					if (vertexID==30) cout<<vertexID<<" Label: "<<iID<<endl;
					if (vertexID==1) cout<<vertexID<<" Label: "<<iID<<endl;
					if (vertexID==18) cout<<vertexID<<" Label: "<<iID<<endl;
					if (vertexID==26) cout<<vertexID<<" Label: "<<iID<<endl;
					if (vertexID==38) cout<<vertexID<<" Label: "<<iID<<endl;
					if (vertexID==48) cout<<vertexID<<" Label: "<<iID<<endl;
					*/
				} else if (strcmp(nodeName,"edge")== 0) {
					if (e == 0) g->SetPoints(pts);

					csource =(char *)xmlGetProp(nodeLevel2,source);
					isource = atoi(csource);
					ctarget =(char *)xmlGetProp(nodeLevel2,target);
					itarget = atoi(ctarget);
					cweight =(char *)xmlGetProp(nodeLevel2,weight);
					iweight = atof(cweight);

					// If a distance is given we need to select the edges that are shorter 
					// than the distance	
					if (distance == -1) {//distance is not given
						edgeID=g->AddEdge(vertexNodeLabelToID[isource],vertexNodeLabelToID[itarget]);
						edgeWeightsArr->InsertNextValue(iweight);	
						key= pair<int,int>(vertexNodeLabelToID[isource], vertexNodeLabelToID[itarget]);
						edgeSourgeTargetToID[key]=edgeID;
						edgeColors->InsertValue(e,0);
						e++;
					} else if (iweight < distance) {
							//Add this edge to the graph	
							edgeID=g->AddEdge(vertexNodeLabelToID[isource],vertexNodeLabelToID[itarget]);
							edgeWeightsArr->InsertNextValue(iweight);	
							key= pair<int,int>(vertexNodeLabelToID[isource], vertexNodeLabelToID[itarget]);
							edgeSourgeTargetToID[key]=edgeID;
							edgeColors->InsertValue(e,0);
							e++;
					}
					
				} else if (strcmp(nodeName,"text")== 0) 
							{ //Do Nothing 
							}
				   else {
						cout<<"Unrecognized Node Element!"<<nodeName<<endl;
						return false;
					} //end:inner-if
			}//for
		} else {cout<<"Incorrect Graph File Format!"<<nodeName<<endl;
				return false;
		}
	} //for

	//if (e > 0) avg = avg/e;  //Average edge length for pyramidial region computation
	//cout<<"iste averaj "<<avg<<endl;
	//xmlSaveFile("xmlfile_copy.xml", doc);
	xmlFreeDoc(doc);
	//Ad vertex ids
	g->GetVertexData()->AddArray(labels);
	//Add object colors
	g->GetVertexData()->AddArray(colorCodeArr);
	//Ad edge weights
	g->GetEdgeData()->AddArray(edgeWeightsArr);
	//Add Edge colors
	g->GetEdgeData()->AddArray(edgeColors);

	return true;
};


/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: TestGraphAlgorithms.cxx,v $


  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.
 

=========================================================================*/
/*-------------------------------------------------------------------------
  Copyright 2008 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
  the U.S. Government retains certain rights in this software.
-------------------------------------------------------------------------*/

//PerformAlgorithm(ren, degree,0, 0, "ChannelColors", 0, 4,"EdgeColors");
void PerformAlgorithm(vtkRenderer* ren, vtkAlgorithm* alg, 
  double xoffset, double yoffset, 
  const char* vertColorArray, double vertMin, double vertMax, 
  const char* edgeColorArray, double edgeMin = 0, double edgeMax = 0)
{
  VTK_CREATE(vtkGraphToPolyData, graphToPoly);
  graphToPoly->SetInputConnection(alg->GetOutputPort());

  VTK_CREATE(vtkGlyphSource2D, glyph);
  glyph->SetGlyphTypeToVertex();
  
/*  If we want vertices as only 2D circles
	glyph->SetGlyphTypeToCircle();
	glyph->FilledOn();
	glyph->CrossOff();
	glyph->SetScale(25.0);
*/

  VTK_CREATE(vtkGlyph3D, vertexGlyph);
  vertexGlyph->SetInputConnection(0, graphToPoly->GetOutputPort());
  vertexGlyph->SetInputConnection(1, glyph->GetOutputPort());
  
  // Vertices as 3D spheres 
  VTK_CREATE(vtkSphereSource, sphere);
  sphere->SetThetaResolution(8); sphere->SetPhiResolution(8);
  sphere->SetRadius(4);
  vertexGlyph->SetInputConnection(1, sphere->GetOutputPort());


  VTK_CREATE(vtkPolyDataMapper, vertexMapper);
  vertexMapper->SetInputConnection(vertexGlyph->GetOutputPort());
  vertexMapper->SetScalarModeToUsePointFieldData();

  // Create a color map and color the vertices
  vtkLookupTable * lut = vtkLookupTable::New();
  lut->SetTableRange(0, 5);
  lut->SetNumberOfColors(6);
  lut->Build();
  lut->SetTableValue(0, 1.0, 1.0, 1.0, 1); //Edge color white
		//lut->SetTableValue(0, 0.0, 0.0, 0.0, 1); //Edge color black
		//lut->SetTableValue(0, 0.0, 0.0, 1.0, 1); //Edge color blue
  lut->SetTableValue(1, 1.0, 1.0, 0.0, 1);  //yellow
  lut->SetTableValue(2,  0.4980, 1.0000, 0.8314 , 1); //cyan
  lut->SetTableValue(3, 1.0, 0.0, 0.0, 1); //Red
  lut->SetTableValue(4, 0.6275, 0.1255, 0.9412, 1); //purple
 //lut->SetTableValue(5, 0.5765, 0.4392, 0.8588, 1); //medium-purple for pyramidial
 lut->SetTableValue(5, 1.0, 1.0, 1.0, 1);
  if (vertColorArray)
    {
	    vertexMapper->SetLookupTable(lut);
		//mapper->SetScalarModeToUsePointFieldData();
		vertexMapper->SetColorModeToMapScalars();
		//mapper->SelectColorArray("Color"); 
		vertexMapper->SelectColorArray(vertColorArray);
		vertexMapper->SetScalarRange(vertMin, vertMax);	
    }

  // Set vertex shape to circle:
  //vertexMapper->


  VTK_CREATE(vtkActor, vertexActor);
  vertexActor->GetProperty()->SetPointSize(10.0);
  vertexActor->SetPosition(xoffset, yoffset, 0.001);
  //vertexActor->GetPickable()


//Load in the texture map. A texture is any unsigned char image. If it
//is not of this type, you will have to map it through a lookup table
//or by using vtkImageShiftScale.

	//VTK_CREATE(vtkBMPReader,bmpReader);
	//bmpReader->SetFileName("p1.bmp");
	//bmpReader->Update();
	//VTK_CREATE(vtkTexture, atext);
	//atext->SetInputConnection(bmpReader->GetOutputPort());
	//atext->InterpolateOn();
	
	//vertexActor->SetTexture(atext);
 vertexActor->SetMapper(vertexMapper);


  VTK_CREATE(vtkPolyDataMapper, edgeMapper);
  edgeMapper->SetInputConnection(graphToPoly->GetOutputPort());
  edgeMapper->SetScalarModeToUseCellFieldData();
  if (edgeColorArray)
    {
		edgeMapper->SetLookupTable(lut);
		//mapper->SetScalarModeToUsePointFieldData();
		edgeMapper->SetColorModeToMapScalars();

    edgeMapper->SelectColorArray(edgeColorArray);
    edgeMapper->SetScalarRange(edgeMin, edgeMax);
    }
  VTK_CREATE(vtkActor, edgeActor);
  edgeActor->SetMapper(edgeMapper);
  edgeActor->SetPosition(xoffset, yoffset, 0);

  ren->AddActor(vertexActor);
  ren->AddActor(edgeActor);
}

// vtk does not provide edge remove.
// This does not delete out edges yet!
vtkGraph* RemoveEdge(vtkGraph* original_graph, vtkIdType u,vtkIdType v)
{
	// first find the edge id using source and target vertices.
	vtkEdgeType id_of_edge_to_remove;
	map<pair<int,int>, vtkEdgeType>::iterator ste;
	ste = edgeSourgeTargetToID.find(pair<int,int>((int)u,(int)v));
	if (ste != edgeSourgeTargetToID.end())
		id_of_edge_to_remove= ste->second;
	else {
			cout<<"No Edge from "<<u<<" to "<<v<<endl;
			return original_graph;
	};

	//An edgeId has been found. Delete the edge
	vtkSelectionSource* ss = vtkSelectionSource::New();
	ss->AddID(-1, id_of_edge_to_remove.Id);
	ss->SetFieldType(vtkSelectionNode::EDGE);
	ss->SetContentType(vtkSelectionNode::INDICES);
	ss->SetInverse(1);

	vtkExtractSelectedGraph* e = vtkExtractSelectedGraph::New();
	e->SetInput(0, original_graph);
	e->SetInputConnection(1, ss->GetOutputPort());
	e->Update();

	vtkGraph* subgraph = e->GetOutput();
	return subgraph;

};

/**************************************************************************************/
// vtk does not have RemoveEdge.
/**************************************************************************************/
template<class EdgeType>
vtkGraph* RemoveEdge(vtkGraph* original_graph, EdgeType id_of_edge_to_remove)
{	
	vtkSelectionSource* ss = vtkSelectionSource::New();
	ss->AddID(-1, id_of_edge_to_remove.Id);
	ss->SetFieldType(vtkSelectionNode::EDGE);
	ss->SetContentType(vtkSelectionNode::INDICES);
	ss->SetInverse(1);

	vtkExtractSelectedGraph* e = vtkExtractSelectedGraph::New();
	e->SetInput(0, original_graph);
	e->SetInputConnection(1, ss->GetOutputPort());
	e->Update();

	vtkGraph* subgraph = e->GetOutput();
	return subgraph;

};

/**************************************************************************************/
// vtk does not have RemoveVertex
/**************************************************************************************/
vtkGraph* RemoveVertex(vtkGraph* original_graph, vtkIdType id_of_vertex_to_remove)
{

	// We first need to check if the vertex to be deleted has any adjacent vertices.
	// If there are, the edges betwen this and the adjacent vertices have to be deleted

	/*
	vtkIdType v;
	VTK_CREATE(vtkAdjacentVertexIterator, adjacent);
	original_graph->GetAdjacentVertices(id_of_vertex_to_remove,adjacent);
	while (adjacent->HasNext()) {
		v = adjacent->Next();
		//Remove Edge(id_of_vertex_to_remove, v)
		original_graph=RemoveEdge(original_graph,id_of_vertex_to_remove,v);
	};
*/


	// First delete Out edges
	VTK_CREATE(vtkOutEdgeIterator, outEdges);
	original_graph->GetOutEdges(id_of_vertex_to_remove, outEdges);
	//vtkOutEdgeType out_edge_to_be_deleted;
	while (outEdges->HasNext()) {
		vtkOutEdgeType out_edge_to_be_deleted = outEdges->Next();
		original_graph = RemoveEdge(original_graph, out_edge_to_be_deleted);
	};

	// Now delete In edges
    VTK_CREATE(vtkInEdgeIterator, inEdges);
    vtkInEdgeType in_edge_to_be_deleted;
	original_graph->GetInEdges(id_of_vertex_to_remove, inEdges);
	while (inEdges->HasNext()) {
		in_edge_to_be_deleted = inEdges->Next();
		original_graph = RemoveEdge(original_graph, in_edge_to_be_deleted);
	};

	// Now delete the vertex

	vtkSelectionSource* ss = vtkSelectionSource::New();
	ss->AddID(-1, id_of_vertex_to_remove);
	ss->SetFieldType(vtkSelectionNode::VERTEX);
	ss->SetContentType(vtkSelectionNode::INDICES);
	ss->SetInverse(1);

	vtkExtractSelectedGraph* e = vtkExtractSelectedGraph::New();
	e->SetInput(0, original_graph);
	e->SetInputConnection(1, ss->GetOutputPort());
	e->Update();

	vtkGraph* subgraph = e->GetOutput();
	return subgraph;

};

// Records the displays in order to get an animation
void recordImage (vtkRenderer* ren, vtkRenderWindow* win) {

	for(int counter = 0; counter< 1; counter++)
	{
		//ren->GetActiveCamera()->Azimuth();
		ren->GetActiveCamera()->Roll(-90);
		ren->ResetCameraClippingRange();
		ren->Render();
		vtkWindowToImageFilter* img = vtkWindowToImageFilter::New();
		img->SetInput(win);
		vtkTIFFWriter *writer = vtkTIFFWriter::New();
		writer->SetInput(img->GetOutput());
		char buff[1024];
		sprintf(buff,"C:\\Render\\render_%d.tif",counter);
		writer->SetFileName(buff);
		writer->Update();
	}

	/*
   for(int counter = 0; counter< 37; counter++)
	{
		ren->GetActiveCamera()->Zoom(1.01);
		//ren->GetActiveCamera()->Roll(2);
		ren->ResetCameraClippingRange();
		ren->Render();
		vtkWindowToImageFilter* img = vtkWindowToImageFilter::New();
		img->SetInput(win);
		vtkTIFFWriter *writer = vtkTIFFWriter::New();
		writer->SetInput(img->GetOutput());
		char buff[1024];
		sprintf(buff,"C:\\Render\\render0_%d.tif",counter);
		writer->SetFileName(buff);
		writer->Update();
	}

	for(int counter = 0; counter< 30; counter++)
	{
		ren->GetActiveCamera()->Azimuth(1);
		ren->GetActiveCamera()->Roll(-5);
		ren->ResetCameraClippingRange();
		ren->Render();
		vtkWindowToImageFilter* img = vtkWindowToImageFilter::New();
		img->SetInput(win);
		vtkTIFFWriter *writer = vtkTIFFWriter::New();
		writer->SetInput(img->GetOutput());
		char buff[1024];
		sprintf(buff,"C:\\Render\\render1_%d.tif",counter);
		writer->SetFileName(buff);
		writer->Update();
	}
	

	for(int counter = 0; counter<45; counter++)
	{
		ren->GetActiveCamera()->Azimuth(1);
		ren->ResetCameraClippingRange();
		ren->Render();
		vtkWindowToImageFilter* img = vtkWindowToImageFilter::New();
		img->SetInput(win);
		vtkTIFFWriter *writer = vtkTIFFWriter::New();
		writer->SetInput(img->GetOutput());
		char buff[1024];
		sprintf(buff,"C:\\Render\\render2_%d.tif",counter);
		writer->SetFileName(buff);
		writer->Update();
	}

	for(int counter = 0; counter< 35; counter++)
	{
		ren->GetActiveCamera()->Zoom(1.01);
		ren->GetActiveCamera()->Roll(1);
		ren->ResetCameraClippingRange();
		ren->Render();
		vtkWindowToImageFilter* img = vtkWindowToImageFilter::New();
		img->SetInput(win);
		vtkTIFFWriter *writer = vtkTIFFWriter::New();
		writer->SetInput(img->GetOutput());
		char buff[1024];
		sprintf(buff,"C:\\Render\\render3_%d.tif",counter);
		writer->SetFileName(buff);
		writer->Update();
	}

	for(int counter = 0; counter< 7; counter++)
	{
		ren->GetActiveCamera()->Zoom(0.9);
		//ren->GetActiveCamera()->Roll(2);
		ren->ResetCameraClippingRange();
		ren->Render();
		vtkWindowToImageFilter* img = vtkWindowToImageFilter::New();
		img->SetInput(win);
		vtkTIFFWriter *writer = vtkTIFFWriter::New();
		writer->SetInput(img->GetOutput());
		char buff[1024];
		sprintf(buff,"C:\\Render\\render4_%d.tif",counter);
		writer->SetFileName(buff);
		writer->Update();
	}

	for(int counter = 0; counter<90; counter++)
	{
		ren->GetActiveCamera()->Azimuth(1);
		ren->ResetCameraClippingRange();
		ren->Render();
		vtkWindowToImageFilter* img = vtkWindowToImageFilter::New();
		img->SetInput(win);
		vtkTIFFWriter *writer = vtkTIFFWriter::New();
		writer->SetInput(img->GetOutput());
		char buff[1024];
		sprintf(buff,"C:\\Render\\render5_%d.tif",counter);
		writer->SetFileName(buff);
		writer->Update();
	}

*/
	/*
	for(int counter = 0; counter< 100; counter++)
	{
		//iren->FlyTo(ren,coords[counter][0],coords[counter][1],coords[counter][2]);
		//ren->GetActiveCamera()->Azimuth(2);
		ren->GetActiveCamera()->Roll(0);
		ren->GetActiveCamera()->Azimuth(0);
		//ren->GetActiveCamera()->Zoom(1.01);
		ren->ResetCameraClippingRange();
		ren->Render();
		vtkWindowToImageFilter* img = vtkWindowToImageFilter::New();
		img->SetInput(win);
		vtkTIFFWriter *writer = vtkTIFFWriter::New();
		writer->SetInput(img->GetOutput());
		char buff[1024];
		sprintf(buff,"C:\\Render\\render_%d.tif",counter);
		writer->SetFileName(buff);
		writer->Update();
	}
*/
}

int main(int argc, char *argv[])
{
	cout<<endl<<endl;
	if (argc <4) {
			cout<<"Incorrect Usage! Check Parameters" <<endl;
			return 0;
	};
	VTK_CREATE(vtkRenderer, ren);
	vtkMutableUndirectedGraph* g=vtkMutableUndirectedGraph::New();
	//Read the vertices and edges into graph g
	//argv1= The name of the input network
	// Python Arguments
	//argv2= distance 
	//argv3= an integer that shows either the graph is standard or a minimum spanning tree
	//argv4= cutoff for pyramidial
	float distance=atof(argv[2]);
	int stdMst=atoi(argv[3]);
	float cutoff=atof(argv[4]);

	//cout<<"Here is what you sent me distance:"<<distance<<" stdMst: "<<stdMst<<" cutoff: "<<cutoff<<endl;

	readXML(argv[1],distance,g);

 	if (stdMst == 2) pyramidial(g,cutoff);
	if (stdMst == 1) {
		VTK_CREATE(vtkBoostKruskalMinimumSpanningTree,mst);
		mst->SetInput(g);
		mst->SetEdgeWeightArrayName("EdgeWeights");
		mst->Update();

		/*
		Kruskal's MST outputs a vtkSelection, which selects the edges that 
		are in the MST. This selection itself cannot be displayed in a graph 
		view, it is simply a list of edge indices. We need to use 
		vtkExtractSelectedGraph, with the original graph as the first 
		input and the vtkSelection output of MST as the second input. This will 
		extract the MST from the graph. The output of vtkExtractSelectedGraph 
		can be sent to the graph view.
		*/
		VTK_CREATE(vtkExtractSelectedGraph, mstGraph);
		mstGraph->SetInput(0,g);
		mstGraph->SetInput(1,mst->GetOutput());
		PerformAlgorithm(ren, mstGraph, 0, 0, "ChannelColors", 0, 5,"EdgeColors");
	} else {
		// Test vertex degree
		VTK_CREATE(vtkVertexDegree, degree);
		degree->SetInput(g);
		PerformAlgorithm(ren, degree,0, 0, "ChannelColors", 0, 5,"EdgeColors");
	}

	//SET VERTEX LABELS
	vtkGraphToPolyData* poly = vtkGraphToPolyData::New();
	poly->SetInput(g);
	VTK_CREATE(vtkLabeledDataMapper, labelMapper);
	labelMapper->SetInputConnection( poly->GetOutputPort() );
	labelMapper->SetLabelModeToLabelFieldData();
	labelMapper->SetInputArrayToProcess( 0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "id" );
	labelMapper->GetLabelTextProperty()->SetColor(1, 1, 1);
	VTK_CREATE(vtkActor2D, labelActor);
	labelActor->SetMapper( labelMapper );	
	labelMapper->GetLabelTextProperty()->SetFontFamilyToCourier();
	labelMapper->GetLabelTextProperty()->SetFontSize(10);
	labelMapper->GetLabelTextProperty()->SetShadow(0);
	ren->AddActor(labelActor );

	//labelMapper->GetLabelTextProperty()->SetFontSize(10);



	//Update the Graph interactively
	//Delete a vertex
	// g=vtkMutableUndirectedGraph::SafeDownCast(RemoveVertex(g,394));
	//Delete an edge
	//g=vtkMutableUndirectedGraph::SafeDownCast(RemoveEdge(g,213,387));
	// This one has many edges
	// g=vtkMutableUndirectedGraph::SafeDownCast(RemoveVertex(g,213));

	//Set the background color - olive
	//ren->SetBackground(0.2300, 0.3700, 0.1700 );
	//We need purple for icons
	//ren->SetBackground(0.6275, 0.1255, 0.9412); 

	//Set the background color - white
	VTK_CREATE(TestRenderViewUpdater, updater);
	ren->SetBackground(0, 0, 0 );
	//grey
	//ren->SetBackground(0.7529, 0.7529, 0.7529);
	//ren->SetBackground(0.8275, 0.8275, 0.8275);
	VTK_CREATE(vtkRenderWindowInteractor, iren);
	VTK_CREATE(vtkRenderWindow, win);
	win->SetSize(1200,1200);	
	win->AddRenderer( ren );
	win->SetInteractor( iren );
	//win->Render();

	/*
	VTK_CREATE(vtkRenderView, view);
  view->SetupRenderWindow(win);
updater->AddView(view);
	
  view->GetRenderer()->ResetCamera();
  view->Update();
*/
	iren->Initialize();
	win->Render();
	recordImage (ren, win);
	iren->Start();

//	iren->FlyTo(ren,100,100,100);

	//Clean up
	/*
	g->Delete();
	iren->Delete();
	win->Delete();
	ren->Delete();
	labelActor->Delete();
	labelMapper->Delete();
	poly->Delete();
	*/
	return 0;


}
