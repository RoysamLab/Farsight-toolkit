#include "tissuenets.h"

// Description:
// Sets color codes for biological objects. A color array for objects in a graph has   /
// to be constructed and attached to the graph by using the color codes here           /
map<char,int> BioNet::SetColorCode() {
	map<char,int> channelColorsTable;
	// Returned Scalars will be mapped to a color during the display
	channelColorsTable['M']=1; //yellow for Microglia
	channelColorsTable['E']=2; //cyan for Endothelials
	channelColorsTable['A']=3; //Red for Astrocytes
	channelColorsTable['N']=4; //Purple for Neurons
	channelColorsTable['1']=3; //1 has been used for Astrocytes
	//channelColorsTable['W']=5;  // white for selected vertices						
		//5 if the vertex is a member of out-pyramidial region
	return channelColorsTable;
};

// Description:
// A helper-function that is used for detecting cycles. If the vertex v is already in  /
// vertices, it means there is a back-edge from it causing a cycle.
bool BioNet::CycleDetected(set<int>* vertices, vtkIdType v) {
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
void BioNet::Averages() {

  //numEdges = 0;
	int inNodes=0, outNodes=0;
  vtkIdType neighborCount =0;
  float subtotal=0, avg=0;

  ofstream myfile ("pyramidal.txt");

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
		if (!CycleDetected(targetvertices, e.Target)) {
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
  cout<<"Averages for Pyramidal Cell Layer has been written to Pyramidal.txt"<<endl<<endl;
  cout<<"Number of nodes in Pyramidal Cell Layer: "<<inNodes<<endl;  
  cout<<"Number of nodes located out of Pyramidal Cell Layer: "<<outNodes<<endl;

  } else cout << "Unable to open file to write averages for Pyramidal Cell Layer";

}

// Description:
// Constructs a graph from a given XGMML file
bool BioNet::ReadXGMML(char* graphFileName, float n) {
	VTK_CREATE(vtkPoints,pts);

   TiXmlDocument doc(graphFileName);
   doc.LoadFile();
   TiXmlHandle docHandle( &doc );
   TiXmlElement* levelOneElement =
   docHandle.FirstChild("graph").Element();
   TiXmlElement *levelTwoElement, *levelThreeElement;

   
//	xmlDocPtr doc;
//	xmlNodePtr nodeLevel1;
//	xmlNodePtr nodeLevel2;
//	xmlChar* attr=NULL;
	const char *nodeName, *label;
	//Node attributes
//		xmlChar *id=(xmlChar*)"id";
//		xmlChar* label=(xmlChar*)"label";
//		xmlChar *weight=(xmlChar*)"weight";
		//Attributes for graphics element
//				xmlChar *type=(xmlChar*)"type"; //ex:circle
//				//coordinates of nodes
//				xmlChar* X=(xmlChar*)"x"; 
//				xmlChar* Y=(xmlChar*)"y";
//				xmlChar* Z=(xmlChar*)"z";


	//Attributes for Edge elements
//		xmlChar *source=(xmlChar*)"source";
//		xmlChar *target=(xmlChar*)"target";	
	vtkIdType iID,iX,iY,iZ; // used for graph coordinates
	vtkIdType isource, itarget; //used for defining graph edges
	double iweight; 
	int e=0;    
	vtkIdType vertexID; 
	vtkEdgeType edgeID;	
	vtkIntArray* labels = vtkIntArray::New();
	labels->SetName("Label");

	//vtkDoubleArray* edgeIDs = vtkDoubleArray::New();
	//edgeIDs->SetName("EdgeIDs");

	//////////////////
	vtkStringArray* edgeIDs = vtkStringArray::New();
	edgeIDs->SetName("EdgeIDs");
	vtkStringArray* edgeIDs2 = vtkStringArray::New();
	edgeIDs->SetName("EdgeIDs2");

	vtkDoubleArray* vertexIDs = vtkDoubleArray::New();
	vertexIDs->SetName("vertexIDs");

	vtkDoubleArray* vertexIDs2 = vtkDoubleArray::New();
	vertexIDs->SetName("vertexIDs2");

	//vtkEdgeType* edgeIDs = new vtkEdgeType();

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
	channelColors = SetColorCode();
	//doc = xmlParseFile(graphFileName);
   while(levelOneElement) {
	   nodeName = levelOneElement->Value();	   

	//for(	nodeLevel1 = doc->children;
	//	nodeLevel1 != NULL;
	//	nodeLevel1 = nodeLevel1->next)
	//{
	//	nodeName = (char*)nodeLevel1->name;
		//Check if this is a graph
		//begin:outer-if
	   //Check if this is a graph
    if (strcmp(nodeName,"graph") == 0)
      {
      levelTwoElement = levelOneElement->FirstChildElement();
      while(levelTwoElement)
        {
        nodeName = (char*)levelTwoElement->Value();
        //Check if this is a node
        if (strcmp(nodeName,"node") == 0)			
          {
          if(levelTwoElement->QueryIntAttribute("id", (int *)&iID) != TIXML_SUCCESS)
            {
            cerr << "ERROR: encountered a node with no ID" << endl;
            return false;
            }
          levelThreeElement = levelTwoElement->FirstChildElement();
          if(!levelThreeElement)
            {
            cerr << "ERROR: encountered a node with no child element" << endl;
            }
          if(levelThreeElement->QueryIntAttribute("x", (int *)&iX) != TIXML_SUCCESS)
            {
            cerr << "ERROR: graphics element with no X value" << endl;
            return false;
            }
          if(levelThreeElement->QueryIntAttribute("y", (int *)&iY) != TIXML_SUCCESS)
            {
            cerr << "ERROR: graphics element with no Y value" << endl;
            return false;
            }
          if(levelThreeElement->QueryIntAttribute("z", (int *)&iZ) != TIXML_SUCCESS)
            {
            cerr << "ERROR: graphics element with no Z value" << endl;
            return false;
            }


//		if (strcmp(nodeName,"graph") == 0) {
//			for(nodeLevel2 = nodeLevel1->children;
//				nodeLevel2 != NULL;
//				nodeLevel2 = nodeLevel2->next)
//			{
//				nodeName = (char*)nodeLevel2->name;
				//Check if this is a node
				//begin:inner-if
				//if (strcmp(nodeName,"node")== 0) {
				//	cnodeValue = (char *)(xmlGetProp(nodeLevel2,id));
				//	iID = atoi(cnodeValue);		
				//	cnodeValue = (char *)(xmlGetProp(nodeLevel2->children->next,X));
				//	iX = atoi(cnodeValue);
				//	cnodeValue = (char *)(xmlGetProp(nodeLevel2->children->next,Y));
				//	iY = atoi(cnodeValue);
				//	cnodeValue = (char *)(xmlGetProp(nodeLevel2->children->next,Z));
				//	iZ = atoi(cnodeValue);
			
					
					//Add one vertex for this point
					vertexID=g->AddVertex();
					vertexNodeLabelToID[iID] = vertexID;
					vertexIDs->InsertNextValue(i);
					vertexIDs2->InsertNextValue(i);
					
					//Show each node id in the graph
					labels->InsertValue(vertexID, iID);
		            label = levelTwoElement->Attribute("label");
					if(label == NULL){
						cerr << "ERROR: encountered a node with no label" << endl;
						return false;}


					//Set the colors of vertices
					colorCodeArr->InsertValue(i,channelColors[label[0]]);
					i++;
					pts->InsertPoint(vertexID,iX,iY,iZ);
				} else if (strcmp(nodeName,"edge")== 0) {
					if (e == 0) g->SetPoints(pts);
					if(levelTwoElement->QueryIntAttribute("source", (int *)&isource) != TIXML_SUCCESS) 
					{
						cerr << "ERROR: encountered an edge with no source" << endl;
						return false;
					}
					if(levelTwoElement->QueryIntAttribute("target", (int *)&itarget) != TIXML_SUCCESS)
					{
						cerr << "ERROR: encountered an edge with no target" << endl;
						return false;
					}
					if(levelTwoElement->QueryDoubleAttribute("weight", (double *)&iweight) != TIXML_SUCCESS)
					{
						cerr << "ERROR: encountered an edge with no weight" << endl;
						return false;
					}


					//csource =(char *)xmlGetProp(nodeLevel2,source);
					//isource = atoi(csource);
					//ctarget =(char *)xmlGetProp(nodeLevel2,target);
					//itarget = atoi(ctarget);
					//cweight =(char *)xmlGetProp(nodeLevel2,weight);
					//iweight = atof(cweight);

					// Use the following if edge lengths are required for histograms
					// cout<<iweight<<endl;

					if (n == -1) {
						//Add this edge to the graph	
						edgeID=g->AddEdge(vertexNodeLabelToID[isource],vertexNodeLabelToID[itarget]);
						//Edge ids have to be string since we later use vtkSurfaceRepresentation for spheres.
						// During event handling this causes problems. We loose the types of selections.
						// But we need a mechanism for separating vertices from edges. If pedigreeIds have an 'e'
						// as the first letter, it represents and edge otherwise it has to be a vertex.
						stringstream ss;
						string s;
						ss << e;
						s = ss.str();
					
						edgeIDs->InsertNextValue(s);
						edgeIDs2->InsertNextValue(s);
						edgeWeightsArr->InsertNextValue(iweight);
						key= pair<int,int>(vertexNodeLabelToID[isource], vertexNodeLabelToID[itarget]);
						edgeSourgeTargetToID[key]=edgeID;
						edgeColors->InsertValue(e,0);
						avg = avg + iweight;  // average of edge weights is needed for pyramidial region
						e++;
					} else if (iweight < n) { // n is distance which is the first parameter to the program
						
						edgeID=g->AddEdge(vertexNodeLabelToID[isource],vertexNodeLabelToID[itarget]);
						//Edge ids have to be string since we later use vtkSurfaceRepresentation for spheres.
						// During event handling this causes problems. We loose the types of selections.
						// But we need a mechanism for separating vertices from edges. If pedigreeIds have an 'e'
						// as the first letter, it represents and edge otherwise it has to be a vertex.
						stringstream ss;
						string s;
						ss << e;
						s = ss.str();
					
						edgeIDs->InsertNextValue(s);
						edgeIDs2->InsertNextValue(s);
						edgeWeightsArr->InsertNextValue(iweight);
						key= pair<int,int>(vertexNodeLabelToID[isource], vertexNodeLabelToID[itarget]);
						edgeSourgeTargetToID[key]=edgeID;
						edgeColors->InsertValue(e,0);
						avg = avg + iweight;  // average of edge weights is needed for pyramidial region
						e++;}





				} else if (strcmp(nodeName,"text")== 0) 
							{ //Do Nothing 
							}
				   else {
						cout<<"Cannot recognize Node Element!"<<nodeName<<endl;
						return false;
					} //end:inner-if
			   levelTwoElement = levelTwoElement->NextSiblingElement();
			}//while
		} else {cout<<"Incorrect Graph File Format!"<<nodeName<<endl;
				return false;
		}

		levelOneElement=levelOneElement->NextSiblingElement();
	} //while

	//If we are here then the network has some vertices/edges
	setNetworkStatusOn();
	//if (e > 0) avg = avg/e;  //Average edge length for pyramidial region computation
	//cout<<"iste averaj "<<avg<<endl;
	//xmlSaveFile("xmlfile_copy.xml", doc);
	//xmlFreeDoc(doc);
	//Ad vertex ids
	g->GetVertexData()->AddArray(labels);
	//Add object colors
	g->GetVertexData()->AddArray(colorCodeArr);
	//Ad edge weights
	g->GetEdgeData()->AddArray(edgeWeightsArr);
	//Add Edge colors
	g->GetEdgeData()->AddArray(edgeColors);

	// Edge Ids are required for displaying and updating networks
	g->GetEdgeData()->AddArray(edgeIDs);
	g->GetVertexData()->AddArray(vertexIDs);

	g->GetVertexData()->SetPedigreeIds(vertexIDs2);
	g->GetEdgeData()->SetPedigreeIds(edgeIDs2);
	//edgeIDs->SetName("EdgeIDs");
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
void BioNet::Display(vtkRenderView* ren, vtkAlgorithm* alg, 
  double xoffset, double yoffset, 
  const char* vertColorArray, double vertMin, double vertMax, 
  const char* edgeColorArray, double edgeMin = 0, double edgeMax = 0)
{
  VTK_CREATE(vtkGraphToPolyData, graphToPoly);
  graphToPoly->SetInputConnection(alg->GetOutputPort());

  VTK_CREATE(vtkGlyphSource2D, glyph);
  glyph->SetGlyphTypeToVertex();
  

 /*
//  If we want vertices as only 2D circles
	glyph->SetGlyphTypeToCircle();
	glyph->FilledOn();
	glyph->CrossOff();
	glyph->SetScale(5.0);
*/

  VTK_CREATE(vtkGlyph3D, vertexGlyph);
  vertexGlyph->SetInputConnection(0, graphToPoly->GetOutputPort());
  vertexGlyph->SetInputConnection(1, glyph->GetOutputPort());
 
  

 
  // Vertices as 3D spheres 
  VTK_CREATE(vtkSphereSource, sphere);
  sphere->SetThetaResolution(8); sphere->SetPhiResolution(8);
  sphere->SetRadius(4);
  vertexGlyph->SetInputConnection(1, sphere->GetOutputPort());
  //


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

  //ren->GetRenderer()->AddActor(vertexActor);
  //ren->GetRenderer()->AddActor(edgeActor);

  //vtkGeometryFilter* vv=vtkGeometryFilter::New();


  VTK_CREATE(vtkSurfaceRepresentation, vertexRep);
  vertexGlyph->Update();
  //vertexGlyph->GetOutput()->Print(cerr);
  vertexRep->SetInputConnection(vertexGlyph->GetOutputPort());
  //vertexRep->SetColorModeToMapScalars();
  vertexRep->SetCellColorArrayName(vertColorArray);
  vertexRep->SetCellColorLUT(lut);
  vertexRep->SetScalarRange(vertMin, vertMax);
  ren->AddRepresentation(vertexRep);

  // We want to pass both the selections and the instance of BioNet class.
  // But we have only one clientdata to be set for this.
	BioNetAndSelectionLink* bnsl=new BioNetAndSelectionLink();
	VTK_CREATE(vtkCallbackCommand,vertexPicked);
	//Jeff Baumes:vertexPicked->SetCallback(this->ProcessVertexSelection);
	
	vertexPicked->SetCallback(this->ProcessEvents);
	bnsl->bioObj=this;
	bnsl->edgeLink=vertexRep->GetSelectionLink();
	vertexPicked->SetClientData((void *)static_cast<BioNetAndSelectionLink*>(bnsl));
	vertexRep->AddObserver(vtkCommand::SelectionChangedEvent, vertexPicked);
	//vertexRep->AddObserver(vtkCommand::KeyPressEvent, vertexPicked);
	//vertexRep->AddObserver(vtkCommand::CharEvent, vertexPicked);
}

// vtk does not provide edge remove.
// This does not delete out edges yet!
/*
void BioNet::RemoveEdge(vtkGraph* original_graph, vtkIdType u,vtkIdType v)
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
*/

// Description:
// vtk does not have RemoveEdge.
template<class EdgeType>
void BioNet::RemoveEdge(EdgeType id_of_edge_to_remove)
{
	//Store the current edge pedigree ids
	//vtkAbstractArray* VpedigreeIDs= g->GetVertexData()->GetPedigreeIds();
	vtkAbstractArray* EpedigreeIDs= g->GetEdgeData()->GetPedigreeIds();

	vtkSelectionSource* ss = vtkSelectionSource::New();
	ss->AddID(-1, id_of_edge_to_remove.Id);
	ss->SetFieldType(vtkSelectionNode::EDGE);
	ss->SetContentType(vtkSelectionNode::PEDIGREEIDS);
	ss->SetInverse(1);

	vtkExtractSelectedGraph* e = vtkExtractSelectedGraph::New();
	e->SetInput(0, g); //Pass the original graph first
	e->SetInputConnection(1, ss->GetOutputPort());
	e->Update();
	
	//Update the original graph
	g=vtkMutableUndirectedGraph::SafeDownCast(e->GetOutput()); 
	//Take care of the edge pedigree Ids. Restore them
	g->GetEdgeData()->SetPedigreeIds(EpedigreeIDs);
	//this->UpdateView();

};

// Description:
// vtk does not have RemoveVertex
void BioNet::RemoveVertex(vtkIdType id_of_vertex_to_remove)
{
	//Store the current vertex pedigree ids
	vtkAbstractArray* VpedigreeIDs= g->GetVertexData()->GetPedigreeIds();
	//vtkAbstractArray* EpedigreeIDs= g->GetEdgeData()->GetPedigreeIds();

	//vtkDataArray* Vp= g->GetVertexData()->GetPedigreeIds();
	//Vp->RemoveTuple(id_of_vertex_to_remove);

	vtkSelectionSource* ss = vtkSelectionSource::New();
	ss->AddID(-1, id_of_vertex_to_remove);
	ss->SetFieldType(vtkSelectionNode::VERTEX);
	ss->SetContentType(vtkSelectionNode::PEDIGREEIDS);
	ss->SetInverse(1);

	vtkExtractSelectedGraph* e = vtkExtractSelectedGraph::New();
	e->SetInput(0, g);
	e->SetInputConnection(1, ss->GetOutputPort());
	e->Update();

	//Update the original graph
	g = vtkMutableUndirectedGraph::SafeDownCast(e->GetOutput()); 
	//Take care of the vertex pedigree Ids. Restore them
	g->GetVertexData()->SetPedigreeIds(VpedigreeIDs);

	//cout<<"Testing ... !"<<endl;
	//cout<<"pedigree array size: "<<VpedigreeIDs->GetNumberOfTuples()<<endl;
	//cout<<"Num og vertices in the graph: "<<g->GetNumberOfVertices()<<endl;
};

// Description:
// Records the displays in order to get an animation
void BioNet::RecordImage (vtkRenderer* ren, vtkRenderWindow* win) {

	ren->GetActiveCamera()->Zoom(1);
	for(int counter = 0; counter< 9; counter++)
	{
		//ren->GetActiveCamera()->Azimuth(10);		
		ren->GetActiveCamera()->Roll(-10);
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
		
	ren->GetActiveCamera()->Elevation(-89);
	for(int counter = -90; counter< 90; counter++)
	{
		printf("%d\n",counter);
		ren->GetActiveCamera()->Elevation(1);
		//ren->GetActiveCamera()->Roll(-10);
		ren->ResetCameraClippingRange();
		ren->Render();
		vtkWindowToImageFilter* img = vtkWindowToImageFilter::New();
		img->SetInput(win);
		vtkTIFFWriter *writer = vtkTIFFWriter::New();
		writer->SetInput(img->GetOutput());
		char buff[10024];
		sprintf(buff,"C:\\Render\\render_%d.tif",counter+90);
		writer->SetFileName(buff);
		writer->Update();
	}

/*	for(int counter = 0; counter< 90; counter++)
	{
		printf("%d\n",counter);
		ren->GetActiveCamera()->Elevation(1);
		//ren->GetActiveCamera()->Roll(-10);
		ren->ResetCameraClippingRange();
		ren->Render();
		vtkWindowToImageFilter* img = vtkWindowToImageFilter::New();
		img->SetInput(win);
		vtkTIFFWriter *writer = vtkTIFFWriter::New();
		writer->SetInput(img->GetOutput());
		char buff[10024];
		sprintf(buff,"C:\\Render\\render_%d.tif",counter);
		writer->SetFileName(buff);
		writer->Update();
	}*/

	/*	for(int counter = 0; counter< 9; counter++)
	{
		ren->GetActiveCamera()->Elevation(10);
		//ren->GetActiveCamera()->Roll(-10);
		ren->ResetCameraClippingRange();
		ren->Render();
		vtkWindowToImageFilter* img = vtkWindowToImageFilter::New();
		img->SetInput(win);
		vtkTIFFWriter *writer = vtkTIFFWriter::New();
		writer->SetInput(img->GetOutput());
		char buff[10024];
		sprintf(buff,"C:\\Render\\render_%d.tif",counter);
		writer->SetFileName(buff);
		writer->Update();
	}*/

}

// Description:
//
template<class T> 
void ClearSelections(set<T>* s1) {
	s1->erase(s1->begin(),s1->end());
};

// Description:
//
template<class T2> 
void CopySets(set<T2> s1, set<T2>* s2) {
	set<T2>::iterator i1;
	for(i1 = s1.begin(); i1 != s1.end(); i1++) {
		s2->insert(*i1);
	}

};

void BioNet::UpdateView(){
	this->degree->SetInput(this->g);
	this->view->RemoveAllRepresentations();
	this->view->GetRenderer()->RemoveAllViewProps();

	this->Display(this->view, this->degree,0, 0, "ChannelColors", 0, 5,"EdgeColors");	

	this->AddLabels();
	BioNetAndSelectionLink* bnsl=new BioNetAndSelectionLink();
	bnsl->bioObj=this;

	edgePicked=vtkCallbackCommand::New();
	selections=vtkSurfaceRepresentation::New();
	edgePicked->SetCallback(this->ProcessEvents);

	edgePicked->SetClientData((void *)static_cast<BioNetAndSelectionLink*>(bnsl));
	selections->AddObserver(vtkCommand::SelectionChangedEvent, edgePicked);

	//selections=vtkSurfaceRepresentation::New();	
	view->SetSelectionType(vtkSelectionNode::PEDIGREEIDS);
	selections->SetInputConnection(this->ToPoly()->GetOutputPort());
	bnsl->edgeLink = vtkSelectionLink::New();
	selections->SetSelectionLink(bnsl->edgeLink);
	view->AddRepresentation(selections);
	view->Update();
}

void BioNet::ProcessEvents(vtkObject* caller, 
				   unsigned long event,
				   void* clientdata, 
				   void* callerdata) 
{
	vtkIdType id;
	set <double>::iterator i;
	set <string>::iterator j;
	//Get the BioNet object
	BioNetAndSelectionLink* bnsl=(BioNetAndSelectionLink *)clientdata;
	//Switch-BEGIN
    switch( event )
    {
    case vtkCommand::SelectionChangedEvent:
		if ((bnsl->bioObj->sve->vertices.size() == 0) && (bnsl->bioObj->sve->edges.size() == 0)) {
			// Delete old selections from BioNet object first
			ClearSelections(&bnsl->bioObj->finalSelections->vertices);
			ClearSelections(&bnsl->bioObj->finalSelections->edges);
			bnsl->bioObj->sve=bnsl->bioObj->GetSelections(bnsl->edgeLink);
			//  Store what we read initially
			CopySets(bnsl->bioObj->sve->vertices, &bnsl->bioObj->finalSelections->vertices);
			CopySets(bnsl->bioObj->sve->edges, &bnsl->bioObj->finalSelections->edges);
			bnsl->bioObj->numOfPasses++; //We read either the vertices or edges now
//			cout<<"numOfPasses = "<<bnsl->bioObj->numOfPasses<<endl;
		};

		if (bnsl->bioObj->sve->vertices.size() == 0) {
				set<string> tempStrSet;
				//(bnsl->bioObj->sve->edges);
				CopySets(bnsl->bioObj->sve->edges, &tempStrSet);
				//set<string> tempStrMap(bnsl->bioObj->sve->edges);
				//Store the old vertex selections
				//cout<<"11111"<<endl;
				copy(bnsl->bioObj->sve->edges.begin(),bnsl->bioObj->sve->edges.end(), tempStrSet.begin());
				//Get the new selections - new vertices
				bnsl->bioObj->GetSelections(bnsl->edgeLink);
				//copy(bnsl->bioObj->sve->edges.begin(),bnsl->bioObj->sve->edges.end(), ostream_iterator<string>(cout, " "));
				//cout<<endl;
				//copy(bnsl->bioObj->sve->vertices.begin(),bnsl->bioObj->sve->vertices.end(), ostream_iterator<double>(cout, " "));
				cout<<endl;
				copy(tempStrSet.begin(),tempStrSet.end(), bnsl->bioObj->sve->edges.begin());
				bnsl->bioObj->numOfPasses++; //The second pass is now completed
				//cout<<"numOfPasses1 = "<<bnsl->bioObj->numOfPasses<<endl;
		}

		//cout<<"bbbbbbbbbbbbb"<<endl;

		if (bnsl->bioObj->sve->edges.size() == 0) {
				set<double> tempDblSet(bnsl->bioObj->sve->vertices);

				copy(bnsl->bioObj->sve->vertices.begin(),bnsl->bioObj->sve->vertices.end(), tempDblSet.begin());	
			
//				copy(bnsl->bioObj->sve->vertices.begin(),bnsl->bioObj->sve->vertices.end(), ostream_iterator<double>(cout, " "));

				//Store the old vertex selections
				//copy(tempDblMap.begin(),tempDblMap.end(), ostream_iterator<double>(cout, " "));

				//Get the new selections - new edges
				bnsl->bioObj->sve=bnsl->bioObj->GetSelections(bnsl->edgeLink);//Fills out new edges
//				copy(bnsl->bioObj->sve->edges.begin(),bnsl->bioObj->sve->edges.end(), ostream_iterator<string>(cout, " "));
//				cout<<endl;

				//Update only the vertices with old ones
				CopySets(tempDblSet,&bnsl->bioObj->sve->vertices);
				bnsl->bioObj->numOfPasses++;
//				cout<<"numOfPasses2 = "<<bnsl->bioObj->numOfPasses<<endl;				
			}


	if (bnsl->bioObj->numOfPasses >= 3) //We have both vertices and edges now
		{
			// Take care of a VTK Bug first
			// If only an edge is picked, an additional vertex is included with the selection
			// Remove this vertex! Or don't copy it at all!
			//if ((bnsl->bioObj->sve->edges.size() == 1) && (bnsl->bioObj->sve->vertices.size() == 1)) {
				//Store only edges
			//	CopySets(bnsl->bioObj->sve->edges, &bnsl->bioObj->finalSelections->edges);
			//}
			//else {
				// Store the selections 
				CopySets(bnsl->bioObj->sve->vertices, &bnsl->bioObj->finalSelections->vertices);
				CopySets(bnsl->bioObj->sve->edges, &bnsl->bioObj->finalSelections->edges);
			//}			
			//Delete the selected vertices and edges 
			ClearSelections(&bnsl->bioObj->sve->vertices);
			ClearSelections(&bnsl->bioObj->sve->edges);
			bnsl->bioObj->numOfPasses=0; //Prepare for the next round	
			cout<<"Number of selected vertices: "<<bnsl->bioObj->finalSelections->vertices.size()<<endl;
			cout<<"Number of selected edges:: "<<bnsl->bioObj->finalSelections->edges.size()<<endl;	

			cout<<"Selected vertices: "<<endl;
			for(i = bnsl->bioObj->finalSelections->vertices.begin(); 
			i != bnsl->bioObj->finalSelections->vertices.end(); i++)
				cout<<"Id : "<<*i<<endl;				
		
			cout<<"Selected edges: "<<endl;
			for(j = bnsl->bioObj->finalSelections->edges.begin(); 
			j != bnsl->bioObj->finalSelections->edges.end(); j++)
				cout<<"Id : "<<*j<<endl;

			//bnsl->bioObj->g->GetVertexData()->GetArray("ChannelColors")->RemoveTuple(435);

			//bnsl->bioObj->UpdateView();

			/*
			//VTK_CREATE(vtkVertexDegree, degree); //Required for displaying graphs
	
			bnsl->bioObj->degree->SetInput(bnsl->bioObj->g);
			//degree->SetInput(e->GetOutput());
	
			bnsl->bioObj->view->RemoveAllRepresentations();

			bnsl->bioObj->Display(bnsl->bioObj->view, bnsl->bioObj->degree,0, 0, "ChannelColors", 0, 5,"EdgeColors");	
			//selections->AddObserver(vtkCommand::SelectionChangedEvent, edgePicked);
			bnsl->bioObj->view->AddRepresentation(bnsl->bioObj->selections);
			//view->Update();
			//bnsl->bioObj->view->GetRenderer()->ResetCamera();
			bnsl->bioObj->view->Update();
*/


		}
        break;

    case vtkCommand::KeyPressEvent:

		// Rotate by x-axis
		if(string(bnsl->bioObj->iren->GetKeySym())=="x"){		
			bnsl->bioObj->view->GetRenderer()->GetActiveCamera()->Elevation(10);
			bnsl->bioObj->view->GetRenderer()->ResetCameraClippingRange();
			bnsl->bioObj->view->GetRenderer()->Render();
			bnsl->bioObj->UpdateView();
		};

		// Rotate by x-axis - opposite direction
		if(string(bnsl->bioObj->iren->GetKeySym())=="X"){		
			bnsl->bioObj->view->GetRenderer()->GetActiveCamera()->Elevation(-10);
			bnsl->bioObj->view->GetRenderer()->ResetCameraClippingRange();
			bnsl->bioObj->view->GetRenderer()->Render();
			bnsl->bioObj->UpdateView();
		};

		// Rotate by y-axis
		if(string(bnsl->bioObj->iren->GetKeySym())=="y"){
			bnsl->bioObj->view->GetRenderer()->GetActiveCamera()->Azimuth(10);
			bnsl->bioObj->view->GetRenderer()->ResetCameraClippingRange();
			bnsl->bioObj->view->GetRenderer()->Render();
			bnsl->bioObj->UpdateView();
		};

		// Rotate by y-axis - opposite direction
		if(string(bnsl->bioObj->iren->GetKeySym())=="Y"){
			bnsl->bioObj->view->GetRenderer()->GetActiveCamera()->Azimuth(-10);
			bnsl->bioObj->view->GetRenderer()->ResetCameraClippingRange();
			bnsl->bioObj->view->GetRenderer()->Render();
			bnsl->bioObj->UpdateView();
		};

		// Rotate by z-axis
		if(string(bnsl->bioObj->iren->GetKeySym())=="z"){
			bnsl->bioObj->view->GetRenderer()->GetActiveCamera()->Roll(10);
			bnsl->bioObj->view->GetRenderer()->ResetCameraClippingRange();
			bnsl->bioObj->view->GetRenderer()->Render();
			bnsl->bioObj->UpdateView();
		};

		// Rotate by z-axis - opposite direction
		if(string(bnsl->bioObj->iren->GetKeySym())=="Z"){
			bnsl->bioObj->view->GetRenderer()->GetActiveCamera()->Roll(-10);
			bnsl->bioObj->view->GetRenderer()->ResetCameraClippingRange();
			bnsl->bioObj->view->GetRenderer()->Render();
			bnsl->bioObj->UpdateView();
		};

		
		if(string(bnsl->bioObj->iren->GetKeySym())=="Delete"){
			// EDGE DELETION
			/*
			string spid;
			vtkIdType id1;
			set <string>::reverse_iterator s;
			for(s = bnsl->bioObj->finalSelections->edges.rbegin(); s != bnsl->bioObj->finalSelections->edges.rend(); s++) {
				cout<<"Id : "<<*s<<endl;
		
				//Edge ids are type of string. Convert string to int first
				istringstream buffer(*s);
				int eid;
				buffer >> eid;
				cout<<"iste eid: "<<eid<<endl;
				id1 = (vtkIdType) eid;
				spid = bnsl->bioObj->GetNetwork()->GetEdgeData()->GetPedigreeIds()->GetVariantValue(id1).ToString();
				cout<<"Iste spid : "<<spid<<endl;


				istringstream buffer2(spid);
				buffer>>eid;

				cout<<"iste eid2: "<<eid<<endl;

				vtkEdgeType eeid;
				eeid.Id=(vtkIdType) eid;

				bnsl->bioObj->RemoveEdge(eeid);
				cout<<"AS-edge"<<endl;
				}
				*/
			//			cout<<"EDGE DELETION ENDS "<<endl;

			//reverse(bnsl->bioObj->finalSelections->vertices.begin(),bnsl->bioObj->finalSelections->vertices.end());

			// VERTEX DELETION
			double pid;
			//Delete vertices
			set <double>::reverse_iterator r;
			for(r = bnsl->bioObj->finalSelections->vertices.rbegin(); r != bnsl->bioObj->finalSelections->vertices.rend(); r++) {
				id = (vtkIdType) *r;
				pid = bnsl->bioObj->GetNetwork()->GetVertexData()->GetPedigreeIds()->GetVariantValue(id).ToDouble();
				bnsl->bioObj->RemoveVertex(pid);
			}

			//Now Update the view
			bnsl->bioObj->UpdateView();
			bnsl->bioObj->numOfPasses=0;
			ClearSelections(&bnsl->bioObj->sve->vertices);
			ClearSelections(&bnsl->bioObj->sve->edges);
		}
			
			
		/*
				VTK_CREATE(vtkVertexDegree, degree); //Required for displaying graphs
				degree->SetInput(bnsl->bioObj->GetNetwork());				
				bnsl->bioObj->view->RemoveAllRepresentations();
				bnsl->bioObj->Display(bnsl->bioObj->GetRenderView(), degree,0, 0, "ChannelColors", 0, 5,"EdgeColors");
				bnsl->bioObj->view->Update();
		*/


		/*
			//Delete vertices
			for(i = bnsl->bioObj->finalSelections->vertices.begin(); i != bnsl->bioObj->finalSelections->vertices.end(); i++) {
				cout<<"Id : "<<*i<<endl;
				id = (vtkIdType) *i;
				bnsl->bioObj->RemoveVertex(id);
				cout<<"AS"<<endl;
				//bnsl->bioObj->RemoveVertex(id);
			}
			
			// Deletion has been completed. Prepare for the next selections
			bnsl->bioObj->numOfPasses=0;
			bnsl->bioObj->visited=false;
			ClearSelections(&bnsl->bioObj->sve->vertices);
			ClearSelections(&bnsl->bioObj->sve->edges);
			ClearSelections(&bnsl->bioObj->finalSelections->vertices);
			ClearSelections(&bnsl->bioObj->finalSelections->edges);

			}
*/

        break;
    default:
		break;
    } //Switch-END

}

			/*
				cout<<"Id : "<<435<<endl;
				id = (vtkIdType) 435;
				bnsl->bioObj->RemoveVertex(id);
				cout<<"Id : "<<435<<"silindi"<<endl;

				cout<<"Id : "<<426<<endl;
				id = (vtkIdType) 426;
				bnsl->bioObj->RemoveVertex(id);
				cout<<"Id : "<<426<<"silindi"<<endl;

				cout<<"Id : "<<394<<endl;
				id = (vtkIdType) 394;
				bnsl->bioObj->RemoveVertex(id);
				cout<<"Id : "<<394<<"silindi"<<endl;
				cout<<"AS"<<endl;
			*/

// Description:
// A helper-function that is used for detecting cycles. If the vertex v is already in  /
// vertices, it means there is a back-edge from it causing a cycle.
void BioNet::InsertIntoMap(set<double>* vertices, double d) {
	set<double>::iterator pos;
	pos=vertices->find(d);
	if (pos==vertices->end()) {
		vertices->insert(d);
	};
}	

		/*
		if ((bnsl->bioObj->numOfPasses == 2) && (bnsl->bioObj->visited==0)) {
			cout<<"Geldim"<<endl;
					// Store the selections 
					CopySets(bnsl->bioObj->sve->vertices, &bnsl->bioObj->finalSelections->vertices);
					CopySets(bnsl->bioObj->sve->edges, &bnsl->bioObj->finalSelections->edges);
			
					//Empty the temporary locations
					//Delete the selected vertices and edges 
					ClearSelections(&bnsl->bioObj->sve->vertices);
					ClearSelections(&bnsl->bioObj->sve->edges);
					bnsl->bioObj->numOfPasses=0; //Prepare for the next round
		}
*/

		/*
		if ((bnsl->bioObj->numOfPasses >= 2) && (bnsl->bioObj->visited==true)) //We have both vertices and edges
		{
			cout<<"ALLLLLLL"<<endl;
			// Take care of a VTK Bug first
			// If only an edge is picked, an additional vertex is included with the selection
			// Remove this vertex! Or don't copy it at all!
			if ((bnsl->bioObj->sve->edges.size() == 1) && (bnsl->bioObj->sve->vertices.size() == 1)) {
				//Store only edges
				CopySets(bnsl->bioObj->sve->edges, &bnsl->bioObj->finalSelections->edges);
				cout<<"ALLLLLLL222222"<<endl;
			}
			else {
				// Store the selections 
				CopySets(bnsl->bioObj->sve->vertices, &bnsl->bioObj->finalSelections->vertices);
				CopySets(bnsl->bioObj->sve->edges, &bnsl->bioObj->finalSelections->edges);
			}
				//Empty the temporary locations
				//Delete the selected vertices and edges 
			ClearSelections(&bnsl->bioObj->sve->vertices);
			ClearSelections(&bnsl->bioObj->sve->edges);
			bnsl->bioObj->numOfPasses=0; //Prepare for the next round	
	cout<<"Final edge size soyle: "<<bnsl->bioObj->finalSelections->edges.size()<<" tamam"<<endl;
		cout<<"Final vertex size soyle: "<<bnsl->bioObj->finalSelections->vertices.size()<<" tamam"<<endl;	
		}
*/

/************************
		bnsl=(BioNetAndSelectionLink *)clientdata;
		bnsl->bioObj->sve=bnsl->bioObj->GetSelections(bnsl->edgeLink);
		cout<<"key "<<bnsl->bioObj->iren->GetKeySym()<<endl;
		if(string(bnsl->bioObj->iren->GetKeySym())=="Delete"){
			cout<<"Num of vertices for sve3 is: "<<bnsl->bioObj->sve->vertices.size()<<endl;
			cout<<"Num of edges for sve3 is: "<<bnsl->bioObj->sve->vertices.size()<<endl;

			for(i = bnsl->bioObj->sve->vertices.begin(); i != bnsl->bioObj->sve->vertices.end(); i++) {
				cout<<"Id : "<<*i<<endl;
				id = (vtkIdType) *i;
				bnsl->bioObj->RemoveVertex(id);
			}
		} else if(string(bnsl->bioObj->iren->GetKeySym())=="a"){
				cout<<"Vertex Addition Starts "<<endl;
		}

		//bnsl->bioObj->view->Update();
		cout<<"LBR"<<endl;
************************/


// Description:
//
SelectedVerticesAndEdges* BioNet::GetSelections(vtkSelectionLink* sel) {

	SelectedVerticesAndEdges* sve = new SelectedVerticesAndEdges();
	vtkSelection *p = sel->GetSelection();
	
	for(int counter=0; counter < p->GetNumberOfNodes(); counter++)
	{
		vtkSelectionNode *n = p->GetNode(counter);
		vtkAbstractArray* arra = n->GetSelectionList();
		 for (int i=0;i<n->GetSelectionList()->GetNumberOfTuples();i++) {
			 vtkVariant aaa=arra->GetVariantValue(i);
			// Decide if this is a vertex or an edge
			 if (aaa.GetTypeAsString() == "double" ) //VERTEX ids have type of double 
				 //vtkIdType id = (vtkIdType) aaa.ToInt();
				 //vtkIdType id = (vtkIdType) this->GetNetwork()->GetVertexData()->GetPedigreeIds()->GetVariantValue(aaa).ToDouble()
				 //this->InsertIntoMap(&sve->vertices, this->GetNetwork()->GetVertexData()->GetPedigreeIds()->GetVariantValue(id).ToDouble());				
		 
				 this->InsertIntoMap(&sve->vertices, (aaa.ToDouble()));
			 else if (aaa.GetTypeAsString() == "string" ) //edge ids have type of string
				  sve->edges.insert(aaa.ToString());
		 }
		 cout<<endl<<endl;
	}

	return sve;

}

// Description:
// Constructor for BioNet
BioNet::BioNet() {
	g=vtkMutableUndirectedGraph::New();
	emptyNetwork=true;	/*The network is initially has no vertices and edges*/
	view=vtkRenderView::New();
	//SelectedVerticesAndEdges*
	sve=new SelectedVerticesAndEdges();
	finalSelections = new SelectedVerticesAndEdges(); // Initialize selections - nothing selected yet
	numOfPasses=0; //Event handler does not return anything yet!
	visited=false;
	degree=vtkVertexDegree::New();

}  

// Description:
// Destructor for BioNet
BioNet::~BioNet() {

	view->Delete();
	selections->Delete();
	edgePicked->Delete();
	g->Delete();
	iren->Delete();
	win->Delete();	
	//labelActor->Delete();
	//labelMapper->Delete();
	//poly->Delete();
}


vtkRenderView* BioNet::GetRenderView(){
	return view;
};

// Description:
//
vtkMutableUndirectedGraph* BioNet::GetNetwork() {
	return g;
};

// Description:
// Converts a given network to GraphToPolyData form
vtkGraphToPolyData* BioNet::ToPoly() {
	vtkGraphToPolyData* poly = vtkGraphToPolyData::New();
	poly->SetInput(GetNetwork());
	return poly;
}

// Description:
// Attaches Labels to vertices
void BioNet::AddLabels() {
	
	vtkGraphToPolyData* poly = ToPoly();
	//Add Labels	
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
	this->GetRenderView()->GetRenderer()->AddActor(labelActor);
	//view->GetRenderer()->AddActor(labelActor );
};

// Description:
// Creates the interactor so picking can occur
void BioNet::Interact()
{	
	BioNetAndSelectionLink* bnsl=new BioNetAndSelectionLink();
	bnsl->bioObj=this;

	edgePicked=vtkCallbackCommand::New();
	selections=vtkSurfaceRepresentation::New();	
	view->SetupRenderWindow(win);
	view->GetRenderer()->ResetCamera();
	view->Update();

	view->SetSelectionType(vtkSelectionNode::PEDIGREEIDS);
	selections->SetInputConnection(this->ToPoly()->GetOutputPort());
	selections->SetSelectionLink(bnsl->edgeLink);
	view->AddRepresentation(selections);
	view->Update();
    view->GetRenderer()->ResetCamera();
	edgePicked->SetCallback(this->ProcessEvents);

	edgePicked->SetClientData((void *)static_cast<BioNetAndSelectionLink*>(bnsl));
	selections->AddObserver(vtkCommand::SelectionChangedEvent, edgePicked);
	
	win->GetInteractor()->GetInteractorStyle()->AddObserver(vtkCommand::KeyPressEvent,edgePicked);
}

// Description:
// create the standard render window and interactor
void BioNet::RenderWin(){
	iren=vtkRenderWindowInteractor::New();
	win=vtkRenderWindow::New();
	win->SetSize(900, 900);
	win->SetInteractor(iren);
	//win->GetInteractor()->GetInteractorStyle()->AddObserver(vtkCommand::KeyPressEvent,edgePicked);
	VTK_CREATE(vtkInteractorStyleTrackballCamera,vC);
	iren->SetInteractorStyle(vC);
}

// Description:
// Constructor for BioNetAndSelectionLink
BioNetAndSelectionLink::BioNetAndSelectionLink() {
	bioObj=new BioNet(); 
	//bioObj->sve->edges
	edgeLink=vtkSelectionLink::New();
};


void BioNet::Kruskalmst() {

	VTK_CREATE(vtkBoostKruskalMinimumSpanningTree,mst);
	mst=vtkBoostKruskalMinimumSpanningTree::New();
	mst->SetInput(this->g);
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
	mstGraph=vtkExtractSelectedGraph::New();
	//VTK_CREATE(vtkExtractSelectedGraph, mstGraph);
	mstGraph->SetInput(0,this->g);
	mstGraph->SetInput(1,mst->GetOutput());
};

int main(int argc, char *argv[])
{
	VTK_CREATE(vtkVertexDegree, degree); //Required for displaying graphs

	cout<<endl<<endl;
	if (argc <4) {
			cout<<"Incorrect Usage! Check Parameters" <<endl;
			return 0;
	};

	//Read the vertices and edges into graph g
	//argv1= The name of the input network
	// Python Arguments
	//argv2= distance 
	//argv3= an integer that shows either the graph is standard or a minimum spanning tree
	//argv4= cutoff for pyramidial
	float distance=atof(argv[2]);
	int stdMst=atoi(argv[3]);
	float cutoff=atof(argv[4]);
	
	//Create an empty network
	BioNet* network=new BioNet();
	network->cutoff=cutoff;

	//Read the vertices and edges into the network
	network->ReadXGMML(argv[1],distance);
	// This is required to be done for the Display function
	degree->SetInput(network->GetNetwork());
	//Add Labels to vertices
	network->AddLabels();

	if (stdMst == 2) {
		network->cutoff=cutoff;
		network->Averages();
		network->Display(network->GetRenderView(), degree,0, 0, "ChannelColors", 0, 5,"EdgeColors");
	} else if (stdMst == 1) {
			network->Kruskalmst();
			network->Display(network->GetRenderView(), network->mstGraph,0, 0, "ChannelColors", 0, 5,"EdgeColors");			
	}
	else { 
		//Display the constructed network
		if (!network->NetworkEmpty())
			network->Display(network->GetRenderView(), degree,0, 0, "ChannelColors", 0, 5,"EdgeColors");
		else {cout<<"Network has no nodes!"<<endl;
			return 0;
		};
	};

	//Start the renderer
	network->RenderWin();
	//Start interactions
	network->Interact();
	//network->RecordImage(network->view->GetRenderer(),network->win);
	//Show the view
	network->iren->Initialize();
	network->iren->Start();

	return 0;


}
