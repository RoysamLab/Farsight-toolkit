//=======================================================================
// Copyright 2001 Jeremy G. Siek, Andrew Lumsdaine, Lie-Quan Lee, 
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//=======================================================================
#include "FTKgraph.h"


//***********************************************************************************************************
//  WRITE/READ IMAGES
//***********************************************************************************************************
template <typename T>
typename T::Pointer readImage(const char *filename)
{
    printf("Reading %s ... ",filename);
    typedef typename itk::ImageFileReader<T> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    
    ReaderType::GlobalWarningDisplayOff();
    reader->SetFileName(filename);
    try
    {
        reader->Update();
    }
    catch(itk::ExceptionObject &err)
    {
        std::cout << "ExceptionObject caught!" <<std::endl;
        std::cout << err << std::endl;
        //return EXIT_FAILURE;
    }
    printf("Done.\n");
    return reader->GetOutput();
}

template <typename T>
int writeImage(typename T::Pointer im, const char* filename)
{
    printf("Writing %s ... ",filename);
    typedef typename itk::ImageFileWriter<T> WriterType;
    
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(filename);
    writer->SetInput(im);
    try
    {
        writer->Update();
    }
    catch(itk::ExceptionObject &err)
    {
        std::cout << "ExceptionObject caught!" <<std::endl;
        std::cout << err << std::endl;
        return EXIT_FAILURE;
    }
    printf("Done.\n");
    return EXIT_SUCCESS;
}




FTKgraph::FTKgraph()
{
}


vtkSmartPointer<vtkTable> FTKgraph::AdjacencyGraph_All(InputImageType::Pointer ip_image, OutputImageType::Pointer op_image, bool CytoImage)
{
	cyto_image = CytoImage;
	vtkSmartPointer<vtkTable> graphtable = vtkSmartPointer<vtkTable>::New();
	graphtable->Initialize();
	vtkSmartPointer<vtkStringArray> column1 = vtkSmartPointer<vtkStringArray>::New();
    column1->SetName( "Source" );
    graphtable->AddColumn(column1);
    vtkSmartPointer<vtkStringArray> column2 = vtkSmartPointer<vtkStringArray>::New();
    column2->SetName( "Target" );
	graphtable->AddColumn(column2);

	graphtable = constructGraphTable_All(ip_image, op_image, graphtable, cyto_image);

	return graphtable;

}

vtkSmartPointer<vtkTable> FTKgraph::constructGraphTable_All(InputImageType::Pointer ip_image, OutputImageType::Pointer op_image, vtkSmartPointer<vtkTable> graphtable, bool CytoImage)
{
    ftkgnt * g1 = new ftkgnt();
	g1->runLabFilter(ip_image, op_image, CytoImage);
	int check;
    // Property accessors

	std::vector< FeatureCalcType::LabelPixelType > labels = g1->labelFilter->GetLabels();
	for (int i=0; i<(int)labels.size(); ++i)
    {
		FeatureCalcType::LabelPixelType id = labels.at(i);
		if(id == 0) 
			continue;
		
		check = 1;
		for(int i=0; i<(int)graphtable->GetNumberOfRows(); ++i)
	    {
		    if((static_cast<vtkVariant>(id)==graphtable->GetValue(i,0)) || (static_cast<vtkVariant>(id)==graphtable->GetValue(i,1))) 
		    {
			    check = 0;
				break;
			}
	    }	
		if (check == 1)
		{
			g1->RAG = g1->BuildRAG(id);
			
			index_map = get(boost::vertex_index, g1->RAG);
			name_map = get(boost::vertex_name, g1->RAG); 

			graphtable = BuildGraphTable( g1->RAG, graphtable);
		}
    }

    return graphtable;
}



vtkSmartPointer<vtkTable> FTKgraph::AdjacencyGraph_ID(unsigned short id, InputImageType::Pointer ip_image, OutputImageType::Pointer op_image, bool CytoImage) 
{
	cyto_image = CytoImage;
	vtkSmartPointer<vtkTable> graphtable = vtkSmartPointer<vtkTable>::New();
	graphtable->Initialize();
	vtkSmartPointer<vtkStringArray> column1 = vtkSmartPointer<vtkStringArray>::New();
    column1->SetName( "Source" );
    graphtable->AddColumn(column1);
    vtkSmartPointer<vtkStringArray> column2 = vtkSmartPointer<vtkStringArray>::New();
    column2->SetName( "Target" );
	graphtable->AddColumn(column2);

	graphtable = constructGraphTable_ID(id, ip_image, op_image, graphtable, cyto_image);

	return graphtable;
}

vtkSmartPointer<vtkTable> FTKgraph::constructGraphTable_ID(unsigned short id, InputImageType::Pointer ip_image, OutputImageType::Pointer op_image, vtkSmartPointer<vtkTable> graphtable, bool CytoImage)
{
    ftkgnt * g1 = new ftkgnt();
	g1->runLabFilter(ip_image, op_image, CytoImage);
	
    // Property accessors
	
	g1->RAG = g1->BuildRAG(id);
    
    index_map = get(boost::vertex_index, g1->RAG);
    name_map = get(boost::vertex_name, g1->RAG); 
    
    graphtable = BuildGraphTable(g1->RAG, graphtable);
  
    return graphtable;
}



vtkSmartPointer<vtkTable> FTKgraph::BuildGraphTable( GraphType g, vtkSmartPointer<vtkTable> table)
{
   for (boost::tie(vi,vi_end) = vertices(g) ; vi != vi_end ; ++vi)
   {
	   for (boost::tie(ai,ai_end) = adjacent_vertices(*vi, g) ; ai != ai_end ; ++ai)
	   {
           flag = 1;
	       for(int i=0; i<(int)table->GetNumberOfRows(); ++i)
	       {
		       //std::cout<<i<<std::endl;
		       if(((static_cast<vtkVariant>(name_map[*vi])==table->GetValue(i,0))&&(static_cast<vtkVariant>(name_map[*ai])==table->GetValue(i,1))) || 
			       ((static_cast<vtkVariant>(name_map[*ai])==table->GetValue(i,0))&&(static_cast<vtkVariant>(name_map[*vi])==table->GetValue(i,1))))
		       {
			       flag = 0;
			       break;
		       }
	       }
	       if(flag == 1)
	       {
	           vtkSmartPointer<vtkVariantArray> model_data1 = vtkSmartPointer<vtkVariantArray>::New();
     	       model_data1->InsertNextValue(static_cast<vtkVariant>(name_map[*vi]));
		       model_data1->InsertNextValue(static_cast<vtkVariant>(name_map[*ai]));
		       table->InsertNextRow(model_data1);		  
	       }	   
       }
   }
   return table;
}



void FTKgraph::DisplayGraph(vtkSmartPointer<vtkTable> graphtable)
{
    
    if(!graphtable)
		return;
	vtkSmartPointer<vtkTableToGraph> TTG = vtkSmartPointer<vtkTableToGraph>::New();
    TTG->SetInputData(0, graphtable);
    TTG->AddLinkVertex("Source", "Vertex", false);
    TTG->AddLinkVertex("Target", "Vertex", false);
    TTG->AddLinkEdge("Source", "Target");
	
	/* VTK_CREATE(vtkFast2DLayoutStrategy, vertexStrategy);

	VTK_CREATE(vtkGraphLayout, vertexLayout);
	vertexLayout->SetLayoutStrategy(vertexStrategy);
	vertexLayout->SetInputConnection(TTG->GetOutputPort());

	VTK_CREATE(vtkRenderer, Ren1);

    VTK_CREATE(vtkGraphToGlyphs, vertexGlyph);
    vertexGlyph->SetInputConnection(vertexLayout->GetOutputPort());
	vertexGlyph->SetGlyphType(9);
	vertexGlyph->FilledOn();
	vertexGlyph->SetRenderer(Ren1);
	vertexGlyph->SetScreenSize(0.1);

	VTK_CREATE(vtkPolyDataMapper, vertexMapper);
    vertexMapper->SetInputConnection(vertexGlyph->GetOutputPort());
    vertexMapper->SetScalarModeToUsePointFieldData();

	VTK_CREATE(vtkActor, vertexActor);
    vertexActor->GetProperty()->SetPointSize(5.0);
    vertexActor->SetMapper(vertexMapper);
	vertexActor->GetProperty()->SetColor(1.0, 0.0, 0.0);
	vertexActor->GetProperty()->SetOpacity(1.0);

	VTK_CREATE(vtkArcParallelEdgeStrategy, edgeStrategy);

	VTK_CREATE(vtkEdgeLayout, edgeLayout);
	edgeLayout->SetLayoutStrategy(edgeStrategy);
	edgeLayout->SetInputConnection(vertexLayout->GetOutputPort());

	VTK_CREATE(vtkGraphToPolyData, graphToPoly);
	graphToPoly->SetInputConnection(edgeLayout->GetOutputPort());

	VTK_CREATE(vtkPolyDataMapper, edgeMapper);
    edgeMapper->SetInputConnection(graphToPoly->GetOutputPort());
	edgeMapper->SetScalarModeToUseCellFieldData();

	VTK_CREATE(vtkActor, edgeActor);
	edgeActor->SetMapper(edgeMapper);
	edgeActor->GetProperty()->SetColor(1.0, 1.0, 0.0);
	edgeActor->GetProperty()->SetOpacity(0.8);

	
	VTK_CREATE(vtkRenderer, Ren2);
	Ren2->AddActor(edgeActor);
	VTK_CREATE(vtkRenderer, Ren3);
	Ren3->AddActor(vertexActor);

	VTK_CREATE(vtkRenderWindow, RenWin);
	RenWin->AddRenderer(Ren2);
	RenWin->AddRenderer(Ren3);
	RenWin->SetSize(1000,1000);


	VTK_CREATE(vtkRenderWindowInteractor, iRen);
	iRen->SetRenderWindow(RenWin);

	RenWin->Render();
	iRen->Start();*/






   vtkSmartPointer<vtkViewTheme> theme;
   theme.TakeReference(vtkViewTheme::CreateMellowTheme());
   theme->SetLineWidth(5);
   theme->SetCellOpacity(0.9);
   theme->SetCellAlphaRange(0.5,0.5);
   theme->SetPointSize(10);
   theme->SetSelectedCellColor(1,0,1);
   theme->SetSelectedPointColor(1,0,1); 

   vtkSmartPointer<vtkGraphLayoutView> view = vtkSmartPointer<vtkGraphLayoutView>::New();
   view->AddRepresentationFromInputConnection(TTG->GetOutputPort());
	view->SetEdgeLabelVisibility(true);
	view->SetEdgeLabelArrayName("Distance");
   view->SetLayoutStrategyToSimple2D();
   view->SetVertexLabelArrayName("label");
   view->VertexLabelVisibilityOn();
   view->SetVertexLabelFontSize(20);
  
   view->GetRenderWindow()->SetSize(600, 600);
   view->ResetCamera();
   view->Render();
   view->GetInteractor()->Start();

   return ;
}


