#include "Image3DView.h"

Image3DView::Image3DView(ftk::Image::Pointer image,std::vector<std::map<int, ftk::Object::Point> > CenterMapVector, LabelImageViewQT * imview, ObjectSelection * sels)
{
	if(!image)return;

	this->LabelImageData = image;
	this->ImageView = imview;
	this->CenterMap = ImageView->GetCenterMapVector();// needs to be fixed later;
	this->bBoxMap = ImageView->GetBoxMapVector();
	this->Selection = sels;
	double before = clock( ) ;
	this->CreateActors();
	double after = clock( ) ;
	double time_elapsed = (after - before) ;
	std::cout << "Time elapsed is " << time_elapsed <<std::endl;	

	// Initialize Flags:
	this->labelsVisible = true;
	this->stacksVisible = true;

	// Setup the render window, interactor and render:
	this->SetupRenderWindow();
	this->CreateInteractorStyle();
	this->Render3DStack(0);

	// Connect signals and slot
	connect(ImageView, SIGNAL(emitTimeChanged()), this, SLOT(RefreshTimeActors()));
	connect(Selection, SIGNAL(changed()), this, SLOT(RefreshTrackActors()));

}
Image3DView::~Image3DView()
{
}

void Image3DView::CreateActors()
{ 	
	for (int i= 0; i<(int)CenterMap.size() ; i++)
	{
		this->CreateBoundaryActors(i);
		this->CreateCentroidActor(i);
		this->CreateLabelActor(i);
	}
}

void Image3DView::CreateCentroidActor(int currentT)
{
	std::map<int, ftk::Object::Point> centermap = CenterMap[currentT];
	std::map<int, ftk::Object::Point>::iterator cmap_iter;
	float col = (float)currentT/100.0;
	for(cmap_iter = centermap.begin();cmap_iter!= centermap.end(); cmap_iter++)
	{
		vtkSmartPointer<vtkSphereSource> centroidSphere =  vtkSmartPointer<vtkSphereSource>::New();
		centroidSphere->SetCenter((double)(*cmap_iter).second.x,
								  (double)(*cmap_iter).second.y,
								  (double)(*cmap_iter).second.z);
		centroidSphere->SetRadius(0.3);	
		// Setup Mapper and Actor:
		vtkSmartPointer<vtkPolyDataMapper> centroidMapper =  vtkSmartPointer<vtkPolyDataMapper>::New();
		centroidMapper->SetInputConnection(centroidSphere->GetOutputPort()); 
		vtkSmartPointer<vtkActor> centroidactor = vtkSmartPointer<vtkActor>::New();
		centroidactor->SetMapper(centroidMapper);
		centroidactor->GetProperty()->SetOpacity(0.7);
		centroidactor->GetProperty()->SetColor(1,0,0);
		CellActorsMap.at(currentT)[(*cmap_iter).first].centroidActor = centroidactor;
	}

}
void Image3DView::CreateBoundaryActors(int currentT)
{
	if(!LabelImageData) return;

	int ch = 0;
	int v, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15, v16, v17, v18, v19, v20, v21, v22, v23, v24, v25, v26;
	std::map<int, ftk::Object::Box> bbox =  bBoxMap[currentT];
	std::map<int, ftk::Object::Point> centermap = CenterMap[currentT];
	std::map<int, CellActors> actorsmap;
	std::map<int, ftk::Object::Box>::iterator bbox_iter;


	
	for(bbox_iter = bbox.begin();bbox_iter!=bbox.end(); bbox_iter++)
	{
		double xc = (double) centermap[(*bbox_iter).first].x;
		double yc = (double) centermap[(*bbox_iter).first].y;
		double zc = (double) centermap[(*bbox_iter).first].z;

		int xmin = (*bbox_iter).second.min.x;
		int xmax = (*bbox_iter).second.max.x;
		
		int ymin = (*bbox_iter).second.min.y;
		int ymax = (*bbox_iter).second.max.y;
		
		int zmin = (*bbox_iter).second.min.z;
		int zmax = (*bbox_iter).second.max.z;

		//create a vtkimage i.e a small image of the cell (segmented) bounded by it's bounding box:		
		vtkSmartPointer<vtkImageData> imdata = vtkSmartPointer<vtkImageData>::New();
		imdata->SetSpacing(1,1,1);
		imdata->SetOrigin(0,0,0);
		imdata->SetDimensions(xmax-xmin+1,ymax-ymin+1,zmax-zmin+1);
		imdata->SetScalarTypeToInt();
		imdata->SetNumberOfScalarComponents(1);
		imdata->AllocateScalars();
		int *imdataptr = (int *)imdata->GetScalarPointer();
		
		for(int k = zmin; k <= zmax; k++)
		{
			for(int i= ymin; i <= ymax; i++)
			{
				for(int j= xmin; j <= xmax; j++)
				{
					v = (int)LabelImageData->GetPixel(currentT, ch, k, i, j);
					*imdataptr = 0;
					if (v > 0)														// could be unnecessary
					{
						v1 = (int)LabelImageData->GetPixel(currentT, ch, k-1, i, j);
						v2 = (int)LabelImageData->GetPixel(currentT, ch, k+1, i, j);
						v3 = (int)LabelImageData->GetPixel(currentT, ch, k, i+1, j);
						v4 = (int)LabelImageData->GetPixel(currentT, ch, k, i-1, j);
						v5 = (int)LabelImageData->GetPixel(currentT, ch, k, i, j-1);
						v6 = (int)LabelImageData->GetPixel(currentT, ch, k, i, j+1);

						v7 = (int)LabelImageData->GetPixel(currentT, ch, k, i+1, j+1);
						v8 = (int)LabelImageData->GetPixel(currentT, ch, k, i+1, j-1);
						v9 = (int)LabelImageData->GetPixel(currentT, ch, k, i-1, j+1);
						v10 = (int)LabelImageData->GetPixel(currentT, ch, k, i-1, j-1);

						v11 = (int)LabelImageData->GetPixel(currentT, ch, k+1, i, j+1);
						v12 = (int)LabelImageData->GetPixel(currentT, ch, k+1, i, j-1);
						v13 = (int)LabelImageData->GetPixel(currentT, ch, k-1, i, j+1);
						v14 = (int)LabelImageData->GetPixel(currentT, ch, k-1, i, j-1);

						v15 = (int)LabelImageData->GetPixel(currentT, ch, k+1, i+1, j);
						v16 = (int)LabelImageData->GetPixel(currentT, ch, k+1, i-1, j);
						v17 = (int)LabelImageData->GetPixel(currentT, ch, k-1, i+1, j);
						v18 = (int)LabelImageData->GetPixel(currentT, ch, k-1, i-1, j);

						v19 = (int)LabelImageData->GetPixel(currentT, ch, k+1, i+1, j+1);
						v20 = (int)LabelImageData->GetPixel(currentT, ch, k+1, i+1, j-1);
						v21 = (int)LabelImageData->GetPixel(currentT, ch, k+1, i-1, j+1);
						v22 = (int)LabelImageData->GetPixel(currentT, ch, k+1, i-1, j-1);

						v23 = (int)LabelImageData->GetPixel(currentT, ch, k-1, i-1, j+1);
						v24 = (int)LabelImageData->GetPixel(currentT, ch, k-1, i-1, j-1);
						v25 = (int)LabelImageData->GetPixel(currentT, ch, k-1, i+1, j+1);
						v26 = (int)LabelImageData->GetPixel(currentT, ch, k-1, i+1, j-1);

						//if(v!=v1 || v!=v2 || v!=v3 || v!=v4 || v!=v5 || v!=v6)
						if(v!=v1 || v!=v2 || v!=v3 || v!=v4 || v!=v5 || v!=v6 
							|| v!=v7 || v!=v8 || v!=v9 || v!=v10 || v!=v11|| v!=v12 
							|| v!=v13 || v!=v14 || v!=v15 || v!=v16 || v!=v17 || v!=v18
							|| v!=v19 || v!=v20 || v!=v21 || v!=v22 || v!=v23 || v!=v24
							|| v!=v25 || v!=v26)
						{
							*imdataptr = 1;
						}
					}
					++imdataptr;
				}
			}
		}

		// Contour the cell:
		vtkSmartPointer<vtkMarchingCubes> contourf = vtkSmartPointer<vtkMarchingCubes>::New();
		contourf->ComputeNormalsOff();
		contourf->ComputeScalarsOff();
		contourf->ComputeGradientsOff();
		contourf->SetInput(imdata);
		contourf->SetNumberOfContours(1);
		contourf->SetValue(0,0.9);
		contourf->Update();
		// Color the contour:
		vtkSmartPointer<vtkFloatArray> cellcolor = vtkSmartPointer<vtkFloatArray>::New();
		cellcolor->SetNumberOfComponents(1);
		cellcolor->SetName("colors");
		float col = (float)(*bbox_iter).first/40.0;
		for( int i=0; i<(int) contourf->GetOutput()->GetNumberOfCells(); ++i) cellcolor->InsertNextTuple1(col);
		contourf->GetOutput()->GetCellData()->SetScalars(cellcolor);	
		// Smooth the contour:
		//vtkSmartPointer<vtkSmoothPolyDataFilter> smoothf = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
		//smoothf->SetInput(contourf->GetOutput());
		//smoothf->SetRelaxationFactor(0.1);
		//smoothf->SetNumberOfIterations(20);
		//smoothf->Update();
		// Create the Mapper and the Actor:
		vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInput(contourf->GetOutput());
		mapper->Update();
		vtkSmartPointer<vtkActor> contourActor = vtkSmartPointer<vtkActor>::New();
		contourActor->SetMapper(mapper);
		contourActor->GetProperty()->SetOpacity(0.8);
		contourActor->SetPosition((double)xmin,(double)ymin,(double)zmin);
		CellActors act;
		act.boundaryActor = contourActor;
		actorsmap.insert(std::pair<int,CellActors>((*bbox_iter).first,act));
	}
	CellActorsMap.push_back(actorsmap);
}
void Image3DView::CreateLabelActor(int currentT)
{
	std::map<int, ftk::Object::Point> centermap = CenterMap[currentT];
	std::map<int, ftk::Object::Point>::iterator cmap_iter;

	vtkSmartPointer<vtkPolyData> labelpolydata = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPoints> labelpoints = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkStringArray> labels = vtkSmartPointer<vtkStringArray>::New();	
	vtkSmartPointer<vtkIntArray> sizes = vtkSmartPointer<vtkIntArray>::New();

	labels->SetNumberOfValues(centermap.size());
	labels->SetName("labels");
	sizes->SetNumberOfValues(centermap.size());
	sizes->SetName("sizes");

	double labeloffset = 2.0;
	for(cmap_iter = centermap.begin();cmap_iter!= centermap.end(); cmap_iter++)
	{
		int return_id = labelpoints->InsertNextPoint((double)(*cmap_iter).second.x+labeloffset, 
													 (double)(*cmap_iter).second.y+labeloffset, 
													 (double)(*cmap_iter).second.z+labeloffset);
		//Labels:
		stringstream time_id;
		time_id <<"("<<(*cmap_iter).first<<","<<currentT<<")";
		labels->SetValue(return_id, time_id.str());
		sizes->SetValue(return_id, (*cmap_iter).first);
	}

	// Add the labels to the centroid locations:
	labelpolydata->SetPoints(labelpoints);
	labelpolydata->GetPointData()->AddArray(labels);
	labelpolydata->GetPointData()->AddArray(sizes);

	// Generate the label hierarchy.
	vtkSmartPointer<vtkPointSetToLabelHierarchy> pointSetToLabelHierarchyFilter =vtkSmartPointer<vtkPointSetToLabelHierarchy>::New();
	pointSetToLabelHierarchyFilter->SetInputConnection( labelpolydata->GetProducerPort());
	pointSetToLabelHierarchyFilter->SetLabelArrayName("labels");
	pointSetToLabelHierarchyFilter->SetPriorityArrayName("sizes");
	pointSetToLabelHierarchyFilter->Update();

	// Create a mapper and actor for the labels.
	vtkSmartPointer<vtkLabelPlacementMapper> labelMapper = vtkSmartPointer<vtkLabelPlacementMapper>::New();
	labelMapper->SetInputConnection(pointSetToLabelHierarchyFilter->GetOutputPort());
	vtkSmartPointer<vtkActor2D> labelactor = vtkSmartPointer<vtkActor2D>::New();
	labelactor->SetMapper(labelMapper);

	ImageLabelsVector.push_back(labelactor);
}

void Image3DView::SetupRenderWindow()
{
	MainWindow = new QMainWindow();
	MainWindow->setWindowTitle(tr("Contour View:"));

	QVTKView = new QVTKWidget(MainWindow);
	MainWindow->setCentralWidget(QVTKView);

	Renderer = vtkSmartPointer<vtkRenderer>::New();
	Renderer->BackingStoreOff();
	RenderWindow  = QVTKView->GetRenderWindow();
	QVTKView->GetRenderWindow()->AddRenderer(Renderer);
	RenderWindow->Render();

	mode = static_cast<ftk::Image::PtrMode>(2); 
	VTKImage = LabelImageData->GetVtkPtr(0,0,mode); 
	Renderer->ResetCamera(VTKImage->GetBounds());  // change image bounds.
	MainWindow->show();
	QVTKView->show();
	QVTKView->setMinimumSize(LabelImageData->GetImageInfo()->numColumns,LabelImageData->GetImageInfo()->numRows);
}
void Image3DView::Render3DStack(int currentT)
{
	std::map<int, CellActors> actorsmap = CellActorsMap.at(currentT) ;
	std::map<int, CellActors>::iterator act_iter;

	for (act_iter = actorsmap.begin(); act_iter != actorsmap.end() ; act_iter++)
	{
		Renderer->AddActor((*act_iter).second.boundaryActor);					// render the contours
		Renderer->AddActor((*act_iter).second.centroidActor);					// render the centroids
	}
	if(labelsVisible)
		Renderer->AddActor(ImageLabelsVector.at(currentT));						// render the labels
}
void Image3DView::RefreshTimeActors(void)
{
	
	if(stacksVisible)
	{
		this->Renderer->RemoveAllViewProps();
		this->Render3DStack(ImageView->GetCurrentTimeVal());
		this->Renderer->Modified();
		this->RenderWindow->Render();
	}
}
void Image3DView::CreateInteractorStyle(void)
{
	this->Interactor = this->QVTKView->GetRenderWindow()->GetInteractor();
	//keep mouse command observers, but change the key ones
	this->keyPress = vtkSmartPointer<vtkCallbackCommand>::New();
	this->keyPress->SetCallback(HandleKeyPress);
	this->keyPress->SetClientData(this);

	this->Interactor->RemoveObservers(vtkCommand::KeyPressEvent);
	this->Interactor->RemoveObservers(vtkCommand::KeyReleaseEvent);
	this->Interactor->RemoveObservers(vtkCommand::CharEvent);
	this->Interactor->AddObserver(vtkCommand::KeyPressEvent, this->keyPress);
}
void Image3DView::HandleKeyPress(vtkObject* caller, unsigned long event, void* clientdata, void* callerdata)
{
	Image3DView* view = (Image3DView*)clientdata;
	char key = view->Interactor->GetKeyCode();
	switch (key)
	{
	//case 'a':
	//	view->AutomationDock->toggleViewAction();
	//	break;

	case 'l':
		view->ToggleLabelVisibility();
		break;

	case 'v':
		view->Toggle3DStackTrackView();
		break;
	//case 's':
	//	view->SplitTraces();
	//	break;

	//case 'f':
	//	view->FlipTraces();
	//	break;

	//case 'w':
	//	view->SaveToFile();
	//	break;

	//case 't':
	//	view->SelectTrees();
	//	break;
	//case 'q':
	//	view->updateSelectionHighlights();
	//	break;
	//case 'z':
	//	view->HandleHippocampalDataset();
	//	break;
	//default:
	//	break;
	}

}
void Image3DView::ToggleLabelVisibility()
{
	int currentT = ImageView->GetCurrentTimeVal();
	if(labelsVisible) // I only have two actors
	{
		this->Renderer->RemoveActor2D(ImageLabelsVector.at(currentT));
		this->labelsVisible = false;
	}
	else
	{
		this->Renderer->AddActor2D(ImageLabelsVector.at(currentT));
		this->labelsVisible = true;
	}
	this->Renderer->Modified();
	this->RenderWindow->Render();
}
void Image3DView::Toggle3DStackTrackView()
{
	if(stacksVisible)
	{
		this->stacksVisible = false;
		this->RefreshTrackActors();
	}
	else
	{
		this->stacksVisible = true;
		this->RefreshTimeActors();
	}
}
void Image3DView::RefreshTrackActors(void)
{
	if(stacksVisible) return;
	std::set<long int> selected_ids = Selection->getSelections();
	
	if (selected_ids.size()!=0)
	{
		this->Renderer->RemoveAllViewProps();
		std::map<int, CellActors>::iterator act_iter;  
		std::set<long int>::iterator id_iter;

		for(id_iter = selected_ids.begin(); id_iter!=selected_ids.end(); id_iter++)
		{
			for(int time=0 ; time<(int)CellActorsMap.size(); time++)
			{
				std::map<int, CellActors>  cellactors = CellActorsMap.at(time);
				act_iter = cellactors.find((*id_iter));
				if (act_iter != cellactors.end())
				{	
					Renderer->AddActor((*act_iter).second.boundaryActor);					// render the contours
					Renderer->AddActor((*act_iter).second.centroidActor);					// render the centroids
				}
			}
		}
		this->Renderer->Modified();
		this->RenderWindow->Render();
	}
}