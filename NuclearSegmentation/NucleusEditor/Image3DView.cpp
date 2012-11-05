#include "Image3DView.h"

Image3DView::Image3DView(ftk::Image::Pointer image, ftk::Image::Pointer labimage, LabelImageViewQT * imview, ObjectSelection * sels)
{
	if(!image)return;

	this->LabelImageData = labimage;
	this->ImageData = image;
	this->ImageView = imview;
	this->CenterMap = ImageView->GetCenterMapVector();// needs to be fixed later;
	this->bBoxMap = ImageView->GetBoxMapVector();
	this->Selection = sels;

	//Generate render data:
	double before = clock( ) ;
	this->CreateVolumes();
	this->CreateActors();
	double after = clock( ) ;
	double time_elapsed = (after - before) ;
	std::cout << "Time elapsed is " << time_elapsed <<std::endl;	
	
	// Setup the render window, interactor and render:
	this->SetupRenderWindow();
	this->CreateInteractorStyle();
	
	this->GenerateTracks();
	this->CreateTrackPoints();
	this->CreateTracks();
//	this->CreateBoundingBox();
	
	// Initialize Flags:
	this->labelsVisible = true;
	this->stacksVisible = true;
	this->volumesVisible = true;
	this->countoursVisible = true;
	this->tracksVisible = true;
	// Render:
	this->Render3DStack(ImageView->GetCurrentTimeVal());

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
		//centroidSphere->SetCenter((double)(*cmap_iter).second.x,
		//						  (double)(*cmap_iter).second.y,
		//						 // (double)(*cmap_iter).second.z);
		//						  (double)((*cmap_iter).second.z*Z_SPACING));
		centroidSphere->SetRadius(0.3);	
		// Setup Mapper and Actor:
		vtkSmartPointer<vtkPolyDataMapper> centroidMapper =  vtkSmartPointer<vtkPolyDataMapper>::New();
		centroidMapper->SetInputConnection(centroidSphere->GetOutputPort()); 
		vtkSmartPointer<vtkActor> centroidactor = vtkSmartPointer<vtkActor>::New();
		centroidactor->SetMapper(centroidMapper);
		centroidactor->GetProperty()->SetOpacity(1);
		centroidactor->GetProperty()->SetColor(1,0,0);
		centroidactor->SetPosition((double)(*cmap_iter).second.x,
								  (double)(*cmap_iter).second.y,
								 // (double)(*cmap_iter).second.z);
								  (double)((*cmap_iter).second.z*Z_SPACING));
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
		vtkSmartPointer<vtkUnsignedCharArray> cellcolor = vtkSmartPointer<vtkUnsignedCharArray>::New();
		cellcolor->SetNumberOfComponents(3);
		cellcolor->SetName("colors");
		//vtkSmartPointer<vtkFloatArray> cellcolor = vtkSmartPointer<vtkFloatArray>::New();
		//cellcolor->SetNumberOfComponents(1);
		//cellcolor->SetName("colors");
		//float col = (float)(*bbox_iter).first/40.0;
		unsigned char color[3];
		this->GetBoudaryActorColor((*bbox_iter).first,color);
//		for( int i=0; i<(int) contourf->GetOutput()->GetNumberOfCells(); ++i) cellcolor->InsertNextTuple1(col);
		for( int i=0; i<(int) contourf->GetOutput()->GetNumberOfCells(); ++i) cellcolor->InsertNextTupleValue(color);
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
//		contourActor->SetPosition((double)xmin,(double)ymin,(double)zmin);
		contourActor->SetPosition((double)xmin,(double)ymin,(double)(zmin*Z_SPACING));
//		contourActor->SetPosition((double)xc,(double)yc,(double)(zc*Z_SPACING));
		//printf("my bbx: %lf, %lf\n",(double)xmin,(double)xmax);
		//printf("my bbx: %lf, %lf\n",(double)ymin,(double)ymax);
		//printf("my bbx: %lf, %lf\n",(double)zmin,(double)zmax);

		double b[6];
		contourActor->GetBounds(b);
		//printf("their bbx: %lf, %lf\n",(double)b[0],(double)b[1]);
		//printf("their bbx: %lf, %lf\n",(double)b[2],(double)b[3]);
		//printf("their bbx: %lf, %lf\n",(double)b[4],(double)b[5]);

	
		//double pos[3];
		//contourActor->GetPosition(pos);
		//printf("x pos set: %lf\n",xc);
		//printf("y pos set: %lf\n",yc);
		//printf("z pos set: %lf\n",(double)zc*Z_SPACING);
		//printf("x pos: %lf\n",pos[0]);
		//printf("y pos: %lf\n",pos[1]);
		//printf("z pos: %lf\n",pos[2]);

		CellActors act;
		act.boundaryActor = contourActor;
		actorsmap.insert(std::pair<int,CellActors>((*bbox_iter).first,act));
	}
	CellActorsMap.push_back(actorsmap);
}
void Image3DView::GetBoudaryActorColor(int id, unsigned char color[3])
{
	int color_index = id%6;
	switch(color_index)
	{
	//Cyan
	case 0:
	    	color[0]=0;
			color[1]=255;
			color[2]=255;
			
					 break;
	//Royal Blue 	65-105-225
	case 1:
	    	color[0]=65;
			color[1]=105;
			color[2]=255;
			
					 break;
	//Red
	case 2:
	    	color[0]=255;
			color[1]=0;
			color[2]=0;
			
					 break;
	//Blue
	case 3:
	    	color[0]=0;
			color[1]=0;
			color[2]=255;
			
					 break;
	
	//Orange 	255-165-0
	case 4:
	    	color[0]=255;
			color[1]=165;
			color[2]=0;
			
					 break;
	//Violet 	238-130-238
	case 5:
	    	color[0]=255;
			color[1]=255;
			color[2]=0;
			 
					 break;
	}

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
													// (double)(*cmap_iter).second.z+labeloffset);
													 (double)((*cmap_iter).second.z*Z_SPACING+labeloffset));
		//Labels:
		std::stringstream time_id;
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
	pointSetToLabelHierarchyFilter->GetTextProperty()->SetColor(1.0,1.0,1.0);
	//pointSetToLabelHierarchyFilter->GetTextProperty()->SetFontSize (15);
	//pointSetToLabelHierarchyFilter->GetTextProperty()->SetBold(2);
	pointSetToLabelHierarchyFilter->Update();

	// Create a mapper and actor for the labels.
	vtkSmartPointer<vtkLabelPlacementMapper> labelMapper = vtkSmartPointer<vtkLabelPlacementMapper>::New();
	labelMapper->SetInputConnection(pointSetToLabelHierarchyFilter->GetOutputPort());
	vtkSmartPointer<vtkActor2D> labelactor = vtkSmartPointer<vtkActor2D>::New();
	labelactor->SetMapper(labelMapper);

	ImageLabelsVector.push_back(labelactor);
}
void Image3DView::CreateVolumes(void)
{
	float colors[3]={1.0,1.0,1.0};
	mode = static_cast<ftk::Image::PtrMode>(2); 

	for(int T = 0; T< ImageData->GetImageInfo()->numTSlices; ++T)
	{
		vtkSmartPointer<vtkImageData> vtkimage = vtkSmartPointer<vtkImageData>::New();
		vtkimage = ImageData->GetVtkPtr(T,0,mode);									// first channel is assumed to be nuclear channel
		vtkimage->SetSpacing(1.0,1.0,Z_SPACING);
		vtkimage->Update();
		vtkSmartPointer<vtkPiecewiseFunction> opacityTransferFunction = vtkSmartPointer<vtkPiecewiseFunction>::New();
		opacityTransferFunction->AddPoint(2,0.0);
		opacityTransferFunction->AddPoint(50,0.8);
		//opacityTransferFunction->AddPoint(55,0.0);
		//opacityTransferFunction->AddPoint(255,1);

		vtkSmartPointer<vtkVolumeProperty> volumeProperty = vtkSmartPointer<vtkVolumeProperty>::New();
		vtkSmartPointer<vtkColorTransferFunction> colorTransferFunction = vtkSmartPointer<vtkColorTransferFunction>::New();
		//colorTransferFunction->AddRGBPoint(5.0,0.0,0.0,0.0);
		colorTransferFunction->AddRGBPoint(0.0,0.0,0.0,0.0);
		//colorTransferFunction->AddRGBPoint(50.0,colors[0],colors[1],colors[2]);
		colorTransferFunction->AddRGBPoint(255.0,colors[0],colors[1],colors[2]);
		volumeProperty->SetColor(colorTransferFunction);
		volumeProperty->SetScalarOpacity(opacityTransferFunction);
		volumeProperty->ShadeOff();
		volumeProperty->SetInterpolationTypeToLinear();
		volumeProperty->DisableGradientOpacityOn();
#ifndef WIN32
		vtkSmartPointer<vtkOpenGLVolumeTextureMapper2D> volumeMapper = vtkSmartPointer<vtkOpenGLVolumeTextureMapper2D>::New();
		volumeMapper->SetMaximumNumberOfPlanes(50);
#else
		//vtkSmartPointer<vtkOpenGLVolumeTextureMapper3D> volumeMapper = vtkSmartPointer<vtkOpenGLVolumeTextureMapper3D>::New();
		vtkSmartPointer<vtkGPUVolumeRayCastMapper> volumeMapper = vtkSmartPointer<vtkGPUVolumeRayCastMapper>::New();
		
		//volumeMapper->SetPreferredMethodToNVidia();
		volumeMapper->SetSampleDistance(1);
#endif
		volumeMapper->SetInput(vtkimage);
		volumeMapper->Update();
		vtkSmartPointer<vtkVolume> volume = vtkSmartPointer<vtkVolume>::New();
		volume->SetMapper(volumeMapper);
		volume->SetProperty(volumeProperty);
		volume->Update();
		printf("Creating Volume at T: %d\n",T);
		ImageVolumes.push_back(volume);
	}
}
void Image3DView::SetupRenderWindow(void)
{
	MainWindow = new QMainWindow();
	MainWindow->setWindowTitle(tr("5D Viewer:"));

	QVTKView = new QVTKWidget(MainWindow);
	MainWindow->setCentralWidget(QVTKView);
	// Opacity Slider:
	DockWidget = new QDockWidget(MainWindow);
	DockWidget->setWindowTitle(tr(""));
	DockWidget->setAllowedAreas(Qt::RightDockWidgetArea\
								|Qt::LeftDockWidgetArea|Qt::TopDockWidgetArea|Qt::BottomDockWidgetArea);

	opacitySlider = new QSlider();
	opacitySlider->setOrientation(Qt::Horizontal);
	opacitySlider->setMaximum(20);
	opacitySlider->setSingleStep(1);
	opacitySlider->setPageStep(1);
	opacitySlider->setTickInterval(1);
	opacitySlider->setTickPosition(QSlider::TicksBelow);

	opacityLabel = new QLabel("Opacity");
	opacityLabel->setDisabled(false);

	brightnessSlider = new QSlider();
	brightnessSlider->setOrientation(Qt::Horizontal);
	brightnessSlider->setMaximum(25);
	brightnessSlider->setSingleStep(1);
	brightnessSlider->setPageStep(1);
	brightnessSlider->setTickInterval(1);
	brightnessSlider->setTickPosition(QSlider::TicksBelow);

	brightnessLabel = new QLabel("Brightness");
	brightnessLabel->setDisabled(false);

	QHBoxLayout * sliderLayout = new QHBoxLayout;
	sliderLayout->addWidget(opacityLabel);
	sliderLayout->addWidget(opacitySlider);
	sliderLayout->addWidget(brightnessLabel);
	sliderLayout->addWidget(brightnessSlider);

	QWidget * widget = new QWidget;
    widget->setLayout(sliderLayout);
	DockWidget->setWidget(widget);
	MainWindow->addDockWidget(Qt::BottomDockWidgetArea,DockWidget);

	connect(opacitySlider, SIGNAL(valueChanged(int)), this, SLOT(ChangeOpacity(void)));
	connect(brightnessSlider, SIGNAL(valueChanged(int)), this, SLOT(ChangeBrightness(void)));


	// Renderer and Render Window:
	Renderer = vtkSmartPointer<vtkRenderer>::New();
	Renderer->BackingStoreOff();
	Renderer->SetBackground(0.0,0.0,0.0);
//	Renderer->SetBackground(0.8,0.8,0.8);
	RenderWindow  = QVTKView->GetRenderWindow();
	QVTKView->GetRenderWindow()->AddRenderer(Renderer);
	RenderWindow->Render();
   	Renderer->ResetCamera(ImageVolumes.at(0)->GetBounds());  // change image bounds.
	// Show everything:
	MainWindow->show();
	//opacitySlider->show();
	QVTKView->show();
	QVTKView->setMinimumSize(ImageData->GetImageInfo()->numColumns,ImageData->GetImageInfo()->numRows);

}
void Image3DView::ChangeOpacity(void)
{
	//int opacityValue = opacitySlider->value();
	//opacityValue = (opacityValue*10)+55;
	double opacityValue = (double)opacitySlider->value();
	opacityValue /= 20.0;
	vtkSmartPointer<vtkPiecewiseFunction> opacityTransferFunction = vtkSmartPointer<vtkPiecewiseFunction>::New();
	//opacityTransferFunction->AddPoint(50,0.0);
	//opacityTransferFunction->AddPoint(255,opacityValue);
	opacityTransferFunction->AddPoint(2,0.0);
	opacityTransferFunction->AddPoint(50,opacityValue);
	for(int i=0; i<ImageVolumes.size(); ++i)
	{
		ImageVolumes.at(i)->GetProperty()->SetScalarOpacity(opacityTransferFunction);
		ImageVolumes.at(i)->Modified();
	}
	this->Renderer->Modified();
	this->RenderWindow->Render();
}
void Image3DView::ChangeBrightness(void)
{
	//int brightnessValue = brightnessSlider->value();
	//brightnessValue = 255-((brightnessValue*10)+5);
	double brightnessValue = (double)brightnessSlider->value();
	brightnessValue /=25.0;
	vtkColorTransferFunction *colorTransferFunction = vtkColorTransferFunction::New();
	colorTransferFunction->AddRGBPoint(0.0, 0.0, 0.0, 0.0);
	colorTransferFunction->AddRGBPoint(255*(1-brightnessValue),1,1,1);
//	colorTransferFunction->AddRGBPoint(brightnessValue,1,1,1);

	for(int i=0; i<ImageVolumes.size(); ++i)
	{
		ImageVolumes.at(i)->GetProperty()->SetColor(colorTransferFunction);
		ImageVolumes.at(i)->Modified();
	}
	this->Renderer->Modified();
	this->RenderWindow->Render();
}

void Image3DView::CreateBoundingBox(void)
{
	vtkSmartPointer<vtkPolyData> boxpolydata = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPoints> boxpoints = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> boxcells = vtkSmartPointer<vtkCellArray>::New();

	double bbox[6];
	this->ImageVolumes.at(0)->GetBounds(bbox);
	for(int z = 4; z < 6; ++z)
	{
		int x = 0;
		int y = x+2;
		unsigned int id1,id2;
		id1 = boxpoints->InsertNextPoint(bbox[x],bbox[y],bbox[z]);
		id2 = boxpoints->InsertNextPoint(bbox[x],bbox[y],bbox[z]);
		boxcells->InsertNextCell(2);
		boxcells->InsertCellPoint(id1);
		boxcells->InsertCellPoint(id2);
		id2 = boxpoints->InsertNextPoint(bbox[x],bbox[y+1],bbox[z]);			
		boxcells->InsertNextCell(2);
		boxcells->InsertCellPoint(id1);
		boxcells->InsertCellPoint(id2);
		id1 = boxpoints->InsertNextPoint(bbox[x+1],bbox[y+1],bbox[z]);
		boxcells->InsertNextCell(2);
		boxcells->InsertCellPoint(id1);
		boxcells->InsertCellPoint(id2);

		x += 1;
		y = x+2;
		id1 = boxpoints->InsertNextPoint(bbox[x],bbox[y],bbox[z]);
		id2 = boxpoints->InsertNextPoint(bbox[x],bbox[y+1],bbox[z]);
		boxcells->InsertNextCell(2);
		boxcells->InsertCellPoint(id1);
		boxcells->InsertCellPoint(id2);
	}
	boxpolydata->SetPoints(boxpoints);
	boxpolydata->SetLines(boxcells);
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInput(boxpolydata);
	boxActor = vtkSmartPointer<vtkActor>::New();
	boxActor->SetMapper(mapper);
	boxActor->GetProperty()->SetLineWidth(3.0);
	boxActor->GetProperty()->SetColor(1.0,0.0,0.0);
	this->Renderer->AddActor(boxActor);
	this->Renderer->Modified();
	this->RenderWindow->Render();	

}
void Image3DView::GenerateTracks(void)
{
	printf("Started CreateTracks\n");
	int max_track_num = 0;
	for(int i = 0; i< CenterMap.size(); ++i)			// iterate through time
	{
		std::map<int,ftk::Object::Point>::iterator id_iter; 
		for(id_iter = CenterMap.at(i).begin(); id_iter != CenterMap.at(i).end();++id_iter)
		{
			max_track_num = MAX(id_iter->first,max_track_num);
		}
	}
	printf("Computed max_track_num = %d\n",max_track_num);
	printf("About to start creating tobj\n");
	TraceData = new TraceObject();
	for(int id = 0; id<= max_track_num; ++id)
	{
		TraceLine *tline = new TraceLine();
		bool once = false;
		for(int i = 0; i< CenterMap.size(); ++i)
		{
			std::map<int,ftk::Object::Point>::iterator id_iter = CenterMap[i].find(id); 
			if(id_iter!=CenterMap[i].end())
			{
				TraceBit tbit;
				tbit.x = id_iter->second.x;
				tbit.y = id_iter->second.y;
				tbit.z = id_iter->second.z*Z_SPACING;
				tbit.time = i;
					
				tbit.id = id_iter->first;
				printf("About to add TraceBit\n");
				tline->AddTraceBit(tbit);
				once = true;
			}
		}

		if(once)
		{
			TraceData->GetTraceLinesPointer()->push_back(tline);
			printf("Added one tline\n");
		}
		else
			delete tline;
	}
	printf("Finished CreateTracks() \n");
	return;
}
void Image3DView::CreateTrackPoints(void)
{
	std::vector<TraceBit> tbits = TraceData->CollectTraceBits();
	vtkSmartPointer<vtkPolyData> point_poly = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();

	for(int counter=0; counter<tbits.size(); counter++)
	{
		int return_id = points->InsertNextPoint(tbits[counter].x,tbits[counter].y,tbits[counter].z);
		cells->InsertNextCell(1);
		cells->InsertCellPoint(return_id);
	}

	printf("About to create poly\n");
	point_poly->SetPoints(points);
	point_poly->SetVerts(cells);

   // Mapper:
	vtkSmartPointer<vtkPolyDataMapper> cubemap = vtkSmartPointer<vtkPolyDataMapper>::New();
	cubemap->SetInput(point_poly);
	cubemap->GlobalImmediateModeRenderingOn();

	// Actor(OpenGL):
	bitsActor = vtkSmartPointer<vtkActor>::New();
	bitsActor->SetMapper(cubemap);
	bitsActor->PickableOff();
	bitsActor->GetProperty()->SetPointSize(1);
	bitsActor->GetProperty()->SetOpacity(1);
	bitsActor->GetProperty()->SetColor(1,1,1);
	this->Renderer->AddActor(bitsActor);
	this->Renderer->Modified();
	this->RenderWindow->Render();	

}

void Image3DView::CreateTracks(void)
{
	TrackPoly =  TraceData->GetVTKPolyData();
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInput(TrackPoly);
	tracksActor = vtkSmartPointer<vtkActor>::New();
	tracksActor->SetMapper(mapper);
	tracksActor->GetProperty()->SetLineWidth(3.0);
	this->Renderer->AddActor(tracksActor);
	this->Renderer->Modified();
	this->RenderWindow->Render();	
}

void Image3DView::Render3DStack(int currentT)
{
	if(countoursVisible)
	{
		std::map<int, CellActors> actorsmap = CellActorsMap.at(currentT) ;
		std::map<int, CellActors>::iterator act_iter;

		for (act_iter = actorsmap.begin(); act_iter != actorsmap.end() ; act_iter++)
		{
			Renderer->AddActor((*act_iter).second.boundaryActor);					// render the contours
		//	//Renderer->AddActor((*act_iter).second.centroidActor);					// render the centroids
		//	
		}
	}
	if(labelsVisible)
		Renderer->AddActor(ImageLabelsVector.at(currentT));						// render the labels
	if(volumesVisible)
		Renderer->AddVolume(ImageVolumes.at(currentT));

}
void Image3DView::RefreshTimeActors(void)
{
	
	if(stacksVisible)
	{
		this->Renderer->RemoveAllViewProps();
		this->Render3DStack(ImageView->GetCurrentTimeVal());
		if(tracksVisible)
		{
			this->Renderer->AddActor(bitsActor);
			this->Renderer->AddActor(tracksActor);
		}
		//this->Renderer->AddActor(boxActor);
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
	case 'l':
		view->ToggleLabelVisibility();
		break;
	case 'd':
		view->Toggle3DStackTrackView();
		break;
	case 'v':
		view->ToggleVolumeVisibility();
		break;
	case 's':
		view->ToggleContourVisibility();
		break;
	case 't':
		view->ToggleTrackVisibility();
		break;
	}

}
void Image3DView::ToggleTrackVisibility()
{
	if(tracksVisible) 
	{
		this->Renderer->RemoveActor(bitsActor);
		this->Renderer->RemoveActor(tracksActor);
		this->tracksVisible = false;
	}
	else
	{
		this->Renderer->AddActor(bitsActor);
		this->Renderer->AddActor(tracksActor);
		this->tracksVisible = true;
	}
	this->Renderer->Modified();
	this->RenderWindow->Render();
}

void Image3DView::ToggleContourVisibility()
{
	int currentT = ImageView->GetCurrentTimeVal();
	std::map<int, CellActors> actorsmap = CellActorsMap.at(currentT) ;
	std::map<int, CellActors>::iterator act_iter;
	if(countoursVisible) 
	{
		for (act_iter = actorsmap.begin(); act_iter != actorsmap.end() ; act_iter++)
		{
			Renderer->RemoveActor((*act_iter).second.boundaryActor);					// render the contours
		}
		this->countoursVisible = false;
	}
	else
	{
		for (act_iter = actorsmap.begin(); act_iter != actorsmap.end() ; act_iter++)
		{
			Renderer->AddActor((*act_iter).second.boundaryActor);					// render the contours
		}
		this->countoursVisible = true;
	}
	this->Renderer->Modified();
	this->RenderWindow->Render();
}

void Image3DView::ToggleVolumeVisibility()
{
	int currentT = ImageView->GetCurrentTimeVal();
	if(volumesVisible) 
	{
		this->Renderer->RemoveVolume(ImageVolumes.at(currentT));
		this->volumesVisible = false;
	}
	else
	{
		this->Renderer->AddVolume(ImageVolumes.at(currentT));
		this->volumesVisible = true;
	}
	this->Renderer->Modified();
	this->RenderWindow->Render();
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
