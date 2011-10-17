#include "KymoGraphView.h"

void TrackingKymoView::GenerateImages()
{
	float col[3]={1.0,1.0,1.0};
	m_vtkim = vtkSmartPointer<vtkImageData>::New();
//	InputImageType::SizeType imsize = m_model->getRawImagePointer(0)->GetLargestPossibleRegion().GetSize();
	unsigned short xsize = (my4DImg->GetImageInfo()->numColumns);
	unsigned short ysize = (my4DImg->GetImageInfo()->numRows);
	unsigned short tsize = (my4DImg->GetImageInfo()->numTSlices);

	m_vtkim->SetExtent(0,xsize-1,0,ysize-1,0,tsize-1);
	m_vtkim->SetScalarTypeToUnsignedChar();
	m_vtkim->AllocateScalars();
	m_vtkim->SetSpacing(1.0,1.0,5.0);

	//Loop Constants:
	unsigned short xy = xsize*ysize;
	ftk::Image::PtrMode mode;
	mode = static_cast<ftk::Image::PtrMode>(0); 

/*	for(unsigned short t=0; t< tsize; t++)
	{
		Input2DImageType::Pointer my2DImage = getProjection(my4DImg->GetItkPtr<InputPixelType>(t,0,mode));
		memcpy(static_cast<unsigned char*>(m_vtkim->GetScalarPointer())+xsize*ysize*t,my2DImage->GetBufferPointer(),sizeof(unsigned char)*xsize*ysize);
	}*/

	Input2DImageType::Pointer my2DImage = Input2DImageType::New();
	Input2DImageType::RegionType region;
	Input2DImageType::SizeType size;
	Input2DImageType::IndexType index;
	index[0]=0;
	index[1]=0;
	size[0] = xsize;
	size[1] = ysize;
	region.SetSize(size);
	region.SetIndex(index);
	my2DImage->SetRegions(region);
	my2DImage->Allocate();

	for(unsigned short t=0; t< tsize; t++)
	{
		SliceIteratorType inputIt(my4DImg->GetItkPtr<InputPixelType>(t,0,mode),my4DImg->GetItkPtr<InputPixelType>(t,0,mode)->GetLargestPossibleRegion());
		LinearIteratorType outputIt(my2DImage,my2DImage->GetLargestPossibleRegion());

		inputIt.SetFirstDirection(0);
		inputIt.SetSecondDirection(1);
		outputIt.SetDirection(0);
		outputIt.GoToBegin();
		while ( ! outputIt.IsAtEnd() )
		{
			while ( ! outputIt.IsAtEndOfLine() )
			{
				outputIt.Set( itk::NumericTraits<unsigned short>::NonpositiveMin() );
				++outputIt;
			}
			outputIt.NextLine();
		}

		inputIt.GoToBegin();
		outputIt.GoToBegin();
		while( !inputIt.IsAtEnd() )
		{
			while ( !inputIt.IsAtEndOfSlice() )
			{
				while ( !inputIt.IsAtEndOfLine() )
				{
					outputIt.Set( MAX( outputIt.Get(), inputIt.Get() ));
					++inputIt;
					++outputIt;
				}
				outputIt.NextLine();
				inputIt.NextLine();
			}
			outputIt.GoToBegin();
			inputIt.NextSlice();
		}

		memcpy(static_cast<unsigned char*>(m_vtkim->GetScalarPointer())+xsize*ysize*t,my2DImage->GetBufferPointer(),sizeof(unsigned char)*xsize*ysize);
	}
	
	m_vtkvolume = getOneVTKVolume(m_vtkim,col);
	m_vtkvolume->Print(std::cout);
	m_vtkrenderer->AddVolume(m_vtkvolume);
	m_vtkrenderer->ResetCamera(m_vtkim->GetBounds());
}

/*Input2DImageType::Pointer TrackingKymoView::getProjection(InputImageType::Pointer im)
{
	Input2DImageType::Pointer output = Input2DImageType::New();
	Input2DImageType::RegionType region;
	Input2DImageType::SizeType size;
	Input2DImageType::IndexType index;
	index[0]=0;
	index[1]=0;
	size[0] = im->GetLargestPossibleRegion().GetSize()[0];
	size[1] = im->GetLargestPossibleRegion().GetSize()[1];
	region.SetSize(size);
	region.SetIndex(index);
	output->SetRegions(region);
	output->Allocate();

	SliceIteratorType inputIt(im,im->GetLargestPossibleRegion());
	LinearIteratorType outputIt(output,output->GetLargestPossibleRegion());

	inputIt.SetFirstDirection(0);
	inputIt.SetSecondDirection(1);
	outputIt.SetDirection(0);



	outputIt.GoToBegin();
	while ( ! outputIt.IsAtEnd() )
	{
		while ( ! outputIt.IsAtEndOfLine() )
		{
			outputIt.Set( itk::NumericTraits<unsigned short>::NonpositiveMin() );
			++outputIt;
		}
		outputIt.NextLine();
	}

	inputIt.GoToBegin();
	outputIt.GoToBegin();
	while( !inputIt.IsAtEnd() )
	{
		while ( !inputIt.IsAtEndOfSlice() )
		{
			while ( !inputIt.IsAtEndOfLine() )
			{
				outputIt.Set( MAX( outputIt.Get(), inputIt.Get() ));
				++inputIt;
				++outputIt;
			}
			outputIt.NextLine();
			inputIt.NextLine();
		}
		outputIt.GoToBegin();
		inputIt.NextSlice();

	}
	return output;
}*/
vtkSmartPointer<vtkVolume> TrackingKymoView::getOneVTKVolume(vtkSmartPointer<vtkImageData> vtkim, float colors[3])
{

	vtkSmartPointer<vtkPiecewiseFunction> opacityTransferFunction = vtkSmartPointer<vtkPiecewiseFunction>::New();
	opacityTransferFunction->AddPoint(2,0.0);
	opacityTransferFunction->AddPoint(50,0.8);

	vtkSmartPointer<vtkVolumeProperty> volumeProperty = vtkSmartPointer<vtkVolumeProperty>::New();

	
	vtkSmartPointer<vtkColorTransferFunction> colorTransferFunction = vtkSmartPointer<vtkColorTransferFunction>::New();
	colorTransferFunction->AddRGBPoint(0.0,0.0,0.0,0.0);
	colorTransferFunction->AddRGBPoint(50.0,colors[0],colors[1],colors[2]);
	volumeProperty->SetColor(colorTransferFunction);
	volumeProperty->SetScalarOpacity(opacityTransferFunction);
	volumeProperty->ShadeOff();
	volumeProperty->SetInterpolationTypeToLinear();
	volumeProperty->DisableGradientOpacityOn();
	


#ifndef WIN32
	vtkSmartPointer<vtkOpenGLVolumeTextureMapper2D> vMapper = vtkSmartPointer<vtkOpenGLVolumeTextureMapper2D>::New();
	vMapper->SetMaximumNumberOfPlanes(50);
#else
	/*vtkSmartPointer<vtkVolumeRayCastMapper> vMapper = vtkSmartPointer<vtkVolumeRayCastMapper>::New();
	vtkSmartPointer<vtkVolumeRayCastCompositeFunction> volume_ray_cast = vtkSmartPointer<vtkVolumeRayCastCompositeFunction>::New();
	vMapper->SetVolumeRayCastFunction(volume_ray_cast);
	//vMapper->SetBlendModeToMaximumIntensity();
	//vMapper->SetMaximumImageSampleDistance(20);
	//vMapper->SetMinimumImageSampleDistance(1);*/
    
	vtkSmartPointer<vtkOpenGLVolumeTextureMapper3D> vMapper = vtkSmartPointer<vtkOpenGLVolumeTextureMapper3D>::New();
	vMapper->SetPreferredMethodToNVidia();
	vMapper->SetSampleDistance(1);
	
#endif
	vMapper->SetInput(vtkim);
	vMapper->Update();
	vtkSmartPointer<vtkVolume> volume = vtkSmartPointer<vtkVolume>::New();
	volume->SetMapper(vMapper);
	volume->SetProperty(volumeProperty);
	volume->Update();
	volume->Print(std::cout);

	printf("______________________________________\n");

	return volume;
}
void TrackingKymoView::AddSliders()
{
  vtkSliderRepresentation2D *sliderRep = vtkSliderRepresentation2D::New();
  sliderRep->SetValue(0.8);
  sliderRep->SetTitleText("Opacity");
  sliderRep->GetPoint1Coordinate()->SetCoordinateSystemToNormalizedDisplay();
  sliderRep->GetPoint1Coordinate()->SetValue(0.2,0.1);
  sliderRep->GetPoint2Coordinate()->SetCoordinateSystemToNormalizedDisplay();
  sliderRep->GetPoint2Coordinate()->SetValue(0.8,0.1);
  sliderRep->SetSliderLength(0.02);
  sliderRep->SetSliderWidth(0.03);
  sliderRep->SetEndCapLength(0.01);
  sliderRep->SetEndCapWidth(0.03);
  sliderRep->SetTubeWidth(0.005);
  sliderRep->SetMinimumValue(0.0);
  sliderRep->SetMaximumValue(1.0);

  vtkSliderWidget *sliderWidget = vtkSliderWidget::New();
  sliderWidget->SetInteractor(m_imageview->GetRenderWindow()->GetInteractor());
  sliderWidget->SetRepresentation(sliderRep);
  sliderWidget->SetAnimationModeToAnimate();

  m_callback_brightness = vtkSlider2DKymoCallbackBrightness::New();
  m_callback_brightness->volume = m_vtkvolume;
  sliderWidget->AddObserver(vtkCommand::InteractionEvent,m_callback_brightness);
  sliderWidget->EnabledOn();
  

// slider 2

  vtkSliderRepresentation2D *sliderRep2 = vtkSliderRepresentation2D::New();
  sliderRep2->SetValue(0.8);
  sliderRep2->SetTitleText("Brightness");
  sliderRep2->GetPoint1Coordinate()->SetCoordinateSystemToNormalizedDisplay();
  sliderRep2->GetPoint1Coordinate()->SetValue(0.2,0.9);
  sliderRep2->GetPoint2Coordinate()->SetCoordinateSystemToNormalizedDisplay();
  sliderRep2->GetPoint2Coordinate()->SetValue(0.8,0.9);
  sliderRep2->SetSliderLength(0.02);
  sliderRep2->SetSliderWidth(0.03);
  sliderRep2->SetEndCapLength(0.01);
  sliderRep2->SetEndCapWidth(0.03);
  sliderRep2->SetTubeWidth(0.005);
  sliderRep2->SetMinimumValue(0.0);
  sliderRep2->SetMaximumValue(1.0);

  vtkSliderWidget *sliderWidget2 = vtkSliderWidget::New();
  sliderWidget2->SetInteractor(m_imageview->GetRenderWindow()->GetInteractor());
  sliderWidget2->SetRepresentation(sliderRep2);
  sliderWidget2->SetAnimationModeToAnimate();

  m_callback_contrast = vtkSlider2DKymoCallbackContrast::New();
  m_callback_contrast->volume = m_vtkvolume;
  sliderWidget2->AddObserver(vtkCommand::InteractionEvent,m_callback_contrast);
  sliderWidget2->EnabledOn();

}
void TrackingKymoView::GenerateTracks()
{
	printf("Started GenerateTracks\n");
//	std::vector<std::vector<ftk::LabelImageFeatures>> *f = m_model->getFeatures();
	std::vector< std::vector<ftk::IntrinsicFeatures> > *f = &myfeatures;
	int max_track_num = 0;
	for(int cox=0; cox< f->size(); cox++)
	{
		for(int coy =0; coy < (*f)[cox].size(); coy ++)
		{
			max_track_num = MAX((*f)[cox][coy].num,max_track_num);
		}
	}
	printf("Computed max_track_num = %d\n",max_track_num);
	double time_offset  = 5;//(m_model->getRawImagePointer(0)->GetLargestPossibleRegion().GetSize()[2])*3;//*((*f)[0][0].spacing[2]);


	printf("About to start creating tobj\n");
	m_tobj = new TraceObject();
	for(int counter= 0; counter<= max_track_num; counter++)
	{
		TraceLine *tline = new TraceLine();
		bool once = false;
		for(int cox=0; cox< f->size(); cox++)
		{
			for(int coy =0; coy < (*f)[cox].size(); coy ++)
			{
				if((*f)[cox][coy].num == counter)
				{
					TraceBit tbit;
					tbit.x = (*f)[cox][coy].Centroid[0]*1.0;
					tbit.y = (*f)[cox][coy].Centroid[1]*1.0;
					tbit.z = ((*f)[cox][coy].time)*(time_offset);
					tbit.time = (*f)[cox][coy].time;
					tbit.id = (*f)[cox][coy].num;
					printf("About to add TraceBit\n");
					tline->AddTraceBit(tbit);
					once = true;
				}
			}
		}
		if(once)
		{
			m_tobj->GetTraceLinesPointer()->push_back(tline);
			printf("Added one tline\n");
		}
		else
		{
			delete tline;
		}
	}
	printf("Finished GenerateTracks() \n");
	return;
} 

vtkSmartPointer<vtkActor> TrackingKymoView::getTrackPoints(std::vector<TraceBit> vec)
{
	vtkSmartPointer<vtkPolyData> point_poly = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkStringArray> labels = vtkSmartPointer<vtkStringArray>::New();	// Labels Array
	labels->SetNumberOfValues(vec.size());
	labels->SetName("labels");
	vtkSmartPointer<vtkIntArray> sizes = vtkSmartPointer<vtkIntArray>::New();			// Priority Arrays used by vtkPointSetToLabelHierarchy
	sizes->SetNumberOfValues(vec.size());
	sizes->SetName("sizes");

	for(int counter=0; counter<vec.size(); counter++)
	{
		int return_id = points->InsertNextPoint(vec[counter].x,vec[counter].y,vec[counter].z);
		cells->InsertNextCell(1);
		cells->InsertCellPoint(return_id);
		//Labels:
		stringstream strId;
		strId <<"("<<vec[counter].id<<","<<vec[counter].time<<")";
		labels->SetValue(return_id, strId.str());
		sizes->SetValue(return_id, vec[counter].id);
	}

	printf("About to create poly\n");
	point_poly->SetPoints(points);
	point_poly->SetVerts(cells);
	point_poly->GetPointData()->AddArray(labels);
	point_poly->GetPointData()->AddArray(sizes);

   // Mapper:
	vtkSmartPointer<vtkPolyDataMapper> cubemap = vtkSmartPointer<vtkPolyDataMapper>::New();
	cubemap->SetInput(point_poly);
	cubemap->GlobalImmediateModeRenderingOn();

	// Actor(OpenGL):
	vtkSmartPointer<vtkActor> cubeact = vtkSmartPointer<vtkActor>::New();
	cubeact->SetMapper(cubemap);
	//cubeact->SetPickable(0);
	cubeact->PickableOn();
	cubeact->GetProperty()->SetPointSize(5);
	cubeact->GetProperty()->SetOpacity(.5);

	// Generate the label hierarchy.
	vtkSmartPointer<vtkPointSetToLabelHierarchy> pointSetToLabelHierarchyFilter =vtkSmartPointer<vtkPointSetToLabelHierarchy>::New();
	pointSetToLabelHierarchyFilter->SetInputConnection( point_poly->GetProducerPort());
	pointSetToLabelHierarchyFilter->SetLabelArrayName("labels");
	pointSetToLabelHierarchyFilter->SetPriorityArrayName("sizes");
	pointSetToLabelHierarchyFilter->Update();

	// Create a mapper and actor for the labels.
	vtkSmartPointer<vtkLabelPlacementMapper> labelMapper = vtkSmartPointer<vtkLabelPlacementMapper>::New();
	labelMapper->SetInputConnection(pointSetToLabelHierarchyFilter->GetOutputPort());
	labelActor = vtkSmartPointer<vtkActor2D>::New();
	labelActor->SetMapper(labelMapper);


	// Add axes:
	double bounds[6];
	point_poly->GetBounds(bounds);
	
	cubeAxesActor = vtkSmartPointer<vtkCubeAxesActor2D>::New();
	cubeAxesActor->SetBounds(bounds);
	cubeAxesActor->SetCamera(m_vtkrenderer->GetActiveCamera());
	cubeAxesActor->SetFontFactor(0.5);
	cubeAxesActor->ScalingOff();
	cubeAxesActor->SetXLabel("x"); 
	cubeAxesActor->SetYLabel("y"); 
	cubeAxesActor->SetZLabel("t");
	cubeAxesActor->PickableOn();
	cubeAxesActor->GetXAxisActor2D()->LabelVisibilityOff();
	cubeAxesActor->GetYAxisActor2D()->LabelVisibilityOff();
	cubeAxesActor->GetZAxisActor2D()->LabelVisibilityOff();


	m_vtkrenderer->AddActor(cubeAxesActor);
	m_vtkrenderer->AddActor(labelActor);
	labelsVisible = true;

	return cubeact;

}
void TrackingKymoView::SaveAnimation()
{
	vtkCamera *vcam = m_vtkrenderer->GetActiveCamera();
	for(int co = 0; co < 100; co++)
	{
		vcam->Roll(1.0);
		vcam->Elevation(0.5);
		vcam->Azimuth(0.2);
		m_vtkrenderer->Render();
	}
}
void TrackingKymoView::refreshImages()
{



}
void TrackingKymoView::refreshSelection(void)
{
	 
	currentTime = ImageView->GetCurrentTimeVal();
	std::set<long int> selected_ids = Selection->getSelections();
	std::set<long int>::iterator it;
	for(it=selected_ids.begin() ; it != selected_ids.end(); it++ )
		std::cout<<*it<<std::endl;
	it=selected_ids.begin();
	if (selected_ids.size()==1)
	{
		if(!singleTracksVisible)
		{
			m_vtkrenderer->AddActor(m_trackactor);
			m_vtkrenderer->RemoveActor(currentTrackActor);

			for(int counter=0; counter< TraceBitVector.size(); counter++)
			{
				if ( TraceBitVector[counter].time == currentTime && TraceBitVector[counter].id == (int) *it)
				{
					this->ResetOriginalColors();
					double bounds[6];
					pointerWidget3d->SetPosition(TraceBitVector[counter].x,TraceBitVector[counter].y,TraceBitVector[counter].z);
					pointerWidget3d->SetEnabled(1);
					bounds[0] = TraceBitVector[counter].x - 70;
					bounds[1] = TraceBitVector[counter].x + 70;
					bounds[2] = TraceBitVector[counter].y - 70;
					bounds[3] = TraceBitVector[counter].y + 70;
					bounds[4] = TraceBitVector[counter].z - 70;
					bounds[5] = TraceBitVector[counter].z + 70;
					m_trackpoly->GetCellData()->GetScalars("colors")->SetTuple3(TraceBitVector[counter].track_marker,255.0,255.0,255.0);
					m_trackpoly->Modified();
					m_vtkrenderer->ResetCamera(bounds);				
					m_imageview->GetRenderWindow()->Render();

				}
			}
		}
		else
		{
			m_vtkrenderer->RemoveActor(m_trackactor);
			m_vtkrenderer->RemoveActor(currentTrackActor);
			std::vector<TraceLine*>* m_tline_pointer = m_tobj->GetTraceLinesPointer();

			for(int counter=0; counter< m_tline_pointer->size(); counter++)
				{
					if (m_tline_pointer->at(counter)->GetId() == (int) *it)
					{
					//	this->ResetOriginalColors();
						printf("Started creating vtkPolyData for rendering purposes ... ");
						vtkSmartPointer<vtkPolyData> poly_traces=vtkSmartPointer<vtkPolyData>::New();
						vtkSmartPointer<vtkUnsignedCharArray> point_scalars = vtkSmartPointer<vtkUnsignedCharArray>::New();	//responsible for coloring the traces
						point_scalars->SetNumberOfComponents(3);
						point_scalars->SetName("colors");
						vtkSmartPointer<vtkPoints> line_points=vtkSmartPointer<vtkPoints>::New();
						line_points->SetDataTypeToDouble();
						vtkSmartPointer<vtkCellArray> line_cells=vtkSmartPointer<vtkCellArray>::New();
						printf("Starting CreatePolyDataRecursive\n");
						m_tobj->CreatePolyDataRecursive(m_tline_pointer->at(counter),point_scalars,line_points,line_cells);
						printf("Finished CreatePolyDataRecursive\n");
						poly_traces->SetPoints(line_points);
						poly_traces->SetLines(line_cells);	
						poly_traces->GetCellData()->SetScalars(point_scalars);

						vtkSmartPointer<vtkPolyDataMapper> trackmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
						trackmapper->SetInput(poly_traces);
						currentTrackActor->SetMapper(trackmapper);
						currentTrackActor->GetProperty()->SetLineWidth(3.0);
						m_vtkrenderer->AddActor(currentTrackActor);
						m_vtkrenderer->Modified();
						m_imageview->GetRenderWindow()->Render();
					}
				}

		}
	}
		
}
void TrackingKymoView::CreatePointer3D()
{ 
	pointerWidget3d = vtkSmartPointer<vtkPointWidget>::New();
	double sceneBounds[6];
	m_vtkrenderer->ComputeVisiblePropBounds(sceneBounds);
	pointerWidget3d->PlaceWidget(sceneBounds);
	pointerWidget3d->SetInteractor(m_imageview->GetInteractor());
	pointerWidget3d->AllOff();
	double pos[3];
	pos[0]= (sceneBounds[1]-sceneBounds[0])/2.0;
	pos[1]= (sceneBounds[3]-sceneBounds[2])/2.0;
	pos[2]= (sceneBounds[5]-sceneBounds[4])/2.0;

	pointerWidget3d->SetPosition(pos);
	pointerWidget3d->SetEnabled(1);
}

// this function needs to be called before rerendering:
void TrackingKymoView::ResetOriginalColors()
{ 
	vtkSmartPointer<vtkUnsignedCharArray> Original =  vtkSmartPointer<vtkUnsignedCharArray>::New();
	Original->SetNumberOfComponents(3);
	Original->SetName("colors");
	Original->DeepCopy(m_tobj->OriginalColors);
	this->UpdateViewColors(Original);
}
void TrackingKymoView::SetTraceBitClasses(int n)
{
	this->numClasses = n;
	for(int i=0; i<(int)tableVector.size() ; ++i) // iterate over time
	{
		for (int row = 0; (int)row < tableVector.at(i)->GetNumberOfRows(); ++row)
		{
			vtkVariant id = tableVector.at(i)->GetValueByName(row,"ID");
			vtkVariant class_id = tableVector.at(i)->GetValueByName(row,"prediction_active");
			for( int j=0; j<(int)TraceBitVector.size(); ++j)
			{
				if(TraceBitVector[j].time == i && TraceBitVector[j].id == id.ToInt())
				{
					TraceBitVector[j].class_id = class_id.ToInt();
				}
			}
		}
	}
	
}
void TrackingKymoView::SetTrackClasses(int n)
{
	this->numClasses = n;						// number of track classes
	for (int row = 0; (int)row < TFTable->GetNumberOfRows(); ++row)
	{
		vtkVariant id = TFTable->GetValueByName(row,"ID");
		vtkVariant class_id = TFTable->GetValueByName(row,"prediction_active");
		for( int j=0; j<(int)TraceBitVector.size(); ++j)
		{
			if(TraceBitVector[j].id == id.ToInt())
			{
				TraceBitVector[j].class_id = class_id.ToInt();
			}
		}
	}

	
}

std::vector<std::vector<unsigned char> >  TrackingKymoView::GetClassColors()
{
	std::vector<std::vector<unsigned char> > classColors;
	std::vector<unsigned char> classcolor;
	for (int clas = 1; clas<= numClasses; ++clas)
	{
	
		switch(clas)
		{
			//Blue
			case 1:
					classcolor.push_back(0);
					classcolor.push_back(0);
					classcolor.push_back(255);


					
							 break;
			//Red
			case 2:
	    			classcolor.push_back(255);
					classcolor.push_back(0);
					classcolor.push_back(0);					
							 break;
			//Green
			case 3:
	    			classcolor.push_back(0);
					classcolor.push_back(255);
					classcolor.push_back(0);
					
							 break;
			//Cyan
			case 4:
	    			classcolor.push_back(0); 
					classcolor.push_back(255);
					classcolor.push_back(255);
					
							 break;
			
			//Orange 	255-165-0
			case 5:
	    			classcolor.push_back(255);
					classcolor.push_back(165);
					classcolor.push_back(0);
					
							 break;
			//Violet 	238-130-238
			case 6:
	    			classcolor.push_back(238);
					classcolor.push_back(130);
					classcolor.push_back(238);
					 
							 break;
			//Yellow 	255-255-0
			case 7:
	    			classcolor.push_back(255);
					classcolor.push_back(255);
					classcolor.push_back(0);
					
							 break;
			//Dark Green 	0-100-0
			case 8:
	    			classcolor.push_back(0);
					classcolor.push_back(100);
					classcolor.push_back(0);
					
							 break;
			//Royal Blue 	65-105-225
			case 9:
	    			classcolor.push_back(65);
					classcolor.push_back(105);
					classcolor.push_back(225);
					
							 break;
			//Gray 	190-190-190
			case 10: 
	    			classcolor.push_back(190);
					classcolor.push_back(190);
					classcolor.push_back(190);
					
							 break;
	  }
		classColors.push_back(classcolor);
		classcolor.clear();
	}
  return classColors;
}
vtkSmartPointer<vtkUnsignedCharArray> TrackingKymoView::CreateClassColorArray()
{
	std::vector<std::vector<unsigned char> > my_colors = this->GetClassColors();
	vtkSmartPointer<vtkUnsignedCharArray> ClassColorsVtk =  vtkSmartPointer<vtkUnsignedCharArray>::New();
	ClassColorsVtk->SetNumberOfComponents(3);
	ClassColorsVtk->SetName("colors");
//	ClassColorsVtk->SetNumberOfValues((int)TraceBitVector.size());
	


	for (int c=0 ;c < my_colors.size() ; ++c)
	{
		for( int j=0; j<(int)TraceBitVector.size(); ++j)
		{
			if( TraceBitVector[j].class_id == c+1)
			{
				unsigned char color[3];
				color[0] = my_colors.at(c).at(0);
				color[1] = my_colors.at(c).at(1);
				color[2] = my_colors.at(c).at(2);
				ClassColorsVtk->InsertTupleValue(TraceBitVector[j].track_marker, color);
			}
		}
	}
	return  ClassColorsVtk;
}

void TrackingKymoView::UpdateViewColors(vtkSmartPointer<vtkUnsignedCharArray> ColorArray)
{
	m_vtkrenderer->AddActor(m_trackactor);
	m_vtkrenderer->Modified();
	m_trackpoly->GetCellData()->RemoveArray("colors");
	m_trackpoly->GetCellData()->SetScalars(ColorArray);
	m_trackpoly->Modified();
	m_imageview->GetRenderWindow()->Render();

}
void TrackingKymoView::UpdateTraceBitClassesView()
{
	vtkSmartPointer<vtkUnsignedCharArray> ColorArray = vtkSmartPointer<vtkUnsignedCharArray>::New();
	ColorArray = this->CreateClassColorArray();
	this->UpdateViewColors(ColorArray);
	
}
void TrackingKymoView::HandleKeyPress(vtkObject* caller, unsigned long event, void* clientdata, void* callerdata)
{
	TrackingKymoView* view = (TrackingKymoView*)clientdata;
	char key = view->Interactor->GetKeyCode();
	switch (key)
	{
	case 'v':
		view->ToggleSelectionMode();
		break;

	case 'l':
		view->ToggleLabelVisibility();
		break;
	}

}
void TrackingKymoView::ToggleLabelVisibility()
{
	if(labelsVisible) // I only have two actors
	{
		m_vtkrenderer->RemoveActor2D(labelActor);
		labelsVisible = false;
	}
	else
	{
		m_vtkrenderer->AddActor2D(labelActor);
		labelsVisible = true;
	}
	m_vtkrenderer->Modified();
	m_imageview->GetRenderWindow()->Render();

}
void TrackingKymoView::ToggleSelectionMode()
{
	if(singleTracksVisible)
	{
		this->ResetOriginalColors();
		singleTracksVisible = false;
	}
	else
		singleTracksVisible = true;
}

void TrackingKymoView::CreateInteractorStyle()
{
	this->Interactor = this->m_imageview->GetRenderWindow()->GetInteractor();
	//keep mouse command observers, but change the key ones
	this->keyPress = vtkSmartPointer<vtkCallbackCommand>::New();
	this->keyPress->SetCallback(HandleKeyPress);
	this->keyPress->SetClientData(this);

	this->Interactor->RemoveObservers(vtkCommand::KeyPressEvent);
	this->Interactor->RemoveObservers(vtkCommand::KeyReleaseEvent);
	this->Interactor->RemoveObservers(vtkCommand::CharEvent);
	this->Interactor->AddObserver(vtkCommand::KeyPressEvent, this->keyPress);
}