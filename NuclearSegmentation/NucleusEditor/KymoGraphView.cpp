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
	std::vector<std::vector<ftk::IntrinsicFeatures>> *f = &myfeatures;
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
	vtkSmartPointer<vtkPoints> points=vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> cells=vtkSmartPointer<vtkCellArray>::New();
	for(int counter=0; counter<vec.size(); counter++)
	{
		int return_id = points->InsertNextPoint(vec[counter].x,vec[counter].y,vec[counter].z);
		cells->InsertNextCell(1);
		cells->InsertCellPoint(return_id);
	}
	printf("About to create poly\n");
	point_poly->SetPoints(points);
	point_poly->SetVerts(cells);

	vtkSmartPointer<vtkPolyDataMapper> cubemap = vtkSmartPointer<vtkPolyDataMapper>::New();
	cubemap->SetInput(point_poly);
	cubemap->GlobalImmediateModeRenderingOn();
	vtkSmartPointer<vtkActor> cubeact = vtkSmartPointer<vtkActor>::New();
	cubeact->SetMapper(cubemap);
	cubeact->SetPickable(0);
	cubeact->GetProperty()->SetPointSize(5);
	cubeact->GetProperty()->SetOpacity(.5);
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
}