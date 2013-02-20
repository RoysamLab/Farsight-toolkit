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
	//m_vtkvolume->Print(std::cout);
	m_vtkrenderer->AddVolume(m_vtkvolume);
	m_vtkrenderer->ResetCamera(m_vtkim->GetBounds());

}

void TrackingKymoView::SaveMovie(void)
{
/*	std::vector<TraceLine*>* m_tline_pointer = m_tobj->GetTraceLinesPointer();

	float col[3]={1.0,1.0,1.0};
	vtkSmartPointer<vtkImageData> vtkimdata = vtkSmartPointer<vtkImageData>::New();
//	InputImageType::SizeType imsize = m_model->getRawImagePointer(0)->GetLargestPossibleRegion().GetSize();
	unsigned short xsize = (my4DImg->GetImageInfo()->numColumns);
	unsigned short ysize = (my4DImg->GetImageInfo()->numRows);
	unsigned short tsize = (my4DImg->GetImageInfo()->numTSlices);

	vtkimdata->SetExtent(0,xsize-1,0,ysize-1,0,tsize-1);
	vtkimdata->SetScalarTypeToUnsignedChar();
	vtkimdata->AllocateScalars();
	vtkimdata->SetSpacing(1.0,1.0,1.0);

	//Loop Constants:
	unsigned short xy = xsize*ysize;
	ftk::Image::PtrMode mode;
	mode = static_cast<ftk::Image::PtrMode>(0); 
	Input2DImageType::Pointer MaxIntProjImage = Input2DImageType::New();
	Input2DImageType::RegionType region;
	Input2DImageType::SizeType size;
	Input2DImageType::IndexType index;
	index[0]=0;
	index[1]=0;
	size[0] = xsize;
	size[1] = ysize;
	region.SetSize(size);
	region.SetIndex(index);
	MaxIntProjImage->SetRegions(region);
	MaxIntProjImage->Allocate();



	for(int t=0; t< (int)tsize; t++)
	{
		// Get the 2-D max intensity projection
		SliceIteratorType inputIt(my4DImg->GetItkPtr<InputPixelType>((unsigned short)t,0,mode),my4DImg->GetItkPtr<InputPixelType>(t,0,mode)->GetLargestPossibleRegion());
		LinearIteratorType outputIt(MaxIntProjImage,MaxIntProjImage->GetLargestPossibleRegion());
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

		RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
		rescaleFilter->SetInput(MaxIntProjImage);
		rescaleFilter->SetOutputMinimum(0);
		rescaleFilter->SetOutputMaximum(255);
		rescaleFilter->Update();

		RGB2DFilterType::Pointer rgbfilter = RGB2DFilterType::New();
		rgbfilter->SetInput(rescaleFilter->GetOutput());
		rgbfilter->SetColormap( RGB2DFilterType::Hot );
		rgbfilter->Update();
		//Connector2DType::Pointer connector2d = Connector2DType::New();
		//connector2d->SetInput(MaxIntProjImage);
		//connector2d->Update();
		


		RGBConnector2DType::Pointer connector2d = RGBConnector2DType::New();
		connector2d->SetInput(rgbfilter->GetOutput());
		connector2d->Update();

		// Get the the text image (Number Image)for this time slice:
		vtkSmartPointer<vtkImageData> LabelImage = vtkSmartPointer<vtkImageData>::New();
		LabelImage->SetExtent(0,xsize-1,0,ysize-1,0,0);
		LabelImage->SetNumberOfScalarComponents(3);
		LabelImage->SetScalarTypeToUnsignedChar();
		LabelImage->AllocateScalars();


		// Create an image of text

		for(int i=0;i<(int)myfeatures[t].size();++i)
		 {
			 stringstream ss;
			 ss<< myfeatures[t][i].num;
			 std::cout<< i<<std::endl;

			vtkSmartPointer<vtkFreeTypeUtilities> freeType = vtkSmartPointer<vtkFreeTypeUtilities>::New();
			vtkSmartPointer<vtkTextProperty> textProperty = vtkSmartPointer<vtkTextProperty>::New();
			textProperty->SetColor( 1.0,0.0,0.0 );
			textProperty->SetFontSize(10);

			 vtkSmartPointer<vtkImageData> textImage = vtkSmartPointer<vtkImageData>::New();
			 textImage->SetNumberOfScalarComponents(3);
			 freeType->RenderString(textProperty,ss.str().c_str(), textImage);
			 std::cout<<ss.str()<<std::endl;

			 //vtkSmartPointer<vtkImageBlend> blend = vtkSmartPointer<vtkImageBlend>::New();
			 //blend->AddInputConnection(drawing->GetOutputPort());
			 //blend->AddInputConnection(textImage->GetProducerPort());
			 //blend->SetOpacity(0,.0);
			 //blend->SetOpacity(1,1);
			 //blend->Update();
			 float bbox[4];
			 bbox[0] = myfeatures[t][i].BoundingBox[0];
			 bbox[1] = myfeatures[t][i].BoundingBox[1];
			 bbox[2] = myfeatures[t][i].BoundingBox[2];
			 bbox[3] = myfeatures[t][i].BoundingBox[3];

			 this->AddLabelToVTKImage(LabelImage,textImage,bbox);
		}

		vtkSmartPointer<vtkImageBlend> newblend = vtkSmartPointer<vtkImageBlend>::New();
		newblend->AddInput(connector2d->GetOutput());
		newblend->AddInputConnection(LabelImage->GetProducerPort());
	    newblend->SetOpacity(0,.6);
		newblend->SetOpacity(1,1);
		newblend->Update();

		 std::cout<<t<<std::endl;
//		 vtkSmartPointer<vtkTIFFWriter> tiffWriter = vtkSmartPointer<vtkTIFFWriter>::New();
		 vtkSmartPointer<vtkPNGWriter> tiffWriter = vtkSmartPointer<vtkPNGWriter>::New();
		 stringstream ss;
		 ss<<t;
//		 std::string filename = "C:\\Users\\amerouan\\Desktop\\FeaturesTests\\movie_t"+ss.str()+".tif";
		 std::string filename = "C:\\Users\\amerouan\\Desktop\\FeaturesTests\\movie_t"+ss.str()+".png";
		 tiffWriter->SetFileName(filename.c_str());
	//	 tiffWriter->SetInput(connector2d->GetOutput());
		 tiffWriter->SetInput(newblend->GetOutput());
		 tiffWriter->Write();
	}	
*/
}
void TrackingKymoView::AddLabelToVTKImage(vtkSmartPointer<vtkImageData> labelImage, vtkSmartPointer<vtkImageData> newlabelImage,float bbox[])
{
	int* labdim = newlabelImage->GetDimensions();
	std::cout << "Dims: " << " x: " << labdim[0] << " y: " << labdim[1] << " z: " << labdim[2] << std::endl;
	int* dims = labelImage->GetDimensions();
	std::cout << "Dims: " << " x: " << dims[0] << " y: " << dims[1] << " z: " << dims[2] << std::endl;
	// Get the iteration range :
	int xbegin = (int)bbox[0];
	int xend = xbegin+(int)labdim[0]-1;
	int ybegin = (int)bbox[2];
	int yend = ybegin+(int)labdim[1]-1;
	
	// might crash for large size of labdim or small size of labelImage(needs more checking)
	if(xend>=dims[0])
	{
		xbegin = dims[0]-labdim[0]-1;
		xend = dims[0] - 1;
	}
	if(yend>=dims[1])
	{
		ybegin = dims[1]-labdim[1]-1;
		yend = dims[1] - 1;
	}

	int xlab = 0;
	int ylab = 0;
	for (int y = ybegin; y < yend; y++)
	{
		for (int x = xbegin; x < xend; x++)
		{
			unsigned char* pixel1 = static_cast<unsigned char*>(labelImage->GetScalarPointer(x,y,0));
			unsigned char* pixel2 = static_cast<unsigned char*>(newlabelImage->GetScalarPointer(xlab,ylab,0));
			pixel1[0] = pixel2[0];
			pixel1[1] = pixel2[1];
			pixel1[2] = pixel2[2];
			++xlab;
		}
		++ylab;
		xlab = 0;
	}


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
	point_poly = vtkSmartPointer<vtkPolyData>::New();
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
		std::stringstream strId;
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
	cubeact = vtkSmartPointer<vtkActor>::New();
	cubeact->SetMapper(cubemap);
	//cubeact->SetPickable(0);
	cubeact->PickableOff();
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
	case 'd':
		view->Delete();
		break;
	case 'p':
		view->TogglePickMode();
		break;
	case 'c':
		view->ConnectNodes();
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

	this->SetupInteractorStyle();
	this->myCellPicker = vtkSmartPointer<vtkCellPicker>::New();
	this->myCellPicker->SetTolerance(0.004);
	this->Interactor->SetPicker(this->myCellPicker);
	this->isPicked = vtkSmartPointer<vtkCallbackCommand>::New();
	this->isPicked->SetCallback(PickCell);

	//isPicked caller allows observer to intepret click 
	this->isPicked->SetClientData(this);            
	this->Interactor->AddObserver(vtkCommand::RightButtonPressEvent,isPicked);

}
void TrackingKymoView::PickCell(vtkObject* caller, unsigned long event, void* clientdata, void* callerdata)
{
	TrackingKymoView* view = (TrackingKymoView*)clientdata;
	int *pos = view->Interactor->GetEventPosition();
	view->Interactor->GetPicker()->Pick(pos[0],pos[1],0.0,view->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer());
	vtkCellPicker *cell_picker = (vtkCellPicker *)view->Interactor->GetPicker();
	if(!view->Interactor->GetControlKey())
	{
		view->TSelection.clear();
		view->NSelection.clear();
		view->ResetOriginalColors();
	}
	if(view->Interactor->GetControlKey()&& view->NSelection.size()== 2) // make sure that not more than 2 nodes are selected at once
	{
		view->NSelection.clear();
		view->RemoveSphereActors();
	}

	if (cell_picker->GetCellId() == -1) 
			std::cout<<"nothing has been selected\n";
	else if(cell_picker->GetViewProp()!=NULL)
	{
		
		std::cout<<"selected:"<<cell_picker->GetCellId()<<"\n";
		if(view->cubeact->GetPickable()==1)
		{
			view->NSelection.insert(cell_picker->GetCellId());
			double pick_pos[3];
			cell_picker->GetPickPosition(pick_pos);
			view->DrawSphere(pick_pos);
			TraceBit tbit = (reinterpret_cast<TraceLine*>(view->m_tobj->hashp[cell_picker->GetCellId()]))->points_hash[cell_picker->GetCellId()];
			view->Selection->SetCurrentTime(tbit.time);
			view->Selection->select(tbit.id);
		}
		else
		{
			view->TSelection.insert(cell_picker->GetCellId());
			view->m_trackpoly->GetCellData()->GetScalars("colors")->SetTuple3(cell_picker->GetCellId(),255.0,255.0,255.0);
			view->m_trackpoly->Modified();
			view->m_imageview->GetRenderWindow()->Render();
		}
	}
	
}
void TrackingKymoView::RemoveSphereActors(void)
{
	m_vtkrenderer->RemoveActor(sphereActor1);
	m_vtkrenderer->RemoveActor(sphereActor2);
	m_vtkrenderer->Modified();
	m_imageview->GetRenderWindow()->Render();
}

void TrackingKymoView::CreateSphereActors(void)
{
	vtkSmartPointer<vtkSphereSource> sphereSource1 = vtkSmartPointer<vtkSphereSource>::New();
	//sphereSource1->SetCenter(0.0, 0.0, 0.0);
	sphereSource1->SetRadius(0.6);

	vtkSmartPointer<vtkSphereSource> sphereSource2 = vtkSmartPointer<vtkSphereSource>::New();
	//sphereSource2->SetCenter(0.0, 0.0, 0.0);
	sphereSource2->SetRadius(0.6);

	vtkSmartPointer<vtkPolyDataMapper> mapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper1->SetInputConnection(sphereSource1->GetOutputPort());
	vtkSmartPointer<vtkPolyDataMapper> mapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper2->SetInputConnection(sphereSource2->GetOutputPort());

	sphereActor1 = vtkSmartPointer<vtkActor>::New();
	sphereActor2 = vtkSmartPointer<vtkActor>::New();
	sphereActor1->SetMapper(mapper1);
	sphereActor2->SetMapper(mapper2);
	sphereActor1->GetProperty()->SetColor(1,1,1);
	sphereActor2->GetProperty()->SetColor(1,1,1);
}

void TrackingKymoView::DrawSphere(double pos[])
{
	if(NSelection.empty()) return;
	// get point position
	if(NSelection.size() == 1)
	{
		sphereActor1->SetPosition(pos);
		m_vtkrenderer->AddActor(sphereActor1);
	}
	else if (NSelection.size() == 2)
	{
		sphereActor2->SetPosition(pos);
		m_vtkrenderer->AddActor(sphereActor2);
	}
	m_vtkrenderer->Modified();
	m_imageview->GetRenderWindow()->Render();

}
void TrackingKymoView::SetupInteractorStyle(void)
{
	vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
	this->Interactor->SetInteractorStyle(style);
	this->m_imageview->GetRenderWindow()->Render();
}
void TrackingKymoView::Delete(void)
{
	if(TSelection.empty()) return;
	this->DeleteAndRelabelData();
	this->UpdateTracePolyData();
	this->UpdateLabels();
	this->m_imageview->GetRenderWindow()->Render();
	TSelection.clear();		// clear selection
}
void TrackingKymoView::UpdateTracePolyData()
{
	//printf("I am in  UpdateTracePolyData\n");
	m_trackpoly = m_tobj->GetVTKPolyData();
	m_trackmapper->SetInput(m_trackpoly);
	m_trackactor->SetMapper(m_trackmapper);
	m_trackactor->GetProperty()->SetLineWidth(3.0);
	m_trackpoly->Modified();
}

void TrackingKymoView::ConnectNodes(void)
{
	if(NSelection.empty()|| NSelection.size()!=2) return;
	std::vector<TraceLine*> selected_tline;
	std::set<unsigned int>::iterator node_vtk_ids_iter;
	// store the selected trace lines:
	for(node_vtk_ids_iter = NSelection.begin(); node_vtk_ids_iter!= NSelection.end(); ++node_vtk_ids_iter)
	{
		TraceLine *tline = reinterpret_cast<TraceLine*>(m_tobj->hashp[*node_vtk_ids_iter]);
		selected_tline.push_back(tline);
	}
	// Check if Nodes Can be Connected:
	if(!IsValidNodeConnection())
	{
		printf("-----------------------------------------------------------------------\n ");
		printf("Node connection not allowed for this case, try to delete tracks first..\n ");
		printf("-----------------------------------------------------------------------\n ");
		NSelection.clear();		// clear selection
	    this->RemoveSphereActors(); // clear selected nodes from the view
		return;
	}
	// find the trace lines that need to be connected and delete them from the trace line data:
	std::vector<TraceLine*> * trace_lines = m_tobj->GetTraceLinesPointer();
	std::vector<TraceLine*>::iterator trace_lines_iter;
	for(int i=0; i<selected_tline.size(); ++i)
	{
		for(trace_lines_iter = trace_lines->begin() ;trace_lines_iter!= trace_lines->end(); ++trace_lines_iter)
		{
			if(selected_tline.at(i)->GetId() == (*trace_lines_iter)->GetId())
			{
				trace_lines->erase(trace_lines_iter);
				break;
			}
		}
	}
	// make a new trace line by connecting the other two trace lines:
	ObjectSelection::Point point;
	std::vector<ObjectSelection::Point> selected_points;
	int new_id = this->GetMaxId()+1;
	TraceLine * new_tline = new TraceLine();
	for(int i = 0; i<selected_tline.size();++i)
	{
		TraceLine::TraceBitsType::iterator tbit_iter = selected_tline.at(i)->GetTraceBitIteratorBegin();
		while(tbit_iter!=selected_tline.at(i)->GetTraceBitIteratorEnd())
		{
			point.id = tbit_iter->id;
			point.time = tbit_iter->time;
			point.new_id = new_id;
			selected_points.push_back(point);
			tbit_iter->id = new_id;			
			new_tline->AddTraceBit(*tbit_iter);
			++tbit_iter;
		}
	}
	new_tline->SetId(new_id);
	trace_lines->push_back(new_tline);
	this->UpdateTracePolyData();
	this->UpdateLabels();
	this->m_imageview->GetRenderWindow()->Render();
	NSelection.clear();		// clear selection
	this->RemoveSphereActors(); // clear selected nodes from the view
	Selection->SelectPoints(selected_points);

}	

bool TrackingKymoView::IsValidNodeConnection(void)
{
	if(NSelection.empty()|| NSelection.size()!=2) return false;
	std::set<unsigned int>::iterator node_vtk_ids_iter;
	std::vector<unsigned int> node_ids; 
	for(node_vtk_ids_iter = NSelection.begin(); node_vtk_ids_iter!= NSelection.end(); ++node_vtk_ids_iter)
	{
			TraceLine *tline = reinterpret_cast<TraceLine*>(m_tobj->hashp[*node_vtk_ids_iter]);
			int tbit_time = tline->points_hash[*node_vtk_ids_iter].time;
			TraceLine::TraceBitsType::iterator tbit_iter = tline->GetTraceBitIteratorBegin();
			bool endNode = true;
			bool beginNode = true;
			std::vector<bool> tmp_vector;
			//bool singleNode = false;
			while(tbit_iter!=tline->GetTraceBitIteratorEnd())
			{
				if(tbit_iter->time > tbit_time)
					endNode = false;
				else if (tbit_iter->time < tbit_time)
					beginNode = false;
				++tbit_iter;
			}
			 tline->points_hash[*node_vtk_ids_iter].end = endNode;
			 tline->points_hash[*node_vtk_ids_iter].begin = beginNode;
			 node_ids.push_back(*node_vtk_ids_iter);
	}
	// now compare the two nodes:
	if(node_ids.size()!=2)
	{
		printf("something is up with the selection\n");
		return false;
	}
	TraceBit tbit1 = (reinterpret_cast<TraceLine*>(m_tobj->hashp[node_ids.at(0)]))->points_hash[node_ids.at(0)];
	TraceBit tbit2 = (reinterpret_cast<TraceLine*>(m_tobj->hashp[node_ids.at(1)]))->points_hash[node_ids.at(1)];
	bool isvalid = false;
	// one is at end and the other is at begin:
	if(tbit1.end == true && tbit1.begin == false && tbit2.end == false && tbit2.begin == true)
	{
		int timediff = (tbit2.time - tbit1.time);
		if( timediff >0 && timediff< 3)
			isvalid = true;
	}
	else if (tbit1.end == false && tbit1.begin == true && tbit2.end == true && tbit2.begin == false)
	{
		int timediff = (tbit1.time - tbit2.time);
		if( timediff >0 && timediff< 4)
			isvalid = true;
	}
	// one is at end or begin and the other is a single node:
	else if (tbit1.end == true && tbit1.begin == true && tbit2.end == true && tbit2.begin == false)
	{
		int timediff = (tbit1.time - tbit2.time);
		if( timediff >0 && timediff< 3)
			isvalid = true;
	}
	else if (tbit1.end == true && tbit1.begin == true && tbit2.end == false && tbit2.begin == true)
	{
		int timediff = (tbit2.time - tbit1.time);
		if( timediff >0 && timediff< 3)
			isvalid = true;
	}

	else if (tbit1.end == false && tbit1.begin == true && tbit2.end == true && tbit2.begin == true)
	{
		int timediff = (tbit1.time - tbit2.time);
		if( timediff >0 && timediff< 3)
			isvalid = true;
	}
	else if (tbit1.end == true && tbit1.begin == false && tbit2.end == true && tbit2.begin == true)
	{
		int timediff = (tbit2.time - tbit1.time);
		if( timediff >0 && timediff< 3)
			isvalid = true;
	}
	return isvalid;
}
void TrackingKymoView::TogglePickMode(void)
{
	if(cubeact->GetPickable()==1)
	{
		cubeact->PickableOff();
		m_trackactor->PickableOn();
	}
	else
	{
		m_trackactor->PickableOff();
		cubeact->PickableOn();
	}
}

void TrackingKymoView::DeleteAndRelabelData(void)
{
	if(TSelection.empty()) return;
	
	std::vector<unsigned int> tmp_selection;
	std::set<unsigned int> samet_selection;	// list of sub-tracks selected within the same track
	//std::set<unsigned int> difft_selection;
	std::map<int, std::set<unsigned int> > selection_map;

	std::set<unsigned int>::iterator sel_iter;

	for(sel_iter = TSelection.begin(); sel_iter!= TSelection.end(); ++sel_iter)
		tmp_selection.push_back(*sel_iter);

	// figure out the number of different selected tracks:
	std::set<int> selected_ids;
	for(int i = 0; i< tmp_selection.size(); ++i)
	{
		TraceLine *tline1 = reinterpret_cast<TraceLine*>(m_tobj->hashc[tmp_selection.at(i)]);
		selected_ids.insert(tline1->GetId());
	}
	std::set<int>::iterator sel_ids_iter;
	for(sel_ids_iter = selected_ids.begin(); sel_ids_iter!= selected_ids.end(); ++sel_ids_iter)
	{
		std::set<unsigned int> samet_selection;
		for(int i = 0; i< tmp_selection.size(); ++i)
		{
			TraceLine *tline = reinterpret_cast<TraceLine*>(m_tobj->hashc[tmp_selection.at(i)]);
			int id = tline->GetId();
			if(tline->GetId() == *sel_ids_iter)
				samet_selection.insert(tmp_selection.at(i));
		}
		selection_map[*sel_ids_iter] = samet_selection;
	}
	std::map<int, std::set<unsigned int> >::iterator select_map_iter;
	for(select_map_iter = selection_map.begin(); select_map_iter!= selection_map.end(); ++select_map_iter)
	{
		if(!(*select_map_iter).second.empty()&& (*select_map_iter).second.size()>1)
			this->DeleteSameTData((*select_map_iter).second);
		else if(!(*select_map_iter).second.empty()&& (*select_map_iter).second.size()==1)
			this->DeleteDiffTData((*select_map_iter).second);
	}
	if(!points_from_delete.empty())
	{
		Selection->SelectPoints(points_from_delete);
		points_from_delete.clear();
	}

	
}
void  TrackingKymoView::DeleteSameTData(std::set<unsigned int> samet_selection)
{
	std::set<unsigned int>::iterator sel_iter = samet_selection.begin();
	TraceLine *tline = reinterpret_cast<TraceLine*>(m_tobj->hashc[*sel_iter]);
	int new_id = this->GetMaxId()+1;

	// first find the trace line in the data to erase it:
	std::vector<TraceLine*>* trace_lines  = m_tobj->GetTraceLinesPointer();
	std::vector<TraceLine*>::iterator trace_lines_iter;
	for(trace_lines_iter = trace_lines->begin() ;trace_lines_iter!= trace_lines->end(); ++trace_lines_iter)
	{
		if(tline->GetId() == (*trace_lines_iter)->GetId())
		{
			trace_lines->erase(trace_lines_iter);
			break;
		}
	}

	printf("I am here\n");

	TraceLine * tline_tmp = new TraceLine();
	tline_tmp = tline;
	sel_iter = samet_selection.begin();
	while(sel_iter != samet_selection.end())
	{
		TraceBit tb1 = tline->subtrace_hash[*sel_iter].at(0);
		TraceBit tb2 = tline->subtrace_hash[*sel_iter].at(1);
		int maxtime = MAX(tb1.time,tb2.time);
		tline_tmp = this->DeleteTlineRecursive(tline_tmp,*sel_iter,new_id,maxtime);
		++new_id;
		++sel_iter;
	}
	trace_lines->push_back(tline_tmp);

}
TraceLine * TrackingKymoView::DeleteTlineRecursive(TraceLine * tline, unsigned int vtk_cell_id, int new_id,int maxtime)
{
		std::vector<TraceLine*>* trace_lines  = m_tobj->GetTraceLinesPointer();
		TraceLine::TraceBitsType::iterator tbit_iter = tline->GetTraceBitIteratorBegin();
		TraceLine * new_tline1 = new TraceLine();
		TraceLine * new_tline2 = new TraceLine();
		ObjectSelection::Point point;

		while(tbit_iter != tline->GetTraceBitIteratorEnd())
		{
			if((*tbit_iter).time < maxtime)
			{
				new_tline1->AddTraceBit(*tbit_iter);
			}
			else 
			{
				point.id = tbit_iter->id;
				point.time = tbit_iter->time;
				(*tbit_iter).id = new_id;
				point.new_id = new_id;
				new_tline2->AddTraceBit(*tbit_iter);
				points_from_delete.push_back(point);
			}
			++tbit_iter;
		}
		trace_lines->push_back(new_tline1);
		new_tline2->SetId(new_id);
	return new_tline2;
}

void  TrackingKymoView::DeleteDiffTData(std::set<unsigned int> difft_selection)
{
	std::set<unsigned int>::iterator sel_iter;
	for(sel_iter = difft_selection.begin(); sel_iter!= difft_selection.end(); ++sel_iter)
	{
		unsigned int vtk_cell_id = *sel_iter;
		
		TraceLine *tline = reinterpret_cast<TraceLine*>(m_tobj->hashc[vtk_cell_id]);
		std::cout<<"deleting cell:"<<vtk_cell_id<<" with id: "<<tline->GetId()<<std::endl;
		TraceBit tb1 = tline->subtrace_hash[vtk_cell_id].at(0);
		TraceBit tb2 = tline->subtrace_hash[vtk_cell_id].at(1);
		int maxtime = MAX(tb1.time,tb2.time);

		std::vector<TraceLine*>* trace_lines  = m_tobj->GetTraceLinesPointer();
		std::vector<TraceLine*>::iterator trace_lines_iter;
		for(trace_lines_iter = trace_lines->begin() ;trace_lines_iter!= trace_lines->end(); ++trace_lines_iter)
		{
			if(tline->GetId() == (*trace_lines_iter)->GetId())
			{
				int new_id = this->GetMaxId()+1;
				trace_lines->erase(trace_lines_iter);
				TraceLine::TraceBitsType::iterator tbit_iter = tline->GetTraceBitIteratorBegin();
				TraceLine * new_tline1 = new TraceLine();
				TraceLine * new_tline2 = new TraceLine();
				ObjectSelection::Point point;
				
				while(tbit_iter != tline->GetTraceBitIteratorEnd())
				{
					if((*tbit_iter).time < maxtime)
					{
						new_tline1->AddTraceBit(*tbit_iter);
					}
					else 
					{
						point.id = tbit_iter->id;
						point.time = tbit_iter->time;
						(*tbit_iter).id = new_id;
						point.new_id = new_id;
						points_from_delete.push_back(point);
						new_tline2->AddTraceBit(*tbit_iter);
					}
					++tbit_iter;
				}
				trace_lines->push_back(new_tline1);
				new_tline2->SetId(new_id);
				trace_lines->push_back(new_tline2);
				break;
			}
		}
	}

}
int TrackingKymoView::GetMaxId(void)
{
	int max = 0;

	std::vector<TraceLine*>* trace_lines  = m_tobj->GetTraceLinesPointer();
		//std::cout<<"tline size at max id: "<<trace_lines->size()<<"\n";

	for(int i=0; i<trace_lines->size();++i)
		max = MAX(trace_lines->at(i)->GetId(),max);
	return max;
}

void TrackingKymoView::UpdateLabels(void)
{
	//printf("entered UpdateLabels\n");
	m_vtkrenderer->RemoveActor(labelActor);
	m_vtkrenderer->Modified();
	point_poly->GetPointData()->RemoveArray("labels");
	point_poly->GetPointData()->RemoveArray("sizes");
	point_poly->Modified();


	std::vector<TraceLine*>* trace_lines  = m_tobj->GetTraceLinesPointer();
	//std::cout<<"tline size at UpdateLabels: "<<trace_lines->size()<<"\n";

	std::vector<TraceLine*>::iterator trace_lines_iter;

	//vtkSmartPointer<vtkPolyData> new_point_poly = vtkSmartPointer<vtkPolyData>::New();

	vtkSmartPointer<vtkStringArray> labels = vtkSmartPointer<vtkStringArray>::New();	// Labels Array
	//labels->SetNumberOfValues(TraceBitVector.size());
	labels->SetName("labels");
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkIntArray> sizes = vtkSmartPointer<vtkIntArray>::New();			// Priority Arrays used by vtkPointSetToLabelHierarchy
	//sizes->SetNumberOfValues(vec.size());
	sizes->SetName("sizes");

	int id = this->GetMaxId();
	//printf("max id: %d\n",id);
	for(trace_lines_iter = trace_lines->begin() ;trace_lines_iter!= trace_lines->end(); ++trace_lines_iter)
	{
		TraceLine::TraceBitsType::iterator tbit_iter = (*trace_lines_iter)->GetTraceBitIteratorBegin();
		//	std::cout<<"I am in trace: "<<(*trace_lines_iter)->GetId()<<"\n";
		while(tbit_iter != (*trace_lines_iter)->GetTraceBitIteratorEnd())
		{
			if(tbit_iter->id == id)
			{
			//	std::cout<<"I am in trace bit: "<<tbit_iter->id<<", "<< tbit_iter->time<<"\n";
			//	scanf("%d");
			}
			int return_id = points->InsertNextPoint(tbit_iter->x,tbit_iter->y,tbit_iter->z);
			//Labels:
			std::stringstream strId;
			strId <<"("<<tbit_iter->id<<","<<tbit_iter->time<<")";
			labels->InsertValue(return_id, strId.str());
			sizes->InsertValue(return_id, tbit_iter->id);
			++tbit_iter;
		}
	}
	//printf("finished part of UpdateLabels\n");
	point_poly->SetPoints(points);
	point_poly->GetPointData()->AddArray(labels);
	point_poly->GetPointData()->AddArray(sizes);
	point_poly->Modified();
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
	m_vtkrenderer->AddActor(labelActor);
	m_vtkrenderer->Modified();
}


		
			


