/*=========================================================================
Copyright 2009 Rensselaer Polytechnic Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 
=========================================================================*/

#include "TrackingVolumeView.h"

void TrackingVolumeView::setTime(int t)
{
	printf("Entered setTime\n");
	//m_vtkrenderer->Render();
	m_vtkvolumes[m_currenttime]->VisibilityOff();
	printf("Removed old one\n");
	m_currenttime = t;
	m_currentvol = m_vtkvolumes[m_currenttime];
	m_vtkvolumes[m_currenttime]->VisibilityOn();
	double time_offset  = 0;//(m_model->getRawImagePointer(0)->GetLargestPossibleRegion().GetSize()[2])*3;
	m_vtkvolumes[m_currenttime]->SetPosition(0,0,m_currenttime*time_offset);
	printf("made the new volume visible\n");
	//m_vtkrenderer->Render();
	//m_imageview->markCachedImageAsDirty();
//	m_callback_contrast->volume = m_currentvol;
//	m_callback_brightness->volume = m_currentvol;
	m_imageview->repaint();
	printf("Repainted \n");
	
}
void TrackingVolumeView::GenerateTracks()
{
	printf("Started GenerateTracks\n");
	std::vector<std::vector<ftk::LabelImageFeatures>> *f = m_model->getFeatures();
	int max_track_num = 0;
	for(int cox=0; cox< f->size(); cox++)
	{
		for(int coy =0; coy < (*f)[cox].size(); coy ++)
		{
			max_track_num = MAX((*f)[cox][coy].num,max_track_num);
		}
	}
	printf("Computed max_track_num = %d\n",max_track_num);
	double time_offset  = 0;//(m_model->getRawImagePointer(0)->GetLargestPossibleRegion().GetSize()[2])*3;//*((*f)[0][0].spacing[2]);
	
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
					tbit.x = (*f)[cox][coy].centroid[0]*1.0;
					tbit.y = (*f)[cox][coy].centroid[1]*1.0;
					tbit.z = (*f)[cox][coy].centroid[2]*3.0+(*f)[cox][coy].time*(time_offset);
					tbit.id = (*f)[cox][coy].time;
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



vtkSmartPointer<vtkActor> TrackingVolumeView::getTrackPoints(std::vector<TraceBit> vec)
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

void TrackingVolumeView::AddSliders()
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

  m_callback_brightness = vtkSlider2DCallbackBrightness::New();
  m_callback_brightness->volume = &m_vtkvolumes;
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

  m_callback_contrast = vtkSlider2DCallbackContrast::New();
  m_callback_contrast->volume = &m_vtkvolumes;
  sliderWidget2->AddObserver(vtkCommand::InteractionEvent,m_callback_contrast);
  sliderWidget2->EnabledOn();

}
void TrackingVolumeView::refreshImages(int t)
{

}

void TrackingVolumeView::refreshSelection()
{

}

void TrackingVolumeView::GenerateImages()
{

	float col[3]={1.0,1.0,1.0};
	/*float org[3];
	org[0] = 0; org[1]=0; org[2] = 0;
	m_model->getRawImagePointer(m_currenttime)->SetOrigin(org);;*/
	printf("Creating Volumes ... ");
	for(int counter=0; counter<m_model->getNumTimePoints(); counter++)
	{
		printf(" %d");
		m_volconnects.push_back(ConnectorType::New());
		m_volconnects[counter]->SetInput(m_model->getRawImagePointer(counter));
		m_volconnects[counter]->Update();
		m_vtkims.push_back(vtkSmartPointer<vtkImageData>::New());
		
		m_vtkims[counter]->DeepCopy(m_volconnects[counter]->GetOutput());
		m_vtkims[counter]->SetSpacing(1.0,1.0,3.0);

		//m_vtkims.push_back(m_volconnects[counter]->GetOutput());
		//m_vtkims[counter]->SetSpacing(1.0,1.0,5.0);
		m_vtkvolumes.push_back(getOneVTKVolume(m_vtkims[counter],col));	
		m_vtkrenderer->AddVolume(m_vtkvolumes[counter]);
		m_vtkvolumes[counter]->VisibilityOff();
	}
	printf("\n");
	m_vtkvolumes[m_currenttime]->VisibilityOn();
	m_currentvol = m_vtkvolumes[m_currenttime];
	m_vtkrenderer->ResetCamera(m_vtkims[m_currenttime]->GetBounds());
	
	//begin cheating
	//vtkSmartPointer<vtkCubeSource> sqsource= vtkSmartPointer<vtkCubeSource>::New();
	//sqsource->SetBounds(0,480,0,400,-0.01,0.01);
	//vtkSmartPointer<vtkPolyDataMapper> sqmap = vtkSmartPointer<vtkPolyDataMapper>::New();
	//sqmap->SetInput(sqsource->GetOutput());
	//vtkSmartPointer<vtkActor> sqact = vtkSmartPointer<vtkActor>::New();
	//sqact->SetMapper(sqmap);
	//sqact->SetPickable(1);
	//sqact->GetProperty()->SetOpacity(0.01);
	//m_vtkrenderer->AddActor(sqact);
	//end cheating

	for(int counter=m_model->getNumTimePoints()-1; counter>=0; counter--)
	{
		setTime(counter);
	}
	for(int counter=m_model->getNumTimePoints()-1; counter>=0; counter--)
	{
		setTime(counter);
	}
	//m_vtkrenderer->Render();
	//m_imageview->repaint();
}