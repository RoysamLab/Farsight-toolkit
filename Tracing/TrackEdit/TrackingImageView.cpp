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

#include "TrackingImageView.h"


void TrackingImageView::InitializeSliders()
{
  sliderRep1 = vtkSmartPointer<vtkSliderRepresentation2D>::New();
  sliderRep1->SetValue(0.3);
  sliderRep1->SetTitleText("Brightness");
  sliderRep1->GetPoint1Coordinate()->SetCoordinateSystemToNormalizedDisplay();
  sliderRep1->GetPoint1Coordinate()->SetValue(0.2,0.9);
  sliderRep1->GetPoint2Coordinate()->SetCoordinateSystemToNormalizedDisplay();
  sliderRep1->GetPoint2Coordinate()->SetValue(0.8,0.9);
  sliderRep1->SetSliderLength(0.01);
  sliderRep1->SetSliderWidth(0.03);
  sliderRep1->SetEndCapLength(0.005);
  sliderRep1->SetEndCapWidth(0.03);
  sliderRep1->SetTubeWidth(0.005);
  sliderRep1->SetMinimumValue(0.0);
  sliderRep1->SetMaximumValue(1.0);

  sliderWidget1 = vtkSmartPointer<vtkSliderWidget>::New();
  sliderWidget1->SetInteractor(m_imageview->GetInteractor());
  sliderWidget1->SetRepresentation(sliderRep1);
  sliderWidget1->SetAnimationModeToAnimate();

  m_brcallback = vtkSmartPointer<vtkCallbackCommand>::New();
  m_brcallback->SetCallback(changeBrightness);
  m_brcallback->SetClientData(static_cast<vtkImageShiftScale*>(m_vtkimageshiftscale));
  sliderWidget1->AddObserver(vtkCommand::InteractionEvent,m_brcallback);
  sliderWidget1->EnabledOn();
}
void TrackingImageView::InitializeCallbacks()
{
	m_imcallback = vtkSmartPointer<vtkCallbackCommand>::New();
	m_imcallback->SetCallback(pickPixel);
	m_imcallback->SetClientData(this);
	m_imageview->GetInteractor()->AddObserver(vtkCommand::LeftButtonPressEvent,m_imcallback);

	m_kycallback = vtkSmartPointer<vtkCallbackCommand>::New();
	m_kycallback->SetCallback(keyPress);
	m_kycallback->SetClientData(this);
	m_imageview->GetInteractor()->AddObserver(vtkCommand::KeyReleaseEvent,m_kycallback);


	//begin cheating: cheat vtk by drawing a cube around the region - used to pick the image
	vtkSmartPointer<vtkCubeSource> sqsource= vtkSmartPointer<vtkCubeSource>::New();
	sqsource->SetBounds(0,480,0,400,-0.01,0.01);
	vtkSmartPointer<vtkPolyDataMapper> sqmap = vtkSmartPointer<vtkPolyDataMapper>::New();
	sqmap->SetInput(sqsource->GetOutput());
	vtkSmartPointer<vtkActor> sqact = vtkSmartPointer<vtkActor>::New();
	sqact->SetMapper(sqmap);
	sqact->SetPickable(1);
	sqact->GetProperty()->SetOpacity(0.01);
	m_vtkrenderer->AddActor(sqact);
	//end cheating

	m_picker = vtkSmartPointer<vtkCellPicker>::New();
	m_picker->SetTolerance(0.0001);
	m_imageview->GetInteractor()->SetPicker(m_picker);

	vtkSmartPointer<vtkInteractorStyleImage> isi = vtkSmartPointer<vtkInteractorStyleImage>::New();
	m_imageview->GetInteractor()->SetInteractorStyle(isi);


	InitializeSliders();
}


void TrackingImageView::changeBrightness(vtkObject *caller, unsigned long eid, void* client_data, void* caller_data)
{
	vtkSliderWidget *sliderWidget = reinterpret_cast<vtkSliderWidget*>(caller);
	vtkImageShiftScale * ss = reinterpret_cast<vtkImageShiftScale*>(client_data);
	double value = static_cast<vtkSliderRepresentation *>(sliderWidget->GetRepresentation())->GetValue();
	ss->SetScale(value*10);
}


void TrackingImageView::keyPress(vtkObject * object, unsigned long eid, void *clientdata, void * callerdata)
{
	TrackingImageView * timv = (TrackingImageView*) clientdata;
	int keypressed = timv->m_imageview->GetInteractor()->GetKeyCode();
	DEBUG1("Keyinput = %d\n",timv->m_imageview->GetInteractor()->GetKeyCode());
	DEBUG1("Key repeat count = %d\n",timv->m_imageview->GetInteractor()->GetRepeatCount());

	//'m/M' for merge
	//'d/D' for delete

	switch(keypressed)
	{
	case 'm': case 'M': 
		if(timv->m_select->get(-1).size()>0)//Need more than 1 to merge but renaming can be accomplished by merge
		{
			printf("Enter the trackid to merge the segments to: ");
			int tid;
			scanf("%d",&tid);
			DEBUG2("Merging the selected cells...\n");
			timv->m_model->mergeTracks(timv->m_select,tid);
			timv->m_select->clear();
		}
		break;
	case 'd': case 'D':
		if(timv->m_select->get(-1).size()>0)//Need more than 0 to delete
		{
			printf("Deleting the selected cells...\n");
			timv->m_model->deleteTracks(timv->m_select);
			timv->m_select->clear();
		}
		break;
	case 's': case 'S':
		char buff[1024];
		fgets(buff,1024,stdin);
		printf("Saving all the labeled images back to disk. Are you sure you want to continue (y/n)?\n");
		buff[0]=0;
		fgets(buff,1024,stdin);
		char ch = buff[0];
		if(ch=='y')
		{
			timv->m_model->saveLabels();
		}
		else
		{
			printf("Aborted\n");
		}
		break;


	}
	//timv->setTime(timv->m_currenttime);
}
void TrackingImageView::pickPixel(vtkObject * object, unsigned long  eid, void * clientdata, void * callerdata)
{
	TrackingImageView *timv = (TrackingImageView*)clientdata;
	int * eventpos = timv->m_imageview->GetInteractor()->GetEventPosition();
	int shiftkey = timv->m_imageview->GetInteractor()->GetShiftKey();
	DEBUG2("Shift key = %d\n",shiftkey);
	timv->m_picker->Pick(eventpos[0],eventpos[1],0.0,timv->m_vtkrenderer);
	DEBUG2("Event at %d %d\n",eventpos[0],eventpos[1]);
	double PickPos[3];
	
	if(timv->m_picker->GetViewProp()!=NULL)
	{
		timv->m_picker->GetPickPosition(PickPos);
		DEBUG2("Real World pos = %0.2lf %0.2lf %0.2lf\n",PickPos[0],PickPos[1],PickPos[2]);
		std::vector<int> labels = timv->m_model->GetLabelsAlongZ(PickPos[0],PickPos[1],timv->m_currenttime);
		for(unsigned int counter = 0; counter < labels.size(); counter++)
		{
			if(shiftkey==0)
			{
				if(!timv->m_select->isSelected(labels[counter],timv->m_currenttime))
				{
					timv->m_select->add(labels[counter],timv->m_currenttime);
				}
				else
				{
					timv->m_select->remove(labels[counter],timv->m_currenttime);
				}
			}
			else
			{
				if(!timv->m_select->isSelected(labels[counter],timv->m_currenttime))
				{
					timv->m_select->add(labels[counter],-1);
				}
				else
				{
					timv->m_select->remove(labels[counter],-1);
				}
			}
		}
	//	std::copy(labels.begin(), labels.end(), std::ostream_iterator<int>(std::cout));
	}
}

void TrackingImageView::refreshImages(int t)
{
	DEBUG3("Entered refreshImages(%d)\n",t);
		// update the borders and the features.
		//showborders[t]->Delete();
		Connector2DType::Pointer conn = Connector2DType::New();
		conn->SetInput(get2DMaskedImage(m_model->getRawImagePointer(t),m_model->getLabelImagePointer(t)));
		conn->Update();
		vtkSmartPointer<vtkImageData> vtkim2d = vtkSmartPointer<vtkImageData>::New();
		vtkim2d->DeepCopy(conn->GetOutput());
		showimages[t]=vtkim2d;
		showborders[t] = getActorForPolyData(get2DBoundary(m_model->getLabelImagePointer(t)));
		showborders[t]->GetProperty()->SetColor(1,0,0);
		DEBUG2("GotActorForPolyData\n");
		std::vector<vtkSmartPointer<vtkTextActor>> acts;
		for(unsigned int fc = 0; fc<(*features)[t].size(); fc++)
		{
			acts.push_back(getActorForFeature((*features)[t][fc]));
		}
		if(t==m_currenttime)
		{
			for(unsigned int counter=0; counter<shownumbers[m_currenttime].size(); counter++)
			{
				m_vtkrenderer->RemoveActor(shownumbers[m_currenttime][counter]);
			}
		}
		shownumbers[t] = acts;
		if(t==m_currenttime)
		{
			for(unsigned int counter=0; counter<shownumbers[m_currenttime].size(); counter++)
			{
				m_vtkrenderer->AddActor2D(shownumbers[m_currenttime][counter]);
			}
		}
		setTime(m_currenttime);
		//refreshSelection();
	DEBUG3("Leaving refreshImages(%d)\n",t);
}
void TrackingImageView::GenerateImages()
{

	DEBUG3("Entered GenerateImages\n");
	int time_points = m_model->getNumTimePoints();
	showimages.clear();
	showborders.clear();
	for(int counter=0; counter<time_points; counter++)
	{
		Connector2DType::Pointer conn = Connector2DType::New();
		conn->SetInput(get2DMaskedImage(m_model->getRawImagePointer(counter),m_model->getLabelImagePointer(counter)));
		conn->Update();
		vtkSmartPointer<vtkImageData> vtkim2d = vtkSmartPointer<vtkImageData>::New();
		vtkim2d->DeepCopy(conn->GetOutput());
		//vtkim2d=conn->GetOutput();
		//conn->Print(std::cout);
		showimages.push_back(vtkim2d);
		showborders.push_back(getActorForPolyData(get2DBoundary(m_model->getLabelImagePointer(counter))));
		(showborders.back())->GetProperty()->SetColor(0,1,0);
		double r,g,b;
		showborders.back()->GetProperty()->GetColor(r,g,b);
		DEBUG1("I got color as [%0.1lf %0.1lf %0.1lf]\n",r,g,b);
		std::vector<vtkSmartPointer<vtkTextActor>> acts;
		for(unsigned int fc = 0; fc<(*features)[counter].size(); fc++)
		{
			acts.push_back(getActorForFeature((*features)[counter][fc]));
		}
		shownumbers.push_back(acts);
		if(counter == m_currenttime)
		{
			for(unsigned int fc = 0; fc < (*features)[counter].size(); fc++)
			{
			//	acts[fc]->Print(std::cout);
				m_vtkrenderer->AddActor2D(acts[fc]);
			}
		}
	}

	DEBUG3("Exiting GenerateImages\n");
}

void TrackingImageView::setTime(int t)
{
	DEBUG3("got t = %d\n",t);
	
	
	m_vtkrenderer->RemoveActor(m_vtkcontouractor);
//	m_vtkrenderer->RemoveActor(m_vtkimageactor);
	
	for(unsigned int counter=0; counter<shownumbers[m_currenttime].size(); counter++)
	{
		m_vtkrenderer->RemoveActor(shownumbers[m_currenttime][counter]);
	}
	m_currenttime = t;
	//m_vtkrenderer->RemoveAllViewProps();
	

	//vtkSmartPointer<vtkCubeSource> sqsource= vtkSmartPointer<vtkCubeSource>::New();
	//sqsource->SetBounds(0,480,0,400,-0.01,0.01);
	//vtkSmartPointer<vtkPolyDataMapper> sqmap = vtkSmartPointer<vtkPolyDataMapper>::New();
	//sqmap->SetInput(sqsource->GetOutput());
	//vtkSmartPointer<vtkActor> sqact = vtkSmartPointer<vtkActor>::New();
	//sqact->SetMapper(sqmap);
	//sqact->SetPickable(1);
	//sqact->GetProperty()->SetOpacity(0.01);
	//m_vtkrenderer->AddActor(sqact);

	//m_vtkimageactor->SetInput(showimages[t]);
	//m_vtkrenderer->RemoveActor(m_vtkcontouractor);
	//m_vtkcontouractor->Delete();
	m_vtkcontouractor = showborders[t];
	(showborders[t])->GetProperty()->SetColor(1,0,0);
	m_vtkimageshiftscale->SetInput(showimages[t]);
	m_vtkrenderer->AddActor(m_vtkcontouractor);
//	m_vtkrenderer->AddActor(m_vtkimageactor);
	for(unsigned int counter=0; counter<shownumbers[t].size(); counter++)
	{
		m_vtkrenderer->AddActor2D(shownumbers[t][counter]);
	}

	refreshSelection();
	m_imageview->GetRenderWindow()->Render();

}


void TrackingImageView::refreshSelection()
{
	DEBUG3("Entered refreshSelection\n");
	std::set<int> selections = m_select->get(m_currenttime);
	std::set<int>::iterator iter = selections.begin();


	bool created_flag = false;
	
	m_vtkrenderer->RemoveActor(m_selectionactor);
	if(selections.size()>0)
	{
		

		std::vector<int> box;
		DEBUG2("got the feature vector for %d\n",m_currenttime);
		vtkSmartPointer<vtkAppendPolyData> apppoly = vtkSmartPointer<vtkAppendPolyData>::New();
		while(iter!=selections.end())
		{
			int val = *iter;
		std::vector<ftk::LabelImageFeatures>::iterator fiter = (*features)[m_currenttime].begin();
		std::vector<ftk::LabelImageFeatures>::iterator fend = (*features)[m_currenttime].end();

		while(fiter!=fend)
		{
	/*		while(fiter->num < val)
			{
				++fiter;
				if(fiter == fend)
					break;
			}
			if(fiter== fend)
				break;*/

			if(fiter->num==val)
			{
				box = fiter->boundingbox;
				apppoly->AddInput(getRectangle(box[0]-2,box[2]-2,box[1]+3,box[3]+3));
	//			++fiter;
				created_flag = true;
				break;
			}
			else
			{
	//			++fiter;
			}
			++fiter;
			if(fiter==fend)
				break;
		}
			++iter;


			// get bbox from features[m_curenttime] 
			// generate square polydata of bbox and collect them together
		}

		DEBUG2("Got the poly for all the rectangles\n");
		if(created_flag == true)
		{
			printf("Created_flag = true\n");
			m_selectionactor = vtkSmartPointer<vtkActor>::New();
			vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			mapper->SetInput(apppoly->GetOutput());
			m_selectionactor->SetMapper(mapper);
			m_selectionactor->GetProperty()->SetColor(1,1,0);
			m_vtkrenderer->AddActor(m_selectionactor);
			m_imageview->repaint();
		}

	}
	

	// create actor for that and assign m_selectionactor to that
	// add m_selectionactor to m_vtkrenderer
}