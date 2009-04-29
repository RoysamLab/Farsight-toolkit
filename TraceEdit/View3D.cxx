/*
  Farsight ToolKit 3D Viewer: 
  v1: implements render window functionality
    render window creation in "View3d::RenderWin"
    line objects mapped and actors created in "View3d::LineAct"
    actor added and window created in "View3d::addAct"
  v2: add picking functionality:
    access to observing a pick event is in "View3d::interact"
    the interacter calls for "View3d::PickPoint"
  v3: include contourFilter and rayCast renderers
*/
#include "View3D.h"

View3d::View3d()
{
  this->lineAct = vtkActor::New();
  this->lineMap = vtkPolyDataMapper::New();
  this->ren = vtkRenderer::New();
  this->renWin = vtkRenderWindow::New();
  this->iren = vtkRenderWindowInteractor::New();
  this->tobj = new TraceObject;
	gapTol = 2;
	gapMax = 30;
	smallLine = 5;
}

View3d::~View3d()
{
  this->lineAct->Delete();
  this->lineMap->Delete();
  this->ren->Delete();
  this->renWin->Delete();
  this->iren->Delete();
  delete this->tobj;
}

void View3d::Initialize(int argc, char **argv)
{
	int num_loaded = 0;
  // load as many files as possible. Provide offset for differentiating types
	for(int counter=1; counter<argc; counter++)
	{
		int len = strlen(argv[counter]);
		if(strcmp(argv[counter]+len-3,"swc")==0)
		{
			printf("I detected swc\n");
			this->tobj->ReadFromSWCFile(argv[counter]);
		}
		else if (strcmp(argv[counter]+len-3,"xml")==0)
		{
			printf("I detected xml\n");
			this->tobj->ReadFromRPIXMLFile(argv[counter]);
		}
		else if( strcmp(argv[counter]+len-3,"tks")==0)
		{
			printf("I detected tks\n");
			this->tobj->ReadFromFeatureTracksFile(argv[counter],num_loaded);
		}
		else if( strcmp(argv[counter]+len-3,"tif")==0 ||
             strcmp(argv[counter]+len-3, "pic")==0)
		{
			printf("I detected a tif file\n");
			this->rayCast(argv[counter]);
			this->AddVolumeSliders();
		}
		num_loaded++;
	}
}

/*  set up the standard render window */
void View3d::RenderWin()
{
  this->ren->SetBackground(0,0,0);
  this->renWin->AddRenderer(ren);
  this->renWin->SetSize(640, 480);          //window sixe can be dragged to fit, good start
  this->renWin->SetInteractor(iren);
  vtkInteractorStyleTrackballCamera* style = vtkInteractorStyleTrackballCamera::New();
  this->iren->SetInteractorStyle( style );        //trackball control
}

/*  update trace data and return a vtkactor  */
vtkActor* View3d::LineAct()
{
  this->poly_line_data = this->tobj->GetVTKPolyData();
 // this->poly_line_data = this->tobj->GetVTKPolyData();
  this->lineMap->SetInput(this->poly_line_data);
  this->lineAct->SetMapper(lineMap);
  this->lineAct->GetProperty()->SetColor(0,1,0);
  printf("Point size %f\n",lineAct->GetProperty()->GetPointSize());
  this->lineAct->GetProperty()->SetPointSize(2);
  printf("later: Point size %f\n",lineAct->GetProperty()->GetPointSize());
  this->lineAct->GetProperty()->SetLineWidth(2);
  return lineAct;
}
vtkActor* View3d::AddBranchIllustrators()
{
	this->poly = tobj->generateBranchIllustrator();
	this->polymap = vtkSmartPointer<vtkPolyDataMapper>::New();
	polymap->SetInput(this->poly);
	this->bactor = vtkSmartPointer<vtkActor>::New();
	this->bactor->SetMapper(this->polymap);
	this->bactor->SetPickable(0);
	return bactor;
	//ren->AddActor(bactor);
	//bactor->Print(std::cout);

}
void View3d::AddPointsAsPoints(std::vector<TraceBit> vec)
{
	vtkSmartPointer<vtkCubeSource> cube_src = vtkSmartPointer<vtkCubeSource>::New();
	cube_src->SetBounds(-0.2,0.2,-0.2,0.2,-0.2,0.2);
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
  vtkSmartPointer<vtkGlyph3D> glyphs = vtkSmartPointer<vtkGlyph3D>::New();
  glyphs->SetSource(cube_src->GetOutput());
  glyphs->SetInput(point_poly);
  vtkSmartPointer<vtkPolyDataMapper> cubemap = vtkSmartPointer<vtkPolyDataMapper>::New();
  cubemap->SetInput(glyphs->GetOutput());
  cubemap->GlobalImmediateModeRenderingOn();
  vtkSmartPointer<vtkActor> cubeact = vtkSmartPointer<vtkActor>::New();
  cubeact->SetMapper(cubemap);
  cubeact->SetPickable(0);
  cubeact->GetProperty()->SetPointSize(5);
  cubeact->GetProperty()->SetOpacity(.5);
  ren->AddActor(cubeact);

}

void View3d::addAct(vtkActor *Actor)
{ /*  add the actor to render ren */
  ren->AddActor(Actor);
  std::cout << "Added actors "<< std::endl;
 // renWin->Render();
}

void View3d::interact()
{ /*  creates the interactor so picking can occur*/
  sphere = vtkSphereSource::New();
  sphere->SetRadius(3);
  sphereMap = vtkPolyDataMapper::New();
  sphereMap->SetInput(sphere->GetOutput());
  sphereMap->GlobalImmediateModeRenderingOn();

  sphereAct = vtkActor::New();
  sphereAct->SetMapper(sphereMap);
  sphereAct->GetProperty()->SetOpacity(.3);
  sphereAct->VisibilityOff();
  sphereAct->SetPickable(0);            //dont want to pick the sphere itself
  ren->AddActor(sphereAct);           //sphere is used to mark the picks


  
  cell_picker = vtkCellPicker::New();
  //double tolerance = (.003*renWin->GetScreenSize())/ renWin->GetSize();
  //cell_picker->SetTolerance(tolerance);
  cell_picker->SetTolerance(0.004);
  iren->SetPicker(cell_picker);

  isPicked = vtkCallbackCommand::New();
  //isPicked->SetCallback(PickPoint);
  isPicked->SetCallback(PickCell);
  isPicked->SetClientData(this);          //isPicked caller alows observer to intepret click 
  
  keyPress = vtkCallbackCommand::New();
  keyPress->SetCallback(SetMode);
  keyPress->SetClientData(this);

  renWin->Render();
  iren->AddObserver(vtkCommand::RightButtonPressEvent,isPicked);
  iren->AddObserver(vtkCommand::KeyPressEvent,keyPress);
  iren->Start();
  
}
void View3d::SetMode(vtkObject* caller, unsigned long event, void* clientdata, void* callerdata)
{
	bool changes=false;
  View3d* view = (View3d*)clientdata;
  char key = view->iren->GetKeyCode();
  switch (key)
  {
    case 'l':
    {
		if (view->IDList.size()<= 0)
		{
			std::cout<<  "Nothing Selected \n";
		}
		else
		{
		  std::cout<<"These lines\t";
		  for (int i = 0; i < view->IDList.size(); i++)
		  {
			std::cout<<  "\t"<<view->IDList[i];   
		  } 
		  std::cout<< " \t are selected \n" ;
		}
    }
	return;
    break;
	case 'c':
		{
			if (view->IDList.size()<= 0)
			{
				std::cout<<  "Nothing Selected \n";
				return;
			}
			else
			{
				view->IDList.clear();
				std::cout<< "cleared list\n";
				changes = true;
			}
			
		}
		break;
    case 'd':
    {
		if(view->IDList.size()>=1)
		{
		  std::cout<<"selected lines \n";
		  for (int i = 0; i < view->IDList.size(); i++)
		  {
			std::cout<<  "\t"<<view->IDList[i];
			view->deleteTrace(view, reinterpret_cast<TraceLine*>(view->tobj->hashc[view->IDList[i]])); 
		  } 
		  std::cout<< " \t deleted" <<std::endl;
		  changes = true;
		}
		else
		{
			std::cout<<  "Nothing to Delete \n";
			return;
		}
    }
    break;
	case 'm':
    {
		if(view->IDList.size()>=1)
		{		  
			std::sort(view->IDList.begin(), view->IDList.end());
			std::reverse(view->IDList.begin(), view->IDList.end());
			int numTrace = view->IDList.size();		
			int i,j, s=0, exist =0;
			std::vector<TraceLine*> traceList;
			std::cout<< "elements passed \t" << numTrace << std::endl;
			for (i = 0;i<numTrace; i++)
			{
				//std::cout<< "added new trace" <<std::endl;
				traceList.push_back( reinterpret_cast<TraceLine*>(view->tobj->hashc[view->IDList[i]]));
				s=traceList.size()-1;
				exist =0;
				//std::cout<<"trace size "<<  traceList.size() << std::endl;
				if (traceList.size()>0)
				{
					j = 0;
					while ( (j < s)&&(exist==0))
					{	
						//std::cout<<"check trace "<<  s <<" to " <<j << std::endl;						
						if (traceList[s]->GetId()== traceList[j]->GetId())
						{
							traceList.pop_back();	
							//std::cout<<"duplicate trace not added\n";	
							exist=1;
						}
						j++;			
					}
				}
			}
		  view->MinEndPoints(view,traceList);
		  changes = true;
		}
		else
		{
			std::cout<<  "Nothing to merge \n";
			return;
		}
    }
    break;
  case 's':
	  {
		  if(view->IDList.size()>=1)
		  {
		  std::unique( view->IDList.begin(), view->IDList.end());
		  for (int i = 0; i < view->IDList.size(); i++)
			{
			view->tobj->splitTrace(view->IDList[i]);
			} 
		  changes = true;
		  /*view->IDList.clear(); 
		  view->sphereAct->VisibilityOff();
		  view->LineAct();
		  view->renWin->Render();*/
		  }
		  else
		  {
			std::cout<<  "Nothing to split\n";
			return;
		  }
	  }
    break;
  case 'f':
	  {
	  if(view->IDList.size()>=1)
	  {
		  for(int i=0; i< view->IDList.size(); i++)
		  {
			  view->tobj->ReverseSegment(reinterpret_cast<TraceLine*>(view->tobj->hashc[view->IDList[i]]));
		  }
		  changes = true;
	  }
	  }
	  break;
	  case 'w':
		{
			std::string fileName;
			std::cout<<"new filename: \n";
			std::cin>>fileName;
			view->tobj->WriteToSWCFile((char *)fileName.c_str());	
		}
	  break;
  }
  if (changes == true)
  {
	  view->sphereAct->VisibilityOff();
	  view->IDList.clear();
	  //view->ren->RemoveAllViewProps();
	  //view->ren->RemoveViewProp()
	  view->ren->RemoveActor(view->bactor);
	  //view->ren->RemoveActor(view->AddBranchIllustrators());	//branch illistrators as a vtk actor
	  view->LineAct();
	  view->addAct(view->AddBranchIllustrators());
	  //view->AddBranchIllustrators();
	  view->renWin->Render();
	  changes = false;
  }
}

void View3d::PickCell(vtkObject* caller, unsigned long event, void* clientdata, void* callerdata)
{ /*  PickPoint allows fot the point id and coordinates to be returned 
  as well as adding a marker on the last picked point
  R_click to select point on line  */
  //printf("Entered pickCell\n");
  View3d* view = (View3d*)clientdata;       //acess to view3d class so view-> is used to call functions

  int *pos = view->iren->GetEventPosition();
  view->iren->GetPicker()->Pick(pos[0],pos[1],0.0,view->iren->GetRenderWindow()->GetRenderers()->GetFirstRenderer());
  vtkCellPicker *cell_picker = (vtkCellPicker *)view->iren->GetPicker();
  if (cell_picker->GetCellId() == NULL) 
  {
    view->sphereAct->VisibilityOff();     //not working quite yet but sphere will move
  }
  else if(cell_picker->GetViewProp()!=NULL) 
  {
    //double selPt[3];
    double pickPos[3];
    view->cell_picker->GetPickPosition(pickPos);    //this is the coordinates of the pick
    //view->cell_picker->GetSelectionPoint(selPt);    //not so usefull but kept
    unsigned int cell_id = cell_picker->GetCellId();  
    view->IDList.push_back(cell_id);
    TraceLine *tline = reinterpret_cast<TraceLine*>(view->tobj->hashc[cell_id]);
	int t =1;
	view->HighlightSelected(tline, t);
    tline->Getstats();              //prints the id and end coordinates to the command prompt 
    view->sphereAct->SetPosition(pickPos);    //sets the selector to new point
    view->sphereAct->VisibilityOn();      //deleteTrace can turn it off 
    view->poly_line_data->Modified();
	
  }
  view->renWin->Render();             //update the render window
}
void View3d::MinEndPoints(View3d* view,std::vector<TraceLine*> traceList)
{
	int i,j, s=0, exist =0;

	std::vector<compTrace> compList;

	
	for (i=0;i<traceList.size()-1; i++)
	{
		for (j=i+1; j<traceList.size(); j++)
		{
			compTrace newComp;
			newComp.Trace1= traceList[i];
			newComp.Trace2= traceList[j];
			newComp.Trace1->EndPtDist(newComp.Trace2,newComp.endPT1, newComp.endPT2, newComp.dist);
			if (!(newComp.dist>newComp.Trace1->GetSize()/gapTol)&&!(newComp.dist>newComp.Trace2->GetSize()/gapTol))		//&&(newComp.dist<gapMax)
			{
				std::cout<<"added comparison\n";
				compList.push_back(newComp);
				
			}
			/*else
			{
				std::cout<<"distance "
				<< newComp.Trace1->GetId()
				<< " and " << newComp.Trace2->GetId()
				<<" is too large \n";
			}*/			
		}
	}
	if (compList.size()>=1)
	{
		i=0, j = 0;
		while (i<compList.size()-1 )
		{
			exist = 0;
			while ((exist == 0)&&(j<compList.size()-1))
			{
				j++;
				if (compList[i].Trace1->GetId()==compList[j].Trace1->GetId())
				{
					if (compList[i].endPT1==compList[j].endPT1)
					{
						std::cout<<"Conflict "<<compList[i].Trace1->GetId()<<" to "<<compList[i].Trace2->GetId()
							<<" and "<<compList[j].Trace1->GetId()<<" to "<<compList[j].Trace2->GetId()<<std::endl;
						exist=1;}
				}
				else if(compList[i].Trace1->GetId()==compList[j].Trace2->GetId())
				{
					if (compList[i].endPT1==compList[j].endPT2)
					{
						std::cout<<"Conflict "<<compList[i].Trace1->GetId()<<" to "<<compList[i].Trace2->GetId()
							<<" and "<<compList[j].Trace1->GetId()<<" to "<<compList[j].Trace2->GetId()<<std::endl;
						exist=1;}
				}
				else if (compList[i].Trace2->GetId()==compList[j].Trace1->GetId())
				{
					if (compList[i].endPT2==compList[j].endPT1)
					{
						std::cout<<"Conflict "<<compList[i].Trace1->GetId()<<" to "<<compList[i].Trace2->GetId()
							<<" and "<<compList[j].Trace1->GetId()<<" to "<<compList[j].Trace2->GetId()<<std::endl;
						exist=1;}
				}
				else if(compList[i].Trace2->GetId()==compList[j].Trace2->GetId())
				{
					if (compList[i].endPT2==compList[j].endPT2)
					{
						std::cout<<"Conflict "<<compList[i].Trace1->GetId()<<" to "<<compList[i].Trace2->GetId()
							<<" and "<<compList[j].Trace1->GetId()<<" to "<<compList[j].Trace2->GetId()<<std::endl;
						exist=1;}
				}
				
			}
			if (exist==1)
			{
				if (compList[i].dist<compList[j].dist)
				{compList.erase(compList.begin()+j);}
				else
				{compList.erase(compList.begin()+i);}
				std::cout<<"Conflict resolved"<< std::endl;
				j=i;
			}
			else
			{
				i++;
				j=i;
			}
		}
		std::cout<<"trace size "<<  traceList.size() << "\tNumber of computed distances\t" << compList.size()<<std::endl;
		for (i=0;i<compList.size(); i++)
		{
			std::cout<<"Trace\t"<<compList[i].Trace1->GetId()<< "\t compaired to trace\t"<<compList[i].Trace2->GetId() 
				<<" gap size of:\t"<<compList[i].dist<< " endpts "<< compList[i].endPT1<< " and "<< compList[i].endPT2<<std::endl;
			if (compList[i].dist<= gapMax/gapTol)
			{
				tobj->mergeTraces(compList[i].endPT1,compList[i].endPT2);
			}
			else if(compList[i].dist<= gapMax)
			{
				char ans;
				std::cout<<"Distance of: "<< compList[i].dist <<" is greater than: " << gapMax/gapTol<< " Merge y\\n\?";
				std::cin>>ans;
				if (ans=='y')
				{
					tobj->mergeTraces(compList[i].endPT1,compList[i].endPT2);
				}
			}
		}//send to merge
	}
}
void View3d::HighlightSelected(TraceLine* tline, int t)
{
	TraceLine::TraceBitsType::iterator iter = tline->GetTraceBitIteratorBegin();
	TraceLine::TraceBitsType::iterator iterend = tline->GetTraceBitIteratorEnd();

  while(iter!=iterend)
  {
	  poly_line_data->GetPointData()->GetScalars()->SetTuple1(iter->marker,1/t);
	  ++iter;
  }
	

}
void View3d::deleteTrace(View3d* view,TraceLine *tline/*int id*/)
{
  
  /*TraceLine *tline = ;*/
  std::vector<unsigned int> * vtk_cell_ids = tline->GetMarkers();

  vtkIdType ncells; vtkIdType *pts;
  for(int counter=0; counter<vtk_cell_ids->size(); counter++)
  {
    view->poly_line_data->GetCellPoints((vtkIdType)(*vtk_cell_ids)[counter],ncells,pts);
    pts[1]=pts[0];
  }
  std::vector<TraceLine*> *children = tline->GetBranchPointer();
  if(children->size()!=0)
  {
    
    for(int counter=0; counter<children->size(); counter++)
    {
      view->poly_line_data->GetCellPoints((*(*children)[counter]->GetMarkers())[0],ncells,pts);
      pts[1]=pts[0];
      view->tobj->GetTraceLinesPointer()->push_back((*children)[counter]);  
      (*children)[counter]->SetParent(NULL);
    }
    // remove the children now
    children->clear();
  }         //finds and removes children
  // remove from parent
  std::vector<TraceLine*>* siblings;
  if(tline->GetParent()!=NULL)
  {
    siblings=tline->GetParent()->GetBranchPointer();
	if(siblings->size()==2)
	{
		// its not a branch point anymore
		TraceLine *tother1;
		if(tline==(*siblings)[0])
			tother1 = (*siblings)[1];
		else
			tother1 = (*siblings)[0];
		tother1->SetParent(NULL);
		siblings->clear();
		TraceLine::TraceBitsType::iterator iter1,iter2;
		iter1= tline->GetParent()->GetTraceBitIteratorEnd();
		iter2 = tother1->GetTraceBitIteratorBegin();
		iter1--;
		
		view->tobj->mergeTraces((*iter1).marker,(*iter2).marker);
		tline->SetParent(NULL);
		//view->tobj->RemoveTraceLine(tline);// FIXME: No idea why I'm doing this, but only this makes it work.
		delete tline;
		return;
	}
  }
  else
  {
    siblings = view->tobj->GetTraceLinesPointer();
  }
  std::vector<TraceLine*>::iterator iter = siblings->begin();
  std::vector<TraceLine*>::iterator iterend = siblings->end();
  while(iter != iterend)
  {
    if(*iter== tline)
    {
      siblings->erase(iter);
      break;
    }
    ++iter;
  }
  tline->SetParent(NULL);

  
}

class vtkSlider2DCallbackBrightness : public vtkCommand
{
public:
  static vtkSlider2DCallbackBrightness *New() 
    { return new vtkSlider2DCallbackBrightness; }
  virtual void Execute(vtkObject *caller, unsigned long, void*)
    {
      vtkSliderWidget *sliderWidget = 
        reinterpret_cast<vtkSliderWidget*>(caller);
      double value = static_cast<vtkSliderRepresentation *>(sliderWidget->GetRepresentation())->GetValue();

      vtkSmartPointer<vtkPiecewiseFunction> opacityTransferFunction = vtkSmartPointer<vtkPiecewiseFunction>::New();
    opacityTransferFunction->AddPoint(2,0.0);
    //opacityTransferFunction->AddPoint(2+256*(1-value),0.2);
    opacityTransferFunction->AddPoint(50,value);
     this->volume->GetProperty()->SetScalarOpacity(opacityTransferFunction);
    }
  vtkSlider2DCallbackBrightness() {

  }

  vtkSmartPointer<vtkVolume> volume;
};

class vtkSlider2DCallbackContrast : public vtkCommand
{
public:
  static vtkSlider2DCallbackContrast *New() 
    { return new vtkSlider2DCallbackContrast; }
  virtual void Execute(vtkObject *caller, unsigned long, void*)
    {
      vtkSliderWidget *sliderWidget = 
        reinterpret_cast<vtkSliderWidget*>(caller);
      double value = static_cast<vtkSliderRepresentation *>(sliderWidget->GetRepresentation())->GetValue();

      /*vtkSmartPointer<vtkPiecewiseFunction> colorTransferFunction = vtkSmartPointer<vtkPiecewiseFunction>::New();
    colorTransferFunction->AddPoint(1.0,0.0f,0.0f,0.0f);
    colorTransferFunction->AddPoint(255*(1-value),0.5f,0.5f,0.0f);*/
    vtkColorTransferFunction *colorTransferFunction = vtkColorTransferFunction::New();
    colorTransferFunction->AddRGBPoint(0.0, 0.0, 0.0, 0.0);
    colorTransferFunction->AddRGBPoint(255*(1-value),0.5,0.5,0);
     this->volume->GetProperty()->SetColor(colorTransferFunction);
    }
  vtkSlider2DCallbackContrast() {

  }

  vtkSmartPointer<vtkVolume> volume;
};

class vtkSlider2DCallbackContourThreshold : public vtkCommand
{
public:
  static vtkSlider2DCallbackContourThreshold *New() 
    { return new vtkSlider2DCallbackContourThreshold; }
  virtual void Execute(vtkObject *caller, unsigned long, void*)
    {
      vtkSliderWidget *sliderWidget = 
        reinterpret_cast<vtkSliderWidget*>(caller);
      double value = static_cast<vtkSliderRepresentation *>(sliderWidget->GetRepresentation())->GetValue();
    printf("value of threshold = %d\n",int(value*255.0));
    cfilter->Print(std::cout);
    cfilter->GetInput()->Print(std::cout);
    double curr_value = cfilter->GetValue(0);
    if(fabs(curr_value-value*255)>10)
    {
    cfilter->SetValue(0,255*value);
    //ping me if this is really necessary and I'll show you a cross-platform
    //way to do it.  -Zack
    //Sleep(1000);
    }

    }
  vtkSlider2DCallbackContourThreshold() {

  }

  vtkSmartPointer<vtkContourFilter> cfilter;
};
void View3d::AddContourThresholdSliders()
{
  vtkSliderRepresentation2D *sliderRep2 = vtkSliderRepresentation2D::New();
  sliderRep2->SetValue(0.8);
  sliderRep2->SetTitleText("Threshold");
  sliderRep2->GetPoint1Coordinate()->SetCoordinateSystemToNormalizedDisplay();
  sliderRep2->GetPoint1Coordinate()->SetValue(0.2,0.8);
  sliderRep2->GetPoint2Coordinate()->SetCoordinateSystemToNormalizedDisplay();
  sliderRep2->GetPoint2Coordinate()->SetValue(0.8,0.8);
  sliderRep2->SetSliderLength(0.02);
  sliderRep2->SetSliderWidth(0.03);
  sliderRep2->SetEndCapLength(0.01);
  sliderRep2->SetEndCapWidth(0.03);
  sliderRep2->SetTubeWidth(0.005);
  sliderRep2->SetMinimumValue(0.0);
  sliderRep2->SetMaximumValue(1.0);

  vtkSliderWidget *sliderWidget2 = vtkSliderWidget::New();
  sliderWidget2->SetInteractor(iren);
  sliderWidget2->SetRepresentation(sliderRep2);
  sliderWidget2->SetAnimationModeToAnimate();

  vtkSlider2DCallbackContourThreshold *callback_contour = vtkSlider2DCallbackContourThreshold::New();
  callback_contour->cfilter = cfilt;
  sliderWidget2->AddObserver(vtkCommand::InteractionEvent,callback_contour);
  sliderWidget2->EnabledOn();

}


class vtkSubRep: public vtkPlaybackRepresentation
{
public:
  static vtkSubRep * New(){ return new vtkSubRep;}
  virtual void Play() 
  {
    slice_counter++; 
    slice_counter = slice_counter % (im_pointer.size());
    vmapper->SetInput((im_pointer)[slice_counter]);
  }
  virtual void Stop() { 
    cout<< "stopped"; 
    slice_counter=0;
    vmapper->SetInput((im_pointer)[slice_counter]);
  }
  virtual void ForwardOneFrame() 
  {
    slice_counter++; 
    slice_counter = slice_counter % (im_pointer.size());
    vmapper->SetInput((im_pointer)[slice_counter]);
    cout << "forward one frame\n";
  }
  virtual void BackwardOneFrame() 
  {
    slice_counter--; 
    slice_counter = slice_counter+im_pointer.size();
    slice_counter = slice_counter % (im_pointer.size());
    vmapper->SetInput((im_pointer)[slice_counter]);
    cout << "backward one frame\n";
  }
  virtual void JumpToBeginning() 
  {
    slice_counter=0;
    vmapper->SetInput((im_pointer)[slice_counter]);
    cout << "jump to beginning\n";
  }
  virtual void JumpToEnd() 
  {
    slice_counter= im_pointer.size()-1;
    vmapper->SetInput((im_pointer)[slice_counter]);
    cout << "jump to end\n";
  }
  int slice_counter;
  std::vector<vtkSmartPointer<vtkImageData> > im_pointer;
  vtkSmartPointer<vtkVolumeMapper> vmapper;
};
void View3d::AddPlaybackWidget(char *filename)
{

  vtkSubRep *playbackrep = vtkSubRep::New();
  playbackrep->slice_counter=0;
  playbackrep->GetPositionCoordinate()->SetCoordinateSystemToNormalizedDisplay();
  playbackrep->GetPositionCoordinate()->SetValue(0.2,0.1);
  playbackrep->GetPosition2Coordinate()->SetCoordinateSystemToNormalizedDisplay();
  playbackrep->GetPosition2Coordinate()->SetValue(0.8,0.1);
  

  vtkSmartPointer<vtkPlaybackWidget> pbwidget = vtkSmartPointer<vtkPlaybackWidget>::New();
  pbwidget->SetRepresentation(playbackrep);
  pbwidget->SetInteractor(iren);

  const unsigned int Dimension = 3;
  typedef unsigned char  PixelType;
  typedef itk::Image< PixelType, Dimension >   ImageType;
  typedef itk::ImageFileReader< ImageType >    ReaderType;
  ReaderType::Pointer imreader = ReaderType::New();

  imreader->SetFileName(filename);
  imreader->Update();

  
  ImageType::SizeType size = imreader->GetOutput()->GetLargestPossibleRegion().GetSize();

  std::vector<vtkSmartPointer<vtkImageData> > &vtkimarray = playbackrep->im_pointer;
  vtkimarray.reserve(size[2]-1);
  printf("About to create vtkimarray contents in a loop\n");
  for(int counter=0; counter<size[2]-1; counter++)
  {
    vtkimarray.push_back(vtkSmartPointer<vtkImageData>::New());
    vtkimarray[counter]->SetScalarTypeToUnsignedChar();
    vtkimarray[counter]->SetDimensions(size[0],size[1],2);
    vtkimarray[counter]->SetNumberOfScalarComponents(1);
    vtkimarray[counter]->AllocateScalars();
    vtkimarray[counter]->SetSpacing(1/2.776,1/2.776,1);
    memcpy(vtkimarray[counter]->GetScalarPointer(),imreader->GetOutput()->GetBufferPointer()+counter*size[0]*size[1]*sizeof(unsigned char),size[0]*size[1]*2*sizeof(unsigned char));
  }

  printf("finished memcopy in playback widget\n");
  vtkPiecewiseFunction *opacityTransferFunction = vtkPiecewiseFunction::New();
  opacityTransferFunction->AddPoint(2,0.0);
  opacityTransferFunction->AddPoint(20,0.2);

  vtkColorTransferFunction *colorTransferFunction = vtkColorTransferFunction::New();
  colorTransferFunction->AddRGBPoint(0.0, 0.0, 0.0, 0.0);
  colorTransferFunction->AddRGBPoint(20.0,1.0,1.0,1.0);

  
  vtkVolumeProperty *volumeProperty = vtkVolumeProperty::New();
    volumeProperty->SetColor(colorTransferFunction);
    volumeProperty->SetScalarOpacity(opacityTransferFunction);
    volumeProperty->SetInterpolationTypeToLinear();

  vtkSmartPointer<vtkOpenGLVolumeTextureMapper3D> volumeMapper = vtkSmartPointer<vtkOpenGLVolumeTextureMapper3D>::New();
  volumeMapper->SetSampleDistance(0.5);
  volumeMapper->SetInput(vtkimarray[playbackrep->slice_counter]);
  
  vtkSmartPointer<vtkVolume> volume = vtkSmartPointer<vtkVolume>::New();
    volume->SetMapper(volumeMapper);
    volume->SetProperty(volumeProperty);
  volume->SetPickable(0);
  ren->AddVolume(volume);
  printf("Added volume in playback widget");
  vol = volume;
  playbackrep->vmapper=volumeMapper;
  /*playbackrep->im_pointer = &vtkimarray;*/
  playbackrep->Print(std::cout);
  
  renWin->Render();
}
void View3d::AddVolumeSliders()
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
  sliderWidget->SetInteractor(iren);
  sliderWidget->SetRepresentation(sliderRep);
  sliderWidget->SetAnimationModeToAnimate();

  vtkSlider2DCallbackBrightness *callback_brightness = vtkSlider2DCallbackBrightness::New();
  callback_brightness->volume = vol;
  sliderWidget->AddObserver(vtkCommand::InteractionEvent,callback_brightness);
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
  sliderWidget2->SetInteractor(iren);
  sliderWidget2->SetRepresentation(sliderRep2);
  sliderWidget2->SetAnimationModeToAnimate();

  vtkSlider2DCallbackContrast *callback_contrast = vtkSlider2DCallbackContrast::New();
  callback_contrast->volume = vol;
  sliderWidget2->AddObserver(vtkCommand::InteractionEvent,callback_contrast);
  sliderWidget2->EnabledOn();

}
void View3d::readImg(char* sourceFile)
{ /*vtkImageReader2 * readMe = vtkImageReader2::New();
  readMe->SetDataOrigin(0.0,0.0,0.0);
  readMe->SetDataExtent(0,255,0,255,0,255);
  readMe->SetDataSpacing(1,1,1);
  readMe->SetDataByteOrderToLittleEndian();
  readMe->SetFileName(sourceFile);
  std::cout<< sourceFile << '\n';
  readMe->Update();*/
  
  const unsigned int Dimension = 3;
  typedef unsigned char  PixelType;
  typedef itk::Image< PixelType, Dimension >   ImageType;
  typedef itk::ImageFileReader< ImageType >    ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( sourceFile );
  try
  {
    reader->Update();
    }
  catch( itk::ExceptionObject & exp )
    {
    std::cerr << "Exception thrown while reading the input file " << std::endl;
    std::cerr << exp << std::endl;
    //return EXIT_FAILURE;
    }
  std::cout << "Image Read " << std::endl;
//itk-vtk connector
  typedef itk::ImageToVTKImageFilter<ImageType> ConnectorType;
  ConnectorType::Pointer connector= ConnectorType::New();
  connector->SetInput( reader->GetOutput() );

  
  cfilt = vtkSmartPointer<vtkContourFilter>::New();
  cfilt->SetInput(connector->GetOutput());
  cfilt->SetValue(0,10);            // this affects render
  //cfilt = vtkSmartPointer<vtkContourFilter>::New();
  //cfilt = contour;


  
  //vtkSmartPointer<vtkPolyDataNormals> normals = vtkSmartPointer<vtkPolyDataNormals>::New();
  //normals->SetInput(cfilt->GetOutput());
  //normals->SetFeatureAngle(15);

  volMap = vtkSmartPointer<vtkPolyDataMapper>::New();
  volMap->SetInput(cfilt->GetOutput());
  
  volAct = vtkSmartPointer<vtkActor>::New();
  volAct->SetMapper(volMap);
  //volAct->GetProperty()->SetOpacity(.3);
  volAct->GetProperty()->SetRepresentationToWireframe();
  volAct->GetProperty()->SetColor(0.5,0.5,0.5);
  volAct->SetPickable(0);
  //volAct->SetScale(1/2.776,1/2.776,1);
  
  ren->AddActor(volAct);
  renWin->Render();
  std::cout << "contour rendered \n";
}

void View3d::rayCast(char *raySource)
{
  const unsigned int Dimension = 3;
  typedef unsigned char  PixelType;
  typedef itk::Image< PixelType, Dimension >   ImageType;
  typedef itk::ImageFileReader< ImageType >    ReaderType;
  ReaderType::Pointer i2spReader = ReaderType::New();
  i2spReader->SetFileName( raySource );
  try
  {
    i2spReader->Update();
    }
  catch( itk::ExceptionObject & exp )
    {
    std::cerr << "Exception thrown while reading the input file " << std::endl;
    std::cerr << exp << std::endl;
    //return EXIT_FAILURE;
    }
  std::cout << "Image Read " << std::endl;
//itk-vtk connector
  typedef itk::ImageToVTKImageFilter<ImageType> ConnectorType;
  ConnectorType::Pointer connector= ConnectorType::New();
  connector->SetInput( i2spReader->GetOutput() );
  vtkImageToStructuredPoints *i2sp = vtkImageToStructuredPoints::New();
  i2sp->SetInput(connector->GetOutput());


  
  ImageType::SizeType size = i2spReader->GetOutput()->GetLargestPossibleRegion().GetSize();
  vtkSmartPointer<vtkImageData> vtkim = vtkSmartPointer<vtkImageData>::New();
  vtkim->SetScalarTypeToUnsignedChar();
  vtkim->SetDimensions(size[0],size[1],size[2]);
  vtkim->SetNumberOfScalarComponents(1);
  vtkim->AllocateScalars();

  memcpy(vtkim->GetScalarPointer(),i2spReader->GetOutput()->GetBufferPointer(),size[0]*size[1]*size[2]*sizeof(unsigned char));

// Create transfer mapping scalar value to opacity
  vtkPiecewiseFunction *opacityTransferFunction = vtkPiecewiseFunction::New();

  opacityTransferFunction->AddPoint(2,0.0);
  opacityTransferFunction->AddPoint(20,0.1);
  opacityTransferFunction->AddPoint(40,0.1);
  // Create transfer mapping scalar value to color
  // Play around with the values in the following lines to better vizualize data
  vtkColorTransferFunction *colorTransferFunction = vtkColorTransferFunction::New();
    colorTransferFunction->AddRGBPoint(0.0, 0.0, 0.0, 0.0);
  colorTransferFunction->AddRGBPoint(50.0,0.5,0.5,0);

  // The property describes how the data will look
  vtkVolumeProperty *volumeProperty = vtkVolumeProperty::New();
    volumeProperty->SetColor(colorTransferFunction);
    volumeProperty->SetScalarOpacity(opacityTransferFunction);
  //  volumeProperty->ShadeOn();
    volumeProperty->SetInterpolationTypeToLinear();

  vtkSmartPointer<vtkOpenGLVolumeTextureMapper3D> volumeMapper = vtkSmartPointer<vtkOpenGLVolumeTextureMapper3D>::New();
  volumeMapper->SetSampleDistance(0.75);
  volumeMapper->SetInput(vtkim);

  // The volume holds the mapper and the property and
  // can be used to position/orient the volume
  vtkVolume *volume = vtkVolume::New();
    volume->SetMapper(volumeMapper);
    volume->SetProperty(volumeProperty);
  volume->SetPickable(0);
  ren->AddVolume(volume);
  vol = volume;
  renWin->Render();
  std::cout << "RayCast rendered \n";
}
