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

/*
  Farsight ToolKit 3D Viewer: 
  v1: implements render window functionality
    render window creation in "View3D::RenderWin"
    line objects mapped and actors created in "View3D::LineAct"
    actor added and window created in "View3D::addAct"
  v2: add picking functionality:
    access to observing a pick event is in "View3D::interact"
    the interacter calls for "View3D::PickPoint"
  v3: include contourFilter and rayCast renderers
  v4: converted to Qt, member names changed to fit "VTK style" more closely.
  v5: automated functions implemented structure in place for "PACE".
*/

#include "ftkGUI/PlotWindow.h"
#include "ftkGUI/HistoWindow.h"

#include "itkImageFileReader.h"
#include "itkImageToVTKImageFilter.h"

#include <QAction>
#include <QtGui>
#include <QVTKWidget.h>

#include "vtkActor.h"
#include "vtkCallbackCommand.h"
#include "vtkCellArray.h"
#include "vtkCellPicker.h"
#include "vtkColorTransferFunction.h"
#include "vtkCommand.h"
#include "vtkContourFilter.h"
#include "vtkCubeSource.h"
#include "vtkFloatArray.h"
#include "vtkGlyph3D.h"
#include "vtkImageData.h"
#include "vtkImageToStructuredPoints.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkLODActor.h"
#include "vtkOpenGLVolumeTextureMapper3D.h"
#include "vtkPiecewiseFunction.h"
#include "vtkPlaybackRepresentation.h"
#include "vtkPlaybackWidget.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkRendererCollection.h"
#include "vtkRenderWindow.h"
#include "vtkSliderRepresentation2D.h"
#include "vtkSliderWidget.h"
#include "vtkSphereSource.h"
#include "vtkVolume.h"
#include "vtkVolumeProperty.h"

#include "TraceBit.h"
#include "TraceGap.h"
#include "TraceLine.h"
#include "TraceObject.h"
#include "MergeModel.h"
#include "View3DHelperClasses.h"
#include "View3D.h"

View3D::View3D(int argc, char **argv)
{
  this->tobj = new TraceObject;
  int num_loaded = 0;
  this->Volume=0;
  this->somaopacity = 1.0;

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
             strcmp(argv[counter]+len-4,"tiff")==0 ||
			 strcmp(argv[counter]+len-3, "pic")==0||
			 strcmp(argv[counter]+len-3, "PIC")==0)
		  {
			printf("I detected a 3d image file\n");
			this->rayCast(argv[counter]);
			//this->AddVolumeSliders();
		  }
		num_loaded++;
	  }

this->QVTK = 0;
this->tobj->gapTol = .5;
this->tobj->gapMax = 10;
this->smallLine = 5;
this->SelectColor =.1;
this->lineWidth= 2;
this->Initialize();
this->tobj->setSmallLineColor(.25);
this->tobj->setMergeLineColor(.4);
this->Ascending = Qt::AscendingOrder;
}

View3D::~View3D()
{
  if(this->QVTK)
    {
    delete this->QVTK;
    }
  delete this->tobj;
}

void View3D::Initialize()
{
  this->CreateGUIObjects();
  this->CreateLayout();
  this->CreateInteractorStyle();
  this->CreateActors();
  this->resize(640, 480);
  this->setWindowTitle(tr("Trace editor"));
  this->QVTK->GetRenderWindow()->Render();
  this->setupLinkedSpace();
}
void View3D::setupLinkedSpace()
{  
	this->tobj->gapList.clear();
	this->plot = NULL;
	this->histo = NULL;
	this->table = new QTableView();
	this->TreeTable =new QTableView();
	this->MergeGaps = new MergeModel(this->tobj->gapList);
	connect(MergeGaps,SIGNAL(modelChanged()),this, SLOT(updateSelectionHighlights()));
}

/*Set up the components of the interface */
void View3D::CreateGUIObjects()
{
  //Set up the main window's central widget
  this->CentralWidget = new QWidget(this);
  this->setCentralWidget(this->CentralWidget);

  //Set up a QVTK Widget for embedding a VTK render window in Qt.
  this->QVTK = new QVTKWidget(this->CentralWidget);
  this->Renderer = vtkRenderer::New();
  this->QVTK->GetRenderWindow()->AddRenderer(this->Renderer);

  //Set up the menu bar
  this->fileMenu = this->menuBar()->addMenu(tr("&File"));
  this->exitAction = fileMenu->addAction(tr("&Exit"));

  this->somaMenu = this->menuBar()->addMenu(tr("&Somas"));
  this->loadSoma = somaMenu->addAction(tr("&Load Soma Data"));
  this->loadSeed = somaMenu->addAction(tr("&Load Seed Data"));
  this->somaSettings = somaMenu->addAction(tr("&Settings"));

  this->viewSomas = somaMenu->addAction(tr("&Show Somas"));
  this->viewSomas->setCheckable(true);
  this->viewSomas->setDisabled(true);

  this->viewSeeds = somaMenu->addAction(tr("&Show Seeds"));
  this->viewSeeds->setCheckable(true);
  this->viewSeeds->setDisabled(true);

  //Set up the buttons that the user will use to interact with this program. 
  this->ListButton = new QPushButton("&list selections", this->CentralWidget);
  this->ClearButton = new QPushButton("&clear selection", this->CentralWidget);
  this->DeleteButton = new QPushButton("&delete trace", this->CentralWidget);
  this->MergeButton = new QPushButton("&merge traces", this->CentralWidget);
  this->SplitButton = new QPushButton("&split trace", this->CentralWidget);
  this->FlipButton = new QPushButton("&flip trace direction", this->CentralWidget);
  this->WriteButton = new QPushButton("&write to .swc file", this->CentralWidget);
  this->SettingsButton = new QPushButton("&edit settings", this->CentralWidget);
  this->AutomateButton = new QPushButton("&Automatic Selection", this->CentralWidget);

  //Setup the tolerance settings editing window
  this->SettingsWidget = new QWidget();
  QIntValidator *intValidator = new QIntValidator(1, 100, this->SettingsWidget);
  QDoubleValidator *doubleValidator =
    new QDoubleValidator(0, 100, 2, this->SettingsWidget);
  this->MaxGapField = new QLineEdit(this->SettingsWidget);
  this->MaxGapField->setValidator(intValidator);
  this->GapToleranceField = new QLineEdit(this->SettingsWidget);
  this->GapToleranceField->setValidator(intValidator);
  this->LineLengthField = new QLineEdit(this->SettingsWidget);
  this->LineLengthField->setValidator(doubleValidator);
  this->ColorValueField = new QLineEdit(this->SettingsWidget);
  this->ColorValueField->setValidator(doubleValidator);
  this->LineWidthField = new QLineEdit(this->SettingsWidget);
  this->LineWidthField->setValidator(doubleValidator);
  this->ApplySettingsButton = new QPushButton("&Apply", this->SettingsWidget);
  this->CancelSettingsButton = new QPushButton("&Cancel", this->SettingsWidget);

  //Load soma window setup
  this->LoadSomaWidget = new QWidget();
  this->LoadSomaWidget->setWindowTitle(tr("Soma File Load"));
  //this->somaRegex->setPattern(".+\\.([a-z]+)"); //Set up regular expression for soma file field
  //this->somaRegex->setCaseSensitivity(Qt::CaseInsensitive);
  this->somaRegex.setPattern(".+\\.([a-z]+)"); 
  this->somaRegex.setCaseSensitivity(Qt::CaseInsensitive);
  this->SomaFileField = new QLineEdit(this->LoadSomaWidget);

  this->OpenSomaButton = new QPushButton("&Open", this->LoadSomaWidget);
  this->CancelSomaButton = new QPushButton("&Cancel", this->LoadSomaWidget);
  this->BrowseSomaButton = new QPushButton("&Browse", this->LoadSomaWidget);

  //Load seed window setup
  this->LoadSeedWidget = new QWidget();
  this->LoadSeedWidget->setWindowTitle(tr("Seed File Load"));
  this->seedRegex.setPattern(".+\\.([a-z]+)");
  this->seedRegex.setCaseSensitivity(Qt::CaseInsensitive);
  this->SeedFileField = new QLineEdit(this->LoadSeedWidget);

  this->OpenSeedButton = new QPushButton("&Open", this->LoadSeedWidget);
  this->CancelSeedButton = new QPushButton("&Cancel", this->LoadSeedWidget);
  this->BrowseSeedButton = new QPushButton("&Browse", this->LoadSeedWidget);

  //Soma settings window setup
  this->SomaSettingsWidget = new QWidget();
  this->SomaSettingsWidget->setWindowTitle(tr("Soma Settings"));

  this->SomaOpacityField = new QLineEdit(this->SomaSettingsWidget);
  this->SomaOpacityField->setValidator(doubleValidator);

  this->ApplySomaSettingsButton = new QPushButton("&Apply", this->SomaSettingsWidget);
  this->CancelSomaSettingsButton = new QPushButton("&Cancel", this->SomaSettingsWidget);
  

  //Setup connections for the menu bar
  connect(this->exitAction, SIGNAL(triggered()), this, SLOT(close()));
  connect(this->loadSoma, SIGNAL(triggered()), this, SLOT(ShowLoadSomaWindow()));
  connect(this->loadSeed, SIGNAL(triggered()), this, SLOT(ShowLoadSeedWindow()));
  connect(this->viewSomas, SIGNAL(triggered()), this, SLOT(ToggleSomas()));
  connect(this->viewSeeds, SIGNAL(triggered()), this, SLOT(ToggleSeeds()));
  connect(this->somaSettings, SIGNAL(triggered()), this, SLOT(ShowSomaSettingsWindow()));

  //Setup connections for the buttons
  connect(this->ListButton, SIGNAL(clicked()), this, SLOT(ListSelections()));
  connect(this->ClearButton, SIGNAL(clicked()), this, SLOT(ClearSelection()));
  connect(this->DeleteButton, SIGNAL(clicked()), this, SLOT(DeleteTraces()));
  connect(this->MergeButton, SIGNAL(clicked()), this, SLOT(MergeTraces()));
  connect(this->SplitButton, SIGNAL(clicked()), this, SLOT(SplitTraces()));
  connect(this->FlipButton, SIGNAL(clicked()), this, SLOT(FlipTraces()));
  connect(this->WriteButton, SIGNAL(clicked()), this, SLOT(WriteToSWCFile()));
  connect(this->AutomateButton, SIGNAL(clicked()), this, SLOT(SLine()));
  connect(this->SettingsButton, SIGNAL(clicked()), this,
    SLOT(ShowSettingsWindow()));
  connect(this->ApplySettingsButton, SIGNAL(clicked()), this,
    SLOT(ApplyNewSettings()));
  connect(this->CancelSettingsButton, SIGNAL(clicked()), this,
    SLOT(HideSettingsWindow()));

  //Set up connections for the Load Soma Data window
  connect(this->OpenSomaButton, SIGNAL(clicked()), this, SLOT(GetSomaFile()));
  connect(this->CancelSomaButton, SIGNAL(clicked()), this, SLOT(HideLoadSomaWindow()));
  connect(this->BrowseSomaButton, SIGNAL(clicked()), this, SLOT(GetSomaPath()));

  //Set up connections for the Load Seed Data window
  connect(this->OpenSeedButton, SIGNAL(clicked()), this, SLOT(GetSeedFile()));
  connect(this->CancelSeedButton, SIGNAL(clicked()), this, SLOT(HideLoadSeedWindow()));
  connect(this->BrowseSeedButton, SIGNAL(clicked()), this, SLOT(GetSeedPath()));

  //Set up connections for the Soma Settings window
  connect(this->CancelSomaSettingsButton, SIGNAL(clicked()), this, SLOT(HideSomaSettingsWindow()));
  connect(this->ApplySomaSettingsButton, SIGNAL(clicked()), this, SLOT(ApplySomaSettings()));
}

void View3D::CreateLayout()
{
  //layout for the main window
	QGridLayout *buttonLayout = new QGridLayout();
	buttonLayout->addWidget(this->AutomateButton, 1, 0);
	buttonLayout->addWidget(this->ListButton, 2, 0);
	buttonLayout->addWidget(this->ClearButton, 3, 0);
	buttonLayout->addWidget(this->DeleteButton, 4, 0);
	buttonLayout->addWidget(this->MergeButton, 5, 0);
	buttonLayout->addWidget(this->SplitButton, 6, 0);
	buttonLayout->addWidget(this->FlipButton, 7, 0);
	buttonLayout->addWidget(this->SettingsButton, 9, 0);
	buttonLayout->addWidget(this->WriteButton, 10, 0);
	buttonLayout->setSpacing(10);
	QGridLayout *viewerLayout = new QGridLayout(this->CentralWidget);
	viewerLayout->addWidget(this->QVTK, 0, 1);
	viewerLayout->addLayout(buttonLayout, 0, 0);

  //layout for the settings window
  QGridLayout *settingsLayout = new QGridLayout(this->SettingsWidget);
  QLabel *maxGapLabel = new QLabel(tr("Maximum gap length:"));
  settingsLayout->addWidget(maxGapLabel, 0, 0);
  settingsLayout->addWidget(this->MaxGapField, 0, 1);
  QLabel *gapToleranceLabel = new QLabel(tr("Gap length tolerance:"));
  settingsLayout->addWidget(gapToleranceLabel, 1, 0);
  settingsLayout->addWidget(this->GapToleranceField, 1, 1);
  QLabel *lineLengthLabel = new QLabel(tr("Small line length:"));
  settingsLayout->addWidget(lineLengthLabel, 2, 0);
  settingsLayout->addWidget(this->LineLengthField, 2, 1);
  QLabel *colorValueLabel = new QLabel(tr("Color value RGB scalar 0 to 1:"));
  settingsLayout->addWidget(colorValueLabel, 3, 0);
  settingsLayout->addWidget(this->ColorValueField, 3, 1);
  QLabel *lineWidthLabel = new QLabel(tr("Line width:"));
  settingsLayout->addWidget(lineWidthLabel, 4, 0);
  settingsLayout->addWidget(this->LineWidthField, 4, 1);
  settingsLayout->addWidget(this->ApplySettingsButton, 5, 0);
  settingsLayout->addWidget(this->CancelSettingsButton, 5, 1); 

  //Layout for the Load Soma Window
  QGridLayout *LoadSomaLayout = new QGridLayout(this->LoadSomaWidget);
  QLabel *SomaFileLabel = new QLabel(tr("Soma Image File:"));
  LoadSomaLayout->addWidget(SomaFileLabel, 0, 0);
  LoadSomaLayout->addWidget(this->SomaFileField, 0, 1);
  LoadSomaLayout->addWidget(this->BrowseSomaButton, 0, 2);
  LoadSomaLayout->addWidget(this->OpenSomaButton, 1, 0);
  LoadSomaLayout->addWidget(this->CancelSomaButton, 1, 1);

  //Layout for the Load Seed Window
  QGridLayout *LoadSeedLayout = new QGridLayout(this->LoadSeedWidget);
  QLabel *SeedFileLabel = new QLabel(tr("Seed Point File: (.txt)"));
  LoadSeedLayout->addWidget(SeedFileLabel, 0, 0);
  LoadSeedLayout->addWidget(this->SeedFileField, 0, 1);
  LoadSeedLayout->addWidget(this->BrowseSeedButton, 0, 2);
  LoadSeedLayout->addWidget(this->OpenSeedButton, 1, 0);
  LoadSeedLayout->addWidget(this->CancelSeedButton, 1, 1);

  //Layout for the Soma Settings Window
  QGridLayout *SomaSettingsLayout = new QGridLayout(this->SomaSettingsWidget);
  QLabel *SomaOpacityLabel = new QLabel(tr("Soma Opacity (0 to 1):"));
  SomaSettingsLayout->addWidget(SomaOpacityLabel, 0, 0);
  SomaSettingsLayout->addWidget(this->SomaOpacityField, 0, 1);
  SomaSettingsLayout->addWidget(this->ApplySomaSettingsButton, 1, 0);
  SomaSettingsLayout->addWidget(this->CancelSomaSettingsButton, 1, 1);
}

void View3D::CreateInteractorStyle()
{
  this->Interactor = this->QVTK->GetRenderWindow()->GetInteractor();
	//keep mouse command observers, but change the key ones
  this->keyPress = vtkCallbackCommand::New();
  this->keyPress->SetCallback(HandleKeyPress);
  this->keyPress->SetClientData(this);

  this->Interactor->RemoveObservers(vtkCommand::KeyPressEvent);
	this->Interactor->RemoveObservers(vtkCommand::KeyReleaseEvent);
	this->Interactor->RemoveObservers(vtkCommand::CharEvent);
	this->Interactor->AddObserver(vtkCommand::KeyPressEvent, this->keyPress);

  //use trackball control for mouse commands
  vtkSmartPointer<vtkInteractorStyleTrackballCamera> style =
    vtkInteractorStyleTrackballCamera::New();
  this->Interactor->SetInteractorStyle(style);
  this->CellPicker = vtkCellPicker::New();
  this->CellPicker->SetTolerance(0.004);
  this->Interactor->SetPicker(this->CellPicker);
  this->isPicked = vtkCallbackCommand::New();
  this->isPicked->SetCallback(PickCell);

  //isPicked caller allows observer to intepret click 
  this->isPicked->SetClientData(this);            
  this->Interactor->AddObserver(vtkCommand::RightButtonPressEvent,isPicked);
}

void View3D::CreateActors()
{
  this->LineActor = vtkActor::New();
  this->LineMapper = vtkPolyDataMapper::New();
  this->UpdateLineActor();
  this->LineActor->SetPickable(1);
  this->Renderer->AddActor(this->LineActor);

  this->UpdateBranchActor();
  this->Renderer->AddActor(this->BranchActor);
 
  if(this->Volume!=NULL)
  {
	  this->Renderer->AddVolume(this->Volume);
	  this->AddVolumeSliders();
  }

  //this->Renderer->AddActor(this->VolumeActor);

  //sphere is used to mark the picks
  this->CreateSphereActor();
  Renderer->AddActor(this->SphereActor);
}

void View3D::CreateSphereActor()
{
  this->Sphere = vtkSphereSource::New();
  this->Sphere->SetRadius(3);
  this->SphereMapper = vtkPolyDataMapper::New();
  this->SphereMapper->SetInput(this->Sphere->GetOutput());
  this->SphereMapper->GlobalImmediateModeRenderingOn();

  this->SphereActor = vtkActor::New();
  this->SphereActor->SetMapper(this->SphereMapper);
  this->SphereActor->GetProperty()->SetOpacity(.3);
  this->SphereActor->VisibilityOff();
  this->SphereActor->SetPickable(0);            //dont want to pick the sphere itself
}

/* update settings */
void View3D::ShowSettingsWindow()
{
  //make sure the values in the input fields are up-to-date
  this->MaxGapField->setText(QString::number(this->tobj->gapMax));
  this->GapToleranceField->setText(QString::number(this->tobj->gapTol));
  this->LineLengthField->setText(QString::number(this->smallLine));
  this->ColorValueField->setText(QString::number(this->SelectColor));
  this->LineWidthField->setText(QString::number(this->lineWidth));
  this->SettingsWidget->show();
}

void View3D::ApplyNewSettings()
{
  this->tobj->gapMax = this->MaxGapField->text().toInt();
  this->tobj->gapTol = this->GapToleranceField->text().toDouble();
  this->smallLine = this->LineLengthField->text().toFloat();
  this->SelectColor = this->ColorValueField->text().toDouble();
  this->lineWidth = this->LineWidthField->text().toFloat();
  //this->Rerender();
  this->SettingsWidget->hide();
}

void View3D::HideSettingsWindow()
{
  this->SettingsWidget->hide();
}
bool View3D::setTol()
{
	bool change = false;
	char select=0;
	std::cout<<"Settings Configuration:\n gap (t)olerance:\t" <<this->tobj->gapTol
		<<"\ngap (m)ax:\t"<<this->tobj->gapMax
		<<"\n(s)mall line:\t"<< smallLine
		<<"\nselection (c)olor:\t"<<SelectColor
		<<"\nline (w)idth:\t"<<lineWidth
		<<"\n(e)nd edit settings\n";
	while (select !='e')
	{
		std::cout<< "select option:\t"; 
		std::cin>>select;
		switch(select)
		{
		case 'm':
			{
				int newMax;
				std::cout<< "maximum gap length\n";
				std::cin>>newMax;
				if (newMax!=this->tobj->gapMax)
				{
					this->tobj->gapMax=newMax;
					change= true;
				}
				break;
			}//end of 'm'
		case 't':
			{
				double newTol;
				std::cout<< "gap length tolerance\n";
				std::cin>>newTol;
				if (newTol!=this->tobj->gapTol)
				{
					this->tobj->gapTol=newTol;
					change= true;
				}
				break;
			}//end of 't'
		case 's':
			{
				float newSmall;
				std::cout<< "small line length\n";
				std::cin>>newSmall;
				if (newSmall!=smallLine)
				{
					smallLine=newSmall;
					change= true;
				}
				break;
			}// end of 's'
		case 'c':
			{
				double newColor;
				std::cout<< "color value RGB scalar 0 to 1\n";
				std::cin>>newColor;
				if (newColor!=SelectColor)
				{
					SelectColor=newColor;
					change= true;
				}
				break;
			}//end of 'c'
		case 'w':
			{
				float newWidth;
				std::cout<<"line Width\n";
				std::cin>>newWidth;
				if (newWidth!=lineWidth)
				{
					lineWidth=newWidth;
					change= true;
				}
				break;
			}
		}//end of switch
	}// end of while
	if (change== true)
	{std::cout<<"Settings Configuration are now:\n gap tollerance:\t" <<this->tobj->gapTol
		<<"\ngap max:\t"<<this->tobj->gapMax
		<<"\nsmall line:\t"<< smallLine
		<<"\nselection color:\t"<<SelectColor
		<<"\nline width:\t"<<lineWidth;
	}
	return change;
}
/*	picking	*/
void View3D::PickCell(vtkObject* caller, unsigned long event, void* clientdata, void* callerdata)
{ /*  PickPoint allows fot the point id and coordinates to be returned 
  as well as adding a marker on the last picked point
  R_click to select point on line  */
  //printf("Entered pickCell\n");
  View3D* view = (View3D*)clientdata;       //acess to view3d class so view-> is used to call functions

  int *pos = view->Interactor->GetEventPosition();
  view->Interactor->GetPicker()->Pick(pos[0],pos[1],0.0,view->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer());
  vtkCellPicker *cell_picker = (vtkCellPicker *)view->Interactor->GetPicker();
  if (cell_picker->GetCellId() == -1) 
  {
    view->SphereActor->VisibilityOff();     //not working quite yet but sphere will move
  }
  else if(cell_picker->GetViewProp()!=NULL) 
  {
    //double selPt[3];
    double pickPos[3];
    view->CellPicker->GetPickPosition(pickPos);    //this is the coordinates of the pick
    //view->cell_picker->GetSelectionPoint(selPt);    //not so usefull but kept
    unsigned int cell_id = cell_picker->GetCellId();  
    view->IDList.push_back(cell_id);
    TraceLine *tline = reinterpret_cast<TraceLine*>(view->tobj->hashc[cell_id]);
	
	view->HighlightSelected(tline, view->SelectColor);
    tline->Getstats();              //prints the id and end coordinates to the command prompt 
    view->SphereActor->SetPosition(pickPos);    //sets the selector to new point
    view->SphereActor->VisibilityOn();      //deleteTrace can turn it off 
    view->poly_line_data->Modified();

    //update the head Qt view here too...
  }
  view->QVTK->GetRenderWindow()->Render();             //update the render window
}


void View3D::HighlightSelected(TraceLine* tline, double color)
{
	TraceLine::TraceBitsType::iterator iter = tline->GetTraceBitIteratorBegin();
	TraceLine::TraceBitsType::iterator iterend = tline->GetTraceBitIteratorEnd();
	if (color == -1)
	{
		color = tline->getTraceColor();
	}
  while(iter!=iterend)
  {
	  //poly_line_data->GetPointData()->GetScalars()->SetTuple1(iter->marker,1/t);
	  poly_line_data->GetPointData()->GetScalars()->SetTuple1(iter->marker,color);
	  ++iter;
  }
	

}

void View3D::Rerender()
{
  this->SphereActor->VisibilityOff();
  this->IDList.clear();
  this->Renderer->RemoveActor(this->BranchActor);
  //this->Renderer->RemoveActor(this->VolumeActor);
  this->UpdateLineActor();
  this->UpdateBranchActor();
  this->Renderer->AddActor(this->BranchActor);
  //this->Renderer->AddActor(this->VolumeActor);
  this->QVTK->GetRenderWindow()->Render();
}

void View3D::UpdateLineActor()
{
  this->poly_line_data = this->tobj->GetVTKPolyData();
  this->LineMapper->SetInput(this->poly_line_data);
  this->LineActor->SetMapper(this->LineMapper);
  this->LineActor->GetProperty()->SetColor(0,1,0);
  this->LineActor->GetProperty()->SetPointSize(2);
  this->LineActor->GetProperty()->SetLineWidth(lineWidth);
}

void View3D::UpdateBranchActor()
{
	this->poly = tobj->generateBranchIllustrator();
	this->polymap = vtkSmartPointer<vtkPolyDataMapper>::New();
	polymap->SetInput(this->poly);
	this->BranchActor = vtkSmartPointer<vtkActor>::New();
	this->BranchActor->SetMapper(this->polymap);
	this->BranchActor->SetPickable(0);
	//Renderer->AddActor(BranchActor);
	//BranchActor->Print(std::cout);

}
void View3D::AddPointsAsPoints(std::vector<TraceBit> vec)
{
	vtkSmartPointer<vtkCubeSource> cube_src = vtkSmartPointer<vtkCubeSource>::New();
	cube_src->SetBounds(-0.2,0.2,-0.2,0.2,-0.2,0.2);
  vtkSmartPointer<vtkPolyData> point_poly = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPoints> points=vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> cells=vtkSmartPointer<vtkCellArray>::New();
  for(unsigned int counter=0; counter<vec.size(); counter++)
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
  Renderer->AddActor(cubeact);

}

void View3D::HandleKeyPress(vtkObject* caller, unsigned long event,
                            void* clientdata, void* callerdata)
{
  View3D* view = (View3D*)clientdata;
  char key = view->Interactor->GetKeyCode();
  switch (key)
    {
    case 'l':
      view->ListSelections();
      break;

    case 'c':
      view->ClearSelection();
      break;

    case 'd':
		view->tobj->gapList.clear();
      view->DeleteTraces();
      break;
    
	  case 'm':
      view->MergeTraces();
      break;
    
    case 's':
      view->SplitTraces();
      break;

    case 'f':
      view->FlipTraces();
      break;

	  case 'w':
      view->WriteToSWCFile();
	  break;

	  case 't':
      view->ShowSettingsWindow();
      break;

	  case 'a':
		  view->SLine();
		  /*std::cout<<"select small lines\n";
		  view->tobj->FindMinLines(view->smallLine);
		  view->Rerender();*/
	    break;

	  case 'q':
		  view->updateSelectionHighlights();
		  break;

    case '-':
      if(view->IDList.size()>=1)
        {
        TraceLine* tline=reinterpret_cast<TraceLine*>(
          view->tobj->hashc[view->IDList[view->IDList.size()-1]]);
        view->HighlightSelected(tline, tline->getTraceColor()-.25);
        view->IDList.pop_back();
        cout<< " These lines: ";
        for (unsigned int i = 0; i < view->IDList.size(); i++)
          {
          cout<<  "\t"<<view->IDList[i];   
          } 
        cout<< " \t are selected \n" ;   
        }
      break;

    default:
      break;
    }
}
/*	Actions		*/
void View3D::SLine()
{
	int numLines, i;
	this->tobj->FindMinLines(this->smallLine);
	numLines= this->tobj->SmallLines.size();
	this->Rerender();
	QMessageBox Myquestion;
	Myquestion.setText("Number of selected small lines:  " 
		+ QString::number(numLines));
	Myquestion.setInformativeText("Delete these small lines?" );
	Myquestion.setStandardButtons(QMessageBox::Yes|QMessageBox::No);
	Myquestion.setDefaultButton(QMessageBox::Yes);
	int ret = Myquestion.exec();
	switch (ret) 
	{	
	case QMessageBox::Yes:
	{
		for (i=0;i<numLines;i++)
		{	//std::cout << "Deleted line:" << i<< std::endl;					
			//this->tobj->SmallLines[i]->Getstats();
			this->tobj->RemoveTraceLine(this->tobj->SmallLines[i]);
		}
		this->tobj->SmallLines.clear();
	}
	break;
	case QMessageBox::No:
	 {
		 this->tobj->SmallLines.clear();
	 }
	 break;
	}
	this->Rerender();
}


void View3D::ListSelections()
{
  QMessageBox *selectionInfo = new QMessageBox;
  QString listText;
  QString selectedText;
  if (this->IDList.size()<= 0)
    {
    listText = tr("No traces selected");
    }
  else
    {
    listText += QString::number(this->IDList.size()) + " lines are selected\n";
    for (unsigned int i = 0; i < this->IDList.size(); i++)
      {
      selectedText += QString::number(this->IDList[i]) + "\n";   
      } 
    }
  selectionInfo->setText(listText);
  selectionInfo->setDetailedText(selectedText);
  selectionInfo->show();
}

void View3D::ClearSelection()
{
	QMessageBox *selectionInfo = new QMessageBox();
	QString selectText;
	if (this->IDList.size()<= 0)
    {
		selectText=tr("\t\tNothing Selected\t\t");
    }
	else
    {
		this->IDList.clear();
		selectText=tr("\t\tcleared list\t\t");
		this->Rerender();
    }
	selectionInfo->setText(selectText);
	selectionInfo->show();
}

/*	delete traces functions	*/
void View3D::DeleteTraces()
{
  if(this->IDList.size()>=1)
    {
    std::cout<<"selected lines \n";
    for (unsigned int i = 0; i < this->IDList.size(); i++)
      {
      std::cout<<  "\t"<<this->IDList[i];
      this->DeleteTrace(reinterpret_cast<TraceLine*>(this->tobj->hashc[this->IDList[i]])); 
      } 
    std::cout<< " \t deleted" <<std::endl;
    this->Rerender();
    }
  else
    {
    std::cout<<  "Nothing to Delete \n";
    }
}

void View3D::DeleteTrace(TraceLine *tline)
{
  std::vector<unsigned int> * vtk_cell_ids = tline->GetMarkers();

  vtkIdType ncells; vtkIdType *pts;
  for(unsigned int counter=0; counter<vtk_cell_ids->size(); counter++)
    {
    this->poly_line_data->GetCellPoints((vtkIdType)(*vtk_cell_ids)[counter],ncells,pts);
    pts[1]=pts[0];
    }
  std::vector<TraceLine*> *children = tline->GetBranchPointer();
  if(children->size()!=0)
    {
    for(unsigned int counter=0; counter<children->size(); counter++)
      {
      this->poly_line_data->GetCellPoints((*(*children)[counter]->GetMarkers())[0],ncells,pts);
      pts[1]=pts[0];
      this->tobj->GetTraceLinesPointer()->push_back((*children)[counter]);  
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
        { 
        tother1 = (*siblings)[1];
        }
      else
        {
        tother1 = (*siblings)[0];
        }
      tother1->SetParent(NULL);
      siblings->clear();
      TraceLine::TraceBitsType::iterator iter1,iter2;
      iter1= tline->GetParent()->GetTraceBitIteratorEnd();
      iter2 = tother1->GetTraceBitIteratorBegin();
      iter1--;
		
      this->tobj->mergeTraces((*iter1).marker,(*iter2).marker);
      tline->SetParent(NULL);
      delete tline;
      return;
	    }
    }
  else
    {
    siblings = this->tobj->GetTraceLinesPointer();
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

/*	merging functions	*/
void View3D::MergeTraces()
{
	if (this->tobj->gapList.size() > 1)
	{
		if(this->plot)
		{
			this->plot->close();
			//delete plot;
		}
		if (this->table)
		{
			this->table->close();
		}
		this->MergeSelectedTraces();
		//this->Rerender();
		this->tobj->gapList.clear();
	}
	else
	{
		int conflict =0;
	  if(this->IDList.size()>=1)
		{		  
		std::sort(this->IDList.begin(), this->IDList.end());
		std::reverse(this->IDList.begin(), this->IDList.end());
		int numTrace = this->IDList.size();		
		int i,j, s=0, exist =0;
		std::vector<TraceLine*> traceList;
		//std::cout<< "elements passed \t" << numTrace << std::endl;
		for (i = 0;i<numTrace; i++)
		  {
		  traceList.push_back( reinterpret_cast<TraceLine*>(
			this->tobj->hashc[this->IDList[i]]));
		  s=traceList.size()-1;
		  exist = 0;
		  if (traceList.size()>0)
			{
			j = 0;
			while ( (j < s)&&(exist==0))
			  {	
			  if (traceList[s]->GetId()== traceList[j]->GetId())
				{
				traceList.pop_back();	
				exist = 1;
				}
			  j++;			
			  }
			}
		  }
		conflict = this->tobj->createGapLists(traceList);
		}
	  else
		{
			conflict = this->tobj->createGapLists(this->tobj->GetTraceLines());
		}	
	  unsigned int i; 
	  QMessageBox *MergeInfo = new QMessageBox;
		if (this->tobj->gapList.size() > 1)
		{
			for (i=0;i<this->tobj->gapList.size(); i++)
			{
				this->tobj->gapList[i]->compID = i;
				//this->HighlightSelected(this->tobj->gapList[i]->Trace1, .25);
				//this->HighlightSelected(this->tobj->gapList[i]->Trace2, .25);
			}
			MergeInfo->setText("\nNumber of computed distances:\t" 
				+ QString::number(this->tobj->gapList.size())
				+"\nConflicts resolved:\t" + QString::number(conflict)
				+"\nEdit selection or press merge again");
			this->ShowMergeStats();
		}//end if this->tobj->gapList size > 1
		else 
		{
			if (this->tobj->gapList.size() ==1)
			{		
				tobj->mergeTraces(this->tobj->gapList[0]->endPT1,this->tobj->gapList[0]->endPT2);
				this->Rerender();
				MergeInfo->setText(this->myText + "\nOne Trace merged");
			}
			else
			{
				this->Rerender();
				MergeInfo->setText("\nNo merges possible, set higher tolerances\n"); 
			}	
		}		
		MergeInfo->show();
		this->myText.clear();
		this->poly_line_data->Modified();
		this->QVTK->GetRenderWindow()->Render();
	}//end else size
}

void View3D::ShowMergeStats()
{
	this->MergeGaps->SetTraceGaps(this->tobj->gapList);
	this->table->setModel(this->MergeGaps->GetModel());
	this->table->setSelectionModel(this->MergeGaps->GetSelectionModel()); 
	this->table->setSelectionBehavior(QAbstractItemView::SelectRows);
	this->table->sortByColumn(8,Ascending);
	this->table->update();
	this->table->show();
	this->plot = new PlotWindow(this->MergeGaps->GetSelectionModel());
	this->plot->show();
}
void View3D::updateSelectionHighlights()
{
	bool selected = false;
	int curID;
	std::vector<int> GapIDs = this->MergeGaps->GetSelectedGapIDs();
	for (int i = 0; i < this->tobj->gapList.size(); i++)
	{
		curID = this->tobj->gapList[i]->compID;
		selected = false;	int j =0;
		while(!selected && j < GapIDs.size())
		{
			if ( curID == GapIDs[j])
			{
				selected = true;
			}
			else
			{
				j++;
			}
		}//end gapids size j
		if (selected == true)
		{
			this->HighlightSelected(this->tobj->gapList[i]->Trace1, .15);
			this->HighlightSelected(this->tobj->gapList[i]->Trace2, .15);
		}//end true
		else
		{
			this->HighlightSelected(this->tobj->gapList[i]->Trace1, .25);
			this->HighlightSelected(this->tobj->gapList[i]->Trace2, .25);
		}//end else true
	}//end for i
	this->poly_line_data->Modified();
	this->QVTK->GetRenderWindow()->Render();
}
void View3D::MergeSelectedTraces()
{
	this->updateSelectionHighlights(); //if it didnt
	std::vector<int> GapIDs = this->MergeGaps->GetSelectedGapIDs();
	bool selected = false;
	int curID;
	QMessageBox *MergeInfo = new QMessageBox;
	double currentAngle=0;
	QPushButton *mergeAll;
	QPushButton *mergeNone;
	unsigned int i=0, j=0,mergeCount=0;
	MergeInfo->setText("Merge Function");
	if (GapIDs.size() > 1)
	{
		for (i = 0; i < this->tobj->gapList.size(); i++)
		{
			curID = this->tobj->gapList[i]->compID;
			selected = false;	
			j =0;
			while(!selected && j < GapIDs.size())
			{
				if ( curID == GapIDs[j])
				{
					selected = true;
				}
				else
				{
					j++;
				}
			}//end gapids size j
			if (selected == true)
			{
				this->tobj->mergeTraces(this->tobj->gapList[i]->endPT1,this->tobj->gapList[i]->endPT2);
			}
		}	
		MergeInfo->setText("merged " + QString::number(GapIDs.size()) + " traces.");	
		MergeInfo->show();
		MergeInfo->exec();
	}
	else
	{
		for (i=0;i<this->tobj->gapList.size(); i++)
		{
			currentAngle=fabs(this->tobj->gapList[i]->angle); 
			//if (currentAngle>=PI/2) 
			//{
			//	currentAngle= PI-currentAngle;
			//}
			if 	((this->tobj->gapList[i]->dist<= this->tobj->gapMax*this->tobj->gapTol&& (currentAngle < 1.1))||
				(this->tobj->gapList[i]->dist<= this->tobj->gapMax && (currentAngle < .6)))
	  		{
				this->dtext+= "\nTrace " + QString::number(this->tobj->gapList[i]->Trace1->GetId());
				this->dtext+= " and "+ QString::number(this->tobj->gapList[i]->Trace2->GetId() );
				this->dtext+="\tgap size of:" + QString::number(this->tobj->gapList[i]->dist); 				
				this->dtext+="\tangle of" + QString::number(currentAngle*180/PI);	//this->tobj->gapList[i]->angle
				tobj->mergeTraces(this->tobj->gapList[i]->endPT1,this->tobj->gapList[i]->endPT2);
				++mergeCount;
	  		}	//end of if
			else if (this->tobj->gapList[i]->dist<= this->tobj->gapMax*(1+this->tobj->gapTol)&& (currentAngle < .3))
	  		{
				this->HighlightSelected(this->tobj->gapList[i]->Trace1, .125);
				this->HighlightSelected(this->tobj->gapList[i]->Trace2, .125);
				this->candidateGaps.push_back( this->tobj->gapList[i]);
				this->grayText+="\nAngle of:\t" + QString::number(currentAngle*180/PI);	//this->tobj->gapList[i]->angle
				this->grayText+="\tWith a distance of:\t" + QString::number(this->tobj->gapList[i]->dist);
				
			 } //end of else
		}//end of for merge
		myText+="\tNumber of merged lines:\t" + QString::number(mergeCount);
		if (this->candidateGaps.size()>=1)
		{
			this->poly_line_data->Modified();
			this->QVTK->GetRenderWindow()->Render();
			myText+="\nNumber of further possible lines:\t" + QString::number(this->candidateGaps.size());
			mergeAll = MergeInfo->addButton("Merge All", QMessageBox::YesRole);
			mergeNone = MergeInfo->addButton("Merge None", QMessageBox::NoRole);
			MergeInfo->setDetailedText(grayText);
		}		//end if graylist size
		else
		{
			this->Rerender();
			MergeInfo->setDetailedText(dtext);
		}		//end else graylist size
		MergeInfo->setText(myText);	
		MergeInfo->show();
		MergeInfo->exec();
		if(MergeInfo->clickedButton()==mergeAll)
		{
			for (j=0; j<this->candidateGaps.size();j++)
			{
				tobj->mergeTraces(this->candidateGaps[j]->endPT1,this->candidateGaps[j]->endPT2);
			}
			this->Rerender();
		}
		else if(MergeInfo->clickedButton()==mergeNone)
		{
			this->candidateGaps.clear();
		}
	}
	this->Rerender();
	this->tobj->gapList.clear();
	this->candidateGaps.clear();
	this->myText.clear();
	this->dtext.clear();
}
/*	other trace modifiers	*/
void View3D::SplitTraces()
{
  if(this->IDList.size()>=1)
    {
    std::unique(this->IDList.begin(), this->IDList.end());
    for(unsigned int i = 0; i < this->IDList.size(); i++)
      {
      this->tobj->splitTrace(this->IDList[i]);
      } 
    this->Rerender();
    }
  else
    {
    cout << "Nothing to split" << endl;
    }
}

void View3D::FlipTraces()
{
  if(this->IDList.size()>=1)
    {
    for(unsigned int i=0; i< this->IDList.size(); i++)
      {
      this->tobj->ReverseSegment(reinterpret_cast<TraceLine*>(
        this->tobj->hashc[this->IDList[i]]));
      }
    this->Rerender();
    }
}

void View3D::WriteToSWCFile()
{
  QString fileName = QFileDialog::getSaveFileName(
    this,
    tr("Save File"),
    "",
    tr("SWC Images (*.swc)"));
	this->tobj->WriteToSWCFile(fileName.toStdString().c_str());	
}


/*	Soma display stuff	*/
void View3D::ShowLoadSomaWindow()
{
	this->LoadSomaWidget->show();
	//this->OpenSomaButton->setDisabled(true);
}

void View3D::GetSomaPath()
{
	//Open file browser to select a file, takes file path and inserts it into the file field text box
	QString path;
	path = QFileDialog::getOpenFileName(this->LoadSomaWidget, "choose a file to load", QString::null, QString::null);
	this->SomaFileField->setText(path);
}
void View3D::GetSomaFile()
{

	//Checks that file has the the correct format
	if(this->somaRegex.exactMatch(this->SomaFileField->text()))
	{
		//Takes input from file field text box and casts it as an std string
		std::string fileName;
		std::string fileEx;
		fileName = this->SomaFileField->text().toStdString();

		//Takes captured file extension from the regex and stores it as an std string
		fileEx = this->somaRegex.cap(1).toStdString();
		for(unsigned int i = 0; i < fileEx.size(); i++)
		{
			fileEx[i] = tolower(fileEx[i]);
		}
			
		//Checks for proper file extension and accepts
		if(fileEx == "tiff" || fileEx == "tif")
		{
			//Checks that file exists
			if(this->CheckFileExists(fileName.c_str()) == 1)
			{
				this->SomaFile = fileName;
				this->viewSomas->setEnabled(true);
				this->LoadSomaWidget->hide();
				this->readImg(this->SomaFile);
				return;
			}
			else
			{
				std::cerr << "TIFF file doesn't exist." << std::endl;
				return;
			}
		}
		else
		{
			std::cerr << "Invalid File Type: Please load correct soma TIFF image." << std::endl;
			return;
		}
	}
	else
	{
		std::cerr << "Invalid File Type: Please load correct soma TIFF image." << std::endl;
	}
}

void View3D::HideLoadSomaWindow()
{
	this->LoadSomaWidget->hide();
}

void View3D::ShowLoadSeedWindow()
{
	this->LoadSeedWidget->show();
}

void View3D::GetSeedPath()
{
	//Gets seed file path from the file browser
	QString path;
	path = QFileDialog::getOpenFileName(this->LoadSeedWidget, "choose a file to load", QString::null, QString::null);
	this->SeedFileField->setText(path);
}

void View3D::GetSeedFile()
{
	//Checks that file has the the correct format
	if(this->seedRegex.exactMatch(this->SeedFileField->text()))
	{
		//Takes input from file field text box and casts it as an std string
		std::string fileName;
		std::string fileEx;
		fileName = this->SeedFileField->text().toStdString();

		//Takes captured file extension from the regex and stores it as an std string
		fileEx = this->seedRegex.cap(1).toStdString();
		for(unsigned int i = 0; i < fileEx.size(); i++)
		{
			fileEx[i] = tolower(fileEx[i]);
		}
			
		//Checks for proper file extension and accepts
		if(fileEx == "txt")
		{
			//Checks that file exists
			if(this->CheckFileExists(fileName.c_str()) == 1)
			{
				this->SeedFile = fileName;
				this->viewSeeds->setEnabled(true);
				this->LoadSeedWidget->hide();
				this->generateSeeds();
				return;
			}
			else
			{
				std::cerr << "TXT file doesn't exist." << std::endl;
				return;
			}
		}
		else
		{
			std::cerr << "Invalid File Type: Please load correct seed point TXT image." << std::endl;
			return;
		}
	}
	else
	{
		std::cerr << "Invalid File Type: Please load correct seed point TXT image." << std::endl;
	}
}

void View3D::HideLoadSeedWindow()
{
	this->LoadSeedWidget->hide();
}

void View3D::ShowSomaSettingsWindow()
{
	this->SomaOpacityField->setText(QString::number(this->somaopacity));
	this->SomaSettingsWidget->show();
}

void View3D::HideSomaSettingsWindow()
{
	this->SomaSettingsWidget->hide();
}

void View3D::ToggleSomas()
{
	//If view somas options is checked, render somas
	if(this->viewSomas->isChecked())
	{
		this->Renderer->AddActor(this->VolumeActor);
		this->QVTK->GetRenderWindow()->Render();
		std::cout << "Somas Rendered" << std::endl;
		return;
	}
	//If unchecked, remove them
	else
	{
		this->Renderer->RemoveActor(this->VolumeActor);
		this->QVTK->GetRenderWindow()->Render();
		std::cout << "Somas Removed" << std::endl;
		return;
	}
}

void View3D::ToggleSeeds()
{
	//If view seeds option is checked, render seeds
	if(this->viewSeeds->isChecked())
	{
		this->Renderer->AddActor(this->glyphActor);
		this->QVTK->GetRenderWindow()->Render();
		std::cout << "Seed Points Rendered" << std::endl;
		return;
	}
	//If unchecked, remove them
	else
	{
		this->Renderer->RemoveActor(this->glyphActor);
		this->QVTK->GetRenderWindow()->Render();
		std::cout << "Seeds Removed" << std::endl;
		return;
	}
}

void View3D::ApplySomaSettings()
{
	//If somas are currently on screen, change opacity and rerender
	if(this->viewSomas->isChecked())
	{
		
		this->somaopacity = this->SomaOpacityField->text().toDouble();
		this->VolumeActor->GetProperty()->SetOpacity(this->somaopacity);
		this->SomaSettingsWidget->hide();
		this->QVTK->GetRenderWindow()->Render();
		return;
	}
	//Else just change the soma opacity variable
	else
	{
		this->somaopacity = this->SomaOpacityField->text().toDouble();
		this->SomaSettingsWidget->hide();
	}
}

void View3D::AddContourThresholdSliders()
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
  sliderWidget2->SetInteractor(Interactor);
  sliderWidget2->SetRepresentation(sliderRep2);
  sliderWidget2->SetAnimationModeToAnimate();

  vtkSlider2DCallbackContourThreshold *callback_contour = vtkSlider2DCallbackContourThreshold::New();
  callback_contour->cfilter = this->ContourFilter;
  sliderWidget2->AddObserver(vtkCommand::InteractionEvent,callback_contour);
  sliderWidget2->EnabledOn();

}

void View3D::AddPlaybackWidget(char *filename)
{

  vtkSubRep *playbackrep = vtkSubRep::New();
  playbackrep->slice_counter=0;
  playbackrep->GetPositionCoordinate()->SetCoordinateSystemToNormalizedDisplay();
  playbackrep->GetPositionCoordinate()->SetValue(0.2,0.1);
  playbackrep->GetPosition2Coordinate()->SetCoordinateSystemToNormalizedDisplay();
  playbackrep->GetPosition2Coordinate()->SetValue(0.8,0.1);
  

  vtkSmartPointer<vtkPlaybackWidget> pbwidget = vtkSmartPointer<vtkPlaybackWidget>::New();
  pbwidget->SetRepresentation(playbackrep);
  pbwidget->SetInteractor(Interactor);

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
  for(unsigned int counter=0; counter<size[2]-1; counter++)
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
  colorTransferFunction->AddRGBPoint(20.0,1,0,0);

  
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
  Renderer->AddVolume(volume);
  printf("Added volume in playback widget");
  this->Volume = volume;
  playbackrep->vmapper=volumeMapper;
  /*playbackrep->im_pointer = &vtkimarray;*/
  playbackrep->Print(std::cout);
  
  this->QVTK->GetRenderWindow()->Render();
}
void View3D::AddVolumeSliders()
{
  vtkSliderRepresentation2D *sliderRep = vtkSliderRepresentation2D::New();
  sliderRep->SetValue(0.1);
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
  sliderWidget->SetInteractor(Interactor);
  sliderWidget->SetRepresentation(sliderRep);
  sliderWidget->SetAnimationModeToAnimate();

  vtkSlider2DCallbackBrightness *callback_brightness = vtkSlider2DCallbackBrightness::New();
  callback_brightness->volume = this->Volume;
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
  sliderWidget2->SetInteractor(Interactor);
  sliderWidget2->SetRepresentation(sliderRep2);
  sliderWidget2->SetAnimationModeToAnimate();

  vtkSlider2DCallbackContrast *callback_contrast = vtkSlider2DCallbackContrast::New();
  callback_contrast->volume = this->Volume;
  sliderWidget2->AddObserver(vtkCommand::InteractionEvent,callback_contrast);
  sliderWidget2->EnabledOn();

}
void View3D::readImg(std::string sourceFile)
{
  
  //Set up the itk image reader
  const unsigned int Dimension = 3;
  typedef unsigned char  PixelType;
  typedef itk::Image< PixelType, Dimension >   ImageType;
  typedef itk::ImageFileReader< ImageType >    ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( sourceFile );

  //Test opening and reading the input file
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
  //Itk image to vtkImageData connector
  typedef itk::ImageToVTKImageFilter<ImageType> ConnectorType;
  ConnectorType::Pointer connector= ConnectorType::New();
  connector->SetInput( reader->GetOutput() );

  //Route vtkImageData to the contour filter
  this->ContourFilter = vtkContourFilter::New();
  this->ContourFilter->SetInput(connector->GetOutput());
  this->ContourFilter->SetValue(0,10);            // this affects render

  this->ContourFilter->Update();

  //Route contour filter output to the mapper
  this->VolumeMapper = vtkPolyDataMapper::New();
  this->VolumeMapper->SetInput(this->ContourFilter->GetOutput());

  //Declare actor and set properties
  this->VolumeActor = vtkActor::New();
  this->VolumeActor->SetMapper(this->VolumeMapper);

  //this->VolumeActor->GetProperty()->SetRepresentationToWireframe();
  this->VolumeActor->GetProperty()->SetOpacity(this->somaopacity);
  this->VolumeActor->GetProperty()->SetColor(0.5,0.5,0.5);
  this->VolumeActor->SetPickable(0);

  std::cout << "Nuclei contour produced" << std::endl;
}

//Function to generate the seed point glyphs
void View3D::generateSeeds()
{
	//Call readSeeds function to process seedpoint file and return vector of points
	this->sdPts = readSeeds(this->SeedFile.c_str());

	//Set up float array
	this->pcoords = vtkSmartPointer<vtkFloatArray>::New();
	this->pcoords->SetNumberOfComponents(3);
	this->pcoords->SetNumberOfTuples(this->sdPts.size());

	//Read in vector of points into the float array
	for(unsigned int i = 0; i < this->sdPts.size(); i++)
	{
		float pts[3] = {this->sdPts[i].x, this->sdPts[i].y, this->sdPts[i].z};
		this->pcoords->SetTuple(i, pts);
	}

	//Set up vtkPoints from the float array
	this->points = vtkSmartPointer<vtkPoints>::New();
	this->points->SetData(this->pcoords);

	//Set up vtkPolyData from vtkPoints
	this->seedPoly = vtkSmartPointer<vtkPolyData>::New();
	this->seedPoly->SetPoints(this->points);

	//Create seed glyphs from sphere source, located at points in the polydata
	this->sphereSource = vtkSmartPointer<vtkSphereSource>::New();
	this->seedGlyph = vtkSmartPointer<vtkGlyph3D>::New();
	this->seedGlyph->SetInput(this->seedPoly);
	this->seedGlyph->SetSource(sphereSource->GetOutput());
	this->seedGlyph->SetVectorModeToUseNormal();
	this->seedGlyph->SetScaleModeToScaleByVector();
	this->seedGlyph->SetScaleFactor(20);
	this->seedGlyph->GeneratePointIdsOn();

	//Set up mapper and actor for the seeds
	this->glyphMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	this->glyphMapper->SetInput(this->seedGlyph->GetOutput());
	this->glyphActor = vtkSmartPointer<vtkLODActor>::New();
	this->glyphActor->SetMapper(this->glyphMapper);

	std::cout << "Seed Points Processed" << std::endl;
}

//Function to read in seed point coordinates from a .txt file
std::vector<point> View3D::readSeeds(std::string seedSource)
{
	point p;
	std::vector<point> pts;
	pts.clear();
	//int numLines = 0;
	
	//std::cout << seedSource << " filename" <<  std::endl;
	ifstream input(seedSource.c_str());

	if(!input)
	{
		std::cerr << "ERROR: Seed file can't be opened" << std::endl;
	}

	//Read in coordinates until the end of file, adding each set of point coordinates into the point vector
	while(true)
	{
		if(input.eof())
		{
			break;
		}	
		if(!(input >> p.x >> p.y >> p.z))
		{
			//std::cerr << "Invalid data point" << std::endl;
			input.clear();
			input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		}
		else
		{
			pts.push_back(p);
		}
	}
	return pts;
}

void View3D::rayCast(char *raySource)
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
  opacityTransferFunction->AddPoint(50,0.1);
 // opacityTransferFunction->AddPoint(40,0.1);
  // Create transfer mapping scalar value to color
  // Play around with the values in the following lines to better vizualize data
  vtkColorTransferFunction *colorTransferFunction = vtkColorTransferFunction::New();
    colorTransferFunction->AddRGBPoint(0.0, 0.0, 0.0, 0.0);
  colorTransferFunction->AddRGBPoint(50.0,1,0,0);

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
//  Renderer->AddVolume(volume);
  this->Volume = volume;
 // this->QVTK->GetRenderWindow()->Render();
  std::cout << "RayCast generated \n";
}

void View3D::closeEvent(QCloseEvent *event)
{
  event->accept();
}

bool View3D::CheckFileExists(const char *filename)
{
	FILE *fp = fopen(filename, "r");
	if(fp == NULL)
	{
		return false;
	}
	else
	{
		fclose(fp);
		return true;
	}
}
