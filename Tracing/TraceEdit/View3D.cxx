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
  v6: "ALISA" implemented
  v7: new GUI and file control
*/

#include "ftkGUI/PlotWindow.h"
//#include "ftkGUI/HistoWindow.h"
#include "ftkGUI/TableWindow.h"

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
#include "branchPT.h"
#include "TraceModel.h"
#include "MergeModel.h"
#include "ImageActors.h"
#include "TraceProjectManager.h"
#include "View3DHelperClasses.h"
#include "View3D.h"

View3D::View3D(int argc, char **argv)
{
	this->tobj = new TraceObject;
	int num_loaded = 0;
	this->Volume=0;
	bool tracesLoaded = false;
	this->translateImages = true;	//this is for testing a switch is needed
	this->Date.currentDate();
	this->Time.currentTime();
	this->Image.clear();
	this->TraceFiles.clear();
	this->SomaFile.clear();
	this->tempTraceFile.clear();
	this->ImageActors = new ImageRenderActors();
	this->EditLogDisplay = new QTextEdit();
	this->EditLogDisplay->setReadOnly(true);
	this->EditLogDisplay->setLineWrapMode(QTextEdit::NoWrap);
	this->EditLogDisplay->append("Farsight Trace Editor Started at: \nDate: \t" + this->Date.currentDate().toString( "ddd MMMM d yy" ) );
	this->EditLogDisplay->append("Time: \t" + this->Time.currentTime().toString( "h:m:s ap" ) );
	
  // load as many files as possible. Provide offset for differentiating types
  for(int counter=1; counter<argc; counter++)
    {
    int len = strlen(argv[counter]);
    if(strcmp(argv[counter]+len-3,"swc")==0)
      {
      printf("I detected swc\n");
      this->tobj->ReadFromSWCFile(argv[counter]);
	  tracesLoaded = true;
	  this->TraceFiles.append( QString(argv[counter]));
	  //this->TraceFiles
       }
    else if (strcmp(argv[counter]+len-3,"xml")==0)
      {
      printf("I detected xml\n");
      this->tobj->ReadFromRPIXMLFile(argv[counter]);
	  this->TraceFiles.append( QString(argv[counter]));
	  tracesLoaded = true;
      }
    else if (strcmp(argv[counter]+len-3,"vtk")==0)
      {
      printf("I detected vtk\n");
      this->tobj->ReadFromVTKFile(argv[counter]);
	  this->TraceFiles.append( QString(argv[counter]));
	  tracesLoaded = true;
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
	  this->Image.append( QString(argv[counter]));
	  this->ImageActors->loadImage(QString(argv[counter]).toStdString(), "Image");
      }
    num_loaded++;
    }
		this->CreateBootLoader();
	//if (num_loaded < 1)
	//{	
	//	this->hide();
	//	this->CreateBootLoader();
	//	return;
	//}//end load
	//else
	//{		
	//	this->show();
	//	if (tracesLoaded)
	//	{
	//		this->ShowTreeData();
	//	}
	//	this->statusBar()->showMessage(tr("Ready"));
	//}
}
View3D::View3D(TraceObject *Traces)
{
	this->tobj = Traces;
	this->Volume=0;
	this->Initialize();
	this->ShowTreeData();
	this->statusBar()->showMessage(tr("Trace Editor Started"));
}
void View3D::CreateBootLoader()
{
	// Create a window that allows files to be loaded
	this->bootLoadFiles = new QWidget;
	this->BootTrace = new QPushButton("Trace",this->bootLoadFiles);
	connect(this->BootTrace, SIGNAL(clicked()), this, SLOT(getTraceFile()));
	this->BootImage = new QPushButton("Image",this->bootLoadFiles);
	connect(this->BootImage, SIGNAL(clicked()), this, SLOT(getImageFile()));
	this->BootSoma = new QPushButton("Somas", this->bootLoadFiles);
	connect(this->BootSoma, SIGNAL(clicked()), this, SLOT(getSomaFile()));
	this->GetAUserName =  new QComboBox(this->bootLoadFiles);
	this->GetAUserName->setEditable(true);
	this->GetAUserName->setInsertPolicy(QComboBox::InsertAtCurrent);
	this->GetAUserName->addItems(this->TraceEditSettings.value("boot/userName", "default user").toStringList());
	this->GetLab = new QComboBox(this->bootLoadFiles);
	this->GetLab->setEditable(true);
	this->GetLab->setInsertPolicy(QComboBox::InsertAtCurrent);
	this->GetLab->addItems(this->TraceEditSettings.value("boot/LabName", "Roysam Lab").toStringList());
	this->scale = new QDoubleSpinBox(this->bootLoadFiles);
	this->scale->setValue(this->TraceEditSettings.value("boot/scale", 1).toDouble());
	this->scale->setSingleStep(.01);
	this->okBoot = new QPushButton("Start",this->bootLoadFiles);
	connect(this->okBoot, SIGNAL(clicked()), this, SLOT(OkToBoot()));
	this->Reload = new QPushButton("Reload", this->bootLoadFiles);
	connect(this->Reload, SIGNAL(clicked()), this, SLOT(ReloadState()));
	this->BootProject = new QPushButton("Project", this->bootLoadFiles);
	connect(this->BootProject, SIGNAL(clicked()), this, SLOT(LoadProject()));
	QFormLayout *LoadLayout = new QFormLayout(this->bootLoadFiles);
	LoadLayout->addRow(tr("User Name "), this->GetAUserName);
	LoadLayout->addRow(tr("Lab Name "), this->GetLab);
	LoadLayout->addRow(tr("Trace File"), this->BootTrace);
	LoadLayout->addRow(tr("Image File"), this->BootImage);
	LoadLayout->addRow(tr("Somas File"), this->BootSoma);
	LoadLayout->addRow(tr("uM Per Voxel"), this->scale);
	LoadLayout->addRow(tr("Run Trace Editor"), this->okBoot);
	LoadLayout->addRow(tr("Reload Previous Session"), this->Reload);
	LoadLayout->addRow(tr("Project"), this->BootProject);
	this->bootLoadFiles->show();
	this->bootLoadFiles->move(this->TraceEditSettings.value("boot/pos",QPoint(40, 59)).toPoint());
}
void View3D::ReloadState()
{
	int i;
	this->Image = this->TraceEditSettings.value("lastOpen/Image").toStringList();
	if(!this->Image.isEmpty()) 
	{
		for (i = 0; i < this->Image.size(); i ++)
		{
			this->ImageActors->loadImage(this->Image.at(i).toStdString(),"Image");
		}
	}
	
	this->tempTraceFile = this->TraceEditSettings.value("lastOpen/Temp").toStringList();
	if (!this->tempTraceFile.isEmpty())
	{	//opens last saved file 
		this->TraceFiles = this->tempTraceFile;
		// if saved file find log file
		QString logFileName = this->TraceEditSettings.value("lastOpen/Log").toString();
		if (!logFileName.isEmpty())
		{	//read and add old log file to new one
			QFile readLog(logFileName);
			if (readLog.open(QIODevice::ReadOnly | QIODevice::Text))
			{
				QTextStream in(&readLog);
				this->EditLogDisplay->append(in.readAll());
			}
		}//end of has log
	}
	else
	{	//if the last file was never saved open this
		this->TraceFiles = this->TraceEditSettings.value("lastOpen/Trace").toStringList();
	}
	if(!this->TraceFiles.isEmpty()) 
	{		
		for (i = 0; i < this->TraceFiles.size(); i++)
		{
			std::string traceFile = this->TraceFiles.at(i).toStdString();
			if(this->TraceFiles.at(i).endsWith("swc"))
			{
				this->tobj->ReadFromSWCFile((char*)traceFile.c_str());
			}
			else if(this->TraceFiles.at(i).endsWith("xml"))
			{
				this->tobj->ReadFromRPIXMLFile((char*)traceFile.c_str());
			}
			else if (this->TraceFiles.at(i).endsWith("vtk"))
			{
				this->tobj->ReadFromVTKFile((char*)traceFile.c_str());
			}
		}//end of sending trace file to proper reader
	}// end of trace files
	this->SomaFile = this->TraceEditSettings.value("lastOpen/Soma").toStringList();
	if(!this->SomaFile.isEmpty()) 
	{
		for (i = 0; i < this->SomaFile.size(); i++)
		{
			this->ImageActors->loadImage(this->SomaFile.at(i).toStdString(), "Soma");
		}
	}
	this->OkToBoot();
}
void View3D::OkToBoot()
{
	if(!this->TraceFiles.isEmpty() || !this->Image.isEmpty() || !this->SomaFile.isEmpty())
	{
		int i = 0;
		//this->show();
		this->uMperVoxel = this->scale->value();
		this->TraceEditSettings.setValue("boot/scale", this->uMperVoxel);
		this->Initialize();
		if (!this->TraceFiles.isEmpty() )
		{	
			this->EditLogDisplay->append("Trace file:");
			for (i = 0; i < this->TraceFiles.size(); i++)
			{
				this->EditLogDisplay->append("\t" + this->TraceFiles.at(i));
			}
			this->ShowTreeData();
		}
		if (!this->Image.isEmpty())
		{
			this->EditLogDisplay->append("Image file:");
			for (i=0; i < this->Image.size(); i++)
			{
				this->EditLogDisplay->append("\t" + this->Image.at(i));
			}
		}
		if (!this->SomaFile.isEmpty())
		{
			this->EditLogDisplay->append("Soma file: ");
			for (i=0; i < this->SomaFile.size(); i++)
			{
				this->EditLogDisplay->append("\t" + this->SomaFile.at(i));
			}
		}
		this->UserName = this->GetAUserName->currentText();
		QStringList allUsers = this->TraceEditSettings.value("boot/userName").toStringList();
		allUsers.removeAll(this->UserName);
		allUsers.prepend(this->UserName);
		this->TraceEditSettings.setValue("boot/userName", allUsers);
		this->EditLogDisplay->append("User: \t" + this->UserName);

		this->LabName = this->GetLab->currentText();
		QStringList allLabs = this->TraceEditSettings.value("boot/LabName").toStringList();
		allLabs.removeAll(this->LabName);
		allLabs.prepend(this->LabName);
		this->TraceEditSettings.setValue("boot/LabName", allLabs);
		this->EditLogDisplay->append("Lab: \t" + this->LabName);
		this->TraceEditSettings.setValue("boot/pos", this->bootLoadFiles->pos());
		this->bootLoadFiles->close();
		this->show();
//std::vector<std::string> LoadedImages = this->ImageActors->GetImageList();
//		for (unsigned int i =0; i < LoadedImages.size(); i++)
//		{
//			std::cout << LoadedImages.at(i)<< "\t";
//		}
	}
	else
	{
		QMessageBox *bootFailed = new QMessageBox;
		bootFailed->setText("There are no files to open. Please select a file to continue.");
		bootFailed->show();
	}
}
QString View3D::getSomaFile()
{
	QString somaFiles = QFileDialog::getOpenFileName(this , "Choose a Soma file to load", ".", 
	  tr("Image File ( *.tiff *.tif *.pic *.PIC ) "));
	if(!somaFiles.isEmpty())
	{
		this->SomaFile.append( somaFiles);
		this->ImageActors->loadImage(somaFiles.toStdString(), "Soma");
	}
	return somaFiles.section('/',-1);
}
QString View3D::getTraceFile()
{	std::string traceFile;
	QString trace = QFileDialog::getOpenFileName(this , "Load Trace Data", ".",
		tr(" TraceFile ( *.xml *.swc *.vtk " ));
	if (!trace.isEmpty())
	{
		this->TraceFiles.append( trace);
		traceFile = trace.toStdString();
		if(trace.endsWith("swc"))
		{
			this->tobj->ReadFromSWCFile((char*)traceFile.c_str());
		}
		else if(trace.endsWith("xml"))
		{
			this->tobj->ReadFromRPIXMLFile((char*)traceFile.c_str());
		}
		else if (trace.endsWith("vtk"))
		{
			this->tobj->ReadFromVTKFile((char*)traceFile.c_str());
		}
	}
	return trace.section('/',-1);
}
QString View3D::getImageFile()
{
	QString NewImageFile = QFileDialog::getOpenFileName(this , "Load Trace Image Data", ".",
		tr("Trace Image ( *.tiff *.tif *.pic *.PIC *.mhd" ));
	if (!NewImageFile.isEmpty())
	{
		this->Image.append( NewImageFile);
		int imgNum = this->ImageActors->loadImage(NewImageFile.toStdString(), "Image");
		if (this->translateImages)
		{
			//get x
			bool ok = false;
			double xx, yy, zz;
			while (!ok)
			{
				xx = QInputDialog::getDouble(this, tr("set x"), tr("Shift x by"), 0, -60000, 60000, 1, &ok);
			}
			ok = false;
			while (!ok)
			{
				yy = QInputDialog::getDouble(this, tr("set y"), tr("Shift y by"), 0, -60000, 60000, 1, &ok);
			}
			ok = false;
			while (!ok)
			{
				zz = QInputDialog::getDouble(this, tr("set z"), tr("Shift z by"), 0, -60000, 60000, 1, &ok);
			}
			this->ImageActors->ShiftImage(imgNum, xx, yy, zz );
		}
	}
	return NewImageFile.section('/',-1);
}
void View3D::LoadTraces()
{
	QString trace = this->getTraceFile();
	if (!trace.isEmpty())
	{
		this->statusBar()->showMessage(tr("Loading Trace") + trace);
		this->EditLogDisplay->append("Trace file: \t" + this->TraceFiles.last());
		if (this->tobj->FeatureHeaders.size() >=1)
		  {
			 this->TreeModel = new TraceModel(this->tobj->GetTraceLines(), this->tobj->FeatureHeaders);
		  }
		else
		  {
			 this->TreeModel->SetTraces(this->tobj->GetTraceLines());
		  }
		this->TreeModel->scaleFactor = this->uMperVoxel;
		this->ShowTreeData();
		this->UpdateLineActor();
		this->UpdateBranchActor();
		this->QVTK->GetRenderWindow()->Render();
		this->TreeModel->SetTraces(this->tobj->GetTraceLines());
	}	
	else
	{
		this->statusBar()->showMessage("Please select a Trace Data file");
	}
}

void View3D::LoadImageData()
{
	QString newImage = this->getImageFile();
	if (!newImage.isEmpty())
	{
		this->statusBar()->showMessage("Loading Image file" + newImage);
		this->EditLogDisplay->append("Image file: \t" + this->Image.last());
		
		this->Renderer->AddVolume(this->ImageActors->RayCastVolume(-1));
		this->QVTK->GetRenderWindow()->Render();
		//this->Rerender();
		this->statusBar()->showMessage("Image File Rendered");
	}
	else
	{
		this->statusBar()->showMessage("Please select an Image file");
	}
}
void View3D::LoadSomaFile()
{
	QString somaFile = this->getSomaFile();
	if(!somaFile.isNull())
	{
		this->Renderer->AddActor(this->ImageActors->ContourActor(-1));
		this->QVTK->GetRenderWindow()->Render();
		this->EditLogDisplay->append("Soma file: \t" + this->SomaFile.last());
		this->statusBar()->showMessage("Somas Rendered");
	}
}
void View3D::LoadProject()
{
	QString projectFile= QFileDialog::getOpenFileName(this , "Load Trace Project File", ".",
		tr("project ( *.xml" ));
	if (!projectFile.isEmpty())
	{
		ProjectManager * project = new ProjectManager((char*)projectFile.toStdString().c_str());
		for (unsigned int i = 0; i < project->size(); i++)
		{ 
			std::string FileName = project->GetFileName(i);
			QString type = QString(project->GetFileType(i).c_str());
			if ((type == "Image")||(type == "Soma"))
			{
				if (type == "Image")
				{
					this->Image.append(QString(FileName.c_str()));
		//this->EditLogDisplay->append("Image file: \t" + this->Image.last());
				}
				else {
					this->SomaFile.append(QString(FileName.c_str()));
		//this->EditLogDisplay->append("Soma file: \t" + this->SomaFile.last());
				}
				this->ImageActors->loadImage(FileName, project->GetFileType(i), 
					project->GetTranslationX(i),project->GetTranslationY(i),project->GetTranslationZ(i));

			}//end type image
			else if (type == "Trace")
			{
				this->tobj->SetTraceOffset(project->GetTranslationX(i),
					project->GetTranslationY(i),project->GetTranslationZ(i));
				QString trace = QString(FileName.c_str());
				if(trace.endsWith("swc"))
				{
					this->tobj->ReadFromSWCFile((char*)FileName.c_str());
				}
				else if(trace.endsWith("xml"))
				{
					this->tobj->ReadFromRPIXMLFile((char*)FileName.c_str());
				}
				else if (trace.endsWith("vtk"))
				{
					this->tobj->ReadFromVTKFile((char*)FileName.c_str());
				}
				this->TraceFiles.append(QString(FileName.c_str()));
		//this->EditLogDisplay->append("Trace file: \t" + this->TraceFiles.last());
			}//end type trace
			else if (type == "Log")
			{
			}//end type log
		}//end of for project size
		this->OkToBoot();
	}// end of project !empty
}
void View3D::SetImgInt()
{
	if (this->ImageActors->NumberOfImages()>=1)
	{	// image intensity values at each Trace Bit of trace line
		this->tobj->ImageIntensity(this->ImageActors->GetImageData(-1));
	}
}
void View3D::TraceBitImageIntensity(int ImgID)
{
	if (this->ImageActors->NumberOfImages()>=1)
	{	// image intensity values at each Trace Bit of trace line
		this->tobj->ImageIntensity(this->ImageActors->GetImageData(ImgID));
	}
}
View3D::~View3D()
{
  if(this->QVTK)
    {
    delete this->QVTK;
    }
  if(this->GapsPlotView)
    {
    delete this->GapsPlotView;
    }
  if(this->OpacitySlider)
    {
    this->OpacitySlider->Delete();
    }
  if(this->BrightnessSlider)
    {
    this->BrightnessSlider->Delete();
    }
  delete this->tobj;
  delete this->undoBuff;
  //the various Qt objects should be getting deleted by closeEvent and
  //parent/child relationships...
}

void View3D::Initialize()
{
	this->QVTK = 0;
	this->OpacitySlider = 0;
	this->BrightnessSlider = 0;
	this->tobj->gapTol = this->TraceEditSettings.value("mainWin/gapTol", .5).toDouble() ;
	this->tobj->gapMax = this->TraceEditSettings.value("mainWin/gapMax", 10).toInt();
	this->SmallLineLength = this->TraceEditSettings.value("mainWin/smallLine", 10).toInt();
	this->SelectColor =this->TraceEditSettings.value("mainWin/selectColor", .1).toDouble();
	this->lineWidth= this->TraceEditSettings.value("mainWin/LineWidth", 2).toDouble();
	this->GapsPlotView = NULL;
	this->TreePlot = NULL;
	this->FTKTable = NULL;
	this->GapsTableView = NULL;

	this->tobj->setSmallLineColor(.25);
	this->tobj->setMergeLineColor(.4);
	this->Ascending = Qt::AscendingOrder;

	this->undoBuff = new bufferType;
	this->numDeleted = 0;
	this->numMerged = 0;
	this->numSplit = 0;
	this->CreateGUIObjects();
	this->CreateLayout();
	this->CreateInteractorStyle();
	this->CreateActors();
	this->resize(this->TraceEditSettings.value("mainWin/size",QSize(850, 480)).toSize());
	this->move(this->TraceEditSettings.value("mainWin/pos",QPoint(40, 59)).toPoint());
	if (!this->TraceFiles.isEmpty())
	{
		this->setWindowTitle(tr("Trace Editor: ")+ this->TraceFiles.last().section('/',-1));
	}
	else
	{
		this->setWindowTitle(tr("Trace Editor"));
	}
	this->QVTK->GetRenderWindow()->Render();
	this->setupLinkedSpace();
}

void View3D::setupLinkedSpace()
{  
  this->tobj->Gaps.clear();
  this->MergeGaps = new MergeModel(this->tobj->Gaps);
  this->MergeGaps->setParent(this);
  if (this->tobj->FeatureHeaders.size() >=1)
  {
	  this->TreeModel = new TraceModel(this->tobj->GetTraceLines(), this->tobj->FeatureHeaders);
  }
  else
  {
	  this->TreeModel = new TraceModel(this->tobj->GetTraceLines());
  }
  this->TreeModel->setParent(this);
  this->connect(this->MergeGaps->GetObjectSelection(), SIGNAL(changed()), 
	  this,SLOT(updateSelectionHighlights()));
  this->connect(this->TreeModel->GetObjectSelection(), SIGNAL(changed()), 
	  this, SLOT(updateTraceSelectionHighlights()));

}

/*Set up the components of the interface */
void View3D::CreateGUIObjects()
{
  //Set up the main window's central widget
  this->CentralWidget = new QWidget(this);
  this->setCentralWidget(this->CentralWidget);

  //Set up a QVTK Widget for embedding a VTK render window in Qt.
  this->QVTK = new QVTKWidget(this->CentralWidget);
  this->Renderer = vtkSmartPointer<vtkRenderer>::New();
  this->QVTK->GetRenderWindow()->AddRenderer(this->Renderer);

  //Set up the menu bar
 
  this->saveAction = new QAction(tr("&Save as..."), this->CentralWidget);
    connect(this->saveAction, SIGNAL(triggered()), this, SLOT(SaveToFile()));
	this->saveAction->setStatusTip("Save results to file");
  this->exitAction = new QAction(tr("&Exit"), this->CentralWidget);
	connect(this->exitAction, SIGNAL(triggered()), this, SLOT(close()));
	this->exitAction->setStatusTip("Exit the Trace Editor");
  this->loadTraceAction = new QAction("Load Trace", this->CentralWidget);
    connect(this->loadTraceAction, SIGNAL(triggered()), this, SLOT(LoadTraces()));
    this->loadTraceAction->setStatusTip("Load traces from .xml or .swc file");
  this->loadTraceImage = new QAction("Load Image", this->CentralWidget);
	connect (this->loadTraceImage, SIGNAL(triggered()), this, SLOT(LoadImageData()));
	this->loadTraceImage->setStatusTip("Load an Image to RayCast Rendering");
//Loading soma data
  this->loadSoma = new QAction("Load Somas", this->CentralWidget);
   connect(this->loadSoma, SIGNAL(triggered()), this, SLOT(LoadSomaFile()));

 //Set up the buttons that the user will use to interact with this program. 
  this->ListButton = new QAction("List", this->CentralWidget);
	connect(this->ListButton, SIGNAL(triggered()), this, SLOT(ListSelections()));
	this->ListButton->setStatusTip("List all selections");
  this->ClearButton = new QAction("Clear", this->CentralWidget); 
	connect(this->ClearButton, SIGNAL(triggered()), this, SLOT(ClearSelection()));
	this->ClearButton->setStatusTip("Clear all selections");
  this->SelectTreeAction = new QAction("Select Tree", this->CentralWidget); 
	connect(this->SelectTreeAction, SIGNAL(triggered()), this, SLOT(SelectTrees()));
	this->SelectTreeAction->setStatusTip("Select the entire tree");
  this->DeleteButton = new QAction("Delete", this->CentralWidget);
	connect(this->DeleteButton, SIGNAL(triggered()), this, SLOT(DeleteTraces()));
	this->DeleteButton->setStatusTip("Delete all selected traces");
  this->MergeButton = new QAction("Merge", this->CentralWidget);
	connect(this->MergeButton, SIGNAL(triggered()), this, SLOT(MergeTraces()));
	this->MergeButton->setStatusTip("Start Merge on selected traces");
  this->BranchButton = new QAction("Branch", this->CentralWidget);
	connect(this->BranchButton, SIGNAL(triggered()), this, SLOT(AddNewBranches()));
	this->BranchButton->setStatusTip("Add branches to trunk");
  this->SplitButton = new QAction("Split", this->CentralWidget); 
	connect(this->SplitButton, SIGNAL(triggered()), this, SLOT(SplitTraces()));
	this->SplitButton->setStatusTip("Split traces at point where selected");
  this->FlipButton = new QAction("Flip", this->CentralWidget);
	connect(this->FlipButton, SIGNAL(triggered()), this, SLOT(FlipTraces()));
	this->FlipButton->setStatusTip("Flip trace direction");
  this->SettingsButton = new QAction("Settings", this->CentralWidget);
	connect(this->SettingsButton, SIGNAL(triggered()), this,
		SLOT(ShowSettingsWindow()));
	this->SettingsButton->setStatusTip("edit the display and tolerance settings");
  this->AutomateButton = new QAction("Small Lines", this->CentralWidget);
	connect(this->AutomateButton, SIGNAL(triggered()), this, SLOT(SLine()));
	this->AutomateButton->setStatusTip("Automatic selection of all small lines");
  this->root = new QAction("Set Root", this->CentralWidget);
	connect(this->root, SIGNAL(triggered()), this, SLOT(SetRoots()));
	this->root->setStatusTip("Solve Branch order by defining Root Trace Lines");
  this->explodeTree = new QAction("Break", this->CentralWidget);
  connect(this->explodeTree, SIGNAL(triggered()), this, SLOT( ExplodeTree()));
  this->explodeTree->setStatusTip("Break tree into segments,aka Explode. Tree can be rebuilt using set root");

  this->UndoButton = new QAction("&Undo", this->CentralWidget);  
	connect(this->UndoButton, SIGNAL(triggered()), this, SLOT(UndoAction()));
  this->RedoButton = new QAction("&Redo", this->CentralWidget);
	connect(this->RedoButton, SIGNAL(triggered()), this, SLOT(RedoAction()));
	this->ImageIntensity = new QAction("Intensity", this->CentralWidget);
	connect(this->ImageIntensity, SIGNAL(triggered()), this, SLOT(SetImgInt()));
  //Setup the tolerance settings editing window
  this->SettingsWidget = new QWidget();
  QIntValidator *intValidator = new QIntValidator(1, 100, this->SettingsWidget);
  //QDoubleValidator *doubleValidator =
  //  new QDoubleValidator(0, 100, 2, this->SettingsWidget);
  this->MaxGapField = new QLineEdit(this->SettingsWidget);
  this->MaxGapField->setValidator(intValidator);
  this->GapToleranceField = new QLineEdit(this->SettingsWidget);
  this->GapToleranceField->setValidator(intValidator);
  this->LineLengthField = new QSpinBox(this->SettingsWidget);
  this->LineLengthField->setRange(0,100);
  this->ColorValueField = new QSpinBox(this->SettingsWidget);
  this->ColorValueField->setRange(0,100);
  this->LineWidthField = new QSpinBox(this->SettingsWidget);
  this->LineWidthField->setRange(1,5);
  this->ApplySettingsButton = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
  //this->CancelSettingsButton = new QPushButton("&Cancel", this->SettingsWidget);
  connect(this->ApplySettingsButton, SIGNAL(accepted()), this, SLOT(ApplyNewSettings()));
  connect(this->ApplySettingsButton, SIGNAL(rejected()), this, SLOT(HideSettingsWindow()));

	QStringList types;
	types <<"0 = undefined" << "1 = soma" <<"2 = axon" <<"3 = dendrite" 
		<<"4 = apical dendrite" <<"5 = fork point" <<"6 = end point" <<"7 = custom";
	this->typeCombo = new QComboBox;
	this->typeCombo->addItems(types);
	connect(this->typeCombo, SIGNAL(activated( int )), this, SLOT(SetTraceType(int )));
	
// Lables for the status bar to show edit counts
	this->SplitLabel = new QLabel(this);
	this->SplitLabel->setText(QString::number(this->numSplit));
	this->MergeLabel = new QLabel(this);
	this->MergeLabel->setText(QString::number(this->numMerged));
	this->DeleteLabel = new QLabel(this);
	this->DeleteLabel->setText(QString::number(this->numDeleted));
}

void View3D::CreateLayout()
{
	this->fileMenu = this->menuBar()->addMenu(tr("&File"));
	this->fileMenu->addAction(this->loadTraceAction);
	this->fileMenu->addAction(this->loadTraceImage);
	this->fileMenu->addAction(this->loadSoma);
	this->fileMenu->addAction(this->saveAction);
	this->fileMenu->addAction(this->exitAction);

	//this->ShowToolBars = this->menuBar()->addMenu(tr("Tool Bars"));

  this->EditsToolBar = addToolBar(tr("Edit Toolbar"));
  this->EditsToolBar->setToolTip("EditToolBar");
  this->menuBar()->addAction(this->EditsToolBar->toggleViewAction());
  //this->ShowToolBars->addAction(this->EditsToolBar->toggleViewAction());
  /*this->EditsToolBar->addAction(this->saveAction);
  this->EditsToolBar->addAction(this->exitAction);
  this->EditsToolBar->addSeparator();*/
  this->EditsToolBar->addAction(this->AutomateButton);
  this->EditsToolBar->addAction(this->ListButton);
  this->EditsToolBar->addAction(this->ClearButton);
  this->EditsToolBar->addAction(this->SelectTreeAction);
  this->EditsToolBar->addSeparator();
  this->EditsToolBar->addAction(this->DeleteButton);
  this->EditsToolBar->addAction(this->MergeButton);
  this->EditsToolBar->addAction(this->SplitButton);
  this->EditsToolBar->addAction(this->FlipButton);
  this->EditsToolBar->addWidget(this->typeCombo);
  this->EditsToolBar->addSeparator();
  this->EditsToolBar->addAction(this->loadSoma);
  this->EditsToolBar->addAction(this->SettingsButton);

  this->BranchToolBar = addToolBar(tr("Branch Toolbar"));
  this->BranchToolBar->setToolTip("Branch Toolbar");
  this->menuBar()->addAction(this->BranchToolBar->toggleViewAction());
  //this->ShowToolBars->addAction(this->BranchToolBar->toggleViewAction());
  this->BranchToolBar->addAction(this->explodeTree);
  this->BranchToolBar->addAction(this->BranchButton);
  this->BranchToolBar->addAction(this->root);
  this->BranchToolBar->addAction(this->ImageIntensity);

  QGridLayout *viewerLayout = new QGridLayout(this->CentralWidget);
  viewerLayout->addWidget(this->QVTK, 0, 0);
  //may add a tree view here
   //layout for the settings window
  QFormLayout *settingsLayout = new QFormLayout(this->SettingsWidget);
  settingsLayout->addRow(tr("Maximum gap length:"), this->MaxGapField);
  settingsLayout->addRow(tr("Gap length tolerance:"),this->GapToleranceField);
  settingsLayout->addRow(tr("Small line length:"),this->LineLengthField);
  settingsLayout->addRow(tr("Color value RGB scalar 0 to 1:"),this->ColorValueField);
  settingsLayout->addRow(tr("Line width:"),this->LineWidthField);
  settingsLayout->addRow(this->ApplySettingsButton);
 

  this->statusBar()->addPermanentWidget(new QLabel("Statistics: Split: ", this));
  this->statusBar()->addPermanentWidget(this->SplitLabel,0);
  this->statusBar()->addPermanentWidget(new QLabel(" Merged: ", this));
  this->statusBar()->addPermanentWidget(this->MergeLabel,0);
  this->statusBar()->addPermanentWidget(new QLabel(" Deleted: ", this));
  this->statusBar()->addPermanentWidget(this->DeleteLabel,0);

  this->InformationDisplays = new QDockWidget("Edit Log Information", this);
  this->InformationDisplays->setWidget(this->EditLogDisplay);
  this->addDockWidget(Qt::LeftDockWidgetArea, this->InformationDisplays);
  this->menuBar()->addAction(this->InformationDisplays->toggleViewAction());
  this->InformationDisplays->hide();
  //this->ShowToolBars->addAction(this->InformationDisplays->toggleViewAction());
}

void View3D::CreateInteractorStyle()
{
  this->Interactor = this->QVTK->GetRenderWindow()->GetInteractor();
  //keep mouse command observers, but change the key ones
  this->keyPress = vtkSmartPointer<vtkCallbackCommand>::New();
  this->keyPress->SetCallback(HandleKeyPress);
  this->keyPress->SetClientData(this);

  this->Interactor->RemoveObservers(vtkCommand::KeyPressEvent);
  this->Interactor->RemoveObservers(vtkCommand::KeyReleaseEvent);
  this->Interactor->RemoveObservers(vtkCommand::CharEvent);
  this->Interactor->AddObserver(vtkCommand::KeyPressEvent, this->keyPress);

  //use trackball control for mouse commands
  vtkSmartPointer<vtkInteractorStyleTrackballCamera> style =
    vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
  this->Interactor->SetInteractorStyle(style);
  this->CellPicker = vtkSmartPointer<vtkCellPicker>::New();
  this->CellPicker->SetTolerance(0.004);
  this->Interactor->SetPicker(this->CellPicker);
  this->isPicked = vtkSmartPointer<vtkCallbackCommand>::New();
  this->isPicked->SetCallback(PickCell);

  //isPicked caller allows observer to intepret click 
  this->isPicked->SetClientData(this);            
  this->Interactor->AddObserver(vtkCommand::RightButtonPressEvent,isPicked);
}

void View3D::CreateActors()
{
  this->LineActor = vtkSmartPointer<vtkActor>::New();
  this->LineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  this->UpdateLineActor();
  this->LineActor->SetPickable(1);
  this->Renderer->AddActor(this->LineActor);

  this->UpdateBranchActor();
  this->Renderer->AddActor(this->BranchActor);
  for (unsigned int i = 0; i < this->ImageActors->NumberOfImages(); i++)
  {
	  if (this->ImageActors->isRayCast(i))
	  {
		  this->Renderer->AddVolume(this->ImageActors->RayCastVolume(i));
		  // test image access functions
		 /* if (!this->TraceFiles.isEmpty())
		  {
			  this->tobj->ImageIntensity(this->ImageActors->GetImageData(i));
		  }*/
		  /*std::vector<double>  dimensions = this->ImageActors->GetImageSize(i);
		  std::cout<< "\nTest of image accessors \nx\t"<< dimensions[0] << "\ty\t"<< dimensions[1]<< "\tz\t"<< dimensions[2]<< "\n";
		  double intensity = this->ImageActors->pointData(i, 646,215,5);
		  std::cout<< "x:646 y:215 z:5\t" << intensity<< "\n";*/
		  //These functions are not required
	  }
	  else
	  {
		  this->Renderer->AddActor(this->ImageActors->ContourActor(i));
	  }
  }
  this->QVTK->GetRenderWindow()->Render();
  //sphere is used to mark the picks
  this->CreateSphereActor();
  Renderer->AddActor(this->SphereActor);
}

void View3D::CreateSphereActor()
{
  this->Sphere = vtkSmartPointer<vtkSphereSource>::New();
  this->Sphere->SetRadius(3);
  this->SphereMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  this->SphereMapper->SetInput(this->Sphere->GetOutput());
  this->SphereMapper->GlobalImmediateModeRenderingOn();

  this->SphereActor = vtkSmartPointer<vtkActor>::New();
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
  this->LineLengthField->setValue(this->SmallLineLength);
  this->ColorValueField->setValue(this->SelectColor*100);
  this->LineWidthField->setValue(this->lineWidth);
  this->SettingsWidget->show();
}

void View3D::ApplyNewSettings()
{
  this->tobj->gapMax = this->MaxGapField->text().toInt();
  this->tobj->gapTol = this->GapToleranceField->text().toDouble();
  this->SmallLineLength = this->LineLengthField->text().toFloat();
  this->SelectColor = this->ColorValueField->text().toDouble()/100;
  this->lineWidth = this->LineWidthField->text().toFloat();
  //this->Rerender();
  this->SettingsWidget->hide();
  this->statusBar()->showMessage(tr("Applying new settings"),3000);
}

void View3D::HideSettingsWindow()
{
  this->SettingsWidget->hide();
}
/*  picking */
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
    view->SelectedTraceIDs.push_back(cell_id);
    TraceLine *tline = reinterpret_cast<TraceLine*>(view->tobj->hashc[cell_id]);
	if (view->tobj->Gaps.size() < 1)
	{
		view->TreeModel->SelectByIDs(tline->GetId());
		//view->HighlightSelected(tline, view->SelectColor);
		//tline->Getstats();              //prints the id and end coordinates to the command prompt 
		view->SphereActor->SetPosition(pickPos);    //sets the selector to new point
		view->SphereActor->VisibilityOn();      //deleteTrace can turn it off 
		view->poly_line_data->Modified();
	}
	else
	{ //int size = view->tobj->Gaps.size();
		int id = tline->GetId();
		view->MergeGaps->SelectbyTraceID(id);
		view->poly_line_data->Modified();
	}
    //update the head Qt view here too...
  }// end if pick
  view->QVTK->GetRenderWindow()->Render();             //update the render window
}

void View3D::updateTraceSelectionHighlights()
{
	this->UpdateLineActor();
	std::vector<TraceLine*> Selections = this->TreeModel->GetSelectedTraces();
	for (unsigned int i = 0; i < Selections.size(); i++)
	{
		this->HighlightSelected(Selections[i],this->SelectColor);
	}
	this->poly_line_data->Modified();
	this->QVTK->GetRenderWindow()->Render();
	this->statusBar()->showMessage(tr("Selected\t")
		+ QString::number(Selections.size()) +tr("\ttraces"));
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
  this->statusBar()->showMessage(tr("Rerender Image"));
  this->SphereActor->VisibilityOff();
  this->SelectedTraceIDs.clear();
  /*this->MergeGaps->GetSelectionModel()->clearSelection();*/
  this->Renderer->RemoveActor(this->BranchActor);
  this->UpdateLineActor();
  this->UpdateBranchActor();
  this->Renderer->AddActor(this->BranchActor); 
  this->TreeModel->SetTraces(this->tobj->GetTraceLines()); 
  this->QVTK->GetRenderWindow()->Render();
  this->FTKTable->setModels(this->TreeModel->getDataTable(), this->TreeModel->GetObjectSelection());
  this->FTKTable->update();
  this->TreePlot->setModels(this->TreeModel->getDataTable(), this->TreeModel->GetObjectSelection());
  this->TreePlot->update();
	this->SplitLabel->setText(QString::number(this->numSplit));
	this->MergeLabel->setText(QString::number(this->numMerged));
	this->DeleteLabel->setText(QString::number(this->numDeleted));
  this->statusBar()->showMessage(tr("Finished Rerendering Image"));
}

void View3D::UpdateLineActor()
{
	
  this->poly_line_data = this->tobj->GetVTKPolyData();
  this->poly_line_data->Modified();
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
    //view->tobj->Gaps.clear();
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
      view->SaveToFile();
    break;

    case 't':
      view->ShowSettingsWindow();
      break;

    case 'a':
      view->SLine();
      break;

    case 'q':
      view->updateSelectionHighlights();
      break;

    //case '-':
    //  if(view->SelectedTraceIDs.size()>=1)
    //    {
    //    TraceLine* tline=reinterpret_cast<TraceLine*>(
    //      view->tobj->hashc[view->SelectedTraceIDs[view->SelectedTraceIDs.size()-1]]);
    //    view->HighlightSelected(tline, tline->getTraceColor()-.25);
    //    view->SelectedTraceIDs.pop_back();
    //    cout<< " These lines: ";
    //    for (unsigned int i = 0; i < view->SelectedTraceIDs.size(); i++)
    //      {
    //      cout<<  "\t"<<view->SelectedTraceIDs[i];   
    //      } 
    //    cout<< " \t are selected \n" ;   
    //    }
    //  break;

	case 'b':
		view->AddNewBranches();
		break;

    default:
      break;
    }
}
/*  Actions   */
void View3D::SLine()
{
  int numLines;
  this->tobj->FindMinLines(this->SmallLineLength);
  numLines= this->tobj->SmallLines.size();
  this->TreeModel->SelectByIDs(this->tobj->SmallLines);
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
	this->DeleteTraces();
    this->tobj->SmallLines.clear();
  }
  break;
  case QMessageBox::No:
   {
     this->tobj->SmallLines.clear();
	 this->TreeModel->GetObjectSelection()->clear();
   }
   break;
  }
  this->Rerender();
}


void View3D::ListSelections()
{
	std::vector<TraceLine* > IDs = this->TreeModel->GetSelectedTraces();
	QMessageBox *selectionInfo = new QMessageBox;
	QString listText;
	QString selectedText;
	if(this->tobj->BranchPoints.size()> 0)
	{
		listText += QString::number(this->tobj->BranchPoints.size())
			+ " Branches are unsolved/n";
	}
	else if (IDs.size()<= 0)
	{
	listText = tr("No traces selected");
	}
	else
	{
	listText += QString::number(IDs.size()) + " lines are selected\n";
	for (unsigned int i = 0; i < IDs.size(); i++)
	  {
		  selectedText += QString(IDs[i]->stats().c_str()) + "\n";   
	  } 
	}
	this->statusBar()->showMessage(listText);
	selectionInfo->setText(listText);
	selectionInfo->setDetailedText(selectedText);
	selectionInfo->show();
}
void View3D::ShowTreeData()
{
	if (this->FTKTable)
	{
		this->FTKTable->close();
	}
	if (this->TreePlot)
	{
		this->TreePlot->close();
	}
	this->FTKTable = new TableWindow();
	this->FTKTable->setModels(this->TreeModel->getDataTable(), this->TreeModel->GetObjectSelection());
	this->FTKTable->setWindowTitle("Trace Object Features Table");
	this->FTKTable->move(this->TraceEditSettings.value("TraceTable/pos",QPoint(32, 561)).toPoint());
	this->FTKTable->resize(this->TraceEditSettings.value("TraceTable/size",QSize(600, 480)).toSize());
	this->FTKTable->show();

	this->TreePlot = new PlotWindow();
	this->TreePlot->setModels(this->TreeModel->getDataTable(), this->TreeModel->GetObjectSelection());
	this->TreePlot->setWindowTitle("Trace Object Features Plot");
	this->TreePlot->move(this->TraceEditSettings.value("TracePlot/pos",QPoint(890, 59)).toPoint());
	this->TreePlot->show();
}

void View3D::ClearSelection()
{
	if(this->GapsPlotView)
	{
	  this->GapsPlotView->close();
	}
	if (this->GapsTableView)
	{
	  this->GapsTableView->close();
	}
	this->tobj->Gaps.clear();
	this->candidateGaps.clear();	
	this->tobj->BranchPoints.clear();
	this->myText.clear();
	this->dtext.clear();
	this->Rerender();
	this->statusBar()->showMessage("All Clear", 4000);

	cout << this->TreePlot->pos().x() << ", " << this->TreePlot->pos().y() << endl;
}

void View3D::SelectTrees()
{
	std::vector<TraceLine*> roots = this->TreeModel->getRoots();
	if (roots.size() > 0)
	{
		//this->ClearSelection();
		std::vector<int> ids;
		ids = this->tobj->GetTreeIDs(roots);
		this->TreeModel->SelectByIDs(ids);
	}//end root size
}
/*  delete traces functions */
void View3D::DeleteTraces()
{
	unsigned int i;
	this->statusBar()->showMessage(tr("Deleting"));
	std::vector<TraceLine*> traceList = this->TreeModel->GetSelectedTraces();
	if (traceList.size() >=1)
	{
		this->EditLogDisplay->append(tr("Deleted\t") + QString::number(traceList.size()) + tr("\ttraces"));
		this->EditLogDisplay->append( "\tID\tType\tSize\tLength\tEuclidian Length\tRadii\tFragmentation Smoothness\tParent ID");
		for (i = 0; i < traceList.size(); i++)
		{				
			this->EditLogDisplay->append( QString(traceList[i]->stats().c_str()));
			this->DeleteTrace(traceList[i]); 
		}
		this->numDeleted += (int) traceList.size();
		this->ClearSelection();
		this->statusBar()->showMessage(tr("Deleted\t") + QString::number(traceList.size()) + tr("\ttraces"));
	}
	else
	{
		this->statusBar()->showMessage(tr("Nothing to Delete \n"));
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

/*	branching functions	*/
void View3D::SetRoots()
{
	std::vector<int> ids = this->TreeModel->GetSelecectedIDs();
	this->EditLogDisplay->append("Setting roots");
	int numToSolve= this->tobj->solveParents(ids);
	this->tobj->cleanTree();
	this->Rerender();
	this->TreeModel->SetTraces(this->tobj->GetTraceLines());
	this->statusBar()->showMessage(QString::number(numToSolve)+ " Remaining Branches");
  //this->Rerender();
}
void View3D::AddNewBranches()
{
	TraceLine* trunk;
	std::vector <TraceLine*> newChildren;
	std::vector <TraceLine*> selected = this->TreeModel->GetSelectedTraces();
	if (selected.size() > 1)
	{
		trunk = selected[0];
		if (trunk->Orient(selected[2]))
		{
			this->tobj->ReverseSegment(trunk);
		}
		for (unsigned int i = 1; i < selected.size(); i ++)
		{
			if (!selected[i]->Orient(trunk))
			{
				this->tobj->ReverseSegment(selected[i]);
			}
			newChildren.push_back(selected[i]);
		}
		this->AddChildren(trunk, newChildren);
		this->ClearSelection();
		this->statusBar()->showMessage(tr("Update Tree Plots"));
		this->TreeModel->SetTraces(this->tobj->GetTraceLines());
		this->statusBar()->showMessage(tr("Branching complete"));
	}
}
void View3D::ExplodeTree()
{
	this->tobj->BranchPoints.clear();
	std::vector<TraceLine*> roots = this->TreeModel->getRoots();
	for (unsigned int i = 0; i < roots.size(); i++)
	{
		this->tobj->explode(roots.at(i));
	}
	//this->tobj->cleanTree();
	this->Rerender();
	this->TreeModel->SetTraces(this->tobj->GetTraceLines());
}
void View3D::AddChildren(TraceLine *trunk, std::vector<TraceLine*> childTraces)
{
	for (unsigned int i = 0; i < childTraces.size(); i++)
	{
		if (childTraces[i]->GetParentID() == -1)
		{
			childTraces[i]->SetParent(trunk);
			trunk->AddBranch(childTraces[i]);
		}//end check/make connection 
	}//end child trace size loop
}
/*  merging functions */
void View3D::MergeTraces()
{
  if (this->tobj->Gaps.size() > 1)
    {
    if(this->GapsPlotView)
      {
      this->GapsPlotView->close();
      }
    if (this->GapsTableView)
      {
      this->GapsTableView->close();
      }
    this->MergeSelectedTraces();  	
    }
  else
    {
    int conflict = 0;
	std::vector<TraceLine*> traceList = this->TreeModel->GetSelectedTraces();
	if (traceList.size() >=1)
	{
		conflict = this->tobj->createGapLists(traceList);
     }
    else
    {
		conflict = this->tobj->createGapLists(this->tobj->GetTraceLines());
    }       
	if(conflict == -1)
    {
        //user aborted
        this->tobj->Gaps.clear();
        return;
    }
    unsigned int i; 
    QMessageBox MergeInfo;
    if (this->tobj->Gaps.size() > 1)
    {
      for (i=0;i<this->tobj->Gaps.size(); i++)
      {
        this->tobj->Gaps[i]->compID = i;
        //this->HighlightSelected(this->tobj->Gaps[i]->Trace1, .25);
        //this->HighlightSelected(this->tobj->Gaps[i]->Trace2, .25);
      }
      MergeInfo.setText("\nNumber of computed distances:\t" 
        + QString::number(this->tobj->Gaps.size())
        +"\nConflicts resolved:\t" + QString::number(conflict)
        +"\nEdit selection or press merge again");
      this->ShowMergeStats();
    }//end if this->tobj->Gaps size > 1
    else 
      {
      if (this->tobj->Gaps.size() ==1)
        {   
        tobj->mergeTraces(this->tobj->Gaps[0]->endPT1,this->tobj->Gaps[0]->endPT2);
		this->EditLogDisplay->append("Merged Trace:\t"  + QString(this->tobj->Gaps[0]->stats().c_str()));/*
			+ QString::number(this->tobj->Gaps[0]->Trace1->GetId()) + "\tto\t" 
			+ QString::number(this->tobj->Gaps[0]->Trace2->GetId()));*/
		this->numMerged++;
		this->ClearSelection();
        MergeInfo.setText(this->myText + "\nOne Trace merged");
		this->statusBar()->showMessage(tr("Update Tree Plots"));
		this->TreeModel->SetTraces(this->tobj->GetTraceLines());
		this->statusBar()->showMessage(tr("Done With Merge"));
        }
      else
        {
        this->Rerender();
        MergeInfo.setText("\nNo merges possible, set higher tolerances\n"); 
        } 
      }   
    MergeInfo.exec();
    this->myText.clear();
    this->poly_line_data->Modified();
    this->QVTK->GetRenderWindow()->Render();
    }//end else size
}

void View3D::ShowMergeStats()
{
  this->MergeGaps->SetTraceGaps(this->tobj->Gaps);
  this->GapsTableView = new TableWindow(); 
  this->GapsTableView->setModels(this->MergeGaps->getDataTable(), this->MergeGaps->GetObjectSelection());
  this->GapsTableView->setWindowTitle("Computed Features for Merge");
  this->GapsTableView->show();
  this->GapsPlotView = new PlotWindow();
  this->GapsPlotView->setModels(this->MergeGaps->getDataTable(), this->MergeGaps->GetObjectSelection());
  this->GapsPlotView->setWindowTitle("Computed Features for Merge");
  this->GapsPlotView->show();
}

void View3D::updateSelectionHighlights()
{
  bool selected = false;
  int curID;
  std::vector<int> GapIDs = this->MergeGaps->GetSelectedGapIDs();
  for (unsigned int i = 0; i < this->tobj->Gaps.size(); i++)
  {
    curID = this->tobj->Gaps[i]->compID;
    selected = false;
    unsigned int j =0;
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
      this->HighlightSelected(this->tobj->Gaps[i]->Trace1, .05);
      this->HighlightSelected(this->tobj->Gaps[i]->Trace2, .05);
    }//end true
    else
    {
      this->HighlightSelected(this->tobj->Gaps[i]->Trace1, .3);
      this->HighlightSelected(this->tobj->Gaps[i]->Trace2, .3);
    }//end else true
  }//end for i
  this->poly_line_data->Modified();
  this->QVTK->GetRenderWindow()->Render();
  this->statusBar()->showMessage(tr("Done"));
}
void View3D::MergeSelectedTraces()
{
	this->statusBar()->showMessage(tr("Merging"));
  this->updateSelectionHighlights(); //if it didnt
  std::vector<int> GapIDs = this->MergeGaps->GetSelectedGapIDs();
  bool selected = false;
  int curID;
  QMessageBox MergeInfo;
  double aveCost=0;
  QPushButton *mergeAll;
  QPushButton *mergeNone;
  unsigned int i=0, j=0;
  MergeInfo.setText("Merge Function");
  if (GapIDs.size() > 1)
  {
    for (i = 0; i < this->tobj->Gaps.size(); i++)
    {
      curID = this->tobj->Gaps[i]->compID;
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
		  aveCost += this->tobj->Gaps[i]->cost;
        this->tobj->mergeTraces(this->tobj->Gaps[i]->endPT1,this->tobj->Gaps[i]->endPT2);
		this->dtext += this->tobj->Gaps[i]->stats().c_str();
      }
    } 
    MergeInfo.setText("merged " + QString::number(GapIDs.size()) + " traces.");  
	this->EditLogDisplay->append("merged " + QString::number(GapIDs.size()) + " traces.");  
	this->numMerged += (int)GapIDs.size();
	this->dtext+="\nAverage Cost of Merge:\t" 
		+ QString::number( aveCost / GapIDs.size() );
	this->dtext+= "\nSum of Merge costs " + QString::number(aveCost);
	this->EditLogDisplay->append(this->dtext);
	MergeInfo.setDetailedText(this->dtext);
    MergeInfo.exec();
  }
  else
  {
    for (i=0;i<this->tobj->Gaps.size(); i++)
    {
	  if  (this->tobj->Gaps[i]->cost <=5)
      {
		  this->dtext+= this->tobj->Gaps[i]->stats().c_str();/*"\nTrace " + QString::number(this->tobj->Gaps[i]->Trace1->GetId());
        this->dtext+= " and "+ QString::number(this->tobj->Gaps[i]->Trace2->GetId() );
		this->dtext+="\tcost of:" + QString::number(this->tobj->Gaps[i]->cost); */
        this->HighlightSelected(this->tobj->Gaps[i]->Trace1, .125);
        this->HighlightSelected(this->tobj->Gaps[i]->Trace2, .125);
        this->candidateGaps.push_back( this->tobj->Gaps[i]);
	  } //end of cost<5
    }//end of for merge    
  if (this->candidateGaps.size()>=1)
  {      
      myText+="\nNumber of possible lines:\t" + QString::number(this->candidateGaps.size());
      mergeAll = MergeInfo.addButton("Merge All", QMessageBox::YesRole);
      mergeNone = MergeInfo.addButton("Merge None", QMessageBox::NoRole);
	  MergeInfo.setText(myText); 
	  MergeInfo.setDetailedText(this->dtext);
	  this->poly_line_data->Modified();
	  this->QVTK->GetRenderWindow()->Render();
	  MergeInfo.exec();    
	  if(MergeInfo.clickedButton()==mergeAll)
	  {
		  this->EditLogDisplay->append("Merged " + 
			  QString::number(this->candidateGaps.size()) + " traces");
		  this->EditLogDisplay->append(this->dtext);
		  this->numMerged += (int)this->candidateGaps.size();
		  for (j=0; j<this->candidateGaps.size();j++)
		  {
			tobj->mergeTraces(this->candidateGaps[j]->endPT1,this->candidateGaps[j]->endPT2);
		  }
	  }
	  else if(MergeInfo.clickedButton()==mergeNone)
	  {
		  this->candidateGaps.clear();
		  this->statusBar()->showMessage("Merge Canceled");
	  }
    }   //end if graylist size
    else
    {
		this->statusBar()->showMessage("nothing to merge");
    }   //end else   
  }
  this->ClearSelection();
  this->statusBar()->showMessage(tr("Update Tree Plots"));
  this->TreeModel->SetTraces(this->tobj->GetTraceLines());
  this->statusBar()->showMessage(tr("Done With Merge"));
}
/*  other trace modifiers */
void View3D::SplitTraces()
{
  if(this->SelectedTraceIDs.size()>=1)
    {
    std::unique(this->SelectedTraceIDs.begin(), this->SelectedTraceIDs.end());
    for(unsigned int i = 0; i < this->SelectedTraceIDs.size(); i++)
	{
	  this->tobj->splitTrace(this->SelectedTraceIDs[i]);
	} 
	this->numSplit += (int) this->SelectedTraceIDs.size();
	this->EditLogDisplay->append("split " + QString::number(this->SelectedTraceIDs.size()) +" Traces");
	this->ClearSelection();
	this->statusBar()->showMessage(tr("Update Tree Plots"));
	this->TreeModel->SetTraces(this->tobj->GetTraceLines());
    }
  else
    {
    this->statusBar()->showMessage(tr("Nothing to split"), 1000); 
    }
}

void View3D::FlipTraces()
{
	std::vector<TraceLine*> traceList = this->TreeModel->GetSelectedTraces();
	if (traceList.size() >=1)
	{
		this->EditLogDisplay->append("Fipped " +  QString::number(traceList.size()) + " Traces");
		for (unsigned int i = 0; i < traceList.size(); i++)
		{
			this->EditLogDisplay->append("\tTrace\t" + QString::number(traceList[i]->GetId()));
			this->tobj->ReverseSegment(traceList[i]);
		}		
		this->ClearSelection();
		this->statusBar()->showMessage(tr("Reversed Selected"));
	}
	else
	{
		this->statusBar()->showMessage(tr("Nothing Selected"));
	}
}
void View3D::SetTraceType(int newType)
{
	std::vector<TraceLine*> traceList = this->TreeModel->GetSelectedTraces();
	if (traceList.size() >=1)
	{
		this->EditLogDisplay->append( QString::number(traceList.size())
			+ " Traces set to type: " +  QString::number( newType));
		for (unsigned int i = 0; i < traceList.size(); i++)
		{
			this->EditLogDisplay->append( "\tTrace\t" + QString::number(traceList[i]->GetId()));
			traceList[i]->SetType((unsigned char) newType);
			traceList[i]->setTraceColor(this->tobj->getTraceLUT((unsigned char) newType));
		}
		this->ClearSelection();
		this->TreeModel->SetTraces(this->tobj->GetTraceLines());
		this->statusBar()->showMessage(QString::number(traceList.size())
			+ tr(" traces set to type")+ QString::number(newType));
	}
	else
	{
		this->statusBar()->showMessage(tr("Nothing Selected"));
	}
}
void View3D::SaveToFile()
{
  //display a save file dialog
  QString fileName = QFileDialog::getSaveFileName(
    this,
    tr("Save File"),
    "",
    tr("SWC Images (*.swc);;VTK files (*.vtk)"));

  //if the user pressed cancel, bail out now.
  if(fileName.isNull())
    {
    return;
    }
  this->tempTraceFile.append( fileName);
  //make sure the user supplied an appropriate output file format
  if(!fileName.endsWith(".vtk") && !fileName.endsWith(".swc"))
    {
    QMessageBox::critical(this, tr("Unsupported file format"),
      tr("Trace editor only supports output to .swc and .vtk file formats"));
    this->SaveToFile();
    return;
    }

  if(fileName.endsWith(".swc"))
    {
    this->tobj->WriteToSWCFile(fileName.toStdString().c_str()); 
    }
  else if(fileName.endsWith(".vtk"))
    {
    this->tobj->WriteToVTKFile(fileName.toStdString().c_str()); 
    }
  this->statusBar()->showMessage("File saved as:\t" + fileName.section('/',-2));
  this->EditLogDisplay->append(QString("File saved as: %1  at time: %2").arg(fileName) 
	   .arg(this->Time.currentTime().toString( "h:m:s ap" )));
//Edit Log written to file
  QString logFileName = this->tempTraceFile.last().section('.',0,-1);
  logFileName.append("_log.txt");
  QFile logFile(logFileName);
  if (!logFile.open(QIODevice::WriteOnly | QIODevice::Text))
  {
	  return;
  }
  else
  {
	  this->TraceEditSettings.setValue("lastOpen/Log", logFileName);
	  QTextStream out(&logFile);
	  out << this->EditLogDisplay->toPlainText();
  }
}


/*  Soma display stuff  */
void View3D::AddContourThresholdSliders()
{
  vtkSliderRepresentation2D *sliderRep2 =
    vtkSmartPointer<vtkSliderRepresentation2D>::New();
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

  vtkSliderWidget *sliderWidget2 = vtkSmartPointer<vtkSliderWidget>::New();
  sliderWidget2->SetInteractor(Interactor);
  sliderWidget2->SetRepresentation(sliderRep2);
  sliderWidget2->SetAnimationModeToAnimate();

  vtkSlider2DCallbackContourThreshold *callback_contour =
    vtkSmartPointer<vtkSlider2DCallbackContourThreshold>::New();
  callback_contour->cfilter = this->ContourFilter;
  sliderWidget2->AddObserver(vtkCommand::InteractionEvent,callback_contour);
  sliderWidget2->EnabledOn();

}


void View3D::AddVolumeSliders()
{
  vtkSliderRepresentation2D *opacitySliderRep =
    vtkSliderRepresentation2D::New();
  opacitySliderRep->SetValue(0.1);
  opacitySliderRep->SetTitleText("Opacity");
  opacitySliderRep->GetPoint1Coordinate()->
    SetCoordinateSystemToNormalizedDisplay();
  opacitySliderRep->GetPoint1Coordinate()->SetValue(0.2,0.1);
  opacitySliderRep->GetPoint2Coordinate()->
    SetCoordinateSystemToNormalizedDisplay();
  opacitySliderRep->GetPoint2Coordinate()->SetValue(0.8,0.1);
  opacitySliderRep->SetSliderLength(0.02);
  opacitySliderRep->SetSliderWidth(0.03);
  opacitySliderRep->SetEndCapLength(0.01);
  opacitySliderRep->SetEndCapWidth(0.03);
  opacitySliderRep->SetTubeWidth(0.005);
  opacitySliderRep->SetMinimumValue(0.0);
  opacitySliderRep->SetMaximumValue(1.0);

  this->OpacitySlider = vtkSliderWidget::New();
  this->OpacitySlider->SetInteractor(this->Interactor);
  this->OpacitySlider->SetRepresentation(opacitySliderRep);
  this->OpacitySlider->SetAnimationModeToAnimate();

  vtkSlider2DCallbackOpacity *callback_opacity =
    vtkSlider2DCallbackOpacity::New();
  callback_opacity->volume = this->Volume;
  this->OpacitySlider->AddObserver(vtkCommand::InteractionEvent,callback_opacity);
  this->OpacitySlider->EnabledOn();
  opacitySliderRep->Delete();
  callback_opacity->Delete();
  

// slider 2

  vtkSliderRepresentation2D *brightnessSliderRep =
    vtkSliderRepresentation2D::New();
  brightnessSliderRep->SetValue(0.8);
  brightnessSliderRep->SetTitleText("Brightness");
  brightnessSliderRep->GetPoint1Coordinate()->
    SetCoordinateSystemToNormalizedDisplay();
  brightnessSliderRep->GetPoint1Coordinate()->SetValue(0.2,0.9);
  brightnessSliderRep->GetPoint2Coordinate()->
    SetCoordinateSystemToNormalizedDisplay();
  brightnessSliderRep->GetPoint2Coordinate()->SetValue(0.8,0.9);
  brightnessSliderRep->SetSliderLength(0.02);
  brightnessSliderRep->SetSliderWidth(0.03);
  brightnessSliderRep->SetEndCapLength(0.01);
  brightnessSliderRep->SetEndCapWidth(0.03);
  brightnessSliderRep->SetTubeWidth(0.005);
  brightnessSliderRep->SetMinimumValue(0.0);
  brightnessSliderRep->SetMaximumValue(1.0);

  this->BrightnessSlider = vtkSliderWidget::New();
  this->BrightnessSlider->SetInteractor(this->Interactor);
  this->BrightnessSlider->SetRepresentation(brightnessSliderRep);
  this->BrightnessSlider->SetAnimationModeToAnimate();

  vtkSlider2DCallbackBrightness *callback_brightness =
    vtkSlider2DCallbackBrightness::New();
  callback_brightness->volume = this->Volume;
  this->BrightnessSlider->AddObserver(vtkCommand::InteractionEvent,callback_brightness);
  this->BrightnessSlider->EnabledOn();
  brightnessSliderRep->Delete();
  callback_brightness->Delete();
}
 
void View3D::UndoAction()
{
	
	if(!(this->undoBuff->UndoOrRedo(0)))
	{
		return;
	}
	else
	{
		std::pair<std::string, TraceObject> undostate = this->undoBuff->getState();
		TraceObject newstate = undostate.second;
		*(this->tobj) = newstate;
		Rerender();
	}
	
}


void View3D::RedoAction()
{
	/*
	if(!(this->undoBuff->UndoOrRedo(1)))
	{
		return;
	}
	else
	{
		std::pair<std::string, TraceObject> undostate = this->undoBuff->getState();
		TraceObject newstate = undostate.second;
		*(this->tobj) = newstate;
		Rerender();
	}
	*/
}


void View3D::closeEvent(QCloseEvent *event)
{	
	this->TraceEditSettings.setValue("mainWin/size", this->size());
	this->TraceEditSettings.setValue("mainWin/pos",	this->pos());
	this->TraceEditSettings.setValue("mainWin/gapTol", this->tobj->gapTol ) ;
	this->TraceEditSettings.setValue("mainWin/gapMax", this->tobj->gapMax);
	this->TraceEditSettings.setValue("mainWin/smallLine", this->SmallLineLength);
	this->TraceEditSettings.setValue("mainWin/selectColor", this->SelectColor);
	this->TraceEditSettings.setValue("mainWin/LineWidth", this->lineWidth);
	this->TraceEditSettings.setValue("lastOpen/Image", this->Image);
	this->TraceEditSettings.setValue("lastOpen/Trace",this->TraceFiles);
	this->TraceEditSettings.setValue("lastOpen/Soma", this->SomaFile);
	this->TraceEditSettings.setValue("lastOpen/Temp", this->tempTraceFile);
 if(this->FTKTable)
  {
	  this->TraceEditSettings.setValue("TraceTable/pos", this->FTKTable->pos());
	  this->TraceEditSettings.setValue("TraceTable/size", this->FTKTable->size());
	  this->FTKTable->close();
  }	
 if(this->TreePlot)
  {
	  this->TraceEditSettings.setValue("TracePlot/pos", this->TreePlot->pos());
	  this->TreePlot->close();
  }	
  this->TraceEditSettings.sync();
  if(this->GapsPlotView)
    {
    this->GapsPlotView->close();
    }
  if(this->GapsTableView)
    {
    this->GapsTableView->close();
    }
  event->accept();
}

