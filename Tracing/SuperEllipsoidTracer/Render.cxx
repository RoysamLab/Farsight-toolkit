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

/**
 \brief Code to render the traces . 
 \author $ Author: Amit Mukherjee $
 \version $Revision 1.0 $
 \date May 05 2009
*/
// Copyright Farsight Toolkit 2009, Rensselaer Polytechnic institute Troy NY 12180.

# include <iostream>
#include "itkCommand.h"
#include "itkVector.h"
#include "itkImage.h"
#include "itkVTKImageExport.h"
#include "itkVTKImageImport.h"
#include "itkCurvatureFlowImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkRGBPixel.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "vtkImageData.h"

#include "vtkImageImport.h"
#include "vtkImageExport.h"
#include "vtkImageActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleImage.h" 

#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkCamera.h"
#include "vtkProperty.h"

#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkPointPicker.h"
#include "vtkCellPicker.h"
#include "vtkCommand.h"
#include "vtkRendererCollection.h"

#include "vtkLineSource.h"
#include "vtkTubeFilter.h"
#include "vtkPolyDataMapper.h"
#include "vtkPoints.h"
#include "vtkFloatArray.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkCubeSource.h"
#include "vtkCallbackCommand.h"
#include <stdio.h>
#include "vtkType.h"
#include "vtkPointData.h"
#include "vtkImageData.h"
#include "vtkImageToStructuredPoints.h"
#include "vtkPiecewiseFunction.h"
#include "vtkColorTransferFunction.h"
#include "vtkVolumeProperty.h"
#include "vtkVolumeRayCastMapper.h"
#include "vtkVolumeRayCastCompositeFunction.h"
#include "vtkVolume.h"
#include "vtkOpenGLVolumeTextureMapper3D.h"
#include "vtkSliderWidget.h"
#include "vtkSliderRepresentation2D.h"

#include "TraceNode.h"
#include "tinyxml.h"

#define VTK_CREATE(type, var) vtkSmartPointer<type> var = vtkSmartPointer<type>::New()


typedef itk::Image< unsigned short, 3 > InputImageType;
typedef itk::Image< unsigned char, 3 > ImageType;
typedef itk::VTKImageExport< ImageType > ExportFilterType;

bool ReadNodeXMLFile(const char*, std::vector<TraceNode*>& );
void AddVolumeSliders();

void ConnectPipelines(ExportFilterType::Pointer exporter, vtkImageImport*& importer)
{
  importer->SetUpdateInformationCallback(exporter->GetUpdateInformationCallback());
  importer->SetPipelineModifiedCallback(exporter->GetPipelineModifiedCallback());
  importer->SetWholeExtentCallback(exporter->GetWholeExtentCallback());
  importer->SetSpacingCallback(exporter->GetSpacingCallback());
  importer->SetOriginCallback(exporter->GetOriginCallback());
  importer->SetScalarTypeCallback(exporter->GetScalarTypeCallback());
  importer->SetNumberOfComponentsCallback(exporter->GetNumberOfComponentsCallback());
  importer->SetPropagateUpdateExtentCallback(exporter->GetPropagateUpdateExtentCallback());
  importer->SetUpdateDataCallback(exporter->GetUpdateDataCallback());
  importer->SetDataExtentCallback(exporter->GetDataExtentCallback());
  importer->SetBufferPointerCallback(exporter->GetBufferPointerCallback());
  importer->SetCallbackUserData(exporter->GetCallbackUserData());
}

class vtkSlider2DCallbackBrightness : public vtkCommand
{
public:
  static vtkSlider2DCallbackBrightness *New() 
    { return new vtkSlider2DCallbackBrightness; }
  virtual void Execute(vtkObject *caller, unsigned long f, void* g)
    {
      vtkSliderWidget *sliderWidget = 
        reinterpret_cast<vtkSliderWidget*>(caller);
      double value = static_cast<vtkSliderRepresentation *>(sliderWidget->GetRepresentation())->GetValue();

      vtkSmartPointer<vtkPiecewiseFunction> opacityTransferFunction = vtkSmartPointer<vtkPiecewiseFunction>::New();
    opacityTransferFunction->AddPoint(1,0.0);
    //opacityTransferFunction->AddPoint(2+256*(1-value),0.2);
    opacityTransferFunction->AddPoint(255,value);
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
  virtual void Execute(vtkObject *caller, unsigned long f, void* g)
    {
      vtkSliderWidget *sliderWidget = 
        reinterpret_cast<vtkSliderWidget*>(caller);
      double value = static_cast<vtkSliderRepresentation *>(sliderWidget->GetRepresentation())->GetValue();

    vtkColorTransferFunction *colorTransferFunction = vtkColorTransferFunction::New();
    colorTransferFunction->AddRGBPoint(0.0, 0.0, 0.0, 0.0);
    colorTransferFunction->AddRGBPoint(255*(1-value),1.0,1.0,1.0);
     this->volume->GetProperty()->SetColor(colorTransferFunction);
    }
  vtkSlider2DCallbackContrast() {

  }

  vtkSmartPointer<vtkVolume> volume;
};


double getEuclideanDistance(const std::vector<TraceNode*>& NodeContainer, long nID, TraceNode* nd)	{
	std::vector<TraceNode*>::const_iterator itn = NodeContainer.begin();
	while (itn < NodeContainer.end())	{
		if ((*itn)->ID == nID)	{
			double dist = vcl_sqrt( (nd->loc[0] - (*itn)->loc[0]) * (nd->loc[0] - (*itn)->loc[0]) + 
						(nd->loc[1] - (*itn)->loc[1]) * (nd->loc[1] - (*itn)->loc[1]) + 
						(nd->loc[2] - (*itn)->loc[2]) * (nd->loc[2] - (*itn)->loc[2]) );
			return dist;
		}
	}
	return 0.0;
}
/*
void PrintStatistics(const std::vector<TraceNode*>& NodeContainer)	{
	std::vector<double> length;
	length.resize(NodeContainer.size());
	std::vector<double>::iterator itl = length.begin();
	while (itl < length.end())	{
		(*itl) = 0.0;
		++itl;
	}
	
	std::vector<TraceNode*>::const_iterator itn = NodeContainer.begin();

	while(itn < NodeContainer.end())	{
		for (int i = 0; i < (*itn)->nbrID.size() ; i++) {
			long tID = (*itn)->TraceID;
			if ((*itn)->nbrID[i] < (*itn)->ID)	{
				length[tID] += getEuclideanDistance(NodeContainer, (*itn)->nbrID[i], (*itn));
			}
		}
		++itn;
	}

}
*/
int main(const int argc, char** argv)	{
	if (argc != 3)	{
		std::cout << "Usage: "<<argv[0] << " [ImageFile] [TraceFile]" <<std::endl;
		return EXIT_FAILURE;
	}
	
	std::cout << "Reading " << argv[1] << " ..." ;
    typedef itk::ImageFileReader< InputImageType > ReaderType;
    ReaderType::Pointer reader  = ReaderType::New();
    reader->SetFileName( argv[1] );

	typedef itk::RescaleIntensityImageFilter<InputImageType, ImageType> RescalerType;
	RescalerType::Pointer rescaler = RescalerType::New();
	rescaler->SetOutputMaximum(255);
	rescaler->SetOutputMinimum(0);
	rescaler->SetInput(reader->GetOutput());
    ImageType::Pointer im2 = rescaler->GetOutput();
	try {
		im2->Update();
	}
	catch (itk::ExceptionObject &e)	{
		std::cerr << "Could not read file "<< std::endl << e << std::endl;
	}

	std::cout << "done! " <<std::endl << "Reading " << argv[2] << " ...";

	double spacing[] = {1.0, 1.0, 1.0};
	im2->SetSpacing(spacing);

	ExportFilterType::Pointer itkExporter = ExportFilterType::New();
	itkExporter->SetInput( im2 );

    vtkImageImport* vtkImporter = vtkImageImport::New();  
    ConnectPipelines(itkExporter, vtkImporter);
    
	VTK_CREATE(vtkPiecewiseFunction, opacityTransferFunction);
	opacityTransferFunction->AddPoint(1,0.0);
	//opacityTransferFunction->AddPoint(20,1.0);
	opacityTransferFunction->AddPoint(255,1.0);

	VTK_CREATE(vtkColorTransferFunction, colorTransferFunction);
    colorTransferFunction->AddRGBPoint(0.0, 0.0, 0.0, 0.0);
	colorTransferFunction->AddRGBPoint(20.0,1.0,1.0,1.0);

	VTK_CREATE(vtkVolumeProperty, volumeProperty);
    volumeProperty->SetColor(colorTransferFunction);
    volumeProperty->SetScalarOpacity(opacityTransferFunction);
	volumeProperty->SetInterpolationTypeToLinear();

	VTK_CREATE(vtkOpenGLVolumeTextureMapper3D, volumeMapper);
	volumeMapper->SetSampleDistance(0.75);
	volumeMapper->SetInput(vtkImporter->GetOutput());

	VTK_CREATE(vtkVolume, volume);
    volume->SetMapper(volumeMapper);
    volume->SetProperty(volumeProperty);
	volume->SetPickable(0);
		
	/// display traces
	VTK_CREATE(vtkPolyData, traces);
	VTK_CREATE(vtkFloatArray, scalars);
	scalars->SetNumberOfComponents(1);
	VTK_CREATE(vtkPoints, points);
	points->SetDataTypeToDouble();
	VTK_CREATE(vtkCellArray, cells);

	VTK_CREATE(vtkPolyDataMapper, traceMapper);
	VTK_CREATE(vtkActor, traceActor);

	std::vector<vtkActor*> branchpoints;
	/// read the traces
	std::vector<TraceNode*> NodeContainer;
	if (ReadNodeXMLFile(argv[2], NodeContainer))	{
		
		std::cout << "done! " <<std::endl << "Rendering .... wait!! " << std::endl;
		// create the nodes
		std::vector<TraceNode*>::iterator itr;
		// build a map<ID, PointID>
		std::map<long, vtkIdType> IDmap;
		for(itr = NodeContainer.begin(); itr != NodeContainer.end(); ++itr)	{
			IDmap[(*itr)->ID] = points->InsertNextPoint((*itr)->loc.GetDataPointer());
			scalars->InsertNextTuple1(0.5);
		}
		
		// create the cells
		for(itr = NodeContainer.begin(); itr != NodeContainer.end(); ++itr)	{

			for (int i = 0; i < (*itr)->nbrID.size() ; i++) {
				if ((*itr)->nbrID[i] < (*itr)->ID)	{
					cells->InsertNextCell(2);
					cells->InsertCellPoint(IDmap[(*itr)->ID]);
					cells->InsertCellPoint(IDmap[(*itr)->nbrID[i]]);
				}
			}
			if ((*itr)->nbrID.size() > 2)	{
				vtkSphereSource *sphere = vtkSphereSource::New();
				sphere->SetRadius(3);
				vtkPolyDataMapper* sphereMap = vtkPolyDataMapper::New();
				sphereMap->SetInput(sphere->GetOutput());
				vtkActor* sphereAct = vtkActor::New();
				sphereAct->SetMapper(sphereMap);
				sphereAct->GetProperty()->SetOpacity(1);
				sphereAct->SetPosition( (*itr)->loc[0], (*itr)->loc[1], (*itr)->loc[2] );
				sphereAct->VisibilityOn();
				branchpoints.push_back(sphereAct);
			}
		}
		traces->SetPoints(points);
		traces->SetLines(cells);
		traces->GetPointData()->SetScalars(scalars);
		traceMapper->SetInput(traces);
		traceActor->SetMapper(traceMapper);
		traceActor->GetProperty()->SetPointSize(2.0);
		traceActor->GetProperty()->SetLineWidth(2.0);

		//PrintStatistics(NodeContainer);
	}
	else {
		std::cout << "Failed !! " << std::endl << "Rendering the volume only " << std::endl;
	}
	
	VTK_CREATE(vtkRenderer, renderer);
    renderer->AddVolume(volume);
	if (NodeContainer.size() > 0)	{
		renderer->AddActor(traceActor);
	}
	//add branch point actors
	std::vector<vtkActor*>::iterator bpit = branchpoints.begin();
	for( ; bpit < branchpoints.end(); ++bpit)	{
		renderer->AddActor(*bpit);
	}

    renderer->SetBackground(0.4392, 0.5020, 0.5647);

	VTK_CREATE(vtkRenderWindow, renWin);
	VTK_CREATE(vtkRenderWindowInteractor, iren);

	renWin->SetSize(500, 500);
    renWin->AddRenderer(renderer);
	renWin->SetInteractor(iren);

	VTK_CREATE(vtkInteractorStyleTrackballCamera, style);
    iren->SetRenderWindow(renWin);
    iren->SetInteractorStyle(style);

	
	//Add volume sliders
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
	callback_brightness->volume = volume;
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
	callback_contrast->volume = volume;
	sliderWidget2->AddObserver(vtkCommand::InteractionEvent,callback_contrast);
	sliderWidget2->EnabledOn();

    // Bring up the render window and begin interaction.
    renWin->Render();
    iren->Start();

}

bool ReadNodeXMLFile(const char* xmlfname, std::vector<TraceNode*>& NodeContainer) {
	NodeContainer.reserve(10000);
	TiXmlDocument doc(xmlfname);
	if (!doc.LoadFile()) {
		return false;
	}

	//scan each Superellipse
	TiXmlNode* xmlnode; 

	for ( xmlnode = doc.FirstChild(); xmlnode != 0; xmlnode = xmlnode->NextSibling()) 	{

		//verify if the xmlnode is a type element
		if (xmlnode->Type()!=TiXmlNode::ELEMENT)	{
			continue;
		}

		//verify if the xmlnode is a superellipse, if not 
		if (strcmp(xmlnode->Value(),"Superellipse"))	{
			continue;
		}

		TraceNode *n = new TraceNode();
		TiXmlAttribute* pAttrib = xmlnode->ToElement()->FirstAttribute();
		while (pAttrib)	{
			if (!strcmp(pAttrib->Name(),"ID"))	{
				int temp = -1;
				if (pAttrib->QueryIntValue(&temp)==TIXML_SUCCESS)
					n->ID = temp;
			}
			else if (!strcmp(pAttrib->Name(),"TraceID"))	{
				int temp = -1;
				if (pAttrib->QueryIntValue(&temp)==TIXML_SUCCESS)
					n->TraceID = temp;
			}
			
			else if (!strcmp(pAttrib->Name(),"x"))	{
				double temp = -1.0;
				if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)
					n->loc[0] = temp;
			}

			else if (!strcmp(pAttrib->Name(),"y"))	{
				double temp = -1.0;
				if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)
					n->loc[1] = temp;
			}
			
			else if (!strcmp(pAttrib->Name(),"z"))	{
				double temp = -1.0;
				if (pAttrib->QueryDoubleValue(&temp)==TIXML_SUCCESS)
					n->loc[2] = temp;
			}
			
			pAttrib=pAttrib->Next();
		}

		
		TiXmlNode* nbr; 
		for ( nbr = xmlnode->FirstChild(); nbr != 0; nbr = nbr->NextSibling())		{
			TiXmlAttribute* nAttr = nbr->ToElement()->FirstAttribute();
			if (!strcmp(nAttr->Name(),"ID"))	{
				int temp = -1;
				if (nAttr->QueryIntValue(&temp)==TIXML_SUCCESS)
					n->nbrID.push_back(temp);
			}
		}
		
		//store in container
		NodeContainer.push_back(n);
	}
	return true;

}
