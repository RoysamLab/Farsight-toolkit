/**
 \brief Code to render the traces . 
 \author $ Author: Amit Mukherjee, Arunachalam Narayanaswamy $
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
	im2->Update();
	std::cout << "done! " <<std::endl << "Reading " << argv[2] << " ...";

	ExportFilterType::Pointer itkExporter = ExportFilterType::New();
	itkExporter->SetInput( im2 );

    vtkImageImport* vtkImporter = vtkImageImport::New();  
    ConnectPipelines(itkExporter, vtkImporter);
    
	VTK_CREATE(vtkPiecewiseFunction, opacityTransferFunction);
	opacityTransferFunction->AddPoint(2,0.0);
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
			for (unsigned int i = 0; i < (*itr)->nbrID.size() ; i++) {
				if ((*itr)->nbrID[i] < (*itr)->ID)	{
					cells->InsertNextCell(2);
					cells->InsertCellPoint(IDmap[(*itr)->ID]);
					cells->InsertCellPoint(IDmap[(*itr)->nbrID[i]]);
				}
			}
		}
		traces->SetPoints(points);
		traces->SetLines(cells);
		traces->GetPointData()->SetScalars(scalars);
		traceMapper->SetInput(traces);
		traceActor->SetMapper(traceMapper);
		traceActor->GetProperty()->SetPointSize(2.0);
		traceActor->GetProperty()->SetLineWidth(2.0);
	}
	else {
		std::cout << "Failed !! " << std::endl << "Rendering the volume only " << std::endl;
	}

	
	VTK_CREATE(vtkRenderer, renderer);
    renderer->AddVolume(volume);
	if (NodeContainer.size() > 0)	{
		renderer->AddActor(traceActor);
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

    // Bring up the render window and begin interaction.
    renWin->Render();
    iren->Start();

}

bool ReadNodeXMLFile(const char* xmlfname, std::vector<TraceNode*>& NodeContainer) {
	NodeContainer.reserve(1000);
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
