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
#if defined(_MSC_VER)
#pragma warning(disable : 4996)
#endif
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkVTKImageExport.h>
#include <itkVTKImageImport.h>
#include <vtkImageImport.h>
#include <vtkRenderLargeImage.h>
#include <itkImage.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageSliceConstIteratorWithIndex.h>
#include <itkImageLinearIteratorWithIndex.h>


#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkProperty.h>
#include <vtkOBJReader.h>
#include <vtkAssembly.h>
#include <vtkLODActor.h>
#include <vtkStripper.h>
#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkImageData.h>
#include <vtkImageReader.h>
#include <vtkImageWriter.h>
#include <vector>
#include <vtkStructuredPointsWriter.h>
#include <vtkXMLWriter.h>
#include <vtkTIFFWriter.h>
#include <vtkContourFilter.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkSTLWriter.h>
#include <vtkSTLReader.h>
#include <vtkMarchingCubes.h>
#include <vtkCamera.h>
#include <vtkOutlineFilter.h>
#include <vtkRenderLargeImage.h>
#include <vtkPLYWriter.h>
#include <vtkPLYReader.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkDecimatePro.h>
#include <vtkDiscreteMarchingCubes.h>
#include <vtkQuadricDecimation.h>
#include <vtkQuadricClustering.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkAppendPolyData.h>
#include <vtkPointData.h>

#include "TraceEdit/TraceObject.h"

//#pragma warning(disable:4996)
//#pragma warning(disable:4101)
//#pragma warning(disable:4018)

#define REFLECT 1450
#define SKIP 20

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))
typedef unsigned char InputPixelType;
typedef unsigned char OutputPixelType;

typedef itk::Vector<unsigned char, 3> VectorPixelType;
typedef itk::Image<VectorPixelType, 3> ColorImageType;
typedef itk::Image<VectorPixelType, 2> Color2DImageType;

typedef itk::Image<InputPixelType,3> InputImageType;
typedef itk::Image<OutputPixelType,3> OutputImageType;
typedef itk::Image<short int,3> LabelImageType;
typedef itk::Image<int,3> IntImageType;

typedef itk::Image<InputPixelType,2> Input2DImageType;
typedef itk::Image<InputPixelType,2> Output2DImageType;

typedef itk::ImageRegionConstIterator<InputImageType> ConstIteratorType;
typedef itk::ImageRegionIterator<InputImageType> IteratorType;

typedef itk::ImageRegionConstIterator<LabelImageType> ConstLabelIteratorType;
typedef itk::ImageRegionIterator<LabelImageType> LabelIteratorType;


typedef itk::ImageRegionConstIterator<Input2DImageType> Const2DIteratorType;
typedef itk::ImageRegionIterator<Input2DImageType> twoDIteratorType;

typedef itk::ImageRegionConstIterator<ColorImageType> ConstColorIteratorType;
typedef itk::ImageRegionIterator<ColorImageType> ColorIteratorType;

typedef itk::ImageLinearIteratorWithIndex< Input2DImageType > LinearIteratorType;
typedef itk::ImageSliceConstIteratorWithIndex< InputImageType > SliceIteratorType;

typedef itk::ImageLinearIteratorWithIndex< Color2DImageType > LinearColorIteratorType;
typedef itk::ImageSliceConstIteratorWithIndex< ColorImageType > SliceColorIteratorType;


typedef itk::VTKImageExport<InputImageType> ExportFilterType;


template <typename ITK_Exporter, typename VTK_Importer>
void ConnectPipelines(ITK_Exporter exporter, VTK_Importer* importer)
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

struct cubecoord{
	unsigned short sx,sy,sz;
	unsigned short ex,ey,ez;
};

InputImageType::Pointer getEmpty(int s1,int s2, int s3)
{
	InputImageType::Pointer p = InputImageType::New();
	InputImageType::SizeType size;
	InputImageType::IndexType index;
	InputImageType::RegionType region;
	size[0] = s1; size[1] = s2; size[2] = s3;
	index.Fill(0);
	region.SetSize(size);
	region.SetIndex(index);
	p->SetRegions(region);
	p->Allocate();
	return p;
}
	template <typename T>
typename T::Pointer readImage(const char *filename)
{
	printf("Reading %s ... ",filename);
	typedef typename itk::ImageFileReader<T> ReaderType;
	typename ReaderType::Pointer reader = ReaderType::New();

	ReaderType::GlobalWarningDisplayOff();
	reader->SetFileName(filename);
	try
	{
		reader->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "ExceptionObject caught!" <<std::endl;
		std::cerr << err << std::endl;
		//return EXIT_FAILURE;
	}
	printf("Done.\n");
	return reader->GetOutput();

}



struct Vec3f{
	double x,y,z;
};


void getSmoothed(std::vector<Vec3f> &line)
{
	int pc = 0;
	//Vec3f temp;
	while(pc <2)
	{
		for(unsigned int counter=1; counter<line.size()-1; counter++)
		{	
			line[counter].x=(line[counter-1].x+line[counter+1].x)/2.0;
			line[counter].y=(line[counter-1].y+line[counter+1].y)/2.0;
			line[counter].z=(line[counter-1].z+line[counter+1].z)/2.0;
		}
		pc++;
	}
}

vtkSmartPointer<vtkPolyData> GetLines(char*filename_microgliatrace)
{
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> cellarray = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();
	{
		FILE *fp = fopen(filename_microgliatrace,"r");
		if(fp == NULL)
		{
			printf("Couldnt open %s file for reading\n",filename_microgliatrace);
			return 0;
		}
		printf("Reading file %s\n",filename_microgliatrace);
		float x,y,z;
		bool start = true;
		std::vector<Vec3f> line;
		line.clear();

		while(fscanf(fp,"%f %f %f ",&x,&y,&z)>0)
		{
			if(x==-1)
			{
				getSmoothed(line);

			//	printf("Starting to add trace\n");
				
				for(unsigned int counter=0; counter< line.size(); counter++)
				{
					int ret = points->InsertNextPoint(line[counter].x,line[counter].y,line[counter].z);
					if(counter<line.size()-1)
					{
						cellarray->InsertNextCell(2);
						cellarray->InsertCellPoint(ret);
						cellarray->InsertCellPoint(ret+1);
					}
				}
				line.clear();
				start = true;
				continue;
			}
      //fscanf(fp,"%*f %*f %*d %*d %*d");
			if( fscanf(fp,"%*f %*f") == EOF )
        {
        cerr << "End of file reached within fscanf" << endl;
        }
			start = false;
			Vec3f temp;
			temp.x=x/SKIP;
			temp.y=(y)/SKIP;
			temp.z=z;
			line.push_back(temp);
		}
		poly->SetPoints(points);
		poly->SetLines(cellarray);
		fclose(fp);
	}
	poly->Print(cout);
	printf("Returning Polydata\n");
	return poly;
}

void WriteTraceToPLY(char *filename_trace, char*filename_write)
{
	TraceObject *tobj = new TraceObject();
	tobj->ReadFromRPIXMLFile(filename_trace);
	//tobj->Print(std::cout);
	vtkSmartPointer<vtkPolyData> poly = tobj->GetVTKPolyData();
	poly->GetPointData()->SetScalars(NULL);
	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
	writer->SetFileName(filename_write);
	writer->SetFileTypeToBinary();
	writer->SetInputData(poly);
	writer->Write();
	delete tobj;
	

}

void WriteSWCToPLY(char *filename_trace, char*filename_write)
{
	TraceObject *tobj = new TraceObject();
	tobj->ReadFromSWCFile(filename_trace);
	//tobj->Print(std::cout);
	vtkSmartPointer<vtkPolyData> poly = tobj->GetVTKPolyData();
	poly->GetPointData()->SetScalars(NULL);
	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
	writer->SetFileName(filename_write);
	writer->SetFileTypeToBinary();
	writer->SetInputData(poly);
	writer->Write();
	delete tobj;
	

}

void Binarize(InputImageType::Pointer im)
{
	IteratorType iter(im,im->GetLargestPossibleRegion());
	for(iter.GoToBegin(); !iter.IsAtEnd(); ++iter)
	{
		iter.Set((iter.Get()!=0)?255:0);
	}
}
void WriteLines(char*filename_microgliatrace, char *filename_write)
{

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> cellarray = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();
	{
		FILE *fp = fopen(filename_microgliatrace,"r");
		if(fp == NULL)
		{
			printf("Couldnt open %s file for reading\n",filename_microgliatrace);
			return ;
		}
		printf("Reading file %s\n",filename_microgliatrace);
		float x,y,z;
		bool start = true;
		std::vector<Vec3f> line;
		line.clear();

		while(fscanf(fp,"%f %f %f ",&x,&y,&z)>0)
		{
			if(x==-1)
			{
				//getSmoothed(line);

			//	printf("Starting to add trace\n");
				
				for(unsigned int counter=0; counter< line.size(); counter++)
				{
					int ret = points->InsertNextPoint(line[counter].x,line[counter].y,line[counter].z);
					if(counter<line.size()-1)
					{
						cellarray->InsertNextCell(2);
						cellarray->InsertCellPoint(ret);
						cellarray->InsertCellPoint(ret+1);
					}
				}
				line.clear();
				start = true;
				continue;
			}
      //fscanf(fp,"%*f %*f %*d %*d %*d");
			if( fscanf(fp,"%*f %*f %*f %*f %*f") == EOF )
        {
        cerr << "End of file reached within fscanf" << endl;
        }
			start = false;
			Vec3f temp;
			temp.x=x/SKIP;
			temp.y=(y)/SKIP;
			temp.z=z/SKIP;
			line.push_back(temp);
		}
		poly->SetPoints(points);
		poly->SetLines(cellarray);
		fclose(fp);
	}
	//poly->Print(cout);

	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
	writer->SetFileName(filename_write);
	writer->SetFileTypeToBinary();
	writer->SetInputData(poly);
	writer->Write();


	printf("Returning Polydata\n");
}

#define IM(a,b,c) p[(c)*slicesize+(b)*linesize+(a)]
vtkSmartPointer<vtkImageData> allocateImage(int dimensions[3])
{
	vtkSmartPointer<vtkImageData> imdata = vtkSmartPointer<vtkImageData>::New();
	imdata->SetSpacing(1,1,1);
	imdata->SetOrigin(0,0,0);
	imdata->SetDimensions(dimensions[0],dimensions[1],dimensions[2]);
	imdata->SetNumberOfScalarComponents(1, NULL);
	imdata->AllocateScalars(VTK_UNSIGNED_CHAR, 1);
	printf("Begin memset\n");
	memset(imdata->GetScalarPointer(),0,dimensions[0]*dimensions[1]*dimensions[2]);
	printf("End memset\n");
	return imdata;
}
void loadImageFromFile(char*filename,vtkSmartPointer<vtkImageData> imdata, int dimensions[3])
{
	

	FILE *fp = fopen(filename,"r");
	if(fp==NULL)
	{
		
		printf("file read error\n");
		return;
	}

	double x,y,z;
	//double l1,l2;
	//int slicesize=dimensions[0]*dimensions[1];
	//int linesize=dimensions[1];
	unsigned char *p = (unsigned char*)imdata->GetScalarPointer();
	if(p==NULL)
		printf("got null!\n");
	else
		printf("not null!\n");
	int n=0;
	while(1)
	  {
		if( fscanf(fp,"%lf %lf %lf %*f %*f",&z,&y,&x) == EOF )
      {
      cerr << "End of file reached within fscanf" << endl;
      }
		//fscanf(fp,"%lf %lf %lf %*lf %*lf",&z,&y,&x);
		n++;
		if(n%10000==0)
			printf("%d\r",n);
		if(feof(fp))
			break;
		//for(int counter=0; counter<1; counter++)
	//	{
		*((unsigned char*)(imdata->GetScalarPointer(x/SKIP,y/SKIP,z)))=255;
	//	}
		//IM((int)x/SKIP,(int)y,(int)z)=255;
	  }
	fclose(fp);
	printf("Done\n");
	//imdata->Print(cout);
	//scanf("%*d");
}
void generate_stl(const char *filenames_format_string, int min_n, int max_n,int step_n, int dimensions[], const char * stl_filename)
{
	vtkSmartPointer<vtkImageData> im=allocateImage(dimensions);

	for(int co=min_n; co<=max_n; co+=step_n)
	{
		for(char ch = 'E'; ch <= 'H'; ch++)
		{
			char buff[1024];
			sprintf(buff, filenames_format_string ,co,ch);
			loadImageFromFile(buff,im,dimensions);
		}
	}

	vtkSmartPointer<vtkPolyDataMapper> mapper;
	vtkSmartPointer<vtkLODActor> actor;
	vtkSmartPointer<vtkMarchingCubes> contourf = vtkSmartPointer<vtkMarchingCubes>::New();
	contourf->SetInputData(im);
	contourf->SetValue(0,127);
	contourf->ComputeNormalsOff();
	contourf->ComputeScalarsOff();
	contourf->ComputeGradientsOff();
//	contourf->SetInputMemoryLimit(500000);
	printf("About to update contourf");
	contourf->Update();
	printf("Contour Generated");
	vtkSmartPointer<vtkSmoothPolyDataFilter> smoothf = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
	smoothf->SetInputData(contourf->GetOutput());
	smoothf->SetRelaxationFactor(0.1);
	smoothf->SetNumberOfIterations(20);
	smoothf->Update();
	printf("Contour smoothed\n");
    mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData(smoothf->GetOutput());
	printf("Begin updating mapper\n");
	mapper->Update();
	vtkSmartPointer<vtkPLYWriter> stlwriter = vtkSmartPointer<vtkPLYWriter>::New();
	stlwriter->SetInputData(smoothf->GetOutput());
	stlwriter->SetFileName(stl_filename);
	stlwriter->SetFileTypeToBinary();
	printf("Begin Writing stl file %s\n", stl_filename);
	stlwriter->Write();
	printf("Done writing to stl file");
	//im->Delete();
}


void render_stl_file_new(char filename[][256], int n, double colors[500][3], bool mask[],float scale[])
{
	vtkRenderer* ren1 = vtkRenderer::New();
	double val[6] ;
	for(int counter=0; counter< n; counter++)
	{
		if(mask[counter]==0)
			continue;
		{
		//vtkSmartPointer<vtkPLYReader> reader=vtkSmartPointer<vtkPLYReader>::New();
		vtkSmartPointer<vtkPolyDataReader> reader=vtkSmartPointer<vtkPolyDataReader>::New();
	//	printf("strlen = %d\n",strlen(filename[counter]));
/*
		if(filename[counter][strlen(filename[counter])-1]=='l')//lame way of detecting stl files
			reader = vtkSmartPointer<vtkSTLReader>::New();
		else
			reader = vtkSmartPointer<vtkPLYReader>::New();
			*/
	//	printf("Loading %s\n",filename[counter]);
		reader->SetFileName(filename[counter]);
		vtkSmartPointer<vtkPolyDataNormals> polynorm = vtkSmartPointer<vtkPolyDataNormals>::New();
		polynorm->SetInputData(reader->GetOutput());
		
		/*vtkDecimatePro* decimater = vtkDecimatePro::New();
		decimater->SetInput(reader->GetOutput());
		decimater->PreserveTopologyOn();
		decimater->SetTargetReduction(0.7);*/
		vtkSmartPointer<vtkPolyDataMapper> mapper=vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputData(polynorm->GetOutput());
		mapper->ImmediateModeRenderingOff();

		mapper->Update();

		vtkSmartPointer<vtkActor> actor=vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		//actor->SetColor(0,1,0);
		
		actor->GetProperty()->SetColor(colors[counter][0],colors[counter][1],colors[counter][2]);
		printf("Using color %lf %lf %lf\n",colors[counter][0],colors[counter][1],colors[counter][2]);
	//	actor->RotateZ(90);
		actor->SetScale(scale[0],scale[1],scale[2]);
	
	//	actor->Print(cout);
		actor->GetBounds(val);
		printf("bounds are %lf %lf %lf %lf %lf %lf\n",val[0],val[1],val[2],val[3],val[4],val[5]);
		vtkSmartPointer<vtkOutlineFilter> outlinef = vtkSmartPointer<vtkOutlineFilter>::New();
		outlinef->SetInputData(reader->GetOutput());
		vtkSmartPointer<vtkPolyDataMapper> outlinemap = vtkSmartPointer<vtkPolyDataMapper>::New();
		outlinemap->SetInputData(outlinef->GetOutput());
		vtkSmartPointer<vtkActor> outlineact = vtkSmartPointer<vtkActor>::New();
		outlineact->SetMapper(outlinemap);
		outlineact->GetProperty()->SetColor(1,1,1);
		//outlineact->RotateZ(90);
		actor->GetProperty()->SetLineWidth(1);
		actor->GetProperty()->SetInterpolationToFlat();
		ren1->AddActor(actor);
	//	ren1->AddActor(outlineact);
	//	reader->Delete();
		}
		
	}
	ren1->SetBackground(0,0,0);
	

	vtkRenderWindow* rwin = vtkRenderWindow::New();
	rwin->AddRenderer(ren1);
	rwin->SetSize(1024,1024);
	
	ren1->ResetCamera();
	/*vtkActorCollection* actcol = ren1->GetActors();
	actcol->Print(cout);
	vtkActor * a=actcol->GetLastActor();
	if(a==NULL)
		printf("a is null\n");
	do
	{
		a->AddPosition(0,0,5000);
		a->GetBounds(val);
		printf("bounds are %lf %lf %lf %lf %lf %lf\n",val[0],val[1],val[2],val[3],val[4],val[5]);
	}while(a!=actcol->GetLastActor());*/
	
	//val[3]=0;
	//val[5]=0;
	//ren1->ResetCamera(val);
	/*ren1->GetActiveCamera()->SetFocalPoint(0,0,-5000);
	ren1->GetActiveCamera()->SetPosition(val[0]/2,val[3]/2,5000);
	ren1->GetActiveCamera()->ComputeViewPlaneNormal();
	ren1->GetActiveCamera()->SetViewUp(0,1,0);
	ren1->GetActiveCamera()->SetClippingRange(1,10000);
	*/
	ren1->GetActiveCamera()->Zoom(1.3);
	printf("Distance to focal point %lf\n",ren1->GetActiveCamera()->GetDistance());
	double *pos = ren1->GetActiveCamera()->GetPosition();
	printf("Camera position %lf %lf %lf\n",pos[0],pos[1],pos[2]);
	
	//pos[2]=2355;
	//ren1->GetActiveCamera()->SetPosition(pos);
	//pos = ren1->GetActiveCamera()->GetPosition();
	//printf("Camera position %lf %lf %lf\n",pos[0],pos[1],pos[2]);
	//printf("Distance to focal point %lf\n",ren1->GetActiveCamera()->GetDistance());

	vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::New();
	iren->SetRenderWindow(rwin);
	vtkInteractorStyleTrackballCamera* style = vtkInteractorStyleTrackballCamera::New();
  iren->SetInteractorStyle( style );

	rwin->Render();
	vtkRenderLargeImage *renderLarge=vtkRenderLargeImage::New();
	renderLarge->SetInput(ren1);
	renderLarge->SetMagnification(1);
	vtkTIFFWriter *writer = vtkTIFFWriter::New();
	writer->SetInputData(renderLarge->GetOutput());
	writer->SetFileName("astrocyte_traces.tif");
	//writer->Write();
	iren->Start();
}
void render_stl_file(char filename[][256], int n, double colors[500][3],bool mask[])
{


	vtkRenderer* ren1 = vtkRenderer::New();
	double val[6] ;
	for(int counter=0; counter< n; counter++)
	{
		if(mask[counter]==0)
			continue;
		if(counter>4)
		{
		vtkPolyDataReader* reader=vtkPolyDataReader::New();
		printf("strlen = %d\n",(int)strlen(filename[counter]));
/*
		if(filename[counter][strlen(filename[counter])-1]=='l')//lame way of detecting stl files
			reader = vtkSmartPointer<vtkSTLReader>::New();
		else
			reader = vtkSmartPointer<vtkPLYReader>::New();
			*/
		printf("Loading %s\n",filename[counter]);
		reader->SetFileName(filename[counter]);
		/*vtkDecimatePro* decimater = vtkDecimatePro::New();
		decimater->SetInput(reader->GetOutput());
		decimater->PreserveTopologyOn();
		decimater->SetTargetReduction(0.7);*/
		vtkSmartPointer<vtkPolyDataMapper> mapper=vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputData(reader->GetOutput());
		mapper->ImmediateModeRenderingOn();
		mapper->Update();
		vtkSmartPointer<vtkLODActor> actor=vtkSmartPointer<vtkLODActor>::New();
		actor->SetMapper(mapper);
		//actor->SetColor(0,1,0);
		
		actor->GetProperty()->SetColor(colors[counter][0],colors[counter][1],colors[counter][2]);
	//	actor->RotateZ(90);
		actor->SetScale(1,1,3);
	//	actor->Print(cout);
		actor->GetBounds(val);
		printf("bounds are %lf %lf %lf %lf %lf %lf\n",val[0],val[1],val[2],val[3],val[4],val[5]);
		vtkSmartPointer<vtkOutlineFilter> outlinef = vtkSmartPointer<vtkOutlineFilter>::New();
		outlinef->SetInputData(reader->GetOutput());
		vtkSmartPointer<vtkPolyDataMapper> outlinemap = vtkSmartPointer<vtkPolyDataMapper>::New();
		outlinemap->SetInputData(outlinef->GetOutput());
		vtkSmartPointer<vtkActor> outlineact = vtkSmartPointer<vtkActor>::New();
		outlineact->SetMapper(outlinemap);
		outlineact->GetProperty()->SetColor(1,1,1);
		//outlineact->RotateZ(90);
		actor->GetProperty()->SetLineWidth(1);
		ren1->AddActor(actor);
	//	ren1->AddActor(outlineact);
		reader->Delete();
		}
		else
		{
		vtkPLYReader* reader=vtkPLYReader::New();
		printf("strlen = %d\n",(int)strlen(filename[counter]));
/*
		if(filename[counter][strlen(filename[counter])-1]=='l')//lame way of detecting stl files
			reader = vtkSmartPointer<vtkSTLReader>::New();
		else
			reader = vtkSmartPointer<vtkPLYReader>::New();
			*/
		printf("Loading %s\n",filename[counter]);
		reader->SetFileName(filename[counter]);
		/*vtkDecimatePro* decimater = vtkDecimatePro::New();
		decimater->SetInput(reader->GetOutput());
		decimater->PreserveTopologyOn();
		decimater->SetTargetReduction(0.7);*/
		vtkSmartPointer<vtkPolyDataMapper> mapper=vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputData(reader->GetOutput());
		mapper->ImmediateModeRenderingOn();
		mapper->Update();
		vtkSmartPointer<vtkLODActor> actor=vtkSmartPointer<vtkLODActor>::New();
		actor->SetMapper(mapper);
		//actor->SetColor(0,1,0);
		
		actor->GetProperty()->SetColor(colors[counter][0],colors[counter][1],colors[counter][2]);
		//actor->RotateZ(90);
		actor->SetScale(1,1,3);
	//	actor->Print(cout);
		actor->GetBounds(val);
		printf("bounds are %lf %lf %lf %lf %lf %lf\n",val[0],val[1],val[2],val[3],val[4],val[5]);
		vtkSmartPointer<vtkOutlineFilter> outlinef = vtkSmartPointer<vtkOutlineFilter>::New();
		outlinef->SetInputData(reader->GetOutput());
		vtkSmartPointer<vtkPolyDataMapper> outlinemap = vtkSmartPointer<vtkPolyDataMapper>::New();
		outlinemap->SetInputData(outlinef->GetOutput());
		vtkSmartPointer<vtkActor> outlineact = vtkSmartPointer<vtkActor>::New();
		outlineact->SetMapper(outlinemap);
		outlineact->GetProperty()->SetColor(1,1,1);
		//outlineact->RotateZ(90);
		
		ren1->AddActor(actor);
	//	ren1->AddActor(outlineact);
		reader->Delete();
		}
		
	}
	ren1->SetBackground(0,0,0);
	

	vtkRenderWindow* rwin = vtkRenderWindow::New();
	rwin->AddRenderer(ren1);
	rwin->SetSize(1024,1024);
	
	ren1->ResetCamera();
	/*vtkActorCollection* actcol = ren1->GetActors();
	actcol->Print(cout);
	vtkActor * a=actcol->GetLastActor();
	if(a==NULL)
		printf("a is null\n");
	do
	{
		a->AddPosition(0,0,5000);
		a->GetBounds(val);
		printf("bounds are %lf %lf %lf %lf %lf %lf\n",val[0],val[1],val[2],val[3],val[4],val[5]);
	}while(a!=actcol->GetLastActor());*/
	
	//val[3]=0;
	//val[5]=0;
	//ren1->ResetCamera(val);
	/*ren1->GetActiveCamera()->SetFocalPoint(0,0,-5000);
	ren1->GetActiveCamera()->SetPosition(val[0]/2,val[3]/2,5000);
	ren1->GetActiveCamera()->ComputeViewPlaneNormal();
	ren1->GetActiveCamera()->SetViewUp(0,1,0);
	ren1->GetActiveCamera()->SetClippingRange(1,10000);
	*/
	ren1->GetActiveCamera()->Zoom(1.3);
	printf("Distance to focal point %lf\n",ren1->GetActiveCamera()->GetDistance());
	double *pos = ren1->GetActiveCamera()->GetPosition();
	printf("Camera position %lf %lf %lf\n",pos[0],pos[1],pos[2]);
	
	//pos[2]=2355;
	//ren1->GetActiveCamera()->SetPosition(pos);
	//pos = ren1->GetActiveCamera()->GetPosition();
	//printf("Camera position %lf %lf %lf\n",pos[0],pos[1],pos[2]);
	//printf("Distance to focal point %lf\n",ren1->GetActiveCamera()->GetDistance());

	vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::New();
	iren->SetRenderWindow(rwin);

	rwin->Render();
	vtkRenderLargeImage *renderLarge=vtkRenderLargeImage::New();
	renderLarge->SetInput(ren1);
	renderLarge->SetMagnification(1);
	vtkTIFFWriter *writer = vtkTIFFWriter::New();
	writer->SetInputData(renderLarge->GetOutput());
	writer->SetFileName("astrocyte_traces.tif");
	writer->Write();
	iren->Start();
	
}

void render_single_stl_file(char filename[],double colors[][3])
{

	vtkRenderer* ren1 = vtkRenderer::New();
	double val[6] ;
	int counter = 0;
	{
		{
		vtkPolyDataReader* reader=vtkPolyDataReader::New();
	//	printf("strlen = %d\n",strlen(filename[counter]));
/*
		if(filename[counter][strlen(filename[counter])-1]=='l')//lame way of detecting stl files
			reader = vtkSmartPointer<vtkSTLReader>::New();
		else
			reader = vtkSmartPointer<vtkPLYReader>::New();
			*/
		printf("Loading %s\n",filename);
		reader->SetFileName(filename);
		/*vtkDecimatePro* decimater = vtkDecimatePro::New();
		decimater->SetInput(reader->GetOutput());
		decimater->PreserveTopologyOn();
		decimater->SetTargetReduction(0.7);*/
		vtkSmartPointer<vtkPolyDataMapper> mapper=vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputData(reader->GetOutput());
		mapper->ImmediateModeRenderingOff();
		mapper->Update();
		vtkSmartPointer<vtkLODActor> actor=vtkSmartPointer<vtkLODActor>::New();
		actor->SetMapper(mapper);
		//actor->SetColor(0,1,0);
		
		actor->GetProperty()->SetColor(colors[counter][0],colors[counter][1],colors[counter][2]);
		actor->RotateZ(90);
		actor->SetScale(1,1,3);
	//	actor->Print(cout);
		actor->GetBounds(val);
		printf("bounds are %lf %lf %lf %lf %lf %lf\n",val[0],val[1],val[2],val[3],val[4],val[5]);
		vtkSmartPointer<vtkOutlineFilter> outlinef = vtkSmartPointer<vtkOutlineFilter>::New();
		outlinef->SetInputData(reader->GetOutput());
		vtkSmartPointer<vtkPolyDataMapper> outlinemap = vtkSmartPointer<vtkPolyDataMapper>::New();
		outlinemap->SetInputData(outlinef->GetOutput());
		vtkSmartPointer<vtkActor> outlineact = vtkSmartPointer<vtkActor>::New();
		outlineact->SetMapper(outlinemap);
		outlineact->GetProperty()->SetColor(1,1,1);
		//outlineact->RotateZ(90);
		
		ren1->AddActor(actor);
	//	ren1->AddActor(outlineact);
		reader->Delete();
		}
		
	}
	ren1->SetBackground(0.25,0.25,0.25);
	

	vtkRenderWindow* rwin = vtkRenderWindow::New();
	rwin->AddRenderer(ren1);
	rwin->SetSize(500,500);
	
	ren1->ResetCamera();
	/*vtkActorCollection* actcol = ren1->GetActors();
	actcol->Print(cout);
	vtkActor * a=actcol->GetLastActor();
	if(a==NULL)
		printf("a is null\n");
	do
	{
		a->AddPosition(0,0,5000);
		a->GetBounds(val);
		printf("bounds are %lf %lf %lf %lf %lf %lf\n",val[0],val[1],val[2],val[3],val[4],val[5]);
	}while(a!=actcol->GetLastActor());*/
	
	//val[3]=0;
	//val[5]=0;
	//ren1->ResetCamera(val);
	/*ren1->GetActiveCamera()->SetFocalPoint(0,0,-5000);
	ren1->GetActiveCamera()->SetPosition(val[0]/2,val[3]/2,5000);
	ren1->GetActiveCamera()->ComputeViewPlaneNormal();
	ren1->GetActiveCamera()->SetViewUp(0,1,0);
	ren1->GetActiveCamera()->SetClippingRange(1,10000);
	*/
//	ren1->GetActiveCamera()->Zoom(3.3);
	printf("Distance to focal point %lf\n",ren1->GetActiveCamera()->GetDistance());
	double *pos = ren1->GetActiveCamera()->GetPosition();
	printf("Camera position %lf %lf %lf\n",pos[0],pos[1],pos[2]);
	
	//pos[2]=2355;
	//ren1->GetActiveCamera()->SetPosition(pos);
	//pos = ren1->GetActiveCamera()->GetPosition();
	//printf("Camera position %lf %lf %lf\n",pos[0],pos[1],pos[2]);
	//printf("Distance to focal point %lf\n",ren1->GetActiveCamera()->GetDistance());

	vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::New();
	iren->SetRenderWindow(rwin);

	rwin->Render();
	vtkRenderLargeImage *renderLarge=vtkRenderLargeImage::New();
	renderLarge->SetInput(ren1);
	renderLarge->SetMagnification(12);
	//vtkTIFFWriter *writer = vtkTIFFWriter::New();
//	writer->SetInput(renderLarge->GetOutput());
//	writer->SetFileName("SiteVisit123_30000X9000.tif");
	//writer->Write();
	iren->Start();
	
}

vtkSmartPointer<vtkPolyData> getVTKPolyDataPrecise(LabelImageType::Pointer label)
{
	LabelIteratorType liter = LabelIteratorType(label,label->GetLargestPossibleRegion());
	liter.GoToBegin();

	//find the maximum number of cells
	unsigned short max1 = 0;
	for(liter.GoToBegin();!liter.IsAtEnd();++liter)
		max1 = MAX(max1,liter.Get());

	//find all the cubes in which cells lie
	cubecoord* carray = new cubecoord[max1+1];
	for(int counter=0; counter<=max1; counter++)
	{
		carray[counter].sx=60000; carray[counter].sy=60000;carray[counter].sz=60000;
		carray[counter].ex=0;carray[counter].ey=0;carray[counter].ez=0;
	}

	typedef itk::ImageRegionConstIteratorWithIndex<LabelImageType> ConstLabelIteratorWithIndex;
	ConstLabelIteratorWithIndex cliter = ConstLabelIteratorWithIndex(label,label->GetLargestPossibleRegion());
	InputImageType::IndexType index;
	for(cliter.GoToBegin();!cliter.IsAtEnd();++cliter)
	{
		int cur = cliter.Get();
		if(cur!=0)
		{
			index = cliter.GetIndex();
			carray[cur].sx= MIN(index[0],carray[cur].sx);
			carray[cur].sy= MIN(index[1],carray[cur].sy);
			carray[cur].sz= MIN(index[2],carray[cur].sz);
			carray[cur].ex= MAX(index[0],carray[cur].ex);
			carray[cur].ey= MAX(index[1],carray[cur].ey);
			carray[cur].ez= MAX(index[2],carray[cur].ez);
		}
	}

	//find the largest image size we need
	unsigned short wx=0,wy=0,wz=0;
	for(int counter=1; counter<=max1; counter++)
	{
		wx = MAX(carray[counter].ex-carray[counter].sx+1,wx);
		wy = MAX(carray[counter].ey-carray[counter].sy+1,wy);
		wz = MAX(carray[counter].ez-carray[counter].sz+1,wz);
	}
	// accommodate padding
	wx = wx+2;wy = wy +2; wz = wz+2;
	printf("wx wy wz %u %u %u\n",wx,wy,wz);
	// create a tiny image of maximum size


	//appendfilter->UserManagedInputsOn();
	//appendfilter->SetNumberOfInputs(max1);
	vtkSmartPointer<vtkAppendPolyData> appendfilter = vtkSmartPointer<vtkAppendPolyData>::New();

/**************/
			ExportFilterType::Pointer itkexporter = ExportFilterType::New();
			
			vtkSmartPointer<vtkImageImport> vtkimporter = vtkSmartPointer<vtkImageImport>::New();
			ConnectPipelines(itkexporter,(vtkImageImport *)vtkimporter);
			vtkSmartPointer<vtkMarchingCubes> contourf = vtkSmartPointer<vtkMarchingCubes>::New();
			contourf->SetInputData(vtkimporter->GetOutput());
			contourf->SetValue(0,127);
			contourf->ComputeNormalsOff();
			contourf->ComputeScalarsOff();
			contourf->ComputeGradientsOff();
			vtkSmartPointer<vtkSmoothPolyDataFilter> smoothf = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
			smoothf->SetInputData(contourf->GetOutput());
			smoothf->SetRelaxationFactor(0.3);
			smoothf->SetNumberOfIterations(20);

			vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
			transform->PostMultiply();
	
			transform->Identity();
			vtkSmartPointer<vtkTransformPolyDataFilter> tf = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
			tf->SetTransform(transform);
			tf->SetInputData(smoothf->GetOutput());
		
/******************/



InputImageType::Pointer t = getEmpty(wx,wy,wz);
	for(int counter=1; counter<=max1; counter++)
	{

		
			
			//printf("Maximum tiny image size I need is [%d %d %d]\n",wx,wy,wz);

			if(carray[counter].sx > 59999)
				continue;
			printf("Working..\n");
	//		scanf("%*d");
			InputImageType::SizeType size;
			InputImageType::RegionType region;
			index.Fill(1);

			region.SetIndex(index);
			region.SetSize(size);


			LabelImageType::SizeType lsize;
			LabelImageType::IndexType lindex;
			LabelImageType::RegionType lregion;

			itkexporter->SetInput(t);

			
			
		

			t->FillBuffer(0);
			lsize[0] = carray[counter].ex-carray[counter].sx+1;
			lsize[1] = carray[counter].ey-carray[counter].sy+1;
			lsize[2] = carray[counter].ez-carray[counter].sz+1;

			lindex[0] = carray[counter].sx;
			lindex[1] = carray[counter].sy;
			lindex[2] = carray[counter].sz;

			lregion.SetIndex(lindex);
			lregion.SetSize(lsize);
			LabelIteratorType localiter = LabelIteratorType(label,lregion);

			size = lsize;
			region.SetSize(size);
			IteratorType iter = IteratorType(t,region);
			for(localiter.GoToBegin(),iter.GoToBegin();!localiter.IsAtEnd();++localiter,++iter)
			{
				if(localiter.Get()==counter)
				{
					iter.Set(255);
				}
			}
			t->Modified();
			vtkimporter->Modified();

			transform->Identity();	
			transform->Translate(carray[counter].sx-1,carray[counter].sy-1,carray[counter].sz-1);
			tf->SetTransform(transform);
			tf->Update();
		vtkSmartPointer<vtkPolyData> pol=vtkSmartPointer<vtkPolyData>::New();
		pol->DeepCopy(tf->GetOutput());
	//	tf->GetOutput()->Print(std::cout);
		
		appendfilter->AddInputData(pol);
		//appendfilter->Update();
	
		//appendfilter->SetInputByNumber(counter-1,tf->GetOutput());
	//	appendfilter->Update();
	//	appendfilter->GetOutput()->Print(std::cout);
		//if(counter>500)
		//	break;
		printf("Completed %d/%d\r",counter,max1);
	//	scanf("%*d");
	}

	appendfilter->Update();
	vtkSmartPointer<vtkDecimatePro> decimate = vtkSmartPointer<vtkDecimatePro>::New();
	decimate->SetInputData(appendfilter->GetOutput());
	decimate->SetTargetReduction(0.1);
	//decimate->SetNumberOfDivisions(32,32,32);
	printf("Decimating the contours...");
	decimate->Update();
	printf("Done\n");
	printf("Smoothing the contours after decimation...");
	vtkSmartPointer<vtkSmoothPolyDataFilter> smoothfinal = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
	smoothfinal->SetRelaxationFactor(0.2);
	smoothfinal->SetInputData(decimate->GetOutput());
	smoothfinal->SetNumberOfIterations(0);
	smoothfinal->Update();
	printf("Done\n");
	delete [] carray;
	vtkSmartPointer<vtkPolyData> out = smoothfinal->GetOutput();
	return out;
}
void generate_stl_from_tif(char *filename_tif,char *stl_filename)
{
	vtkSmartPointer<vtkImageData> im;	
	ExportFilterType::Pointer itkexporter = ExportFilterType::New();
	LabelImageType::Pointer itkim = readImage<LabelImageType>(filename_tif);
	//Binarize(itkim);
	vtkSmartPointer<vtkPolyData> poly = getVTKPolyDataPrecise(itkim);
	//itkexporter->SetInput(itkim);
	//vtkSmartPointer<vtkImageImport> vtkimporter = vtkSmartPointer<vtkImageImport>::New();
	//ConnectPipelines(itkexporter,(vtkImageImport *)vtkimporter);
	//vtkimporter->Update();
	//im = vtkimporter->GetOutput();

	vtkSmartPointer<vtkPolyDataMapper> mapper;
//	vtkSmartPointer<vtkLODActor> actor;
//	vtkSmartPointer<vtkMarchingCubes> contourf = vtkSmartPointer<vtkMarchingCubes>::New();
//	contourf->SetInput(im);
//	contourf->SetValue(0,127);
//	contourf->ComputeNormalsOff();
//	contourf->ComputeScalarsOff();
//	contourf->ComputeGradientsOff();
////	contourf->SetInputMemoryLimit(500000);
//	printf("About to update contourf\n");
//	contourf->Update();
//	printf("Contour generated\n");
//	vtkSmartPointer<vtkSmoothPolyDataFilter> smoothf = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
//	smoothf->SetInput(contourf->GetOutput());
//	smoothf->SetRelaxationFactor(0.2);
//	smoothf->SetNumberOfIterations(10);
//	smoothf->Update();
//	printf("Contour smoothed\n");
//	vtkSmartPointer<vtkDecimatePro> decimate = vtkSmartPointer<vtkDecimatePro>::New();
//	decimate->SetInput(smoothf->GetOutput());
//	decimate->SetTargetReduction(0.9);
//	decimate->Update();
//	printf("Contour decimated\n");
    mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData(poly);
	printf("Begin updating mapper\n");
	mapper->Update();
	vtkSmartPointer<vtkPolyDataWriter> stlwriter = vtkSmartPointer<vtkPolyDataWriter>::New();
	stlwriter->SetInputData(poly);
	stlwriter->SetFileName(stl_filename);
	stlwriter->SetFileTypeToBinary();
	printf("Begin Writing stl file %s\n", stl_filename);
	stlwriter->Write();
	printf("Done writing to stl file");
	return;	
}

bool file_exists(char *filename)
{
	FILE * fp = fopen(filename,"r");
	if(fp==NULL)
	{
		return false;
	}
	else
	{
		fclose(fp);
		return true;
	}
}
int main(int argc, char**argv)
{
	if(argc<2)
	{
		printf("Usage: exe input_parameter_file\n");
		return 1;
	}
	FILE * fp = fopen(argv[1],"r");
	if(fp==NULL)
	{
		printf("Couldn't load the parameters file\n");
		return 1;
	}
	char buff[1024];
	char input_filenames[500][256];
	char stl_filenames[500][256];
	double colors[500][3];
	bool mask[500];
	int counter=0;
	float scale[3];
	char bufff[1024];
	fgets(bufff,1023,fp);
	if( sscanf(bufff,"%f %f %f",scale, scale+1,scale+2) == EOF )
    {
    cerr << "End of file reached within fscanf" << endl;
    }
	printf("scale = %f %f %f\n",scale[0],scale[1],scale[2]);
	while(fgets(buff,1000,fp)!=NULL)
	{
		char filename[1024];

		sscanf(buff,"%s %lf %lf %lf",filename,colors[counter],colors[counter]+1,colors[counter]+2);
		printf(" I read filename = |%s| colors[] = |%lf| |%lf| |%lf|\n",filename, colors[counter][0],colors[counter][1], colors[counter][2]);
		if(strcmp(&filename[strlen(filename)-4],".xml")==0)
		{
			//trace
			strcpy(input_filenames[counter],filename);
			sprintf(stl_filenames[counter],"cache/cache_%s.vtk",filename);
			if(!file_exists(stl_filenames[counter]))
				WriteTraceToPLY(input_filenames[counter],stl_filenames[counter]);
		}
		else if(strcmp(&filename[strlen(filename)-4],".swc")==0)
		{
			strcpy(input_filenames[counter],filename);
			sprintf(stl_filenames[counter],"cache/cache_%s.vtk",filename);
			if(!file_exists(stl_filenames[counter]))
				WriteSWCToPLY(input_filenames[counter],stl_filenames[counter]);
		}
		else
		{
			//image
			strcpy(input_filenames[counter],filename);
			sprintf(stl_filenames[counter],"cache/cache_%s.vtk",filename);
			if(!file_exists(stl_filenames[counter]))
				generate_stl_from_tif(input_filenames[counter],stl_filenames[counter]);
		}
		mask[counter] = 1;
		counter++;
		//fscanf(fp,"%*c");
	}
	fclose(fp);
	render_stl_file_new(&stl_filenames[0],counter,colors,mask,scale);
	//render_single_stl_file(stl_filenames[0],&colors[0]);
	return 0;

}
int main_old(int argc, char**argv)
{
	vtkSphereSource* sphere = vtkSphereSource::New();
	sphere->SetRadius(1.0);
	sphere->SetThetaResolution(18);
	sphere->SetPhiResolution(18);

	/*const char *filenames[]={"w44.obj","w43.obj","w42.obj","w41.obj",
		"w34.obj","w33.obj","w32.obj","w31.obj",
		"w24.obj","w23.obj","w22.obj","w21.obj",
		"w14.obj","w13.obj","w12.obj","w11.obj",
		"vessel1.obj","vessel2.obj","vessel3.obj","vessel4.obj"};*/

	/*const char *filename_tracefiles[]={"transformed_100um2percent25x1unmixed02TracedPoints.txt",
		"transformed_100um2percent25x1unmixed03TracedPoints.txt",
		"transformed_100um2percent25x2unmixed02TracedPoints.txt",
		"transformed_100um2percent25x2unmixed03TracedPoints.txt",
		"transformed_100um2percent25x3unmixed02TracedPoints.txt",
		"transformed_100um2percent25x3unmixed03TracedPoints.txt",
		"transformed_100um2percent25x4unmixed02TracedPoints.txt",
		"transformed_100um2percent25x4unmixed03TracedPoints.txt"};*/

	/*double colors[][3]={ {1,0,0},{1,0,1},{1,1,0},{0,1,0},
		{1,0,0},{1,0,1},{1,1,0},{0,1,0},
		{1,0,0},{1,0,1},{1,1,0},{0,1,0},
		{1,0,0},{1,0,1},{1,1,0},{0,1,0},
		{0,1,1},{0,1,1},{0,1,1},{0,1,1} };*/

	/*double trace_colors[][3]={ {1,0,0},{1,1,0} };*/
	vtkSmartPointer<vtkOBJReader> reader;
	vtkSmartPointer<vtkPolyDataMapper> mapper;
	vtkSmartPointer<vtkLODActor> actor;
	//vtkStripper *stripper;
	//vtkRenderer* ren1 = vtkRenderer::New();
	//vtkAssembly * assembly = vtkAssembly::New();
	vtkSmartPointer<vtkPolyData> polydata;
	
	int dimensions[3]={400,1400,70};
	
	//const char *stl_filenames[]={"transformed_to_Montage8Hunmixed_fused.tiff/class1.ply","transformed_to_Montage8Hunmixed_fused.tiff/class2.ply","transformed_to_Montage8Hunmixed_fused.tiff/class3.ply","transformed_to_Montage8Hunmixed_fused.tiff/class4.ply","transformed_to_Montage8Hunmixed_fused.tiff/vessel_montage.ply","test_write_lines.ply","vessel.stl","nuclei_class1.stl","nuclei_class2.stl","nuclei_class3.stl","nuclei_class4.stl"};
	//generate_stl("transformed_100um2percent25x%dunmixed04.npts",1,4,1,dimensions,"vessel.stl");
 //   generate_stl("Nuclei_%d.xml_class1.txt",1,4,1,dimensions,"nuclei_class1.stl");
//	generate_stl("Nuclei_%d.xml_class2.txt",1,4,1,dimensions,"nuclei_class2.stl");
//	generate_stl("Nuclei_%d.xml_class3.txt",1,4,1,dimensions,"nuclei_class3.stl");
//	generate_stl("Nuclei_%d.xml_class4.txt",1,4,1,dimensions,"nuclei_class4.stl");
	//render_stl_file("auto_generate_vessel.stl");
	/*double col[][3]={
    {0,1,0},
		{1,1,0},
		{1,0,1},
		{1,0,0},
		{0,1,1}
	};*/

	//const char *filename_traces_format_input  = "transformed_to_Montage8Hunmixed_fused.tiff/transformed_Montage%d%cunmixed0%dTracedPoints.txt";
	const char *filename_traces_format_output = "transformed_to_Montage8Hunmixed_fused.tiff/transformed_Montage%d%cunmixed0%dTracedPoints.ply";
	const char *filename_vessel_input = "transformed_to_Montage8Hunmixed_fused.tiff/transformed_Montage%d%cunmixed01.npts";
	/*const char *filename_cells_input[4]={"transformed_to_Montage8Hunmixed_fused.tiff/Montage%d%cunmixed02.xml_class1.txt",
								   "transformed_to_Montage8Hunmixed_fused.tiff/Montage%d%cunmixed02.xml_class2.txt",
								   "transformed_to_Montage8Hunmixed_fused.tiff/Montage%d%cunmixed02.xml_class3.txt",
								   "transformed_to_Montage8Hunmixed_fused.tiff/Montage%d%cunmixed02.xml_class4.txt"};*/

	generate_stl(filename_vessel_input,1,15,1,dimensions,"transformed_to_Montage8Hunmixed_fused.tiff/vessel_montage.ply");
	
	//generate_stl(filename_cells_input[0],1,15,1,dimensions,stl_filenames[0]);

	//generate_stl(filename_cells_input[1],1,15,1,dimensions,stl_filenames[1]);

	//generate_stl(filename_cells_input[2],1,15,1,dimensions,stl_filenames[2]);

	//generate_stl(filename_cells_input[3],1,15,1,dimensions,stl_filenames[3]);

	
	//return 0;
	//char buff_in[1024];
	//char buff_out[1024];
	



	//for(int counter=1; counter<16; counter++)
	//{
	//	for(char ch = 'E'; ch <='H'; ch++)
	//	{
	//		for( int co = 3; co <=5; co+=2)
	//		{
	//			sprintf(buff_in,filename_traces_format_input,counter,ch,co);
	//			sprintf(buff_out,filename_traces_format_output,counter,ch,co);
	//			printf("Processing %s...",buff_in);
	//			WriteLines(buff_in,buff_out);
	//			printf("Done\n");
	//		}
	//	}
	//}

	char *filenames_traces_ply[125];
	double color_t[125][3];


	int pc = 0;
	for(int counter=1; counter<16; counter++)
	{
		for(char ch = 'E'; ch <='H'; ch++)
		{
			for( int co = 3; co <=5; co+=2)
			{
				if(co==5)
				{
					color_t[pc][0]=1;color_t[pc][1]=0;color_t[pc][2]=0;
				}
				else
				{
					color_t[pc][0]=1;color_t[pc][1]=1;color_t[pc][2]=0;
				}
				filenames_traces_ply[pc]= new char[1024];
				sprintf(filenames_traces_ply[pc++],filename_traces_format_output,counter,ch,co);
			}
		}
	}
	for(int c= 0; c<5; c++)
		filenames_traces_ply[pc++]=new char[1024];


	strcpy(filenames_traces_ply[120],"transformed_to_Montage8Hunmixed_fused.tiff/class1.ply");
	strcpy(filenames_traces_ply[121],"transformed_to_Montage8Hunmixed_fused.tiff/class2.ply");
	strcpy(filenames_traces_ply[122],"transformed_to_Montage8Hunmixed_fused.tiff/class3.ply");
	strcpy(filenames_traces_ply[123],"transformed_to_Montage8Hunmixed_fused.tiff/class4.ply");
	strcpy(filenames_traces_ply[124],"transformed_to_Montage8Hunmixed_fused.tiff/vessel_montage.ply");

	color_t[120][0]=0;color_t[120][1]=1;color_t[120][2]=0;
	color_t[121][0]=1;color_t[121][1]=1;color_t[121][2]=0;
	color_t[122][0]=1;color_t[122][1]=0;color_t[122][2]=1;
	color_t[123][0]=1;color_t[123][1]=0;color_t[123][2]=0;
	color_t[124][0]=0;color_t[124][1]=1;color_t[124][2]=1;

	//render_stl_file(filenames_traces_ply,125,color_t);

	for(int counter=0; counter<pc; counter++)
	{
		delete [] filenames_traces_ply[counter];
	}

//	render_stl_file(stl_filenames,1,col);
	
	//vtkSmartPointer<vtkImageData> im=allocateImage(dimensions);
	//
	//for(int co=0; co<4; co++)
	//{
	//	char buff[1024];
	//	sprintf(buff, "transformed_100um2percent25x%dunmixed04.npts",co+1);
	//	loadImageFromFile(buff,im,dimensions);
	//}
	////vtkImageWriter* writer = vtkImageWriter::New();
	////writer->SetFilePattern("output%d.bmp");
	//
	//vtkSmartPointer<vtkContourFilter> contourf = vtkSmartPointer<vtkContourFilter>::New();
	//contourf->SetInput(im);
	//contourf->SetValue(0,0.5);
	//vtkSmartPointer<vtkSmoothPolyDataFilter> smoothf = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
	//smoothf->SetInput(contourf->GetOutput());
	//smoothf->SetRelaxationFactor(0.1);
	//smoothf->SetNumberOfIterations(20);
 //   mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	//mapper->SetInput(smoothf->GetOutput());
	//printf("Begin updating mapper\n");
	//mapper->Update();
	//vtkSmartPointer<vtkSTLWriter> stlwriter = vtkSmartPointer<vtkSTLWriter>::New();
	//stlwriter->SetInput(smoothf->GetOutput());
	//stlwriter->SetFileName("vessel.stl");
	//stlwriter->SetFileTypeToBinary();
	//printf("Begin Writing stl file\n");
	//stlwriter->Write();
	//actor= vtkSmartPointer<vtkLODActor>::New();
	//actor->SetMapper(mapper);
	////actor->SetScale(1,1,3);
	//ren1->AddActor(actor);
	//ren1->SetBackground(0,0,0);

	//vtkRenderWindow* rwin = vtkRenderWindow::New();
	//rwin->AddRenderer(ren1);

	//vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::New();
	//iren->SetRenderWindow(rwin);

	//rwin->Render();

	//iren->Start();

	//writer->SetFileName("output_montage.raw");
	//
	//writer->SetInput(im);
	//writer->SetFileDimensionality(3);
	//writer->Write();
	//writer->Delete();
	
	return 0;

	//for(int counter=0; counter<20; counter++)
	//{
	//	reader= vtkSmartPointer<vtkOBJReader>::New();
	//	reader->SetFileName(filenames[counter]);
	//	reader->Update();
	//	mapper= vtkSmartPointer<vtkPolyDataMapper>::New();
	//	mapper->SetInput(reader->GetOutput());
	//	mapper->Update();
	//	printf("I updated the mapper\n");
	//	actor= vtkSmartPointer<vtkLODActor>::New();
	//	actor->SetMapper(mapper);
	//	printf("Mapper assigned to the actor\n");
	//	actor->GetProperty()->SetColor(colors[counter][0],colors[counter][1],colors[counter][2]);
	//	assembly->AddPart(actor);
	//	printf("Added File %s\n",filenames[counter]);
	//}

	//for(int counter=0; counter<8; counter++)
	//{
	//	polydata = GetLines(filename_tracefiles[counter]);
	//	printf("Counter = %d\n",counter);
	//	mapper= vtkSmartPointer<vtkPolyDataMapper>::New();
	//	printf("About to assign polydata to the mapper\n");
	//	mapper->SetInput(polydata);
	//	mapper->Update();
	//	printf("Update the mapper\n");
	//	actor= vtkSmartPointer<vtkLODActor>::New();
	//	actor->SetMapper(mapper);
	//	actor->GetProperty()->SetColor(trace_colors[counter%2][0],trace_colors[counter%2][1],trace_colors[counter%2][2]);
	//	assembly->AddPart(actor);
	//}

	//printf("created the assembly\n");
	//ren1->AddActor(assembly);
	//ren1->SetBackground(0,0,0);

	//vtkRenderWindow* rwin = vtkRenderWindow::New();
	//rwin->AddRenderer(ren1);

	//vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::New();
	//iren->SetRenderWindow(rwin);

	//rwin->Render();

	//iren->Start();
	//sphere->Delete();

	//ren1->Delete();
	//rwin->Delete();
	//iren->Delete();

	//return 0;
}
