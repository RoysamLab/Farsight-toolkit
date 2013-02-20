#include "helpers.h"


#include <itkVTKImageExport.h>
#include <itkVTKImageImport.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkProperty.h>
#include <vtkMarchingCubes.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkImageImport.h>
#include <vtkImageConstantPad.h>
#include <vtkRenderLargeImage.h>
#include <vtkTIFFWriter.h>
#include <vtkVectorText.h>
#include <vtkFollower.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkAppendPolyData.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkVolume.h>

#include <vtkVolumeProperty.h>
#include <vtkPiecewiseFunction.h>
#include <vtkColorTransferFunction.h>
#include <vtkOpenGLVolumeTextureMapper2D.h>
#include <vtkOpenGLVolumeTextureMapper3D.h>
#include <vtkVolumeRayCastMapper.h>
#include <vtkVolumeRayCastCompositeFunction.h>
#include <vtkFixedPointVolumeRayCastMapper.h>
#include <vtkLightCollection.h>
#include <vtkLight.h>
#include <vtkCamera.h>
#include <vtkCubeSource.h>
#include <vtkPolyDataNormals.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper2D.h>
#include <vtkProperty2D.h>

#define OFFSET 70

using namespace helpers;
#if defined(_MSC_VER)
#pragma warning(disable: 4996)
#endif
double start_t,end_t,diff_t;
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

#define MAX_TIME 200
#define MAX_TAGS 4
#define MAX_LABEL 10000
#define PAUSE {printf("%d:>",__LINE__);scanf("%*d");}

bool file_exists(char *filename)
{
	FILE * fp = fopen(filename,"r");
	if(fp!=NULL)
	{
		fclose(fp);
		return true;
	}
	return false;
}


template <typename T>
typename T::Pointer readImage(const char *filename)
{
	printf("Reading %s ... \n",filename);
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
	printf("Done\n");
	return reader->GetOutput();

}
template <typename T>
int writeImage(typename T::Pointer im, const char* filename)
{
	printf("Writing %s ... \n",filename);
	typedef typename itk::ImageFileWriter<T> WriterType;

	typename WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(filename);
	writer->SetInput(im);
	try
	{
		writer->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "ExceptionObject caught!" <<std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
void getArrayFromStdVector(std::vector<FeaturesType> &f, FeaturesType	*&farray)
{
	farray = new FeaturesType[f.size()];
	for(unsigned int counter=0; counter<f.size(); counter++)
	{
		farray[counter]=f[counter];
	}
}

void getStdVectorFromArray(FeaturesType *farray, int n,std::vector<FeaturesType> &f)
{
	f.reserve(n);
	for(int counter=0; counter<n; counter++)
	{
		f.push_back(farray[counter]);
	}
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
		carray[counter].sx=60000;carray[counter].sy=60000;carray[counter].sz=60000;
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
	// create a tiny image of maximum size


	//appendfilter->UserManagedInputsOn();
	//appendfilter->SetNumberOfInputs(max1);
	vtkSmartPointer<vtkAppendPolyData> appendfilter = vtkSmartPointer<vtkAppendPolyData>::New();

	for(int counter=1; counter<=max1; counter++)
	{

			if(carray[counter].sx>=carray[counter].ex)
				continue;
			if(counter==169)
				printf("%d %d %d %d %d %d\n",carray[counter].sx,carray[counter].ex,carray[counter].sy,carray[counter].ey,carray[counter].sz,carray[counter].ez);
			InputImageType::Pointer t = getEmpty(wx,wy,wz);
		//	printf("Maximum tiny image size I need is [%d %d %d]\n",wx,wy,wz);

			InputImageType::SizeType size;
			InputImageType::RegionType region;
			index.Fill(1);

			region.SetIndex(index);
			region.SetSize(size);


			LabelImageType::SizeType lsize;
			LabelImageType::IndexType lindex;
			LabelImageType::RegionType lregion;


			ExportFilterType::Pointer itkexporter = ExportFilterType::New();
			itkexporter->SetInput(t);
			vtkSmartPointer<vtkImageImport> vtkimporter = vtkSmartPointer<vtkImageImport>::New();
			ConnectPipelines(itkexporter,(vtkImageImport *)vtkimporter);
			vtkSmartPointer<vtkMarchingCubes> contourf = vtkSmartPointer<vtkMarchingCubes>::New();
			contourf->SetInput(vtkimporter->GetOutput());
			contourf->SetValue(0,127);
			contourf->ComputeNormalsOff();
			contourf->ComputeScalarsOff();
			contourf->ComputeGradientsOff();
			vtkSmartPointer<vtkSmoothPolyDataFilter> smoothf = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
			smoothf->SetInput(contourf->GetOutput());
			smoothf->SetRelaxationFactor(0.2);
			smoothf->SetNumberOfIterations(20);



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

				vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
	transform->PostMultiply();
	vtkSmartPointer<vtkTransformPolyDataFilter> tf = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	tf->SetTransform(transform);
			transform->Identity();
			transform->Translate(carray[counter].sx-1,carray[counter].sy-1,carray[counter].sz-1);

			tf->SetInput(smoothf->GetOutput());
		
		tf->Update();
	//	tf->GetOutput()->Print(std::cout);
		
		appendfilter->AddInput(tf->GetOutput());
	
		//appendfilter->SetInputByNumber(counter-1,tf->GetOutput());
	//	appendfilter->Update();
	//	appendfilter->GetOutput()->Print(std::cout);
		printf("Completed %d/%d\r",counter,max1);
	//	scanf("%*d");
	}

	appendfilter->Update();
	delete [] carray;
	vtkSmartPointer<vtkPolyData> out = appendfilter->GetOutput();
	return out;
}

std::vector<vtkSmartPointer<vtkPolyData> > getIndividualPolyData(LabelImageType::Pointer label)
{
	LabelIteratorType liter = LabelIteratorType(label,label->GetLargestPossibleRegion());
	liter.GoToBegin();

	std::vector<vtkSmartPointer<vtkPolyData> > cells;
	//find the maximum number of cells
	unsigned short max1 = 0;
	for(liter.GoToBegin();!liter.IsAtEnd();++liter)
		max1 = MAX(max1,liter.Get());

	//find all the cubes in which cells lie
	cubecoord* carray = new cubecoord[max1+1];
	for(int counter=0; counter<=max1; counter++)
	{
		carray[counter].sx=60000;carray[counter].sy=60000;carray[counter].sz=60000;
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
	// create a tiny image of maximum size


	//appendfilter->UserManagedInputsOn();
	//appendfilter->SetNumberOfInputs(max1);
	//vtkSmartPointer<vtkAppendPolyData> appendfilter = vtkSmartPointer<vtkAppendPolyData>::New();

	for(int counter=1; counter<=max1; counter++)
	{

			if(carray[counter].sx>=carray[counter].ex)
				continue;
			InputImageType::Pointer t = getEmpty(wx,wy,wz);
		//	printf("Maximum tiny image size I need is [%d %d %d]\n",wx,wy,wz);

			InputImageType::SizeType size;
			InputImageType::RegionType region;
			index.Fill(1);

			region.SetIndex(index);
			region.SetSize(size);


			LabelImageType::SizeType lsize;
			LabelImageType::IndexType lindex;
			LabelImageType::RegionType lregion;


			ExportFilterType::Pointer itkexporter = ExportFilterType::New();
			itkexporter->SetInput(t);
			vtkSmartPointer<vtkImageImport> vtkimporter = vtkSmartPointer<vtkImageImport>::New();
			ConnectPipelines(itkexporter,(vtkImageImport *)vtkimporter);
			vtkSmartPointer<vtkMarchingCubes> contourf = vtkSmartPointer<vtkMarchingCubes>::New();
			contourf->SetInput(vtkimporter->GetOutput());
			contourf->SetValue(0,127);
			contourf->ComputeNormalsOff();
			contourf->ComputeScalarsOff();
			contourf->ComputeGradientsOff();
			vtkSmartPointer<vtkSmoothPolyDataFilter> smoothf = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
			smoothf->SetInput(contourf->GetOutput());
			smoothf->SetRelaxationFactor(0.2);
			smoothf->SetNumberOfIterations(20);



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

				vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
	transform->PostMultiply();
	vtkSmartPointer<vtkTransformPolyDataFilter> tf = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	tf->SetTransform(transform);
			transform->Identity();
			transform->Translate(carray[counter].sx-1,carray[counter].sy-1,carray[counter].sz-1);

			tf->SetInput(smoothf->GetOutput());
		
		tf->Update();
	//	tf->GetOutput()->Print(std::cout);
		vtkSmartPointer<vtkPolyData> out1 = vtkSmartPointer<vtkPolyData>::New();
		out1->DeepCopy(tf->GetOutput());
		cells.push_back(out1);
		//appendfilter->AddInput(tf->GetOutput());
	
		//appendfilter->SetInputByNumber(counter-1,tf->GetOutput());
	//	appendfilter->Update();
	//	appendfilter->GetOutput()->Print(std::cout);
		printf("Completed %d/%d\r",counter,max1);
	//	scanf("%*d");
	}

	//appendfilter->Update();
	delete [] carray;
	//vtkSmartPointer<vtkPolyData> out = appendfilter->GetOutput();
	return cells;
}
std::vector<vtkSmartPointer<vtkTextActor> > getTextActors(std::vector<FeaturesType> f[][MAX_TAGS],const int current_time)
{
	std::vector<vtkSmartPointer<vtkTextActor> > avec;
	char buff[100];
	//FIXME
	for(int tag_counter =1; tag_counter<=2; tag_counter++)
	{
		for(int time_counter=current_time;time_counter<=current_time; time_counter++)
		{
			for(int counter=0; counter< f[time_counter-1][tag_counter].size(); counter++)
			{

				vtkSmartPointer<vtkTextActor> textact = vtkSmartPointer<vtkTextActor>::New();
			//	textact->ScaledTextOn();
				sprintf(buff,"%d",f[time_counter-1][tag_counter][counter].num);
				textact->SetInput(buff);
				textact->GetPositionCoordinate()->SetCoordinateSystemToWorld();
				textact->GetPositionCoordinate()->SetValue(f[time_counter-1][tag_counter][counter].Centroid[0],f[time_counter-1][tag_counter][counter].Centroid[1],f[time_counter-1][tag_counter][counter].Centroid[2]*5);
				textact->GetPosition2Coordinate()->SetCoordinateSystemToWorld();
				textact->GetPosition2Coordinate()->SetValue(f[time_counter-1][tag_counter][counter].Centroid[0]+10,f[time_counter-1][tag_counter][counter].Centroid[1]+10,f[time_counter-1][tag_counter][counter].Centroid[2]*5);
				if(tag_counter==1)
				{
					textact->GetTextProperty()->SetColor(1,1,0);
				}
				else
				{
					textact->GetTextProperty()->SetColor(0,1,1);
				}
				textact->GetTextProperty()->SetJustificationToCentered();
				textact->GetTextProperty()->SetFontSize(16);
				avec.push_back(textact);

			//	vtkSmartPointer<vtkVectorText> atext = vtkSmartPointer<vtkVectorText>::New();
			//	sprintf(buff,"%d",f[time_counter-1][tag_counter][counter].num);
			//	atext->SetText(buff);

				//vtkSmartPointer<vtkPolyDataMapper> amapper = vtkSmartPointer<vtkPolyDataMapper>::New();
				//amapper->SetInput(atext->GetOutput());

				//vtkSmartPointer<vtkFollower> textActor = vtkSmartPointer<vtkFollower>::New();
				//textActor->SetMapper(amapper);
				//textActor->SetScale(8,8,8);
				//textActor->AddPosition(f[time_counter-1][tag_counter][counter].x,f[time_counter-1][tag_counter][counter].y,f[time_counter-1][tag_counter][counter].z*2+5);
				//textActor->GetProperty()->SetColor(1,0,0);
				//avec.push_back(textActor);
			}
		}
	}
	return avec;
}


//#define CACHE_PREFIX "D:/ucb dataset/output/ena/cache"

bool compare(FeaturesType a, FeaturesType b)
{
	return a.time<b.time;
}


typedef vtkSmartPointer<vtkPolyData> SP_PDM;
void renderPolyData(std::vector<SP_PDM> vec,std::vector<vtkSmartPointer<vtkTextActor> > text_actors,int timenum)
{
	printf("Entered renderPolyData\n");
	vtkSmartPointer<vtkRenderer> ren = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renwin = vtkSmartPointer<vtkRenderWindow>::New();

	//float vtkcolor[][3]={0,1,0,
	//	0,154/255.0,25/255.0,
	//	207/255.0,141/255.0,0,
	//	1,0,0};
 float vtkcolor[4][3]={0,1,0,
		0,154/255.0,25/255.0,
		207/255.0,141/255.0,0,
		1,0,0};


	//vec[0]->Print(cout);
	
	for(int counter=0; counter<vec.size(); counter++)
	{
		vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInput(vec[counter]);
		mapper->ImmediateModeRenderingOff();
		vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		actor->GetProperty()->SetColor(vtkcolor[counter][0],vtkcolor[counter][1],vtkcolor[counter][2]);
		actor->SetScale(0.96,0.96,4);
		ren->AddActor(actor);
	}

	printf("Out of for loop\n");
	for(int counter=0; counter<text_actors.size(); counter++)
	{
		ren->AddActor(text_actors[counter]);
		
	}

	
	//vtkSmartPointer<vtkTextActor> textact = vtkSmartPointer<vtkTextActor>::New();
	////textact->ScaledTextOn();
	//textact->SetInput("123");

	//textact->GetPositionCoordinate()->SetCoordinateSystemToWorld();
	//textact->GetPositionCoordinate()->SetValue(240,200,15);
	//textact->GetPosition2Coordinate()->SetCoordinateSystemToWorld();
	//textact->GetPosition2Coordinate()->SetValue(240+5,200+5,15);
	//ren->AddActor(textact);

	ren->SetBackground(1,1,1);

	renwin->AddRenderer(ren);
	renwin->SetSize(700,700);

	
	vtkSmartPointer<vtkRenderWindowInteractor> renwininteract = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renwininteract->SetRenderWindow(renwin);
	renwin->Render();

	vtkRenderLargeImage *renderLarge=vtkRenderLargeImage::New();
	renderLarge->SetInput(ren);
	renderLarge->SetMagnification(1);
	vtkTIFFWriter *writer = vtkTIFFWriter::New();
	writer->SetInput(renderLarge->GetOutput());
	char buf_rand[1000];
	//printf("rendering Rendering_%d.tif\n",timenum+OFFSET);
	//sprintf(buf_rand,"Rendering_%d.tif",timenum+OFFSET);
	//writer->SetFileName(buf_rand);
	//writer->Write();

	/*for(int counter=0; counter<text_actors.size(); counter++)
	{
		text_actors[counter]->SetCamera(ren->GetActiveCamera());
	}*/
	
	
	

     renwininteract->Start();

}

void hsv_to_rgb(float hsv[], float rgb[])
{

	float h =  hsv[0];
	float s = hsv[1];
	float v = hsv[2];
  int h_i  = (h*6);
  float f = h*6 - h_i;
  float p = v * (1 - s);
  float q = v * (1 - f*s);
  float t = v * (1 - (1 - f) * s);
  switch(h_i)
  {
  case 0:
	  rgb[0] = v; rgb[1] = t; rgb[2] = p;break;
  case 1:
	  rgb[0] = q; rgb[1] = v; rgb[2] = p;break;
  case 2:
	  rgb[0] = p; rgb[1] = v; rgb[2] = t;break;
  case 3:
	  rgb[0] = p; rgb[1] = q; rgb[2] = v;break;
  case 4:
	  rgb[0] = t; rgb[1] = p; rgb[2] = v;break;
  case 5:
  default:
	  rgb[0] = v; rgb[1] = p; rgb[2] = q;break;
  }
}

//#define MAX_TRACKS_NUM 2000
void get_rgb_for_num(int num,float rgb[3])
{
	float golden_conjugate = 0.61803;
	//num = num % MAX_TRACKS_NUM;
	float hsv[3];
	srand(747);
	hsv[0] = fmod((num)*213*golden_conjugate,1);
	srand(1024);
	for(int counter = 0; counter< num; counter++)
		hsv[1] = rand()*1.0/RAND_MAX*0.3+0.7;
	//hsv[2] = fmod(num*1231*golden_conjugate,1)*0.45+0.5;
	srand(2343);
	for(int counter = 0; counter < num; counter++)
		hsv[2] = rand()*1.0/RAND_MAX*0.45+0.5;
	hsv_to_rgb(hsv,rgb);
}

void renderTrackingData(std::vector<SP_PDM> vec,std::vector<FeaturesType> &fv, int timenum, int channel)
{
	printf("Entered renderPolyData\n");
	vtkSmartPointer<vtkRenderer> ren = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renwin = vtkSmartPointer<vtkRenderWindow>::New();

	//float vtkcolor[][3]={0,1,0,
	//	0,154/255.0,25/255.0,
	//	207/255.0,141/255.0,0,
	//	1,0,0};
 float vtkcolor[4][3]={0,1,0,
		0,154/255.0,25/255.0,
		207/255.0,141/255.0,0,
		1,0,0};


	//vec[0]->Print(cout);
	
	double golden_conjugate = 0.618033988749895;
	double h = 5*rand()*1.0/RAND_MAX;
	float s = 1.0;
	float v = rand()*1.0/RAND_MAX;
	float rgb[3],hsv[3];
	srand(1023);
	if(vec.size()!=fv.size())
	{
		printf("Something is wrong...\n");
		scanf("%*d");
	}
	for(int counter=0; counter<vec.size(); counter++)
	{
		vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		vtkSmartPointer<vtkPolyDataNormals> polynorm = vtkSmartPointer<vtkPolyDataNormals>::New();
		polynorm->SetInput(vec[counter]);
		mapper->SetInput(polynorm->GetOutput());
		//mapper->SetInput(vec[counter]);
		mapper->ImmediateModeRenderingOff();
		vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		//h = h + golden_conjugate;
		//h = fmod(h,1);
		//hsv[0] = h;
		//hsv[0] = rand()*1.0/RAND_MAX;
		//s = rand()*1.0/RAND_MAX*0.5 + 0.5;
		//hsv[1] = s;
		//v = v + golden_conjugate;
		//v = fmod(v,1);
		//v = rand()*1.0/RAND_MAX;
		//hsv[2] = v*0.45+0.5;
		//printf("%0.3f %0.3f %0.3f\n",hsv[0],hsv[1],hsv[2]);
		//hsv_to_rgb(hsv,rgb);
		get_rgb_for_num(fv[counter].num,rgb);
		printf("%d = %f %f %f\n",fv[counter].num,rgb[0],rgb[1],rgb[2]);
		actor->GetProperty()->SetColor(rgb[0],rgb[1],rgb[2]);
		//actor->GetProperty()->Print(std::cout);
		actor->GetProperty()->SetAmbient(0.4);
		//actor->GetProperty()->SetDiffuse(0);
		//actor->GetProperty()->SetAmbientColor(1,1,1);
		//actor->GetProperty()->SetDiffuseColor(1,1,1);
		actor->SetScale(0.96,0.96,4);
		ren->AddActor(actor);
	}

	printf("Out of for loop\n");
	vtkSmartPointer<vtkCubeSource> cbsource = vtkSmartPointer<vtkCubeSource>::New();
	cbsource->SetBounds(-1,513,-1,513,-1,45);
	vtkSmartPointer<vtkPolyDataMapper> polymap = vtkSmartPointer<vtkPolyDataMapper>::New();
	polymap->SetInput(cbsource->GetOutput());
	vtkSmartPointer <vtkActor> cbact = vtkSmartPointer<vtkActor>::New();
	cbact->SetMapper(polymap);
	cbact->GetProperty()->SetOpacity(0);
	ren->AddActor(cbact);
	
	//vtkSmartPointer<vtkTextActor> textact = vtkSmartPointer<vtkTextActor>::New();
	////textact->ScaledTextOn();
	//textact->SetInput("123");

	//textact->GetPositionCoordinate()->SetCoordinateSystemToWorld();
	//textact->GetPositionCoordinate()->SetValue(240,200,15);
	//textact->GetPosition2Coordinate()->SetCoordinateSystemToWorld();
	//textact->GetPosition2Coordinate()->SetValue(240+5,200+5,15);
	//ren->AddActor(textact);

	ren->SetBackground(0,0,0);
	//ren->TwoSidedLightingOn();
	//ren->MakeLight();
	//ren->Print(std::cout);
	//ren->GetActiveCamera()->Zoom(1.001);
	//ren->GetActiveCamera()->ParallelProjectionOn();
	ren->GetLights()->InitTraversal();
	
	//ren->GetLights()->GetNextItem()->Print(std::cout);

	renwin->AddRenderer(ren);
	//renwin->SetSize(700,700);

	renwin->FullScreenOn();
	vtkSmartPointer<vtkRenderWindowInteractor> renwininteract = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renwininteract->SetRenderWindow(renwin);
	renwin->Render();

	vtkRenderLargeImage *renderLarge=vtkRenderLargeImage::New();
	renderLarge->SetInput(ren);
	renderLarge->SetMagnification(3);
	vtkTIFFWriter *writer = vtkTIFFWriter::New();
	writer->SetInput(renderLarge->GetOutput());
	char buf_rand[1000];
	//printf("rendering Rendering_%d.tif\n",timenum+OFFSET);
	sprintf(buf_rand,"L:\\thesis_figures\\02102009_1455-624_harvard_tracking_renderings\\testTSeries-02102009-1455-624_ch%d_cycle%03d.tif",channel,timenum+1);
	writer->SetFileName(buf_rand);
	writer->Write();
	

	/*for(int counter=0; counter<text_actors.size(); counter++)
	{
		text_actors[counter]->SetCamera(ren->GetActiveCamera());
	}*/
	
	
	

    // renwininteract->Start();

}

void renderTrackingData_withSegmentation(std::vector<SP_PDM> vec1,std::vector<FeaturesType> &fv1,int ch1,std::vector<SP_PDM> vec2,std::vector<FeaturesType> &fv2, int ch2,int timenum)
{
	printf("Entered renderPolyData\n");
	vtkSmartPointer<vtkRenderer> ren = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renwin = vtkSmartPointer<vtkRenderWindow>::New();

	//float vtkcolor[][3]={0,1,0,
	//	0,154/255.0,25/255.0,
	//	207/255.0,141/255.0,0,
	//	1,0,0};
 float vtkcolor[4][3]={0,1,0,
		0,154/255.0,25/255.0,
		207/255.0,141/255.0,0,
		1,0,0};


	//vec[0]->Print(cout);
 float sp[3] = {1.15,1.15,4.0};
	double golden_conjugate = 0.618033988749895;
	double h = 5*rand()*1.0/RAND_MAX;
	float s = 1.0;
	float v = rand()*1.0/RAND_MAX;
	float rgb[3],hsv[3];
	srand(1023);
	if(vec1.size()!=fv1.size())
	{
		printf("Something is wrong...\n");
		scanf("%*d");
	}
	for(int counter=0; counter<vec1.size(); counter++)
	{
		vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		vtkSmartPointer<vtkPolyDataNormals> polynorm = vtkSmartPointer<vtkPolyDataNormals>::New();
		polynorm->SetInput(vec1[counter]);
		mapper->SetInput(polynorm->GetOutput());
		//mapper->SetInput(vec[counter]);
		mapper->ImmediateModeRenderingOff();
		vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		//h = h + golden_conjugate;
		//h = fmod(h,1);
		//hsv[0] = h;
		//hsv[0] = rand()*1.0/RAND_MAX;
		//s = rand()*1.0/RAND_MAX*0.5 + 0.5;
		//hsv[1] = s;
		//v = v + golden_conjugate;
		//v = fmod(v,1);
		//v = rand()*1.0/RAND_MAX;
		//hsv[2] = v*0.45+0.5;
		//printf("%0.3f %0.3f %0.3f\n",hsv[0],hsv[1],hsv[2]);
		//hsv_to_rgb(hsv,rgb);
		get_rgb_for_num(fv1[counter].num,rgb);
		printf("%d = %f %f %f\n",fv1[counter].num,rgb[0],rgb[1],rgb[2]);
		actor->GetProperty()->SetColor(rgb[0],rgb[1],rgb[2]);
		//actor->GetProperty()->Print(std::cout);
		actor->GetProperty()->SetAmbient(0.4);
		//actor->GetProperty()->SetDiffuse(0);
		//actor->GetProperty()->SetAmbientColor(1,1,1);
		//actor->GetProperty()->SetDiffuseColor(1,1,1);
		actor->SetScale(sp[0],sp[1],sp[2]);
		ren->AddActor(actor);
	}
	for(int counter=0; counter<vec2.size(); counter++)
	{
		vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		vtkSmartPointer<vtkPolyDataNormals> polynorm = vtkSmartPointer<vtkPolyDataNormals>::New();
		polynorm->SetInput(vec2[counter]);
		mapper->SetInput(polynorm->GetOutput());
		//mapper->SetInput(vec[counter]);
		mapper->ImmediateModeRenderingOff();
		vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		h = h + golden_conjugate;
		h = fmod(h,1);
		hsv[0] = h;
		//hsv[0] = rand()*1.0/RAND_MAX;
		s = rand()*1.0/RAND_MAX*0.5 + 0.5;
		hsv[1] = s;
		//v = v + golden_conjugate;
		//v = fmod(v,1);
		v = rand()*1.0/RAND_MAX;
		hsv[2] = v*0.45+0.5;
		//printf("%0.3f %0.3f %0.3f\n",hsv[0],hsv[1],hsv[2]);
		//hsv_to_rgb(hsv,rgb);
		get_rgb_for_num(fv2[counter].num,rgb);
		//printf("%d = %f %f %f\n",fv[counter].num,rgb[0],rgb[1],rgb[2]);
		//actor->GetProperty()->SetColor(rgb[0],rgb[1],rgb[2]);
		actor->GetProperty()->SetColor(0.7,0.7,0.7);
		//actor->GetProperty()->Print(std::cout);
		actor->GetProperty()->SetAmbient(0.4);
		//actor->GetProperty()->SetDiffuse(0);
		//actor->GetProperty()->SetAmbientColor(1,1,1);
		//actor->GetProperty()->SetDiffuseColor(1,1,1);
		actor->SetScale(sp[0],sp[1],sp[2]);

		//vtkSmartPointer<vtkSphereSource> spsource = vtkSmartPointer<vtkSphereSource>::New();
		//spsource->SetRadius(4);
		//spsource->SetPhiResolution(20);
		//spsource->SetThetaResolution(20);
		////spsource->SetCenter(fv[counter].Centroid[0],fv[counter].Centroid[1],40);
		//vtkSmartPointer<vtkPolyDataMapper2D> spmap = vtkSmartPointer<vtkPolyDataMapper2D>::New();
		//spmap->SetInput(spsource->GetOutput());
		//vtkSmartPointer<vtkActor2D> spact = vtkSmartPointer<vtkActor2D>::New();
		//spact->SetMapper(spmap);
		//spact->GetProperty()->SetColor(1,0,0);
		//spact->GetPositionCoordinate()->SetCoordinateSystemToWorld();
		//spact->GetPositionCoordinate()->SetValue(fv2[counter].Centroid[0]*sp[0],fv2[counter].Centroid[1]*sp[1],fv2[counter].Centroid[2]*sp[2]);
		//spact->GetPosition2Coordinate()->SetCoordinateSystemToWorld();
		//spact->GetPosition2Coordinate()->SetValue(fv2[counter].Centroid[0]*sp[0]+10,fv2[counter].Centroid[1]*sp[1]+10,fv2[counter].Centroid[2]*sp[2]);
		////spact->SetScale(0.96,0.96);
		//ren->AddActor(spact);
		ren->AddActor(actor);
	}


	printf("Out of for loop\n");
	vtkSmartPointer<vtkCubeSource> cbsource = vtkSmartPointer<vtkCubeSource>::New();
	cbsource->SetBounds(-1,513,-1,513,-1,45);
	vtkSmartPointer<vtkPolyDataMapper> polymap = vtkSmartPointer<vtkPolyDataMapper>::New();
	polymap->SetInput(cbsource->GetOutput());
	vtkSmartPointer <vtkActor> cbact = vtkSmartPointer<vtkActor>::New();
	cbact->SetMapper(polymap);
	cbact->GetProperty()->SetOpacity(0);
	ren->AddActor(cbact);
	
	//vtkSmartPointer<vtkTextActor> textact = vtkSmartPointer<vtkTextActor>::New();
	////textact->ScaledTextOn();
	//textact->SetInput("123");

	//textact->GetPositionCoordinate()->SetCoordinateSystemToWorld();
	//textact->GetPositionCoordinate()->SetValue(240,200,15);
	//textact->GetPosition2Coordinate()->SetCoordinateSystemToWorld();
	//textact->GetPosition2Coordinate()->SetValue(240+5,200+5,15);
	//ren->AddActor(textact);

	ren->SetBackground(0,0,0);
	//ren->TwoSidedLightingOn();
	//ren->MakeLight();
	//ren->Print(std::cout);
	//ren->GetActiveCamera()->Zoom(1.001);
	//ren->GetActiveCamera()->ParallelProjectionOn();
	ren->GetLights()->InitTraversal();
	
	//ren->GetLights()->GetNextItem()->Print(std::cout);

	renwin->AddRenderer(ren);
	//renwin->SetSize(700,700);

	renwin->FullScreenOn();
	vtkSmartPointer<vtkRenderWindowInteractor> renwininteract = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renwininteract->SetRenderWindow(renwin);
	renwin->Render();

	vtkRenderLargeImage *renderLarge=vtkRenderLargeImage::New();
	renderLarge->SetInput(ren);
	renderLarge->SetMagnification(3);
	vtkTIFFWriter *writer = vtkTIFFWriter::New();
	writer->SetInput(renderLarge->GetOutput());
	char buf_rand[1000];
	//printf("rendering Rendering_%d.tif\n",timenum+OFFSET);
	sprintf(buf_rand,"L:\\thesis_figures\\harvard_tracking_renderings\\testTSeries-11172010-1053-017_cycle%03d.tif",timenum+1);
	writer->SetFileName(buf_rand);
	writer->Write();
	

	/*for(int counter=0; counter<text_actors.size(); counter++)
	{
		text_actors[counter]->SetCamera(ren->GetActiveCamera());
	}*/
	
	
	

     //renwininteract->Start();

}
void renderSegmentationData(std::vector<SP_PDM> vec,std::vector<FeaturesType> &fv, int timenum, int channel)
{
	printf("Entered renderPolyData\n");
	vtkSmartPointer<vtkRenderer> ren = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renwin = vtkSmartPointer<vtkRenderWindow>::New();

	//float vtkcolor[][3]={0,1,0,
	//	0,154/255.0,25/255.0,
	//	207/255.0,141/255.0,0,
	//	1,0,0};
 float vtkcolor[4][3]={0,1,0,
		0,154/255.0,25/255.0,
		207/255.0,141/255.0,0,
		1,0,0};


	//vec[0]->Print(cout);
	
	double golden_conjugate = 0.618033988749895;
	double h = rand()*1.0/RAND_MAX;
	float s = 1.0;
	float v = rand()*1.0/RAND_MAX;
	float rgb[3],hsv[3];
	if(vec.size()!=fv.size())
	{
		printf("Something is wrong...\n");
		scanf("%*d");
	}
	for(int counter=0; counter<vec.size(); counter++)
	{
		vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		vtkSmartPointer<vtkPolyDataNormals> polynorm = vtkSmartPointer<vtkPolyDataNormals>::New();
		polynorm->SetInput(vec[counter]);
		mapper->SetInput(polynorm->GetOutput());
		//mapper->SetInput(vec[counter]);
		mapper->ImmediateModeRenderingOff();
		vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		h = h + golden_conjugate;
		h = fmod(h,1);
		hsv[0] = h;
		//hsv[0] = rand()*1.0/RAND_MAX;
		s = rand()*1.0/RAND_MAX*0.5 + 0.5;
		hsv[1] = s;
		//v = v + golden_conjugate;
		//v = fmod(v,1);
		v = rand()*1.0/RAND_MAX;
		hsv[2] = v*0.45+0.5;
		//printf("%0.3f %0.3f %0.3f\n",hsv[0],hsv[1],hsv[2]);
		//hsv_to_rgb(hsv,rgb);
		get_rgb_for_num(fv[counter].num,rgb);
		//printf("%d = %f %f %f\n",fv[counter].num,rgb[0],rgb[1],rgb[2]);
		//actor->GetProperty()->SetColor(rgb[0],rgb[1],rgb[2]);
		actor->GetProperty()->SetColor(0,1,0);
		//actor->GetProperty()->Print(std::cout);
		actor->GetProperty()->SetAmbient(0.4);
		//actor->GetProperty()->SetDiffuse(0);
		//actor->GetProperty()->SetAmbientColor(1,1,1);
		//actor->GetProperty()->SetDiffuseColor(1,1,1);
		actor->SetScale(0.96,0.96,4);

		vtkSmartPointer<vtkSphereSource> spsource = vtkSmartPointer<vtkSphereSource>::New();
		spsource->SetRadius(4);
		spsource->SetPhiResolution(20);
		spsource->SetThetaResolution(20);
		//spsource->SetCenter(fv[counter].Centroid[0],fv[counter].Centroid[1],40);
		vtkSmartPointer<vtkPolyDataMapper2D> spmap = vtkSmartPointer<vtkPolyDataMapper2D>::New();
		spmap->SetInput(spsource->GetOutput());
		vtkSmartPointer<vtkActor2D> spact = vtkSmartPointer<vtkActor2D>::New();
		spact->SetMapper(spmap);
		spact->GetProperty()->SetColor(1,0,0);
		spact->GetPositionCoordinate()->SetCoordinateSystemToWorld();
		spact->GetPositionCoordinate()->SetValue(fv[counter].Centroid[0]*0.96,fv[counter].Centroid[1]*0.96,fv[counter].Centroid[2]*4);
		spact->GetPosition2Coordinate()->SetCoordinateSystemToWorld();
		spact->GetPosition2Coordinate()->SetValue(fv[counter].Centroid[0]*0.96+10,fv[counter].Centroid[1]*0.96+10,fv[counter].Centroid[2]*4);
		//spact->SetScale(0.96,0.96);
		ren->AddActor(spact);
		ren->AddActor(actor);
	}

	printf("Out of for loop\n");
	vtkSmartPointer<vtkCubeSource> cbsource = vtkSmartPointer<vtkCubeSource>::New();
	cbsource->SetBounds(-1,513,-1,513,-1,45);
	vtkSmartPointer<vtkPolyDataMapper> polymap = vtkSmartPointer<vtkPolyDataMapper>::New();
	polymap->SetInput(cbsource->GetOutput());
	vtkSmartPointer <vtkActor> cbact = vtkSmartPointer<vtkActor>::New();
	cbact->SetMapper(polymap);
	cbact->GetProperty()->SetOpacity(0);
	ren->AddActor(cbact);
	
	//vtkSmartPointer<vtkTextActor> textact = vtkSmartPointer<vtkTextActor>::New();
	////textact->ScaledTextOn();
	//textact->SetInput("123");

	//textact->GetPositionCoordinate()->SetCoordinateSystemToWorld();
	//textact->GetPositionCoordinate()->SetValue(240,200,15);
	//textact->GetPosition2Coordinate()->SetCoordinateSystemToWorld();
	//textact->GetPosition2Coordinate()->SetValue(240+5,200+5,15);
	//ren->AddActor(textact);

	ren->SetBackground(0,0,0);
	//ren->TwoSidedLightingOn();
	//ren->MakeLight();
	//ren->Print(std::cout);
	//ren->GetActiveCamera()->Zoom(1.001);
	//ren->GetActiveCamera()->ParallelProjectionOn();
	ren->GetLights()->InitTraversal();
	
	//ren->GetLights()->GetNextItem()->Print(std::cout);

	renwin->AddRenderer(ren);
	//renwin->SetSize(700,700);

	renwin->FullScreenOn();
	vtkSmartPointer<vtkRenderWindowInteractor> renwininteract = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renwininteract->SetRenderWindow(renwin);
	renwin->Render();

	vtkRenderLargeImage *renderLarge=vtkRenderLargeImage::New();
	renderLarge->SetInput(ren);
	renderLarge->SetMagnification(3);
	vtkTIFFWriter *writer = vtkTIFFWriter::New();
	writer->SetInput(renderLarge->GetOutput());
	char buf_rand[1000];
	//printf("rendering Rendering_%d.tif\n",timenum+OFFSET);
	sprintf(buf_rand,"L:\\thesis_figures\\02102009_1455-624_harvard_tracking_renderings\\Seg_testTSeries-02102009-1455-624_ch%d_cycle%03d.tif",channel,timenum+1);
	writer->SetFileName(buf_rand);
	writer->Write();
	

	/*for(int counter=0; counter<text_actors.size(); counter++)
	{
		text_actors[counter]->SetCamera(ren->GetActiveCamera());
	}*/
	
	
	

     //renwininteract->Start();

}
LabelImageType::Pointer getFlippedUD(LabelImageType::Pointer label)
{
	LabelImageType::Pointer out = LabelImageType::New();
	out->SetRegions(label->GetLargestPossibleRegion());
	out->Allocate();

	typedef itk::ImageRegionIteratorWithIndex<LabelImageType> IRWI;
	IRWI iter1(label,label->GetLargestPossibleRegion());

	LabelImageType::SizeType size = out->GetLargestPossibleRegion().GetSize();
	LabelImageType::IndexType index;
	for(iter1.GoToBegin();!iter1.IsAtEnd();++iter1)
	{
		index = iter1.GetIndex();
		index[1] = size[1] - index[1] - 1;
		out->SetPixel(index,iter1.Get());
	}
	return out;
}
void createTrackFeatures(std::vector<FeaturesType> fvector[MAX_TIME][MAX_TAGS], std::vector<ftk::TrackFeatures> &tfs, int c,int num_t)
{
	int max_track_num = 0;
	for(int t = 0; t< num_t; t++)
	{
		for(unsigned int counter=0; counter< fvector[t][c-1].size(); counter++)
		{
			max_track_num = MAX(max_track_num,fvector[t][c-1][counter].num);
		}
	}

	for(int counter=1; counter <= max_track_num; counter++)
	{
		ftk::TrackFeatures trackf;
		trackf.intrinsic_features.clear();
		for(int t = 0; t< num_t;t++)
		{
		for(unsigned int counter1 = 0; counter1 < fvector[t][c-1].size(); counter1++)
		{
			if(fvector[t][c-1][counter1].num == counter)
			{
				trackf.intrinsic_features.push_back(fvector[t][c-1][counter1]);
			}
		}
		}
		std::sort(trackf.intrinsic_features.begin(),trackf.intrinsic_features.end(),compare);
		tfs.push_back(trackf);
		//PRINTF("Added %d elements to tfs\n",counter);
	}
}
bool mysortpredicate(FeaturesType a, FeaturesType b)
{
	return a.num<b.num;
}
int main(int argc, char**argv)
{
	printf("Started\n");
	//int counter = 0;
  int pc = 1;
  int num_t = 1; //atoi(argv[pc++]);
  int num_channels = 1;//atoi(argv[pc++]);
 
//int c = 4;//atoi(argv[2]);

  //for (int counter = 0; counter< num_channels; counter++)
  //{
	 // sscanf(argv[pc++],"%f,%f,%f",&colors[counter][0],&colors[counter][1],&colors[counter][2]);
  //}

  //std::ifstream infile;

  //infile.open(argv[pc],std::ifstream::in);
  //std::vector<std::string> fnames;
  //while(1)
  //{
	 // char buff[1024];
	 //infile.getline(buff,1023);
	 //if(infile.eof())
		// break;
	 //fnames.push_back(buff);
	 //std::cout<<buff;
  //}
  //return 0;
  //int num_other_channels = 0;

  for(int c = 2; c<=2; c++)
  {
	  LabelImageType::Pointer segmented; // FIXME
	  InputImageType::Pointer images;
	  std::vector<FeaturesType> fvector;
	  int max_time = 30;
	  //#define BASE "C:\\Users\\Arun\\Research\\Tracking\\berkeley\\cache\\wF5p120507m1s5\\"
#define BASE "L:\\Tracking\\cache\\testTSeries-11172010-1053-017\\"
	  char buff[1024];

	  //	if(c<2)
	  //		sprintf(buff,BASE"labeled_tracks_wF5p120507m1s5_w%d_t%d.tif",c+1,counter+1);
	 /* for(int counter =0; counter < max_time ; counter++)
	  {
		  fvector.clear();
		  sprintf(buff,BASE"labeled_tracks_TSeries-11172010-1053-017_Cycle%03d_CurrentSettings_Ch%d.tif",counter+1,c);
		  segmented = getFlippedUD(readImage<LabelImageType>(buff));
		  InputImageType::Pointer tempim = readImage<InputImageType>(BASE"smoothed_TSeries-02102009-1455-624_Cycle001_CurrentSettings_Ch4.tif");
		  std::vector<SP_PDM> v;
		  getFeatureVectorsFarsight(segmented,tempim,fvector,counter,c);
		  v = getIndividualPolyData(segmented);
		  sort(fvector.begin(),fvector.end(),mysortpredicate);
		  renderTrackingData(v,fvector,counter,c);
	  }
	  for(int counter =0; counter < max_time ; counter++)
	  {
		  fvector.clear();
		  sprintf(buff,BASE"clabeled_TSeries-02102009-1455-624_Cycle%03d_CurrentSettings_Ch%d.tif",counter+1,c);
		  segmented = getFlippedUD(readImage<LabelImageType>(buff));
		  InputImageType::Pointer tempim = readImage<InputImageType>(BASE"smoothed_TSeries-02102009-1455-624_Cycle001_CurrentSettings_Ch4.tif");
		  std::vector<SP_PDM> v;
		  getFeatureVectorsFarsight(segmented,tempim,fvector,counter,c);
		  v = getIndividualPolyData(segmented);
		  sort(fvector.begin(),fvector.end(),mysortpredicate);
		  renderSegmentationData(v,fvector,counter,c);
	  }*/
  }
  // For simultaneous rendering of segmented DC and tracked T cells
  {
	  LabelImageType::Pointer segmented1,segmented2; 
	  InputImageType::Pointer images;
	  std::vector<FeaturesType> fvector1,fvector2;
	  int max_time = 50;
	  int track_c = 2;
	  int dc_c = 3;
	  char buff[1024];
	  for(int counter =0; counter < max_time ; counter++)
	  {
		  fvector1.clear();
		  fvector2.clear();
		  sprintf(buff,BASE"labeled_tracks_TSeries-11172010-1053-017_Cycle%03d_CurrentSettings_Ch%d.tif",counter+1,track_c);
		  segmented1 = getFlippedUD(readImage<LabelImageType>(buff));
		  sprintf(buff,BASE"labeled_TSeries-11172010-1053-017_Cycle%03d_CurrentSettings_Ch%d.tif",counter+1,dc_c);
		  segmented2 = getFlippedUD(readImage<LabelImageType>(buff));
		  InputImageType::Pointer tempim = readImage<InputImageType>(BASE"unmixed_TSeries-11172010-1053-017_Cycle001_CurrentSettings_Ch2.tif");
		  std::vector<SP_PDM> v1,v2;
		  getFeatureVectorsFarsight(segmented1,tempim,fvector1,counter,track_c);
		  getFeatureVectorsFarsight(segmented2,tempim,fvector2,counter,dc_c);
		  v1 = getIndividualPolyData(segmented1);
		  v2 = getIndividualPolyData(segmented2);
		  sort(fvector1.begin(),fvector1.end(),mysortpredicate);
		  sort(fvector2.begin(),fvector2.end(),mysortpredicate);
		  renderTrackingData_withSegmentation(v1,fvector1,track_c,v2,fvector2,dc_c,counter);
	  }
  }
  return 0;
}
int old_main(int argc, char **argv)
{
	//ST();

	printf("Started\n");
	int c=1;
	//int counter = 0;
  int pc = 1;
  int num_t = 1; //atoi(argv[pc++]);
  int num_channels = 4;//atoi(argv[pc++]);
 


  //for (int counter = 0; counter< num_channels; counter++)
  //{
	 // sscanf(argv[pc++],"%f,%f,%f",&colors[counter][0],&colors[counter][1],&colors[counter][2]);
  //}

  //std::ifstream infile;

  //infile.open(argv[pc],std::ifstream::in);
  //std::vector<std::string> fnames;
  //while(1)
  //{
	 // char buff[1024];
	 //infile.getline(buff,1023);
	 //if(infile.eof())
		// break;
	 //fnames.push_back(buff);
	 //std::cout<<buff;
  //}
  //return 0;
  //int num_other_channels = 0;
 
	LabelImageType::Pointer segmented[MAX_TIME][MAX_TAGS]; // FIXME
	InputImageType::Pointer images[MAX_TIME][MAX_TAGS];
	std::vector<FeaturesType> fvector[MAX_TIME][MAX_TAGS];
	#define BASE "C:\\Users\\Arun\\Research\\Tracking\\berkeley\\cache\\wF5p120507m1s5\\"
	for(int c = 0; c< num_channels; c++)
	{
		for(int counter=0; counter< num_t; counter++)
		{
			char buff[1024];
			if(c<2)
				sprintf(buff,BASE"labeled_tracks_wF5p120507m1s5_w%d_t%d.tif",c+1,counter+1);
			else
				sprintf(buff,BASE"labeled_wF5p120507m1s5_w%d_t%d.tif",c+1,counter+1);
			segmented[counter][c] = getFlippedUD(readImage<LabelImageType>(buff));
		}
	}
	InputImageType::Pointer tempim = readImage<InputImageType>(BASE"unmixed_wF5p120507m1s5_w1_t1.tif");
	for(int counter=0; counter< num_t; counter++)
	{
		std::vector<SP_PDM> v;
		for(int c = 1; c<=num_channels; c++)
		{
			getFeatureVectorsFarsight(segmented[counter][c-1],tempim,fvector[counter][c-1],counter,c);
			v.push_back(getVTKPolyDataPrecise(segmented[counter][c-1]));
		}
		renderPolyData(v,getTextActors(fvector,counter),counter);
	}
	//for(int c=1;c<(num_channels+1);c++)
	
		//int c=1;
	//for(int t = 0; t< num_t ; t++)
	//{
	//	images[t][c-1]=readImage<InputImageType>(argv[pc++]);
	//}
 

  return 0;
  }
