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

std::vector<vtkSmartPointer<vtkTextActor>> getTextActors(std::vector<FeaturesType> f[][MAX_TAGS],const int current_time)
{
	std::vector<vtkSmartPointer<vtkTextActor>> avec;
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
void renderPolyData(std::vector<SP_PDM> vec,std::vector<vtkSmartPointer<vtkTextActor>> text_actors,int timenum)
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
		actor->SetScale(1,1,4);
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

	ren->SetBackground(0,0,0);

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
	printf("rendering Rendering_%d.tif\n",timenum+OFFSET);
	sprintf(buf_rand,"Rendering_%d.tif",timenum+OFFSET);
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
int main(int argc, char **argv)
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
