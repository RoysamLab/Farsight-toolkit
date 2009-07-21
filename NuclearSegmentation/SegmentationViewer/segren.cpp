#include "ftkNuclearSegmentation.h"
#include <direct.h>

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
#include <vtkPolyDataSource.h>
#include <vtkPolyDataNormals.h>
#include <vtkDecimatePro.h>
#include <vtkDiscreteMarchingCubes.h>
#include <vtkQuadricDecimation.h>
#include <vtkQuadricClustering.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkAppendPolyData.h>
#include <vtkColorTransferFunction.h>
#include <vtkLookupTable.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#define REFLECT 1450
#define SKIP 20

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))
typedef unsigned char inputPixelType; //want this
typedef unsigned char OutputPixelType;

typedef itk::Vector<unsigned char, 3> VectorPixelType;
typedef itk::Image<VectorPixelType, 3> ColorImageType;
typedef itk::Image<VectorPixelType, 2> Color2DImageType;

typedef itk::Image<inputPixelType,3> inputImageType; //want this
typedef itk::Image<OutputPixelType,3> OutputImageType;
typedef itk::Image<short int,3> LabelImageType;
typedef itk::Image<int,3> IntImageType;

typedef itk::Image<inputPixelType,2> Input2DImageType;
typedef itk::Image<inputPixelType,2> Output2DImageType;

typedef itk::ImageRegionConstIterator<inputImageType> ConstIteratorType;
typedef itk::ImageRegionIterator<inputImageType> IteratorType;

typedef itk::ImageRegionConstIterator<LabelImageType> ConstLabelIteratorType;
typedef itk::ImageRegionIterator<LabelImageType> LabelIteratorType;


typedef itk::ImageRegionConstIterator<Input2DImageType> Const2DIteratorType;
typedef itk::ImageRegionIterator<Input2DImageType> twoDIteratorType;

typedef itk::ImageRegionConstIterator<ColorImageType> ConstColorIteratorType;
typedef itk::ImageRegionIterator<ColorImageType> ColorIteratorType;

typedef itk::ImageLinearIteratorWithIndex< Input2DImageType > LinearIteratorType;
typedef itk::ImageSliceConstIteratorWithIndex< inputImageType > SliceIteratorType;

typedef itk::ImageLinearIteratorWithIndex< Color2DImageType > LinearColorIteratorType;
typedef itk::ImageSliceConstIteratorWithIndex< ColorImageType > SliceColorIteratorType;


typedef itk::VTKImageExport<inputImageType> ExportFilterType; //want this

struct cubeCoord
{
	unsigned short sx, sy, sz;
	unsigned short ex, ey, ez;
};

inputImageType::Pointer getEmpty(int s1, int s2, int s3)
{
	inputImageType::Pointer p = inputImageType::New();
	inputImageType::SizeType size;
	inputImageType::IndexType index;
	inputImageType::RegionType region;
	size[0] = s1;
	size[1] = s2;
	size[2] = s3;
	index.Fill(0);
	region.SetSize(size);
	region.SetIndex(index);
	p->SetRegions(region);
	p->Allocate();
	return p;
}

template <typename T>
typename T::Pointer readImage(const std::string filename)
{
	std::cout << "Reading " << filename << " ..." << std::endl;
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
		std::cerr << "Exception object caught!" << std::endl;
		std::cerr << err << std::endl;
	}
	std::cout << "Done." << std::endl;
	return reader->GetOutput();
}

template <typename ITK_Exporter, typename VTK_Importer>
void connectPipelines(ITK_Exporter exporter, VTK_Importer *importer)
//void ConnectPipelines(ITK_Exporter exporter, VTK_Importer* importer)
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

/*
vtkSmartPointer<vtkPolyData> getVTKPolyDataPrecise(LabelImageType::Pointer label)
{
	LabelIteratorType labelItr = LabelIteratorType(label, label->GetLargestPossibleRegion());
	labelItr.GoToBegin();

	//Find the maximum number of cells
	int max = 0;
	for(labelItr.GoToBegin();!labelItr.IsAtEnd();++labelItr)
	{
		max = MAX(max, labelItr.Get());
	}

	//Find all the cubes in which cells lie
	cubeCoord* cubeArray = new cubeCoord[max + 1];
	for(int i = 0; i <= max; i++)
	{
		cubeArray[i].sx = 6000;
		cubeArray[i].sy = 6000;
		cubeArray[i].sz = 6000;
		cubeArray[i].ex = 0;
		cubeArray[i].ey = 0;
		cubeArray[i].ez = 0;
	}
    typedef itk::ImageRegionConstIteratorWithIndex<LabelImageType> ConstLabelIteratorWithIndex;
	ConstLabelIteratorWithIndex labItr = ConstLabelIteratorWithIndex(label,label->GetLargestPossibleRegion());
	inputImageType::IndexType index;
	for(labItr.GoToBegin();!labItr.IsAtEnd();++labItr)
	{
		int current = labItr.Get();
		if(current != 0)
		{
			index = labItr.GetIndex();
			cubeArray[current].sx = MIN(index[0], cubeArray[current].sx);
			cubeArray[current].sy = MIN(index[1], cubeArray[current].sy);
			cubeArray[current].sx = MIN(index[2], cubeArray[current].sz);
			cubeArray[current].ex = MIN(index[0], cubeArray[current].ex);
			cubeArray[current].ey = MIN(index[1], cubeArray[current].ey);
			cubeArray[current].ez = MIN(index[2], cubeArray[current].ez);
		}
	}
	//Find the largest image size necessary
	unsigned short wx = 0, wy = 0, wz = 0;
	for(int j = 1; j <= max; j++)
	{
		wx = MAX(cubeArray[j].ex-cubeArray[j].sx+1,wx);
		wy = MAX(cubeArray[j].ey-cubeArray[j].sy+1,wy);
		wz = MAX(cubeArray[j].ez-cubeArray[j].sz+1,wz);
	}

	//Accommodate padding
	wx += 2;
	wy += 2;
	wz += 2;
	std::cout << "wx: " << wx << " wy: " << wy << " wz: " << wz << std::endl;

	vtkSmartPointer<vtkAppendPolyData> appendfilter = vtkSmartPointer<vtkAppendPolyData>::New();

	ExportFilterType::Pointer itkexporter = ExportFilterType::New();
	vtkSmartPointer<vtkImageImport> vtkimporter = vtkSmartPointer<vtkImageImport>::New();
	connectPipelines(itkexporter,(vtkImageImport *)vtkimporter);

	std::cout << "here" << std::endl;
	vtkSmartPointer<vtkMarchingCubes> contourFilter = vtkSmartPointer<vtkMarchingCubes>::New();
	contourFilter->SetInput(vtkimporter->GetOutput());
	contourFilter->SetValue(0,127);
	contourFilter->ComputeNormalsOff();
	contourFilter->ComputeScalarsOff();
	contourFilter->ComputeGradientsOff();

	vtkSmartPointer<vtkSmoothPolyDataFilter> smoothFilter = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
	smoothFilter->SetInput(contourFilter->GetOutput());
	smoothFilter->SetRelaxationFactor(0.3);
	smoothFilter->SetNumberOfIterations(20);

	vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
	transform->PostMultiply();
	transform->Identity();
	
	vtkSmartPointer<vtkTransformPolyDataFilter> transFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	transFilter->SetTransform(transform);
	transFilter->SetInput(smoothFilter->GetOutput());

	inputImageType::Pointer t = getEmpty(wx,wy,wz);
	for(int k = 1; k <= max; k++)
	{
		inputImageType::SizeType size;
		inputImageType::RegionType region;
		index.Fill(1);

		region.SetIndex(index);
		region.SetSize(size);

		LabelImageType::SizeType lsize;
		LabelImageType::IndexType lindex;
		LabelImageType::RegionType lregion;
		
		itkexporter->SetInput(t);

		t->FillBuffer(0);
		lsize[0] = cubeArray[k].ex-cubeArray[k].sx+1;
		lsize[1] = cubeArray[k].ey-cubeArray[k].sy+1;
		lsize[2] = cubeArray[k].ez-cubeArray[k].sz+1;

		lindex[0] = cubeArray[k].sx;
		lindex[1] = cubeArray[k].sy;
		lindex[2] = cubeArray[k].sz;

		lregion.SetIndex(lindex);
		lregion.SetSize(lsize);
		LabelIteratorType localItr = LabelIteratorType(label,lregion);

		size = lsize;
		region.SetSize(size);
		IteratorType itr = IteratorType(t,region);
		for(localItr.GoToBegin(),itr.GoToBegin();!localItr.IsAtEnd();++localItr, ++itr)
		{
			if(localItr.Get() == k)
			{
				itr.Set(255);
			}
		}
		t->Modified();
		vtkimporter->Modified();

		transform->Identity();
		transform->Translate(cubeArray[k].sx-1,cubeArray[k].sy-1,cubeArray[k].sz-1);
		transFilter->SetTransform(transform);
		transFilter->Update();

		//vtkSmartPointer<vtkPolyData> pol = vtkSmartPointer<vtkPolyData>::New();
		vtkPolyData *pol = vtkPolyData::New();
		pol->DeepCopy(transFilter->GetOutput());

		appendfilter->AddInput(pol);
	}

	appendfilter->Update();
	vtkSmartPointer<vtkDecimatePro> decimate = vtkSmartPointer<vtkDecimatePro>::New();
	decimate->SetInput(appendfilter->GetOutput());
	decimate->SetTargetReduction(0.75);
	std::cout << "Decimating the contours..." << std::endl;
	decimate->Update();
	std::cout << "Done." << std::endl;

	std::cout << "Smoothing the contours after decimation..." << std::endl;
	vtkSmartPointer<vtkSmoothPolyDataFilter> smoothFinal = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
	smoothFinal->SetRelaxationFactor(0.2);
	smoothFinal->SetInput(decimate->GetOutput());
	smoothFinal->SetNumberOfIterations(0);
	smoothFinal->Update();
	vtkSmartPointer<vtkPolyData> output = smoothFinal->GetOutput();

	delete [] cubeArray;
	return output;
}
*/

vtkSmartPointer<vtkPolyData> getVTKPolyDataPrecise(LabelImageType::Pointer label)
{
	LabelIteratorType liter = LabelIteratorType(label,label->GetLargestPossibleRegion());
	liter.GoToBegin();

	//find the maximum number of cells
	unsigned short max1 = 0;
	for(liter.GoToBegin();!liter.IsAtEnd();++liter)
	{
		max1 = MAX(max1,liter.Get());
	}

	//find all the cubes in which cells lie
	cubeCoord* carray = new cubeCoord[max1+1];
	for(int counter=0; counter<=max1; counter++)
	{
		carray[counter].sx=6000;carray[counter].sy=6000;carray[counter].sz=6000;
		carray[counter].ex=0;carray[counter].ey=0;carray[counter].ez=0;
	}

	typedef itk::ImageRegionConstIteratorWithIndex<LabelImageType> ConstLabelIteratorWithIndex;
	ConstLabelIteratorWithIndex cliter = ConstLabelIteratorWithIndex(label,label->GetLargestPossibleRegion());
	inputImageType::IndexType index;
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

	ExportFilterType::Pointer itkexporter = ExportFilterType::New();
			
	vtkSmartPointer<vtkImageImport> vtkimporter = vtkSmartPointer<vtkImageImport>::New();
	connectPipelines(itkexporter,(vtkImageImport *)vtkimporter);
	vtkSmartPointer<vtkMarchingCubes> contourf = vtkSmartPointer<vtkMarchingCubes>::New();


	int a = vtkimporter->GetNumberOfScalarComponents();
	std::cout << a << " components" << std::endl;

	char* scalars = vtkimporter->GetScalarArrayName();


	contourf->SetInput(vtkimporter->GetOutput());
	contourf->SetValue(0,127);
	contourf->ComputeNormalsOff();
	contourf->ComputeScalarsOff();
	contourf->ComputeGradientsOff();
	vtkSmartPointer<vtkSmoothPolyDataFilter> smoothf = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
	smoothf->SetInput(contourf->GetOutput());
	smoothf->SetRelaxationFactor(0.3);
	smoothf->SetNumberOfIterations(20);

	vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
	transform->PostMultiply();
	
	transform->Identity();
	vtkSmartPointer<vtkTransformPolyDataFilter> tf = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	tf->SetTransform(transform);
	tf->SetInput(smoothf->GetOutput());
	inputImageType::Pointer t = getEmpty(wx,wy,wz);
	for(int counter=1; counter<=max1; counter++)
	{
		//	printf("Maximum tiny image size I need is [%d %d %d]\n",wx,wy,wz);

		inputImageType::SizeType size;
		inputImageType::RegionType region;
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
		
		appendfilter->AddInput(pol);
		printf("Completed %d/%d\r",counter,max1);
	}

	appendfilter->Update();
	vtkSmartPointer<vtkDecimatePro> decimate = vtkSmartPointer<vtkDecimatePro>::New();
	decimate->SetInput(appendfilter->GetOutput());
	decimate->SetTargetReduction(0.75);
	//decimate->SetNumberOfDivisions(32,32,32);
	printf("Decimating the contours...");
	decimate->Update();
	printf("Done\n");
	printf("Smoothing the contours after decimation...");
	vtkSmartPointer<vtkSmoothPolyDataFilter> smoothfinal = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
	smoothfinal->SetRelaxationFactor(0.2);
	smoothfinal->SetInput(decimate->GetOutput());
	smoothfinal->SetNumberOfIterations(0);
	smoothfinal->Update();
	printf("Done\n");
	vtkSmartPointer<vtkPolyData> out = smoothfinal->GetOutput();
	return out;
}

void generate(const std::string filename, std::string ply, std::vector<vtkSmartPointer<vtkPolyData> > &data)
{
	vtkSmartPointer<vtkImageData> img;
	LabelImageType::Pointer itkImg = readImage<LabelImageType>(filename);
	vtkSmartPointer<vtkPolyData> polyDat = getVTKPolyDataPrecise(itkImg);

	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInput(polyDat);
	mapper->Update();
	data.push_back(polyDat);

	vtkSmartPointer<vtkPolyDataWriter> plyWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
	//double *range = new double;
	//range = polyDat->GetScalarRange();
	//std::cout << *range << std::endl;
	//plyWriter->WriteScalarData();
	plyWriter->SetInput(polyDat);

	ofstream plyoutput(ply.c_str());
	//plyWriter->SetFileName(ply.c_str());
	plyWriter->WriteToOutputStringOn();
	plyWriter->Write();
	//std::cout << plyWriter->GetOutputString() << std::endl;
	plyoutput << plyWriter->GetOutputString();
	plyoutput.close();
}

void renderData(int n, float scale[], std::vector<vtkSmartPointer<vtkPolyData> > &data, std::vector<std::string> plys)
{
	double colors[10][3];
	colors[0][0] = 1; colors[0][1] = 0; colors[0][2] = 0;
	colors[1][0] = 0; colors[1][1] = 1; colors[1][2] = 0;
	colors[2][0] = 0; colors[2][1] = 0; colors[2][2] = 1;
	colors[3][0] = 1; colors[3][1] = 1; colors[3][2] = 0;
	colors[4][0] = 1; colors[4][1] = 0; colors[4][2] = 1;
	colors[5][0] = 0; colors[5][1] = 1; colors[5][2] = 1;
	colors[6][0] = .5; colors[6][1] = 0; colors[7][1]= 0;
	colors[7][0] = 0; colors[7][1] = .5; colors[7][2] = 0;
	colors[8][0] = 0; colors[8][1] = 0; colors[8][2] = .5;
	colors[9][0] = .5; colors[9][1] = .5; colors[9][2] = 0;

	vtkRenderer *renderer = vtkRenderer::New();

	double val[6];
	for(int i = 0; i < n; i++)
	{
		int clr = i % 10;
		std::cout << clr << std::endl;
		vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
		reader->SetFileName(plys[i].c_str());
		//vtkSmartPointer<vtkPolyDataNormals> polyNorm = vtkSmartPointer<vtkPolyDataNormals>::New();
		//polyNorm->SetInput(reader->GetOutput());

		vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		//mapper->SetInput(polyNorm->GetOutput());
		mapper->SetInput(reader->GetOutput());
		mapper->ImmediateModeRenderingOff();

		mapper->Update();
		//std::cout << "here" << std::endl;
		/*
		vtkColorTransferFunction *colorFun = vtkColorTransferFunction::New();
		colorFun->AddRGBPoint(0.0,0.0,0.0,0.0);
		colorFun->AddRGBPoint(150,0.5,0.1,0.0);
		*/


		vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		actor->GetProperty()->SetColor(colors[clr][0], colors[clr][1], colors[clr][2]);
		std::cout << colors[clr][0] << colors[clr][1] << colors[clr][2];
		//actor->GetProperty()->SetColor((clrs[i])[0], (clrs[i])[1], (clrs[i])[2]);
		//std::cout << clrs[n][0] << clrs[n][1] << clrs[n][2] << std::endl;
		//actor->RotateZ(90);
		actor->SetScale(scale[0], scale[1], scale[2]);

		actor->GetBounds(val);
		std::cout << "The bounds are: ";
		for(int m = 0; m < 6; m++)
		{
			std::cout << val[m] << " ";
		}
		std::cout << std::endl;

		renderer->AddActor(actor);
	}
	renderer->SetBackground(0,0,0);

	vtkRenderWindow *renwin = vtkRenderWindow::New();
	renwin->AddRenderer(renderer);
	renwin->SetSize(500,500);
	
	renderer->ResetCamera();
	renderer->GetActiveCamera()->Zoom(1.3);
	std::cout << "Distace to focal point: " << renderer->GetActiveCamera()->GetDistance() << std::endl;
	double *pos = renderer->GetActiveCamera()->GetPosition();
	std::cout << "Camera position: " << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;

	vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
	iren->SetRenderWindow(renwin);
	vtkInteractorStyleTrackballCamera *style = vtkInteractorStyleTrackballCamera::New();
	iren->SetInteractorStyle(style);

	renwin->Render();

	iren->Start();
}

bool file_exists(std::string filename)
{
	ifstream file(filename.c_str());
	if(!file)
	{
		return 0;
	}
	return 1;
}

int main(int argc, char *argv[])
{
	if(argc != 2)
	{
		std::cerr << "Usage: Executable Input-File" << std::endl;
		return -1;
	}

	//std::cout << argv[1] << " " << std::endl;
	std::vector<vtkSmartPointer<vtkPolyData> > polyVect;
	ifstream segParam(argv[1]);
	float scale[3];
	std::vector<std::string> fileNames;
	std::vector<std::string> plyFiles;
	int counter = 0;
	std::string fileName;
	std::string colorFile;
	std::string outputFile;
	std::string fileExt;
	std::string plyFile;

	if(!segParam)
	{
		std::cerr << "Error opening segmentation parameter file" << std::endl;
		return -1;
	}
 
	segParam >> scale[0] >> scale[1] >> scale[2];

	segParam >> fileName;
	if(fileName != " " && fileName != "")
	{
		for(unsigned int i = fileName.size() - 1; i >= 0; i--)
		{
			if(fileName[i] == '.')
			{
				fileExt = fileName.substr(i + 1, (fileName.size() - i - 1));
				outputFile = fileName.substr(0, i);
				break;
			}
		}
		if(fileExt == "tif" || fileExt == "tiff")
		{ 
			//mkdir("cache");
			_mkdir("cache");
			outputFile.append("_out.tif");

			std::string output = "cache/";
			output.append(outputFile);
			std::cout << output << std::endl;
			ftk::NuclearSegmentation *segmentation = new ftk::NuclearSegmentation();
			segmentation->RunGraphColoring(fileName, output);
			//std::cout << segmentation->colorImages.size() << std::endl;

			counter = segmentation->colorImages.size();
			for(unsigned int i = 0; i < segmentation->colorImages.size(); i++)
			{
				colorFile = segmentation->colorImages[i];
				plyFile = segmentation->colorImages[i].append(".ply");
				plyFiles.push_back(plyFile);
				if(!file_exists(plyFiles[i]))
				{
					generate(colorFile, plyFiles[i], polyVect);
				}

			}
			delete segmentation;
		}
		else
		{
			std::cerr << "Invalid file type for segmentation output" << std::endl;
			return -1;
		}
	}
	else
	{
		std::cerr << "Invalid filename in parameters, load valid parameter file" << std::endl;
		return -1;
	}

	/*
	while(!segParam.eof())
	{
		std::string fileName;
		std::string plyFile;
		std::string fileExt;
		std::vector<double> color;
		double x, y, z;
		segParam >> fileName >> x >> y >> z;
		if(fileName != " " && fileName != "")
		{
			for(unsigned int i = fileName.size() - 1; i >= 0; i--)
			{
				if(fileName[i] == '.')
				{
					fileExt = fileName.substr(i+1, (fileName.size() - 1 - i));
					break;
				}
			}
			if(fileExt == "tif" || fileExt == "tiff")
			{
				color.clear();
				color.push_back(x); color.push_back(y); color.push_back(z);
				//std::cout << color[0] << color[1] << color[2] << std::endl;
				
				colors.push_back(color);
				fileNames.push_back(fileName);
				plyFile = fileName.append(".ply");
				plyFiles.push_back(plyFile);
				if(!file_exists(plyFiles[counter]))
				{
					generate(fileNames[counter], plyFiles[counter], polyVect);
				}
			}
			else
			{
				std::cerr << "Improper file type for segmentation output" << std::endl;
				return -1;
			}
			counter++;
		}
	}
	*/

	segParam.close();
	renderData(counter, scale, polyVect, plyFiles);
	return 0;
}