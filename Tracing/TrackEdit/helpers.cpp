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

#include "helpers.h"

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
	template <typename T>
int writeImage(typename T::Pointer im, const char* filename)
{
	printf("Writing %s ... ",filename);
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
	printf("Done.\n");
	return EXIT_SUCCESS;
}

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

Input2DImageType::Pointer get2DEmpty(int s1, int s2)
{
	Input2DImageType::Pointer p = Input2DImageType::New();
	Input2DImageType::SizeType size;
	Input2DImageType::IndexType index;
	Input2DImageType::RegionType region;
	size[0] = s1;size[1] = s2;
	index.Fill(0);
	region.SetSize(size);
	region.SetIndex(index);
	p->SetRegions(region);
	p->Allocate();
	return p;
}

vtkSmartPointer<vtkPolyData> getVTKContourFromITKImage(Input2DImageType::Pointer im2d)
{
	vtkSmartPointer<vtkContourFilter> cofilt = vtkSmartPointer<vtkContourFilter>::New();
	Connector2DType::Pointer conn = Connector2DType::New();
	conn->SetInput(im2d);
	cofilt->SetInput(conn->GetOutput());
	cofilt->SetValue(0,254);
	cofilt->Update();
	vtkSmartPointer<vtkPolyData> out= vtkSmartPointer<vtkPolyData>::New();
	out->DeepCopy(cofilt->GetOutput());// to disconnect pipelines;
	return out;

}
vtkSmartPointer<vtkPolyData> ShiftPolyData(vtkSmartPointer<vtkPolyData> p, double x, double y, double z)
{
	vtkSmartPointer<vtkTransform> trans = vtkSmartPointer<vtkTransform>::New();
	trans->Translate(x,y,z);
	vtkSmartPointer<vtkTransformPolyDataFilter> polytrans = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	polytrans->SetInput(p);
	polytrans->SetTransform(trans);
	polytrans->Update();
	return polytrans->GetOutput();

}
Input2DImageType::Pointer getProjection(InputImageType::Pointer im)
{
//	printf("Entered getProjection\n");
	//im->Print(std::cout);
	//printf("finished printing input\n");
	Input2DImageType::Pointer output = Input2DImageType::New();
	Input2DImageType::RegionType region;
	Input2DImageType::SizeType size;
	Input2DImageType::IndexType index;
	index[0]=0;
	index[1]=0;
	size[0] = im->GetLargestPossibleRegion().GetSize()[0];
	size[1] = im->GetLargestPossibleRegion().GetSize()[1];
	region.SetSize(size);
	region.SetIndex(index);
	output->SetRegions(region);
	//output->Print(std::cout);
//	printf("About to allocate output\n");
	output->Allocate();

//	output->Print(std::cout);
	SliceIteratorType inputIt(im,im->GetLargestPossibleRegion());
	LinearIteratorType outputIt(output,output->GetLargestPossibleRegion());

	inputIt.SetFirstDirection(0);
	inputIt.SetSecondDirection(1);
	outputIt.SetDirection(0);



	outputIt.GoToBegin();
//	printf("Begin setting default values to NonPositiveMin()\n");
	while ( ! outputIt.IsAtEnd() )
	{
		while ( ! outputIt.IsAtEndOfLine() )
		{
			outputIt.Set( itk::NumericTraits<unsigned short>::NonpositiveMin() );
			++outputIt;
		}
		outputIt.NextLine();
	}

	inputIt.GoToBegin();
	outputIt.GoToBegin();
//	printf("Begin getting the projections\n");
	while( !inputIt.IsAtEnd() )
	{
		while ( !inputIt.IsAtEndOfSlice() )
		{
			while ( !inputIt.IsAtEndOfLine() )
			{
				outputIt.Set( MAX( outputIt.Get(), inputIt.Get() ));
				++inputIt;
				++outputIt;
			}
			outputIt.NextLine();
			inputIt.NextLine();
		}
		outputIt.GoToBegin();
		inputIt.NextSlice();

	}
//	printf("Exiting getProjection\n");
	return output;
}
vtkSmartPointer<vtkPolyData> get2DBoundary(LabelImageType::Pointer label)
{
	DEBUG3("Entered get2DBoundary\n");
	LabelIteratorType liter = LabelIteratorType(label,label->GetLargestPossibleRegion());
	liter.GoToBegin();

	//find the maximum number of cells
	unsigned short max1 = 0;
	for(liter.GoToBegin();!liter.IsAtEnd();++liter)
		max1 = MAX(max1,liter.Get());

	//find all the cubes in which cells lie
	std::vector<cubecoord> carray(max1+1);
	for(int counter=0; counter<=max1; counter++)
	{
		carray[counter].sx=60000;carray[counter].sy=60000;carray[counter].sz=60000;
		carray[counter].ex=0;carray[counter].ey=0;carray[counter].ez=0;
	}

	DEBUG2("Found bboxes\n");
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
		if(carray[counter].sx>50000)
			continue;
		wx = MAX(carray[counter].ex-carray[counter].sx+1,wx);
		wy = MAX(carray[counter].ey-carray[counter].sy+1,wy);
		wz = MAX(carray[counter].ez-carray[counter].sz+1,wz);
	}
	// accommodate padding
	wx = wx+2;wy = wy +2; wz = wz+2;
	// create a tiny image of maximum size

	LabelImageType::SizeType globalsize = label->GetLargestPossibleRegion().GetSize();
	Input2DImageType::Pointer t2d = get2DEmpty(wx,wy);
	Input2DImageType::Pointer bound_im = get2DEmpty(globalsize[0],globalsize[1]);
	
	vtkSmartPointer<vtkAppendPolyData> total_contour = vtkSmartPointer<vtkAppendPolyData>::New();
	bound_im->FillBuffer(0);
	typedef itk::ImageRegionIteratorWithIndex<Input2DImageType> Iterator2DWithIndexType;

	PROGRESS("Generating Borders: [");
	for(int counter=1; counter<=max1; counter++)
	{

		if((counter%(max1/10))==0)
			PROGRESS("#");

		if(carray[counter].sx>50000)
		{
			//printf("Skipped 1 contour\n");
			continue;
		}
		Input2DImageType::SizeType size;
		Input2DImageType::IndexType index2d;
		Input2DImageType::RegionType region;
		index2d.Fill(1);

		region.SetIndex(index2d);

		LabelImageType::SizeType lsize;
		LabelImageType::IndexType lindex;
		LabelImageType::RegionType lregion;

		t2d->FillBuffer(0);
		lsize[0] = carray[counter].ex-carray[counter].sx+1;
		lsize[1] = carray[counter].ey-carray[counter].sy+1;
		lsize[2] = carray[counter].ez-carray[counter].sz+1;
		
		lindex[0] = carray[counter].sx;
		lindex[1] = carray[counter].sy;
		lindex[2] = carray[counter].sz;

		lregion.SetIndex(lindex);
		lregion.SetSize(lsize);

		ConstLabelIteratorWithIndex localiter = ConstLabelIteratorWithIndex(label,lregion);

		size[0] = lsize[0];size[1] = lsize[1];
		region.SetSize(size);
		twoDIteratorType iter = twoDIteratorType(t2d,region);


		t2d->FillBuffer(0);
		for(localiter.GoToBegin(),iter.GoToBegin();!localiter.IsAtEnd();++localiter,++iter)
		{
			if(iter.IsAtEnd())
				iter.GoToBegin();
			if(iter.Get()==0)
			{
				iter.Set(255*(localiter.Get()!=0));
			}
		}
		
		vtkSmartPointer<vtkPolyData> contour = getVTKContourFromITKImage(t2d);

		total_contour->AddInput(ShiftPolyData(contour,static_cast<double>(lindex[0])-1,static_cast<double>(lindex[1])-1,1));
		/*
		Iterator2DWithIndexType iter2d = Iterator2DWithIndexType(t2d,region);
		Input2DImageType::IndexType idx;
		idx[0] = lindex[0];idx[1] = lindex[1];
		region.SetIndex(idx);
		twoDIteratorType global2diter = twoDIteratorType(bound_im,region);
		global2diter.GoToBegin();
		for(iter2d.GoToBegin(); !iter2d.IsAtEnd();++iter2d,++global2diter)
		{
			unsigned char cur = iter2d.Get();
			if(cur==0)
				continue;
			index2d = iter2d.GetIndex();
			//check the four neighbors to see if this pixel is in boundary
			index2d[0]--;
			
			if(t2d->GetPixel(index2d)==0)
			{
				global2diter.Set(255);
				continue;
			}
			index2d[0]++;
			index2d[1]--;
			if(t2d->GetPixel(index2d)==0)
			{
				global2diter.Set(255);
				continue;
			}
			index2d[1]++;
			index2d[1]++;
			if(t2d->GetPixel(index2d)==0)
			{
				global2diter.Set(255);
				continue;
			}
			index2d[1]--;
			index2d[0]++;
			if(t2d->GetPixel(index2d)==0)
			{
				global2diter.Set(255);
				continue;
			}
			index2d[0]--;
		}
*/
	}
	PROGRESS("]\n");
	total_contour->Update();
	//total_contour->GetOutput()->Print(std::cout);
	vtkSmartPointer<vtkPolyDataNormals> poly_norm = vtkSmartPointer<vtkPolyDataNormals>::New();
	poly_norm->SetInput(total_contour->GetOutput());
	poly_norm->Update();
	DEBUG3("Exiting get2DBoundary\n");
	return total_contour->GetOutput();

}

vtkSmartPointer<vtkPolyData> getRectangle(double x1, double y1, double x2, double y2)
{
	DEBUG3("Entered getRectangle\n");
	vtkSmartPointer<vtkPoints> sqp = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> sqc = vtkSmartPointer<vtkCellArray>::New();
	double points[3];unsigned int ids[4];
	points[0] = x1; points[1] = y1; points[2] = 1.0;
	ids[0] = sqp->InsertNextPoint(points);
	points[0] = x1; points[1] = y2; points[2] = 1.0;
	ids[1] = sqp->InsertNextPoint(points);
	points[0] = x2; points[1] = y2; points[2] = 1.0;
	ids[2] = sqp->InsertNextPoint(points);
	points[0] = x2; points[1] = y1; points[2] = 1.0;
	ids[3] = sqp->InsertNextPoint(points);
	sqc->InsertNextCell(2);
	sqc->InsertCellPoint(ids[0]); sqc->InsertCellPoint(ids[1]);
	sqc->InsertNextCell(2);
	sqc->InsertCellPoint(ids[1]); sqc->InsertCellPoint(ids[2]);
	sqc->InsertNextCell(2);
	sqc->InsertCellPoint(ids[2]); sqc->InsertCellPoint(ids[3]);
	sqc->InsertNextCell(2);
	sqc->InsertCellPoint(ids[3]); sqc->InsertCellPoint(ids[0]);
	vtkSmartPointer<vtkPolyData> poly  = vtkSmartPointer<vtkPolyData>::New();
	poly->SetPoints(sqp);
	poly->SetLines(sqc);
	
	DEBUG3("Leaving getRectangle\n");
	return poly;

}

Input2DImageType::Pointer get2DBoundaryImage(LabelImageType::Pointer label)
{
	LabelIteratorType liter = LabelIteratorType(label,label->GetLargestPossibleRegion());
	liter.GoToBegin();

	//find the maximum number of cells
	unsigned short max1 = 0;
	for(liter.GoToBegin();!liter.IsAtEnd();++liter)
		max1 = MAX(max1,liter.Get());

	//find all the cubes in which cells lie
	std::vector<cubecoord> carray(max1+1);
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

	LabelImageType::SizeType globalsize = label->GetLargestPossibleRegion().GetSize();
	Input2DImageType::Pointer t2d = get2DEmpty(wx,wy);
	Input2DImageType::Pointer bound_im = get2DEmpty(globalsize[0],globalsize[1]);
	
	vtkSmartPointer<vtkAppendPolyData> total_contour = vtkSmartPointer<vtkAppendPolyData>::New();
	bound_im->FillBuffer(0);
	typedef itk::ImageRegionIteratorWithIndex<Input2DImageType> Iterator2DWithIndexType;

	for(int counter=1; counter<=max1; counter++)
	{

		Input2DImageType::SizeType size;
		Input2DImageType::IndexType index2d;
		Input2DImageType::RegionType region;
		index2d.Fill(1);

		region.SetIndex(index2d);

		LabelImageType::SizeType lsize;
		LabelImageType::IndexType lindex;
		LabelImageType::RegionType lregion;

		t2d->FillBuffer(0);
		lsize[0] = carray[counter].ex-carray[counter].sx+1;
		lsize[1] = carray[counter].ey-carray[counter].sy+1;
		lsize[2] = carray[counter].ez-carray[counter].sz+1;

		lindex[0] = carray[counter].sx;
		lindex[1] = carray[counter].sy;
		lindex[2] = carray[counter].sz;

		lregion.SetIndex(lindex);
		lregion.SetSize(lsize);

		ConstLabelIteratorWithIndex localiter = ConstLabelIteratorWithIndex(label,lregion);

		size[0] = lsize[0];size[1] = lsize[1];
		region.SetSize(size);
		twoDIteratorType iter = twoDIteratorType(t2d,region);


		t2d->FillBuffer(0);
		for(localiter.GoToBegin(),iter.GoToBegin();!localiter.IsAtEnd();++localiter,++iter)
		{
			if(iter.IsAtEnd())
				iter.GoToBegin();
			if(iter.Get()==0)
			{
				iter.Set(255*(localiter.Get()!=0));
			}
		}
		
		
		Iterator2DWithIndexType iter2d = Iterator2DWithIndexType(t2d,region);
		Input2DImageType::IndexType idx;
		idx[0] = lindex[0];idx[1] = lindex[1];
		region.SetIndex(idx);
		twoDIteratorType global2diter = twoDIteratorType(bound_im,region);
		global2diter.GoToBegin();
		for(iter2d.GoToBegin(); !iter2d.IsAtEnd();++iter2d,++global2diter)
		{
			unsigned char cur = iter2d.Get();
			if(cur==0)
				continue;
			index2d = iter2d.GetIndex();
			//check the four neighbors to see if this pixel is in boundary
			index2d[0]--;
			
			if(t2d->GetPixel(index2d)==0)
			{
				global2diter.Set(255);
				continue;
			}
			index2d[0]++;
			index2d[1]--;
			if(t2d->GetPixel(index2d)==0)
			{
				global2diter.Set(255);
				continue;
			}
			index2d[1]++;
			index2d[1]++;
			if(t2d->GetPixel(index2d)==0)
			{
				global2diter.Set(255);
				continue;
			}
			index2d[1]--;
			index2d[0]++;
			if(t2d->GetPixel(index2d)==0)
			{
				global2diter.Set(255);
				continue;
			}
			index2d[0]--;
		}

	}
	return bound_im;

}
vtkSmartPointer<vtkActor> getActorForPolyData(vtkSmartPointer<vtkPolyData> p)
{
	//p->Print(std::cout);
	vtkSmartPointer<vtkPolyDataMapper> p_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	p_mapper->SetInput(p);
	vtkSmartPointer<vtkActor> p_actor = vtkSmartPointer<vtkActor>::New();
	p_actor->SetMapper(p_mapper);
	p_actor->GetProperty()->SetEdgeColor(1, 0, 0);
	//p_actor->GetProperty()->SetLineWidth(2);
	//p_actor->GetProperty()->SetOpacity(1);
	
	//p_actor->GetProperty()->SetEdgeColor(255, 1, 0);
	return p_actor;
}
//Color2DImageType::Pointer getColorBoundaryImage(LabelImageType::Pointer labelled, InputImageType::Pointer im, int channel)
//{
//	double multiplier;
//	if(channel==1)
//		multiplier = 6;
//	else
//		multiplier = 1.5;
//	unsigned char colorarray[][3]={255,0,0,
//		0,154,25,
//		207,141,0,
//		255,0,0};
//	VectorPixelType colorcodes[4];
//	VectorPixelType black;black[0]=0;black[1]=0;black[2]=0;
//
//	for(int counter=0; counter<4; counter++)
//		colorcodes[counter]=colorarray[counter];
//
//	Input2DImageType::Pointer boundary = get2DBoundary(labelled);
//	Color2DImageType::Pointer cimage = Color2DImageType::New();
//	Color2DImageType::RegionType colregion;
//	Color2DImageType::IndexType colindex;colindex.Fill(0);
//	Color2DImageType::SizeType colsize; 
//	colsize[0]=boundary->GetLargestPossibleRegion().GetSize()[0];
//	colsize[1]=boundary->GetLargestPossibleRegion().GetSize()[1];
//	colregion.SetSize(colsize);
//	colregion.SetIndex(colindex);
//	cimage->SetRegions(colregion);
//	cimage->Allocate();
//
//	Input2DImageType::Pointer im2d = getProjection(im);
//	VectorPixelType pix;
//	typedef itk::ImageRegionIterator<Color2DImageType> Color2DIteratorType;
//	Color2DIteratorType iter = Color2DIteratorType(cimage,cimage->GetLargestPossibleRegion());
//	twoDIteratorType biter = twoDIteratorType(boundary,boundary->GetLargestPossibleRegion());
//	twoDIteratorType inputiter = twoDIteratorType(im2d,im2d->GetLargestPossibleRegion());
//	//generate a mask image 
//	InputImageType::Pointer mask=InputImageType::New();
//	mask->SetRegions(labelled->GetLargestPossibleRegion());
//	mask->Allocate();
//	LabelIteratorType liter = LabelIteratorType(labelled,labelled->GetLargestPossibleRegion());
//	IteratorType maskiter = IteratorType(mask,mask->GetLargestPossibleRegion());
//	for(liter.GoToBegin(),maskiter.GoToBegin();!maskiter.IsAtEnd();++maskiter,++liter)
//	{
//		maskiter.Set(255*(liter.Get()!=0));
//	}
//	Input2DImageType::Pointer mask2d = getProjection(mask);
//	twoDIteratorType mask2diter = twoDIteratorType(mask2d,mask2d->GetLargestPossibleRegion());
//	//initialize all iterators
//	iter.GoToBegin();biter.GoToBegin();inputiter.GoToBegin();mask2diter.GoToBegin();
//	for(;!iter.IsAtEnd();++iter,++biter,++inputiter,++liter,++mask2diter)
//	{
//		if(biter.Get()!=0)
//		{
//			pix = colorcodes[channel-1];
//		}
//		else if(mask2diter.Get()!=0)
//		{
//			pix[0] = MIN(inputiter.Get()*multiplier,255);
//			pix[1] = pix[0];
//			pix[2] = pix[0];
//		}
//		else
//		{
//			pix = black;
//		}
//		iter.Set(pix);
//	}
//
//	return cimage;
//}
bool mergeLabels1(LabelImageType::Pointer im, int merge_to,int n1)
{
	LabelIteratorType iter = LabelIteratorType(im,im->GetLargestPossibleRegion());
	unsigned short val;
	bool flag = false;
	for(iter.GoToBegin();!iter.IsAtEnd();++iter)
	{
		val = iter.Get();
		if(val == n1)
		{
			flag = true;
			iter.Set(merge_to);
		}
	}
	return flag;
}
bool mergeLabels2(LabelImageType::Pointer im, int merge_to,int n1, int n2)
{
	LabelIteratorType iter = LabelIteratorType(im,im->GetLargestPossibleRegion());
	unsigned short val;
	bool flag = false;
	for(iter.GoToBegin();!iter.IsAtEnd();++iter)
	{
		val = iter.Get();
		if(val == n1 || val == n2)
		{
			flag = true;
			iter.Set(merge_to);
		}
	}
	return flag;
}
bool mergeLabels3(LabelImageType::Pointer im, int merge_to,int n1, int n2,int n3)
{
	LabelIteratorType iter = LabelIteratorType(im,im->GetLargestPossibleRegion());
	unsigned short val;
	bool flag = false;
	for(iter.GoToBegin();!iter.IsAtEnd();++iter)
	{
		val = iter.Get();
		if(val == n1 || val == n2 || val == n3)
		{
			flag = true;
			iter.Set(merge_to);
		}
	}
	return flag;
}
bool mergeLabels4(LabelImageType::Pointer im, int merge_to,int n1, int n2, int n3, int n4)
{
	LabelIteratorType iter = LabelIteratorType(im,im->GetLargestPossibleRegion());
	unsigned short val;
	bool flag = false;
	for(iter.GoToBegin();!iter.IsAtEnd();++iter)
	{
		val = iter.Get();
		if(val == n1 || val == n2 || val == n3 || val == n4)
		{
			flag = true;
			iter.Set(merge_to);
		}
	}
	return flag;
}
bool mergeLabels5(LabelImageType::Pointer im, int merge_to,int n1, int n2, int n3, int n4, int n5)
{
	LabelIteratorType iter = LabelIteratorType(im,im->GetLargestPossibleRegion());
	unsigned short val;
	bool flag = false;
	for(iter.GoToBegin();!iter.IsAtEnd();++iter)
	{
		val = iter.Get();
		if(val == n1 || val == n2 || val == n3 || val == n4 || val == n5)
		{
			flag = true;
			iter.Set(merge_to);
		}
	}
	return flag;
}

bool deleteLabels1(LabelImageType::Pointer im, int n1)
{
	LabelIteratorType iter = LabelIteratorType(im,im->GetLargestPossibleRegion());
	unsigned short val;
	bool flag = false;
	for(iter.GoToBegin();!iter.IsAtEnd();++iter)
	{
		val = iter.Get();
		if(val == n1)
		{
			flag = true;
			iter.Set(0);
		}
	}
	return flag;
}

bool deleteLabels2(LabelImageType::Pointer im, int n1, int n2)
{
	LabelIteratorType iter = LabelIteratorType(im,im->GetLargestPossibleRegion());
	unsigned short val;
	bool flag = false;
	for(iter.GoToBegin();!iter.IsAtEnd();++iter)
	{
		val = iter.Get();
		if(val == n1 || val == n2)
		{
			flag = true;
			iter.Set(0);
		}
	}
	return flag;
}

bool deleteLabels3(LabelImageType::Pointer im, int n1, int n2, int n3)
{
	LabelIteratorType iter = LabelIteratorType(im,im->GetLargestPossibleRegion());
	unsigned short val;
	bool flag = false;
	for(iter.GoToBegin();!iter.IsAtEnd();++iter)
	{
		val = iter.Get();
		if(val == n1 || val == n2 || val == n3)
		{
			flag = true;
			iter.Set(0);
		}
	}
	return flag;
}

bool deleteLabels4(LabelImageType::Pointer im, int n1, int n2, int n3, int n4)
{
	LabelIteratorType iter = LabelIteratorType(im,im->GetLargestPossibleRegion());
	unsigned short val;
	bool flag = false;
	for(iter.GoToBegin();!iter.IsAtEnd();++iter)
	{
		val = iter.Get();
		if(val == n1 || val == n2 || val == n3 || val == n4)
		{
			flag = true;
			iter.Set(0);
		}
	}
	return flag;
}

bool deleteLabels5(LabelImageType::Pointer im, int n1, int n2, int n3, int n4, int n5)
{
	LabelIteratorType iter = LabelIteratorType(im,im->GetLargestPossibleRegion());
	unsigned short val;
	bool flag = false;
	for(iter.GoToBegin();!iter.IsAtEnd();++iter)
	{
		val = iter.Get();
		if(val == n1 || val == n2 || val == n3 || val == n4 || val == n5)
		{
			flag = true;
			iter.Set(0);
		}
	}
	return flag;
}

vtkSmartPointer<vtkTextActor> getActorForFeature(ftk::LabelImageFeatures f)
{
	vtkSmartPointer<vtkTextActor> tact = vtkSmartPointer<vtkTextActor>::New();
	char buff[500];
	sprintf(buff,"%d",f.num);
	tact->SetInput(buff);
	tact->GetPositionCoordinate()->SetCoordinateSystemToWorld();
	tact->GetPositionCoordinate()->SetValue(f.centroid[0],f.centroid[1],-0.5);
	tact->GetPosition2Coordinate()->SetCoordinateSystemToWorld();
	tact->GetPosition2Coordinate()->SetValue(f.centroid[0]+10,f.centroid[1]+10,-0.5);
	tact->GetTextProperty()->SetBold(1);
	//tact->GetTextProperty()->SetShadow(1);

	if(f.tag == 0)
	{
		tact->GetTextProperty()->SetColor(0,1,0.2);
	}
	else
	{
		tact->GetTextProperty()->SetColor(1,0,0);
	}
	tact->GetTextProperty()->SetJustificationToCentered();
	tact->GetTextProperty()->SetFontSize(16);
	return tact;
}


void getFeatureVectorsFarsight(LabelImageType::Pointer im, InputImageType::Pointer in_image, std::vector<ftk::LabelImageFeatures> & feature_vector, int time, int tag)
{
	printf("Started feature calculation\n");
	if(im->GetLargestPossibleRegion().GetSize()[2]==1)
	{
		//convert it to a 2D Image
		LabelImageType::SizeType linsize = im->GetLargestPossibleRegion().GetSize();
		Label2DImageType::Pointer l2d = Label2DImageType::New();
		Label2DImageType::SizeType l2dsize; l2dsize[0] = linsize[0]; l2dsize[1] = linsize[1];
		Label2DImageType::IndexType l2dindex; l2dindex.Fill(0);
		Label2DImageType::RegionType l2dregion; l2dregion.SetSize(l2dsize); l2dregion.SetIndex(l2dindex);
		l2d->SetRegions(l2dregion);
		l2d->Allocate();
		memcpy(im->GetBufferPointer(),l2d->GetBufferPointer(),sizeof(LabelImageType::PixelType)*l2dsize[0]*l2dsize[1]);

		Input2DImageType::Pointer i2d = Input2DImageType::New();
		i2d->SetRegions(l2dregion);
		i2d->Allocate();
		memcpy(in_image->GetBufferPointer(),i2d->GetBufferPointer(),sizeof(InputImageType::PixelType)*l2dsize[0]*l2dsize[1]);
	
		typedef ftk::LabelImageToFeatures<Input2DImageType::PixelType, Label2DImageType::PixelType, 2> FeatureCalculator2DType;
		FeatureCalculator2DType::Pointer fc2d = FeatureCalculator2DType::New();
		fc2d->SetImageInputs(i2d,l2d);
		fc2d->SetLevel(1);
		fc2d->Update();
		std::vector<Label2DImageType::PixelType> labels = fc2d->GetLabels();
		for(unsigned int counter=0; counter<labels.size();counter++)
		{
			if(labels[counter]==0)
				continue;
			 feature_vector.push_back(fc2d->GetFeatures(labels[counter]));
			 feature_vector.back().num=labels[counter];
			 feature_vector.back().tag = tag;
			 feature_vector.back().time = time;
			 if(feature_vector.back().spacing.size()==0)
			 {
				feature_vector.back().spacing.push_back(1.0);
				feature_vector.back().spacing.push_back(1.0);
			 }
		}
		
	}
	else
	{
		typedef ftk::LabelImageToFeatures<InputImageType::PixelType,LabelImageType::PixelType,3> FeatureCalculatorType;
		FeatureCalculatorType::Pointer fc = FeatureCalculatorType::New();
		fc->SetImageInputs(in_image,im);
		fc->SetLevel(1);
		fc->Update();
		std::vector<LabelImageType::PixelType> labels = fc->GetLabels();
		for(unsigned int counter=0; counter< labels.size(); counter++)
		{
			if(labels[counter]==0)
				continue;
			feature_vector.push_back(fc->GetFeatures(labels[counter]));
			feature_vector.back().num = labels[counter];
			feature_vector.back().tag = tag;
			feature_vector.back().time = time;
			if(feature_vector.back().spacing.size()==0)
			{
				feature_vector.back().spacing.push_back(1.0);
				feature_vector.back().spacing.push_back(1.0);
				feature_vector.back().spacing.push_back(5.0);
			}
		}
	}
	printf("Ended feature calculation\n");
}

void getFeatureVectorsNew(LabelImageType::Pointer im, InputImageType::Pointer in_image, std::vector<Feature> &feature_vector, int time, int tag)
{
	unsigned int num_objects;
	unsigned int base = feature_vector.size();

	{
		typedef itk::LabelGeometryImageFilter<LabelImageType, InputImageType> LabelGeometryFilterType;
		typedef itk::LabelStatisticsImageFilter< InputImageType, LabelImageType> LabelStatisticsFilterType;

		LabelGeometryFilterType::Pointer lgfilter = LabelGeometryFilterType::New();
		LabelStatisticsFilterType::Pointer lsfilter = LabelStatisticsFilterType::New();

		lgfilter->SetInput(im);
		lgfilter->SetIntensityInput(in_image);
		lgfilter->Update();

		num_objects = lgfilter->GetNumberOfLabels();

		lsfilter->SetInput(in_image);
		lsfilter->SetLabelInput(im);
		lsfilter->Update();

		Feature ftemp;

		std::vector <short> labels = lgfilter->GetLabels();

		for(unsigned int counter=0; counter < num_objects; counter++)
		{
			unsigned short num = labels[counter];

			if(num==0)
				continue;
			ftemp.tag = tag;
			ftemp.t = time;
			ftemp.num = labels[counter];
			ftemp.volume = lgfilter->GetVolume(num);

			LabelGeometryFilterType::LabelPointType centroid = lgfilter->GetCentroid(num);
			ftemp.x = centroid[0];
			ftemp.y = centroid[1];
			ftemp.z = centroid[2];

			LabelGeometryFilterType::BoundingBoxType bbox = lgfilter->GetBoundingBox(num);
			ftemp.bbox.sx = bbox[0];
			ftemp.bbox.ex = bbox[1];
			ftemp.bbox.sy = bbox[2];
			ftemp.bbox.ey = bbox[3];
			ftemp.bbox.sz = bbox[4];
			ftemp.bbox.ez = bbox[5];

			LabelGeometryFilterType::AxesLengthType axes_length =  lgfilter->GetAxesLength(num);
			ftemp.xaxis_length = axes_length[0];
			ftemp.yaxis_length = axes_length[1];
			ftemp.zaxis_length = axes_length[2];

			ftemp.avg_intensity = lsfilter->GetMean(num);
			ftemp.min_intensity = lsfilter->GetMinimum(num);
			ftemp.max_intensity = lsfilter->GetMaximum(num);
			ftemp.sd_intensity = lsfilter->GetSigma(num);

			feature_vector.push_back(ftemp);
		}
	}

	// fill up all the texture features now
	
	{

		typedef itk::Statistics::ScalarImageTextureCalculator<LabelImageType> TextureCalculatorType;
		LabelImageType::Pointer new_im = LabelImageType::New();
		new_im->SetRegions(im->GetLargestPossibleRegion());
		new_im->Allocate();

		LabelIteratorType liter(new_im,new_im->GetLargestPossibleRegion());
		IteratorType iter(in_image,in_image->GetLargestPossibleRegion());
		for(iter.GoToBegin(), liter.GoToBegin(); !iter.IsAtEnd();++iter,++liter)
		{
			liter.Set(static_cast<unsigned short>(iter.Get()));
		}

		TextureCalculatorType::Pointer tcalc = TextureCalculatorType::New();
		tcalc->SetInput(new_im);
		tcalc->SetImageMask(im);
		tcalc->SetPixelValueMinMax(0,255);// we want to calculate it only over the range of unsigned char
		tcalc->SetFastCalculations(1);

		for(int counter=0; counter< int(num_objects)-1; counter++)
		{
			LabelImageType::RegionType region;
			LabelImageType::SizeType size;
			LabelImageType::IndexType index;

			index[0] = feature_vector[counter+base].bbox.sx;
			index[1] = feature_vector[counter+base].bbox.sy;
			index[2] = feature_vector[counter+base].bbox.sz;

			size[0] = feature_vector[counter+base].bbox.ex-feature_vector[counter+base].bbox.sx+1;
			size[1] = feature_vector[counter+base].bbox.ey-feature_vector[counter+base].bbox.sy+1;
			size[2] = feature_vector[counter+base].bbox.ez-feature_vector[counter+base].bbox.sz+1;

			region.SetSize(size);
			region.SetIndex(index);
			im->SetRequestedRegion(region);
			new_im->SetRequestedRegion(region);

			tcalc->SetInsidePixelValue(feature_vector[counter+base].num);
			tcalc->Compute();

			TextureCalculatorType::FeatureValueVector* vec = tcalc->GetFeatureMeans();
			feature_vector[counter+base].texture_energy              = vec->ElementAt(0);
			feature_vector[counter+base].texture_entropy            = vec->ElementAt(1);
			feature_vector[counter+base].texture_inv_diff_moment    = vec->ElementAt(2);
			feature_vector[counter+base].texture_inertia            = vec->ElementAt(3);
			feature_vector[counter+base].texture_cluster_shade      = vec->ElementAt(4);
			feature_vector[counter+base].texture_cluster_prominence = vec->ElementAt(5);
		}
		
	}
	return;
}

void getFeatureVectors(LabelImageType::Pointer im,InputImageType::Pointer in_image,std::vector<Feature> &feature_vector,int time,int tag)
{
	
	LabelIteratorType ltempiter(im,im->GetLargestPossibleRegion());

	int num_objects = 0;
	for(ltempiter.GoToBegin();!ltempiter.IsAtEnd();++ltempiter)
	{
		num_objects=MAX(num_objects,ltempiter.Get());
	}
	if(num_objects==0)
	{
		printf("I got num objects = %d tag = %d time = %d\n",num_objects,tag,time);
	}

	Feature ftemp;
	ftemp.x=ftemp.y=ftemp.z=0,ftemp.t=time;
	ftemp.bbox.sx=1e4;ftemp.bbox.sy=1e4;ftemp.bbox.sz=1e4;
	ftemp.bbox.ex=0;ftemp.bbox.ey=0;ftemp.bbox.ez=0;
	ftemp.volume=0,ftemp.tag=tag;ftemp.avg_intensity =0;
//	printf("DEBUG::Number of components : %d\n",num_objects);
	int base = feature_vector.size();
	for(int counter=0;counter<num_objects; counter++)
	{
		feature_vector.push_back(ftemp);
	}

	typedef itk::ImageRegionIteratorWithIndex<LabelImageType> LabelIteratorType;

	LabelIteratorType liter(im,im->GetLargestPossibleRegion());
	IteratorType iter = IteratorType(in_image,in_image->GetLargestPossibleRegion());
	LabelImageType::IndexType index;
	int num;
	for(liter.GoToBegin(),iter.GoToBegin();!liter.IsAtEnd();++liter,++iter)
	{
		if((num=liter.Get())!=0)
		{
			index = liter.GetIndex();
			num--;
			num+=base;
			feature_vector[num].num = num + 1 - base;
			feature_vector[num].volume++;
			feature_vector[num].x+=index[0];
			feature_vector[num].y+=index[1];
			feature_vector[num].z+=index[2];
			feature_vector[num].bbox.sx = MIN(feature_vector[num].bbox.sx,index[0]);
			feature_vector[num].bbox.sy = MIN(feature_vector[num].bbox.sy,index[1]);
			feature_vector[num].bbox.sz = MIN(feature_vector[num].bbox.sz,index[2]);
			feature_vector[num].bbox.ex = MAX(feature_vector[num].bbox.ex,index[0]);
			feature_vector[num].bbox.ey = MAX(feature_vector[num].bbox.ey,index[1]);
			feature_vector[num].bbox.ez = MAX(feature_vector[num].bbox.ez,index[2]);
			feature_vector[num].avg_intensity+=iter.Get();
		}
	}
	for(int counter=0;counter<num_objects; counter++)
	{
		if(feature_vector[base+counter].volume<1)
			continue;
		int volume=feature_vector[base+counter].volume;
		feature_vector[base+counter].x/=volume;
		feature_vector[base+counter].y/=volume;
		feature_vector[base+counter].z/=volume;
		feature_vector[base+counter].avg_intensity/=volume;
	}
	std::vector<Feature>::iterator fiter = feature_vector.begin();
	while(fiter!=feature_vector.end())
	{
		if(fiter->volume<1)
		{
			fiter = feature_vector.erase(fiter);
			continue;
		}
		++fiter;
	}
	return;
}


Input2DImageType::Pointer get2DMaskedImage(InputImageType::Pointer im, LabelImageType::Pointer l)
{
	Input2DImageType::Pointer output = Input2DImageType::New();
	Input2DImageType::SizeType osize;
	osize[0] = im->GetLargestPossibleRegion().GetSize()[0];
	osize[1] = im->GetLargestPossibleRegion().GetSize()[1];
	Input2DImageType::RegionType oregion;
	Input2DImageType::IndexType oindex;
	oindex.Fill(0);
	oregion.SetSize(osize);
	oregion.SetIndex(oindex);
	output->SetRegions(oregion);
	output->Allocate();
	twoDIteratorType iter(output,output->GetLargestPossibleRegion());
	LabelIteratorType liter(l,l->GetLargestPossibleRegion());
	IteratorType imiter(im,im->GetLargestPossibleRegion());

	iter.GoToBegin();
	liter.GoToBegin();
	imiter.GoToBegin();
	output->FillBuffer(0);
	for(;!liter.IsAtEnd();++liter,++imiter)
	{
		iter.Set(MAX(iter.Get(),imiter.Get()));//*(liter.Get()>0)));
		++iter;
		if(iter.IsAtEnd())
		{
			iter.GoToBegin();
		}
	}
	return output;
}

vtkSmartPointer<vtkVolume> getOneVTKVolume(vtkSmartPointer<vtkImageData> vtkim, float colors[3])
{

	vtkSmartPointer<vtkPiecewiseFunction> opacityTransferFunction = vtkSmartPointer<vtkPiecewiseFunction>::New();
	opacityTransferFunction->AddPoint(2,0.0);
	opacityTransferFunction->AddPoint(50,0.8);

	vtkSmartPointer<vtkVolumeProperty> volumeProperty = vtkSmartPointer<vtkVolumeProperty>::New();

	
	vtkSmartPointer<vtkColorTransferFunction> colorTransferFunction = vtkSmartPointer<vtkColorTransferFunction>::New();
	colorTransferFunction->AddRGBPoint(0.0,0.0,0.0,0.0);
	colorTransferFunction->AddRGBPoint(50.0,colors[0],colors[1],colors[2]);
	volumeProperty->SetColor(colorTransferFunction);
	volumeProperty->SetScalarOpacity(opacityTransferFunction);
	volumeProperty->ShadeOff();
	volumeProperty->SetInterpolationTypeToLinear();
	volumeProperty->DisableGradientOpacityOn();
	


#ifndef WIN32
	vtkSmartPointer<vtkOpenGLVolumeTextureMapper2D> vMapper = vtkSmartPointer<vtkOpenGLVolumeTextureMapper2D>::New();
	vMapper->SetMaximumNumberOfPlanes(50);
#else
	/*vtkSmartPointer<vtkVolumeRayCastMapper> vMapper = vtkSmartPointer<vtkVolumeRayCastMapper>::New();
	vtkSmartPointer<vtkVolumeRayCastCompositeFunction> volume_ray_cast = vtkSmartPointer<vtkVolumeRayCastCompositeFunction>::New();
	vMapper->SetVolumeRayCastFunction(volume_ray_cast);
	//vMapper->SetBlendModeToMaximumIntensity();
	//vMapper->SetMaximumImageSampleDistance(20);
	//vMapper->SetMinimumImageSampleDistance(1);*/
    
	vtkSmartPointer<vtkOpenGLVolumeTextureMapper3D> vMapper = vtkSmartPointer<vtkOpenGLVolumeTextureMapper3D>::New();
	vMapper->SetPreferredMethodToNVidia();
	vMapper->SetSampleDistance(1);
	
#endif
	vMapper->SetInput(vtkim);
	vMapper->Update();
	vtkSmartPointer<vtkVolume> volume = vtkSmartPointer<vtkVolume>::New();
	volume->SetMapper(vMapper);
	volume->SetProperty(volumeProperty);
	volume->Update();
	volume->Print(std::cout);

	printf("______________________________________\n");

	return volume;
}
