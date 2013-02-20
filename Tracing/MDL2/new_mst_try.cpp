
#include <stdio.h>
#include <time.h>
#include <vector>
#include <algorithm>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageSliceConstIteratorWithIndex.h>
#include <itkImageLinearIteratorWithIndex.h>
#include <itkConstNeighborhoodIterator.h>
#include <itkMedianImageFilter.h>
#include <itkMeanImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkBinaryMedianImageFilter.h>
#include <itkOtsuThresholdImageFilter.h>
#include <itkShiftScaleImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkBinaryErodeImageFilter.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkDanielssonDistanceMapImageFilter.h>
#include <itkVector.h>
#include <itkVTKImageExport.h>
#include <itkVTKImageImport.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_hungarian_algorithm.h>
#include <itkFastMarchingImageFilter.h>

#include "itkHessian3DToVesselnessMeasureImageFilter.h"
#include "itkMultiScaleHessianBasedMeasureImageFilter.h"
#include "itkSymmetricSecondRankTensor.h"
#include "itkGrayscaleErodeImageFilter.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "mdlMST.h"
#include "mdlUtils.h"

#define MAX(a,b) (((a) > (b))?(a):(b))
#define MIN(a,b) (((a) < (b))?(a):(b))

//typedefs

typedef unsigned char InputPixelType;
typedef unsigned char OutputPixelType;

typedef itk::Vector<unsigned char, 3> VectorPixelType;
typedef itk::Image<VectorPixelType, 3> ColorImageType;
typedef itk::Image<VectorPixelType, 2> Color2DImageType;

typedef itk::Image<InputPixelType,3> InputImageType;
typedef itk::Image<float,3> FloatImageType;
typedef itk::Image<OutputPixelType,3> OutputImageType;
typedef itk::Image<short int,3> LabelImageType;

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

typedef itk::ConstNeighborhoodIterator<InputImageType> NeighborhoodIteratorType;

typedef itk::MedianImageFilter<InputImageType,InputImageType> MedianFilterType;
typedef itk::BinaryMedianImageFilter<InputImageType,InputImageType> BinaryMedianFilterType;

typedef itk::Image<bool,3> BoolImageType;
typedef itk::BinaryThresholdImageFilter<InputImageType,OutputImageType> ThresholdFilterType;
typedef itk::OtsuThresholdImageFilter<InputImageType,OutputImageType> OtsuThresholdFilterType;

typedef itk::Image<short int,3> DistanceImageType;
typedef itk::DanielssonDistanceMapImageFilter<InputImageType,DistanceImageType> DistanceMapFilterType;
typedef DistanceMapFilterType::VectorImageType OffsetImageType;

typedef itk::ConnectedComponentImageFilter<InputImageType,LabelImageType> ConnectedFilterType;
typedef itk::RelabelComponentImageFilter<LabelImageType,LabelImageType> RelabelFilterType;

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
		std::cout << "ExceptionObject caught!" <<std::endl;
		std::cout << err << std::endl;
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
		std::cout << "ExceptionObject caught!" <<std::endl;
		std::cout << err << std::endl;
		return EXIT_FAILURE;
	}
	printf("Done.\n");
	return EXIT_SUCCESS;
}

void get_tiles(InputImageType::RegionType inputr, int tilesizex, int tilesizey, int tilesizez, int borderx, int bordery, int borderz, std::vector<InputImageType::RegionType> &in1, std::vector<InputImageType::RegionType> &in2, std::vector<InputImageType::RegionType> &out1, std::vector<InputImageType::RegionType> &out2)
{
	
	int xsize = inputr.GetSize()[0];
	int ysize = inputr.GetSize()[1];
	int zsize = inputr.GetSize()[2];

	int kx = 0;int ky = 0;int kz = 0;

	kx = xsize /(tilesizex-borderx);
	ky = ysize /(tilesizey-bordery);
	kz = zsize /(tilesizez-borderz);

	int remx = xsize % (tilesizex-borderx);
	int remy = ysize % (tilesizey-bordery);
	int remz = zsize % (tilesizez-borderz);

	if ( remx > 0 )
		kx ++;
	if ( remy > 0 )
		ky ++;
	if ( remz > 0 )
		kz ++;

	for(int xco = 0; xco < kx; xco++)
	{
		for(int yco = 0; yco < ky; yco++)
		{
			for(int zco = 0; zco < kz; zco++)
			{
				InputImageType::SizeType imsize = inputr.GetSize();
				InputImageType::IndexType index;
				InputImageType::SizeType size;


				index.Fill(0);
				size[0] = MIN((xco)*(tilesizex-borderx)+tilesizex-1,imsize[0]-1) -  xco * (tilesizex-borderx) +1;
				size[1] = MIN((yco)*(tilesizey-bordery)+tilesizey-1,imsize[1]-1) -  yco * (tilesizey-bordery) +1;
				size[2] = MIN((zco)*(tilesizez-borderz)+tilesizez-1,imsize[2]-1) -  zco * (tilesizez-borderz) +1;

				InputImageType::RegionType region;
				region.SetIndex(index);
				region.SetSize(size);
				in2.push_back(region);

				InputImageType::RegionType region1;
				index[0] = xco *(tilesizex-borderx);
				index[1] = yco *(tilesizey-bordery);
				index[2] = zco *(tilesizez-borderz);
				region1.SetIndex(index);
				region1.SetSize(size);
				in1.push_back(region1);


				if(xco != 0)
				{
					size[0] = size[0] - borderx/2;
					index[0] = borderx/2;
				}
				if(xco != kx-1)
				{
					size[0] = size[0] - borderx/2;
				}

				if(yco != 0)
				{
					size[1] = size[1] - bordery/2;
					index[1] = bordery/2;
				}
				if(yco != ky-1)
				{
					size[1] = size[1] - bordery/2;
				}

				if(zco != 0)
				{
					size[2] = size[2] - borderz/2;
					index[2] = borderz/2;
				}
				if(zco != kz-1)
				{
					size[2] = size[2] - borderz/2;
				}
				region.SetIndex(index);
				region.SetSize(size);

				out2.push_back(region);

				if(xco!=0)
				{
					index[0] = xco *(tilesizex-borderx)+borderx/2;
				}
				if(yco!=0)
				{
					index[1] = yco *(tilesizey-bordery)+bordery/2;
				}
				if(zco!=0)
				{
					index[2] = zco *(tilesizez-borderz)+borderz/2;
				}
				region1.SetIndex(index);
				region1.SetSize(size);
				out1.push_back(region1);
			}
		}
	}
}

int transfer_function1(int val)
{
	float factor = 1.0/2/40/40;
	return 255*exp(float(-(val-255)*(val-255)*factor));
}

int main()
{
	InputImageType::Pointer im = readImage<InputImageType>("C:/Users/arun/Research/Farsight/exe/bin/CF_1_inverted_bg_sub.tif");
	FILE *fp = fopen("C:/Users/arun/Research/Farsight/exe/bin/seeds.txt","r");

	//IteratorType initer(im,im->GetLargestPossibleRegion());
	//initer.GoToBegin();
	//
	//for(;!initer.IsAtEnd(); ++initer)
	//{
	//	initer.Set(transfer_function1(initer.Get()));
	//}
	//
	//writeImage<InputImageType>(im,"C:/Users/arun/Research/Farsight/exe/bin/hp2_cropped2_filtered.tif");
	//
	//return 0;
	typedef itk::SymmetricSecondRankTensor<double,3> HessianType;
	typedef itk::Hessian3DToVesselnessMeasureImageFilter<float> MeasureType;
	typedef itk::Image<HessianType,3> HessianImageType;
	typedef itk::MultiScaleHessianBasedMeasureImageFilter< InputImageType, HessianImageType, FloatImageType> VesselnessFilterType;
	
	std::vector<InputImageType::RegionType> in1,in2,out1,out2;

	get_tiles(im->GetLargestPossibleRegion().GetSize(),1500,1500,1500,100,100,10,in1,in2,out1,out2);
	
	InputImageType::Pointer om;
/*
	om = InputImageType::New();
	om->SetRegions(im->GetLargestPossibleRegion());
	om->Allocate();
	for(int counter = 0; counter < in1.size(); counter++)
	{
		InputImageType::Pointer imtile = InputImageType::New();//
		imtile->SetRegions(in2[counter]);
		imtile->Allocate();

		in1[counter].Print(std::cout);
		in2[counter].Print(std::cout);
		IteratorType iter1(im,in1[counter]);
		IteratorType iter2(imtile,in2[counter]);
		for(iter1.GoToBegin(),iter2.GoToBegin();!iter1.IsAtEnd(); ++iter1,++iter2)
		{
			iter2.Set(iter1.Get());
		}

		VesselnessFilterType::Pointer vfilt = VesselnessFilterType::New();
		MeasureType::Superclass::Pointer measure = MeasureType::New();
		
		vfilt->SetInput(imtile);
		vfilt->SetHessianToMeasureFilter((VesselnessFilterType::HessianToMeasureFilterType *)measure);
		vfilt->SetSigmaMinimum(3.0);
		vfilt->SetSigmaMaximum(5.0);
		vfilt->SetNumberOfSigmaSteps(3);
		vfilt->SetSigmaStepMethod(VesselnessFilterType::EquispacedSigmaSteps);
		vfilt->Update();
		FloatImageType::Pointer omtile = vfilt->GetOutput();
		
		typedef itk::ImageRegionIterator<FloatImageType> FloatIteratorType;
		FloatIteratorType iter3;
		iter1 = IteratorType(om,out1[counter]);
		iter3 = FloatIteratorType(omtile,out2[counter]);
		for(iter1.GoToBegin(),iter3.GoToBegin();!iter1.IsAtEnd();++iter1,++iter3)
		{
			iter1.Set(iter3.Get());
		}
	}
	writeImage<InputImageType>(om,"C:/Users/arun/Research/Farsight/exe/bin/vesselnesstest.tif");
*/
	om = readImage<InputImageType>("C:/Users/arun/Research/Farsight/exe/bin/vesselnesstest.tif");
	typedef itk::BinaryBallStructuringElement<InputImageType::PixelType,3> StructElementType;
	typedef itk::GrayscaleDilateImageFilter<InputImageType,InputImageType,StructElementType> FilterType1;
	FilterType1::Pointer minfilt = FilterType1::New();
	minfilt->SetInput(om);
	FilterType1::RadiusType radius;
	radius[0] = 1;
	radius[1] = 1;
	radius[2] = 1;
	StructElementType strel;
	strel.SetRadius(radius);
	minfilt->SetKernel(strel);
	minfilt->Update();

	InputImageType::Pointer seed_out  = InputImageType::New();
	seed_out->SetRegions(om->GetLargestPossibleRegion());
	seed_out->Allocate();
	seed_out->FillBuffer(0);

	int thresh_value = 6;
	int number_of_seeds = 200;
	int tnum_seeds = 0;

	typedef itk::ImageRegionIteratorWithIndex<InputImageType> IndexIteratorType;

	IndexIteratorType it1(minfilt->GetOutput(),minfilt->GetOutput()->GetLargestPossibleRegion());
	IteratorType it2(om,om->GetLargestPossibleRegion());

	for(it2.GoToBegin();!it2.IsAtEnd(); ++it2)
	{
		if(it2.Get()>thresh_value)
			tnum_seeds++;
	}
	printf("tnum_seeds = %d\n",tnum_seeds);
	IteratorType it3(seed_out,seed_out->GetLargestPossibleRegion());
	IteratorType it4(im,im->GetLargestPossibleRegion());

	std::vector<mdl::fPoint3D> seeds;
	seeds.clear();
	/*for(it1.GoToBegin(),it2.GoToBegin(),it3.GoToBegin(),it4.GoToBegin();!it1.IsAtEnd();++it1,++it2,++it3,++it4)
	{
		if(it1.Get()==it2.Get() && it4.Get() > 150)
		{
			it3.Set(255);
			InputImageType::IndexType index = it1.GetIndex();
			mdl::fPoint3D seed1;
			seed1.x = index[0];
			seed1.y = index[1];
			seed1.z = index[2];
			seeds.push_back(seed1);
		}
	}*/

	seeds.clear();
	
	while(!feof(fp))
	{
		mdl::fPoint3D seed1;
		fscanf(fp,"%f %f %f",&seed1.x,&seed1.y,&seed1.z);
		if(feof(fp))
			break;
		seed1.x*=1;
		seed1.y*=1;
		seeds.push_back(seed1);
	}
	fclose(fp);
	printf("Seeds.size = %d\n",seeds.size());
	//scanf("%*d");
	mdl::vtkFileHandler * fhd1 = new mdl::vtkFileHandler();
	fhd1->SetNodes(&seeds);
	std::vector<mdl::pairE> nullpairs;
	fhd1->SetLines(&nullpairs);
	std::string outFilename1 = "C:/Users/arun/Research/Farsight/exe/bin/mst_input.vtk";
	fhd1->Write(outFilename1.c_str());
	delete fhd1;


	int edgeRange = 50;
	int morphStrength = 0;
	mdl::MST *mst = new mdl::MST( im );
	mst->SetDebug(false);
	mst->SetUseVoxelRounding(false);
	mst->SetEdgeRange(edgeRange);
	mst->SetPower(1);
	mst->SetSkeletonPoints( &seeds );
	// can choose different weight
	//mst->CreateGraphAndMST(1);
	mst->CreateGraphAndMST(5);
	mst->ErodeAndDialateNodeDegree(morphStrength);

	std::vector<mdl::fPoint3D> nodes = mst->GetNodes();
	std::vector<mdl::pairE> bbpairs = mst->BackboneExtract();

	delete mst;

	std::cerr << "Saving\n";

	//****************************************************************
	// TREE WRITER
	mdl::vtkFileHandler * fhd2 = new mdl::vtkFileHandler();
	fhd2->SetNodes(&nodes);
	fhd2->SetLines(&bbpairs);
	std::string outFilename2 = "C:/Users/arun/Research/Farsight/exe/bin/mst_tree.vtk";
	fhd2->Write(outFilename2.c_str());
	delete fhd2;
	scanf("%*d");

	writeImage<InputImageType>(seed_out,"C:/Users/arun/Research/Farsight/exe/bin/seedimage.tif");
		
}