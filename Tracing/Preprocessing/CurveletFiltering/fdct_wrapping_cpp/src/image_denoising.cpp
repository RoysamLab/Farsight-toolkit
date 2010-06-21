/*
   Copyright (C) 2004 Caltech
   Written by Lexing Ying
*/

#include "fdct_wrapping.hpp"
#include "fdct_wrapping_inline.hpp"

#include <stdlib.h>
#include <time.h>
using namespace std;
using namespace fdct_wrapping_ns;

//itk includes
#include "itkImage.h"
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

//typedefs
typedef unsigned char InputPixelType;
typedef itk::Image<InputPixelType,3> InputImageType;
typedef itk::Image<InputPixelType,2> Input2DImageType;
typedef itk::Image<float,3> FloatImageType;
typedef itk::Image<float,2> Float2DImageType;
  
int nbscales;
int nbangles_coarse;
int ac;

//function definitions
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



template <typename PixelType>
void copyslice(typename itk::Image<PixelType,2>::Pointer im1, typename itk::Image<PixelType,3>::Pointer im2, int slice)
{
	PixelType* inpointer  = im1->GetBufferPointer();
	PixelType* outpointer = im2->GetBufferPointer();
	typename itk::Image<PixelType,2>::SizeType size = im1->GetLargestPossibleRegion().GetSize();
	int copysize = size[0]*size[1];
	memcpy(outpointer+copysize*slice,inpointer,sizeof(PixelType)*copysize);
}
int optionsCreate(const char* optfile, map<string,string>& options)
{
  options.clear();
  ifstream fin(optfile); assert(fin.good());
  string name;  fin>>name;
  while(fin.good()) {
	 char cont[100];	 fin.getline(cont, 99);
	 options[name] = string(cont);
	 fin>>name;
  }
  fin.close();
  return 0;
}






void getCurvetlets(Input2dImageType::Pointer im,Input2dImageType::Pointer &om,Float2dImageType::Pointer &cosim, Float2dImageType::Pointer &sinim)
{
	clock_t ck0, ck1;


	Input2DImageType::SizeType imsize = im->GetLargestPossibleRegion().GetSize();

	int m,n;
	int m = int(pow(2,ceil(log(imsize[1])/log(2.0)))+0.5);
	int n = int(pow(2,ceil(log(imsize[0])/log(2.0)))+0.5);

	CpxNumMat input(n,m);

	typedef itk::ImageRegionIteratorWithIndex<Input2DImageType> IterType1;
	ck0 = clock();
	IterType1 iter1(im,im->GetLargestPossibleRegion());
	std::complex<double> cpxtemp;
	for(iter1.GoToBegin(); !iter1.IsAtEnd(); ++iter1)
	{
		Input2DImageType::IndexType index;
		index = iter1.GetIndex();
		cpxtemp.real = iter1.Get();
		cpxtemp.imag = 0;
		input(index[1],index[0])=cpxtemp;
	}
	
	//compute the fdct wrapping for the F

	CpxNumMat initF(n,m);
	CpxNumMat F(n,m);
	cpxtemp.real = sqrt(n*m);
	cpxtemp.imag = 0;
	initF(0,0) = cpxtemp;
	fdct_wrapping_fftwshift(initF,F);


	//fdct_wrapping_
	vector< vector<CpxNumMat> > c;  //vector<int> extra;
	fdct_wrapping(m, n, nbscales, nbangles_coarse, ac, F, c);
	ck1 = clock();  cout<<"FDCT_WRAPPING_  takes "<<double(ck1-ck0)/CLOCKS_PER_SEC<<" seconds"<<endl;  ck0 = ck1;

	vector< vector<double> > E;	
	for(int cx = 0; cx< c.size();cx++)
	{
		vector<double> row;
		for(int cy = 0;cy<cx.size(); cy++)
		{
			double sum = 0;
			int size1 = c[cx][cy].m();
			int size2 = c[cx][cy].n();
			for(int counter = 0; counter < size1; counter++)
			{
				for(int counter1 =0; counter1 < size2; counter1++)
				{
					sum = sum + abs(c[cx][cy](counter,counter1));
				}
			}
			row.push_back(sqrt(sum/size1/size2))
		}
		E.push_back(row);
	}


	//ifdct_wrapping_
	CpxNumMat y(x); clear(y);
	ifdct_wrapping(m, n, nbscales, nbangles_coarse, ac, c, y);
	ck1 = clock();  cout<<"IFDCT_WRAPPING_ takes "<<double(ck1-ck0)/CLOCKS_PER_SEC<<" seconds"<<endl;  ck0 = ck1;

	CpxNumMat e(m,n);
	for(int i=0; i<m; i++)
		for(int j=0; j<n; j++)
			e(i,j) = x(i,j) - y(i,j);
	cerr<<"accuracy of inversion "<<sqrt(energy(e)/(m*n))<<endl;

	vector< vector<double> > sx, sy;
	vector< vector<double> > fx, fy;
	vector< vector<int> > nx, ny;
	fdct_wrapping_param(m, n, nbscales, nbangles_coarse, ac, sx, sy, fx, fy, nx, ny);
}
int main(int argc, char** argv)
{

  
  //assert(argc==2);
  //get options
  map<string, string> opts;  optionsCreate("options", opts);
  
  //get input data
  map<string,string>::iterator mi;
  

  //mi = opts.find("-m"); assert(mi!=opts.end());
  //{ istringstream ss((*mi).second); ss>>m; }

  //mi = opts.find("-n"); assert(mi!=opts.end());
  //{ istringstream ss((*mi).second); ss>>n; }
  

  mi = opts.find("-nbscales"); assert(mi!=opts.end());
  { istringstream ss((*mi).second); ss>>nbscales; }
  

  mi = opts.find("-nbangles_coarse"); assert(mi!=opts.end());
  { istringstream ss((*mi).second); ss>>nbangles_coarse; }
  

  mi = opts.find("-ac"); assert(mi!=opts.end());
  { istringstream ss((*mi).second); ss>>ac; }
  

  InputImageType::Pointer inputim = readImage<InputImageType>("Lena.jpg");
  
  int slices = inputim->GetLargestPossibleRegion().GetSize()[2];
  InputImageType::Pointer outputim = InputImageType::New();
  outputim->SetRegions(inputim->GetLargestPossibleRegion());
  outputim->Allocate();
  FloatImageType::Pointer cosim = FloatImageType::New();
  cosim->SetRegions(inputim->GetLargestPossibleRegion());
  cosim->Allocate();
  FloatImageType::Pointer sinim = FloatImageType::New();
  sinim->SetRegions(inputim->GetLargestPossibleRegion());
  sinim->Allocate();

  for(int counter = 0; counter < slices; counter++)
  {
	  Input2DImageType::Pointer im2d = getSlice(inputim,counter);
	  Input2DImageType::Pointer om2d;
	  Float2DImageType::Pointer cosim2d,sinim2d;
	  //call single slice 2-d curvelets function
	  getCurvelets(im2d,om2d,cosim2d,sinim2d);
	  copyslice<InputPixelType>(om2d,outputim,counter);
	  copyslice<float>(cosim2d,cosim,counter);
	  copyslice<float>(sinim2d,sinim,counter);
  }

	writeImage<InputImageType>(outputim,"LenaOut.tif");	

  
  return 0;
}

