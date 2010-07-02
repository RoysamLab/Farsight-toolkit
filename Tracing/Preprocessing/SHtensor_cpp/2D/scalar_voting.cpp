#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <iostream>


using namespace std;


//itk includes
#include "itkImage.h"
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionConstIterator.h>
#include "vnl/vnl_vector.h"
#include "vnl/algo/vnl_convolve.txx"


//typedefs
typedef unsigned char InputPixelType;
typedef std::complex<float> Cpx;
typedef itk::Image<InputPixelType,3> InputImageType;
typedef itk::Image<InputPixelType,2> Input2DImageType;
typedef itk::Image<float,3> FloatImageType;
typedef itk::Image<float,2> Float2DImageType;
typedef itk::Image<Cpx,2> Cpx2DImageType;

float sigma;
float order;

template <typename T>typename T::Pointer readImage(const char *filename){	printf("Reading %s ... \n",filename);	typedef typename itk::ImageFileReader<T> ReaderType;	typename ReaderType::Pointer reader = ReaderType::New();	ReaderType::GlobalWarningDisplayOff();	reader->SetFileName(filename);	try	{		reader->Update();	}	catch(itk::ExceptionObject &err)	{		std::cerr << "ExceptionObject caught!" <<std::endl;		std::cerr << err << std::endl;		_exit(0);		//return EXIT_FAILURE;	}
	printf("Done\n");	return reader->GetOutput();}
template <typename T>int writeImage(typename T::Pointer im, const char* filename){	printf("Writing %s ... \n",filename);	typedef typename itk::ImageFileWriter<T> WriterType;	typename WriterType::Pointer writer = WriterType::New();	writer->SetFileName(filename);	writer->SetInput(im);	try	{		writer->Update();	}	catch(itk::ExceptionObject &err)	{		std::cerr << "ExceptionObject caught!" <<std::endl;		std::cerr << err << std::endl;		return EXIT_FAILURE;	}	return EXIT_SUCCESS;}

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
int main(int argc, char**argv)
{
	clock_t ck0, ck1;
	ck0 = clock();
	if(argc < 2 || argc > 5)
	{
		printf("Usage : shtensor_scalar_voting_2d.exe gradmag_file dircos_file dirsin_file [options_file]");
		return -1;
	}
	//assert(argc==3);
	//get options

	map<string, string> opts;  
	if(argc == 5)
		optionsCreate(argv[4], opts);

	map<string,string>::iterator mi;

	mi = opts.find("-sigma"); 
	if(mi!=opts.end())
	{ istringstream ss((*mi).second); ss>>sigma; }
	else
	{ sigma = 5; printf("Chose sigma = 5 as default\n");}

	mi = opts.find("-order"); 
	if(mi!=opts.end())
	{ istringstream ss((*mi).second); ss>>order; }
	else
	{ order = 4; printf("Chose order = 4 as default\n");}


	Input2DImageType::Pointer gradmag = readImage<Input2DImageType>(argv[1]);
	Float2DImageType::Pointer cosim = readImage<Float2DImageType>(argv[2]);
	Float2DImageType::Pointer sinim = readImage<Float2DImageType>(argv[3]);

	Float2DImageType::Pointer Axx,Axy,Ayy,stickness, orientation;
	Cpx2DImageType::Pointer A0,A2,A2n;

	Axx = Float2DImageType::New();
	Axy = Float2DImageType::New();
	Ayy = Float2DImageType::New();
	stickness = Float2DImageType::New();
	orientation = Float2DImageType::New();
	A0 = Cpx2DImageType::New();
	A2 = Cpx2DImageType::New();
	A2n = Cpx2DImageType::New();

	Axx->SetRegions(gradmag->GetLargestPossibleRegion()); Axx->Allocate();
	Axy->SetRegions(gradmag->GetLargestPossibleRegion()); Axy->Allocate();
	Ayy->SetRegions(gradmag->GetLargestPossibleRegion()); Ayy->Allocate();
	stickness->SetRegions(gradmag->GetLargestPossibleRegion()); stickness->Allocate();
	orientation->SetRegions(gradmag->GetLargestPossibleRegion()); orientation->Allocate();
	A0->SetRegions(gradmag->GetLargestPossibleRegion()); A0->Allocate();
	A2->SetRegions(gradmag->GetLargestPossibleRegion()); A2->Allocate();
	A2n->SetRegions(gradmag->GetLargestPossibleRegion()); A2n->Allocate();

	typedef itk::ImageRegionIterator<Float2DImageType> Float2DIteratorType;
	typedef itk::ImageRegionIterator<Input2DImageType> Input2DIteratorType;

	Input2DIteratorType iter1(gradmag,gradmag->GetLargestPossibleRegion());iter1.GoToBegin();
	Float2DIteratorType iter2(cosim,sinim->GetLargestPossibleRegion());iter2.GoToBegin();
	Float2DIteratorType iter3(sinim,sinim->GetLargestPossibleRegion());iter3.GoToBegin();
	Float2DIteratorType iter4(Axx,Axx->GetLargestPossibleRegion());iter4.GoToBegin();
	Float2DIteratorType iter5(Axy,gradmag->GetLargestPossibleRegion());iter5.GoToBegin();
	Float2DIteratorType iter6(Ayy,gradmag->GetLargestPossibleRegion());iter6.GoToBegin();
	Float2DIteratorType iter7(stickness,gradmag->GetLargestPossibleRegion());iter7.GoToBegin();
	Float2DIteratorType iter8(orientation,gradmag->GetLargestPossibleRegion());iter8.GoToBegin();

	typedef itk::ImageRegionIterator<Cpx2DImageType> Cpx2DIteratorType;

	Cpx2DIteratorType iter9(A0,gradmag->GetLargestPossibleRegion());iter9.GoToBegin();
	Cpx2DIteratorType iter10(A2,gradmag->GetLargestPossibleRegion());iter10.GoToBegin();
	Cpx2DIteratorType iter11(A2n,gradmag->GetLargestPossibleRegion());iter11.GoToBegin();


	for(;!iter1.IsAtEnd();)
	{
		float Ix = iter1.Get()*iter2.Get();
		float Iy = -iter1.Get()*iter3.Get();
		float Axx = Ix*Ix;
		float Axy = Ix*Iy;
		float Ayy = Iy*Iy;

		iter4.Set(Axx);
		iter5.Set(Axy);
		iter6.Set(Ayy);
		iter7.Set(sqrt(float((Axx+Ayy)*(Axx+Ayy) - 4*((Axx*Ayy) - Axy*Axy))));//stickness
		iter8.Set(fmod(atan2(float(iter3.Get()),float(iter2.Get())),2*acos(float(0))));//orienation
		iter9.Set(Cpx(Axx+Ayy,0));
		iter10.Set(Cpx(Axx-Ayy, -2*Axy));
		iter11.Set(Cpx(Axx-Ayy, 2*Axy));
		++iter1;++iter2,++iter3;++iter4;++iter5;++iter6;++iter7;++iter8;++iter9;++iter10;++iter11;
	}

	int K = sigma*4.0;

	vnl_matrix<float> mat(2*K+1,2*K+1);


	for(int xco = -K; xco<=K; xco++)
	{
		for(int yco = -K; yco<=K; yco++)
		{
			float norm = sqrt(xco*xco+yco*yco);
			if(norm!=0)
				mat[K+xco][K+yco] = fmod(atan2(yco,xco),2*acos(float(0)));
			else
				mat[K+xco][K+yco] = 0;
		}
	}

	vnl_vector<int> a(3),b(3);
	a[0] = 1;a[1] =2; a[2] = 1;
	b = a;
	for(int co = 1; co< order -1; co++)
		a = vnl_convolve(a,b,0);
	a = a.extract((a.size()+1)/2,0);
	a = 2*a;
	a[a.size()-1] = a[a.size()-1]/2;
	a.flip();
	std::vector< vnl_matrix< Cpx > > W,C;

	int rows = gradmag->GetLargestPossibleRegion().GetSize[1];
	int cols = gradmag->GetLargestPossibleRegion().GetSize[0];

	vnl_matrix < unsigned char > U0(rows,cols);
	int pc = 0;

	for(int co = 0; co <= 4*order-2; co = co + 2)
	{
		vnl_matrix< Cpx > wvals(2*K+1,2*K+1);
		vnl_matrix< Cpx > cvals(rows,cols);
		for(int xco = -K; xco<=K; xco++)
		{
			for(int yco = -K; yco<=K; yco++)
			{
				float norm = sqrt(xco*xco+yco*yco);
				wvals[xco+K][yco+K] = exp(-norm/2.0/sigma/sigma)*Cpx(cos(-mat[xco+K][yco+K]*co),sin(-mat[xco+K][yco+K]*co));
			}
		}
		typedef itk::ImageRegionIteratorWithIndex<Float2DImageType> IterType1;
		IterType1 iter1(stickness,stickness->GetLargestPossibleRegion());
		Float2DIteratorType iter2(orientation,orientation->GetLargestPossibleRegion());
		iter1.GoToBegin();
		iter2.GoToBegin();
		for(;!iter1.IsAtEnd();++iter1,++iter2)
		{
			Float2DImageType::IndexType index;
			index = iter1.Get();
			cvals(index[1],index[0]) = stickness*Cpx(cos(-co*iter2.Get()),sin(-co*iter2.Get()));
		}
		vnl_matrix<Cpx> Utemp = vnl_convolve(cvals,wvals,1);
		for(int xco = K; xco < rows + K ; xco++)
		{
			for(int yco = K; yco < cols + K; yco++)
			{
				U0[xco-K][yco-K] += Utemp[xco][yco];
			}
		}
		pc++;
		if(pc > a.size()-1)
		{
			break;
		}
		//W.push_back(vals);
	}

	Input2DImageType::Pointer out = Input2DImageType::New();
	out->SetRegions(gradmag->GetLargestPossibleRegion());
	out->Allocate();

	typedef itk::ImageRegionIteratorWithIndex<Input2DImageType> IteratorType2;
	IteratorType2 outiter(out,out->GetLargestPossibleRegion());
	outiter.GoToBegin();
	for(;!outiter.IsAtEnd();++outiter)
	{
		Input2DImageType::IndexType index = outiter.GetIndex();
		outiter.Set(U0[index[1]][index[0]]);
	}
	writeImage<Input2DImageType>(out,"test_out.tif");
	return 0;
}
