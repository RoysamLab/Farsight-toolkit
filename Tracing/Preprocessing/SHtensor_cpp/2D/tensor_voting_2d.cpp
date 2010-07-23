#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <iostream>
#include <iomanip>
#include <fftw3.h>

using namespace std;


//itk includes
#include "itkImage.h"
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionConstIterator.h>
#include "itkCastImageFilter.h"
#include "vnl/vnl_vector.h"
#include "vnl/algo/vnl_convolve.h"
#include "vnl/algo/vnl_fft_2d.h"

#define MAX(a,b) (((a) > (b))? (a) : (b))
#define MIN(a,b) (((a) < (b))? (a) : (b))

//typedefs
typedef unsigned char InputPixelType;
typedef std::complex<float> Cpx;
typedef itk::Image<InputPixelType,3> InputImageType;
typedef itk::Image<InputPixelType,2> Input2DImageType;
typedef itk::Image<float,3> FloatImageType;
typedef itk::Image<float,2> Float2DImageType;
typedef itk::Image< Cpx ,2 > Cpx2DImageType;

template <typename T> int writeImage(typename T::Pointer im, const char* filename){
	printf("Writing %s ... \n",filename);	
	typedef typename itk::ImageFileWriter< T > WriterType;	
	typename WriterType::Pointer writer = WriterType::New();	
	writer->SetFileName(filename);	writer->SetInput(im);	
	try	{
		writer->Update();	
	}	
	catch(itk::ExceptionObject &err)	{	
		std::cerr << "ExceptionObject caught!" <<std::endl;	
		std::cerr << err << std::endl;	
		return EXIT_FAILURE;	
	}	
	return EXIT_SUCCESS;
}


template <typename T> typename T::Pointer readImage(const char *filename)
{	
	printf("Reading %s ... \n",filename);	
	typedef typename itk::ImageFileReader<T> ReaderType;
	typename ReaderType::Pointer reader = ReaderType::New();
	ReaderType::GlobalWarningDisplayOff();	reader->SetFileName(filename);	
	try	{	
		reader->Update();	
	}	
	catch(itk::ExceptionObject &err)	{
		std::cerr << "ExceptionObject caught!" <<std::endl;
		std::cerr << err << std::endl;
		_exit(0);		
		//return EXIT_FAILURE;	
	}
	printf("Done\n");	return reader->GetOutput();
}

int optionsCreate(const char* optfile, map<string,string>& options)
{
	options.clear();
	ifstream fin(optfile); assert(fin.good());
	string name;  fin>>name;
	while(fin.good()) 
	{
		char cont[100];	 fin.getline(cont, 99);
		options[name] = string(cont);
		fin>>name;
	}
	fin.close();
	return 0;
}
/*template <typename T1, typename T2>
typename itk::Image< T1,2>::Pointer vnlmat_to_itkimage(vnl_matrix< typename T2 > mat,float scale = 1.0)
{
	typedef typename itk::Image<T1,2> ImageType;
	ImageType::SizeType size;
	ImageType::IndexType index;
	ImageType::RegionType region;
	index.Fill(0);
	size[0] = mat.cols();
	size[1] = mat.rows();
	region.SetSize(size);
	region.SetIndex(index);
	ImageType::Pointer im = ImageType::New();
	im->SetRegions(region);
	im->Allocate();
	for(int xco = 0; xco<size[1]; xco++)
	{
		for(int yco = 0; yco < size[0]; yco++)
		{
			index[1] = xco;
			index[0] = yco;
			im->SetPixel(index,mat[xco][yco]*scale);
		}
	}
	return im;

}*/

vnl_matrix<float> cpx_to_abs(vnl_matrix< std::complex< float> > mat)
{
	vnl_matrix<float> out(mat.rows(),mat.cols());

	for(int x= 0; x< mat.rows(); x++)
	{
		for(int y = 0; y < mat.cols(); y++)
		{
			out[x][y] = abs(mat[x][y]);
		}
	}
	return out;
}

vnl_matrix< std::complex <float> > myconvolve( vnl_matrix< std::complex<float> > &mat1, vnl_matrix< std::complex <float> > &mat2)
{
	//printf("mat1.size() = %d %d\n",mat1.rows(),mat1.cols());
	//printf("mat2.size() = %d %d\n",mat2.rows(),mat2.cols());


	typedef itk::Image<unsigned short, 2> UShortImageType; 
	//writeImage<UShortImageType>(vnlmat_to_itkimage<unsigned short,float>(cpx_to_abs(mat1)),"testmat1.tif");
	//writeImage<UShortImageType>(vnlmat_to_itkimage<unsigned short,float>(cpx_to_abs(mat2),255),"testmat2.tif");
	int tsize = mat1.rows()*mat1.cols()*sizeof(fftwf_complex);
	fftwf_complex *in = (fftwf_complex*)fftwf_malloc(tsize);
	fftwf_complex *out = (fftwf_complex*)fftwf_malloc(tsize);

	printf("Started planning..\n");
	fftwf_plan plan;
#pragma omp critical
	{
		plan = fftwf_plan_dft_2d(mat1.rows(),mat1.cols(),in, out,FFTW_FORWARD,FFTW_MEASURE);
	}
	printf("Finished planning \n");
	memcpy(in,mat1.data_block(),tsize);
	fftwf_execute_dft(plan,in,out);
	printf("executed fftwf\n");
	//fftwf_destroy_plan(plan);
	memcpy(mat1.data_block(),out,tsize);
	//fftwf_free(in);
	//fftwf_free(out);

	//in = (fftwf_complex*)fftwf_malloc(tsize);
	//out = (fftwf_complex*)fftwf_malloc(tsize);
	//plan = fftwf_plan_dft_2d(mat1.rows(),mat1.cols(),in, out,FFTW_FORWARD,FFTW_MEASURE);
	memcpy(in,mat2.data_block(),tsize);
	fftwf_execute_dft(plan,in,out);
	printf("executed fftwf\n");
	fftwf_destroy_plan(plan);
	memcpy(mat2.data_block(),out,tsize);
	//fftwf_free(in);
	//fftwf_free(out);


	vnl_matrix< std::complex < float > > mat3(mat1.rows(),mat1.cols());
	float factor = mat1.rows()*mat1.cols();

	for(int xco = 0; xco < mat1.rows(); ++xco)
	{
		for(int yco = 0; yco < mat1.cols(); ++yco)
		{
			mat3[xco][yco] = mat1[xco][yco]/factor * mat2[xco][yco];

		}
	}

	//in = (fftwf_complex*)fftwf_malloc(tsize);
	//out = (fftwf_complex*)fftwf_malloc(tsize);
#pragma omp critical
	{
		plan = fftwf_plan_dft_2d(mat1.rows(),mat1.cols(),in, out,FFTW_BACKWARD,FFTW_MEASURE);
	}
	memcpy(in,mat3.data_block(),tsize);
	fftwf_execute_dft(plan,in,out);
	printf("executed fftwf\n");
	fftwf_destroy_plan(plan);
	memcpy(mat3.data_block(),out,tsize);
	//writeImage<UShortImageType>(vnlmat_to_itkimage<unsigned short,float>(cpx_to_abs(mat3)),"testmat3.tif");
	fftwf_free(in);
	fftwf_free(out);

	//std::cout<<mat3;


	return mat3;

}


template <typename T1, typename T2>
typename T1::Pointer castImage(typename T2::Pointer im)
{
	typedef typename itk::CastImageFilter<T2,T1> CastImageFilter;
	typename CastImageFilter::Pointer castf = CastImageFilter::New();
	castf->SetInput(im);
	castf->Update();
	return castf->GetOutput();
}







template <typename PixelType>
typename itk::Image<PixelType, 2>::Pointer getSlice(typename itk::Image<PixelType,3>::Pointer im, int slice)
{
	typedef typename itk::Image<PixelType,2> ImageType2D;
	typedef typename itk::Image<PixelType,2> ImageType3D;
	typename ImageType2D::Pointer out = ImageType3D::New();
	typename ImageType2D::SizeType size;
	size[0] = im->GetLargestPossibleRegion().GetSize()[0];
	size[1] = im->GetLargestPossibleRegion().GetSize()[1];
	typename Input2DImageType::IndexType index;
	index.Fill(0);
	typename ImageType2D::RegionType region;
	region.SetSize(size);
	region.SetIndex(index);
	out->SetRegions(region);
	out->Allocate();
	if(out->GetBufferPointer()==NULL)
		printf("Could not allocate memory -1 ... I'm going to crash any moment now.. \n");
	memcpy(out->GetBufferPointer(),im->GetBufferPointer()+slice*size[0]*size[1],size[0]*size[1]*sizeof(PixelType));
	return out;
}


template <typename PixelType>
void copyslice(typename itk::Image<PixelType,2>::Pointer im1, typename itk::Image<PixelType,3>::Pointer im2, int slice)
{
	printf("Copying slice\n");
	PixelType* inpointer  = im1->GetBufferPointer();
	PixelType* outpointer = im2->GetBufferPointer();
	typename itk::Image<PixelType,2>::SizeType size = im1->GetLargestPossibleRegion().GetSize();
	int copysize = size[0]*size[1];
	memcpy(outpointer+copysize*slice,inpointer,sizeof(PixelType)*copysize);
}
int main(int argc, char**argv)
{
	float sigma;
	float order;
	time_t t1,t2;
	clock_t ck0, ck1;
	ck0 = clock();
	time(&t1);
	if(argc < 2 || argc > 5)
	{
		printf("Usage : tensor_voting_2d.exe gradmag_file dircos_file dirsin_file [options_file]");
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
	int numt;
	mi = opts.find("-num_threads"); 
	if(mi!=opts.end())
	{ istringstream ss((*mi).second); ss>>numt; }
	else
	{ numt = 8; printf("Chose num_threads = 8 as default \n"); }


	InputImageType::Pointer ingradmag = readImage<InputImageType>(argv[1]);
	FloatImageType::Pointer incosim = readImage<FloatImageType>(argv[2]);
	FloatImageType::Pointer insinim = readImage<FloatImageType>(argv[3]);

	InputImageType::Pointer realout = InputImageType::New();
	realout->SetRegions(ingradmag->GetLargestPossibleRegion());
	realout->Allocate();
	int num_slices = ingradmag->GetLargestPossibleRegion().GetSize()[2];
	std::vector<float> factors(num_slices);

#pragma omp parallel for  num_threads(numt)
	for(int sliceco = 0; sliceco < num_slices; sliceco++)
	{

		printf("SLICE %d\n", sliceco);

		Input2DImageType::Pointer gradmag = getSlice<unsigned char>(ingradmag, sliceco);//readImage<Input2DImageType>(argv[1]);
		Float2DImageType::Pointer cosim = getSlice<float>(incosim,sliceco); //readImage<Float2DImageType>(argv[2]);
		Float2DImageType::Pointer sinim = getSlice<float>(insinim,sliceco);//readImage<Float2DImageType>(argv[3]);

		Float2DImageType::Pointer Axx,Axy,Ayy,stickness, orientation;
		Cpx2DImageType::Pointer A0,A2,A2n;

		//	Axx = Float2DImageType::New();
		//	Axy = Float2DImageType::New();
		//	Ayy = Float2DImageType::New();
		stickness = Float2DImageType::New();
		orientation = Float2DImageType::New();
		//	A0 = Cpx2DImageType::New();
		//	A2 = Cpx2DImageType::New();
		//	A2n = Cpx2DImageType::New();

		//	Axx->SetRegions(gradmag->GetLargestPossibleRegion()); Axx->Allocate();
		//	Axy->SetRegions(gradmag->GetLargestPossibleRegion()); Axy->Allocate();
		//	Ayy->SetRegions(gradmag->GetLargestPossibleRegion()); Ayy->Allocate();
		stickness->SetRegions(gradmag->GetLargestPossibleRegion()); stickness->Allocate();
		orientation->SetRegions(gradmag->GetLargestPossibleRegion()); orientation->Allocate();
		//	A0->SetRegions(gradmag->GetLargestPossibleRegion()); A0->Allocate();
		//	A2->SetRegions(gradmag->GetLargestPossibleRegion()); A2->Allocate();
		//	A2n->SetRegions(gradmag->GetLargestPossibleRegion()); A2n->Allocate();

		typedef itk::ImageRegionIterator<Float2DImageType> Float2DIteratorType;
		typedef itk::ImageRegionIterator<Input2DImageType> Input2DIteratorType;

		Input2DIteratorType iter1(gradmag,gradmag->GetLargestPossibleRegion());iter1.GoToBegin();
		Float2DIteratorType iter2(sinim,sinim->GetLargestPossibleRegion());iter2.GoToBegin();
		Float2DIteratorType iter3(cosim,sinim->GetLargestPossibleRegion());iter3.GoToBegin();
		//	Float2DIteratorType iter4(Axx,Axx->GetLargestPossibleRegion());iter4.GoToBegin();
		//	Float2DIteratorType iter5(Axy,gradmag->GetLargestPossibleRegion());iter5.GoToBegin();
		//	Float2DIteratorType iter6(Ayy,gradmag->GetLargestPossibleRegion());iter6.GoToBegin();
		Float2DIteratorType iter7(stickness,gradmag->GetLargestPossibleRegion());iter7.GoToBegin();
		Float2DIteratorType iter8(orientation,gradmag->GetLargestPossibleRegion());iter8.GoToBegin();

		typedef itk::ImageRegionIterator< Cpx2DImageType > Cpx2DIteratorType;

		//	Cpx2DIteratorType iter9(A0,gradmag->GetLargestPossibleRegion());iter9.GoToBegin();
		//	Cpx2DIteratorType iter10(A2,gradmag->GetLargestPossibleRegion());iter10.GoToBegin();
		//	Cpx2DIteratorType iter11(A2n,gradmag->GetLargestPossibleRegion());iter11.GoToBegin();

		//printf("Begin iterators ...\n");
		for(;!iter1.IsAtEnd();)
		{
			float cosv = iter2.Get();
			float sinv = iter3.Get();
			if(cosv != cosv)
				cosv = 1;
			if(sinv != sinv)
				sinv = 0;
			float Ix = iter1.Get()*1.0f*cosv;
			float Iy = -iter1.Get()*1.0f*sinv;
			//printf("Ix = %f Iy = %f\n",Ix,Iy);
			float Axx = Ix*Ix;
			float Axy = Ix*Iy;
			float Ayy = Iy*Iy;

			//iter4.Set(Axx);
			//iter5.Set(Axy);
			//iter6.Set(Ayy);
			iter7.Set(sqrt(float((Axx+Ayy)*(Axx+Ayy) - 4*(Axx*Ayy - Axy*Axy))));//stickness
			iter8.Set(fmod(atan2(float(sinv),float(cosv)),2*acos(float(0))));//orienation
			//iter9.Set(Cpx(Axx+Ayy,0));
			//iter10.Set(Cpx(Axx-Ayy, -2*Axy));
			//iter11.Set(Cpx(Axx-Ayy, 2*Axy));
			++iter1;++iter2,++iter3;
			//++iter4;++iter5;++iter6;
			++iter7;++iter8;
			//++iter9;++iter10;++iter11;
		}
		//printf("Done with copying..\n");

		typedef itk::Image<unsigned short, 2> UShort2DImageType;
		/*writeImage<UShort2DImageType>(castImage<UShort2DImageType,Float2DImageType>(stickness),"stickness.tif");
		Input2DImageType::Pointer out_orient =  Input2DImageType::New();
		out_orient->SetRegions(orientation->GetLargestPossibleRegion());
		out_orient->Allocate();
		Input2DIteratorType iter_o(out_orient,out_orient->GetLargestPossibleRegion());
		Float2DIteratorType iter_(orientation,orientation->GetLargestPossibleRegion());
		for(iter_o.GoToBegin(),iter_.GoToBegin(); !iter_o.IsAtEnd();++iter_o,++iter_)
		{
		iter_o.Set(iter_.Get()*180/2/acos(float(0)));
		}
		writeImage<Input2DImageType>(out_orient,"orientation.tif");
		*/
		int K = sigma*4.0;

		vnl_matrix<float> mat(2*K+1,2*K+1);

		//printf("Creating mat...\n");
		for(int xco = -K; xco<=K; xco++)
		{
			for(int yco = -K; yco<=K; yco++)
			{
				float norm = sqrt(float(xco*xco+yco*yco));
				if(norm!=0)
					mat[K+xco][K+yco] = fmod(atan2(float(yco/norm),float(xco/norm)),2*acos(float(0)));
				else
					mat[K+xco][K+yco] = 0;
			}
		}
		//printf("Done creating mat...\n");
		//printf("mat.size() = %d %d\n",mat.rows(),mat.cols());
		//std::cout<<setprecision(2)<<mat;
		vnl_vector<int> a(3),b(3);
		a[0] = 1;a[1] =2; a[2] = 1;
		b = a;
		for(int co = 1; co< order -1; co++)
			a = vnl_convolve(a,b);
		vnl_vector<int> acopy  = a;
		a = a.extract((a.size()+1)/2,0);
		a = 2*a;
		a[a.size()-1] = a[a.size()-1]/2;
		a.flip();
		std::vector< vnl_matrix< Cpx > > W,C;
		//std::cout<<acopy;
		int rows = gradmag->GetLargestPossibleRegion().GetSize()[1];
		int cols = gradmag->GetLargestPossibleRegion().GetSize()[0];

		vnl_matrix < std::complex<float> > U2n(rows,cols,std::complex<float>(0,0));
		int pc = 0;
		//printf("Entering main loop...\n");
		for(int co = 0; co <= 4*order-2; co = co + 2)
		{
			int xsize_t = 2*K+1+cols+1;
			int ysize_t = 2*K+1+rows+1;
			int xsize = 1; while(xsize< xsize_t) xsize <<=1;
			int ysize = 1; while(ysize< ysize_t) ysize <<=1;
			//xsize = xsize_t;
			//ysize = ysize_t;
			vnl_matrix< Cpx > wvals(ysize,xsize,0);
			vnl_matrix< Cpx > cvals(ysize,xsize,0);

			for(int xco = -K; xco<=K; xco++)
			{
				for(int yco = -K; yco<=K; yco++)
				{
					float norm = sqrt(float(xco*xco+yco*yco));
					float lval = exp(-norm*norm/2.0/sigma/sigma);
					wvals[xco+K][yco+K] = Cpx(lval*cos(-mat[xco+K][yco+K]*co),lval*sin(-mat[xco+K][yco+K]*co));
				}
			}

			typedef itk::ImageRegionIteratorWithIndex< Float2DImageType > IterType1;
			IterType1 iter1(stickness,stickness->GetLargestPossibleRegion());
			Float2DIteratorType iter2(orientation,orientation->GetLargestPossibleRegion());
			iter1.GoToBegin();
			iter2.GoToBegin();
			for(;!iter1.IsAtEnd();++iter1,++iter2)
			{
				Float2DImageType::IndexType index;
				index = iter1.GetIndex();
				cvals[index[1]][index[0]] = iter1.Get()*Cpx(cos(-(co-2)*iter2.Get()),sin(-(co-2)*iter2.Get())); // NOTE THE CHANGE HERE
			}
			//printf("Copied wvals and cvals... entering myconvolve\n");
			//std::cout<<setprecision(2)<<cvals.extract(10,10);
			vnl_matrix< Cpx > Utemp = myconvolve(cvals,wvals);
			//printf("Finished myconvolve\n");
			//std::cout<<setprecision(2)<<Utemp.extract(10,10);
			for(int xco = K; xco < rows + K ; xco++)
			{
				for(int yco = K; yco < cols + K; yco++)
				{
					U2n[xco-K][yco-K] += float(acopy[pc]) * (Utemp[xco][yco]);
				}
			}

			pc++;
			if(pc > acopy.size()-1)
			{
				printf("pc = %d.. breaking\n",pc);
				break;
			}
			printf("Ended one more iteration of loop\n");
		}
		vnl_matrix< float > sticknew(U2n.rows(),U2n.cols());
		for(int xco = 0; xco < U2n.rows(); xco++)
		{
			for(int yco = 0; yco < U2n.cols(); yco++)
			{
				sticknew[xco][yco] = abs(U2n[xco][yco]);
			}
		}


		float min1 = 1e10f;
		float max1 = -1;
		float rescale_factor = 1.0;
		for(int xco = 0; xco < rows ; xco++)
		{
			for(int yco = 0; yco < cols; yco++)
			{
				min1 = MIN(sticknew[xco][yco],min1);
				max1 = MAX(sticknew[xco][yco],max1);

			}
		}
		printf("min1 = %f max1 = %f\n",min1,max1);
		if(max1 > 255)
		{
			rescale_factor = 255.0/max1;
		}
		Input2DImageType::Pointer out = Input2DImageType::New();
		out->SetRegions(gradmag->GetLargestPossibleRegion());
		out->Allocate();

		typedef itk::ImageRegionIteratorWithIndex<Input2DImageType> IteratorType2;
		IteratorType2 outiter(out,out->GetLargestPossibleRegion());
		outiter.GoToBegin();
		factors[sliceco] = rescale_factor;
		for(;!outiter.IsAtEnd();++outiter)
		{
			Input2DImageType::IndexType index = outiter.GetIndex();
			outiter.Set(MAX(MIN((sticknew[index[1]][index[0]]*rescale_factor),255),0));
		}


		copyslice<unsigned char>(out,realout,sliceco);

	}


	float min_rescale_factor = 999999;
	for(int counter = 0; counter < factors.size(); counter++)
	{
		min_rescale_factor = MIN(min_rescale_factor,factors[counter]);
	}


	typedef itk::ImageRegionIterator< InputImageType > InputIteratorType;

	for(int counter = 0; counter < factors.size(); counter++)
	{
		InputImageType::RegionType region = realout->GetLargestPossibleRegion();
		InputImageType::IndexType index; index.Fill(0);
		index[2] = counter;
		InputImageType::SizeType size  = region.GetSize();
		size[2] = 1;
		region.SetSize(size);
		region.SetIndex(index);
		InputIteratorType initer(realout,region);
		float loss = min_rescale_factor/factors[counter];
		for(initer.GoToBegin(); !initer.IsAtEnd(); ++initer)
		{
			initer.Set(initer.Get()*loss);
		}
	}
	argv[1][strlen(argv[1])-4] = 0;
	char buffer[1024];
	sprintf(buffer,"%s_TV2D.tif",argv[1]);
	writeImage<InputImageType>(realout,buffer);
	time(&t2);
	std::cout<<"tensor_voting_2d preprocessing takes "<< t2-t1 << " seconds (time_t calculation)\n";
	//ck1 = clock();  cout<<"tensor_voting_2d preprocessing takes "<<double(ck1-ck0)/CLOCKS_PER_SEC<<" seconds"<<endl;  ck0 = ck1;
	return 0;
}


