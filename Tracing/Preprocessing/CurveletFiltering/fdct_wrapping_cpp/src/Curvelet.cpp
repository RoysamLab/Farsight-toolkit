#include "Curvelet.h"
using namespace std;
using namespace fdct_wrapping_ns;
Curvelet::Curvelet()
{
	this->nbangles_coarse = 8;
	this->nsigmas_coarse = 2.2;
	this->nsigmas_fine = 2.5;
	this->neighb_weight = .5;
	this->tuning_neighb = .6;
	this->sigma_ratio = 0.03;
	this->numt = 16;
	this->tile_size = 4024;
	this->border = 100;
	this->slices = 0;

	//this->InputImage = NULL;
}
Curvelet::~Curvelet()
{
}
void Curvelet::SetSigma(float sigma)
{
	this->sigma_ratio = sigma;
}
InputImageType::Pointer Curvelet::RunOnInputImage(InputImageType::Pointer InputImage)
{
	//InputImage = NewInputImage;
	slices = InputImage->GetLargestPossibleRegion().GetSize()[2];
	InputImageType::Pointer outputim = InputImageType::New();
	outputim->SetRegions(InputImage->GetLargestPossibleRegion());
	outputim->Allocate();
	FloatImageType::Pointer cosim = FloatImageType::New();
	cosim->SetRegions(InputImage->GetLargestPossibleRegion());
	cosim->Allocate();
	FloatImageType::Pointer sinim = FloatImageType::New();
	sinim->SetRegions(InputImage->GetLargestPossibleRegion());
	sinim->Allocate();

	if(outputim->GetBufferPointer() == NULL || cosim->GetBufferPointer() == NULL || sinim->GetBufferPointer() == NULL)
	{
		printf("Couldnt' allocate memory - 3.. going to crash now\n");
	}
	int max_dim = tile_size;


	int xsize = InputImage->GetLargestPossibleRegion().GetSize()[0];
	int ysize = InputImage->GetLargestPossibleRegion().GetSize()[1];

	int kx = 0;int ky = 0;


	kx = xsize /(max_dim-this->border);
	ky = ysize /(max_dim-this->border);

	int remx = xsize % (max_dim-this->border);
	int remy = ysize % (max_dim-this->border);

	if ( remx > 0 )
		kx ++;
	if ( remy > 0 )
		ky ++;

	for(int xco = 0; xco < kx; xco++)
	{
		for(int yco = 0; yco < ky; yco++)
		{

			InputImageType::SizeType imsize = InputImage->GetLargestPossibleRegion().GetSize();
			InputImageType::IndexType index;
			InputImageType::SizeType size;
			InputImageType::RegionType region;

			index.Fill(0);
			size[0] =  MIN((xco)*(max_dim-this->border)+max_dim-1,imsize[0]-1) -  xco * (max_dim-this->border) +1;
			size[1] =  MIN((yco)*(max_dim-this->border)+max_dim-1,imsize[1]-1) -  yco * (max_dim-this->border) +1;
			size[2] = imsize[2];

			InputImageType::Pointer imtile = InputImageType::New();
			region.SetIndex(index);
			region.SetSize(size);
			imtile->SetRegions(region);
			imtile->Allocate();
			if(imtile->GetBufferPointer()==NULL)
				printf("Couldn't allocate memory - 4 .. going to crash now\n");
			InputImageType::RegionType region1;
			index[0] = xco *(max_dim-this->border);
			index[1] = yco *(max_dim-this->border);
			index[2] = 0;
			region1.SetIndex(index);
			region1.SetSize(size);

			typedef itk::ImageRegionIterator<InputImageType> IteratorType;
			IteratorType iter1(InputImage,region1);
			IteratorType iter2(imtile,region);


			//printf("xco = %d yco = %d :\n",xco,yco);
			region1.Print(std::cout);
			region.Print(std::cout);

			iter1.GoToBegin();
			iter2.GoToBegin();
			for(;!iter1.IsAtEnd();++iter1,++iter2)
			{
				iter2.Set(iter1.Get());
			}


			InputImageType::Pointer outputtile = InputImageType::New();
			outputtile->SetRegions(imtile->GetLargestPossibleRegion());
			outputtile->Allocate();
			FloatImageType::Pointer cosimtile = FloatImageType::New();
			cosimtile->SetRegions(imtile->GetLargestPossibleRegion());
			cosimtile->Allocate();
			FloatImageType::Pointer sinimtile = FloatImageType::New();
			sinimtile->SetRegions(imtile->GetLargestPossibleRegion());
			sinimtile->Allocate();
			if(outputtile->GetBufferPointer() == NULL || cosimtile->GetBufferPointer()==NULL || sinimtile->GetBufferPointer() == NULL )
			{
				printf("Couldn't allocate memory - 5 .. going to crash now ..\n");
			}

			{
#pragma omp parallel for shared(cosimtile,imtile,sinimtile,outputtile)  num_threads(numt)
				for(int counter = 0; counter < slices; counter++)
				{
					//printf("Counter = %d\n",counter);
					Input2DImageType::Pointer im2d = getSlice(imtile,counter);
					Input2DImageType::Pointer om2d;
					Float2DImageType::Pointer cosim2d,sinim2d;
					//call single slice 2-d curvelets function
					getCurveletsForOneSlice(im2d,om2d,cosim2d,sinim2d);
					copyslice<InputPixelType>(om2d,outputtile,counter);
					copyslice<float>(cosim2d,cosimtile,counter);
					copyslice<float>(sinim2d,sinimtile,counter);
				}
			}
			
			//printf("copying the tile\n");
			if(xco != 0)
			{
				size[0] = size[0] - border/2;
				index[0] = border/2;
			}
			if(xco != kx-1)
			{
				size[0] = size[0] - border/2;
			}

			if(yco != 0)
			{
				size[1] = size[1] - border/2;
				index[1] = border/2;
			}
			if(yco != ky-1)
			{
				size[1] = size[1] - border/2;
			}
			size[2] = slices;
			index[2] = 0;


			region.SetIndex(index);
			region.SetSize(size);
			
			if(xco!=0)
			{
				index[0] = xco *(max_dim-border)+border/2;
			}
			if(yco!=0)
			{
				index[1] = yco *(max_dim-border)+border/2;
			}
			
			
			index[2] = 0;
			region1.SetSize(size);
			region1.SetIndex(index);

			iter1 = IteratorType(outputim,region1);
			iter2 = IteratorType(outputtile,region);
			typedef itk::ImageRegionIterator<FloatImageType> FIteratorType;
			FIteratorType iter3(cosim,region1);
			FIteratorType iter4(cosimtile,region);
			FIteratorType iter5(sinim,region1);
			FIteratorType iter6(sinimtile, region);

			iter1.GoToBegin();iter2.GoToBegin();
			iter3.GoToBegin();iter4.GoToBegin();
			iter5.GoToBegin();iter6.GoToBegin();

			for(;!iter1.IsAtEnd();++iter1,++iter2,++iter3,++iter4,++iter5,++iter6)
			{
				iter1.Set(iter2.Get());
				iter3.Set(iter4.Get());
				iter5.Set(iter6.Get());
			}
			//printf("Done with copying the tile to full image\n");
		}
	}
	return outputim;
}

template<typename T>
void Curvelet::circshift(T input, int shiftx, int shifty, T &output)
{
	int m,n;
	m = input.m();
	n = input.n();
	if(1)
	{
		//		printf("shiftx %d shifty %d\n",shiftx, shifty);
		//		scanf("%*d");
	}
	if(output.m() != input.m() || output.n() != input.n())
	{
		output.resize(m,n);
	}
	for(int cx = 0; cx < m; cx++)
	{
		for(int cy = 0; cy < n; cy++)
		{
			if(shiftx < -m || shifty < -n)
				printf("Error here shiftx %d, shifty %d\n",shiftx,shifty);
			output((cx+shiftx+5*m)%m,(cy+shifty+5*n)%n) = input(cx,cy);
		}
	}
}
Input2DImageType::Pointer Curvelet::getSlice(InputImageType::Pointer im, int slice)
{
	Input2DImageType::Pointer out = Input2DImageType::New();
	Input2DImageType::SizeType size;
	size[0] = im->GetLargestPossibleRegion().GetSize()[0];
	size[1] = im->GetLargestPossibleRegion().GetSize()[1];
	Input2DImageType::IndexType index;
	index.Fill(0);
	Input2DImageType::RegionType region;
	region.SetSize(size);
	region.SetIndex(index);
	out->SetRegions(region);
	out->Allocate();
	if(out->GetBufferPointer()==NULL)
		printf("Could not allocate memory -1 ... I'm going to crash any moment now.. \n");
	memcpy(out->GetBufferPointer(),im->GetBufferPointer()+slice*size[0]*size[1],size[0]*size[1]*sizeof(unsigned char));
	return out;
}

template<typename PixelType>
void Curvelet::copyslice(typename itk::Image<PixelType,2>::Pointer im1, typename itk::Image<PixelType,3>::Pointer im2, int slice)
{
	//printf("Copying slice\n");
	PixelType* inpointer  = im1->GetBufferPointer();
	PixelType* outpointer = im2->GetBufferPointer();
	typename itk::Image<PixelType,2>::SizeType size = im1->GetLargestPossibleRegion().GetSize();
	int copysize = size[0]*size[1];
	memcpy(outpointer+copysize*slice,inpointer,sizeof(PixelType)*copysize);
}
void Curvelet::getCurveletsForOneSlice(Input2DImageType::Pointer im,Input2DImageType::Pointer &om,Float2DImageType::Pointer &cosim, Float2DImageType::Pointer &sinim)
{
	Input2DImageType::SizeType imsize = im->GetLargestPossibleRegion().GetSize();

	int m=1,n=1;
	while(m<imsize[1])
		m <<= 1;
	while(n < imsize[0])
		n <<= 1;
	//m = int(pow(2,ceil(log(float(imsize[1]))/log(2.0)))+0.5);
	//n = int(pow(2,ceil(log(float(imsize[0]))/log(2.0)))+0.5);
	int nbscales;
	if(m < n )
	{
		nbscales = int(log(float(m))/log(2.0f)+0.5 - 3);	
	}
	else
	{
		nbscales = int(log(float(n))/log(2.0f)+0.5 - 3);
	}
	//nbangles_coarse = 8;
	//nsigmas_coarse = 2.2;
	//nsigmas_fine = 2.5;
	//neighb_weight = 0.5;
	//tuning_neighb = 0.6;
	int is_real = 0;

	float sigma = sigma_ratio*255;
	int finest = 1;
	nshifts = 1;
	ac = 1;



	//compute the fdct wrapping for the F
	cpx cpxtemp;
	CpxNumMat F(m,n);
	cpxtemp = std::complex<double>(sqrt(float(n*m)),0);
	F(m/2,n/2) = cpxtemp;


	//printf("Computing L^2 norms...");
	//fdct_wrapping_
	vector< vector<CpxNumMat> > c;  //vector<int> extra;
	fdct_wrapping(m, n, nbscales, nbangles_coarse, ac, F, c);
	F.resize(0,0);

	//printf("Done\n");
	CpxNumMat input(m,n);
	typedef itk::ImageRegionIteratorWithIndex<Input2DImageType> IterType1;

	IterType1 iter1(im,im->GetLargestPossibleRegion());

	for(iter1.GoToBegin(); !iter1.IsAtEnd(); ++iter1)
	{
		Input2DImageType::IndexType index;
		index = iter1.GetIndex();
		cpxtemp = cpx(iter1.Get(),0);
		input(index[1],index[0])=cpxtemp;
	}


	vector< vector<double> > sx, sy;
	vector< vector<double> > fx, fy;
	vector< vector<int> > nx, ny;
	fdct_wrapping_param(m, n, nbscales, nbangles_coarse, ac, sx, sy, fx, fy, nx, ny);
	/*
	printf("M = %d N= %d nbscales = %d\n",m,n,nbscales);
	//scanf("%*d");
	printf("SX: ");
	for(int co1 = 0; co1 < sx.size(); co1++)
	{
	printf("%d: ",sx[co1].size());
	for(int co2 = 0; co2 < sx[co1].size(); co2++)
	{
	printf("%lf ",sx[co1][co2]);
	}
	printf("\n");
	}

	printf("SY: ");
	for(int co1 = 0; co1 < sy.size(); co1++)
	{
	printf("%d: ",sy[co1].size());
	for(int co2 = 0; co2 < sy[co1].size(); co2++)
	{
	printf("%lf ",sy[co1][co2]);
	}
	printf("\n");
	}

	printf("fX: ");
	for(int co1 = 0; co1 < fx.size(); co1++)
	{
	printf("%d: ",fx[co1].size());
	for(int co2 = 0; co2 < fx[co1].size(); co2++)
	{
	printf("%0.2lf ",fx[co1][co2]);
	}
	printf("\n");
	}

	printf("fY: ");
	for(int co1 = 0; co1 < fy.size(); co1++)
	{
	printf("%d: ",fy[co1].size());
	for(int co2 = 0; co2 < fy[co1].size(); co2++)
	{
	printf("%0.2lf ",fy[co1][co2]);
	}
	printf("\n");
	}

	printf("NX: ");
	for(int co1 = 0; co1 < nx.size(); co1++)
	{
	printf("%d: ",nx[co1].size());
	for(int co2 = 0; co2 < nx[co1].size(); co2++)
	{
	printf("%d ",nx[co1][co2]);
	}
	printf("\n");
	}

	printf("NY: ");
	for(int co1 = 0; co1 < ny.size(); co1++)
	{
	printf("%d: ",ny[co1].size());
	for(int co2 = 0; co2 < ny[co1].size(); co2++)
	{
	printf("%d ",ny[co1][co2]);
	}
	printf("\n");
	}
	*/
	//scanf("%*d");

	vector< vector<double> > E;	
	for(int cx = 0; cx< c.size();cx++)
	{
		vector<double> row;
		for(int cy = 0;cy< c[cx].size(); cy++)
		{
			double sum = 0;
			int size1 = c[cx][cy].m();
			int size2 = c[cx][cy].n();
			for(int counter = 0; counter < size1; counter++)
			{
				for(int counter1 =0; counter1 < size2; counter1++)
				{
					sum = sum + abs(c[cx][cy](counter,counter1))*abs(c[cx][cy](counter,counter1));
				}
			}
			row.push_back(sqrt(sum/size1/size2));
		}
		E.push_back(row);
	}


	//printf("Done creating E\n");

	int ndone = 0;
	CpxNumMat temp_restored(m,n);
	CpxNumMat temp_restored_cos(m,n);
	CpxNumMat temp_restored_sin(m,n);
	//printf("about to enter for loop\n");
	double pi = 2*acos(0.0f);
	for(int c1 = 0; c1<nshifts; c1++)
	{
		for(int c2 = 0; c2<nshifts; c2++)
		{
			//printf("In loop %d\n",c2);
			CpxNumMat shift_img(m,n);

			int xshift = 1;
			int yshift = 1;
			circshift(input, xshift,yshift,shift_img);
			ndone = ndone + 1;

			//printf("Direct transform, shift nr. %d ...\n",ndone);
			vector< vector< CpxNumMat > > out;
			vector< vector< CpxNumMat > > outcos;
			vector< vector< CpxNumMat > > outsin;
			fdct_wrapping(m, n, nbscales, nbangles_coarse, ac, shift_img, out);

			//printf("Thresholding...\n");

			double thresh = nsigmas_coarse * sigma;
			for(int j = 0; j < nbscales; j++)
			{
				vector<CpxNumMat> temp;
				temp.resize(out[j].size());
				outcos.push_back(temp);
				outsin.push_back(temp);

				if(j== nbscales-1)
				{
					thresh = nsigmas_fine * sigma;
				}
				for( int l = 0; l < out[j].size(); l++)
				{
					//printf("In j = %d/%d l = %d/%d\n",j,out.size(),l,out[j].size());
					double thresh_jl = thresh*E[j][l];
					outcos[j][l] = CpxNumMat(out[j][l].m(),out[j][l].n());
					outsin[j][l] = CpxNumMat(out[j][l].m(),out[j][l].n());
					DblNumMat modcjl(out[j][l].m(),out[j][l].n());
					CpxNumMat argcjl(out[j][l].m(),out[j][l].n());
					for(int cx = 0; cx < out[j][l].m(); cx++)
					{
						for(int cy = 0; cy < out[j][l].n(); cy++)
						{
							modcjl(cx,cy) = abs(out[j][l](cx,cy));
							argcjl(cx,cy) = out[j][l](cx,cy)/modcjl(cx,cy);
						}
					}

					double rowstep = m*1.0/nx[j][l];
					double colstep = n*1.0/ny[j][l];


					double factor  = 1.0/sqrt(float(1+4*neighb_weight*tuning_neighb));
					int evenquad = 1 - (int(ceil(float((l+1)*4.0/out[j].size()))+0.5)%2);
					double theta = fmod(pi/4 - pi/out[j].size() - 2*pi/out[j].size()*l + 2*pi,pi);

					DblNumMat shift1(modcjl.m(),modcjl.n());
					DblNumMat shift2(modcjl.m(),modcjl.n());
					DblNumMat shift3(modcjl.m(),modcjl.n());
					DblNumMat shift4(modcjl.m(),modcjl.n());
					if (evenquad)
					{
						double fcolsjl;
						if( j == 0 || (finest == 2 && j == nbscales -1))
						{
							fcolsjl = 1;
						}
						else
						{
							fcolsjl = fy[j][l];
						}
						int rowshift = - int( fx[j][l]/fcolsjl*rowstep+ 0.5);
						//if(rowshift>100 || rowshift < -100)
						//{
						//	printf("rowshift %d fx[j][l] = %lf fcolsjl = %lf %d %d\n",rowshift,fx[j][l],fcolsjl,j,l);scanf("%*d");
						//}
						circshift(modcjl,1,0,shift1);
						circshift(modcjl,-1,0,shift2);
						circshift(modcjl,rowshift,1,shift3);
						circshift(modcjl,-rowshift,-1,shift4);


					}
					else
					{
						double frowsjl;
						if( j == 0 || (finest == 2 && j == nbscales -1))
						{
							frowsjl = 1;
						}
						else
						{
							frowsjl = fx[j][l];
						}
						int colshift = - int( fy[j][l]/frowsjl*colstep+ 0.5);
						//if(colshift > 100 || colshift < -100)
						//{
						//	printf("colshift %d fy[j][l] = %lf frowsjl = %lf %d %d\n",colshift,fy[j][l],frowsjl,j,l);scanf("%*d");
						//}
						circshift(modcjl,0,1,shift1);
						circshift(modcjl,0,-1,shift2);
						circshift(modcjl,1,colshift,shift3);
						circshift(modcjl,-1,-colshift,shift4);
					}
					//printf("Finished if else j = %d l = %d\n",j,l);
					double test;
					for(int cx = 0; cx < modcjl.m(); cx++)
					{
						for(int cy = 0; cy < modcjl.n(); cy++)
						{
							//printf("In cx = %d cy = %d\n",cx,cy);
							test = modcjl(cx,cy)*modcjl(cx,cy);
							test += neighb_weight*(shift1(cx,cy)*shift1(cx,cy) + shift2(cx,cy)*shift2(cx,cy));
							test += neighb_weight*(shift3(cx,cy)*shift3(cx,cy) + shift4(cx,cy)*shift4(cx,cy));
							test = sqrt(test);
							test *= factor;
							modcjl(cx,cy) = modcjl(cx,cy)*(test > thresh_jl);
							out[j][l](cx,cy) = modcjl(cx,cy)*argcjl(cx,cy);
							//printf("hi1\n");
							if(j!=0)
							{
								//printf("hi2\n");
								//printf("outcos[j].size() = %d\n",outcos[j].size());
								//printf("outsin[j].size() = %d\n",outsin[j].size());
								outcos[j][l](cx,cy) = cos(theta*2)*out[j][l](cx,cy);
								outsin[j][l](cx,cy) = sin(theta*2)*out[j][l](cx,cy);
							}
						}
					}
					//printf("finished last for loop\n");
				}
			}

			CpxNumMat temp_restored_t(m,n);
			CpxNumMat temp_restored_cos_t(m,n);
			CpxNumMat temp_restored_sin_t(m,n);

			clear(temp_restored_t);clear(temp_restored_cos_t);clear(temp_restored_sin_t);
			//printf("About to call ifdct wrapping1\n");
			ifdct_wrapping(m,n,nbscales,nbangles_coarse, ac,out,temp_restored_t);
			//printf("About to call ifdct wrapping2\n");
			ifdct_wrapping(m,n,nbscales,nbangles_coarse, ac,outcos,temp_restored_cos_t);
			//printf("About to call ifdct wrapping3\n");
			ifdct_wrapping(m,n,nbscales,nbangles_coarse, ac,outsin,temp_restored_sin_t);

			circshift(temp_restored_t,-xshift,-yshift,temp_restored);
			circshift(temp_restored_cos_t,-xshift,-yshift,temp_restored_cos);
			circshift(temp_restored_sin_t,-xshift,-yshift,temp_restored_sin);

		}
	}


	//printf("About to copy the output to images\n");
	//printf("temp_restored size = %d %d\n",temp_restored.m(),temp_restored.n());
	typedef itk::ImageRegionIteratorWithIndex<Input2DImageType> IterType2;
	typedef itk::ImageRegionIteratorWithIndex<Float2DImageType> IterType3;

	Input2DImageType::IndexType index;
	index.Fill(0);
	Input2DImageType::SizeType size2;
	size2[0] = imsize[0];
	size2[1] = imsize[1];
	Input2DImageType::RegionType region;
	region.SetIndex(index);
	region.SetSize(size2);

	om = Input2DImageType::New();
	cosim = Float2DImageType::New();
	sinim = Float2DImageType::New();

	om->SetRegions(region); om->Allocate();
	cosim->SetRegions(region); cosim->Allocate();
	sinim->SetRegions(region); sinim->Allocate();

	if(om->GetBufferPointer()==NULL || cosim->GetBufferPointer() == NULL || sinim->GetBufferPointer()==NULL)
	{
		printf("Couldn't allocate memory - 2.. going to crash now...\n");
	}
	om->FillBuffer(0);
	cosim->FillBuffer(0);
	sinim->FillBuffer(0);
	IterType2 iter2(om,region);
	IterType3 iter3(cosim,region);
	IterType3 iter4(sinim,region);
	iter2.GoToBegin();
	iter3.GoToBegin();
	iter4.GoToBegin();
	//region.Print(std::cout);
	//printf("Beginning loop..\n");
	for(;!iter2.IsAtEnd();++iter2,++iter3,++iter4)
	{
		index = iter2.GetIndex();
		if(temp_restored(index[1],index[0]).real()>=0)
		{
			if(temp_restored(index[1],index[0]).real() <= 255)
				iter2.Set(temp_restored(index[1],index[0]).real());
			else
				iter2.Set(255);
		}
		double oldcos = temp_restored_cos(index[1],index[0]).real();
		double oldsin = temp_restored_sin(index[1],index[0]).real();
		double sum1 = sqrt(oldcos*oldcos + oldsin*oldsin+0.001);
		oldcos /= sum1;
		oldsin /= sum1;
		double angle = fmod(atan2(oldsin,oldcos)+2*pi,2*pi)/2;
		iter3.Set(cos(angle));
		iter4.Set(sin(angle));
	}


	//printf("Returning..\n");
	////ifdct_wrapping_
	//CpxNumMat y(x); clear(y);
	//ifdct_wrapping(m, n, nbscales, nbangles_coarse, ac, c, y);
	//ck1 = clock();  cout<<"IFDCT_WRAPPING_ takes "<<double(ck1-ck0)/CLOCKS_PER_SEC<<" seconds"<<endl;  ck0 = ck1;

	//CpxNumMat e(m,n);
	//for(int i=0; i<m; i++)
	//	for(int j=0; j<n; j++)
	//		e(i,j) = x(i,j) - y(i,j);
	//cerr<<"accuracy of inversion "<<sqrt(energy(e)/(m*n))<<endl;
}