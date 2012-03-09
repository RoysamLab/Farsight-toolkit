#include "MultipleNeuronTracer.h"
#include "time.h"

#ifdef _OPENMP
#include "omp.h"
#endif

MultipleNeuronTracer::MultipleNeuronTracer()
{
}

MultipleNeuronTracer::~MultipleNeuronTracer()
{
}

void MultipleNeuronTracer::LoadCurvImage(std::string fname, unsigned int pad) 
{
	std::cout << "Reading input file "<< fname << std::endl;
	ReaderType::GlobalWarningDisplayOff();
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(fname);
	ImageType3D::Pointer image = reader->GetOutput();
	image->Update();

	std::cout << "Entering LoadCurvImage" << std::endl;
	LoadCurvImage_1(image, pad);
}

void MultipleNeuronTracer::LoadCurvImage_1(ImageType3D::Pointer &image, unsigned int pad)  
{
	ImageType3D::Pointer CurvImage = image;
	padz = pad;

	RescalerType::Pointer rescaler = RescalerType::New();
	rescaler->SetOutputMinimum(0.0);
	rescaler->SetOutputMaximum(1.0);
	rescaler->SetInput(CurvImage);

	//Median filter
	std::cout << "Running Median Filter" << std::endl;
	MedianFilterType::Pointer medfilt = MedianFilterType::New();
	medfilt->SetNumberOfThreads(16);
	medfilt->SetInput(rescaler->GetOutput());
	ImageType3D::SizeType rad = { {1, 1, 1} };
	medfilt->SetRadius(rad);
	medfilt->Update();
	CurvImage = medfilt->GetOutput();

	//pad z slices
	std::cout << "pad z slices" << std::endl;
	itk::Size<3> isz = CurvImage->GetBufferedRegion().GetSize();
	itk::Size<3> osz = isz;
	osz[2] += 2*padz;
	itk::Index<3> indx, ondx;
	
	PaddedCurvImage = ImageType3D::New();
	PaddedCurvImage->SetRegions(osz);
	PaddedCurvImage->Allocate();
	PaddedCurvImage->SetSpacing(CurvImage->GetSpacing());
	
	for(ondx[2] = 0; ondx[2] < osz[2]; ++ondx[2]) 
	{
		indx[2] = (ondx[2] < padz) ? 0 : ondx[2] - padz;
		indx[2] = (ondx[2] >= osz[2]-padz) ? isz[2]-1 : indx[2];
		for(ondx[1] = 0; ondx[1] < osz[1]; ++ondx[1]) 
		{
			indx[1] = ondx[1];
			for(ondx[0] = 0; ondx[0] < osz[0]; ++ondx[0]) 
			{
				indx[0] = ondx[0];
				PaddedCurvImage->SetPixel(ondx, CurvImage->GetPixel(indx));
			}
		}
	}

	std::cout << "Input file size (after zero padding) is " << PaddedCurvImage->GetBufferedRegion().GetSize() << std::endl;
	size = PaddedCurvImage->GetBufferedRegion().GetSize();
	//CurvImage->Delete();
}

///////////////////////////////////////////////////////////////////////
void MultipleNeuronTracer::ReadStartPoints(std::string fname, unsigned int pad) 
{
	padz = pad;

	std::string temp, num;
	std::ifstream infile;
	infile.open(fname.c_str());
	if(!infile.good())
	  {
	  std::cout << "Error reading seed points" << std::endl;
	  exit(1);
	  }
	size_t x1, x2;
	std::cout << "Reading start points " << std::endl;

	while(!infile.eof()) 
	{
		std::getline(infile,temp);
		if (temp.length() < 1)
			continue;
		
		std::cout<<temp; // Prints our STRING.
		x1 = temp.find_first_of("0123456789.");
		x2 = temp.find_first_not_of("0123456789.",x1);
		if ((x2 - x1) > 10)
			continue;
		
		num = temp.substr(x1,x2-x1);
		float x = atof(num.c_str());

		x1 = temp.find_first_of("0123456789.",x2+1);
		x2 = temp.find_first_not_of("0123456789.",x1);
		if ((x2 - x1) > 10)
			continue;
		
		num = temp.substr(x1,x2-x1);
		float y = atof(num.c_str());

		x1 = temp.find_first_of("0123456789.",x2+1);
		x2 = temp.find_first_not_of("0123456789.",x1);
		if (x2 > temp.length())
			x2 = temp.length();
		
		if ((x2 - x1) > 10)
			continue;
		
		num = temp.substr(x1,x2-x1);
		float z = atof(num.c_str());

		itk::Size<3> osz = size;  //original size padz
		osz[2] = osz[2]-padz;
		std::cout <<" after conversion " << x <<" "<< y <<" "<< z << std::endl;

		if ( (x>=0.0) && (y>=0.0) && (z>=0.0) )
		{
			itk::Index<3> n;
			n[0] = long(x + 0.5); 
			if (n[0] >= (unsigned int)osz[0]) 
				n[0] = osz[0]-1;
			n[1] = long(y + 0.5);
			if (n[1] >= (unsigned int)osz[1])
				n[1] = osz[1]-1;
			n[2] = long(z + 0.5);
			if (n[2] >= (unsigned int)osz[2])
				n[2] = osz[2]-1;
			StartPoints.push_back(n);
			std::cout << " is read as " << n << std::endl;
		}
		else
			std::cout << " is discarded (Recommended format XXX YYY ZZZ , Try removing decimal points, add leading zeros in the input text file)" << std::endl;
	}
	infile.close();
}

void MultipleNeuronTracer::ReadStartPoints_1(std::vector< itk::Index<3> > somaCentroids, unsigned int pad) 
{
	padz = pad;

	std::cout << "Reading start points " << std::endl;
	for(int i=0; i<(int)somaCentroids.size(); ++i)
	{
		float x = (float)somaCentroids.at(i)[0];
		float y = (float)somaCentroids.at(i)[1];
		float z = (float)somaCentroids.at(i)[2];
	
		std::cout << x <<" "<< y <<" "<< z << std::endl;
		itk::Size<3> osz = size;  //original size padz
		osz[2] = osz[2]-padz;
		
		if ( (x>=0.0) && (y>=0.0) && (z>=0.0) )
		{
			itk::Index<3> n;
			n[0] = long(x + 0.5); 
			if (n[0] >= (unsigned int)osz[0]) 
				n[0] = osz[0]-1;
			n[1] = long(y + 0.5);
			if (n[1] >= (unsigned int)osz[1])
				n[1] = osz[1]-1;
			n[2] = long(z + 0.5);
			if (n[2] >= (unsigned int)osz[2])
				n[2] = osz[2]-1;
			StartPoints.push_back(n);
			std::cout << " is read as " << n << std::endl;
		}
		else
			std::cout << " is discarded (Recommended format XXX YYY ZZZ , Try removing decimal points, add leading zeros in the input text file)" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////////
void MultipleNeuronTracer::RunTracing(void)
{
	FeatureMain();			//Nice function here that is easy to miss....

	CurrentID = 1;

	//set up the connection image and swc image
	ConnImage = ImageType3D::New();
	ConnImage->SetRegions(PaddedCurvImage->GetBufferedRegion());
	ConnImage->Allocate();
	ConnImage->FillBuffer(MAXVAL);	//MAXVAL is ... needs to be replaced with std::numeric_limit< float >::max()...

	SWCImage = SWCImageType3D::New(); //major memory
	SWCImage->SetRegions(PaddedCurvImage->GetBufferedRegion());
	SWCImage->Allocate();
	SWCImage->FillBuffer(NULL);

	// fill the SWCImage image with start points
	std::vector<IndexType>::iterator startIt;
	int tID = 1;
	
	clock_t fillSWCImage1_start_time = clock();
	for (startIt = StartPoints.begin(); startIt != StartPoints.end(); ++startIt, ++tID)
	{
		itk::Index<3> startIndex = (*startIt);
		startIndex[2] += padz;													//Convert to padded image index
		SWCNode* start_node = new SWCNode(CurrentID++, -1, tID, startIndex);	//This is the seed points SWCNode
		SWCImage->SetPixel(startIndex,start_node);								//Adding all seed points to the SWCImage
		ConnImage->SetPixel(startIndex,0.0f);									//Set the ConnectedImage to 0.0 at all the seed nodes (remember that the Connected image is all initialized with MAXVAL)... 
		SWCNodeContainer.push_back(start_node);									//Fill the SWCNodeContainer with start points
		HeapNode *h = new HeapNode(start_node->ndx, 0.0);						//Heap nodes hold an (index, value) pair
		PQ.push(h);																//Priority Queue contains the seed nodes now...
	}
	std::cout << "fillSWCImage1 took: " << (clock() - fillSWCImage1_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	clock_t fillSWCImage2_start_time = clock();
	
	long eCounter = 0, TotalePoints;
	itk::ImageRegionConstIterator<ImageType3D> Nit(NDXImage, NDXImage->GetBufferedRegion());
	for (Nit.GoToBegin(); !Nit.IsAtEnd(); ++Nit) 
	{
		if (Nit.Get() > 0)	//Vesselness value is greater than 0
		{
			itk::Index<3> endx = Nit.GetIndex();
			SWCNode* s2 = new SWCNode(0, -1, -1*(++eCounter), endx);	//id = 0, parent_id = -1, tree id = -1 * eCounter, index that this vesselness value is greater than 0
			SWCImage->SetPixel(endx,s2);								//Adding all critical points where vesselness value is greater than 0 to the SWC image
		}
	}
	std::cout << "fillSWCImage2 took: " << (clock() - fillSWCImage2_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	TotalePoints = eCounter;
	std::cout<<"eCounter = "<<eCounter<<std::endl;	//eCounter is just number of nodes that are critical points (but not seed points)
	//std::cout << "No of CTs inserted : " <<  TotalePoints << std::endl;

	//Generating some kind of offset neighborhood... this needs to be done with itkNeighborhoodIterator
	itk::Offset<3> x1 = {{-1, 0 ,0}};
	off.push_back( x1 );
	x1[0] = 1;					// x1 = {{1, 0, 0}}
	off.push_back( x1 );
	x1[0] = 0; 
	x1[1] = -1;					// x1 = {{0, -1, 0}}
	off.push_back( x1 );
	x1[1] = 1;					// x1 = {{0, 1, 0}}
	off.push_back( x1 );
	x1[1] = 0; 
	x1[2] = -1;					// x1 = {{0, 0, -1}}
	off.push_back( x1 );
	x1[2] = 1;					// x1 = {{0, 0, 1}}
	off.push_back( x1 );

	std::vector<OffsetType>::iterator oit;
	bool showMessage = false;
	//std::cout << " Heap size: " << PQ.size() << std::endl;
	float KeyValue;
	
	clock_t PQ_popping_start_time = clock();
	
	while(!PQ.empty())	//For each seed node
	{
		//Take the top HeapNode and remove it from the Priority Queue 
		HeapNode *h = PQ.top();
		PQ.pop();

		//Temporarily store the index and value of the node
		itk::Index<3> ndx = h->ndx;
		KeyValue = h->KeyValue;
		delete h;

		//Don't do anything if the heapnode value is larger than the one in the connected image
		if ( KeyValue > ConnImage->GetPixel(ndx) ) 
			continue;
		

		if ((eCounter <= 0) || (KeyValue > CostThreshold)) 
		{
			if (showMessage == true) 
			{
				std::cout << "NOTE: Exiting the search at cost " << CostThreshold << " However, " << (100*eCounter)/TotalePoints << "%% of the image is still not covered, change cost if necessary!!\r"<< std::endl;
				//std::cout << "Cleaning Heap size: " << PQ.size() << std::endl;
				//std::cout<<"keyvalue = "<<KeyValue<<std::endl;
				showMessage = false;
			}
			
			SWCNode* t  = SWCImage->GetPixel(ndx);
			if ( t != NULL) 
			{
				if (t->TreeID < 0) 
				{
					delete t;
				}
			}
			continue;
		}

		SWCNode* s = SWCImage->GetPixel(ndx);
		if (s != NULL) 
		{
			if (s->TreeID < 0) 
			{
				std::vector<IndexType> Chain;
				SWCNode* L = TBack(ndx, Chain);
				if ( L  != NULL ) 
				{
					float costFactor = GetCostLocal( L , ndx);
					std::vector<IndexType>::reverse_iterator cit;
					SWCNode* par = L;
					for (cit = Chain.rbegin(); cit != Chain.rend(); ++cit) 
					{
						SWCNode* t = SWCImage->GetPixel(*cit);
						if (t == NULL) 
						{
							float val = ConnImage->GetPixel(*cit) * costFactor;
							ConnImage->SetPixel((*cit),val);
							SWCNode* s = new SWCNode(CurrentID++, par, L->TreeID, (*cit));
							SWCImage->SetPixel((*cit),s);
							SWCNodeContainer.push_back(s);
							par->children.push_back(s);
							par = s;
							HeapNode *h = new HeapNode((*cit), val);
							PQ.push(h);
						}
						else 
						{
							if (t->TreeID < 0) 
							{
								delete t;
								eCounter--;
								float val = ConnImage->GetPixel(*cit) * costFactor;
								ConnImage->SetPixel((*cit),val);
								SWCNode* s = new SWCNode(CurrentID++, par, L->TreeID, (*cit));
								SWCImage->SetPixel((*cit),s);
								SWCNodeContainer.push_back(s);
								par->children.push_back(s);
								//std::cout<<"SWCImage Node @ " << (*cit) << "(" << s->ID << ") with parent " << par->ID << "  Cost: " << val << "  " << (100*eCounter)/TotalePoints << "% Remaining.\r";// << std::endl;
								par = s;
								HeapNode *h = new HeapNode((*cit), val);
								PQ.push(h);
							}
						}
					}
				} 
			}
		}

		for (oit = off.begin(); oit < off.end(); ++oit) 
		{
			itk::Index<3> ndx2 = ndx + (*oit);
			if ( (ndx2[0] < 2) || (ndx2[1] < 2) || (ndx2[2] < 2) || (ndx2[0] >= unsigned(size[0] - 2)) || (ndx2[1] >= unsigned(size[1] - 2)) || (ndx2[2] >= unsigned(size[2] - 2)) )  
				continue;
			
			if (SWCImage->GetPixel(ndx2) != NULL) 
			{
				if (SWCImage->GetPixel(ndx2)->TreeID > 0) 
				{
					continue;			
				}
			}
			PixelType P = 1/(PaddedCurvImage->GetPixel(ndx2) + 0.001f);  // consider taking inverse here
			PixelType a1, a2, a3;
			ScanNeighbors(a1,a2,a3, ndx2);
			PixelType aa = Update( a1, a2, a3, P );
			if ( ConnImage->GetPixel(ndx2) > aa )  
			{
				ConnImage->SetPixel(ndx2, aa);
				HeapNode *h = new HeapNode(ndx2, aa);
				PQ.push(h);
			}
		}
	}
	
	
	std::cout << "PQ popping took: " << (clock() - PQ_popping_start_time)/(float) CLOCKS_PER_SEC << std::endl;
	
	clock_t Interpolate1_start_time = clock();
	Interpolate(2.0);
	std::cout << "Interpolate1 took: " << (clock() - Interpolate1_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	clock_t Decimate_start_time = clock();
	Decimate();
	std::cout << "Decimate took: " << (clock() - Decimate_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	clock_t Interpolate2_start_time = clock();
	Interpolate(2.0);	
	std::cout << "Interpolate2 took: " << (clock() - Interpolate2_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	clock_t RemoveIntraSomaNodes_start_time = clock();
	RemoveIntraSomaNodes();
	std::cout << "RemoveIntraSomaNodes took: " << (clock() - RemoveIntraSomaNodes_start_time)/(float) CLOCKS_PER_SEC << std::endl;

}

//////////////////////////////////////////////////////////////////////////////////////
//	INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////

void MultipleNeuronTracer::FeatureMain(void)
{
	time_t FeatureMain_start_time = clock();
	std::cout << std::endl<< "Feature detection 3D" << std::endl;
	NDXImage = ImageType3D::New();
	NDXImage->SetRegions(PaddedCurvImage->GetBufferedRegion());
	NDXImage->Allocate();
	NDXImage->FillBuffer(0.0f);
	
	float sigmas[] =  { 2.0f, 2.8284f, 4.0f, 5.6569f, 8.0f, 11.31f };	//LoG scales
	for (unsigned int i = 0; i < 6; ++i)
	{
		std::cout << "Analysis at " << sigmas[i] << std::endl;
		GetFeature( sigmas[i] );			//I guess this is finding all the critical points and throwing their vesselness values into NDXImage
	}

	//itk::RescaleIntensityImageFilter<ImageType3D, ImageType3D>::Pointer rescaler = itk::RescaleIntensityImageFilter<ImageType3D, ImageType3D>::New();
	//rescaler->SetInput(NDXImage);;
	//rescaler->SetOutputMaximum(255);
	//rescaler->SetOutputMinimum(0);

	//itk::CastImageFilter<ImageType3D, CharImageType3D>::Pointer caster = itk::CastImageFilter<ImageType3D, CharImageType3D>::New();
	//caster->SetInput(rescaler->GetOutput());

	//itk::ImageFileWriter<CharImageType3D>::Pointer writer = itk::ImageFileWriter<CharImageType3D>::New();
	//writer->SetInput(caster->GetOutput());
	//writer->SetFileName("C:\\Data\\Darpa\\TEST_FOR_PIPELINE\\23_2100_4200\\seed_points.tif");
	//writer->Update();
	
	/*ImageType3D::Pointer temp = ImageType3D::New();
	temp->SetRegions(PaddedCurvImage->GetBufferedRegion());
	temp->Allocate();

	itk::ImageRegionConstIterator<ImageType3D> Nit(NDXImage, NDXImage->GetBufferedRegion());
	itk::ImageRegionIterator<ImageType3D> tit(temp, temp->GetBufferedRegion());
	itk::ImageRegionConstIterator<ImageType3D> Cit(PaddedCurvImage, PaddedCurvImage->GetBufferedRegion());
	for (Nit.GoToBegin(), Cit.GoToBegin(), tit.GoToBegin(); !Nit.IsAtEnd(); ++Nit, ++tit, ++Cit)
	{
		if (Nit.Get() > 0) 
		{
			tit.Set(1.0);		
		}
	}*/
}

void MultipleNeuronTracer::GetFeature( float sigma ) 
{
	clock_t LoG_start_time = clock();
	typedef itk::LaplacianRecursiveGaussianImageFilter< ImageType3D , ImageType3D> GFilterType;
	GFilterType::Pointer gauss = GFilterType::New();
	gauss->SetInput( PaddedCurvImage );
	gauss->SetSigma( sigma );
	gauss->SetNormalizeAcrossScale(false);
	//ImageType3D::Pointer smoothedCurvImage = gauss->GetOutput();
	gauss->GetOutput()->Update();
	std::cout << "Laplacian of Gaussian at " << sigma << " took " << (clock() - LoG_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	//itk::RescaleIntensityImageFilter<ImageType3D, ImageType3D>::Pointer rescaler = itk::RescaleIntensityImageFilter<ImageType3D, ImageType3D>::New();
	//rescaler->SetInput(gauss->GetOutput());;
	//rescaler->SetOutputMaximum(255);
	//rescaler->SetOutputMinimum(0);

	//itk::CastImageFilter<ImageType3D, CharImageType3D>::Pointer caster = itk::CastImageFilter<ImageType3D, CharImageType3D>::New();
	//caster->SetInput(rescaler->GetOutput());

	//std::stringstream ss;
	//ss << ceil(sigma);
	//itk::ImageFileWriter<CharImageType3D>::Pointer writer = itk::ImageFileWriter<CharImageType3D>::New();
	//writer->SetInput(caster->GetOutput());
	//writer->SetFileName("C:\\Data\\Darpa\\TEST_FOR_PIPELINE\\23_2100_4200\\LOG_" + ss.str() + ".tif");
	//writer->Update();

	float tot = 0.0f, num = 0.0f;
	itk::ImageRegionIterator<ImageType3D> ittemp(gauss->GetOutput(), gauss->GetOutput()->GetBufferedRegion());
	float gamma = 1.6f;
	float tnorm = vcl_pow(sigma,gamma);

	for(ittemp.GoToBegin(); !ittemp.IsAtEnd(); ++ittemp)
	{
		float q = ittemp.Get()*tnorm;
		ittemp.Set(-1.0f*q);
		tot += q*q;
		num ++;
	}
	//std::cout << "Scale "<< sigma << " had average Energy: " << tot <<std::endl;

	// set the diagonal terms in neighborhood iterator
	itk::Offset<3>
		xp =  {{2 ,  0 ,   0}},
		xn =  {{-2,  0,    0}},
		yp =  {{0,   2,   0}},
		yn =  {{0,  -2,    0}},
		zp =  {{0,   0,    2}},
		zn =  {{0,   0,   -2}};

	itk::Size<3> rad = {{1,1,1}};
	itk::NeighborhoodIterator<ImageType3D> nit(rad , gauss->GetOutput(), gauss->GetOutput()->GetBufferedRegion());
	itk::ImageRegionIterator<ImageType3D> it(gauss->GetOutput(), gauss->GetOutput()->GetBufferedRegion());

	unsigned int
		xy1 =  17, //{ 1 ,   1 ,  0 },
		xy2 =  9,  //{ -1,  -1 ,  0 },
		xy3 =  15, //{ -1,   1 ,  0 },
		xy4 =  11, //{ 1 ,  -1 ,  0 },

		yz1 =  25, //{ 0 ,   1 ,  1 },
		yz2 =  1,  //{ 0 ,  -1 , -1 },
		yz3 =  19, //{ 0 ,  -1 ,  1 },
		yz4 =  7,  //{ 0 ,   1 , -1 },

		xz1 =  23, //{ 1 ,   0 ,  1 },
		xz2 =  3,  //{-1 ,   0 , -1 },
		xz3 =  21, //{-1 ,   0 ,  1 },
		xz4 =  5;  //{ 1 ,   0 , -1 };

	typedef itk::FixedArray< double, 3 > EigenValuesArrayType;
	typedef itk::Matrix< double, 3, 3 > EigenVectorMatrixType;
	typedef itk::SymmetricSecondRankTensor<double,3> TensorType;

	itk::Size<3> sz = PaddedCurvImage->GetBufferedRegion().GetSize();
	sz[0] = sz[0] - 3;
	sz[1] = sz[1] - 3; 
	sz[2] = sz[2] - 3;

	it.GoToBegin();
	nit.GoToBegin();
	itk::Vector<float,3> sp = PaddedCurvImage->GetSpacing();

	long win = long(sigma)/2;
	if (win < 2) 
	{
		win = 2;
	}
	
	long ctCnt = 0;
	while(!nit.IsAtEnd()) 
	{
		itk::Index<3> ndx = it.GetIndex();
		if ( (ndx[0] < 2) || (ndx[1] < 2) || (ndx[2] < 2) ||
			(ndx[0] > (unsigned int)sz[0]) || (ndx[1] > (unsigned int)sz[1]) ||
			(ndx[2] > (unsigned int)sz[2]) )
		{
			++it;
			++nit;
			continue;
		}

		float a1 = 0.0;
		for (unsigned int i=0; i < 13; ++i)
		{
			a1 += vnl_math_max(nit.GetPixel(i), nit.GetPixel(26 - i));
		}
		
		float val = nit.GetPixel(13) ;

		const float thresh1 = 0.005;   // 3% of maximum theshold from Lowe 2004
		const float thresh2 = 0.0003;  // -0.1 percent of range

		if ( ((val - a1/13.0f) > thresh2 ) && ( val > thresh1 ))  
		{
			TensorType h;
			h[0] = gauss->GetOutput()->GetPixel( ndx + xp ) + gauss->GetOutput()->GetPixel( ndx + xn ) - 2*nit.GetPixel( 13 );
			h[3] = gauss->GetOutput()->GetPixel( ndx + yp ) + gauss->GetOutput()->GetPixel( ndx + yn ) - 2*nit.GetPixel( 13 );
			h[5] = gauss->GetOutput()->GetPixel( ndx + zp ) + gauss->GetOutput()->GetPixel( ndx + zn ) - 2*nit.GetPixel( 13 );
			h[1] = nit.GetPixel(xy1) + nit.GetPixel(xy2) - nit.GetPixel(xy3) - nit.GetPixel(xy4);
			h[2] = nit.GetPixel(xz1) + nit.GetPixel(xz2) - nit.GetPixel(xz3) - nit.GetPixel(xz4);
			h[4] = nit.GetPixel(yz1) + nit.GetPixel(yz2) - nit.GetPixel(yz3) - nit.GetPixel(yz4);

			EigenValuesArrayType ev;
			EigenVectorMatrixType em;
			h.ComputeEigenAnalysis (ev, em);

			unsigned int w;
			if (IsPlate(ev, w)) 
			{
				float value = vnl_math_abs(ev[0]) + vnl_math_abs(ev[1]) + vnl_math_abs(ev[2]) - vnl_math_abs(ev[w]);
				if (RegisterIndex(value, ndx, sz, win))	//RegisterIndex returns true if this value is the highest in the neighborhood, otherwise it will return false
				{
					NDXImage->SetPixel(ndx,value);
					ctCnt++;			//CriTical Counter I guess
				}
			}
		}
		++it;
		++nit;
	}
	std::cout <<"Number of CTs at this stage: " << ctCnt <<std::endl;
}

bool MultipleNeuronTracer::IsPlate(const itk::FixedArray<float, 3> &ev, unsigned int &w)  
{
	float L1, L2, L;
	if ( (ev[0] > ev[1]) && (ev[0] > ev[2]) ) 
	{
		w = 0;
		L = ev[0];
		L1 = ev[1]; 
		L2 = ev[2];
		if (ev[1] > ev[2])
		{
			L1 = ev[2]; 
			L2 = ev[1];		
		}
	}

	else if( (ev[1] > ev[0]) && (ev[1] > ev[2]) ) 
	{
		w = 1;
		L = ev[1];
		L1 = ev[0];
		L2 = ev[2];
		if (ev[0] > ev[2]) 
		{
			L1 = ev[2];
			L2 = ev[0];		
		}
	}

	else  
	{
		w = 2;
		L = ev[2];
		L1 = ev[0];
		L2 = ev[1];
		if (ev[0] > ev[1]) 
		{
			L1 = ev[1];
			L2 = ev[0];		
		}
	}

	if /*( (abs(L2)/sqrt(abs(L1*L))<0.25) && (abs(L)+abs(L1)+abs(L2)>0.05) )*/(abs(L2)/sqrt(abs(L1*L))<0.5)//((L - L2) > (L2 - L1) && (L - L2) > vnl_math_abs(L)) 
	{
		return true;
	}
	
	return false;  /// right now this is turned off (Amit)
}


//Searches in some specified window in NDXImage around ndx and see if something larger than value is present
//If there is a value in the neighborhood larger than value then return false, else set the NDXImage at ndx to 0 and return true 
bool MultipleNeuronTracer::RegisterIndex(const float value, itk::Index<3> &ndx, itk::Size<3>& sz, long h = 2) 
{
	itk::Index<3> n;
	bool higherPresent = false;
	for (n[0] = ndx[0]-h; n[0] <= ndx[0]+h; ++n[0]) 
	{
		for (n[1] = ndx[1]-h; n[1] <= ndx[1]+h; ++n[1]) 
		{
			for (n[2] = ndx[2]-h; n[2] <= ndx[2]+h; ++n[2]) 
			{
				if ( (n[0] < 2) || (n[1] < 2) || (n[2] < 2) || (n[0] > (unsigned int)sz[0]) ||
					(n[1] > (unsigned int)sz[1]) || (n[2] > (unsigned int)sz[2]) )
				{
					continue;
				}

				float curval = NDXImage->GetPixel(n);
				if (value > curval) 
				{
					NDXImage->SetPixel(n,0.0f);	//Why do we set this to 0.0? We overwrite with value later anyways...
				}
				else if (value < curval) 
				{
					higherPresent = true;				
				}
			}
		}
	}
	if (higherPresent == true) 
	{
		return false;	
	}
	
	return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////

SWCNode* MultipleNeuronTracer::TBack(itk::Index<3> &ndx, std::vector<IndexType>& Chain)  
{
	SWCNode* Label = NULL;
	itk::Index<3> n;
	itk::Vector<float,3> p, x, d, dold;
	for (int i=0; i<3; i++) 
	{
		p[i] = static_cast<PixelType>(ndx[i]);
		dold[i] = 0.0f;
	}
	bool done = false;
	if (SWCImage->GetPixel(ndx)->TreeID > 0) 
	{
		done = true;
	}
	const float MAXDERV = 10000.0f;

	Chain.push_back(ndx);

	while (done == false) 
	{
		//x
		x = p; 
		x[0]++;
		n.CopyWithRound(x);
		if (n[0] < (unsigned int)size[0])
		{
			d[0] = ConnImage->GetPixel(n);
		}
		else
		{
			d[0] = MAXDERV;
		}

		x = p; 
		x[0]--;
		n.CopyWithRound(x);
		if (n[0] >= 0)    
		{
			d[0] -= ConnImage->GetPixel(n);   
		}
		else 
		{
			d[0] -= MAXDERV; 
		}

		// y
		x = p; 
		x[1]++;
		n.CopyWithRound(x);
		if (n[1] < (unsigned int)size[1]) 
		{
			d[1] = ConnImage->GetPixel(n);
		}
		else
		{
			d[1] = MAXDERV;
		}
		
		x = p; 
		x[1]--;
		n.CopyWithRound(x);
		if (n[1] >= 0)
		{
			d[1] -= ConnImage->GetPixel(n);
		}
		else
		{
			d[1] -= MAXDERV; 
		}

		// z
		x = p; 
		x[2]++;
		n.CopyWithRound(x);
		if (n[2] < (unsigned int)size[2]) 
		{
			d[2] = ConnImage->GetPixel(n); 
		}
		else
		{
			d[2] = MAXDERV;
		}
		
		x = p; 
		x[2]--;
		n.CopyWithRound(x);
		if (n[2] >= 0)
		{
			d[2] -= ConnImage->GetPixel(n);
		}
		else
		{
			d[2] -= MAXDERV;
		}

		double norm2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
		if (norm2>0.001)
		{
			d.Normalize();
			d += dold;
			d.Normalize();
			dold = d;
			d *= 0.5;
			p -= d;
		}

		n.CopyWithRound(p);
		Chain.push_back(n);
		//check termination
		SWCNode *t = SWCImage->GetPixel(n);
		if (t != NULL ) 
		{
			if (t->TreeID > 0) 
			{
				done = true;
				Label = SWCImage->GetPixel(n);
				break;
			}
		}
		if (Chain.size() > 500) 
		{
			//std::cout << "Tree not found for " << ndx << " in 500 steps, exiting!! " << std::endl;
			Chain.clear();
			Label = NULL;
			done = true;
			break;
		}
	}
	return Label;
}

///////////////////////////////////////////////////////////////////////////////////
float MultipleNeuronTracer::GetCost(SWCNode* s, itk::Index<3>& endx ) 
{
	itk::Index<3> base = endx, ndx = s->ndx;
	float cost = 0.0f, angsum = 0.0f, count = 0.01f;
	itk::Vector<float,3> d1, d2 , gd, gd1;
	d1.Filled(0.0);
	gd1.Fill(0.0);
	bool first = true;

	while (count < 500.0f) 
	{
		float d = (ndx[0] - base[0])*(ndx[0] - base[0]) + (ndx[1] - base[1])*(ndx[1] - base[1]) + (ndx[2] - base[2])*(ndx[2] - base[2]) ;
		if ( vcl_sqrt(d) > 6.0f) 
		{
			d2 = d1;
			d1[0] = float(ndx[0] - base[0]);
			d1[1] = float(ndx[1] - base[1]);
			d1[2] = float(ndx[2] - base[2]);
			d1.Normalize();
			if (first == true) 
			{
				first = false;
				gd1 = d1;
			}
			else 
			{
				PixelType w = dot_product(d1.Get_vnl_vector(),d2.Get_vnl_vector());
				if (w < 0.99f) 
				{
					angsum += vcl_acos(vnl_math_abs(w));
				}
				count ++;
			}
			base = ndx;
		}
		s = s->parent;
		if (s == NULL) 
		{
			break;
		}
		ndx = s->ndx;
	}
	gd[0] = float(ndx[0] - endx[0]);
	gd[1] = float(ndx[1] - endx[1]);
	gd[2] = float(ndx[2] - endx[2]);
	gd.Normalize();
	float allowedTurns = 1.0f;
	if (dot_product(gd.Get_vnl_vector(),gd1.Get_vnl_vector() ) >= 0) 
	{
		cost = (angsum/allowedTurns);
		if ( cost > 1.0) 
		{
			cost = 1.0f;		
		}
	}
	else 
	{
		cost = 1.0f;
	}
	
	return cost;
}

/////////////////////////////////////////////////////////////////////////////////////
float MultipleNeuronTracer::GetCostLocal(SWCNode* s, itk::Index<3>& endx ) 
{
	itk::Index<3> base = endx, ndx = s->ndx;
	float cost = 0.0f, count = 0.01f;
	itk::Vector<float,3> d1, d2;
	d2.Filled(0.0);

	d1[0] = float(ndx[0] - base[0]);
	d1[1] = float(ndx[1] - base[1]);
	d1[2] = float(ndx[2] - base[2]);
	d1.Normalize();

	base = ndx;

	while (count < 500.0f) 
	{
		float d = (ndx[0] - base[0])*(ndx[0] - base[0]) + (ndx[1] - base[1])*(ndx[1] - base[1]) + (ndx[2] - base[2])*(ndx[2] - base[2]) ;
		if ( vcl_sqrt(d) > 6.0f) 
		{
			d2[0] = float(ndx[0] - base[0]);
			d2[1] = float(ndx[1] - base[1]);
			d2[2] = float(ndx[2] - base[2]);
			d2.Normalize();

			PixelType w = dot_product(d1.Get_vnl_vector(),d2.Get_vnl_vector());
			if ( w <= 0.0f) 
			{
				cost = 1.0f;
			}
			else if (( w > 0.0f) && (w <= 0.98f)) 
			{
				cost = 1.0 - w;
			}
			else 
			{
				cost = 0.0f;
			}
			break;
		}
		count++;
		s = s->parent;
		if (s == NULL) 
		{
			break;
		}
		ndx = s->ndx;
	}

	return cost;
}

///////////////////////////////////////////////////////////////////////////////////
void MultipleNeuronTracer::ScanNeighbors( PixelType &a1, PixelType &a2, PixelType &a3, itk::Index<3> &ndx) 
{
	a1 = MAXVAL;
	if(ndx[0] > 0)
	{
		a1 = ConnImage->GetPixel(ndx + off.at(0));
	}	
	if (ndx[0] < (unsigned int)size[0]-1) 
	{
		a1 = vnl_math_min(ConnImage->GetPixel(ndx + off.at(1)), a1 );
	}
	
	a2 = MAXVAL;
	if(ndx[1] > 0)  
	{
		a2 = ConnImage->GetPixel(ndx + off.at(2));
	}
	if (ndx[1] < (unsigned int)size[1]-1) 
	{
		a2 = vnl_math_min(ConnImage->GetPixel(ndx + off.at(3)), a2 );
	}
	
	a3 = MAXVAL;
	if(ndx[2] > 0)  
	{
		a3 = ConnImage->GetPixel(ndx + off.at(4));
	}
	if (ndx[2] < (unsigned int)size[2]-1) 
	{
		a3 = vnl_math_min(ConnImage->GetPixel(ndx + off.at(5)), a3 );
	}
}

///////////////////////////////////////////////////////////////////////
PixelType MultipleNeuronTracer::Update( PixelType a1,  PixelType a2,  PixelType a3,  PixelType P )  
{
	if (a1 > a2)  
	{
		PixelType temp = a2;
		a2 = a1;
		a1 = temp;
	}

	if (a2 > a3)  
	{
		PixelType temp = a3;
		a3 = a2;
		a2 = temp;
	}

	if (a1 > a2)  
	{
		PixelType temp = a2;
		a2 = a1;
		a1 = temp;
	}

	PixelType A1 = 0;
	PixelType delta = (a2+a1+a3)*(a2+a1+a3) - 3*(a1*a1 + a2*a2 + a3*a3 - P*P);

	if( delta>=0 ) 
	{
		A1 = ( a2+a1+a3 + vcl_sqrt(delta) ) / 3.0;
	}
	
	if( A1 <= a3 )  
	{
		delta = (a2+a1)*(a2+a1) - 2*(a1*a1 + a2*a2 - P*P);
		A1 = 0;
		if( delta>=0 )  
		{
			A1 = ( a2+a1 + vcl_sqrt(delta) ) / 2.0;
		}
		if( A1 <= a2 ) 
		{
			A1 = a1 + P;
		}
	}
	return A1;
}

///////////////////////////////////////////////////////////////////////////////////

void MultipleNeuronTracer::Decimate() 
{
	//std::cout << "Decimating the tree of size: " << SWCNodeContainer.size() << std::endl;
	std::vector<SWCNode*>::iterator sit;
	for (sit = SWCNodeContainer.begin(); sit != SWCNodeContainer.end(); ++sit) 
	{
		if((*sit)->children.size() >= 2) 
		{
			(*sit)->IsBranch = true;
			(*sit)->IsActive = true;
		}
		else if((*sit)->children.size() == 0) 
		{
			(*sit)->IsLeaf = true;
			(*sit)->IsActive = true;
		}
		else if ((*sit)->parent == NULL) 
		{
			(*sit)->IsActive = true;
		}
		else 
		{
			(*sit)->IsActive = false;
		}
	}

	//std::cout << "Tree labeled: 1" << std::endl;
	for (sit = SWCNodeContainer.begin(); sit != SWCNodeContainer.end(); ++sit) 
	{
		if ((*sit)->IsActive == false) 
		{
			if ((*sit)->parent->IsActive == false)  
			{
				bool chActive = false;
				for (unsigned int i = 0; i < (*sit)->children.size(); ++i) 
				{
					if ((*sit)->children[i]->IsActive == true) 
					{
						chActive = true;					
					}
				}
				if (chActive == false) 
				{
					(*sit)->IsActive = true;				
				}
			}
		}
	}

	const float minOffshootLength = 6;
	//std::cout << "Removing offshoots of length less than " << minOffshootLength  << std::endl;

	for (sit = SWCNodeContainer.begin(); sit != SWCNodeContainer.end(); ++sit) 
	{
		if ((*sit)->IsLeaf == true) 
		{
			//std::cout << "Leaf at" << (*sit)->ndx << " ID:" << (*sit)->ID << " ParentID:" << (*sit)->PID << std::endl;
			SWCNode* par = (*sit)->parent;
			if (par == NULL) 
			{
				continue;
			}
			
			if (par->PID == -1) 
			{
				continue;
			}
			
			itk::Vector<PixelType,3> p1 = (*sit)->pos;
			itk::Vector<PixelType,3> p2 = par->pos;
			itk::Vector<PixelType,3> dp = p1 - p2;
			float d = dp.GetNorm();

			while ( par->IsBranch == false ) 
			{
				p1 = p2;
				par = par->parent;
				if (par == NULL) 
				{
					break;
				}
				p2 = par->pos;
				dp = p1 - p2;
				d += dp.GetNorm();
			}

			if (d < minOffshootLength) 
			{
				SWCNode* n = (*sit);
				while ( n != par ) 
				{
					n->IsActive = false;
					n = n->parent;
				}
				if(par != NULL)
				{
					par->IsBranch = false;
				}
			}
		}
	}

	//std::cout << "Tree labeled: 3" << std::endl;

	std::vector<SWCNode*> NewContainer;
	NewContainer.reserve(SWCNodeContainer.size());

	long newID = 1;
	itk::Array<long> IDLookUp(SWCNodeContainer.size());
	IDLookUp.Fill(0);

	for (unsigned int i=0; i < SWCNodeContainer.size(); ++i) 
	{
		if (SWCNodeContainer[i]->IsActive == true) 
		{
			IDLookUp[i] = newID++;		
		}
	}
	//std::cout << "Lookup generated: " << std::endl;

	for (sit = SWCNodeContainer.begin(); sit != SWCNodeContainer.end(); ++sit) 
	{
		if ((*sit)->IsActive == true) 
		{
			long PID;
			itk::Index<3> ndx;
			for (int i = 0; i < 3; ++i) 
			{
				ndx[i] = long((*sit)->pos[i]);
			}

			SWCNode* par = (*sit)->parent;
			if (par == NULL) 
			{
				PID = -1;
			}
			else 
			{
				while (par->IsActive == false)
				{
					par = par->parent;
					if(par == NULL)
					{
						break;
					}
				}
				if(par == NULL)
				{
					PID = -1;
				}
				
				else
				{
					PID = IDLookUp(par->ID - 1);
					if(PID < 1 || (unsigned int)PID > NewContainer.size())
					{
						continue;
					}
					par = NewContainer[PID-1];
				}
			}

			SWCNode* s = new SWCNode();
			s->ID = IDLookUp((*sit)->ID - 1);
			s->PID = PID;
			s->IsActive = true;
			s->IsBranch = (*sit)->IsBranch;
			s->IsLeaf = (*sit)->IsLeaf;
			s->ndx = ndx;
			s->parent = par;
			s->pos = (*sit)->pos;
			s->TreeID = (*sit)->TreeID;
			if (par != NULL) 
			{
				par->children.push_back(s);
			}		
			NewContainer.push_back(s);
		}
	}
	//std::cout << "NewContainer created: " << std::endl;

	for (sit = SWCNodeContainer.begin(); sit != SWCNodeContainer.end(); ++sit) 
	{
		delete (*sit);
	}

	SWCNodeContainer = NewContainer;
}

///////////////////////////////////////////////////////////////////////////////////
void MultipleNeuronTracer::Interpolate(float sigma) 
{
	//std::cout << "Interpolating the tree: " << std::endl;
	typedef itk::SmoothingRecursiveGaussianImageFilter< ImageType3D , ImageType3D> GFilterType;
	GFilterType::Pointer gauss = GFilterType::New();
	gauss->SetInput( PaddedCurvImage );
	gauss->SetSigma( sigma );
	gauss->SetNormalizeAcrossScale(false);
	//ImageType3D::Pointer smoothedCurvImage = gauss->GetOutput();
	gauss->GetOutput()->Update();

	std::vector<SWCNode*>::iterator sit;
	for (sit = SWCNodeContainer.begin(); sit != SWCNodeContainer.end(); ++sit) 
	{
		float w,x,y,z;
		if (((*sit)->children.size() > 0) && ((*sit)->parent != NULL)) 
		{
			w = vnl_math_max(gauss->GetOutput()->GetPixel((*sit)->ndx), 0.1f);
			x = w * float((*sit)->ndx[0]);
			y = w * float((*sit)->ndx[1]);
			z = w * float((*sit)->ndx[2]);

			if ((*sit)->parent != NULL) 
			{
				float w1 = vnl_math_max(gauss->GetOutput()->GetPixel((*sit)->parent->ndx), 0.1f);
				w += w1;
				x += (w1 * float((*sit)->parent->ndx[0]));
				y += (w1 * float((*sit)->parent->ndx[1]));
				z += (w1 * float((*sit)->parent->ndx[2]));
			}
			for (unsigned int i = 0; i < (*sit)->children.size() ; ++i) 
			{
				float w1 = vnl_math_max(gauss->GetOutput()->GetPixel((*sit)->children[i]->ndx), 0.1f);
				w += w1;
				x += (w1 * float((*sit)->children[i]->ndx[0]));
				y += (w1 * float((*sit)->children[i]->ndx[1]));
				z += (w1 * float((*sit)->children[i]->ndx[2]));
			}
			(*sit)->pos[0] = x / w;
			(*sit)->pos[1] = y / w;
			(*sit)->pos[2] = z / w;
			(*sit)->ndx.CopyWithRound((*sit)->pos);
		}
		else 
		{
			(*sit)->pos[0] = float((*sit)->ndx[0]);
			(*sit)->pos[1] = float((*sit)->ndx[1]);
			(*sit)->pos[2] = float((*sit)->ndx[2]);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////
void MultipleNeuronTracer::LoadSomaImage(std::string somaFileName)
{
	/*typedef itk::ImageFileReader<CharImageType3D> SomaReaderType;
	SomaReaderType::Pointer somaReader = SomaReaderType::New();
	somaReader->SetFileName(somaFileName);
	SomaImage = somaReader->GetOutput();
	somaReader->Update();
	*/

	// Now reading a labeled image for the somas

	typedef itk::ImageFileReader<LabelImageType3D> SomaReaderType;
	SomaReaderType::Pointer somaReader = SomaReaderType::New();
	somaReader->SetFileName(somaFileName);
	SomaImage = somaReader->GetOutput();
	somaReader->Update();
}

void MultipleNeuronTracer::RemoveIntraSomaNodes(void)
{
	std::cout << "Removing nodes that fall inside the somas of the Curvelets Image" << std::endl;

	unsigned int originalSize = SWCNodeContainer.size();
	LabelArrayType somaArray = SomaImage->GetBufferPointer();
	itk::Size<3> im_size = SomaImage->GetBufferedRegion().GetSize();
	int slice_size = im_size[0] * im_size[1];
	int row_size = im_size[0];

	//find the root nodes of each tree
	std::cout << "Finding the root nodes of each tree" << std::endl;
	std::map<long, SWCNode*> treeIDToRootMap;
	std::vector<SWCNode*>::iterator sit;
	for (sit = SWCNodeContainer.begin(); sit != SWCNodeContainer.end(); ++sit)
	{
		//assume that a node with no parent is a centroid
		if( (*sit)->parent == NULL )
		{
			treeIDToRootMap[(*sit)->TreeID] = (*sit);
		}
	}
	
	if(treeIDToRootMap.size() != this->StartPoints.size()){
		std::cout << "Centroids missing!!" << std::endl;
		
		/*std::ofstream cent_out_1("cent1.txt");
		std::ofstream cent_out_2("cent2.txt");
		for(int i = 0; i < this->StartPoints.size(); i++)
			cent_out_1 << this->StartPoints[i][0] << "," << this->StartPoints[i][1] << "," << this->StartPoints[i][2] << std::endl;
		for(int i = 0; i < treeIDToRootMap.size(); i++)
			cent_out_2 << treeIDToRootMap[i]->ndx[0] << "," << treeIDToRootMap[i]->ndx[1] << "," << treeIDToRootMap[i]->ndx[2] << std::endl;
		cent_out_1.close();
		cent_out_2.close();*/
	}
	
	itk::Index<3> dummy_index;
	//Removing nodes
	for (sit = SWCNodeContainer.begin(); sit != SWCNodeContainer.end();)
	{
		//don't check nodes that are outside the extent of the soma image
		if ( !SomaImage->GetLargestPossibleRegion().IsInside( (*sit)->ndx ) )
		{
			++sit;
			continue;
		}

		//don't remove centroid nodes
		if( (*sit)->parent == NULL )
		{
			++sit;
			continue;
		}

		//remove any other node that falls within a soma
		/*if ( SomaImage->GetPixel( (*sit)->ndx ) != 0 )
		{
			delete (*sit);
			sit = SWCNodeContainer.erase(sit);
		}*/
		
		// Removing nodes only lying in the foreground of the current soma
		itk::Index<3> Node_0 = (*sit)->ndx;
		itk::Index<3> Node_1 = treeIDToRootMap[(*sit)->TreeID]->ndx;
		//if ( somaArray[(slice_size * Node_0[2]) + (row_size * Node_0[1]) + Node_0[0]] != 0 )
		if ( SomaImage->GetPixel( (*sit)->ndx ) != 0 )
		{
			//if( somaArray[(slice_size * Node_0[2]) + (row_size * Node_0[1]) + Node_0[0]] == somaArray[(slice_size * Node_1[2]) + (row_size * Node_1[1]) + Node_1[0]])
			if( SomaImage->GetPixel((*sit)->ndx) == SomaImage->GetPixel(treeIDToRootMap[(*sit)->TreeID]->ndx) )
			{
				for(int i = 0; i < this->StartPoints.size(); i++)
				{
					if((*sit)->ndx[0] == this->StartPoints[i][0] && (*sit)->ndx[1] == this->StartPoints[i][1] && (*sit)->ndx[2] == this->StartPoints[i][2])
						std::cout << "Centroid " << (*sit)->ndx[0] << ", " << (*sit)->ndx[1] << ", " << (*sit)->ndx[2] << " deleted!! " <<std::endl;
				}

				delete (*sit);
				sit = SWCNodeContainer.erase(sit);
				//std::cout << "Deleted node. " << std::endl;				
			}	

			else
			{
				SWCNode *parent = (*sit)->parent;
				SWCNode *root = treeIDToRootMap[(*sit)->TreeID];

				if(parent->ndx == root->ndx)
				{
					++sit;
					continue;
				}					
				if(SomaImage->GetPixel(parent->ndx) == SomaImage->GetPixel(root->ndx))
				{
					(*sit)->parent = root;
					(*sit)->PID = root->ID;

					++sit;
				}
				else
				{
					++sit;
					continue;
				}		
			}
		}
		
		//otherwise if its parent lies within a soma reassign it to be a child
		//of the centroid instead.
		else
		{
			SWCNode *parent = (*sit)->parent;
			if ( !SomaImage->GetLargestPossibleRegion().IsInside( parent->ndx ) )
			{
				++sit;
				continue;
			}			

			itk::Index<3> Node_2 = parent->ndx;
			//if( somaArray[(slice_size * Node_2[2]) + (row_size * Node_2[1]) + Node_2[0]] != 0)

			if( SomaImage->GetPixel( parent->ndx ) != 0)
			{					
				if( SomaImage->GetPixel(parent->ndx) == SomaImage->GetPixel(treeIDToRootMap[(*sit)->TreeID]->ndx) )
				{
					(*sit)->parent = treeIDToRootMap[(*sit)->TreeID];
					(*sit)->PID = treeIDToRootMap[(*sit)->TreeID]->ID;
				}					
			}
			++sit;
		}
	}

	size_t newSize = SWCNodeContainer.size();
	std::cout << "Just removed " << originalSize - newSize
		<< " nodes (" << originalSize << " to " << newSize << ")"
		<< std::endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////

void MultipleNeuronTracer::WriteMultipleSWCFiles(std::string fname, unsigned int padz) 
{
	// check number of start points to determine number of files to write, with new filename eachtime
	std::cout << "Total " << SWCNodeContainer.size() << " nodes..." <<std::endl;
	std::vector<SWCNode*>::iterator sit;
	float SCALE = 1.0f;

	for (unsigned int i = 0; i < StartPoints.size(); ++i) 
	{
		std::stringstream ss;
		ss << "_" << i+1 << ".swc";
		std::string fname1 = fname;
		fname1.replace(fname.length()-4,8,ss.str());
		std::cout << "Writing SWCImage file " << fname1 << " \n Tree ID " << i+1 <<std::endl;
		std::ofstream ofile(fname1.c_str());
		//ofile << "#Neuron Tracing Code 3D, RPI" << std::endl;
		//ofile << "#author: AM" << std::endl;
	
		//make the LookUp table
		std::map<long, long> NodeIDToSWCIDMap;
		long ID = 1;
		long rootID = 1;
		for (sit = SWCNodeContainer.begin(); sit != SWCNodeContainer.end(); ++sit) 
		{
			if ((*sit)->TreeID == i+1) 
				NodeIDToSWCIDMap[(*sit)->ID] = ID++;			
		}
		std::cout << ID << " Nodes found  ";

		//create the SWCImage file
		for (sit = SWCNodeContainer.begin(); sit != SWCNodeContainer.end(); ++sit) 
		{
			if ((*sit)->TreeID == i+1) 
			{
				long id = NodeIDToSWCIDMap[(*sit)->ID];
				long pid = -1;
				long type = 3;
				if ((*sit)->PID > 0) 
					pid = NodeIDToSWCIDMap[(*sit)->PID];
				
				if(pid == -1) 
				{
					type = 1;
					rootID = NodeIDToSWCIDMap[(*sit)->ID];
				}

				//hack for when your parent was deleted but you didn't get assigned as
				//a child of the root
				if(pid == 0)
					pid = rootID;
				
				//get radius estimate for this node
				float radius = getRadius((*sit)->pos);

				ofile << id << " " << type << " " << SCALE*(*sit)->pos[0] << " "
					<< SCALE*(*sit)->pos[1] << " " << SCALE*(*sit)->pos[2]-padz
					<< " " <<  " " << radius << " " << pid << std::endl;
			}
		}
		ofile.close();
		std::cout << " file written. " << std::endl;
	}

	for (sit = SWCNodeContainer.begin(); sit != SWCNodeContainer.end(); ++sit) 
		delete (*sit);
	
	std::cout << " done! " << std::endl;
}

float MultipleNeuronTracer::getRadius(itk::Vector<float,3>& pos) 
{
	float r = 2.0f;
	itk::Vector<float,3> m1, m2, m;
	itk::Index<3> ndx;

	for (int iter = 0; iter < 20; ++iter) 
	{
		for( int i = 0; i<3; i++) 
		{
			m1[i] = pos[i] - vnl_math_max(2.0f*r, 5.0f);
			m2[i] = pos[i] + vnl_math_max(2.0f*r, 5.0f);
		}

		std::vector<float> c;
		c.reserve(4*4*int(r*r));
		float i1 = 0.0f, i2 = 0.0f, i1s = 0.0f, i2s = 0.0f;
		for (m[2] = m1[2]; m[2] <= m2[2]; m[2]++) 
		{
			for (m[1] = m1[1]; m[1] <= m2[1]; m[1]++) 
			{
				for (m[0] = m1[0]; m[0] <= m2[0]; m[0]++) 
				{
					ndx.CopyWithRound(m);
					itk::Vector<float,3> mm = pos - m;
					float d = mm.GetNorm();
					if (PaddedCurvImage->GetBufferedRegion().IsInside(ndx)) 
					{
						float val = PaddedCurvImage->GetPixel(ndx);
						if (d < r) 
						{
							i1 += val;
							++i1s;
						}
						else 
						{
							i2 += val;
							++i2s;
						}

						if (vnl_math_abs(d - r) < 0.7f) 
							c.push_back(val);						
					}
				}
			}
		}
		i1 /= i1s;
		i2 /= i2s;
		float dr = 0.0f;
		for (std::vector<float>::iterator it = c.begin(); it < c.end(); ++it) 
			dr += vnl_math_abs((*it) - i1) - vnl_math_abs((*it) - i2);
		
		dr *= 1.0f / float(c.size()); //rate
		dr = vnl_math_max(dr , -1.0f);
		dr = vnl_math_min(dr , 1.0f);
		r -= dr;
		r = vnl_math_max(r , 1.0f) ; 
	}
	return r;
}

///////////////////////////////////////////////////////////////////////////////////
void MultipleNeuronTracer::WriteSWCFile(std::string fname, unsigned int padz) 
{
	std::vector<SWCNode*>::iterator sit;
	std::cout << "Writing SWCImage file " << fname << " with " << SWCNodeContainer.size() << " nodes...";
	std::ofstream ofile(fname.c_str());
	//ofile << "#Neuron Tracing Code 3D, RPI" << std::endl;
	//ofile << "#author: AM" << std::endl;
	for (sit = SWCNodeContainer.begin(); sit != SWCNodeContainer.end(); ++sit) 
	{
		//get radius estimate for this node
		float radius = getRadius((*sit)->pos);
		ofile << (*sit)->ID << " 3 " << (*sit)->pos[0] << " " << (*sit)->pos[1] << " "
			<< (*sit)->pos[2]-padz << " " << radius <<" " << (*sit)->PID << std::endl;
		delete (*sit);
	}
	ofile.close();
	std::cout << " done! " << std::endl;
}


vtkSmartPointer< vtkTable > MultipleNeuronTracer::GetSWCTable(unsigned int padz) 
{
	vtkSmartPointer< vtkTable > SWCTable = vtkSmartPointer< vtkTable >::New();
	SWCTable->Initialize();

	vtkSmartPointer< vtkDoubleArray > column = vtkSmartPointer< vtkDoubleArray >::New();
	column->SetName("A");
	SWCTable->AddColumn(column);
	column = vtkSmartPointer< vtkDoubleArray >::New();
	column->SetName("B");
	SWCTable->AddColumn(column);
	column = vtkSmartPointer< vtkDoubleArray >::New();
	column->SetName("C");
	SWCTable->AddColumn(column);
	column = vtkSmartPointer< vtkDoubleArray >::New();
	column->SetName("D");
	SWCTable->AddColumn(column);
	column = vtkSmartPointer< vtkDoubleArray >::New();
	column->SetName("E");
	SWCTable->AddColumn(column);
	column = vtkSmartPointer< vtkDoubleArray >::New();
	column->SetName("F");
	SWCTable->AddColumn(column);
	column = vtkSmartPointer< vtkDoubleArray >::New();
	column->SetName("G");
	SWCTable->AddColumn(column);

	std::vector<SWCNode*>::iterator sit;
	for (sit = SWCNodeContainer.begin(); sit != SWCNodeContainer.end(); ++sit) 
	{
		//get radius estimate for this node
		float radius = getRadius((*sit)->pos);
		vtkSmartPointer< vtkVariantArray > row = vtkSmartPointer< vtkVariantArray >::New();
		row->InsertNextValue(vtkVariant((*sit)->ID));
		row->InsertNextValue(vtkVariant(3));
		row->InsertNextValue(vtkVariant((*sit)->pos[0]));
		row->InsertNextValue(vtkVariant((*sit)->pos[1]));
		row->InsertNextValue(vtkVariant((*sit)->pos[2]-padz));
		row->InsertNextValue(vtkVariant(radius));
		row->InsertNextValue(vtkVariant((*sit)->PID));
		SWCTable->InsertNextRow(row);
		delete (*sit);
	}

	return SWCTable;	
}

///////////////////////////////////////////////////////////////////////
void MultipleNeuronTracer::GenerateTestImage(void) 
{
	PaddedCurvImage = ImageType3D::New();
	size[0] = 20; 
	size[1] = 20; 
	size[2] = 20;
	PaddedCurvImage->SetRegions(size);
	PaddedCurvImage->Allocate();
	PaddedCurvImage->FillBuffer(0.0);

	itk::Vector<float,3> dir; 
	dir.Fill(1.0f); 
	dir.Normalize();
	itk::Vector<float,3> acc, pos;
	pos.Fill(3.0f);
	itk::Index<3> ndx;
	ndx.CopyWithRound(pos);

	PaddedCurvImage->SetPixel(ndx, 1.0f);

	for (int i=0; i<15; i++) 
	{
		float val = float(rand()%100) / 100.0f;
		for (int j = 0;j<3;j++) 
			acc[j] = (float(rand()%100) / 100.0f) - 0.5f;
		
		dir += acc*0.5;
		dir.Normalize();

		pos += dir;
		ndx.CopyWithRound(pos);
		PaddedCurvImage->SetPixel(ndx,val);
	}

	WriteImage3D(std::string("GeneratedImage.mhd"), PaddedCurvImage);
}

void MultipleNeuronTracer::WriteImage3D(std::string fname, MultipleNeuronTracer::ImageType3D::Pointer image)  
{
	std::cout << "Writing output file "<< fname << std::endl;
	typedef itk::ImageFileWriter<ImageType3D> WriterType;
	WriterType::GlobalWarningDisplayOff();
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(fname);
	writer->SetInput(image);
	writer->Update();
}

///////////////////////////////////////////////////////////////////////////////
void MultipleNeuronTracer::BlackOut(itk::Index<3> &stndx)
{
	for (long z = -3; z <=3 ; ++z) 
	{
		for (long y = -5; y <=5 ; ++y) 
		{
			for (long x = -5; x <=5 ; ++x) 
			{
				itk::Offset<3> off = { {x,y,z} };
				itk::Index<3> n = stndx + off;
				if ( (n[0] < 0) || (n[1] < 0) || (n[2] < 0) ||
					(n[0] >= (unsigned int)size[0]) || (n[1] >= (unsigned int)size[1]) ||
					(n[2] >= (unsigned int)size[2]) )  
				{
						continue;
				}
				PaddedCurvImage->SetPixel(n,1.0f);
				NDXImage->SetPixel(n,0);
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	SWCImage NODE and HEAP NODE
//
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
SWCNode::SWCNode()
{
	this->ID = -1;
	this->PID = -1;
	this->TreeID = -1;
	this->IsLeaf = false;
	this->IsBranch = false;
	this->parent = NULL;
	this->children.reserve(2);
}

SWCNode::SWCNode(long id, long parent_id, long tree_id, itk::Index<3> index)
{
	this->ID = id;
	this->PID = parent_id;
	this->TreeID = tree_id;
	this->ndx = index;
	this->IsLeaf = false;
	this->IsBranch = false;
	this->parent = NULL;
	this->children.reserve(2);
}

SWCNode::SWCNode(long id, SWCNode * parent, long tree_id, itk::Index<3> index)
{
	this->ID = id;
	this->PID = parent->ID;
	this->TreeID = tree_id;
	this->ndx = index;
	this->IsLeaf = false;
	this->IsBranch = false;
	this->parent = parent;
	this->children.reserve(2);
}

HeapNode::HeapNode(itk::Index<3> n1, PixelType d)
{
	ndx = n1;
	KeyValue = d;
}





