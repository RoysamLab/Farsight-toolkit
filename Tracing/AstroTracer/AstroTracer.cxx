#include "AstroTracer.h"
#include "time.h"



#ifdef _OPENMP
#include "omp.h"
#endif

AstroTracer::AstroTracer()
{
}


void AstroTracer::LoadParameters(const char* parametersFileName)
{
	std::map<std::string, std::string> opts;  
	this->optionsCreate(parametersFileName, opts);
	
	std::map<std::string,std::string>::iterator mi;

	mi = opts.find("-intensity_threshold"); 
	if(mi!=opts.end())
	{ std::istringstream ss((*mi).second); ss>>this->intensity_threshold; 
	}
	else
	{ this->intensity_threshold = 0.005; printf("Chose intensity_threshold = 0.005 as default\n");}

	mi = opts.find("-contrast_threshold");
	if(mi!=opts.end())
	{ std::istringstream ss((*mi).second); ss>>this->contrast_threshold; }
	else
	{	  this->contrast_threshold = 0.0003; printf("Chose contrast_threshold = 0.0003 as default\n"); }

	mi = opts.find("-cost_threshold"); 
	if(mi!=opts.end())
	{ std::istringstream ss((*mi).second); ss>>this->cost_threshold; }
	else
	{ this->cost_threshold = 700; printf("Chose cost_threshold = 700 as default\n");}


	mi = opts.find("-offshoot"); 
	if(mi!=opts.end())
	{ std::istringstream ss((*mi).second); ss>>this->offshoot; }
	else
	{ this->offshoot = 10; printf("Chose offshoot = 10 as default\n"); }


	std::cout<<"intensity_threshold="<<this->intensity_threshold<<std::endl;
	std::cout<<"contrast_threshold="<<this->contrast_threshold<<std::endl;
	std::cout<<"cost_threshold="<<this->cost_threshold<<std::endl;
	std::cout<<"offshoot="<<this->offshoot<<std::endl;

}
void AstroTracer::LoadCurvImage(std::string fname, unsigned int pad) 
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

void AstroTracer::LoadCurvImage_1(ImageType3D::Pointer &image, unsigned int pad)  
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
// 	medfilt->SetNumberOfThreads(16);
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
		indx[2] = (ondx[2] <padz) ? 0 : ondx[2] - padz;
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

void AstroTracer::ReadStartPoints(std::string fname, unsigned int pad) 
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
// 			std::cout << " is read as " << n << std::endl;
		}
		else
			std::cout << " is discarded (Recommended format XXX YYY ZZZ , Try removing decimal points, add leading zeros in the input text file)" << std::endl;
	}
	infile.close();
}

void AstroTracer::ReadStartPoints_1(std::vector< itk::Index<3> > somaCentroids, unsigned int pad) 
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

void AstroTracer::RunTracing(void)
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
		SWCNodeContainer.push_back(start_node);									//Fill the _SWCNodeContainer with start points
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


		if ((eCounter <= 0) || (KeyValue > CostThreshold) ) 
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
					// By influencing this cost, atleast one "bad" node will be added to the tree.
					// Consider placing a hard cost-based threshold to avoid crazy nodes being added
					//float costFactor = GetCostLocal( L , ndx);
					float costFactor = GetCostLocal2(L , ndx);

					//if(costFactor > 1.5)
					//	continue;
					
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

void AstroTracer::SetNScales(int nScales){
	
	this->nScales = nScales;
	std::cout << "Number of scales for LoG set to: " << this->nScales << std::endl;
}

void AstroTracer::SetScaleRange(int start_scale, int end_scale){
	
	this->startScale = start_scale;
	this->endScale = end_scale;
	this->nScales = end_scale - start_scale + 1;
	//std::cout << "Scales = " << this->nScales << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////
//	INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////

void AstroTracer::FeatureMainExternal(void)
{
	time_t FeatureMain_start_time = clock();
	std::cout << std::endl<< "Feature detection 3D" << std::endl;
	
	NDXImage = ImageType3D::New();
	NDXImage->SetRegions(PaddedCurvImage->GetBufferedRegion());
	NDXImage->Allocate();
	NDXImage->FillBuffer(0.0f);

	NDXImage2 = ImageType3D::New();///////////////////
	NDXImage2->SetRegions(PaddedCurvImage->GetBufferedRegion());//////////////////////
	NDXImage2->Allocate();/////////////////////////
	NDXImage2->FillBuffer(0.0f);///////////////////////////
		
	this->LoGScaleImage = ImageType3D::New();
	this->LoGScaleImage->SetRegions(this->PaddedCurvImage->GetBufferedRegion());
	this->LoGScaleImage->Allocate();
	this->LoGScaleImage->FillBuffer(0.0f);

	std::cout << "Padded image size: " << this->PaddedCurvImage->GetBufferedRegion().GetSize() << std::endl;

	
	// These are: 2, root2, 4, root4 and so on...	
	float sigmas[] =  { 2.0f, 2.8284f, 4.0f, 5.6569f, 8.0f, 11.31f };

	if(this->nScales == 0)
		this->nScales = 3; // DEFAULT
	if(this->nScales > 6)
		this->nScales = 6; // MAX

	if(this->startScale == 0)
		this->startScale = 1;
	if(this->endScale == 0)
		this->endScale = 3;

	std::vector<float> sigma_vec;
	
	//std::cout << "Sigma vec size: " << sigma_vec.size() << std::endl;
	for(int i = this->startScale-1; i < this->startScale-1 + this->nScales; i++)
		sigma_vec.push_back(sigmas[i]);
	
	std::cout << "Sigma vec size: " << sigma_vec.size() << std::endl;
	
	for(int i = 0; i < sigma_vec.size(); i++){
		std::cout << "Analysis at " << sigma_vec[i] << std::endl;
		GetFeatureExternal(sigma_vec[i], i+1);
	}
	
	std::cout << "Done with GetFeature. " << std::endl;
}

void AstroTracer::OptimizeCoverage(void){

	std::cout << std::endl<< "Optimizing feature coverage for the image." << std::endl;
	
	// Optimizing coverage at single scale	
	float sigma = 5.6569f;
	float sigma_min = sigma;
	float sigma_max = sigma;
	int sigma_intervals = 1;

	StatisticsFilterType::Pointer stats_filter = StatisticsFilterType::New();
	stats_filter->SetInput(this->PaddedCurvImage);
	stats_filter->Update();
	double img_max_val = stats_filter->GetMaximum();
	double img_mean_val = stats_filter->GetMean();

	/*ObjectnessMeasures obj_measures(sigma_min, sigma_max, sigma_intervals, 1);
	obj_measures.alpha = 0.5*img_max_val;
	obj_measures.beta = 0.5*img_max_val;
	obj_measures.gamma = 0.25*img_max_val;

	//std::cout << "Max image value: " << img_max_val << " Mean img value: " << img_mean_val << std::endl;

	//this->ComputeObjectnessImage(obj_measures);

	StatisticsFilterType::Pointer stats_filter2 = StatisticsFilterType::New();
	stats_filter2->SetInput(this->ObjectnessImage);
	stats_filter2->Update();
	double max_objectness = stats_filter2->GetMaximum();

	MultiplyImageFilter::Pointer img_multiplier = MultiplyImageFilter::New();
	img_multiplier->SetInput(this->ObjectnessImage);
	img_multiplier->SetConstant2(1.0/(max_objectness));
	img_multiplier->Update();

	ImageType3D::Pointer normalized_obj_img = img_multiplier->GetOutput();

	StatisticsFilterType::Pointer stats_filter3 = StatisticsFilterType::New();
	stats_filter3->SetInput(normalized_obj_img);
	stats_filter3->Update();
	double mean_objectness = stats_filter3->GetMean();
	*/

	//std::cout << "Max objectness value: " << max_objectness << std::endl << "Mean objectness value: " << mean_objectness << std::endl;

	/*typedef itk::ImageFileWriter<ImageType3D> ImageWriterType;
	ImageWriterType::Pointer image_writer = ImageWriterType::New();
	image_writer->SetFileName("C:\\Prathamesh\\Astrocytes\\CoverageExp\\VesselnessImage.mhd");
	image_writer->SetInput(this->ObjectnessImage);
	image_writer->Update();
	*/


	LoGFilterType::Pointer gauss = LoGFilterType::New();
	gauss->SetInput( PaddedCurvImage );
	gauss->SetSigma( sigma );
	gauss->SetNormalizeAcrossScale(false);
	gauss->GetOutput()->Update();
	
	float tot = 0.0f, num = 0.0f;
	itk::ImageRegionIterator<ImageType3D> ittemp(gauss->GetOutput(), gauss->GetOutput()->GetBufferedRegion());
	float gamma = 1.6f;
	float tnorm = vcl_pow(sigma,gamma);
	for(ittemp.GoToBegin(); !ittemp.IsAtEnd(); ++ittemp){
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

	itk::Vector<float,3> sp = PaddedCurvImage->GetSpacing();


	// This parameter will affect the density og LoG points detected at each scale.
	float win_scale = 2;

	long win = long(sigma)/2;
	//long win = win_scale * long(sigma);
	if (win <2) 
		win = 2;
	
	int opt_iter = 0, max_opt_iter = 10;
	float max_coverage = 0.002;
	float min_coverage = 0.001;
	float thresh1 = 0.01; //0.05; //0.08; //0.005; //0.03   // 3% of maximum theshold from Lowe 2004
	float thresh2 = 0.0001; //0.0009; //0.015; //0.0003;  //0.001 -0.1 percent of range


	while(opt_iter < max_opt_iter){
		
		NDXImage = ImageType3D::New();
		NDXImage->SetRegions(PaddedCurvImage->GetBufferedRegion());
		NDXImage->Allocate();
		NDXImage->FillBuffer(0.0f);

		it.GoToBegin();
		nit.GoToBegin();
		
		long ctCnt = 0;
		long rejectCtCnt = 0;

		float total_fg_objectness = 0.0;
		float total_fg_intensity = 0.0;
		float total_fg_obj_int = 0.0;
		
		while(!nit.IsAtEnd()){

			itk::Index<3> ndx = it.GetIndex();
			if ( (ndx[0] < 2) || (ndx[1] < 2) || (ndx[2] < 2) ||
				(ndx[0] > (unsigned int)sz[0]) || (ndx[1] > (unsigned int)sz[1]) ||
				(ndx[2] > (unsigned int)sz[2]) ){
				++it;
				++nit;
				continue;
			}

			float a1 = 0.0;
			for (unsigned int i=0; i < 13; ++i)
				a1 += vnl_math_max(nit.GetPixel(i), nit.GetPixel(26 - i));
					
			float val = nit.GetPixel(13);
			
			if ( ((val - a1/13.0f) > thresh2 ) && ( val > thresh1 )){

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
				if (IsSeed(ev, w)){

					float value = vnl_math_abs(ev[0]) + vnl_math_abs(ev[1]) + vnl_math_abs(ev[2]) - vnl_math_abs(ev[w]); // How is this score derived?
					if (RegisterIndex(value, ndx, sz, win)){
						

						NDXImage->SetPixel(ndx, value);
						ctCnt++;

			
						//total_fg_objectness += normalized_obj_img->GetPixel(ndx); //this->ObjectnessImage->GetPixel(ndx);
						total_fg_intensity += this->PaddedCurvImage->GetPixel(ndx);
						//total_fg_obj_int += this->ObjectnessImage->GetPixel(ndx) * this->PaddedCurvImage->GetPixel(ndx);
						//total_fg_obj_int += normalized_obj_img->GetPixel(ndx) * this->PaddedCurvImage->GetPixel(ndx);
					}	
					else
						rejectCtCnt++;
				}
			}
			++it;
			++nit;
		}

		double Mi_coverage = (double)ctCnt / img_mean_val;
		//double Mi_coverage2 = (double)ctCnt / mean_objectness;
		//double Mi_coverage3 = (double)ctCnt / (img_mean_val * mean_objectness);
		
		ImageType3D::SizeType im_size = this->PaddedCurvImage->GetBufferedRegion().GetSize();
		//double total_objectness = mean_objectness * (im_size[0]*im_size[1]*im_size[2]);
		double total_intensity = img_mean_val * (im_size[0]*im_size[1]*im_size[2]);

		//double Mi_coverage4 = total_fg_objectness / total_objectness;
		double Mi_coverage5 = total_fg_intensity / total_intensity;
		
		std::cout << std::endl;
		std::cout << "Number of CTs rejected by RegisterIndex() are: " << rejectCtCnt << std::endl;
		std::cout << "Number of CTs at this stage: " << ctCnt <<std::endl;
		
		//std::cout << "Total foreground objectness: " << total_fg_objectness << std::endl;
		//std::cout << "Total foreground intensity: " << total_fg_intensity << std::endl;
		//std::cout << "Total foreground obj*intensity: " << total_fg_obj_int << std::endl;
		//std::cout << "Mi_coverage: " << Mi_coverage << std::endl;
		//std::cout << "Mi_coverage2: " << Mi_coverage2 << std::endl;
		//std::cout << "Mi_coverage3: " << Mi_coverage3 << std::endl;
		//std::cout << "Mi_coverage4: " << Mi_coverage4 << std::endl;
		std::cout << "Iter: " << opt_iter << " Mi_coverage5: " << Mi_coverage5 << " Intensity threshold: " << thresh1 << ", Contrast threshold: " << thresh2 << std::endl;

		opt_iter++;

		if(Mi_coverage5 > max_coverage){
			thresh1 += 0.02;
			thresh2 += 0.0004;
		}
		else if(Mi_coverage5 < min_coverage){
			thresh1 -= 0.02;
			thresh2 -= 0.0004;
		}
		else
			break;		
	}

	// Setting thresholds based on the optimized coverage
	this->intensity_threshold = thresh1;
	this->contrast_threshold = thresh2;

	if(opt_iter == max_opt_iter)
		std::cout << "Coverage might not be correctly optimized. Manual adjustments might be needed. " << std::endl;
	else
		std::cout << "Done with optimizing coverage. " <<  std::endl;

	
	// Code for writing out separate files for each LoG scale
	/*RescalerType::Pointer rescaler2 = RescalerType::New();
	rescaler2->SetInput(NDXImage);
	rescaler2->SetOutputMaximum( 255 );
	rescaler2->SetOutputMinimum( 0 );
	rescaler2->Update();
	itk::CastImageFilter< ImageType3D, CharImageType3D>::Pointer caster2 = itk::CastImageFilter< ImageType3D, CharImageType3D>::New();
	caster2->SetInput(rescaler2->GetOutput());

	itk::ImageFileWriter< CharImageType3D >::Pointer LoGwriter2 = itk::ImageFileWriter< CharImageType3D >::New();

	LoGwriter2->SetFileName("C:\\Prathamesh\\Astrocytes\\CoverageExp\\LoG_Points_1.tif");	
	LoGwriter2->SetInput(caster2->GetOutput());
	LoGwriter2->Update();
	*/		
}

void AstroTracer::FeatureMain(void)
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
		GetFeature( sigmas[i] );			//I guess this is finding all the critical points and throwing their vesselness values into _NDXImage

	}

}

bool HeapNode::operator ==(const HeapNode& h2){

	return this->ndx[0] == h2.ndx[0] && this->ndx[1] == h2.ndx[1] && this->ndx[2] == h2.ndx[2];
}

void AstroTracer::GetFeatureExternal( float sigma , int scale_index){

	clock_t LoG_start_time = clock();
	LoGFilterType::Pointer gauss = LoGFilterType::New();
	gauss->SetInput( PaddedCurvImage );
	gauss->SetSigma( sigma );
	gauss->SetNormalizeAcrossScale(false);
	//ImageType3D::Pointer smoothedCurvImage = gauss->GetOutput();
	gauss->GetOutput()->Update();

	std::cout << "Laplacian of Gaussian at " << sigma << " took " << (clock() - LoG_start_time)/(float) CLOCKS_PER_SEC << std::endl;

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

	//Store LoG images for later
	//this->LoG_Vector.push_back(gauss->GetOutput());

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


	// This parameter will affect the density og LoG points detected at each scale.
	float win_scale = 2;

	long win = long(sigma)/2;
	//long win = win_scale * long(sigma);
	if (win <2) 
	{
		win = 2;
	}

	ImageType3D::Pointer LoGCurrentScaleImage = ImageType3D::New();
	LoGCurrentScaleImage->SetRegions(this->PaddedCurvImage->GetBufferedRegion());
	LoGCurrentScaleImage->Allocate();
	LoGCurrentScaleImage->FillBuffer(0.0f);
	
	std::vector<HeapNode> current_LoG_points;
	long ctCnt = 0;
	long rejectCtCnt = 0;
	
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
		
		float val = nit.GetPixel(13);


		const float thresh1 = 0.05; //0.08; //0.005; //0.03   // 3% of maximum theshold from Lowe 2004
		const float thresh2 = 0.0009; //0.015; //0.0003;  //0.001 -0.1 percent of range
		
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
			if (IsSeed(ev, w)) 
			{
				// What is RegisterIndex used for?
				float value = vnl_math_abs(ev[0]) + vnl_math_abs(ev[1]) + vnl_math_abs(ev[2]) - vnl_math_abs(ev[w]); // How is this score derived?
				if (RegisterIndex(value, ndx, sz, win)) 
				{

					/*itk::Vector<float, 3> pos;
					pos[0] = ndx[0];
					pos[1] = ndx[1];
					pos[2] = ndx[2];
					float radius = getRadius(pos);
					float radius_th = 2.5;

					//std::cout << radius << std::endl;

					if(radius < radius_th)
						continue;
						*/

					NDXImage->SetPixel(ndx, value);
					ctCnt++;

					this->LoGScaleImage->SetPixel(ndx, scale_index-1);
					LoGCurrentScaleImage->SetPixel(ndx, value);
					current_LoG_points.push_back(HeapNode(ndx, value));	

					/*HeapNode current_node(ndx, value);
					
					if(!this->AllLoGPointsVector.empty()){ 
						
						std::vector<HeapNode>::iterator vec_it = std::find(this->AllLoGPointsVector.begin(), this->AllLoGPointsVector.end(), current_node);
						if(vec_it != this->AllLoGPointsVector.end())
							this->AllLoGPointsVector.erase(vec_it);
					}*/

					this->AllLoGPointsVector.push_back(HeapNode(ndx, value));
				}	
				else
					rejectCtCnt++;
			}
		}
		++it;
		++nit;
	}

	std::cout << " Number of CTs rejected by RegisterIndex() are: " << rejectCtCnt << std::endl;
	std::cout << " Number of CTs at this stage: " << ctCnt <<std::endl;
	
	if(!current_LoG_points.empty()){
		
		this->LoGPointsVector.push_back(current_LoG_points);
		this->LoG_Vector.push_back(gauss->GetOutput());
	
	
		// Code for writing out separate files for each LoG scale
		RescalerType::Pointer rescaler2 = RescalerType::New();
		rescaler2->SetInput(LoGCurrentScaleImage);
		rescaler2->SetOutputMaximum( 255 );
		rescaler2->SetOutputMinimum( 0 );
		rescaler2->Update();
		itk::CastImageFilter< ImageType3D, CharImageType3D>::Pointer caster2 = itk::CastImageFilter< ImageType3D, CharImageType3D>::New();
		caster2->SetInput(rescaler2->GetOutput());

		itk::ImageFileWriter< CharImageType3D >::Pointer LoGwriter2 = itk::ImageFileWriter< CharImageType3D >::New();

		if(scale_index == 1)
			LoGwriter2->SetFileName("C:\\Prathamesh\\Astrocytes\\ControlExp\\Testing\\LoG_Points_1.tif");	
		else if(scale_index == 2)
			LoGwriter2->SetFileName("C:\\Prathamesh\\Astrocytes\\DeviceExp\\Testing\\LoG_Points_2.tif");
		else if(scale_index == 3)
			LoGwriter2->SetFileName("C:\\Prathamesh\\Astrocytes\\DeviceExp\\Testing\\LoG_Points_3.tif");
		else if(scale_index == 4)
			LoGwriter2->SetFileName("C:\\Prathamesh\\Astrocytes\\DeviceExp\\Testing\\LoG_Points_4.tif");
		else if(scale_index == 5)
			LoGwriter2->SetFileName("C:\\Prathamesh\\Astrocytes\\DeviceExp\\Testing\\LoG_Points_5.tif");
		else
			LoGwriter2->SetFileName("C:\\Prathamesh\\Astrocytes\\DeviceExp\\Testing\\LoG_Points_6.tif");
		

		LoGwriter2->SetInput(caster2->GetOutput());
		LoGwriter2->Update();

	}
}

void AstroTracer::GetFeature( float sigma ) 
{
	std::cout<<std::endl<<"Get Feature 1";
	
	clock_t LoG_start_time = clock();
	typedef itk::LaplacianRecursiveGaussianImageFilter< ImageType3D , ImageType3D> GFilterType;
	GFilterType::Pointer gauss = GFilterType::New();
	gauss->SetInput( PaddedCurvImage );
	gauss->SetSigma( sigma );
	gauss->SetNormalizeAcrossScale(false);
	gauss->GetOutput()->Update();
	std::cout << "325Laplacian of Gaussian at " << sigma << " took " << (clock() - LoG_start_time)/(float) CLOCKS_PER_SEC << std::endl;

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
	
	const float thresh1 = intensity_threshold;// 0.01;//0.005;   // 3% of maximum theshold from Lowe 2004
	const float thresh2 = contrast_threshold;//0.003;//0.0003;  // -0.1 percent of range

	long ctCnt = 0;
	int inrt = 0; // niclas testing
	while(!nit.IsAtEnd()) 
	{
// 		std::cout<<inrt++<<"\n "<<std::flush;
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

			if (IsSeed(ev, w)) 
			{
				float value = vnl_math_abs(ev[0]) + vnl_math_abs(ev[1]) + vnl_math_abs(ev[2]) - vnl_math_abs(ev[w]);
				if (RegisterIndex(value, ndx, sz, win))	//RegisterIndex returns true if this value is the highest in the neighborhood, otherwise it will return false
				{
					NDXImage->SetPixel(ndx,value);
					ctCnt++;			//CriTical Counter I guess
					//std::cout<<ctCnt<<" ";
				}
			}
		}
		++it;
		++nit;
	}
	std::cout <<"asdfNumber of CTs at this stage: " << ctCnt <<std::endl<<std::flush;
	//out_seeds.close();
}


bool AstroTracer::IsSeed(const itk::FixedArray<float, 3> &ev, unsigned int &w)  
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
	
	// What is the sequence of the sorted list?
	float lambda1 = L2;
	

	//float ballness = abs(L2) / sqrt(abs(L1 * L)); //<2
	//	
	//if ( ballness > 0.1  && (abs(L) + abs(L1) + abs(L2) > 0.05)) //( abs(L-L1)/abs(L2)>0.6 ) //( (abs(L2)/sqrt(abs(L1*L))<0.5) || ( abs(L-L1)/abs(L2)>0.4) )  //( (abs(L)/sqrt(abs(L2*L1))>1.05) || (abs(L1)/sqrt(abs(L2*L))<0.95) ) //(abs(L1)/sqrt(abs(L2*L))<0.9) //((L - L2) > (L2 - L1) && (L - L2) > vnl_math_abs(L)) 
	//{
	//	return true;
	//}
	
	return true;  /// right now this is turned off (Amit)
}

bool AstroTracer::IsBall(const itk::FixedArray<float, 3>& ev, unsigned int& w, double& ballness){

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
	
	// What is the sequence of the sorted list: L1 < L2 < L
	
	ballness = abs(L1) / sqrt(abs(L2) * abs(L)); //<2
		
	//if ( ballness > 0.1  && (abs(L) + abs(L1) + abs(L2) > 0.05)) //( abs(L-L1)/abs(L2)>0.6 ) //( (abs(L2)/sqrt(abs(L1*L))<0.5) || ( abs(L-L1)/abs(L2)>0.4) )  //( (abs(L)/sqrt(abs(L2*L1))>1.05) || (abs(L1)/sqrt(abs(L2*L))<0.95) ) //(abs(L1)/sqrt(abs(L2*L))<0.9) //((L - L2) > (L2 - L1) && (L - L2) > vnl_math_abs(L)) 
	//{
	//	return true;
	//}
	
	return true;  /// right now this is turned off (Amit)

}

bool AstroTracer::RegisterIndex(const float value, itk::Index<3> &ndx, itk::Size<3>& sz, long h = 2) 
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
					NDXImage->SetPixel(n,0.0f);
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
SWCNode* AstroTracer::TBack(itk::Index<3> &ndx, std::vector<IndexType>& Chain)  
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
	//

	float costFactorLabel = 0.0;
	
	if (Chain.size()!=0 )
		costFactorLabel = GetCostLocalLabel( Label , ndx);

	if (costFactorLabel>=10.0)
		{	
			Chain.clear();
			done=true;
			return NULL;}
	else 
		return Label;
	//


	//return Label;
}

///////////////////////////////////////////////////////////////////////////////////
float AstroTracer::GetCost(SWCNode* s, itk::Index<3>& endx ) 
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

float AstroTracer::GetCostLocal2(SWCNode* s, itk::Index<3>& endx){

	itk::Index<3> new_leaf_ndx = endx, leaf_ndx = s->ndx;
	float cost = 1.0f, cost_low_bound = 1.0f, cost_high_bound = 2.0f;
	float current_curvature = 0.0f;
	int ignore_N_nodes = 6.0, most_N_nodes = 500;
	std::vector<float> curvature_vec;
	itk::Vector<float, 3> pos_diff, pos_diff_last;
	pos_diff.Fill(0.0); pos_diff_last.Fill(0.0);

	
	pos_diff[0] = leaf_ndx[0] - new_leaf_ndx[0];
	pos_diff[1] = leaf_ndx[1] - new_leaf_ndx[1];
	pos_diff[2] = leaf_ndx[2] - new_leaf_ndx[2];
	pos_diff.Normalize();

	for(int i = 0; i < most_N_nodes; i++){

		s = s->parent;
		if(s == NULL) 
			break;
		
		new_leaf_ndx = leaf_ndx;
		leaf_ndx = s->ndx;

		pos_diff_last = pos_diff;
		pos_diff[0] = leaf_ndx[0] - new_leaf_ndx[0];
		pos_diff[1] = leaf_ndx[1] - new_leaf_ndx[1];
		pos_diff[2] = leaf_ndx[2] - new_leaf_ndx[2];
		pos_diff.Normalize();
		
		current_curvature = dot_product(pos_diff.Get_vnl_vector(), pos_diff_last.Get_vnl_vector());
		curvature_vec.push_back(current_curvature);
	}

	if(curvature_vec.empty())
		return 0; //cost_low_bound;
	if(curvature_vec.size() < ignore_N_nodes)
		return 0; //cost_low_bound;
	
	float mean_curvature = std::accumulate(curvature_vec.begin(), curvature_vec.end(), 0.0);
	mean_curvature = mean_curvature / (float)curvature_vec.size();
	
	float std_curvature = 0.0;
	for(int i = 0; i < curvature_vec.size(); i++)
		std_curvature = std_curvature + std::pow(curvature_vec[i] - mean_curvature, 2);
	std_curvature = std::sqrt(std_curvature / (float)curvature_vec.size());
	
	cost = std_curvature + cost_low_bound; // Minimum cost is cost_low_bound when std is zero

	if(cost > cost_high_bound)
		cost = cost_high_bound;
			
	return std_curvature;
}

/////////////////////////////////////////////////////////////////////////////////////
float AstroTracer::GetCostLocal(SWCNode* s, itk::Index<3>& endx ) 
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

float AstroTracer::GetCostLocalLabel(SWCNode* s, itk::Index<3>& endx ) //this functions is used in TBack to prevent nodes to join trees, in cases of abrupt directional changes 
{
	itk::Index<3> base = endx, ndx = s->ndx;
	float cost = 0.0f, count = 0.01f, local_count=0.01f;
	itk::Vector<float,3> d1, d2;
	d2.Filled(0.0);

	d1[0] = float(ndx[0] - base[0]); //  leaf-current
	d1[1] = float(ndx[1] - base[1]);
	d1[2] = float(ndx[2] - base[2]);
	d1.Normalize();

	base = ndx;

	float local_abrupt=0.0;//flag for marking if there is an abrupt directional change as we cross the chain towards close ancestors

	while (count < 500.0f) //500 (5)
	{
		float d = (ndx[0] - base[0])*(ndx[0] - base[0]) + (ndx[1] - base[1])*(ndx[1] - base[1]) + (ndx[2] - base[2])*(ndx[2] - base[2]) ;
		if ( vcl_sqrt(d) > 6.0f) //6.0 (0)
		{
			d2[0] = float(ndx[0] - base[0]); //  ancestor-leaf
			d2[1] = float(ndx[1] - base[1]);
			d2[2] = float(ndx[2] - base[2]);
			d2.Normalize();

			PixelType w = dot_product(d1.Get_vnl_vector(),d2.Get_vnl_vector());

			if ( w <= 0.4f) //0.0 (0.2)
			{
				cost = 1.0f; //1.0 (10.0)
				local_abrupt=1.0;
			}
			else if (( w > 0.4f) && (w <= 0.98f))//0.0 && 0.98  (0.2-0.98)
			{
				cost = 1.0 - w;//-
			}
			else 
			{
				cost = 0.0f;//0
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

	if (local_abrupt==1.0)
		cost=10.0;

	return cost;


}

///////////////////////////////////////////////////////////////////////////////////
void AstroTracer::ScanNeighbors( PixelType &a1, PixelType &a2, PixelType &a3, itk::Index<3> &ndx) 
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
PixelType AstroTracer::Update( PixelType a1,  PixelType a2,  PixelType a3,  PixelType P )  
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

void AstroTracer::Decimate() 
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

	const float minOffshootLength = offshoot;
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
void AstroTracer::Interpolate(float sigma) 
{
	//std::cout << "Interpolating the tree: " << std::endl;
	typedef itk::SmoothingRecursiveGaussianImageFilter< ImageType3D , ImageType3D> LoGFilterType;
	LoGFilterType::Pointer gauss = LoGFilterType::New();
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
void AstroTracer::LoadSomaImage(std::string somaFileName)
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

void AstroTracer::RemoveIntraSomaNodes(void)
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
				for(int i = 0; i < this->StartPoints.size(); i++){
					if((*sit)->ndx[0] == this->StartPoints[i][0] && (*sit)->ndx[1] == this->StartPoints[i][1] && (*sit)->ndx[2] == this->StartPoints[i][2])
						std::cout << "Centroid " << (*sit)->ndx[0] << ", " << (*sit)->ndx[1] << ", " << (*sit)->ndx[2] << " deleted!! " <<std::endl;
				}
				
				delete (*sit);
				sit = SWCNodeContainer.erase(sit);
				//std::cout << "Deleted node. " << std::endl;
			}
			else{
				SWCNode *parent = (*sit)->parent;
				SWCNode *root = treeIDToRootMap[(*sit)->TreeID];

				if(parent->ndx == root->ndx){
					++sit;
					continue;
				}

				if(SomaImage->GetPixel(parent->ndx) == SomaImage->GetPixel(root->ndx)){
					(*sit)->parent = root;
					(*sit)->PID = root->ID;

					++sit;
				}
				else{
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

void AstroTracer::WriteMultipleSWCFiles(std::string fname, unsigned int padz) 
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

float AstroTracer::getRadius(itk::Vector<float,3>& pos) 
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

float AstroTracer::getRadiusAndLikelihood(itk::Vector<float,3> &pos, float& likelihood){
	
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

		likelihood = i1 - i2;

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
void AstroTracer::WriteSWCFile(std::string fname, unsigned int padz) 
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

///////////////////////////////////////////////////////////////////////
void AstroTracer::GenerateTestImage(void) 
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

void AstroTracer::WriteImage3D(std::string fname, AstroTracer::ImageType3D::Pointer image)  
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
void AstroTracer::BlackOut(itk::Index<3> &stndx)
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

int AstroTracer::optionsCreate(const char* optfile, std::map<std::string,std::string>& options)
{
	options.clear();
	
	std::ifstream fin(optfile); assert(fin.good());
	std::string name;  fin>>name;
	while(fin.good()) {
		char cont[100];	 fin.getline(cont, 99);
		options[name] = std::string(cont);
		fin>>name;
	}
	fin.close();
	return 0;
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
	ID = -1;
	PID = -1;
	TreeID = -1;
	IsLeaf = false;
	IsBranch = false;
	parent = NULL;
	children.reserve(2);
}

SWCNode::SWCNode(long id, long pid, long tid, itk::Index<3> n)
{
	ID = id;
	PID = pid;
	TreeID = tid;
	ndx = n;
	IsLeaf = false;
	IsBranch = false;
	parent = NULL;
	children.reserve(2);
}

SWCNode::SWCNode(long id, SWCNode * p, long tid, itk::Index<3> n)
{
	ID = id;
	PID = p->ID;
	TreeID = tid;
	ndx = n;
	IsLeaf = false;
	IsBranch = false;
	parent = p;
	children.reserve(2);
}

HeapNode::HeapNode(itk::Index<3> n1, PixelType d)
{
	ndx = n1;
	KeyValue = d;
}

void AstroTracer::CallFeatureMainExternal(void){
	this->FeatureMainExternal();
}

void AstroTracer::ComputeObjectnessImage(ObjectnessMeasures obj_measures){

	//float sigma_min = 2.0f; //0.5f;
	//float sigma_max = 10.0f; //4.0f;
	//int sigma_steps = 5;

	//float alpha = 0.5, beta = 0.5, gamma = 0.25; //5.0;

	//int obj_dim = objectness_type; //1; //0: Blobness, 1: Vesselness, 2: Plateness

	MultiScaleHessianFilterType::Pointer multi_scale_Hessian = MultiScaleHessianFilterType::New();
	multi_scale_Hessian->SetInput(this->PaddedCurvImage);
	multi_scale_Hessian->SetSigmaMin(obj_measures.sigma_min);
	multi_scale_Hessian->SetSigmaMax(obj_measures.sigma_max);
	multi_scale_Hessian->SetNumberOfSigmaSteps(obj_measures.sigma_intervals);

	ObjectnessFilterType::Pointer objectness_filter = ObjectnessFilterType::New();
	objectness_filter->SetScaleObjectnessMeasure(false);
	objectness_filter->SetBrightObject(true);
	objectness_filter->SetAlpha(obj_measures.alpha);
	objectness_filter->SetBeta(obj_measures.beta);
	objectness_filter->SetGamma(obj_measures.gamma);
	objectness_filter->SetObjectDimension(obj_measures.objectness_type);

	//std::cout << obj_measures.alpha << std::endl << obj_measures.beta << std::endl << obj_measures.gamma << std::endl;

	multi_scale_Hessian->Update();
	
	this->ObjectnessImage = multi_scale_Hessian->GetOutput();

	/*typedef itk::ImageFileWriter<ImageType3D> ImageWriterType;
	ImageWriterType::Pointer image_writer = ImageWriterType::New();

	image_writer->SetFileName("C:\\Prathamesh\\Astrocytes\\CoverageExp\\VesselnessImage.mhd");
	image_writer->SetInput(multi_scale_Hessian->GetOutput());
	image_writer->Update();

	ImageWriterType::Pointer image_writer2 = ImageWriterType::New();
	image_writer2->SetFileName("C:\\Prathamesh\\Astrocytes\\CoverageExp\\VesselnessImage_scales.mhd");
	image_writer2->SetInput(multi_scale_Hessian->GetScalesOutput());
	image_writer2->Update();*/
}

void AstroTracer::ComputeFTKObjectnessImage(){

	float sigma = 5.6569f;
	
	// Code for FTK objectness
	LoGFilterType::Pointer gauss = LoGFilterType::New();
	gauss->SetInput( PaddedCurvImage );
	gauss->SetSigma( sigma );
	gauss->SetNormalizeAcrossScale(false);
	gauss->GetOutput()->Update();

	float tot = 0.0f, num = 0.0f;
	itk::ImageRegionIterator<ImageType3D> ittemp(gauss->GetOutput(), gauss->GetOutput()->GetBufferedRegion());
	float gamma = 1.6f;
	float tnorm = vcl_pow(sigma,gamma);

	// Replace by ITK filters, will be faster
	for(ittemp.GoToBegin(); !ittemp.IsAtEnd(); ++ittemp){
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


	// This parameter will affect the density og LoG points detected at each scale.
	float win_scale = 2;

	long win = long(sigma)/2;
	//long win = win_scale * long(sigma);
	if (win <2) 
		win = 2;
		
	const float thresh1 = 0.05; //0.08; //0.005; //0.03   // 3% of maximum theshold from Lowe 2004
	const float thresh2 = 0.0009; //0.015; //0.0003;  //0.001 -0.1 percent of range
	
	
	//Loop over image with the neighborhood for computing objectness
	float alpha = 0.5, beta = 0.5, gamma1 = 0.0025;

	MinMaxCalculatorType::Pointer min_max_calculator = MinMaxCalculatorType::New();
	min_max_calculator->SetImage(this->PaddedCurvImage);
	min_max_calculator->ComputeMaximum();
	double img_max_val = min_max_calculator->GetMaximum();

	ImageType3D::RegionType id_reg;
	ImageType3D::IndexType id_st;
	ImageType3D::SizeType id_sz = PaddedCurvImage->GetBufferedRegion().GetSize();

	id_st[0] = 0;
	id_st[1] = 0;
	id_st[2] = 0;
	
	id_reg.SetSize(id_sz);
	id_reg.SetIndex(id_st);
	
	this->ObjectnessImage = ImageType3D::New();
	this->ObjectnessImage->SetRegions(id_reg);
	this->ObjectnessImage->Allocate();
	this->ObjectnessImage->SetSpacing(PaddedCurvImage->GetSpacing());

	this->ObjectnessImage->FillBuffer(0);
	
	
	while(!nit.IsAtEnd()) 
	{
		itk::Index<3> ndx = it.GetIndex();
		if ( (ndx[0] < 2) || (ndx[1] < 2) || (ndx[2] < 2) ||
			(ndx[0] > (unsigned int)sz[0]) || (ndx[1] > (unsigned int)sz[1]) ||
			(ndx[2] > (unsigned int)sz[2]) ){
			++it;
			++nit;
			continue;
		}
		
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

		ObjectnessMeasures obj_measures(alpha * img_max_val, beta * img_max_val, gamma1 * img_max_val);
		
		this->GetHessianBasedObjectnessMeasures(ev, obj_measures);
		
		this->ObjectnessImage->SetPixel(ndx, obj_measures.vesselness);
		
		++it;
		++nit;
	}
}

ObjectnessMeasures::ObjectnessMeasures(){

	this->alpha = 0.5;
	this->beta = 0.5;
	this->gamma = 0.25;
	
	this->sigma_min = 0.5;
	this->sigma_max = 2.0;
	this->sigma_intervals = 1;
	this->objectness_type = 1;

	this->ballness = 0.0;
	this->plateness = 0.0;
	this->vesselness = 0.0;
	this->noiseness = 0.0;
}

ObjectnessMeasures::ObjectnessMeasures(float alpha, float beta, float gamma){

	this->alpha = alpha;
	this->beta = beta;
	this->gamma = gamma;

	this->sigma_min = 0.5;
	this->sigma_max = 2.0;
	this->sigma_intervals = 1;
	this->objectness_type = 1;

	this->ballness = 0.0;
	this->plateness = 0.0;
	this->vesselness = 0.0;
	this->noiseness = 0.0;
}

ObjectnessMeasures::ObjectnessMeasures(float sigma_min, float sigma_max, float sigma_intervals, int obj_type){

	this->alpha = 0.5;
	this->beta = 0.5;
	this->gamma = 0.25;

	this->sigma_min = sigma_min;
	this->sigma_max = sigma_max;
	this->sigma_intervals = sigma_intervals;
	this->objectness_type = obj_type;

	this->ballness = 0.0;
	this->plateness = 0.0;
	this->vesselness = 0.0;
	this->noiseness = 0.0;
}

void AstroTracer::GetHessianBasedObjectnessMeasures(itk::FixedArray<double, 3>& ev, ObjectnessMeasures& obj_measures){

	bool zeroness = false;
	int lowest_idx = 0;
	itk::FixedArray<double, 3> ev_original;

	ev_original[0] = ev[0];
	ev_original[1] = ev[1];
	ev_original[2] = ev[2];

	ev[0] = std::abs(ev[0]);
	ev[1] = std::abs(ev[1]);
	ev[2] = std::abs(ev[2]);

	float L1, L2, L;
	if ( (ev[0] > ev[1]) && (ev[0] > ev[2]) ) 
	{
		L = ev[0];
		L1 = ev[1]; 
		L2 = ev[2];
		if (ev[1] > ev[2])
		{
			L1 = ev[2]; 
			L2 = ev[1];

			lowest_idx = 2;
		}
		else
			lowest_idx = 1;
	}
	else if( (ev[1] > ev[0]) && (ev[1] > ev[2]) ) 
	{
		L = ev[1];
		L1 = ev[0];
		L2 = ev[2];
		if (ev[0] > ev[2]) 
		{
			L1 = ev[2];
			L2 = ev[0];	

			lowest_idx = 2;
		}
		else
			lowest_idx = 0;
	}
	else  
	{
		L = ev[2];
		L1 = ev[0];
		L2 = ev[1];
		if (ev[0] > ev[1]) 
		{
			L1 = ev[1];
			L2 = ev[0];		

			lowest_idx = 1;
		}
		else
			lowest_idx = 0;
	}

	// What is the sequence of the sorted list: L1 < L2 < L

	if(lowest_idx == 0){
		if(ev_original[1] > 0 || ev_original[2] > 0)
			zeroness = true;
	}
	if(lowest_idx == 1){
		if(ev_original[0] > 0 || ev_original[2] > 0)
			zeroness = true;
	}
	if(lowest_idx == 2){
		if(ev_original[0] > 0 || ev_original[1] > 0)
			zeroness = true;
	}
	
	if(zeroness){
		obj_measures.ballness = 0.0;
		obj_measures.plateness = 0.0;
		obj_measures.vesselness = 0.0;
		obj_measures.noiseness = 0.0;
	}
	else{

		//float alpha = 0.5 * img_max_val; 
		//float beta = 0.5 * img_max_val;
		//float gamma = 0.0025 * img_max_val;

		obj_measures.ballness = L1 / sqrt(L2 * L);
		obj_measures.plateness = L2 / L;
		obj_measures.noiseness = std::sqrt(std::pow(L1, 2) + std::pow(L2, 2) + std::pow(L, 2));

		obj_measures.vesselness = (1 - std::exp(-1 * std::pow(obj_measures.plateness, 2) / (2.0 * std::pow(obj_measures.alpha, 2)))) 
			* std::exp(-1 * std::pow(obj_measures.ballness, 2) / (2.0 * std::pow(obj_measures.beta, 2)))
			* (1 - std::exp(-1 * std::pow(obj_measures.noiseness, 2) / (2.0 * std::pow(obj_measures.gamma, 2))));
	}
}

void AstroTracer::ComputeAstroFeatures(std::string outputFname, std::string IDFname, unsigned int padz, const std::string points_source = "SWC"){

	//Prepare a list of points at which to compute the features
	std::string option1 = "SWC";
	std::string option2 = "LOG";
	std::string option3 = "External";

	std::vector<HeapNode> points_list;
	if(!strcmp(points_source.c_str(), option1.c_str())){
		std::cout << "Computing features around SWC points. " << std::endl;

		std::vector<SWCNode*>::iterator sit;
		for(sit = SWCNodeContainer.begin(); sit != SWCNodeContainer.end(); ++sit){
			itk::Index<3> idx; 
			idx[0] = (*sit)->pos[0]; idx[1] = (*sit)->pos[1]; idx[2] = (*sit)->pos[2];
			points_list.push_back(HeapNode(idx, 0));
		}
	}
	else if(!strcmp(points_source.c_str(), option2.c_str())){
		std::cout << "Computing features around LOG points. " << this->LoGPointsVector.size() << std::endl;

		if(this->AllLoGPointsVector.empty()){
			std::cout << "NO LOG POINTS FOUND. RETURNING. " << std::endl;
			return;
		}

		for(int i = 0; i < this->AllLoGPointsVector.size(); i++){
			//std::vector<HeapNode> currentLoGPoints = this->LoGPointsVector[i];
			//for(int j = 0; j < currentLoGPoints.size(); j++)
			
			points_list.push_back(this->AllLoGPointsVector[i]);
		}
	}
	else
		std::cout << "Incorrect option for calculating the features. Quitting. " << std::endl; 
	
	std::cout << "Computing features around N points. " << points_list.size() << std::endl;

	//Preparing IDImage
	LabelImageType3D::RegionType id_reg;
	LabelImageType3D::IndexType id_st;
	LabelImageType3D::SizeType id_sz = PaddedCurvImage->GetBufferedRegion().GetSize();

	id_st[0] = 0;
	id_st[1] = 0;
	id_st[2] = 0;
	
	id_reg.SetSize(id_sz);
	id_reg.SetIndex(id_st);
	
	IDImage = LabelImageType3D::New();
	IDImage->SetRegions(id_reg);
	IDImage->Allocate();
	IDImage->SetSpacing(PaddedCurvImage->GetSpacing());
	IDImage->FillBuffer(0);


	typedef itk::MinimumMaximumImageCalculator<ImageType3D> MinMaxCalculatorType;
	MinMaxCalculatorType::Pointer min_max_calculator = MinMaxCalculatorType::New();
	min_max_calculator->SetImage(this->PaddedCurvImage);
	min_max_calculator->ComputeMaximum();
	double img_max_val = min_max_calculator->GetMaximum();

	std::vector<double> radius_vec;
	std::vector<double> int_vec, mean_int_vec, var_int_vec, min_int_vec, max_int_vec;

	std::ofstream feature_vector;////////////
	feature_vector.open(outputFname.c_str(), std::ios::out);//////////////////
	//std::cout << "After feature_vector.open"<<std::endl;
	
	unsigned short IDIndex = 0;//ID index
	feature_vector << "ID" << '\t' << "x" << '\t' << "y" << '\t' << "z" << '\t' << "radius" << '\t' << "sph_likelihood" << '\t' << "ball" << '\t' << "plate" << '\t';
	feature_vector << "vessel" << '\t' << "int" << '\t' << "m_int" << '\t' << "v_int" << '\t' << "mx_int" << '\t' << "mn_int" << '\t' << "nucl_dist" << '\t';
	feature_vector << std::endl;

	for(int i = 0; i < points_list.size(); i++){

		//std::cout << " LoG point index: " << IDIndex << std::endl;


		//get radius estimate for this node
		itk::Vector<float, 3> pos;
		pos[0] = points_list[i].ndx[0];
		pos[1] = points_list[i].ndx[1];
		pos[2] = points_list[i].ndx[2];

		float likelihood = 0.0;
		float radius = getRadiusAndLikelihood(pos, likelihood);

		double radiusThresh = 2;
		if(radius < radiusThresh)
			continue;

		ImageType3D::IndexType current_idx;
		current_idx[0] = pos[0];
		current_idx[1] = pos[1];
		current_idx[2] = pos[2]-padz; // PADZ???
		//std::cout << "After current_ide set, LoG_vector_size="<< this->LoG_Vector.size() << std::endl;
		//std::cout << "Current Indices:"<<current_idx[0]<<" "<<current_idx[1]<<" "<<current_idx[2]<<" "<<std::endl;


		////////////////// Code for nearest nuiclei /////////////////////////
		float double_scale_nuclei = 20;

		CharImageType3D::IndexType starting_index_nuclei, end_index_nuclei;
		CharImageType3D::SizeType sub_volume_size_nuclei;
		CharImageType3D::RegionType sub_volume_region_nuclei;
		CharImageType3D::Pointer sub_volume_nuclei;

		starting_index_nuclei[0] = current_idx[0] - double_scale_nuclei; starting_index_nuclei[1] = current_idx[1] - double_scale_nuclei; starting_index_nuclei[2] = current_idx[2] - double_scale_nuclei;
		end_index_nuclei[0] = current_idx[0] + double_scale_nuclei; end_index_nuclei[1] = current_idx[1] + double_scale_nuclei; end_index_nuclei[2] = current_idx[2] + double_scale_nuclei;

		itk::Size<3> sz = PaddedCurvImage->GetBufferedRegion().GetSize();

		//std::cout << "Nuclei: Starting Indices:"<<starting_index_nuclei[0]<<" "<<starting_index_nuclei[1]<<" "<<starting_index_nuclei[2]<<" "<<std::endl;
		//std::cout << "Nuclei: End Indices:"<<end_index_nuclei[0]<<" "<<end_index_nuclei[1]<<" "<<end_index_nuclei[2]<<std::endl;

		if ( (starting_index_nuclei[0] < 0) || (starting_index_nuclei[1] < 0) || (starting_index_nuclei[2] < 0) ||
			(end_index_nuclei[0] > (unsigned int)sz[0]) || (end_index_nuclei[1] > (unsigned int)sz[1]) ||
			(end_index_nuclei[2] > (unsigned int)sz[2]) )
			continue;

		sub_volume_size_nuclei[0] = 2 * double_scale_nuclei; sub_volume_size_nuclei[1] = 2 * double_scale_nuclei; sub_volume_size_nuclei[2] = 2 * double_scale_nuclei;
		
		sub_volume_region_nuclei.SetIndex(starting_index_nuclei);
		sub_volume_region_nuclei.SetSize(sub_volume_size_nuclei);

		typedef itk::BinaryThresholdImageFilter<LabelImageType3D, CharImageType3D> ThresholdFilterType;
		ThresholdFilterType::Pointer threshold_filter = ThresholdFilterType::New();
		threshold_filter->SetLowerThreshold(1);
		threshold_filter->SetInsideValue(255);
		threshold_filter->SetOutsideValue(0);
		threshold_filter->SetInput(this->SomaImage);
		threshold_filter->Update();


		typedef itk::RegionOfInterestImageFilter<CharImageType3D, CharImageType3D> VolumeOfInterestFilterType_nuclei2;
		VolumeOfInterestFilterType_nuclei2::Pointer sub_volume_filter_nuclei = VolumeOfInterestFilterType_nuclei2::New();
		sub_volume_filter_nuclei->SetInput(threshold_filter->GetOutput());
		sub_volume_filter_nuclei->SetRegionOfInterest(sub_volume_region_nuclei);
		sub_volume_filter_nuclei->Update();
		sub_volume_nuclei = sub_volume_filter_nuclei->GetOutput();

		
		/*typedef itk::ExtractImageFilter<LabelImageType3D, CharImageType3D>  ExtractSubVolumeFilterType_nuclei;
		ExtractSubVolumeFilterType_nuclei::Pointer sub_volume_filter_nuclei = ExtractSubVolumeFilterType_nuclei::New();
		sub_volume_filter_nuclei->SetExtractionRegion(sub_volume_region_nuclei);
		sub_volume_filter_nuclei->SetInput(this->SomaImage);
		sub_volume_filter_nuclei->SetDirectionCollapseToIdentity();
		sub_volume_filter_nuclei->Update();
		sub_volume_nuclei = sub_volume_filter_nuclei->GetOutput();*/

		
		// DO NOT USE DANIELSSON FILTER. IT HAS BUGS. SEE: http://www.itk.org/Bug/view.php?id=10757
		/*typedef itk::DanielssonDistanceMapImageFilter< LabelImageType3D, LabelImageType3D > DanielssonDistanceMapImageFilterType;
		DanielssonDistanceMapImageFilterType::Pointer DianielssonFilter = DanielssonDistanceMapImageFilterType::New();

		DianielssonFilter->SetInput(sub_volume_nuclei);
		DianielssonFilter->SetInputIsBinary(false);
		//DianielssonFilter->GetDistanceMap()->Update();
		DianielssonFilter->GetOutput()->Update();

		LabelImageType3D::Pointer sub_volume_distance = DianielssonFilter->GetOutput();*/

		
		SignedMaurerDistanceMapImageFilterType::Pointer MaurerFilter = SignedMaurerDistanceMapImageFilterType::New();
		MaurerFilter->SetInput(sub_volume_nuclei);
		MaurerFilter->SetSquaredDistance(false);
		MaurerFilter->SetUseImageSpacing(false);
		MaurerFilter->SetInsideIsPositive(false);
		MaurerFilter->Update();
		
		ImageType3D::Pointer sub_volume_distance = MaurerFilter->GetOutput();
		
		
		/*typedef itk::ApproximateSignedDistanceMapImageFilter<CharImageType3D, ImageType3D> ApproxSignedDistanceMapImageFilterType;
		ApproxSignedDistanceMapImageFilterType::Pointer ApproxDistFilter = ApproxSignedDistanceMapImageFilterType::New();
		ApproxDistFilter->SetInput(sub_volume_nuclei);
		ApproxDistFilter->SetInsideValue(255);
		ApproxDistFilter->SetOutsideValue(0);
		ApproxDistFilter->Update();

		ImageType3D::Pointer sub_volume_distance = ApproxDistFilter->GetOutput();*/


		// for testing...
		/*if(i == 200){
			
			RescalerType::Pointer rescaler = RescalerType::New();
			rescaler->SetInput(sub_volume_distance);
			rescaler->SetOutputMaximum(255);
			rescaler->SetOutputMinimum(0);
			rescaler->Update();
			itk::CastImageFilter<ImageType3D, CharImageType3D>::Pointer caster = itk::CastImageFilter< ImageType3D, CharImageType3D>::New();
			caster->SetInput(rescaler->GetOutput());

			std::cout << "Writing distance maps to disk. " << std::endl;
			itk::ImageFileWriter< CharImageType3D >::Pointer distance_map_writer = itk::ImageFileWriter< CharImageType3D >::New();
			distance_map_writer->SetFileName("C:\\Users\\msavelon\\Desktop\\Astro\\TrainingWithBill\\distance_map_sub.tif");
			distance_map_writer->SetInput(caster->GetOutput());
			distance_map_writer->Update();

			itk::ImageFileWriter<CharImageType3D>::Pointer distance_map_writer2 = itk::ImageFileWriter<CharImageType3D>::New();
			distance_map_writer2->SetFileName("C:\\Users\\msavelon\\Desktop\\Astro\\TrainingWithBill\\vol_sub.tif");
			distance_map_writer2->SetInput(sub_volume_nuclei);
			distance_map_writer2->Update();

			LabelImageType3D::SizeType sub_size = sub_volume_nuclei->GetLargestPossibleRegion().GetSize();
			LabelImageType3D::SizeType sub_size1 = sub_volume_nuclei->GetBufferedRegion().GetSize();
			
			std::cout << sub_size[0] << ", " << sub_size[1] << ", " << sub_size[2] << std::endl;
			std::cout << sub_size1[0] << ", " << sub_size1[1] << ", " << sub_size1[2] << std::endl;
			std::cout << current_idx[0] << ", " << current_idx[1] << ", " << current_idx[2] << std::endl;
			std::cout << sub_volume_nuclei->GetOrigin() << std::endl;	
		}*/
		
		//itk::Index<3> lndx = current_idx-starting_index;

		LabelImageType3D::IndexType ldx;
		ldx[0] = double_scale_nuclei; //current_idx[0]; 
		ldx[1] = double_scale_nuclei; //current_idx[1]; 
		ldx[2] = double_scale_nuclei; //current_idx[2]; 

		double nucleus_distance = sub_volume_distance->GetPixel(ldx);
		
		//std::cout << "Root: " << i << " Nuc_dist: " << nucleus_distance << std::endl;


		double nucDistThresh = 100*double_scale_nuclei;
		if(nucleus_distance > nucDistThresh)
			continue;


		////////////////////// Code for eigen values based features ////////////////////////////

		float ballness, plateness, vesselness, noiseness;
	
		PixelType scale_index = this->LoGScaleImage->GetPixel(current_idx);
		ImageType3D::Pointer LoG_sc_image = this->LoG_Vector[scale_index];

		itk::Offset<3>
				xp =  {{2 ,  0 ,   0}},
				xn =  {{-2,  0,    0}},
				yp =  {{0,   2,   0}},
				yn =  {{0,  -2,    0}},
				zp =  {{0,   0,    2}},
				zn =  {{0,   0,   -2}};

		itk::Size<3> rad = {{1,1,1}};
		itk::NeighborhoodIterator<ImageType3D> nit(rad , LoG_sc_image, LoG_sc_image->GetBufferedRegion());
		
		nit.SetLocation(current_idx);
		//std::cout << "After nit definition"<<std::endl;

		//itk::ImageRegionIterator<ImageType3D> it(gauss->GetOutput(), gauss->GetOutput()->GetBufferedRegion());

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

		TensorType h;
		h[0] = LoG_sc_image->GetPixel( current_idx + xp ) + LoG_sc_image->GetPixel( current_idx + xn ) - 2*nit.GetPixel( 13 );
		h[3] = LoG_sc_image->GetPixel( current_idx + yp ) + LoG_sc_image->GetPixel( current_idx + yn ) - 2*nit.GetPixel( 13 );
		h[5] = LoG_sc_image->GetPixel( current_idx + zp ) + LoG_sc_image->GetPixel( current_idx + zn ) - 2*nit.GetPixel( 13 );
		h[1] = nit.GetPixel(xy1) + nit.GetPixel(xy2) - nit.GetPixel(xy3) - nit.GetPixel(xy4);
		h[2] = nit.GetPixel(xz1) + nit.GetPixel(xz2) - nit.GetPixel(xz3) - nit.GetPixel(xz4);
		h[4] = nit.GetPixel(yz1) + nit.GetPixel(yz2) - nit.GetPixel(yz3) - nit.GetPixel(yz4);

		EigenValuesArrayType ev, ev_original;
		EigenVectorMatrixType em;
		h.ComputeEigenAnalysis (ev, em);
		float trace = h.GetTrace();

		//std::cout << "After eigen analysis.open"<<std::endl;

		// Amit-style sorting
		
		bool zeroness = false;
		int lowest_idx = 0;
		
		ev_original[0] = ev[0];
		ev_original[1] = ev[1];
		ev_original[2] = ev[2];

		ev[0] = std::abs(ev[0]);
		ev[1] = std::abs(ev[1]);
		ev[2] = std::abs(ev[2]);

		float L1, L2, L;
		if ( (ev[0] > ev[1]) && (ev[0] > ev[2]) ) 
		{
			L = ev[0];
			L1 = ev[1]; 
			L2 = ev[2];
			if (ev[1] > ev[2])
			{
				L1 = ev[2]; 
				L2 = ev[1];

				lowest_idx = 2;
			}
			else
				lowest_idx = 1;
		}
		else if( (ev[1] > ev[0]) && (ev[1] > ev[2]) ) 
		{
			L = ev[1];
			L1 = ev[0];
			L2 = ev[2];
			if (ev[0] > ev[2]) 
			{
				L1 = ev[2];
				L2 = ev[0];	

				lowest_idx = 2;
			}
			else
				lowest_idx = 0;
		}
		else  
		{
			L = ev[2];
			L1 = ev[0];
			L2 = ev[1];
			if (ev[0] > ev[1]) 
			{
				L1 = ev[1];
				L2 = ev[0];		

				lowest_idx = 1;
			}
			else
				lowest_idx = 0;
		}

		// What is the sequence of the sorted list: L1 < L2 < L

		if(lowest_idx == 0){
			if(ev_original[1] > 0 || ev_original[2] > 0)
				zeroness = true;
		}
		if(lowest_idx == 1){
			if(ev_original[0] > 0 || ev_original[2] > 0)
				zeroness = true;
		}
		if(lowest_idx == 2){
			if(ev_original[0] > 0 || ev_original[1] > 0)
				zeroness = true;
		}
		
		if(zeroness){
			ballness = 0.0;
			plateness = 0.0;
			vesselness = 0.0;
		}
		else{

			float alpha = 0.5 * img_max_val; 
			float beta = 0.5 * img_max_val;
			float gamma = 0.0025 * img_max_val;

			ballness = L1 / sqrt(L2 * L);
			plateness = L2 / L;
			noiseness = std::sqrt(std::pow(L1, 2) + std::pow(L2, 2) + std::pow(L, 2));

			vesselness = (1 - std::exp(-1 * std::pow(plateness, 2) / (2.0 * std::pow(alpha, 2)))) 
				* std::exp(-1 * std::pow(ballness, 2) / (2.0 * std::pow(beta, 2)))
				* (1 - std::exp(-1 * std::pow(noiseness, 2) / (2.0 * std::pow(gamma, 2))));
		}
		

		/////////////////////// Code for intensity based features ///////////////////////////////

		float double_scale = 2 * radius; 
		
		//StatisticsFilterType::Pointer stats_filter = StatisticsFilterType::New();
		ImageType3D::IndexType starting_index, end_index;
		ImageType3D::SizeType sub_volume_size;
		ImageType3D::RegionType sub_volume_region;
		ImageType3D::Pointer sub_volume;

		starting_index[0] = current_idx[0] - double_scale; starting_index[1] = current_idx[1] - double_scale; starting_index[2] = current_idx[2] - double_scale;
		end_index[0] = current_idx[0] + double_scale; end_index[1] = current_idx[1] + double_scale; end_index[2] = current_idx[2] + double_scale;

		
		//std::cout << "Intensity: Starting Indices:"<<starting_index[0]<<" "<<starting_index[1]<<" "<<starting_index[2]<<" "<<std::endl;
		//std::cout << "Intensity: End Indices:"<<end_index[0]<<" "<<end_index[1]<<" "<<end_index[2]<<" "<<std::endl;
		//std::cout << "Intensity: Size:"<<sz[0]<<" "<<sz[1]<<" "<<sz[2]<<" "<<std::endl;

		if ( (starting_index[0] < 0) || (starting_index[1] < 0) || (starting_index[2] < 0) ||
			(end_index[0] > (unsigned int)sz[0]) || (end_index[1] > (unsigned int)sz[1]) ||
			(end_index[2] > (unsigned int)sz[2]) )
		continue;

		sub_volume_size[0] = 2 * double_scale; sub_volume_size[1] = 2 * double_scale; sub_volume_size[2] = 2 * double_scale;
		
		sub_volume_region.SetIndex(starting_index);
		sub_volume_region.SetSize(sub_volume_size);
		
		VolumeOfInterestFilterType::Pointer sub_volume_filter = VolumeOfInterestFilterType::New();
		sub_volume_filter->SetInput(this->PaddedCurvImage);
		sub_volume_filter->SetRegionOfInterest(sub_volume_region);
		sub_volume_filter->Update();
		sub_volume = sub_volume_filter->GetOutput();
		//std::cout << "After subvolume_filter"<<std::endl;

		StatisticsFilterType::Pointer stats_filter = StatisticsFilterType::New();
		stats_filter->SetInput(sub_volume);
		stats_filter->Update();

		double intensity = this->PaddedCurvImage->GetPixel(current_idx);
		double mean_intensity = stats_filter->GetMean();
		double var_intensity = stats_filter->GetVariance();
		double max_intensity = stats_filter->GetMaximum();
		double min_intensity = stats_filter->GetMinimum();
		//std::cout << "After statistics"<<std::endl;

		
		//if(feature_vector.good()){

			if(IDImage->GetPixel(current_idx) == 0){

				//std::cout << "Root: " << i << std::endl; //<< " Ballness: " << ballness << " Plateness: " << plateness << " Vesselness: " << vesselness << std::endl;

				IDImage->SetPixel(current_idx, IDIndex+1);

				//feature_vector << (IDIndex+1) << '\t' << pos[0] << '\t' << pos[1] << '\t' << pos[2]-padz << '\t' << radius << '\t' << ballness << '\t' << plateness << '\t';	
				//feature_vector << intensity << '\t' << mean_intensity << '\t' << var_intensity << '\t' << max_intensity << '\t' << min_intensity << '\t' << nucleus_distance << '\t' ;
				//feature_vector << std::endl;

				radius_vec.push_back(radius);
				int_vec.push_back(intensity);
				min_int_vec.push_back(min_intensity);
				max_int_vec.push_back(max_intensity);
				mean_int_vec.push_back(mean_intensity);
				var_int_vec.push_back(var_intensity);


				CandidateRootPoint a_root;
				a_root.featureVector.node = HeapNode(points_list[i]);
				a_root.featureVector.ID = IDIndex;
				a_root.featureVector.ballness = ballness;
				a_root.featureVector.plateness = plateness;
				a_root.featureVector.intensity = intensity;
				a_root.featureVector.maxIntensity = max_intensity;
				a_root.featureVector.minIntensity = min_intensity;
				a_root.featureVector.meanIntensity = mean_intensity;
				a_root.featureVector.varianceIntensity = var_intensity;
				a_root.featureVector.radius = radius;
				a_root.featureVector.sphere_likelihood = likelihood;
				a_root.featureVector.nucleusDistance = nucleus_distance;
				a_root.featureVector.vesselness = vesselness;

				this->AllRootPoints.push_back(a_root);

				IDIndex++;			
			}
		//}
	}

	// Normalizing the features
	double r_min, r_max, int_min, int_max, mean_min, mean_max, var_min, var_max, min_min, min_max, max_min, max_max;

	r_min = *std::min_element(radius_vec.begin(), radius_vec.end());
	r_max = *std::max_element(radius_vec.begin(), radius_vec.end());
	int_min = *std::min_element(int_vec.begin(), int_vec.end());
	int_max = *std::max_element(int_vec.begin(), int_vec.end());
	mean_min = *std::min_element(mean_int_vec.begin(), mean_int_vec.end());
	mean_max = *std::max_element(mean_int_vec.begin(), mean_int_vec.end());
	var_min = *std::min_element(var_int_vec.begin(), var_int_vec.end());
	var_max = *std::max_element(var_int_vec.begin(), var_int_vec.end());
	min_min = *std::min_element(min_int_vec.begin(), min_int_vec.end());
	min_max = *std::max_element(min_int_vec.begin(), min_int_vec.end());
	max_min = *std::min_element(max_int_vec.begin(), max_int_vec.end());
	max_max = *std::max_element(max_int_vec.begin(), max_int_vec.end());


	std::cout << "Writing root feature vector file. " << std::endl;

	if(feature_vector.good()){
	
		for(int i = 0; i < this->AllRootPoints.size(); i++){
			
			radius_vec[i] = (radius_vec[i] - r_min) / (r_max - r_min);
			int_vec[i] = (int_vec[i] - int_min) / (int_max - int_min);
			min_int_vec[i] = (min_int_vec[i] - min_min) / (min_max - min_min);
			max_int_vec[i] = (max_int_vec[i] - max_min) / (max_max - max_min);
			mean_int_vec[i] = (mean_int_vec[i] - mean_min) / (mean_max - mean_min);
			var_int_vec[i] = (var_int_vec[i] - var_min) / (var_max - var_min);

			this->AllRootPoints[i].featureVector.radius = radius_vec[i];
			this->AllRootPoints[i].featureVector.intensity = int_vec[i];
			this->AllRootPoints[i].featureVector.minIntensity = min_int_vec[i];
			this->AllRootPoints[i].featureVector.maxIntensity = max_int_vec[i];
			this->AllRootPoints[i].featureVector.meanIntensity = mean_int_vec[i];
			this->AllRootPoints[i].featureVector.varianceIntensity = var_int_vec[i];

			//std::cout << "Root: " << i << std::endl;
			
			feature_vector << (i+1) << '\t' << this->AllRootPoints[i].featureVector.node.ndx[0] << '\t' << this->AllRootPoints[i].featureVector.node.ndx[1] << '\t' << this->AllRootPoints[i].featureVector.node.ndx[2]-padz << '\t' ;
			feature_vector << this->AllRootPoints[i].featureVector.radius << '\t'  << this->AllRootPoints[i].featureVector.sphere_likelihood << '\t' << this->AllRootPoints[i].featureVector.ballness << '\t' << this->AllRootPoints[i].featureVector.plateness << '\t';	
			feature_vector << this->AllRootPoints[i].featureVector.vesselness << '\t' << this->AllRootPoints[i].featureVector.intensity << '\t' << this->AllRootPoints[i].featureVector.meanIntensity << '\t' << this->AllRootPoints[i].featureVector.varianceIntensity << '\t';
			feature_vector << this->AllRootPoints[i].featureVector.maxIntensity << '\t' << this->AllRootPoints[i].featureVector.minIntensity << '\t' << this->AllRootPoints[i].featureVector.nucleusDistance << '\t';
			feature_vector << std::endl;
		}
	}


	feature_vector.close();/////////////
	std::cout << "Done with feature vector file. " << std::endl;


	//itk::CastImageFilter< ImageType3D, LabelImageType3D>::Pointer caster = itk::CastImageFilter< ImageType3D,  LabelImageType3D>::New();///////
	//caster->SetInput(IDImage);

	itk::ImageFileWriter< LabelImageType3D >::Pointer IDImageWriter = itk::ImageFileWriter< LabelImageType3D >::New();
	IDImageWriter->SetFileName(IDFname.c_str());
	IDImageWriter->SetInput(IDImage);
	try{
		IDImageWriter->Update();
	}
	catch(itk::ExceptionObject e){
		std::cout << e << std::endl;
		//return EXIT_FAILURE;
	}
	std::cout << "Done with RootsImage writing." << std::endl;
}


bool AstroTracer::PopulateLoGImages(std::vector<float> sigma_vec){

	std::cout << std::endl<< "Computing LoG images" << std::endl;

	//float sigmas[] =  { 2.0f, 2.8284f, 4.0f, 5.6569f, 8.0f, 11.31f };
	//for (unsigned int i = 0; i < 6; ++i)//////////
	for(int i = 0; i < sigma_vec.size(); i++){

		std::cout << "Analysis at " << sigma_vec[i] << std::endl;
		
		float sigma = sigma_vec[i];
			
		LoGFilterType::Pointer gauss = LoGFilterType::New();
		gauss->SetInput(PaddedCurvImage);
		gauss->SetSigma(sigma);
		gauss->SetNormalizeAcrossScale(false);
		gauss->GetOutput()->Update();

		std::cout << "Laplacian of Gaussian at " << sigma << std::endl;

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

		//Store LoG images for later
		this->LoG_Vector.push_back(gauss->GetOutput());
	}

	
	if(this->LoG_Vector.empty())
		return false;
	else
		return true;
	
}

void AstroTracer::UseActiveLearningRootsModel(std::string model_path){
}

HeapNode::HeapNode(){


}

RootPointFeatureVector::RootPointFeatureVector(){

	HeapNode();
}

CandidateRootPoint::CandidateRootPoint(){

	RootPointFeatureVector();
}

void AstroTracer::ReadRootPointsExternal(std::string rootPointsFileName){

	// Currently reading the complete feature vector file along with class_ID and storing the features for class 1 (root points) only.

	std::ifstream rootPoints;
	rootPoints.open(rootPointsFileName.c_str(), std::ios::in); 
	std::cout << "Reading root points file. " << std::endl;
	
	std::vector<std::string> str_vec;
	std::string line, str1;
	if(rootPoints.is_open()){
		
		unsigned short line_number = 0;
		CandidateRootPoint root_point;

		while(rootPoints.good()){

			line_number++;
			//std::cout << line_number << std::endl;

			if(!std::getline(rootPoints, line))
				break;

			//Ignore the first line since it is all text
			if(line_number == 1)
				continue;
			
			std::istringstream str_stream(line); 
			while(str_stream.good()){
				
				if(!getline(str_stream, str1, '\t'))
					break;

				str_vec.push_back(str1);
			}

			root_point.featureVector.ID = atof(str_vec[0].c_str());

			itk::Index<3> idx;
			idx[0] = atof(str_vec[1].c_str());
			idx[1] = atof(str_vec[2].c_str());
			idx[2] = atof(str_vec[3].c_str());
			root_point.featureVector.node = HeapNode(idx, 0);

			root_point.featureVector.radius = atof(str_vec[4].c_str());
			root_point.featureVector.sphere_likelihood = atof(str_vec[5].c_str());
			root_point.featureVector.ballness = atof(str_vec[6].c_str());
			root_point.featureVector.plateness = atof(str_vec[7].c_str());
			root_point.featureVector.vesselness = atof(str_vec[8].c_str());
			root_point.featureVector.intensity = atof(str_vec[9].c_str());
			root_point.featureVector.meanIntensity = atof(str_vec[10].c_str());
			root_point.featureVector.varianceIntensity = atof(str_vec[11].c_str());
			root_point.featureVector.maxIntensity = atof(str_vec[12].c_str());
			root_point.featureVector.minIntensity = atof(str_vec[13].c_str());
			root_point.featureVector.nucleusDistance = atof(str_vec[14].c_str());
			
			// Be careful with the last two. These change depending on how the file was saved from nuc_editor.
			root_point.classValue = atof(str_vec[15].c_str());
			root_point.confidenceMeasure = atof(str_vec[16].c_str());

			// ONLY TWO CLASSES OF ROOT POINTS ARE CONSIDERED
			if(root_point.classValue == 1)
				root_point.isRootPoint = true;
			else 
				root_point.isRootPoint = false;

			if(root_point.isRootPoint)
				this->CandidateRootPoints.push_back(root_point);
			
			str_vec.clear();
		}
		rootPoints.close();
	}
	else{
		std::cout << " Could not open root points file. Exiting now. " << std::endl;
		return;
	}

	if(this->CandidateRootPoints.empty()){
		std::cout << " Empty root points file. Quitting. " << std::endl;
		return;
	}
	std::cout << "Root points file read. " << this->CandidateRootPoints.size() << std::endl;
}

void AstroTracer::GetCentroidsForTracing(std::string outputFname, std::string finalIDImageFileName)
{
	int centroid_count=0;
	
	//Creating centroids.txt file with class 1 roots which are closest to a nucleus 

	//Preparing IDImage
	LabelImageType3D::RegionType id_reg;
	LabelImageType3D::IndexType id_st;
	LabelImageType3D::SizeType id_sz = PaddedCurvImage->GetBufferedRegion().GetSize();

	id_st[0] = 0;
	id_st[1] = 0;
	id_st[2] = 0;
	
	id_reg.SetSize(id_sz);
	id_reg.SetIndex(id_st);
	
	this->FinalRootsImage = LabelImageType3D::New();
	this->FinalRootsImage->SetRegions(id_reg);
	this->FinalRootsImage->Allocate();
	this->FinalRootsImage->SetSpacing(PaddedCurvImage->GetSpacing());

	this->FinalRootsImage->FillBuffer(0);

	std::ofstream centroid_points;
	centroid_points.open(outputFname.c_str(), std::ios::out);

	CharImageType3D::IndexType starting_index_nuclei, end_index_nuclei;
	CharImageType3D::SizeType sub_volume_size_nuclei;
	CharImageType3D::RegionType sub_volume_region_nuclei;
	CharImageType3D::Pointer sub_volume_nuclei;
	
	typedef itk::BinaryThresholdImageFilter<LabelImageType3D, CharImageType3D> ThresholdFilterType;
	ThresholdFilterType::Pointer threshold_filter = ThresholdFilterType::New();
	threshold_filter->SetLowerThreshold(1);
	threshold_filter->SetInsideValue(255);
	threshold_filter->SetOutsideValue(0);
	threshold_filter->SetInput(this->SomaImage);
	threshold_filter->Update();
	
	std::cout << "Root points size: " << this->CandidateRootPoints.size() << std::endl;

	//Loop over nuclei
	for(SIZE_T i = 0; i < this->NucleiObjects.size(); i++){

		// Use only astrocyte nuclei
		if(this->NucleiObjects[i].classValue != 1)
			continue;

		// ROI proportional to nuclei scale, assuming spherical nuclei.
		float double_scale_nuclei = 0.5*std::pow((float)this->NucleiObjects[i].intrinsicFeatures.boundingBoxVolume, (float)0.333333);
		
		CharImageType3D::IndexType current_idx;
		current_idx[0] = this->NucleiObjects[i].intrinsicFeatures.centroid.ndx[0];
		current_idx[1] = this->NucleiObjects[i].intrinsicFeatures.centroid.ndx[1];
		current_idx[2] = this->NucleiObjects[i].intrinsicFeatures.centroid.ndx[2] - padz; 

		starting_index_nuclei[0] = current_idx[0] - double_scale_nuclei; starting_index_nuclei[1] = current_idx[1] - double_scale_nuclei; starting_index_nuclei[2] = current_idx[2] - double_scale_nuclei;
		end_index_nuclei[0] = current_idx[0] + double_scale_nuclei; end_index_nuclei[1] = current_idx[1] + double_scale_nuclei; end_index_nuclei[2] = current_idx[2] + double_scale_nuclei;

		LabelImageType3D::SizeType sz = this->SomaImage->GetBufferedRegion().GetSize();

		if ( (starting_index_nuclei[0] < 0) || (starting_index_nuclei[1] < 0) || (starting_index_nuclei[2] < 0) ||
			(end_index_nuclei[0] > (unsigned int)sz[0]) || (end_index_nuclei[1] > (unsigned int)sz[1]) ||
			(end_index_nuclei[2] > (unsigned int)sz[2]) )
			continue;

		std::cout << "Nuclei: "  << i << " Scale: " << double_scale_nuclei << std::endl;

		sub_volume_size_nuclei[0] = 2 * double_scale_nuclei; sub_volume_size_nuclei[1] = 2 * double_scale_nuclei; sub_volume_size_nuclei[2] = 2 * double_scale_nuclei;

		sub_volume_region_nuclei.SetIndex(starting_index_nuclei);
		sub_volume_region_nuclei.SetSize(sub_volume_size_nuclei);


		typedef itk::RegionOfInterestImageFilter<CharImageType3D, CharImageType3D> VolumeOfInterestFilterType_nuclei2;
		VolumeOfInterestFilterType_nuclei2::Pointer sub_volume_filter_nuclei = VolumeOfInterestFilterType_nuclei2::New();
		sub_volume_filter_nuclei->SetInput(threshold_filter->GetOutput());
		sub_volume_filter_nuclei->SetRegionOfInterest(sub_volume_region_nuclei);
		sub_volume_filter_nuclei->Update();
		sub_volume_nuclei = sub_volume_filter_nuclei->GetOutput();


		SignedMaurerDistanceMapImageFilterType::Pointer MaurerFilter = SignedMaurerDistanceMapImageFilterType::New();
		MaurerFilter->SetInput(sub_volume_nuclei);
		MaurerFilter->SetSquaredDistance(false);
		MaurerFilter->SetUseImageSpacing(false);
		MaurerFilter->SetInsideIsPositive(false);
		MaurerFilter->Update();
		
		ImageType3D::Pointer distance_map = MaurerFilter->GetOutput();
				
		double cur_distance;
		double min_distance = 1000.0;
		//double max_distance = -1000.0; 
		//double mean_distance = 0.0;
		//double variance_distance = 0.0;
		//double acc_distance = 0.0;
		//int n_roots = 0;
		
		std::vector<double> distance_array;
		std::multimap<double, HeapNode> distance_idx_map;
		std::multimap<double, HeapNode>::reverse_iterator map_rit;
		typedef std::pair<double, HeapNode> distance_idx_pair_type;
		
		LabelImageType3D::IndexType min_root_idx;
		min_root_idx[0] = this->CandidateRootPoints[0].featureVector.node.ndx[0];
		min_root_idx[1] = this->CandidateRootPoints[0].featureVector.node.ndx[1];
		min_root_idx[2] = this->CandidateRootPoints[0].featureVector.node.ndx[2] - padz;

		for(SIZE_T j = 0; j < this->CandidateRootPoints.size(); j++){

			//std::cout << "Nuclei: "  << i << " Root: " << j << std::endl;
			

			//if(this->CandidateRootPoints[j].isRootPoint){   //No check needed here since the check has already been done before by selecting class 1 points

				LabelImageType3D::IndexType current_root_idx;

				current_root_idx[0] = this->CandidateRootPoints[j].featureVector.node.ndx[0];
				current_root_idx[1] = this->CandidateRootPoints[j].featureVector.node.ndx[1];
				current_root_idx[2] = this->CandidateRootPoints[j].featureVector.node.ndx[2] - padz;

				int offset = 2; //1;
				if(current_root_idx[0] < starting_index_nuclei[0]+offset || current_root_idx[1] < starting_index_nuclei[1]+offset || current_root_idx[2] < starting_index_nuclei[2]+offset ||
					current_root_idx[0] > end_index_nuclei[0]-offset || current_root_idx[1] > end_index_nuclei[1]-offset || current_root_idx[2] > end_index_nuclei[2]-offset)
					continue;
				
				ImageType3D::IndexType relative_root_idx;
				relative_root_idx[0] = current_root_idx[0] - starting_index_nuclei[0];
				relative_root_idx[1] = current_root_idx[1] - starting_index_nuclei[1];
				relative_root_idx[2] = current_root_idx[2] - starting_index_nuclei[2];

				cur_distance = distance_map->GetPixel(relative_root_idx);
				
				if(cur_distance < 0.00000000000001 && cur_distance > 0.0)
					cur_distance = 0.0;

				if(cur_distance > -0.00000000000001 && cur_distance < 0.0)
					cur_distance = 0.0;
				
				if(cur_distance < -1000.0)
					cur_distance = -1000.0;
				
				if(cur_distance > 1000.0)
					cur_distance = 1000.0;
				
				//std::cout << cur_distance << std::endl;

				if(cur_distance < min_distance)
				{
					min_distance = cur_distance;
					min_root_idx[0] = current_root_idx[0];
					min_root_idx[1] = current_root_idx[1];
					min_root_idx[2] = current_root_idx[2] - padz;
				}

				distance_array.push_back(cur_distance);
				distance_idx_map.insert(distance_idx_pair_type(this->CandidateRootPoints[j].featureVector.radius, HeapNode(current_root_idx, 0)));

		}// End of loop over candidate roots	

		if(!distance_idx_map.empty()){

			// Keep the root point with highest scale
			this->CentroidListForTracing.push_back(distance_idx_map.rbegin()->second);
			this->CentroidScales.push_back(distance_idx_map.rbegin()->first);

			/*if(distance_idx_map.size() == 1){

				this->CentroidListForTracing.push_back(distance_idx_map.begin()->second);

				this->RefinedRootImage->SetPixel(min_root_idx, 255);
			}
			else{

				double cur_scale;
				HeapNode max_scale_node;
				std::multimap<double, HeapNode> max_scale_map(distance_idx_map.rbegin(), distance_idx_map.rbegin()+1);
				
				for(map_rit = distance_idx_map.rbegin()+1; map_rit != distance_idx_map.rend(); map_rit++){	
					
					if(max_scale_map.begin()->first - map_rit->first > 0.01)
						max_scale_map.insert(distance_idx_pair_type(map_rit->first, map_rit->second))
				}
				

				this->CentroidListForTracing.push_back(HeapNode(min_root_idx, 0));
			}*/
		}

		
	}//end of loop over nuclei


	
	// Filter the centroids based on density
	int centroid_nhood = 50; //50; 
	std::vector<bool> neighbor_flags(this->CentroidListForTracing.size(), false);
	
	for(int i = 0; i < this->CentroidListForTracing.size(); i++){

		if(neighbor_flags[i])
			continue;
		
		itk::Index<3> cur_idx = this->CentroidListForTracing[i].ndx;
		std::vector<HeapNode> neighbor_points;
		std::multimap<double, HeapNode> neighbor_point_map;
		std::multimap<double, HeapNode>::iterator neighbor_point_map_itr;
		typedef std::pair<double, HeapNode> distance_idx_pair_type;

		for(int j = 0; j < this->CentroidListForTracing.size(); j++){

			if(j == i)
				continue;

			itk::Index<3> neighbor_idx = this->CentroidListForTracing[j].ndx;

			if(std::abs((int)(cur_idx[0] - neighbor_idx[0])) > centroid_nhood || std::abs((int)(cur_idx[1] - neighbor_idx[1])) > centroid_nhood || 
				std::abs((int)(cur_idx[2] - neighbor_idx[2])) > centroid_nhood)
				continue;
			
			neighbor_flags[j] = true;
			neighbor_points.push_back(HeapNode(neighbor_idx, 0));
			neighbor_point_map.insert(distance_idx_pair_type(this->CentroidScales[j], HeapNode(neighbor_idx, 0)));
		}

		if(neighbor_points.empty()){
			neighbor_flags[i] = false;	
			this->DensityFilteredCentroidListForTracing.push_back(HeapNode(cur_idx, 0));
			continue;
		}
		neighbor_flags[i] = true;
		
		//Replace neighboring root points with their centroids
		/*itk::Index<3> centroid_of_centroid;
		int cx = 0, cy = 0, cz = 0;

		for(int k = 0; k < neighbor_points.size(); k++){
			cx = cx + neighbor_points[k].ndx[0];
			cy = cy + neighbor_points[k].ndx[1];
			cz = cz + neighbor_points[k].ndx[2];
		}
		centroid_of_centroid[0] = cx / neighbor_points.size();
		centroid_of_centroid[1] = cy / neighbor_points.size();
		centroid_of_centroid[2] = cz / neighbor_points.size();
		
		this->DensityFilteredCentroidListForTracing.push_back(HeapNode(centroid_of_centroid, 0));
		*/
		
		// Replace neighboring root points with their scale-weighted centoids
		/*itk::Index<3> centroid_of_centroid;
		int cx = 0, cy = 0, cz = 0;
		double weight_sum = 0;

		for(neighbor_point_map_itr = neighbor_point_map.begin(); neighbor_point_map_itr != neighbor_point_map.end(); neighbor_point_map_itr++){
			cx = cx + neighbor_point_map_itr->first * neighbor_point_map_itr->second.ndx[0];
			cy = cy + neighbor_point_map_itr->first * neighbor_point_map_itr->second.ndx[1];
			cz = cz + neighbor_point_map_itr->first * neighbor_point_map_itr->second.ndx[2];

			weight_sum = weight_sum + neighbor_point_map_itr->first;
		}
		centroid_of_centroid[0] = cx / weight_sum;
		centroid_of_centroid[1] = cy / weight_sum;
		centroid_of_centroid[2] = cz / weight_sum;
		
		this->DensityFilteredCentroidListForTracing.push_back(HeapNode(centroid_of_centroid, 0));
		*/

		//Replace neighboring root points with the highest scale point
		this->DensityFilteredCentroidListForTracing.push_back(neighbor_point_map.rbegin()->second);
	}

	if(centroid_points.good()){

		for(int i = 0; i < this->DensityFilteredCentroidListForTracing.size(); i++){

			this->FinalRootsImage->SetPixel(this->DensityFilteredCentroidListForTracing[i].ndx, 255);

			centroid_points << (float)(this->DensityFilteredCentroidListForTracing[i].ndx[0]) << '\t' 
				<< (float)(this->DensityFilteredCentroidListForTracing[i].ndx[1]) << '\t' 
				<< (float)(this->DensityFilteredCentroidListForTracing[i].ndx[2]) << std::endl;

			centroid_count++;
		}

		/*for(int i = 0; i < this->CentroidListForTracing.size(); i++){

			this->FinalRootsImage->SetPixel(this->CentroidListForTracing[i].ndx, 255);

			centroid_points << (float)(this->CentroidListForTracing[i].ndx[0]) << '\t' 
				<< (float)(this->CentroidListForTracing[i].ndx[1]) << '\t' 
				<< (float)(this->CentroidListForTracing[i].ndx[2]) << std::endl;

			centroid_count++;
		}*/

	}

	
	std::cout << "***Points List Size: " << this->CandidateRootPoints.size() << std::endl;
	std::cout << "***Centroids List Size: " << centroid_count << std::endl;


	itk::ImageFileWriter< LabelImageType3D >::Pointer IDImageWriter = itk::ImageFileWriter< LabelImageType3D >::New();
	IDImageWriter->SetFileName(finalIDImageFileName.c_str());
	IDImageWriter->SetInput(this->FinalRootsImage);
	try{
		IDImageWriter->Update();
	}
	catch(itk::ExceptionObject e){
		std::cout << e << std::endl;
		//return EXIT_FAILURE;
	}
	std::cout << " Done with FinalRootImage writing." << std::endl;

	centroid_points.close();
	//End of creating centroids.txt file

}

IntrinsicFeatureVector::IntrinsicFeatureVector(){
}

AssociativeFeatureVector::AssociativeFeatureVector(){

	this->minRootDist = 100000;
	this->maxRootDist = 100000;
	this->meanRootDist = 100000;
	this->varRootDist = 100000;
	this->nRoots = 0;
}

NucleiObject::NucleiObject(){

	IntrinsicFeatureVector();
	AssociativeFeatureVector();
}


void AstroTracer::ReadNucleiFeaturesExternal(std::string nucleiFeaturesFileName){

	std::ifstream nucleiPoints;
	nucleiPoints.open(nucleiFeaturesFileName.c_str(), std::ios::in); 
	std::cout << "Reading nuclei features file. " << std::endl;
	
	std::vector<std::string> str_vec;
	std::string line, str1;
	if(nucleiPoints.is_open()){
		
		unsigned short line_number = 0;
		NucleiObject nuclei_object;

		while(nucleiPoints.good()){

			line_number++;
			//std::cout << line_number << std::endl;

			if(!std::getline(nucleiPoints, line))
				break;

			//Ignore the first line since it is all text
			if(line_number == 1)
				continue;
			
			std::istringstream str_stream(line); 
			while(str_stream.good()){
				
				if(!getline(str_stream, str1, '\t'))
					break;
				
				// Crazy "nan" is replaced by -1
				//if(strcmp(str1.c_str(), "nan") == 0)
				//	str1 = "-1";

				str_vec.push_back(str1);
			}

			nuclei_object.intrinsicFeatures.ID = atof(str_vec[0].c_str());

			itk::Index<3> idx;
			idx[0] = atof(str_vec[1].c_str());
			idx[1] = atof(str_vec[2].c_str());
			idx[2] = atof(str_vec[3].c_str());
			nuclei_object.intrinsicFeatures.centroid = HeapNode(idx, 0);

			nuclei_object.intrinsicFeatures.volume = atof(str_vec[4].c_str());
			nuclei_object.intrinsicFeatures.integratedIntensity = atof(str_vec[5].c_str());
			nuclei_object.intrinsicFeatures.eccentricity = atof(str_vec[6].c_str());
			nuclei_object.intrinsicFeatures.elongation = atof(str_vec[7].c_str());
			nuclei_object.intrinsicFeatures.boundingBoxVolume = atof(str_vec[9].c_str());

			nuclei_object.intrinsicFeatures.meanIntensity = atof(str_vec[11].c_str());
			nuclei_object.intrinsicFeatures.varianceIntensity = atof(str_vec[15].c_str());
			nuclei_object.intrinsicFeatures.meanSurfaceGradient = atof(str_vec[16].c_str());
			nuclei_object.intrinsicFeatures.radiusVariation = atof(str_vec[22].c_str());
			nuclei_object.intrinsicFeatures.shapeMeasure = atof(str_vec[24].c_str());
			
			nuclei_object.intrinsicFeatures.energy = atof(str_vec[26].c_str());
			nuclei_object.intrinsicFeatures.entropy = atof(str_vec[27].c_str());
			nuclei_object.intrinsicFeatures.inverseDiffMoment = atof(str_vec[28].c_str());
			nuclei_object.intrinsicFeatures.inertia = atof(str_vec[29].c_str());
			nuclei_object.intrinsicFeatures.clusterShade = atof(str_vec[30].c_str());
			nuclei_object.intrinsicFeatures.clusterProminence = atof(str_vec[31].c_str());

			nuclei_object.associativeFeatures.astro_total = atof(str_vec[32].c_str());
			nuclei_object.associativeFeatures.astro_avg = atof(str_vec[33].c_str());
			nuclei_object.associativeFeatures.astro_surr = atof(str_vec[34].c_str());
			nuclei_object.associativeFeatures.micro_total = atof(str_vec[35].c_str());
			nuclei_object.associativeFeatures.micro_avg = atof(str_vec[36].c_str());
			nuclei_object.associativeFeatures.micro_surr = atof(str_vec[37].c_str());
			nuclei_object.associativeFeatures.neuro_total = atof(str_vec[38].c_str());
			nuclei_object.associativeFeatures.neuro_avg = atof(str_vec[39].c_str());
			nuclei_object.associativeFeatures.neuro_surr = atof(str_vec[40].c_str());

			this->NucleiObjects.push_back(nuclei_object);

			str_vec.clear();
		}
		nucleiPoints.close();
	}
	else{
		std::cout << " Could not open nuclei features file. Exiting now. " << std::endl;
		return;
	}

	if(this->NucleiObjects.empty()){
		std::cout << " Empty nuclei features file. Quitting. " << std::endl;
		return;
	}
	std::cout << "Nuclei features file read. " << this->NucleiObjects.size() << std::endl;
}

void AstroTracer::ComputeFeaturesFromCandidateRoots(void){

	// Derive this from: nuclei_object.intrinsicFeatures.boundingBoxVolume
	//float double_scale_nuclei = 25;
	

	CharImageType3D::IndexType starting_index_nuclei, end_index_nuclei;
	CharImageType3D::SizeType sub_volume_size_nuclei;
	CharImageType3D::RegionType sub_volume_region_nuclei;
	CharImageType3D::Pointer sub_volume_nuclei;

	/*std::cout << "Writing distance maps to disk. " << std::endl;
	itk::ImageFileWriter< LabelImageType3D >::Pointer distance_map_writer = itk::ImageFileWriter< LabelImageType3D >::New();
	distance_map_writer->SetFileName("C:\\Users\\msavelon\\Desktop\\Astro\\TrainingWithBill\\distance_map.tif");
	distance_map_writer->SetInput(nucleus_distance);
	distance_map_writer->Update();

	itk::ImageFileWriter< LabelImageType3D >::Pointer distance_map_writer3 = itk::ImageFileWriter< LabelImageType3D >::New();
	distance_map_writer3->SetFileName("C:\\Users\\msavelon\\Desktop\\Astro\\TrainingWithBill\\voronoi_map.tif");
	distance_map_writer3->SetInput(voronoi_map);
	distance_map_writer3->Update();*/
	
	
	typedef itk::BinaryThresholdImageFilter<LabelImageType3D, CharImageType3D> ThresholdFilterType;
	ThresholdFilterType::Pointer threshold_filter = ThresholdFilterType::New();
	threshold_filter->SetLowerThreshold(1);
	threshold_filter->SetInsideValue(255);
	threshold_filter->SetOutsideValue(0);
	threshold_filter->SetInput(this->SomaImage);
	threshold_filter->Update();
	
	/*itk::ImageFileWriter< CharImageType3D >::Pointer nuclei_writer = itk::ImageFileWriter<CharImageType3D>::New();
	nuclei_writer->SetFileName("C:\\Users\\msavelon\\Desktop\\Astro\\TrainingWithBill\\Binary_nuclei.tif");
	nuclei_writer->SetInput(threshold_filter->GetOutput());
	nuclei_writer->Update();*/
	
	std::cout << "Root points size: " << this->CandidateRootPoints.size() << std::endl;
	
	std::cout << "Computing nuclei features. " << std::endl;

	//Loop over nuclei
	for(SIZE_T i = 0; i < this->NucleiObjects.size(); i++){

		// ROI proportional to nuclei scale, assuming spherical nuclei.
		float double_scale_nuclei = 0.5*std::pow((float)this->NucleiObjects[i].intrinsicFeatures.boundingBoxVolume, (float)0.333333);

		
		CharImageType3D::IndexType current_idx;
		current_idx[0] = this->NucleiObjects[i].intrinsicFeatures.centroid.ndx[0];
		current_idx[1] = this->NucleiObjects[i].intrinsicFeatures.centroid.ndx[1];
		current_idx[2] = this->NucleiObjects[i].intrinsicFeatures.centroid.ndx[2] - padz; 

		starting_index_nuclei[0] = current_idx[0] - double_scale_nuclei; starting_index_nuclei[1] = current_idx[1] - double_scale_nuclei; starting_index_nuclei[2] = current_idx[2] - double_scale_nuclei;
		end_index_nuclei[0] = current_idx[0] + double_scale_nuclei; end_index_nuclei[1] = current_idx[1] + double_scale_nuclei; end_index_nuclei[2] = current_idx[2] + double_scale_nuclei;

		LabelImageType3D::SizeType sz = this->SomaImage->GetBufferedRegion().GetSize();

		//std::cout << "Nuclei: Starting Indices:"<<starting_index_nuclei[0]<<" "<<starting_index_nuclei[1]<<" "<<starting_index_nuclei[2]<<" "<<std::endl;
		//std::cout << "Nuclei: End Indices:"<<end_index_nuclei[0]<<" "<<end_index_nuclei[1]<<" "<<end_index_nuclei[2]<<std::endl;

	    //std::cout << "Nuclei: "  << i << " Scale: " << double_scale_nuclei << std::endl;

		if ( (starting_index_nuclei[0] < 0) || (starting_index_nuclei[1] < 0) || (starting_index_nuclei[2] < 0) ||
			(end_index_nuclei[0] > (unsigned int)sz[0]) || (end_index_nuclei[1] > (unsigned int)sz[1]) ||
			(end_index_nuclei[2] > (unsigned int)sz[2]) )
			continue;

		//std::cout << "Nuclei: "  << i << " Scale: " << double_scale_nuclei << std::endl;

		sub_volume_size_nuclei[0] = 2 * double_scale_nuclei; sub_volume_size_nuclei[1] = 2 * double_scale_nuclei; sub_volume_size_nuclei[2] = 2 * double_scale_nuclei;

		sub_volume_region_nuclei.SetIndex(starting_index_nuclei);
		sub_volume_region_nuclei.SetSize(sub_volume_size_nuclei);


		typedef itk::RegionOfInterestImageFilter<CharImageType3D, CharImageType3D> VolumeOfInterestFilterType_nuclei2;
		VolumeOfInterestFilterType_nuclei2::Pointer sub_volume_filter_nuclei = VolumeOfInterestFilterType_nuclei2::New();
		sub_volume_filter_nuclei->SetInput(threshold_filter->GetOutput());
		sub_volume_filter_nuclei->SetRegionOfInterest(sub_volume_region_nuclei);
		sub_volume_filter_nuclei->Update();
		sub_volume_nuclei = sub_volume_filter_nuclei->GetOutput();


		SignedMaurerDistanceMapImageFilterType::Pointer MaurerFilter = SignedMaurerDistanceMapImageFilterType::New();
		MaurerFilter->SetInput(sub_volume_nuclei);
		MaurerFilter->SetSquaredDistance(false);
		MaurerFilter->SetUseImageSpacing(false);
		MaurerFilter->SetInsideIsPositive(false);
		MaurerFilter->Update();
		
		ImageType3D::Pointer distance_map = MaurerFilter->GetOutput();
		
				
		double cur_distance;
		double min_distance = 1000.0;
		double max_distance = -1000.0; 
		double mean_distance = 0.0;
		double variance_distance = 0.0;
		double acc_distance = 0.0;
		int n_roots = 0;
		
		//std::cout << std::endl;

		std::vector<double> distance_array;
		for(SIZE_T j = 0; j < this->CandidateRootPoints.size(); j++){

			//std::cout << "Nuclei: "  << i << " Root: " << j << std::endl;
			

			if(CandidateRootPoints[j].isRootPoint){

				LabelImageType3D::IndexType current_root_idx;
				current_root_idx[0] = this->CandidateRootPoints[j].featureVector.node.ndx[0];
				current_root_idx[1] = this->CandidateRootPoints[j].featureVector.node.ndx[1];
				current_root_idx[2] = this->CandidateRootPoints[j].featureVector.node.ndx[2] - padz;

				//if(i > 880)
				//	std::cout << "Nuclei: "  << i << " Root: " << j << std::endl;
				
				int offset = 2; //1;
				if(current_root_idx[0] < starting_index_nuclei[0]+offset || current_root_idx[1] < starting_index_nuclei[1]+offset || current_root_idx[2] < starting_index_nuclei[2]+offset ||
					current_root_idx[0] > end_index_nuclei[0]-offset || current_root_idx[1] > end_index_nuclei[1]-offset || current_root_idx[2] > end_index_nuclei[2]-offset)
					continue;


				//std::cout << "Nuclei: "  << i << " Root: " << j << std::endl;
				
				
				ImageType3D::IndexType relative_root_idx;
				relative_root_idx[0] = current_root_idx[0] - starting_index_nuclei[0];
				relative_root_idx[1] = current_root_idx[1] - starting_index_nuclei[1];
				relative_root_idx[2] = current_root_idx[2] - starting_index_nuclei[2];

				//std::cout << relative_root_idx[0] << ", " << relative_root_idx[1] << ", " << relative_root_idx[2] << std::endl;
				
				cur_distance = distance_map->GetPixel(relative_root_idx);
				
				if(cur_distance < 0.00000000000001 && cur_distance > 0.0)
					cur_distance = 0.0;

				if(cur_distance > -0.00000000000001 && cur_distance < 0.0)
					cur_distance = 0.0;
				
				if(cur_distance < -1000.0)
					cur_distance = -1000.0;
				
				if(cur_distance > 1000.0)
					cur_distance = 1000.0;
				
				//std::cout << cur_distance << std::endl;

				if(cur_distance < min_distance)
					min_distance = cur_distance;
				if(cur_distance > max_distance)
					max_distance = cur_distance;

				distance_array.push_back(cur_distance);
				acc_distance = acc_distance + cur_distance;
				n_roots++;
			}
		}	

		//std::cout << "acc_dist: " << acc_distance << " n_roots: " << n_roots << std::endl;

		if(!distance_array.empty()){

			//std::cout << "acc_dist: " << acc_distance << " n_roots: " << n_roots << std::endl;

			mean_distance = acc_distance / n_roots;
			for(int k = 0; k < distance_array.size(); k++)
				variance_distance = variance_distance + std::pow(distance_array[k] - mean_distance, 2);
		
			variance_distance = variance_distance / n_roots;

			//std::cout << "mean: " << mean_distance << " variance: " << variance_distance << " n_roots: " << n_roots << std::endl;

			this->NucleiObjects[i].associativeFeatures.minRootDist = min_distance; 
			this->NucleiObjects[i].associativeFeatures.maxRootDist = max_distance;
			this->NucleiObjects[i].associativeFeatures.meanRootDist = mean_distance;
			this->NucleiObjects[i].associativeFeatures.varRootDist = variance_distance;
			this->NucleiObjects[i].associativeFeatures.nRoots = n_roots;

			//distance_array.clear();
		}
		else{
			// These features do not exist when no root points are found in the neighborhood
			this->NucleiObjects[i].associativeFeatures.minRootDist = 100000; //double_scale_nuclei;
			this->NucleiObjects[i].associativeFeatures.maxRootDist = 100000; //double_scale_nuclei;
			this->NucleiObjects[i].associativeFeatures.meanRootDist = 100000; //double_scale_nuclei;
			this->NucleiObjects[i].associativeFeatures.varRootDist = 100000; //double_scale_nuclei;
			this->NucleiObjects[i].associativeFeatures.nRoots = 0;

		}//loop over candidate roots

		//distance_array.clear();

	}//end of loop over nuclei
}

void AstroTracer::WriteNucleiFeatures(std::string outputFname){
	
	std::ofstream nuclei_feature_vector;
	nuclei_feature_vector.open(outputFname.c_str(), std::ios::out);
	//std::cout << "After feature_vector.open"<<std::endl;
	
	unsigned short IDIndex = 0;//ID index
	nuclei_feature_vector << "ID" << '\t' << "x" << '\t' << "y" << '\t' << "z" << '\t' << "volume" << '\t' << "sum_int" << '\t' << "mean_int" << '\t';	
	nuclei_feature_vector << "var_int" << '\t' << "eccentricity" << '\t' << "elongation" << '\t' << "mean_surf_gradient" << '\t' << "radius_variation" << '\t'; 
	nuclei_feature_vector << "shape_measure" << '\t' << "energy" << '\t' << "entropy" << '\t' << "inverse_diff_moment" << '\t' << "inertia" << '\t';
	nuclei_feature_vector << "cluster_shade" << '\t' << "cluster_prominence" << '\t' << "Astrocyte_TOTAL" << '\t' << "Astrocyte_AVG" << '\t' << "Astrocyte_SURR" << '\t';
	nuclei_feature_vector << "Microglia_TOTAL" << '\t' << "Microglia_AVG" << '\t' << "Microglia_SURR" << '\t' << "Neurons_TOTAL" << '\t' << "Neurons_AVG" << '\t';
	nuclei_feature_vector << "Neurons_SURR" << '\t' << "min_root_dist" << '\t' << "max_root_dist" << '\t' << "mean_root_dist" << '\t' << "var_root_dist" << '\t' << "n_roots" << '\t';
	nuclei_feature_vector << std::endl;

	for(SIZE_T i = 0; i < this->NucleiObjects.size(); i++){

		NucleiObject nuc = this->NucleiObjects[i];
		
		//if(nuc.associativeFeatures.nRoots != 0){
			nuclei_feature_vector << nuc.intrinsicFeatures.ID << '\t' << nuc.intrinsicFeatures.centroid.ndx[0] << '\t' << nuc.intrinsicFeatures.centroid.ndx[1] << '\t' << nuc.intrinsicFeatures.centroid.ndx[2] << '\t';
			nuclei_feature_vector << nuc.intrinsicFeatures.volume << '\t' << nuc.intrinsicFeatures.integratedIntensity << '\t' << nuc.intrinsicFeatures.meanIntensity << '\t' << nuc.intrinsicFeatures.varianceIntensity << '\t';
			nuclei_feature_vector << nuc.intrinsicFeatures.eccentricity << '\t' << nuc.intrinsicFeatures.elongation << '\t' << nuc.intrinsicFeatures.meanSurfaceGradient << '\t' << nuc.intrinsicFeatures.radiusVariation << '\t';
			nuclei_feature_vector << nuc.intrinsicFeatures.shapeMeasure << '\t' << nuc.intrinsicFeatures.energy << '\t' << nuc.intrinsicFeatures.entropy << '\t' << nuc.intrinsicFeatures.inverseDiffMoment << '\t';
			nuclei_feature_vector << nuc.intrinsicFeatures.inertia << '\t' << nuc.intrinsicFeatures.clusterShade << '\t' << nuc.intrinsicFeatures.clusterProminence << '\t'; 
			nuclei_feature_vector << nuc.associativeFeatures.astro_total << '\t' << nuc.associativeFeatures.astro_avg << '\t' << nuc.associativeFeatures.astro_surr << '\t';
			nuclei_feature_vector << nuc.associativeFeatures.micro_total << '\t' << nuc.associativeFeatures.micro_avg << '\t' << nuc.associativeFeatures.micro_surr << '\t';
			nuclei_feature_vector << nuc.associativeFeatures.neuro_total << '\t' << nuc.associativeFeatures.neuro_avg << '\t' << nuc.associativeFeatures.neuro_surr << '\t';
			nuclei_feature_vector << nuc.associativeFeatures.minRootDist << '\t' << nuc.associativeFeatures.maxRootDist << '\t' << nuc.associativeFeatures.meanRootDist << '\t';
			nuclei_feature_vector << nuc.associativeFeatures.varRootDist << '\t' << nuc.associativeFeatures.nRoots << '\t';
			nuclei_feature_vector << std::endl;
		//}
	}
	nuclei_feature_vector.close();
	std::cout << "Done with nuclei feature vector file." << std::endl;
}

void AstroTracer::ReadFinalNucleiTable(std::string finalNucleiTableFileName){

	std::ifstream nucleiPoints;
	nucleiPoints.open(finalNucleiTableFileName.c_str(), std::ios::in); 
	std::cout << "Reading final nuclei table file. " << std::endl;
	
	std::vector<std::string> str_vec;
	std::string line, str1;
	if(nucleiPoints.is_open()){
		
		unsigned short line_number = 0;
		NucleiObject nuclei_object;

		while(nucleiPoints.good()){

			line_number++;
			//std::cout << line_number << std::endl;

			if(!std::getline(nucleiPoints, line))
				break;

			//Ignore the first line since it is all text
			if(line_number == 1)
				continue;
			
			std::istringstream str_stream(line); 
			while(str_stream.good()){
				
				if(!getline(str_stream, str1, '\t'))
					break;
				
				// Crazy "nan" is replaced by -1
				//if(strcmp(str1.c_str(), "nan") == 0)
				//	str1 = "-1";

				str_vec.push_back(str1);
			}

			nuclei_object.intrinsicFeatures.ID = atof(str_vec[0].c_str());

			itk::Index<3> idx;
			idx[0] = atof(str_vec[1].c_str());
			idx[1] = atof(str_vec[2].c_str());
			idx[2] = atof(str_vec[3].c_str());
			nuclei_object.intrinsicFeatures.centroid = HeapNode(idx, 0);

			nuclei_object.intrinsicFeatures.volume = atof(str_vec[4].c_str());
			nuclei_object.intrinsicFeatures.integratedIntensity = atof(str_vec[5].c_str());
			nuclei_object.intrinsicFeatures.eccentricity = atof(str_vec[6].c_str());
			nuclei_object.intrinsicFeatures.elongation = atof(str_vec[7].c_str());
			nuclei_object.intrinsicFeatures.boundingBoxVolume = atof(str_vec[9].c_str());

			nuclei_object.intrinsicFeatures.meanIntensity = atof(str_vec[11].c_str());
			nuclei_object.intrinsicFeatures.varianceIntensity = atof(str_vec[15].c_str());
			nuclei_object.intrinsicFeatures.meanSurfaceGradient = atof(str_vec[16].c_str());
			nuclei_object.intrinsicFeatures.radiusVariation = atof(str_vec[22].c_str());
			nuclei_object.intrinsicFeatures.shapeMeasure = atof(str_vec[24].c_str());
			
			nuclei_object.intrinsicFeatures.energy = atof(str_vec[26].c_str());
			nuclei_object.intrinsicFeatures.entropy = atof(str_vec[27].c_str());
			nuclei_object.intrinsicFeatures.inverseDiffMoment = atof(str_vec[28].c_str());
			nuclei_object.intrinsicFeatures.inertia = atof(str_vec[29].c_str());
			nuclei_object.intrinsicFeatures.clusterShade = atof(str_vec[30].c_str());
			nuclei_object.intrinsicFeatures.clusterProminence = atof(str_vec[31].c_str());


			nuclei_object.classValue = atof(str_vec[47].c_str());
			nuclei_object.confidenceMeasure = atof(str_vec[48].c_str());


			this->NucleiObjects.push_back(nuclei_object);

			str_vec.clear();
		}
		nucleiPoints.close();
	}
	else{
		std::cout << " Could not open nuclei features file. Exiting now. " << std::endl;
		return;
	}

	if(this->NucleiObjects.empty()){
		std::cout << " Empty nuclei features file. Quitting. " << std::endl;
		return;
	}
	std::cout << "Nuclei features file read. " << this->NucleiObjects.size() << std::endl;
}