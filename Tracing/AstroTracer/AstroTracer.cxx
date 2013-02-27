#include "AstroTracer.h"
#include "time.h"



#ifdef _OPENMP
#include "omp.h"
#endif

AstroTracer::AstroTracer(){
	this->isCoverageOptimized = false;
	this->VBT = new ftkVesselTracer;
}


void AstroTracer::SetInputDataPath(const std::string path){
	this->InputDataPath = path;
}
void AstroTracer::LoadParameters(const char* parametersFileName)
{
	std::map<std::string, std::string> opts;  
	this->optionsCreate(parametersFileName, opts);
	
	std::map<std::string,std::string>::iterator mi;

	mi = opts.find("-intensity_threshold"); 
	if(!this->isCoverageOptimized){
		if(mi!=opts.end())
		{ std::istringstream ss((*mi).second); ss>>this->intensity_threshold; 
		}
		else
		{ this->intensity_threshold = 0.005; printf("Chose intensity_threshold = 0.005 as default\n");}
	}

	mi = opts.find("-contrast_threshold");
	if(!this->isCoverageOptimized){
		if(mi!=opts.end())
		{ std::istringstream ss((*mi).second); ss>>this->contrast_threshold; }
		else
		{	  this->contrast_threshold = 0.0003; printf("Chose contrast_threshold = 0.0003 as default\n"); }
	}

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
void AstroTracer::LoadCurvImageFromPath(std::string fname, unsigned int pad) 
{
	std::cout << "Reading input file "<< fname << std::endl;
	ReaderType::GlobalWarningDisplayOff();
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(fname);
	ImageType3D::Pointer image = reader->GetOutput();
	image->Update();

	std::cout << "Entering LoadCurvImage" << std::endl;
	LoadCurvImage(image, pad);
}

void AstroTracer::LoadCurvImage(ImageType3D::Pointer &image, unsigned int pad)  
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

void AstroTracer::DoPreprocessing(void){
	
	this->VBT->PreprocessData(this->InputDataPath, this->PaddedCurvImage, false, false);
}

void AstroTracer::LoadPreprocessingResults(){

	typedef itk::ImageFileReader<ImageType3D> ReaderType;
	
	ReaderType::Pointer reader1 = ReaderType::New();
	reader1->SetFileName(this->InputDataPath + "_pre.mhd");
	reader1->Update();
	this->PaddedCurvImage = reader1->GetOutput();
	this->VBT->SetInputImage(this->PaddedCurvImage);

	ReaderType::Pointer reader2 = ReaderType::New();
	reader2->SetFileName(this->InputDataPath + "_gx.mhd");
	reader2->Update();
	this->gx = reader2->GetOutput();
	
	ReaderType::Pointer reader3 = ReaderType::New();
	reader3->SetFileName(this->InputDataPath + "_gy.mhd");
	reader3->Update();
	this->gy = reader3->GetOutput();

	ReaderType::Pointer reader4 = ReaderType::New();
	reader4->SetFileName(this->InputDataPath + "_gz.mhd");
	reader4->Update();
	this->gz = reader4->GetOutput();
	
	this->VBT->SetGVFImages(gx, gy, gz);

	ReaderType::Pointer reader5 = ReaderType::New();
	reader5->SetFileName(this->InputDataPath + "_vesselness.mhd");
	reader5->Update();
	this->VesselnessImage = reader5->GetOutput();
	this->VBT->SetVesselnessImage(this->VesselnessImage);

	this->VBT->NormalizeAndRescaleData();
	this->VBT->SphericalBinPreprocess();
	this->VBT->Set_useVesselness(1);
}

void AstroTracer::ReadStartPointsFromPath(std::string fname, unsigned int pad) 
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

void AstroTracer::ReadStartPoints(std::vector< itk::Index<3> > somaCentroids, unsigned int pad) 
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

	////////////////////////////////

	// Compute vesselness image for tracing
	float sigma_min = 1.0;
	float sigma_max = 10.0;
	int sigma_intervals = 10;

	StatisticsFilterType::Pointer stats_filter = StatisticsFilterType::New();
	stats_filter->SetInput(this->PaddedCurvImage);
	stats_filter->Update();
	double img_max_val = stats_filter->GetMaximum();
	
	std::cout << "Max img value: " << img_max_val << std::endl;
	
	ObjectnessMeasures obj_measures(sigma_min, sigma_max, sigma_intervals, 1);
	obj_measures.alpha = 0.5 * img_max_val;
	obj_measures.beta = 0.5 * img_max_val;
	obj_measures.gamma = 0.25 * img_max_val; //0.25 * img_max_val;

	this->ComputeObjectnessImage(obj_measures);

	StatisticsFilterType::Pointer stats_filter2 = StatisticsFilterType::New();
	stats_filter2->SetInput(this->ObjectnessImage);
	stats_filter2->Update();
	double max_objectness = stats_filter2->GetMaximum();

	std::cout << "Max objectness value: " << max_objectness << std::endl;

	MultiplyImageFilter::Pointer img_multiplier = MultiplyImageFilter::New();
	img_multiplier->SetInput(this->ObjectnessImage);
	img_multiplier->SetConstant2(1.0/(max_objectness));
	img_multiplier->Update();

	MultiplyImageFilter::Pointer img_multiplier2 = MultiplyImageFilter::New();
	img_multiplier2->SetInput(this->PaddedCurvImage);
	img_multiplier2->SetConstant2(1.0/(img_max_val));
	img_multiplier2->Update();

	this->PaddedCurvImage = img_multiplier2->GetOutput();
	this->ObjectnessImage = img_multiplier->GetOutput();

	// Create a vesselness-intensity hybrid image
	MultiplyImageFilter::Pointer img_multiplier3 = MultiplyImageFilter::New();
	img_multiplier3->SetInput1(this->PaddedCurvImage);
	img_multiplier3->SetInput2(this->ObjectnessImage);
	img_multiplier3->Update();
	
	this->ObjectnessHybridImage = img_multiplier3->GetOutput();

	// Write the images out
	/*typedef itk::ImageFileWriter<ImageType3D> ImageWriterType;
	ImageWriterType::Pointer image_writer = ImageWriterType::New();
	image_writer->SetFileName("F:\\SfN12_posterChild\\1\\VesselnessImage_astro.mhd");
	image_writer->SetInput(this->ObjectnessImage);
	image_writer->Update();

	ImageWriterType::Pointer image_writer1 = ImageWriterType::New();
	image_writer1->SetFileName("F:\\SfN12_posterChild\\1\\Paddedcurv_astro.mhd");
	image_writer1->SetInput(this->PaddedCurvImage);
	image_writer1->Update();

	ImageWriterType::Pointer image_writer2 = ImageWriterType::New();
	image_writer2->SetFileName("F:\\SfN12_posterChild\\1\\Hybrid_vesselness.mhd");
	image_writer2->SetInput(this->ObjectnessHybridImage);
	image_writer2->Update();
	
	std::cout << "Done with writing vesselness image. " << std::endl;*/

	
	// Do traing on the vesselness image directly!!
	//this->PaddedCurvImage = this->ObjectnessImage;




	////////////////////////////////////////

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
	
	// Initialize containers with start points - call this primary nodes. TreeID is positive for primary nodes.
	clock_t fillSWCImage1_start_time = clock();
	for (startIt = StartPoints.begin(); startIt != StartPoints.end(); ++startIt, ++tID)
	{
		itk::Index<3> startIndex = (*startIt);
		startIndex[2] += padz;													//Convert to padded image index
		SWCNode_astro* start_node = new SWCNode_astro(CurrentID++, -1, tID, startIndex);	//This is the seed points SWCNode_astro
		SWCImage->SetPixel(startIndex,start_node);								//Adding all seed points to the SWCImage
		ConnImage->SetPixel(startIndex,0.0f);									//Set the ConnectedImage to 0.0 at all the seed nodes (remember that the Connected image is all initialized with MAXVAL)... 
		SWCNode_astroContainer.push_back(start_node);									//Fill the _SWCNode_astroContainer with start points
		HeapNode_astro *h = new HeapNode_astro(start_node->ndx, 0.0);						//Heap nodes hold an (index, value) pair
		PQ.push(h);																//Priority Queue contains the seed nodes now...
	}
	std::cout << "fillSWCImage1 took: " << (clock() - fillSWCImage1_start_time)/(float) CLOCKS_PER_SEC << std::endl;

	clock_t fillSWCImage2_start_time = clock();
	
	// Putting all the "interest" points in SWCImage - call this secondary nodes. TreeID is negetive for secondary nodes.
	long eCounter = 0, TotalePoints;
	itk::ImageRegionConstIterator<ImageType3D> Nit(NDXImage, NDXImage->GetBufferedRegion());
	for (Nit.GoToBegin(); !Nit.IsAtEnd(); ++Nit) 
	{
		if (Nit.Get() > 0)	//Vesselness value is greater than 0
		{
			itk::Index<3> endx = Nit.GetIndex();
			SWCNode_astro* s2 = new SWCNode_astro(0, -1, -1*(++eCounter), endx);	//id = 0, parent_id = -1, tree id = -1 * eCounter, index that this vesselness value is greater than 0
			SWCImage->SetPixel(endx,s2);								//Adding all critical points where vesselness value is greater than 0 to the SWC image
		}
	}
	std::cout << "fillSWCImage2 took: " << (clock() - fillSWCImage2_start_time)/(float) CLOCKS_PER_SEC << std::endl;


	TotalePoints = eCounter;
	std::cout << "eCounter = " << eCounter << std::endl;	//eCounter is just number of nodes that are critical points (but not seed points)
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

	// Start tracing by popping from the priority queue. This is initialized with start points only
	while(!PQ.empty())	//For each seed node
	{
		//Take the top HeapNode_astro and remove it from the Priority Queue 
		HeapNode_astro *h = PQ.top();
		PQ.pop();

		//Temporarily store the index and value of the node
		itk::Index<3> ndx = h->ndx;
		KeyValue = h->KeyValue;
		delete h;

		//Don't do anything if the HeapNode_astro value is larger than the one in the connected image
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
			
			SWCNode_astro* t  = SWCImage->GetPixel(ndx);
			if ( t != NULL) 
			{
				if (t->TreeID < 0) 
				{
					delete t;
				}
			}
			continue;
		}

		// Initially SWCImage contains only primary and secondary nodes
		SWCNode_astro* s = SWCImage->GetPixel(ndx);
		// Go in if the current node is either primary or secondary
		if (s != NULL) 
		{
			// Go in if the current node is a secondary node
			if (s->TreeID < 0) 
			{
				std::vector<IndexType> Chain;
				
				// Trace back from a secondary node.
				// The cost function inside TBack is used to decide whether or not the Chain will join the tree
				SWCNode_astro* L = TBack(ndx, Chain);

				if ( L  != NULL ) 
				{
					// By influencing this cost, atleast one "bad" node will be added to the tree.
					// Consider placing a hard cost-based threshold to avoid crazy nodes being added
					// This cost will only slowly discourage traces from going into the high cost areas. This cost is between 0 and 1, 0: will encourage 
					// tracing in this area, 1: will not affect the cost, meaning it will be decided by the intensity cost alone.
					
					float costFactor = GetCost(L, ndx);
					//float costFactor = GetCostOld(L, ndx);
					//float costFactor = GetCostLocal( L , ndx);
					//float costFactor = GetCostLocal2(L , ndx);

					//if(costFactor > 1.5)
					//	continue;
					
					std::vector<IndexType>::reverse_iterator cit;
					SWCNode_astro* par = L;

					for (cit = Chain.rbegin(); cit != Chain.rend(); ++cit) 
					{
						SWCNode_astro* t = SWCImage->GetPixel(*cit);
						if (t == NULL) 
						{
							float val = ConnImage->GetPixel(*cit) * costFactor;
							ConnImage->SetPixel((*cit),val);
							SWCNode_astro* s = new SWCNode_astro(CurrentID++, par, L->TreeID, (*cit));
							SWCImage->SetPixel((*cit),s);
							SWCNode_astroContainer.push_back(s);
							par->children.push_back(s);
							par = s;
							HeapNode_astro *h = new HeapNode_astro((*cit), val);
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
								SWCNode_astro* s = new SWCNode_astro(CurrentID++, par, L->TreeID, (*cit));
								SWCImage->SetPixel((*cit),s);
								SWCNode_astroContainer.push_back(s);
								par->children.push_back(s);
								//std::cout<<"SWCImage Node @ " << (*cit) << "(" << s->ID << ") with parent " << par->ID << "  Cost: " << val << "  " << (100*eCounter)/TotalePoints << "% Remaining.\r";// << std::endl;
								par = s;
								HeapNode_astro *h = new HeapNode_astro((*cit), val);
								PQ.push(h);
							}
						}
					}
				} 
			}
		}
		
		// Consider adding the neighbors of ndx to the queue for tracing. This is done until we hit a primary or secondary node which is already present in SWCImage.
		// This is because of the condition "if(s != NULL)" above.
		for (oit = off.begin(); oit < off.end(); ++oit) 
		{
			itk::Index<3> ndx2 = ndx + (*oit);
			if ( (ndx2[0] < 2) || (ndx2[1] < 2) || (ndx2[2] < 2) || (ndx2[0] >= unsigned(size[0] - 2)) || (ndx2[1] >= unsigned(size[1] - 2)) || (ndx2[2] >= unsigned(size[2] - 2)) )  
				continue;
			
			// If the neighbouring point is already a node (primary/secondary) and if it has been assigned to a tree, then do nothing for it.
			if (SWCImage->GetPixel(ndx2) != NULL) 
			{
				if (SWCImage->GetPixel(ndx2)->TreeID > 0) 
				{
					continue;			
				}
			}

			// This block tries to estimate the cost of pushing a neighboring point to the queue. The cost depends non-linearly on the pixel value at ndx2
			// and the values in the ConnImage at the minimum-cost neighbors in each direction (a1, a2, a3) of ndx2. The initial cost, for the very first
			// point, is just inversely proportional to the pixel value at ndx2 (P).
			PixelType P = 1/(PaddedCurvImage->GetPixel(ndx2) + 0.001f);  // consider taking inverse here
			PixelType a1, a2, a3;

			//a1, a2, a3 are the minimum cost (per ConnImage) neighbors of ndx2 in each direction
			ScanNeighbors(a1,a2,a3, ndx2);

			// This function returns the cost of adding ndx2 to the queue to further consider it for tracing
			PixelType aa = Update( a1, a2, a3, P );

			// Check if the existing cost is greater than the current cost. Go for the lower cost and push ndx2 to the queue
			if ( ConnImage->GetPixel(ndx2) > aa )  
			{
				ConnImage->SetPixel(ndx2, aa);
				HeapNode_astro *h = new HeapNode_astro(ndx2, aa);
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

	//clock_t RemoveIntraSomaNodes_start_time = clock();
	//RemoveIntraSomaNodes();
	//std::cout << "RemoveIntraSomaNodes took: " << (clock() - RemoveIntraSomaNodes_start_time)/(float) CLOCKS_PER_SEC << std::endl;

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

	//NDXImage2 = ImageType3D::New();///////////////////
	//NDXImage2->SetRegions(PaddedCurvImage->GetBufferedRegion());//////////////////////
	//NDXImage2->Allocate();/////////////////////////
	//NDXImage2->FillBuffer(0.0f);///////////////////////////
		
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


void AstroTracer::OptimizeCoverage(std::string coverageFileName, bool writeResult){

	std::cout << std::endl<< "Optimizing feature coverage for the image." << std::endl;

	// Optimizing coverage at single scale	
	float sigma = 5.6569f;
	float sigma_min = sigma;
	float sigma_max = sigma;
	int sigma_intervals = 1;

	double intensity_weight = 0.2; //0.1; //0.4; //0.2;
	double objectness_weight = 1.0 - intensity_weight;

	StatisticsFilterType::Pointer stats_filter = StatisticsFilterType::New();
	stats_filter->SetInput(this->PaddedCurvImage);
	stats_filter->Update();
	double img_max_val = stats_filter->GetMaximum();
	

	ObjectnessMeasures obj_measures(sigma_min, sigma_max, sigma_intervals, 0); // use this for astrocytes
	//ObjectnessMeasures obj_measures(sigma_min, sigma_max, sigma_intervals, 1); // use this for microglia
	obj_measures.alpha = 0.5 * img_max_val;
	obj_measures.beta = 0.5 * img_max_val;
	obj_measures.gamma = 0.25 * img_max_val; //0.25 * img_max_val;

	//std::cout << "Max image value: " << img_max_val << " Mean img value: " << img_mean_val << std::endl;

	this->ComputeObjectnessImage(obj_measures);

	//Write out the vesselenss image
	if(writeResult){
		std::string vesselnessPointsFileName = coverageFileName;
		vesselnessPointsFileName.erase(vesselnessPointsFileName.length()-4, vesselnessPointsFileName.length());
		vesselnessPointsFileName.append("_vesselness.mhd");

		//std::cout << vesselnessPointsFileName << std::endl;

		typedef itk::ImageFileWriter<ImageType3D> ImageWriterType;
		ImageWriterType::Pointer image_writer = ImageWriterType::New();
		image_writer->SetFileName(vesselnessPointsFileName);
		image_writer->SetInput(this->ObjectnessImage);
		image_writer->Update();
	}

	StatisticsFilterType::Pointer stats_filter2 = StatisticsFilterType::New();
	stats_filter2->SetInput(this->ObjectnessImage);
	stats_filter2->Update();
	double max_objectness = stats_filter2->GetMaximum();

	MultiplyImageFilter::Pointer img_multiplier = MultiplyImageFilter::New();
	img_multiplier->SetInput(this->ObjectnessImage);
	img_multiplier->SetConstant2(1.0/(max_objectness));
	img_multiplier->Update();

	MultiplyImageFilter::Pointer img_multiplier2 = MultiplyImageFilter::New();
	img_multiplier2->SetInput(this->PaddedCurvImage);
	img_multiplier2->SetConstant2(1.0/(img_max_val));
	img_multiplier2->Update();

	ImageType3D::Pointer normalized_img = img_multiplier2->GetOutput();

	ImageType3D::Pointer normalized_obj_img = img_multiplier->GetOutput();

	StatisticsFilterType::Pointer stats_filter3 = StatisticsFilterType::New();
	stats_filter3->SetInput(normalized_obj_img);
	stats_filter3->Update();
	double mean_objectness = stats_filter3->GetMean();
	double std_objectness = stats_filter3->GetVariance();

	StatisticsFilterType::Pointer stats_filter4 = StatisticsFilterType::New();
	stats_filter4->SetInput(normalized_img);
	stats_filter4->Update();
	double img_mean_val = stats_filter4->GetMean();
	double img_std_val = stats_filter4->GetVariance();
	

	MultiplyImageFilter::Pointer img_multiplier3 = MultiplyImageFilter::New();
	img_multiplier3->SetInput1(normalized_img);
	img_multiplier3->SetInput2(normalized_obj_img);
	img_multiplier3->Update();
	ImageType3D::Pointer int_obj_prod_img = img_multiplier3->GetOutput();

	MultiplyImageFilter::Pointer img_multiplier4 = MultiplyImageFilter::New();
	img_multiplier4->SetInput(normalized_img);
	img_multiplier4->SetConstant2(intensity_weight);
	img_multiplier4->Update();
	MultiplyImageFilter::Pointer img_multiplier5 = MultiplyImageFilter::New();
	img_multiplier5->SetInput(normalized_obj_img);
	img_multiplier5->SetConstant2(objectness_weight);
	img_multiplier5->Update();
	AddImageFilter::Pointer img_adder = AddImageFilter::New();
	img_adder->SetInput1(img_multiplier4->GetOutput());
	img_adder->SetInput2(img_multiplier5->GetOutput());
	img_adder->Update();
	ImageType3D::Pointer int_obj_sum_img = img_adder->GetOutput();
	
	StatisticsFilterType::Pointer stats_filter5 = StatisticsFilterType::New();
	stats_filter5->SetInput(int_obj_sum_img);
	stats_filter5->Update();
	double added_img_mean_val = stats_filter5->GetMean();
	
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
	
	int opt_iter = 0, max_opt_iter = 10; //20;
	float max_coverage = 0.002; //0.01; //0.0;
	float min_coverage = 0.001; // 0.0;
	float coverage_upper_limit = max_coverage + (0.2*max_coverage); //(0.05*max_coverage);
	float coverage_lower_limit = min_coverage - (0.25*min_coverage); //(0.2*min_coverage); //(0.05*min_coverage);

	std::cout << "Coverage limits: [" << coverage_lower_limit << ", " << coverage_upper_limit << "] " << std::endl;

	float thresh1_step_size = 0.02;
	float thresh2_step_size = 0.0004;

	float thresh1 = 0.01; //0.03; //0.01; //0.05; //0.08; //0.005; //0.03   // 3% of maximum theshold from Lowe 2004
	float thresh2 = 0.0001; //0.0005; //0.0001; //0.0009; //0.015; //0.0003;  //0.001 -0.1 percent of range

	double Mi_coverage4 = 0.0;
	double Mi_coverage5 = 0.0;
	double Mi_coverage6 = 0.0;


	std::ofstream optim_internal_file;
	
	if(writeResult){
		std::string internalOptimFileName = coverageFileName;
		internalOptimFileName.erase(internalOptimFileName.length()-4, internalOptimFileName.length());
		internalOptimFileName.append("_internal.txt");
		optim_internal_file.open(internalOptimFileName.c_str(), std::ios::out);
	}
			
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
		float total_fg_wt_obj_int = 0.0;
		
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

			
						total_fg_objectness += normalized_obj_img->GetPixel(ndx); //this->ObjectnessImage->GetPixel(ndx);
						total_fg_intensity += normalized_img->GetPixel(ndx);  //this->PaddedCurvImage->GetPixel(ndx);
						//total_fg_obj_int += this->ObjectnessImage->GetPixel(ndx) * this->PaddedCurvImage->GetPixel(ndx);
						total_fg_obj_int += normalized_obj_img->GetPixel(ndx) * this->PaddedCurvImage->GetPixel(ndx);
						total_fg_wt_obj_int += int_obj_sum_img->GetPixel(ndx);

					}	
					else
						rejectCtCnt++;
				}
			}
			++it;
			++nit;
		}

		double Mi_coverage = (double)ctCnt / img_mean_val;
		double Mi_coverage2 = (double)ctCnt / mean_objectness;
		double Mi_coverage3 = (double)ctCnt / (img_mean_val * mean_objectness);
		
		ImageType3D::SizeType im_size = this->PaddedCurvImage->GetBufferedRegion().GetSize();
		double total_objectness = mean_objectness * (im_size[0]*im_size[1]*im_size[2]);
		double total_intensity = img_mean_val * (im_size[0]*im_size[1]*im_size[2]);
		double total_int_obj_sum = added_img_mean_val * (im_size[0]*im_size[1]*im_size[2]);

		if(added_img_mean_val < 0.0001)
			total_int_obj_sum = 0.001 * (im_size[0]*im_size[1]*im_size[2]);
	

		Mi_coverage4 = total_fg_objectness / total_objectness;
		Mi_coverage5 = total_fg_intensity / total_intensity;
		Mi_coverage6 = total_fg_wt_obj_int / total_int_obj_sum;

		std::cout << "num: " << total_fg_wt_obj_int << " den: " << added_img_mean_val << std::endl;
		
		std::cout << std::endl;
		std::cout << "Number of CTs at this stage: " << ctCnt <<std::endl;
		std::cout << "Number of CTs rejected by RegisterIndex() are: " << rejectCtCnt << std::endl;
		
		//std::cout << "Total foreground objectness: " << total_fg_objectness << std::endl;
		//std::cout << "Total foreground intensity: " << total_fg_intensity << std::endl;
		//std::cout << "Total foreground obj*intensity: " << total_fg_obj_int << std::endl;
		//std::cout << "Mi_coverage: " << Mi_coverage << std::endl;
		//std::cout << "Mi_coverage2: " << Mi_coverage2 << std::endl;
		//std::cout << "Mi_coverage3: " << Mi_coverage3 << std::endl;
		//std::cout << "Mi_coverage4: " << Mi_coverage4 << std::endl;
		//std::cout << "Mi_coverage5: " << Mi_coverage5 << std::endl;
		//std::cout << "Mi_coverage6: " << Mi_coverage6 << std::endl;

		std::cout << "Iter: " << opt_iter << " Mi_coverage6: " << Mi_coverage6 << " Intensity threshold: " << thresh1 << ", Contrast threshold: " << thresh2 << std::endl;

		if(writeResult){
			if(optim_internal_file.good())
				optim_internal_file << Mi_coverage6 << '\t' << thresh1 << '\t' << thresh2 << '\t' << ctCnt << std::endl;
			else
				std::cout << "Error writing coverage internal file. " << std::endl;
		}

		opt_iter++;
		

		if(Mi_coverage6 > coverage_upper_limit){
			thresh1 += thresh1_step_size;
			thresh2 += thresh2_step_size;
		}
		else if(Mi_coverage6 < coverage_upper_limit && Mi_coverage6 > max_coverage){
			thresh1 += thresh1_step_size / 5.0;
			thresh2 += thresh2_step_size / 5.0;
		}
		else if(Mi_coverage6 < coverage_lower_limit){

			if(thresh1 <= 0.01){
				thresh1 -= thresh1_step_size / 5.0;
				thresh2 -= thresh2_step_size / 5.0;
			}
			else if(thresh1 <= 0.0){
				thresh1 = 0.0;
				thresh2 = 0.0;
				break;
			}
			else{
				thresh1 -= thresh1_step_size;
				thresh2 -= thresh2_step_size;
			}
		}
		else
			break;		
	}

	if(writeResult)
		optim_internal_file.close();

	// Setting thresholds based on the optimized coverage
	this->intensity_threshold = thresh1;
	this->contrast_threshold = thresh2;

	// Printing the results
	if(opt_iter == max_opt_iter){
		std::cout << "Coverage might not be correctly optimized. Manual adjustments might be needed. " << std::endl;
		this->isCoverageOptimized = false;
	}
	else if(thresh1 = 0.0){
		std::cout << "Coverage might not be correctly optimized. Manual adjustments might be needed. " << std::endl;
		this->isCoverageOptimized = false;
	}
	else{
		std::cout << "Done with optimizing coverage. " <<  std::endl;
		this->isCoverageOptimized = true;
	}

	
	if(writeResult){

		std::ofstream coverage_file;
		coverage_file.open(coverageFileName.c_str(), std::ios::out);
		if(coverage_file.good()){
			coverage_file << Mi_coverage6 << std::endl << this->intensity_threshold << std::endl << this->contrast_threshold << std::endl;
			coverage_file << img_mean_val << std::endl << img_std_val << std::endl << mean_objectness << std::endl << std_objectness << std::endl;

			coverage_file.close();
		}
		else
			std::cout << "Error writing coverage file. " << std::endl;
		

		std::string LoGPointsFileName = coverageFileName;
		LoGPointsFileName.erase(LoGPointsFileName.length()-4, LoGPointsFileName.length());
		LoGPointsFileName.append("_optimum_LOG_points.tif");
		
		RescalerType::Pointer rescaler2 = RescalerType::New();
		rescaler2->SetInput(NDXImage);
		rescaler2->SetOutputMaximum( 255 );
		rescaler2->SetOutputMinimum( 0 );
		rescaler2->Update();
		itk::CastImageFilter< ImageType3D, CharImageType3D>::Pointer caster2 = itk::CastImageFilter< ImageType3D, CharImageType3D>::New();
		caster2->SetInput(rescaler2->GetOutput());

		itk::ImageFileWriter< CharImageType3D >::Pointer LoGwriter2 = itk::ImageFileWriter< CharImageType3D >::New();

		LoGwriter2->SetFileName(LoGPointsFileName.c_str());	
		LoGwriter2->SetInput(caster2->GetOutput());
		LoGwriter2->Update();			
	}
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

bool HeapNode_astro::operator ==(const HeapNode_astro& h2){

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
	
	std::vector<HeapNode_astro> current_LoG_points;
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


		const float thresh1 = this->intensity_threshold; //0.05; //0.08; //0.005; //0.03   // 3% of maximum theshold from Lowe 2004
		const float thresh2 = this->contrast_threshold;  //0.0009; //0.015; //0.0003;  //0.001 -0.1 percent of range
		
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
				// What is RegisterIndex used for? For keeping the strongest ridge point across scales
				// value corresponds to strength of being a ridge point?
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
					current_LoG_points.push_back(HeapNode_astro(ndx, value));	

					/*HeapNode_astro current_node(ndx, value);
					
					if(!this->AllLoGPointsVector.empty()){ 
						
						std::vector<HeapNode_astro>::iterator vec_it = std::find(this->AllLoGPointsVector.begin(), this->AllLoGPointsVector.end(), current_node);
						if(vec_it != this->AllLoGPointsVector.end())
							this->AllLoGPointsVector.erase(vec_it);
					}*/

					this->AllLoGPointsVector.push_back(HeapNode_astro(ndx, value));
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
		/*RescalerType::Pointer rescaler2 = RescalerType::New();
		rescaler2->SetInput(LoGCurrentScaleImage);
		rescaler2->SetOutputMaximum( 255 );
		rescaler2->SetOutputMinimum( 0 );
		rescaler2->Update();
		itk::CastImageFilter< ImageType3D, CharImageType3D>::Pointer caster2 = itk::CastImageFilter< ImageType3D, CharImageType3D>::New();
		caster2->SetInput(rescaler2->GetOutput());

		itk::ImageFileWriter< CharImageType3D >::Pointer LoGwriter2 = itk::ImageFileWriter< CharImageType3D >::New();

		if(scale_index == 1)
			LoGwriter2->SetFileName("C:\\Prathamesh\\Astrocytes\\ShearletExp\\Control\\LoG_Points_1.tif");	
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
		*/	
	}
}

void AstroTracer::GetFeature( float sigma ) 
{
	std::cout << "Get Feature 1 " << std::endl;
	
	clock_t LoG_start_time = clock();
	typedef itk::LaplacianRecursiveGaussianImageFilter< ImageType3D , ImageType3D> GFilterType;
	GFilterType::Pointer gauss = GFilterType::New();
	gauss->SetInput( PaddedCurvImage );
	gauss->SetSigma( sigma );
	gauss->SetNormalizeAcrossScale(false);
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
	std::cout << "Number of CTs at this stage: " << ctCnt << std::endl << std::flush;
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
SWCNode_astro* AstroTracer::TBack(itk::Index<3> &ndx, std::vector<IndexType>& Chain)  
{
	// Label is either the top of the chain (a primary node) or NULL
	SWCNode_astro* Label = NULL;

	itk::Index<3> n;
	itk::Vector<float,3> p, x, d, dold;
	for (int i=0; i<3; i++) 
	{
		p[i] = static_cast<PixelType>(ndx[i]);
		dold[i] = 0.0f;
	}
	bool done = false;

	// If this is a primary node, do not trace back, return NULL.
	if (SWCImage->GetPixel(ndx)->TreeID > 0) 
	{
		done = true;
	}
	const float MAXDERV = 10000.0f;

	// Pushed back the cuurent index, ndx, the starting point of the chain.
	Chain.push_back(ndx);

	// Continue tracing in the direction of the derivative of the tracing cost (from ConnImage) until you hit a primary/labeled node or until you are lost!
	while (done == false) 
	{
		// Find the derivative of the ConnImage (cost of tracing) in the three directions
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
		
		// Push back next element index to the chain
		Chain.push_back(n);

		//check termination
		SWCNode_astro *t = SWCImage->GetPixel(n);
		
		// Check if the chain has hit either a valid node (primary or secondary node)
		if (t != NULL ) 
		{
			// If it has hit a primary/labeled node, stop the chain here
			if (t->TreeID > 0) 
			{
				done = true;
				Label = SWCImage->GetPixel(n);
				break;
			}
		}

		// If the chain is growing to be so big, you are misleaded!!
		if (Chain.size() > 500) 
		{
			//std::cout << "Tree not found for " << ndx << " in 500 steps, exiting!! " << std::endl;
			Chain.clear();
			Label = NULL;
			done = true;
			break;
		}
	}
	

	float costFactorLabel = 0.0;
	
	if(Chain.size() != 0 ){

		//this functions is used in TBack to prevent nodes to join trees, in cases of abrupt directional changes 
		costFactorLabel = GetCostLocalLabel(Label , ndx);
	}

	// If the curvature-based cost if too high, then this chain is not added to the tree.
	// This is a good place to add additional contraints for astrocytes?
	if(costFactorLabel >= 10.0){	
		Chain.clear();
		done = true;
		return NULL;
	}
	else 
		return Label;


	//return Label;
}

///////////////////////////////////////////////////////////////////////////////////
float AstroTracer::GetCost(SWCNode_astro* s, itk::Index<3>& endx ) 
{

	// Similar to GetCostLocal, but accumulates the angles along the trace at regular intervals, and returns this as the cost.
	// If the accumulated angle is high, this would indicate abrupt turns.

	// Try to use the std(curvature) in a similar way for cost

	// Try adding vesselness to the cost function

	itk::Index<3> base = endx, ndx = s->ndx;
	float cost = 0.0f, angsum = 0.0f, count = 0.01f, ignore_dist = 0.0;
	itk::Vector<float,3> d1, d2 , gd, gd1;
	d1.Filled(0.0);
	gd1.Fill(0.0);
	bool first = true;
	int max_points = 500; //1000; //500
	std::vector<float> angle_vec, del_scale_vec, del_likelihood_vec, del_vesselness_vec;
	float acc_likelihood = 0.0, acc_del_scale = 0.0, acc_vesselness = 0.0;
	float max_trace_length = 150.0, trace_length = 0.0, trace_length_norm = 0.0;


	while (count < max_points) 
	{
		
		// Length based cost (this is a hard cost, consider giving low weight to this)
		trace_length++;

		if(first)
			ignore_dist = 50.0f;
		else
			ignore_dist = 3.0; //10; //20.0f; //10.0f; // 10 is good for curvature cost

		float d = (ndx[0] - base[0])*(ndx[0] - base[0]) + (ndx[1] - base[1])*(ndx[1] - base[1]) + (ndx[2] - base[2])*(ndx[2] - base[2]) ;
		if ( vcl_sqrt(d) > ignore_dist) //6.0f) 
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
				//count++;

				//del_vesselness_vec.push_back(this->ObjectnessHybridImage->GetPixel(base));
				//itk::Vector<float, 3> pos_base;
				//pos_base[0] = base[0]; pos_base[1] = base[1]; pos_base[2] = base[2];
				//float likelihood_base = 0.0;
				//float scale_base = getRadiusAndLikelihood(pos_base, likelihood_base);
				//del_likelihood_vec.push_back(likelihood_base);
			}
			else 
			{
				// Scale gradient based cost
				itk::Vector<float, 3> pos_ndx, pos_base;
				pos_ndx[0] = ndx[0]; pos_ndx[1] = ndx[1]; pos_ndx[2] = ndx[2];
				pos_base[0] = base[0]; pos_base[1] = base[1]; pos_base[2] = base[2];
				float likelihood_ndx = 0.0, likelihood_base = 0.0;
				float scale_ndx = getRadiusAndLikelihood(pos_ndx, likelihood_ndx);
				float scale_base = getRadiusAndLikelihood(pos_base, likelihood_base);

				float del_scale = scale_ndx - scale_base;
				float del_likelihood = likelihood_ndx * likelihood_base;

				//std::cout << del_scale << ", " << del_likelihood << std::endl;
				
				acc_likelihood += (likelihood_ndx * likelihood_base);
				acc_del_scale += del_scale;
				
				del_scale_vec.push_back(del_scale);
				del_likelihood_vec.push_back(del_likelihood);

				// Vesselness-based cost
				float del_vesselness = this->ObjectnessHybridImage->GetPixel(ndx) * this->ObjectnessHybridImage->GetPixel(base);
				del_vesselness_vec.push_back(del_vesselness);				

				acc_vesselness += del_vesselness;
				
				// Curvature-based cost
				PixelType w = dot_product(d1.Get_vnl_vector(),d2.Get_vnl_vector());
				if (w < 0.99f) 
				{
					angsum += vcl_acos(vnl_math_abs(w));
					angle_vec.push_back(vcl_acos(vnl_math_abs(w)));
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

	if(angle_vec.empty())
		return 0;
	if(angle_vec.size() < ignore_dist)
		return 0;
	
	float mean_curvature = std::accumulate(angle_vec.begin(), angle_vec.end(), 0.0);
	mean_curvature = mean_curvature / (float)angle_vec.size();
	
	float std_curvature = 0.0;
	for(int i = 0; i < angle_vec.size(); i++)
		std_curvature = std_curvature + std::pow(angle_vec[i] - mean_curvature, 2);
	std_curvature = std::sqrt(std_curvature / (float)angle_vec.size());

	
	float mean_del_scale = std::accumulate(del_scale_vec.begin(), del_scale_vec.end(), 0.0);
	mean_del_scale = mean_del_scale / (float)del_scale_vec.size();
	
	float std_del_scale = 0.0;
	for(int i = 0; i < del_scale_vec.size(); i++)
		std_del_scale = std_del_scale + std::pow(del_scale_vec[i] - mean_del_scale, 2);
	std_del_scale = std::sqrt(std_del_scale / (float)del_scale_vec.size());
	

	float mean_del_likelihood = std::accumulate(del_likelihood_vec.begin(), del_likelihood_vec.end(), 0.0);
	mean_del_likelihood = mean_del_likelihood / (float)del_likelihood_vec.size();
	
	float std_del_likelihood = 0.0;
	for(int i = 0; i < del_likelihood_vec.size(); i++)
		std_del_likelihood = std_del_likelihood + std::pow(del_likelihood_vec[i] - mean_del_likelihood, 2);
	std_del_likelihood = std::sqrt(std_del_likelihood / (float)del_likelihood_vec.size());


	float mean_vesselness = std::accumulate(del_vesselness_vec.begin(), del_vesselness_vec.end(), 0.0);
	mean_vesselness = mean_vesselness / (float)del_vesselness_vec.size();
	
	float std_vesselness = 0.0;
	for(int i = 0; i < del_vesselness_vec.size(); i++)
		std_vesselness = std_vesselness + std::pow(del_vesselness_vec[i] - mean_vesselness, 2);
	std_vesselness = std::sqrt(std_vesselness / (float)del_vesselness_vec.size());

	
	//trace_length = count; 

	//std::cout << std_curvature << ", " << std_del_scale << ", " << std_del_likelihood << std::endl;
	//std::cout << mean_del_scale << ", " << mean_del_likelihood << std::endl;
	//std::cout << mean_vesselness << ", " << acc_vesselness << ", " << std_vesselness << std::endl;
	//std::cout << "Trace length: " << trace_length/max_trace_length << std::endl;

	gd[0] = float(ndx[0] - endx[0]);
	gd[1] = float(ndx[1] - endx[1]);
	gd[2] = float(ndx[2] - endx[2]);
	gd.Normalize();

	float allowedTurns = 1.0f; //1.0f; //4.0f;

	// If the trace has taken a U-turn, the dot product will be negetive 
	// Allowing U-turns for astrocytes
	if(true) //(dot_product(gd.Get_vnl_vector(),gd1.Get_vnl_vector() ) >= 0) 
	{
		//cost = (angsum/allowedTurns);
		//cost = std_curvature;
		//cost = mean_curvature;
		//cost = std_del_scale;

		if(mean_del_scale > 1.0)
			mean_del_scale = 1.0;
		if(mean_del_scale < 0.0)
			mean_del_scale = 0.0;

		//cost = 1.0 - mean_del_scale;
		
		if(mean_del_likelihood > 1.0)
			mean_del_likelihood = 1.0;
		if(mean_del_likelihood < 0.0)
			mean_del_likelihood = 0.0;

		//cost = 1.0 - mean_del_likelihood;

		if(acc_del_scale > 1.0)
			acc_del_scale = 1.0;
		if(acc_del_scale < 0.0)
			acc_del_scale = 0.0;
		
		//cost = 1.0 - acc_del_scale;

		if(acc_likelihood > 1.0)
			acc_likelihood = 1.0;
		if(acc_likelihood < 0.0)
			acc_likelihood = 0.0;

		//cost = 1.0 - acc_likelihood;

		if(acc_vesselness > 1.0)
			acc_vesselness = 1.0;
		if(acc_vesselness < 0.0)
			acc_vesselness = 0.0;

		//cost = 1.0 - acc_vesselness;

		if(mean_vesselness > 1.0)
			mean_vesselness = 1.0;
		if(mean_vesselness < 0.0)
			mean_vesselness = 0.0;

		//cost = 1.0 - mean_vesselness;

		trace_length_norm = trace_length / max_trace_length;
		if(trace_length_norm > 1.0)
			trace_length_norm = 1.0;

		//cost = std_vesselness;


		//cost = 10.0 * std_del_likelihood;
		//cost = mean_del_likelihood;

		//cost = 1.0 - acc_likelihood; // try this with scale also
		

		// Combine the costs into "mother of all costs"
		float scale_based_cost = 1.0 - mean_del_scale;
		float likelihood_based_cost = 1.0 - mean_del_likelihood;
		float vesselness_based_cost = 1.0 - mean_vesselness;
		float length_based_cost = trace_length_norm;

		// Equally weighing all the costs right now
		cost = (scale_based_cost + likelihood_based_cost + vesselness_based_cost + length_based_cost)/4.0;

		//std::cout << "Combined cost: " << cost << std::endl;

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

float AstroTracer::GetCostOld(SWCNode_astro* s, itk::Index<3>& endx){

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

float AstroTracer::GetCostLocal2(SWCNode_astro* s, itk::Index<3>& endx){

	// Returns a cost between 0 and 1, based on standard deviation of the curvature in order to avoid abrupt changes in the traces.
	// Curvature is computed as the cosine of successive nodes starting from the current node upto its root.

	itk::Index<3> new_leaf_ndx = endx, leaf_ndx = s->ndx;
	float cost = 1.0f, cost_low_bound = 1.0f, cost_high_bound = 2.0f;
	float current_curvature = 0.0f;
	int ignore_N_nodes = 6.0, most_N_nodes = 1000;
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
		
		current_curvature = vcl_acos(vcl_abs(dot_product(pos_diff.Get_vnl_vector(), pos_diff_last.Get_vnl_vector())));
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
	
	cost = std_curvature; //+ cost_low_bound; // Minimum cost is cost_low_bound when std is zero

	if(cost > 1.0)
		cost = 1.0;
			
	return cost;
}

float AstroTracer::GetCostLocal3(SWCNode_astro* s, itk::Index<3>& endx){


	// Check modifications to the cost function based on gradient of scale, likelihood and intensity

	itk::Index<3> new_leaf_ndx = endx, leaf_ndx = s->ndx;
	float cost = 1.0f, cost_low_bound = 1.0f, cost_high_bound = 2.0f;
	float current_curvature = 0.0f;
	int ignore_N_nodes = 6.0, most_N_nodes = 10000;
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

	return cost;
}
/////////////////////////////////////////////////////////////////////////////////////
float AstroTracer::GetCostLocal(SWCNode_astro* s, itk::Index<3>& endx ) 
{
	// Same as GetCostLocalLabel but does not flag abrupt turns

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

float AstroTracer::GetCostLocalLabel(SWCNode_astro* s, itk::Index<3>& endx ) 
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

		// Make sure you have reached far enough along the parents of the top of the chain to avoid being too local
		if ( vcl_sqrt(d) > 6.0f) //6.0 (0)
		{
			d2[0] = float(ndx[0] - base[0]); //  ancestor-leaf
			d2[1] = float(ndx[1] - base[1]);
			d2[2] = float(ndx[2] - base[2]);
			d2.Normalize();

			PixelType w = dot_product(d1.Get_vnl_vector(),d2.Get_vnl_vector());
			
			// Cost depends on the cosine of the angle between 2 vectors
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

			// The cost is computed only once when you have reached far enough along the parents of the label (i.e s). This could be computed as an accumulated
			// cost as one moves away from the top of the chain (label).
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
	// Sort a1, a2, a3 in ascending order
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
	//std::cout << "Decimating the tree of size: " << SWCNode_astroContainer.size() << std::endl;
	std::vector<SWCNode_astro*>::iterator sit;
	for (sit = SWCNode_astroContainer.begin(); sit != SWCNode_astroContainer.end(); ++sit) 
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
	for (sit = SWCNode_astroContainer.begin(); sit != SWCNode_astroContainer.end(); ++sit) 
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

	for (sit = SWCNode_astroContainer.begin(); sit != SWCNode_astroContainer.end(); ++sit) 
	{
		if ((*sit)->IsLeaf == true) 
		{
			//std::cout << "Leaf at" << (*sit)->ndx << " ID:" << (*sit)->ID << " ParentID:" << (*sit)->PID << std::endl;
			SWCNode_astro* par = (*sit)->parent;
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
				SWCNode_astro* n = (*sit);
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

	std::vector<SWCNode_astro*> NewContainer;
	NewContainer.reserve(SWCNode_astroContainer.size());

	long newID = 1;
	itk::Array<long> IDLookUp(SWCNode_astroContainer.size());
	IDLookUp.Fill(0);

	for (unsigned int i=0; i < SWCNode_astroContainer.size(); ++i) 
	{
		if (SWCNode_astroContainer[i]->IsActive == true) 
		{
			IDLookUp[i] = newID++;		
		}
	}
	//std::cout << "Lookup generated: " << std::endl;

	for (sit = SWCNode_astroContainer.begin(); sit != SWCNode_astroContainer.end(); ++sit) 
	{
		if ((*sit)->IsActive == true) 
		{
			long PID;
			itk::Index<3> ndx;
			for (int i = 0; i < 3; ++i) 
			{
				ndx[i] = long((*sit)->pos[i]);
			}

			SWCNode_astro* par = (*sit)->parent;
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

			SWCNode_astro* s = new SWCNode_astro();
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

	for (sit = SWCNode_astroContainer.begin(); sit != SWCNode_astroContainer.end(); ++sit) 
	{
		delete (*sit);
	}

	SWCNode_astroContainer = NewContainer;
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

	std::vector<SWCNode_astro*>::iterator sit;
	for (sit = SWCNode_astroContainer.begin(); sit != SWCNode_astroContainer.end(); ++sit) 
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

	unsigned int originalSize = SWCNode_astroContainer.size();
	LabelArrayType somaArray = SomaImage->GetBufferPointer();
	itk::Size<3> im_size = SomaImage->GetBufferedRegion().GetSize();
	int slice_size = im_size[0] * im_size[1];
	int row_size = im_size[0];

	//find the root nodes of each tree
	std::cout << "Finding the root nodes of each tree" << std::endl;
	std::map<long, SWCNode_astro*> treeIDToRootMap;
	std::vector<SWCNode_astro*>::iterator sit;
	for (sit = SWCNode_astroContainer.begin(); sit != SWCNode_astroContainer.end(); ++sit)
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
	for (sit = SWCNode_astroContainer.begin(); sit != SWCNode_astroContainer.end();)
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
			sit = SWCNode_astroContainer.erase(sit);
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
				sit = SWCNode_astroContainer.erase(sit);
				//std::cout << "Deleted node. " << std::endl;
			}
			else{
				SWCNode_astro *parent = (*sit)->parent;
				SWCNode_astro *root = treeIDToRootMap[(*sit)->TreeID];

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
			SWCNode_astro *parent = (*sit)->parent;
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

	size_t newSize = SWCNode_astroContainer.size();
	std::cout << "Just removed " << originalSize - newSize
		<< " nodes (" << originalSize << " to " << newSize << ")"
		<< std::endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////
void AstroTracer::WriteMultipleSWCFiles(std::string fname, unsigned int padz) 
{
	// check number of start points to determine number of files to write, with new filename eachtime
	std::cout << "Total " << SWCNode_astroContainer.size() << " nodes..." <<std::endl;
	std::vector<SWCNode_astro*>::iterator sit;
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
		for (sit = SWCNode_astroContainer.begin(); sit != SWCNode_astroContainer.end(); ++sit) 
		{
			if ((*sit)->TreeID == i+1) 
				NodeIDToSWCIDMap[(*sit)->ID] = ID++;			
		}
		std::cout << ID << " Nodes found  ";

		//create the SWCImage file
		for (sit = SWCNode_astroContainer.begin(); sit != SWCNode_astroContainer.end(); ++sit) 
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

	for (sit = SWCNode_astroContainer.begin(); sit != SWCNode_astroContainer.end(); ++sit) 
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
					if (this->PaddedCurvImage->GetBufferedRegion().IsInside(ndx)) 
					{
						float val = this->PaddedCurvImage->GetPixel(ndx);
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
	std::vector<SWCNode_astro*>::iterator sit;
	std::cout << "Writing SWCImage file " << fname << " with " << SWCNode_astroContainer.size() << " nodes...";
	std::ofstream ofile(fname.c_str());
	//ofile << "#Neuron Tracing Code 3D, RPI" << std::endl;
	//ofile << "#author: AM" << std::endl;


	for (sit = SWCNode_astroContainer.begin(); sit != SWCNode_astroContainer.end(); ++sit) 
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
SWCNode_astro::SWCNode_astro()
{
	ID = -1;
	PID = -1;
	TreeID = -1;
	IsLeaf = false;
	IsBranch = false;
	parent = NULL;
	children.reserve(2);
}

SWCNode_astro::SWCNode_astro(long id, long pid, long tid, itk::Index<3> n)
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

SWCNode_astro::SWCNode_astro(long id, SWCNode_astro * p, long tid, itk::Index<3> n)
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

HeapNode_astro::HeapNode_astro(itk::Index<3> n1, PixelType d)
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
	multi_scale_Hessian->SetSigmaMinimum(obj_measures.sigma_min);
	multi_scale_Hessian->SetSigmaMaximum(obj_measures.sigma_max);
	multi_scale_Hessian->SetNumberOfSigmaSteps(obj_measures.sigma_intervals);

	ObjectnessFilterType::Pointer objectness_filter = ObjectnessFilterType::New();
	//ObjectnessFilterType::Pointer objectness_filter = multi_scale_Hessian->GetHessianToMeasureFilter();
	
	objectness_filter->SetScaleObjectnessMeasure(false);
	objectness_filter->SetBrightObject(true);
	objectness_filter->SetAlpha(obj_measures.alpha);
	objectness_filter->SetBeta(obj_measures.beta);
	objectness_filter->SetGamma(obj_measures.gamma);
	objectness_filter->SetObjectDimension(obj_measures.objectness_type);
    
    multi_scale_Hessian->SetHessianToMeasureFilter(objectness_filter);
	
	//std::cout << obj_measures.alpha << std::endl << obj_measures.beta << std::endl << obj_measures.gamma << std::endl;

	try
    {
        multi_scale_Hessian->Update();
	}
    catch (itk::ExceptionObject &err)
    {
        std::cerr << "Error in multiscale Hessian filter" << std::endl;
    }
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

	// This function is depricated - use ComputeAstroFeaturesPipeline
	std::cout << "This function is depricated. Use ComputeAstroFeaturesPipeline. " << std::endl;

	//Prepare a list of points at which to compute the features
	std::string option1 = "SWC";
	std::string option2 = "LOG";
	std::string option3 = "External";

	std::vector<HeapNode_astro> points_list;
	if(!strcmp(points_source.c_str(), option1.c_str())){
		std::cout << "Computing features around SWC points. " << std::endl;

		std::vector<SWCNode_astro*>::iterator sit;
		for(sit = SWCNode_astroContainer.begin(); sit != SWCNode_astroContainer.end(); ++sit){
			itk::Index<3> idx; 
			idx[0] = (*sit)->pos[0]; idx[1] = (*sit)->pos[1]; idx[2] = (*sit)->pos[2];
			points_list.push_back(HeapNode_astro(idx, 0));
		}
	}
	else if(!strcmp(points_source.c_str(), option2.c_str())){
		std::cout << "Computing features around LOG points. " << this->LoGPointsVector.size() << std::endl;

		if(this->AllLoGPointsVector.empty()){
			std::cout << "NO LOG POINTS FOUND. RETURNING. " << std::endl;
			return;
		}

		for(int i = 0; i < this->AllLoGPointsVector.size(); i++){
			//std::vector<HeapNode_astro> currentLoGPoints = this->LoGPointsVector[i];
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

	std::ofstream feature_vector;
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

		double radiusThresh = 2.0;
		if(radius < radiusThresh)
			continue;

		ImageType3D::IndexType current_idx;
		current_idx[0] = pos[0];
		current_idx[1] = pos[1];
		current_idx[2] = pos[2]-padz; // PADZ???
		//std::cout << "After current_ide set, LoG_vector_size="<< this->LoG_Vector.size() << std::endl;
		//std::cout << "Current Indices:"<<current_idx[0]<<" "<<current_idx[1]<<" "<<current_idx[2]<<" "<<std::endl;


		////////////////// Code for nearest nuclei /////////////////////////
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
				a_root.featureVector.node = HeapNode_astro(points_list[i]);
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



void AstroTracer::Set_DistanceMapImage(AstroTracer::ImageType3D::Pointer distance_map_image){
	this->SomaDistanceMapImage = distance_map_image;
}

void AstroTracer::Set_NucleiLabelImage(AstroTracer::LabelImageType3D::Pointer nuc_label_image){
	this->NucleiLabelImage = nuc_label_image;
}

void AstroTracer::ComputeSomaDistanceMap(){
	
	typedef itk::BinaryThresholdImageFilter<LabelImageType3D, CharImageType3D> ThresholdFilterType;
	ThresholdFilterType::Pointer threshold_filter = ThresholdFilterType::New();
	threshold_filter->SetLowerThreshold(1);
	threshold_filter->SetInsideValue(255);
	threshold_filter->SetOutsideValue(0);
	threshold_filter->SetInput(this->SomaImage);
	
	SignedMaurerDistanceMapImageFilterType::Pointer MaurerFilter = SignedMaurerDistanceMapImageFilterType::New();
	MaurerFilter->SetInput(threshold_filter->GetOutput());
	MaurerFilter->SetSquaredDistance(false);
	MaurerFilter->SetUseImageSpacing(false);
	MaurerFilter->SetInsideIsPositive(false);
	MaurerFilter->Update();
	
	this->SomaDistanceMapImage = MaurerFilter->GetOutput();
}

void AstroTracer::ComputeRootPointFeatures(){
	
	std::vector< vtkSmartPointer< vtkTable > > feature_tables;
	feature_tables.resize(1);
		
	std::vector< LabelImageType3D::Pointer > output_images;
	output_images.resize(1);

	std::string featureVectorFileName = this->InputDataPath;
	featureVectorFileName.append("_feature_vector_roots.txt");

	std::string IDImageFileName = this->InputDataPath;
	IDImageFileName.append("_RootsImage.tif");
	
	this->ComputeSomaDistanceMap();
	this->ComputeAstroFeaturesPipeline(featureVectorFileName, IDImageFileName, 0, this->PaddedCurvImage->GetBufferedRegion(), feature_tables, output_images, true); 

	ftk::SaveTable(featureVectorFileName, feature_tables[0]);
}

VBTNode AstroTracer::getVBTFeatures(itk::Index<3> current_idx){

	VBTNode node;
	node.SetLocationFromITKIndex(current_idx);
	
	this->VBT->FitSphereAtVBTNode(node);
	this->VBT->ComputeODFNoPrior(node);
	this->VBT->ComputeODFFeatures(node);
	
	return node;
}

bool AstroTracer::getLocalIntensityFeatures(itk::Index<3> current_idx, itk::Size<3> sz, float radius, float& max_intensity, float& min_intensity, float& mean_intensity, float& var_intensity){

	float double_scale = 2.0 * radius; 
		
	//StatisticsFilterType::Pointer stats_filter = StatisticsFilterType::New();
	ImageType3D::IndexType starting_index, end_index;
	ImageType3D::SizeType sub_volume_size;
	ImageType3D::RegionType sub_volume_region;
	ImageType3D::Pointer sub_volume;

	starting_index[0] = current_idx[0] - double_scale; starting_index[1] = current_idx[1] - double_scale; starting_index[2] = current_idx[2] - double_scale;
	end_index[0] = current_idx[0] + double_scale; end_index[1] = current_idx[1] + double_scale; end_index[2] = current_idx[2] + double_scale;

	if ( (starting_index[0] < 0) || (starting_index[1] < 0) || (starting_index[2] < 0) ||
		(end_index[0] > (unsigned int)sz[0]) || (end_index[1] > (unsigned int)sz[1]) ||
		(end_index[2] > (unsigned int)sz[2]) )
		return false;

	sub_volume_size[0] = 2 * double_scale; sub_volume_size[1] = 2 * double_scale; sub_volume_size[2] = 2 * double_scale;
	
	sub_volume_region.SetIndex(starting_index);
	sub_volume_region.SetSize(sub_volume_size);
	
	VolumeOfInterestFilterType::Pointer sub_volume_filter = VolumeOfInterestFilterType::New();
	sub_volume_filter->SetInput(this->PaddedCurvImage);
	sub_volume_filter->SetRegionOfInterest(sub_volume_region);
	sub_volume_filter->Update();
	sub_volume = sub_volume_filter->GetOutput();
	
	StatisticsFilterType::Pointer stats_filter = StatisticsFilterType::New();
	stats_filter->SetInput(sub_volume);
	stats_filter->Update();

	mean_intensity = stats_filter->GetMean();
	var_intensity = stats_filter->GetVariance();
	max_intensity = stats_filter->GetMaximum();
	min_intensity = stats_filter->GetMinimum();

	return true;
}

void AstroTracer::getHessianEigenFeatures(itk::Index<3> current_idx, double img_max_val, float &ballness, float &plateness, float &vesselness, float &noiseness){

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

		vesselness = (1.0 - std::exp(-1 * std::pow(plateness, 2) / (2.0 * std::pow(alpha, 2)))) 
			* std::exp(-1.0 * std::pow(ballness, 2) / (2.0 * std::pow(beta, 2)))
			* (1.0 - std::exp(-1.0 * std::pow(noiseness, 2) / (2.0 * std::pow(gamma, 2))));
	}
}

void AstroTracer::ComputeAstroFeaturesPipeline(std::string outputFname, std::string IDFname, unsigned int padz, ImageType3D::RegionType regionLocal_inside, 
											   std::vector<vtkSmartPointer<vtkTable> >& features_table_vec, std::vector<LabelImageType3D::Pointer>& out_images, 
											   const bool writeResult){
	
	//Prepare a list of points at which to compute the features
	std::vector<HeapNode_astro> points_list;
	
	std::cout << "Computing features around LOG points for scales: " << this->LoGPointsVector.size() << std::endl;
	if(this->AllLoGPointsVector.empty()){
		std::cout << "NO LOG POINTS FOUND. RETURNING. " << std::endl;
		return;
	}
	if(this->SomaDistanceMapImage.IsNull()){
		std::cout << "Distance map image is NULL. Returning. " << std::endl;
		return;
	}

	for(int i = 0; i < this->AllLoGPointsVector.size(); i++){
		if(regionLocal_inside.IsInside(this->AllLoGPointsVector[i].ndx))
			points_list.push_back(this->AllLoGPointsVector[i]);
	}
	
	std::cout << "Computing features around N points. " << points_list.size() << std::endl;

	if(points_list.size() == 0){
		std::cout << "NO N POINTS FOUND. RETURNING. " << std::endl;
		return;
	}

	if(points_list.size() > 100000){
		std::cout << "\nTOO MANY LOG POINTS\n";
		std::cout << "TILE NO: " << outputFname << "\n";
		return;
	}

	//Preparing IDImage
	LabelImageType3D::RegionType id_reg;
	LabelImageType3D::IndexType id_st;
	LabelImageType3D::SizeType id_sz = this->PaddedCurvImage->GetBufferedRegion().GetSize();
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

	unsigned short IDIndex = 0;//ID index
	itk::Size<3> sz = this->PaddedCurvImage->GetBufferedRegion().GetSize();

	std::vector<std::string> col_names;
	col_names.push_back("ID");
	col_names.push_back("x");
	col_names.push_back("y");
	col_names.push_back("z");
	col_names.push_back("radius");
	col_names.push_back("sph_likelihood");
	col_names.push_back("ball");
	col_names.push_back("plate");
	col_names.push_back("vessel");
	col_names.push_back("int");
	col_names.push_back("m_int");
	col_names.push_back("v_int");
	col_names.push_back("mx_int");
	col_names.push_back("mn_int");
	col_names.push_back("nucl_dist");
	col_names.push_back("radius_VBT");
	col_names.push_back("likelihood_VBT");
	col_names.push_back("odf_std");
	col_names.push_back("odf_mean");
	col_names.push_back("odf_energy");
	col_names.push_back("odf_val1");
	col_names.push_back("odf_val2");
	col_names.push_back("odf_val3");
			
	this->features_table = vtkSmartPointer<vtkTable>::New();
	this->features_table->Initialize();

	for(int i = 0; i < col_names.size(); i++){

		vtkSmartPointer<vtkDoubleArray> table_col = vtkSmartPointer<vtkDoubleArray>::New();		
		table_col->SetName(col_names[i].c_str()); 
		
		this->features_table->AddColumn(table_col);
	}

	std::vector<CandidateRootPoint> root_points(points_list.size(), CandidateRootPoint());

#pragma omp parallel for
	for(int i = 0; i < points_list.size(); i++){

		// current_idx has to be local to the current tile (PaddedCurvImage)
		itk::Vector<float, 3> pos;
		pos[0] = points_list[i].ndx[0];
		pos[1] = points_list[i].ndx[1];
		pos[2] = points_list[i].ndx[2];
		
		ImageType3D::IndexType current_idx;
		current_idx[0] = pos[0];
		current_idx[1] = pos[1];
		current_idx[2] = pos[2]-padz; // PADZ???


		////////// Code for features from ODFs ////////////////////
		VBTNode node = this->getVBTFeatures(current_idx);
		
		////////////////// Code for scale features /////////////////////////
		float likelihood = 0.0;
		float radius = getRadiusAndLikelihood(pos, likelihood);

		double radiusThresh = 2.0;
		if(radius < radiusThresh)
			continue;

		//std::cout << i << std::endl;
		
		////////////////// Code for nearest nuiclei /////////////////////////
		float double_scale_nuclei = 20;
		double nucleus_distance = 0.0;
		nucleus_distance = this->SomaDistanceMapImage->GetPixel(current_idx);

		double nucDistThresh = 10.0 * double_scale_nuclei;
		if(nucleus_distance > nucDistThresh)
			continue;

		////////////////////// Code for eigen values based features ////////////////////////////
		float ballness, plateness, vesselness, noiseness;
		this->getHessianEigenFeatures(current_idx, img_max_val, ballness, plateness, vesselness, noiseness);
	
		/////////////////////// Code for intensity based features ///////////////////////////////
		float max_intensity, min_intensity, mean_intensity, var_intensity, intensity;
		intensity = this->PaddedCurvImage->GetPixel(current_idx);
		
		if(!this->getLocalIntensityFeatures(current_idx, sz, radius, max_intensity, min_intensity, mean_intensity, var_intensity))
			continue;

		if(IDImage->GetPixel(current_idx) == 0){

			//IDImage->SetPixel(current_idx, IDIndex+1);
			//IDImage->SetPixel(current_idx, i);

			CandidateRootPoint a_root;
			a_root.featureVector.node = HeapNode_astro(points_list[i]);
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

			a_root.featureVector.radiusVBT = node.scale;
			a_root.featureVector.likelihoodVBT = node.likelihood;
			a_root.featureVector.odfFeatures.std = node.odfFeatures.std;
			a_root.featureVector.odfFeatures.mean = node.odfFeatures.mean;
			a_root.featureVector.odfFeatures.energy = node.odfFeatures.energy;
			
			a_root.featureVector.odfFeatures.ODFModeVals = node.odfFeatures.ODFModeVals;
			if(node.odfFeatures.ODFModeVals.size() == 1){		
				a_root.featureVector.odfFeatures.ODFModeVals.push_back(0);
				a_root.featureVector.odfFeatures.ODFModeVals.push_back(0);
			}
			else if(node.odfFeatures.ODFModeVals.size() == 2)
				a_root.featureVector.odfFeatures.ODFModeVals.push_back(0);

			root_points[i] = a_root;

			//IDIndex++;			
			
			//std::cout << node.odfFeatures.ODFModeVals[0] << ", " << node.odfFeatures.energy << "," << node.scale <<std::endl;
		}
	}

	for(int i = 0; i < root_points.size(); i++){
		if(root_points[i].featureVector.radius == -1)
			continue;
		
		if(this->IDImage->GetPixel(root_points[i].featureVector.node.ndx) == 0){
			root_points[i].featureVector.ID = i;
			this->IDImage->SetPixel(root_points[i].featureVector.node.ndx, i);
			this->AllRootPoints.push_back(root_points[i]);
		}
	}

	std::cout << "Done with computing root point features: " << this->AllRootPoints.size() << std::endl;

	for(int i = 0; i < this->AllRootPoints.size(); i++){
		
		vtkSmartPointer<vtkVariantArray> table_row = vtkSmartPointer<vtkVariantArray>::New();
		//table_row->InsertNextValue(vtkVariant(i+1));
		table_row->InsertNextValue(vtkVariant(this->AllRootPoints[i].featureVector.ID));
		table_row->InsertNextValue(vtkVariant(this->AllRootPoints[i].featureVector.node.ndx[0]));
		table_row->InsertNextValue(vtkVariant(this->AllRootPoints[i].featureVector.node.ndx[1]));
		table_row->InsertNextValue(vtkVariant(this->AllRootPoints[i].featureVector.node.ndx[2] -  padz));
		table_row->InsertNextValue(vtkVariant(this->AllRootPoints[i].featureVector.radius));
		table_row->InsertNextValue(vtkVariant(this->AllRootPoints[i].featureVector.sphere_likelihood));
		table_row->InsertNextValue(vtkVariant(this->AllRootPoints[i].featureVector.ballness));
		table_row->InsertNextValue(vtkVariant(this->AllRootPoints[i].featureVector.plateness));
		table_row->InsertNextValue(vtkVariant(this->AllRootPoints[i].featureVector.vesselness));
		table_row->InsertNextValue(vtkVariant(this->AllRootPoints[i].featureVector.intensity));
		table_row->InsertNextValue(vtkVariant(this->AllRootPoints[i].featureVector.meanIntensity));
		table_row->InsertNextValue(vtkVariant(this->AllRootPoints[i].featureVector.varianceIntensity));
		table_row->InsertNextValue(vtkVariant(this->AllRootPoints[i].featureVector.maxIntensity));
		table_row->InsertNextValue(vtkVariant(this->AllRootPoints[i].featureVector.minIntensity));
		table_row->InsertNextValue(vtkVariant(this->AllRootPoints[i].featureVector.nucleusDistance));
		table_row->InsertNextValue(vtkVariant(this->AllRootPoints[i].featureVector.radiusVBT));
		table_row->InsertNextValue(vtkVariant(this->AllRootPoints[i].featureVector.likelihoodVBT));
		table_row->InsertNextValue(vtkVariant(this->AllRootPoints[i].featureVector.odfFeatures.std));
		table_row->InsertNextValue(vtkVariant(this->AllRootPoints[i].featureVector.odfFeatures.mean));
		table_row->InsertNextValue(vtkVariant(this->AllRootPoints[i].featureVector.odfFeatures.energy));
		table_row->InsertNextValue(vtkVariant(this->AllRootPoints[i].featureVector.odfFeatures.ODFModeVals[0]));
		table_row->InsertNextValue(vtkVariant(this->AllRootPoints[i].featureVector.odfFeatures.ODFModeVals[1]));
		table_row->InsertNextValue(vtkVariant(this->AllRootPoints[i].featureVector.odfFeatures.ODFModeVals[2]));


		features_table->InsertNextRow(table_row);
	}
	 
	features_table_vec[0] = features_table;
	out_images[0] = IDImage; 
		
	if(writeResult){
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

HeapNode_astro::HeapNode_astro(){
}

RootPointFeatureVector::RootPointFeatureVector(){

	this->radius = -1;
	HeapNode_astro();
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
			root_point.featureVector.node = HeapNode_astro(idx, 0);

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

void AstroTracer::ReadRootPointsPipeline(const std::vector<vtkSmartPointer<vtkTable> > roots_feature_table)
{
	//bool prediction_found = false;
	//for(int i = (int)roots_feature_table[0]->GetNumberOfColumns(); i > 0; i--)
	//{
	//	std::string current_column = std::string(roots_feature_table[0]->GetColumnName(i));
	//	if(current_column.find("prediction") != std::string::npos )
	//	{
	//		prediction_found = true;
	//		break;
	//	}	
	//}

	CandidateRootPoint root_point;
	itk::Index<3> root_idx;
	for(int i = 0; i < roots_feature_table[0]->GetNumberOfRows(); i++)
	{
		root_point.featureVector.ID = roots_feature_table[0]->GetValueByName(i, "ID").ToInt();

		root_idx[0] = roots_feature_table[0]->GetValueByName(i, "x").ToInt(); 
		root_idx[1] = roots_feature_table[0]->GetValueByName(i, "y").ToInt();
		root_idx[2] = roots_feature_table[0]->GetValueByName(i, "z").ToInt();
		root_point.featureVector.node = HeapNode_astro(root_idx, 0);

		root_point.featureVector.radius = roots_feature_table[0]->GetValueByName(i, "radius").ToFloat();
		root_point.featureVector.sphere_likelihood = roots_feature_table[0]->GetValueByName(i, "sph_likelihood").ToFloat();
		root_point.featureVector.ballness = roots_feature_table[0]->GetValueByName(i, "ball").ToFloat();
		root_point.featureVector.plateness = roots_feature_table[0]->GetValueByName(i, "plate").ToFloat();
		root_point.featureVector.vesselness = roots_feature_table[0]->GetValueByName(i, "vessel").ToFloat();
		root_point.featureVector.intensity = roots_feature_table[0]->GetValueByName(i, "int").ToFloat();
		root_point.featureVector.meanIntensity = roots_feature_table[0]->GetValueByName(i, "m_int").ToFloat();
		root_point.featureVector.varianceIntensity = roots_feature_table[0]->GetValueByName(i, "v_int").ToFloat();
		root_point.featureVector.maxIntensity = roots_feature_table[0]->GetValueByName(i, "mx_int").ToFloat();
		root_point.featureVector.minIntensity = roots_feature_table[0]->GetValueByName(i, "mn_int").ToFloat();
		root_point.featureVector.nucleusDistance = roots_feature_table[0]->GetValueByName(i, "nucl_dist").ToFloat();
		
		//if(prediction_found){
		//	root_point.classValue = roots_feature_table[0]->GetValueByName(i, "prediction_active_mg").ToInt();
		//	root_point.confidenceMeasure = roots_feature_table[0]->GetValueByName(i, "confidence_mg").ToFloat();

		//	if(root_point.classValue == 1)
		//		root_point.isRootPoint = true;
		//	else
		//		root_point.isRootPoint = false;

		//	if(root_point.isRootPoint)
		//		this->CandidateRootPoints.push_back(root_point);
		//}
		//else{
		root_point.isRootPoint = true;
		this->CandidateRootPoints.push_back(root_point);
		//}
	}
	
	if(this->CandidateRootPoints.empty()){
		std::cout << "Root points table is empty! Returning. " << std::endl;
		return;
	}
	std::cout << "Root points table read. " << this->CandidateRootPoints.size() << std::endl;
}

void AstroTracer::GetCentroidsForTracing(std::string outputFname, std::string finalIDImageFileName){
	
	int centroid_count=0;
	
	//Creating centroids.txt file with class 1 roots which are closest to a nucleus 

	//Preparing IDImage
	LabelImageType3D::RegionType id_reg;
	LabelImageType3D::IndexType id_st;
	LabelImageType3D::SizeType id_sz = this->PaddedCurvImage->GetBufferedRegion().GetSize();

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

	ImageType3D::IndexType starting_index_nuclei, end_index_nuclei;
	ImageType3D::SizeType sub_volume_size_nuclei;
	ImageType3D::RegionType sub_volume_region_nuclei;
	ImageType3D::Pointer sub_volume_nuclei;
	
	typedef itk::BinaryThresholdImageFilter<LabelImageType3D, CharImageType3D> ThresholdFilterType;
	ThresholdFilterType::Pointer threshold_filter = ThresholdFilterType::New();
	threshold_filter->SetLowerThreshold(1);
	threshold_filter->SetInsideValue(255);
	threshold_filter->SetOutsideValue(0);
	threshold_filter->SetInput(this->SomaImage);
	threshold_filter->Update();

	SignedMaurerDistanceMapImageFilterType::Pointer MaurerFilter = SignedMaurerDistanceMapImageFilterType::New();
	MaurerFilter->SetInput(threshold_filter->GetOutput());
	MaurerFilter->SetSquaredDistance(false);
	MaurerFilter->SetUseImageSpacing(false);
	MaurerFilter->SetInsideIsPositive(false);
	MaurerFilter->Update();

	LabelImageType3D::SizeType sz = this->SomaImage->GetBufferedRegion().GetSize();

	
	LabelGeometryFilterType::Pointer label_geometry_filter = LabelGeometryFilterType::New();
	label_geometry_filter->SetInput(this->SomaImage); 
	label_geometry_filter->Update();
	
	//std::cout << "Root points size: " << this->CandidateRootPoints.size() << std::endl;
	//std::cout << "Nuclei points size: " << this->NucleiObjects.size() << std::endl;

	//Loop over nuclei
	for(size_t i = 0; i < this->NucleiObjects.size(); i++){

		// Use only astrocyte nuclei
		if(this->NucleiObjects[i].classValue != 1)
			continue;

		//std::cout << "Ye!! Asto nuclei found!! " << std::endl;

		// ROI proportional to nuclei scale, assuming spherical nuclei.
		//float double_scale_nuclei = 0.5*std::pow((float)this->NucleiObjects[i].intrinsicFeatures.boundingBoxVolume, (float)0.333333);
		//double_scale_nuclei = 2.0 * double_scale_nuclei; 
		
		CharImageType3D::IndexType current_idx;
		current_idx[0] = this->NucleiObjects[i].intrinsicFeatures.centroid.ndx[0];
		current_idx[1] = this->NucleiObjects[i].intrinsicFeatures.centroid.ndx[1];
		current_idx[2] = this->NucleiObjects[i].intrinsicFeatures.centroid.ndx[2] - padz; 

		
		float sph_rad = 0.5 * label_geometry_filter->GetMajorAxisLength(this->SomaImage->GetPixel(current_idx));
		float double_scale_nuclei = 1.5 * sph_rad;

		starting_index_nuclei[0] = current_idx[0] - double_scale_nuclei; starting_index_nuclei[1] = current_idx[1] - double_scale_nuclei; starting_index_nuclei[2] = current_idx[2] - double_scale_nuclei;
		end_index_nuclei[0] = current_idx[0] + double_scale_nuclei; end_index_nuclei[1] = current_idx[1] + double_scale_nuclei; end_index_nuclei[2] = current_idx[2] + double_scale_nuclei;

		
		//std::cout << starting_index_nuclei[0] << ", " << starting_index_nuclei[1] << ", " << starting_index_nuclei[2] << std::endl;

		if(starting_index_nuclei[0] < 0)
			starting_index_nuclei[0] = 0;
		if(starting_index_nuclei[1] < 0)
			starting_index_nuclei[1] = 0;
		if(starting_index_nuclei[2] < 0)
			starting_index_nuclei[2] = 0;
		if(end_index_nuclei[0] > sz[0])
			end_index_nuclei[0] = sz[0];
		if(end_index_nuclei[1] > sz[1])
			end_index_nuclei[1] = sz[1];
		if(end_index_nuclei[2] > sz[2])
			end_index_nuclei[2] = sz[2];

		//if ( (starting_index_nuclei[0] < 0) || (starting_index_nuclei[1] < 0) || (starting_index_nuclei[2] < 0) ||
		//	(end_index_nuclei[0] > (unsigned int)sz[0]) || (end_index_nuclei[1] > (unsigned int)sz[1]) ||
		//	(end_index_nuclei[2] > (unsigned int)sz[2]) )
		//	continue;



		//std::cout << "Nuclei: "  << i << " Scale: " << double_scale_nuclei << std::endl;

		//sub_volume_size_nuclei[0] = 2 * double_scale_nuclei; sub_volume_size_nuclei[1] = 2 * double_scale_nuclei; sub_volume_size_nuclei[2] = 2 * double_scale_nuclei;
		sub_volume_size_nuclei[0] = end_index_nuclei[0] - starting_index_nuclei[0];
		sub_volume_size_nuclei[1] = end_index_nuclei[1] - starting_index_nuclei[1];
		sub_volume_size_nuclei[2] = end_index_nuclei[2] - starting_index_nuclei[2];
		

		sub_volume_region_nuclei.SetIndex(starting_index_nuclei);
		sub_volume_region_nuclei.SetSize(sub_volume_size_nuclei);

		typedef itk::RegionOfInterestImageFilter<ImageType3D, ImageType3D> VolumeOfInterestFilterType_nuclei2;
		VolumeOfInterestFilterType_nuclei2::Pointer sub_volume_filter_nuclei = VolumeOfInterestFilterType_nuclei2::New();
		sub_volume_filter_nuclei->SetInput(MaurerFilter->GetOutput());
		sub_volume_filter_nuclei->SetRegionOfInterest(sub_volume_region_nuclei);
		sub_volume_filter_nuclei->Update();
		sub_volume_nuclei = sub_volume_filter_nuclei->GetOutput();
		
		ImageType3D::Pointer distance_map = sub_volume_nuclei; //MaurerFilter->GetOutput();
				
		double cur_distance = 0.0;
		double min_distance = 1000.0;
		//double max_distance = -1000.0; 
		//double mean_distance = 0.0;
		//double variance_distance = 0.0;
		//double acc_distance = 0.0;
		//int n_roots = 0;
		
		std::vector<double> distance_array;
		std::multimap<double, HeapNode_astro> distance_idx_map;
		std::multimap<double, HeapNode_astro>::reverse_iterator map_rit;
		typedef std::pair<double, HeapNode_astro> distance_idx_pair_type;
		
		LabelImageType3D::IndexType min_root_idx;
		min_root_idx[0] = this->CandidateRootPoints[0].featureVector.node.ndx[0];
		min_root_idx[1] = this->CandidateRootPoints[0].featureVector.node.ndx[1];
		min_root_idx[2] = this->CandidateRootPoints[0].featureVector.node.ndx[2] - padz;

		for(size_t j = 0; j < this->CandidateRootPoints.size(); j++){

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
				distance_idx_map.insert(distance_idx_pair_type(this->CandidateRootPoints[j].featureVector.radius, HeapNode_astro(current_root_idx, 0)));

		}// End of loop over candidate roots	

		if(!distance_idx_map.empty()){

			// Among the neighboring root points, keep the root point with highest scale 
			this->CentroidListForTracing.push_back(distance_idx_map.rbegin()->second);
			this->CentroidScales.push_back(distance_idx_map.rbegin()->first);

			/*if(distance_idx_map.size() == 1){

				this->CentroidListForTracing.push_back(distance_idx_map.begin()->second);

				this->RefinedRootImage->SetPixel(min_root_idx, 255);
			}
			else{

				double cur_scale;
				HeapNode_astro max_scale_node;
				std::multimap<double, HeapNode_astro> max_scale_map(distance_idx_map.rbegin(), distance_idx_map.rbegin()+1);
				
				for(map_rit = distance_idx_map.rbegin()+1; map_rit != distance_idx_map.rend(); map_rit++){	
					
					if(max_scale_map.begin()->first - map_rit->first > 0.01)
						max_scale_map.insert(distance_idx_pair_type(map_rit->first, map_rit->second))
				}
				

				this->CentroidListForTracing.push_back(HeapNode_astro(min_root_idx, 0));
			}*/
		}
		//else{
			// If the map is empty, it means that the astro nuclei did not have any root points in the neighborhood.
			// In this case, just take the max intensity point in the neighborhood of this nuclei as a root point?
		//}

		
	}//end of loop over nuclei

	std::cout << "Before density filtering: " << this->CentroidListForTracing.size() << std::endl;
	
	// Filter the centroids based on density
	double centroid_nhood = 30.0; //25; //50; 
	std::vector<bool> neighbor_flags(this->CentroidListForTracing.size(), false);
	
	for(int i = 0; i < this->CentroidListForTracing.size(); i++){

		if(neighbor_flags[i])
			continue;
		
		itk::Index<3> cur_idx = this->CentroidListForTracing[i].ndx;
		std::vector<HeapNode_astro> neighbor_points;
		std::multimap<double, HeapNode_astro> neighbor_point_map;
		std::multimap<double, HeapNode_astro>::iterator neighbor_point_map_itr;
		typedef std::pair<double, HeapNode_astro> distance_idx_pair_type;

		for(int j = 0; j < this->CentroidListForTracing.size(); j++){

			if(j == i)
				continue;

			itk::Index<3> neighbor_idx = this->CentroidListForTracing[j].ndx; 
			
			// Please use Eucleidien distance
			//if(std::abs((int)(cur_idx[0] - neighbor_idx[0])) > centroid_nhood || std::abs((int)(cur_idx[1] - neighbor_idx[1])) > centroid_nhood || 
			//	std::abs((int)(cur_idx[2] - neighbor_idx[2])) > centroid_nhood)
			//	continue;
			if(std::sqrt(std::pow((double)(cur_idx[0] - neighbor_idx[0]), 2) + std::pow((double)(cur_idx[1] - neighbor_idx[1]), 2) + std::pow((double)(cur_idx[2] - neighbor_idx[2]), 2)) > centroid_nhood)
				continue;
			
			neighbor_flags[j] = true;
			neighbor_points.push_back(HeapNode_astro(neighbor_idx, 0));
			neighbor_point_map.insert(distance_idx_pair_type(this->CentroidScales[j], HeapNode_astro(neighbor_idx, 0)));
		}

		if(neighbor_points.empty()){
			neighbor_flags[i] = false;	
			this->DensityFilteredCentroidListForTracing.push_back(HeapNode_astro(cur_idx, 0));
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
		
		this->DensityFilteredCentroidListForTracing.push_back(HeapNode_astro(centroid_of_centroid, 0));
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
		
		this->DensityFilteredCentroidListForTracing.push_back(HeapNode_astro(centroid_of_centroid, 0));
		*/

		//Replace neighboring root points with the highest scale point
		this->DensityFilteredCentroidListForTracing.push_back(neighbor_point_map.rbegin()->second);

		//neighbor_point_map.clear();
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

	
	std::cout << "Points list size: " << this->CandidateRootPoints.size() << std::endl;
	std::cout << "Centroids list size: " << centroid_count << std::endl;


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
	std::cout << "Done with FinalRootImage writing: " << finalIDImageFileName << std::endl;

	centroid_points.close();
	//End of creating centroids.txt file
}

void AstroTracer::GetCentroidsForTracingPipeline(std::string outputFname, std::string finalIDImageFileName, unsigned int padz, ImageType3D::RegionType regionLocal_inside, std::vector<vtkSmartPointer<vtkTable> >& centroids_table, std::vector<LabelImageType3D::Pointer>& out_images, const bool writeResult){

	int centroid_count = 0;
	
	//Creating centroids.txt file with class 1 roots which are closest to a nucleus 

	//Preparing IDImage
	LabelImageType3D::RegionType id_reg;
	LabelImageType3D::IndexType id_st;
	LabelImageType3D::SizeType id_sz = this->PaddedCurvImage->GetBufferedRegion().GetSize();

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

	CharImageType3D::IndexType starting_index_nuclei, end_index_nuclei;
	LabelImageType3D::SizeType sz = id_sz; //this->PaddedCurvImage->GetBufferedRegion().GetSize();
	
	std::cout << "Root points size: " << this->CandidateRootPoints.size() << std::endl;
	std::cout << "Nuclei points size: " << this->NucleiObjects.size() << std::endl;

	//Loop over nuclei
	for(size_t i = 0; i < this->NucleiObjects.size(); i++){

		// Use only astrocyte nuclei
		if(this->NucleiObjects[i].classValue != 1)
			continue;

		//std::cout << "Ye!! Asto nuclei found!! " << std::endl;

		// ROI proportional to nuclei scale, assuming spherical nuclei.
		float double_scale_nuclei = 0.5*std::pow((float)this->NucleiObjects[i].intrinsicFeatures.boundingBoxVolume, (float)0.333333);
		double_scale_nuclei = 2.0 * double_scale_nuclei; 
		
		CharImageType3D::IndexType current_idx;
		current_idx[0] = this->NucleiObjects[i].intrinsicFeatures.centroid.ndx[0];
		current_idx[1] = this->NucleiObjects[i].intrinsicFeatures.centroid.ndx[1];
		current_idx[2] = this->NucleiObjects[i].intrinsicFeatures.centroid.ndx[2] - padz; 

		if(!regionLocal_inside.IsInside(current_idx))
			continue;

		starting_index_nuclei[0] = current_idx[0] - double_scale_nuclei; starting_index_nuclei[1] = current_idx[1] - double_scale_nuclei; starting_index_nuclei[2] = current_idx[2] - double_scale_nuclei;
		end_index_nuclei[0] = current_idx[0] + double_scale_nuclei; end_index_nuclei[1] = current_idx[1] + double_scale_nuclei; end_index_nuclei[2] = current_idx[2] + double_scale_nuclei;

		
		//std::cout << starting_index_nuclei[0] << ", " << starting_index_nuclei[1] << ", " << starting_index_nuclei[2] << std::endl;

		if(starting_index_nuclei[0] < 0)
			starting_index_nuclei[0] = 0;
		if(starting_index_nuclei[1] < 0)
			starting_index_nuclei[1] = 0;
		if(starting_index_nuclei[2] < 0)
			starting_index_nuclei[2] = 0;
		if(end_index_nuclei[0] > sz[0])
			end_index_nuclei[0] = sz[0];
		if(end_index_nuclei[1] > sz[1])
			end_index_nuclei[1] = sz[1];
		if(end_index_nuclei[2] > sz[2])
			end_index_nuclei[2] = sz[2];

		//if ( (starting_index_nuclei[0] < 0) || (starting_index_nuclei[1] < 0) || (starting_index_nuclei[2] < 0) ||
		//	(end_index_nuclei[0] > (unsigned int)sz[0]) || (end_index_nuclei[1] > (unsigned int)sz[1]) ||
		//	(end_index_nuclei[2] > (unsigned int)sz[2]) )
		//	continue;


		std::cout << "Nuclei: "  << i << " Scale: " << double_scale_nuclei << std::endl;
				
		double cur_distance = 0.0;
		std::multimap<double, HeapNode_astro> distance_idx_map;
		std::multimap<double, HeapNode_astro>::reverse_iterator map_rit;
		typedef std::pair<double, HeapNode_astro> distance_idx_pair_type;
		
		for(size_t j = 0; j < this->CandidateRootPoints.size(); j++){

			//std::cout << "Nuclei: "  << i << " Root: " << j << std::endl;
		
			LabelImageType3D::IndexType current_root_idx;
			current_root_idx[0] = this->CandidateRootPoints[j].featureVector.node.ndx[0];
			current_root_idx[1] = this->CandidateRootPoints[j].featureVector.node.ndx[1];
			current_root_idx[2] = this->CandidateRootPoints[j].featureVector.node.ndx[2] - padz;

			int offset = 2; //1;
			if(current_root_idx[0] < starting_index_nuclei[0]+offset || current_root_idx[1] < starting_index_nuclei[1]+offset || current_root_idx[2] < starting_index_nuclei[2]+offset ||
				current_root_idx[0] > end_index_nuclei[0]-offset || current_root_idx[1] > end_index_nuclei[1]-offset || current_root_idx[2] > end_index_nuclei[2]-offset)
				continue;
			
			
			distance_idx_map.insert(distance_idx_pair_type(this->CandidateRootPoints[j].featureVector.radius, HeapNode_astro(current_root_idx, 0)));

		}
		if(!distance_idx_map.empty()){

			// Among the neighboring root points, keep the root point with highest scale 
			this->CentroidListForTracing.push_back(distance_idx_map.rbegin()->second);
			this->CentroidScales.push_back(distance_idx_map.rbegin()->first);

		}
		//else{
			// If the map is empty, it means that the astro nuclei did not have any root points in the neighborhood.
			// In this case, just take the max intensity point in the neighborhood of this nuclei as a root point?
		//}

		
	}//end of loop over nuclei

	std::cout << "Before density filtering: " << this->CentroidListForTracing.size() << std::endl;
	
	// Filter the centroids based on density
	double centroid_nhood = 30.0; //25; //50; 
	std::vector<bool> neighbor_flags(this->CentroidListForTracing.size(), false);
	
	for(int i = 0; i < this->CentroidListForTracing.size(); i++){

		if(neighbor_flags[i])
			continue;
		
		itk::Index<3> cur_idx = this->CentroidListForTracing[i].ndx;
		std::vector<HeapNode_astro> neighbor_points;
		std::multimap<double, HeapNode_astro> neighbor_point_map;
		std::multimap<double, HeapNode_astro>::iterator neighbor_point_map_itr;
		typedef std::pair<double, HeapNode_astro> distance_idx_pair_type;

		for(int j = 0; j < this->CentroidListForTracing.size(); j++){

			if(j == i)
				continue;

			itk::Index<3> neighbor_idx = this->CentroidListForTracing[j].ndx; 
			
			// Please use Eucleidien distance
			//if(std::abs((int)(cur_idx[0] - neighbor_idx[0])) > centroid_nhood || std::abs((int)(cur_idx[1] - neighbor_idx[1])) > centroid_nhood || 
			//	std::abs((int)(cur_idx[2] - neighbor_idx[2])) > centroid_nhood)
			//	continue;
			if(std::sqrt(std::pow((double)(cur_idx[0] - neighbor_idx[0]), 2) + std::pow((double)(cur_idx[1] - neighbor_idx[1]), 2) + std::pow((double)(cur_idx[2] - neighbor_idx[2]), 2)) > centroid_nhood)
				continue;
			
			neighbor_flags[j] = true;
			neighbor_points.push_back(HeapNode_astro(neighbor_idx, 0));
			neighbor_point_map.insert(distance_idx_pair_type(this->CentroidScales[j], HeapNode_astro(neighbor_idx, 0)));
		}

		if(neighbor_points.empty()){
			neighbor_flags[i] = false;	
			this->DensityFilteredCentroidListForTracing.push_back(HeapNode_astro(cur_idx, 0));
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
		
		this->DensityFilteredCentroidListForTracing.push_back(HeapNode_astro(centroid_of_centroid, 0));
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
		
		this->DensityFilteredCentroidListForTracing.push_back(HeapNode_astro(centroid_of_centroid, 0));
		*/

		//Replace neighboring root points with the highest scale point
		this->DensityFilteredCentroidListForTracing.push_back(neighbor_point_map.rbegin()->second);
	}

	std::cout << "After density filtering: " << this->DensityFilteredCentroidListForTracing.size() << std::endl;

	for(int i = 0; i < this->DensityFilteredCentroidListForTracing.size(); i++)
		this->FinalRootsImage->SetPixel(this->DensityFilteredCentroidListForTracing[i].ndx, 255);

	out_images[0] = this->FinalRootsImage;

	vtkSmartPointer<vtkDoubleArray> x_array = vtkSmartPointer<vtkDoubleArray>::New();
	vtkSmartPointer<vtkDoubleArray> y_array = vtkSmartPointer<vtkDoubleArray>::New();
	vtkSmartPointer<vtkDoubleArray> z_array = vtkSmartPointer<vtkDoubleArray>::New();
	
	for(int i = 0; i < this->DensityFilteredCentroidListForTracing.size(); i++){
		x_array->InsertNextValue(this->DensityFilteredCentroidListForTracing[i].ndx[0]);
		y_array->InsertNextValue(this->DensityFilteredCentroidListForTracing[i].ndx[1]);
		z_array->InsertNextValue(this->DensityFilteredCentroidListForTracing[i].ndx[2]);		
	}

	x_array->SetName("x");
	y_array->SetName("y");
	z_array->SetName("z");
	
	centroids_table[0]->AddColumn(x_array);
	centroids_table[0]->AddColumn(y_array);
	centroids_table[0]->AddColumn(z_array);
	

	if(writeResult){
		std::ofstream centroid_points;
		centroid_points.open(outputFname.c_str(), std::ios::out);

		if(centroid_points.good()){

			for(int i = 0; i < this->DensityFilteredCentroidListForTracing.size(); i++){
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
		
		std::cout << "Points list size: " << this->CandidateRootPoints.size() << std::endl;
		std::cout << "Centroids list size: " << centroid_count << std::endl;


		itk::ImageFileWriter<LabelImageType3D>::Pointer IDImageWriter = itk::ImageFileWriter<LabelImageType3D>::New();
		IDImageWriter->SetFileName(finalIDImageFileName.c_str());
		IDImageWriter->SetInput(this->FinalRootsImage);
		try{
			IDImageWriter->Update();
		}
		catch(itk::ExceptionObject e){
			std::cout << e << std::endl;
			//return EXIT_FAILURE;
		}
		std::cout << "Done with FinalRootImage writing. " << std::endl;

		centroid_points.close();
	}
}

void AstroTracer::ReadStartPointsInternal(){

	if(this->DensityFilteredCentroidListForTracing.empty()){
		std::cout << "Internal start points list is empty. " << std::endl;
		return;
	}

	//if(!this->StartPoints.empty())
	//	this->StartPoints.clear();
	

	for(int i = 0; i < this->DensityFilteredCentroidListForTracing.size(); i++)
		this->StartPoints.push_back(this->DensityFilteredCentroidListForTracing[i].ndx);
	
	std::cout << "Start points loaded internally: " << this->StartPoints.size();
	
}

void AstroTracer::ReadStartPointsFromVTKTable(const std::vector<vtkSmartPointer<vtkTable> > centroids_table){

	if(centroids_table.empty()){
		std::cout << "Centroids table vector is empty! Returning.";
		return;
	}
	if(centroids_table[0]->GetNumberOfRows() < 1){
		std::cout << "Centroids table is empty! Returning.";
		return;
	}

	ImageType3D::IndexType idx;
	for(int i = 0; i < centroids_table[0]->GetNumberOfRows(); i++){
		idx[0] = centroids_table[0]->GetValueByName(i, "x").ToFloat();
		idx[1] = centroids_table[0]->GetValueByName(i, "y").ToFloat();
		idx[2] = centroids_table[0]->GetValueByName(i, "z").ToFloat();

		this->StartPoints.push_back(idx);
	}
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
			nuclei_object.intrinsicFeatures.centroid = HeapNode_astro(idx, 0);

			nuclei_object.intrinsicFeatures.volume = atof(str_vec[4].c_str());
			nuclei_object.intrinsicFeatures.integratedIntensity = atof(str_vec[5].c_str());
			nuclei_object.intrinsicFeatures.eccentricity = atof(str_vec[6].c_str());
			nuclei_object.intrinsicFeatures.elongation = atof(str_vec[7].c_str());
			nuclei_object.intrinsicFeatures.orientation = atof(str_vec[8].c_str());
			//nuclei_object.intrinsicFeatures.boundingBoxVolume = atof(str_vec[9].c_str());
			nuclei_object.intrinsicFeatures.boundingBoxVolume = nuclei_object.intrinsicFeatures.volume;

			//nuclei_object.intrinsicFeatures.integratedIntensity = atof(str_vec[10].c_str());
			nuclei_object.intrinsicFeatures.meanIntensity = atof(str_vec[11].c_str());
			nuclei_object.intrinsicFeatures.minIntensity = atof(str_vec[12].c_str());
			nuclei_object.intrinsicFeatures.maxIntensity = atof(str_vec[13].c_str());

			nuclei_object.intrinsicFeatures.varianceIntensity = atof(str_vec[15].c_str());
			nuclei_object.intrinsicFeatures.meanSurfaceGradient = atof(str_vec[16].c_str());
			nuclei_object.intrinsicFeatures.interiorGradient = atof(str_vec[17].c_str());
			nuclei_object.intrinsicFeatures.surfaceIntensity = atof(str_vec[18].c_str());
			nuclei_object.intrinsicFeatures.interiorIntensity = atof(str_vec[19].c_str());
			nuclei_object.intrinsicFeatures.intensityRatio = atof(str_vec[20].c_str());
			nuclei_object.intrinsicFeatures.convexity = atof(str_vec[21].c_str());

			nuclei_object.intrinsicFeatures.radiusVariation = atof(str_vec[22].c_str());
			nuclei_object.intrinsicFeatures.surfaceArea = atof(str_vec[23].c_str());
			nuclei_object.intrinsicFeatures.shapeMeasure = atof(str_vec[24].c_str());
			nuclei_object.intrinsicFeatures.sharedBoundary = atof(str_vec[25].c_str());
			
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
			
			// Only for 5 channel analysis
			//nuclei_object.associativeFeatures.vessel_total = atof(str_vec[41].c_str());
			//nuclei_object.associativeFeatures.vessel_avg = atof(str_vec[42].c_str());
			//nuclei_object.associativeFeatures.vessel_surr = atof(str_vec[43].c_str());

			//TEMP: DARPA HACK
			nuclei_object.classValue = atof(str_vec[41].c_str());
			nuclei_object.confidenceMeasure = atof(str_vec[42].c_str());
			nuclei_object.distanceToElectrode = atof(str_vec[43].c_str());

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

void AstroTracer::ReadNucleiFeaturesPipeline(const std::vector<vtkSmartPointer<vtkTable> > nuclei_features_table){
	
	//bool prediction_found = false, root_found = false;
	//for(int i = nuclei_features_table[0]->GetNumberOfColumns(); i > 0; i--){
	//	std::string current_column = nuclei_features_table[0]->GetColumnName(i);
	//	if(current_column.find("prediction") != std::string::npos ){
	//		prediction_found = true;
	//		break;
	//	}
	//}
	//for(int i = nuclei_features_table[0]->GetNumberOfColumns(); i > 0; i--){
	//	std::string current_column = nuclei_features_table[0]->GetColumnName(i);
	//	if(current_column.find("root") != std::string::npos ){
	//		root_found = true;
	//		break;
	//	}
	//}

	NucleiObject nuclei_object;
	itk::Index<3> nuc_idx;
	for(int i = 0; i < nuclei_features_table[0]->GetNumberOfRows(); i++)
	{
		nuclei_object.intrinsicFeatures.ID = nuclei_features_table[0]->GetValueByName(i, "ID").ToInt();
		
		nuc_idx[0] = nuclei_features_table[0]->GetValueByName(i, "centroid_x").ToInt();
		nuc_idx[1] = nuclei_features_table[0]->GetValueByName(i, "centroid_y").ToInt();
		nuc_idx[2] = nuclei_features_table[0]->GetValueByName(i, "centroid_z").ToInt();
		nuclei_object.intrinsicFeatures.centroid = HeapNode_astro(nuc_idx, 0);
		nuclei_object.intrinsicFeatures.volume = nuclei_features_table[0]->GetValueByName(i, "volume").ToFloat();
		//nuclei_object.intrinsicFeatures.integratedIntensity = nuclei_features_table[0]->GetValueByName(i, "sum_int").ToFloat();
		//nuclei_object.intrinsicFeatures.meanIntensity = nuclei_features_table[0]->GetValueByName(i, "mean_int").ToFloat();
		//nuclei_object.intrinsicFeatures.varianceIntensity = nuclei_features_table[0]->GetValueByName(i, "var_int").ToFloat();
		//nuclei_object.intrinsicFeatures.eccentricity = nuclei_features_table[0]->GetValueByName(i, "eccentricity").ToFloat();
		//nuclei_object.intrinsicFeatures.elongation = nuclei_features_table[0]->GetValueByName(i, "elongation").ToFloat();
		nuclei_object.intrinsicFeatures.boundingBoxVolume = nuclei_object.intrinsicFeatures.volume;

		//nuclei_object.intrinsicFeatures.meanSurfaceGradient = nuclei_features_table[0]->GetValueByName(i, "mean_surf_gradient").ToFloat();
		//nuclei_object.intrinsicFeatures.radiusVariation = nuclei_features_table[0]->GetValueByName(i, "radius_variation").ToFloat();
		//nuclei_object.intrinsicFeatures.shapeMeasure = nuclei_features_table[0]->GetValueByName(i, "shape_measure").ToFloat();
		//nuclei_object.intrinsicFeatures.energy = nuclei_features_table[0]->GetValueByName(i, "energy").ToFloat();
		//nuclei_object.intrinsicFeatures.entropy = nuclei_features_table[0]->GetValueByName(i, "entropy").ToFloat();
		//nuclei_object.intrinsicFeatures.inverseDiffMoment = nuclei_features_table[0]->GetValueByName(i, "inverse_diff_moment").ToFloat();
		//nuclei_object.intrinsicFeatures.inertia = nuclei_features_table[0]->GetValueByName(i, "inertia").ToFloat();
		//nuclei_object.intrinsicFeatures.clusterShade = nuclei_features_table[0]->GetValueByName(i, "cluster_shade").ToFloat();
		//nuclei_object.intrinsicFeatures.clusterProminence = nuclei_features_table[0]->GetValueByName(i, "cluster_prominence").ToFloat();
		//nuclei_object.associativeFeatures.astro_total = nuclei_features_table[0]->GetValueByName(i, "TRI_TOTAL").ToFloat();
		//nuclei_object.associativeFeatures.astro_avg = nuclei_features_table[0]->GetValueByName(i, "TRI_AVG").ToFloat();
		//nuclei_object.associativeFeatures.astro_surr = nuclei_features_table[0]->GetValueByName(i, "TRI_SURR").ToFloat();
		//nuclei_object.associativeFeatures.micro_total = nuclei_features_table[0]->GetValueByName(i, "GFP_TOTAL").ToFloat();
		//nuclei_object.associativeFeatures.micro_avg = nuclei_features_table[0]->GetValueByName(i, "GFP_AVG").ToFloat();
		//nuclei_object.associativeFeatures.micro_surr = nuclei_features_table[0]->GetValueByName(i, "GFP_SURR").ToFloat();
		//nuclei_object.associativeFeatures.neuro_total = nuclei_features_table[0]->GetValueByName(i, "Cy5_TOTAL").ToFloat();
		//nuclei_object.associativeFeatures.neuro_avg = nuclei_features_table[0]->GetValueByName(i, "Cy5_AVG").ToFloat();
		//nuclei_object.associativeFeatures.neuro_surr = nuclei_features_table[0]->GetValueByName(i, "Cy5_SURR").ToFloat();

		//if(root_found){
		//	nuclei_object.associativeFeatures.minRootDist = nuclei_features_table[0]->GetValueByName(i, "min_root_dist").ToFloat();
		//	nuclei_object.associativeFeatures.maxRootDist = nuclei_features_table[0]->GetValueByName(i, "max_root_dist").ToFloat();
		//	nuclei_object.associativeFeatures.meanRootDist = nuclei_features_table[0]->GetValueByName(i, "mean_root_dist").ToFloat();
		//	nuclei_object.associativeFeatures.varRootDist = nuclei_features_table[0]->GetValueByName(i, "var_root_dist").ToFloat();
		//	nuclei_object.associativeFeatures.nRoots = nuclei_features_table[0]->GetValueByName(i, "n_roots").ToInt();
		//}
		//if(prediction_found){
		//	nuclei_object.classValue = nuclei_features_table[0]->GetValueByName(i, "prediction_active_mg").ToInt();
		//	nuclei_object.confidenceMeasure = nuclei_features_table[0]->GetValueByName(i, "confidence_mg").ToFloat();
		//}
		
		this->NucleiObjects.push_back(nuclei_object);
	}

	if(this->NucleiObjects.empty()){
		std::cout << "Nuclei feature table is empty! Returning. " << std::endl;
		return; 
	}
	std::cout << "Nuclei features table read. " << this->NucleiObjects.size() << std::endl;
	
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

	LabelGeometryFilterType::Pointer label_geometry_filter = LabelGeometryFilterType::New();
	label_geometry_filter->SetInput(this->SomaImage); 
	label_geometry_filter->Update();
	
	std::cout << "Root points size: " << this->CandidateRootPoints.size() << std::endl;
	
	std::cout << "Computing nuclei features. " << std::endl;

	//Loop over nuclei
	for(size_t i = 0; i < this->NucleiObjects.size(); i++){

		CharImageType3D::IndexType current_idx;
		current_idx[0] = this->NucleiObjects[i].intrinsicFeatures.centroid.ndx[0];
		current_idx[1] = this->NucleiObjects[i].intrinsicFeatures.centroid.ndx[1];
		current_idx[2] = this->NucleiObjects[i].intrinsicFeatures.centroid.ndx[2] - padz; 

		// ROI proportional to nuclei scale, assuming spherical nuclei.
		//float sph_rad = std::pow((double)(0.23877 * this->NucleiObjects[i].intrinsicFeatures.volume), (double)0.333333);
		float sph_rad = 0.5 * label_geometry_filter->GetMajorAxisLength(this->SomaImage->GetPixel(current_idx));
		float double_scale_nuclei = 1.5 * sph_rad;

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

		sub_volume_size_nuclei[0] = 2.0 * double_scale_nuclei; sub_volume_size_nuclei[1] = 2.0 * double_scale_nuclei; sub_volume_size_nuclei[2] = 2.0 * double_scale_nuclei;

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
		for(size_t j = 0; j < this->CandidateRootPoints.size(); j++){

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

			mean_distance = acc_distance / (double)n_roots;
			for(int k = 0; k < distance_array.size(); k++)
				variance_distance = variance_distance + std::pow(distance_array[k] - mean_distance, 2);
		
			variance_distance = variance_distance / (double)n_roots;

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
			this->NucleiObjects[i].associativeFeatures.minRootDist = 100000.0; //double_scale_nuclei;
			this->NucleiObjects[i].associativeFeatures.maxRootDist = 100000.0; //double_scale_nuclei;
			this->NucleiObjects[i].associativeFeatures.meanRootDist = 100000.0; //double_scale_nuclei;
			this->NucleiObjects[i].associativeFeatures.varRootDist = 100000.0; //double_scale_nuclei;
			this->NucleiObjects[i].associativeFeatures.nRoots = 0;

		}//loop over candidate roots

		//distance_array.clear();

	}//end of loop over nuclei
}

void AstroTracer::ComputeFeaturesFromCandidateRootsPipeline(const ImageType3D::RegionType local_region, std::vector<vtkSmartPointer<vtkTable> >& nuclei_features_table, std::string TrainingFileName, std::string outputFname){

	if(this->CandidateRootPoints.empty() || this->NucleiObjects.empty() || this->SomaDistanceMapImage.IsNull()){
		std::cout << "Either root points table or nuclei table or the distance map image is empty. Returning." << std::endl;
		vtkSmartPointer<vtkDoubleArray> min_dist_array = vtkSmartPointer<vtkDoubleArray>::New();
		min_dist_array->SetName("min_root_dist");
		vtkSmartPointer<vtkDoubleArray> max_dist_array = vtkSmartPointer<vtkDoubleArray>::New();
		max_dist_array->SetName("max_root_dist");
		vtkSmartPointer<vtkDoubleArray> mean_dist_array = vtkSmartPointer<vtkDoubleArray>::New();
		mean_dist_array->SetName("mean_root_dist");
		vtkSmartPointer<vtkDoubleArray> var_dist_array = vtkSmartPointer<vtkDoubleArray>::New();
		var_dist_array->SetName("var_root_dist");
		vtkSmartPointer<vtkDoubleArray> n_roots_array = vtkSmartPointer<vtkDoubleArray>::New();
		n_roots_array->SetName("n_roots");

		for(int i = 0; i < this->NucleiObjects.size(); i++){
			min_dist_array->InsertNextValue(100000.0);
			max_dist_array->InsertNextValue(100000.0);
			mean_dist_array->InsertNextValue(100000.0);
			var_dist_array->InsertNextValue(100000.0);
			n_roots_array->InsertNextValue(0);
		}	


		nuclei_features_table[0]->AddColumn(min_dist_array);
		nuclei_features_table[0]->AddColumn(max_dist_array);
		nuclei_features_table[0]->AddColumn(mean_dist_array);
		nuclei_features_table[0]->AddColumn(var_dist_array);
		nuclei_features_table[0]->AddColumn(n_roots_array);
		//return;
	}
	
	else
	{
		ImageType3D::IndexType starting_index_nuclei, end_index_nuclei;
		ImageType3D::SizeType sub_volume_size_nuclei;
		ImageType3D::RegionType sub_volume_region_nuclei;
		ImageType3D::Pointer sub_volume_nuclei;

		//sz is the size of the tile plus the padding
		LabelImageType3D::SizeType sz = this->PaddedCurvImage->GetBufferedRegion().GetSize();
		
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

		LabelGeometryFilterType::Pointer label_geometry_filter = LabelGeometryFilterType::New();
		label_geometry_filter->SetInput(this->SomaImage);
		label_geometry_filter->Update();

		std::cout << "Root points size: " << this->CandidateRootPoints.size() << std::endl;

		std::cout << "Computing nuclei features. " << std::endl;

		//Loop over nuclei
		for(size_t i = 0; i < this->NucleiObjects.size(); i++){

			CharImageType3D::IndexType current_idx;
			current_idx[0] = this->NucleiObjects[i].intrinsicFeatures.centroid.ndx[0];
			current_idx[1] = this->NucleiObjects[i].intrinsicFeatures.centroid.ndx[1];
			current_idx[2] = this->NucleiObjects[i].intrinsicFeatures.centroid.ndx[2] - padz; 

			// ROI proportional to nuclei scale, assuming spherical nuclei.
			//float sph_rad = std::pow((double)(0.23877 * this->NucleiObjects[i].intrinsicFeatures.volume), (double)0.333333);
			float sph_rad = 0.5 * label_geometry_filter->GetMajorAxisLength(this->SomaImage->GetPixel(current_idx));
			float double_scale_nuclei = 1.5 * sph_rad;

			if(!local_region.IsInside(current_idx))
				continue;

			starting_index_nuclei[0] = current_idx[0] - double_scale_nuclei; starting_index_nuclei[1] = current_idx[1] - double_scale_nuclei; starting_index_nuclei[2] = current_idx[2] - double_scale_nuclei;
			end_index_nuclei[0] = current_idx[0] + double_scale_nuclei; end_index_nuclei[1] = current_idx[1] + double_scale_nuclei; end_index_nuclei[2] = current_idx[2] + double_scale_nuclei;


			//std::cout << "Nuclei: Starting Indices:"<<starting_index_nuclei[0]<<" "<<starting_index_nuclei[1]<<" "<<starting_index_nuclei[2]<<" "<<std::endl;
			//std::cout << "Nuclei: End Indices:"<<end_index_nuclei[0]<<" "<<end_index_nuclei[1]<<" "<<end_index_nuclei[2]<<std::endl;

			//std::cout << "Nuclei: "  << i << " Scale: " << double_scale_nuclei << std::endl;

			if ( (starting_index_nuclei[0] < 0) || (starting_index_nuclei[1] < 0) || (starting_index_nuclei[2] < 0) ||
				(end_index_nuclei[0] > (unsigned int)sz[0]) || (end_index_nuclei[1] > (unsigned int)sz[1]) ||
				(end_index_nuclei[2] > (unsigned int)sz[2]) )
				continue;

			//std::cout << "Nuclei: "  << i << " Scale: " << double_scale_nuclei << std::endl;

			sub_volume_size_nuclei[0] = 2.0 * double_scale_nuclei; sub_volume_size_nuclei[1] = 2.0 * double_scale_nuclei; sub_volume_size_nuclei[2] = 2.0 * double_scale_nuclei;

			sub_volume_region_nuclei.SetIndex(starting_index_nuclei);
			sub_volume_region_nuclei.SetSize(sub_volume_size_nuclei);


			typedef itk::RegionOfInterestImageFilter<ImageType3D, ImageType3D> VolumeOfInterestFilterType_nuclei2;
			VolumeOfInterestFilterType_nuclei2::Pointer sub_volume_filter_nuclei = VolumeOfInterestFilterType_nuclei2::New();
			sub_volume_filter_nuclei->SetInput(this->SomaDistanceMapImage);
			sub_volume_filter_nuclei->SetRegionOfInterest(sub_volume_region_nuclei);
			sub_volume_filter_nuclei->Update();
			sub_volume_nuclei = sub_volume_filter_nuclei->GetOutput();

			ImageType3D::Pointer distance_map = sub_volume_nuclei;


			double cur_distance;
			double min_distance = 1000.0;
			double max_distance = -1000.0; 
			double mean_distance = 0.0;
			double variance_distance = 0.0;
			double acc_distance = 0.0;
			int n_roots = 0;

			//std::cout << std::endl;

			std::vector<double> distance_array;
			for(size_t j = 0; j < this->CandidateRootPoints.size(); j++){

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

				mean_distance = acc_distance / (double)n_roots;
				for(int k = 0; k < distance_array.size(); k++)
					variance_distance = variance_distance + std::pow(distance_array[k] - mean_distance, 2);

				variance_distance = variance_distance / (double)n_roots;

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
				this->NucleiObjects[i].associativeFeatures.minRootDist = 100000.0; //double_scale_nuclei;
				this->NucleiObjects[i].associativeFeatures.maxRootDist = 100000.0; //double_scale_nuclei;
				this->NucleiObjects[i].associativeFeatures.meanRootDist = 100000.0; //double_scale_nuclei;
				this->NucleiObjects[i].associativeFeatures.varRootDist = 100000.0; //double_scale_nuclei;
				this->NucleiObjects[i].associativeFeatures.nRoots = 0;

			}//loop over candidate roots

			//distance_array.clear();

		}//end of loop over nuclei

		vtkSmartPointer<vtkDoubleArray> min_dist_array = vtkSmartPointer<vtkDoubleArray>::New();
		min_dist_array->SetName("min_root_dist");
		vtkSmartPointer<vtkDoubleArray> max_dist_array = vtkSmartPointer<vtkDoubleArray>::New();
		max_dist_array->SetName("max_root_dist");
		vtkSmartPointer<vtkDoubleArray> mean_dist_array = vtkSmartPointer<vtkDoubleArray>::New();
		mean_dist_array->SetName("mean_root_dist");
		vtkSmartPointer<vtkDoubleArray> var_dist_array = vtkSmartPointer<vtkDoubleArray>::New();
		var_dist_array->SetName("var_root_dist");
		vtkSmartPointer<vtkDoubleArray> n_roots_array = vtkSmartPointer<vtkDoubleArray>::New();
		n_roots_array->SetName("n_roots");

		for(int i = 0; i < this->NucleiObjects.size(); i++){
			min_dist_array->InsertNextValue(this->NucleiObjects[i].associativeFeatures.minRootDist);
			max_dist_array->InsertNextValue(this->NucleiObjects[i].associativeFeatures.maxRootDist);
			mean_dist_array->InsertNextValue(this->NucleiObjects[i].associativeFeatures.meanRootDist);
			var_dist_array->InsertNextValue(this->NucleiObjects[i].associativeFeatures.varRootDist);
			n_roots_array->InsertNextValue(this->NucleiObjects[i].associativeFeatures.nRoots);
		}	


		nuclei_features_table[0]->AddColumn(min_dist_array);
		nuclei_features_table[0]->AddColumn(max_dist_array);
		nuclei_features_table[0]->AddColumn(mean_dist_array);
		nuclei_features_table[0]->AddColumn(var_dist_array);
		nuclei_features_table[0]->AddColumn(n_roots_array);
	}


	vtkSmartPointer<vtkTable> active_model_table = ftk::LoadTable(TrainingFileName);

	// to generate the Active Learning Matrix
	vnl_matrix<double> act_learn_matrix;
	act_learn_matrix.set_size((int)active_model_table->GetNumberOfColumns() , (int)active_model_table->GetNumberOfRows() - 2);
	for(int row = 2; row<(int)active_model_table->GetNumberOfRows(); ++row)
	{
		for(int col=0; col<(int)active_model_table->GetNumberOfColumns(); ++col)
		{
			act_learn_matrix.put(col, row-2, active_model_table->GetValue(row,col).ToDouble());
		}
	}

	//to generate the std_deviation and the mean vectors
	vnl_vector<double> std_dev_vec, mean_vec; 
	std_dev_vec.set_size((int)active_model_table->GetNumberOfColumns() - 1);
	mean_vec.set_size((int)active_model_table->GetNumberOfColumns() - 1);
	for(int col=1; col<(int)active_model_table->GetNumberOfColumns(); ++col)
	{
		std_dev_vec.put(col-1, active_model_table->GetValue(0,col).ToDouble());
		mean_vec.put(col-1, active_model_table->GetValue(1,col).ToDouble());
	}
	active_model_table->RemoveRow(0);
	active_model_table->RemoveRow(0);
	active_model_table->RemoveColumn(0);

	std::string classification_name = "multi_class";
	double confidence_thresh = 0.25;

	MCLR* mclr = new MCLR();
	//Number of features and classes needed in "add_bias" fuction of MCLR
	mclr->Set_Number_Of_Classes((int)active_model_table->GetNumberOfRows());
	mclr->Set_Number_Of_Features((int)active_model_table->GetNumberOfColumns());

	vtkSmartPointer<vtkTable> test_table  = vtkSmartPointer<vtkTable>::New();
	test_table->Initialize();
	//test_table->SetNumberOfRows(roots_table->GetNumberOfRows());
	for(int col=0; col<(int)active_model_table->GetNumberOfColumns(); ++col)
	{
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName(active_model_table->GetColumnName(col));
		test_table->AddColumn(column);	
	}
	for(int row = 0; row < (int)nuclei_features_table[0]->GetNumberOfRows(); ++row)
	{	
		vtkSmartPointer<vtkVariantArray> model_data1 = vtkSmartPointer<vtkVariantArray>::New();
		for(int col=0; col<(int)active_model_table->GetNumberOfColumns();++col)
		{
			std::string column_name = active_model_table->GetColumnName(col);
			model_data1->InsertNextValue(nuclei_features_table[0]->GetValueByName(row,column_name.c_str()));
		}
		test_table->InsertNextRow(model_data1);
	}	

	////// Final Data  to classify from the model
	vnl_matrix<double> data_classify;
	//if(normalize_from_model)
	//	data_classify =  mclr->Normalize_Feature_Matrix_w(mclr->tableToMatrix_w(test_table), std_dev_vec, mean_vec);
	//else
	//	data_classify =  mclr->Normalize_Feature_Matrix(mclr->tableToMatrix_w(test_table));
	data_classify = mclr->Normalize_Feature_Matrix(mclr->tableToMatrix_w(test_table));
	data_classify = data_classify.transpose();

	vnl_matrix<double> currprob;
	currprob = mclr->Test_Current_Model_w(data_classify, act_learn_matrix);

	std::string prediction_col_name = "prediction_active_" + classification_name;
	std::string confidence_col_name = "confidence_" + classification_name;

	//// Add the Prediction Column 
	std::vector< std::string > prediction_column_names = ftk::GetColumsWithString(prediction_col_name, nuclei_features_table[0] );
	if(prediction_column_names.size() == 0)
	{
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName(prediction_col_name.c_str());
		column->SetNumberOfValues( nuclei_features_table[0]->GetNumberOfRows() );
		nuclei_features_table[0]->AddColumn(column);
	}

	// Add the confidence column
	std::vector< std::string > confidence_column_names = ftk::GetColumsWithString(confidence_col_name, nuclei_features_table[0] );
	if(confidence_column_names.size() == 0)
	{
		vtkSmartPointer<vtkDoubleArray> column_confidence = vtkSmartPointer<vtkDoubleArray>::New();
		column_confidence->SetName(confidence_col_name.c_str());
		column_confidence->SetNumberOfValues( nuclei_features_table[0]->GetNumberOfRows() );
		nuclei_features_table[0]->AddColumn(column_confidence);
	}

	for(int row = 0; row<(int)nuclei_features_table[0]->GetNumberOfRows(); ++row)  
	{
		vnl_vector<double> curr_col = currprob.get_column(row);
		nuclei_features_table[0]->SetValueByName(row, confidence_col_name.c_str(), vtkVariant(curr_col(curr_col.arg_max())));
		if(curr_col(curr_col.arg_max()) > confidence_thresh) 
		{
			nuclei_features_table[0]->SetValueByName(row, prediction_col_name.c_str(), vtkVariant(curr_col.arg_max()+1));						
		}
		else
		{
			nuclei_features_table[0]->SetValueByName(row, prediction_col_name.c_str(), vtkVariant(0));
		}
	}
	//roots_table_vec[0] = roots_table;

	delete mclr;


	ftk::SaveTable(outputFname, nuclei_features_table[0]);

}

void AstroTracer::WriteNucleiFeatures(std::string outputFname){
	
	std::ofstream nuclei_feature_vector;
	nuclei_feature_vector.open(outputFname.c_str(), std::ios::out);
	//std::cout << "After feature_vector.open"<<std::endl;
	
	unsigned short IDIndex = 0;//ID index
	
	/*nuclei_feature_vector << "ID" << '\t' << "x" << '\t' << "y" << '\t' << "z" << '\t' << "volume" << '\t' << "sum_int" << '\t' << "mean_int" << '\t';	
	nuclei_feature_vector << "var_int" << '\t' << "eccentricity" << '\t' << "elongation" << '\t' << "mean_surf_gradient" << '\t' << "radius_variation" << '\t'; 
	nuclei_feature_vector << "shape_measure" << '\t' << "energy" << '\t' << "entropy" << '\t' << "inverse_diff_moment" << '\t' << "inertia" << '\t';
	nuclei_feature_vector << "cluster_shade" << '\t' << "cluster_prominence" << '\t' << "Astro_TOTAL" << '\t' << "Astro_AVG" << '\t' << "Astro_SURR" << '\t';
	nuclei_feature_vector << "Micro_TOTAL" << '\t' << "Micro_AVG" << '\t' << "Micro_SURR" << '\t' << "Neuro_TOTAL" << '\t' << "Neuro_AVG" << '\t';
	nuclei_feature_vector << "Neuro_SURR" << '\t'; //<<  "Vessel_TOTAL" << '\t' << "Vessel_AVG" << '\t' << "Vessel_SURR" << '\t';
	nuclei_feature_vector << "min_root_dist" << '\t' << "max_root_dist" << '\t' << "mean_root_dist" << '\t' << "var_root_dist" << '\t' << "n_roots" << '\t';
	nuclei_feature_vector << std::endl;*/

	nuclei_feature_vector << "ID" << '\t' << "x" << '\t' << "y" << '\t' << "z" << '\t' << "volume" << '\t' << "inegrated_intensity" << '\t';
	nuclei_feature_vector << "eccentricity" << '\t' << "elongation" << '\t' << "orientation" << '\t' << "mean_int" << '\t' << "minimum" << '\t' << "maximum" << '\t' << "variance" << '\t';	
	nuclei_feature_vector << "surface_gradient" << '\t' << "interior_gradient" << '\t' << "surface_intensity" << "\t" << "interior_intensity" << '\t' << "intensity_ratio" << '\t' << "convexity" << '\t';
	nuclei_feature_vector << "radius_variation" << '\t' << "surface_area" << '\t' << "shape" << '\t' << "shared_boundary" << '\t' << "t_energy" << '\t' << "t_entropy" << '\t' << "inverse_diff_moment" << '\t' << "inertia" << '\t';
	nuclei_feature_vector << "cluster_shade" << '\t' << "cluster_prominence" << '\t' << "TRI_TOTAL" << '\t' << "TRI_AVG" << '\t' << "TRI_SURR" << '\t';
	nuclei_feature_vector << "GFP_TOTAL" << '\t' << "GFP_AVG" << '\t' << "GFP_SURR" << '\t' << "Cy5_TOTAL" << '\t' << "Cy5_AVG" << '\t';
	nuclei_feature_vector << "Cy5_SURR" << '\t'; //<<  "Vessel_TOTAL" << '\t' << "Vessel_AVG" << '\t' << "Vessel_SURR" << '\t';
	nuclei_feature_vector << "min_root_dist" << '\t' << "max_root_dist" << '\t' << "mean_root_dist" << '\t' << "var_root_dist" << '\t' << "n_roots" << '\t';
	
	//TEMP: DARPA HACK
	nuclei_feature_vector << "p_neu" << '\t' << "c_neu" << '\t' << "distance_to_electrode" << '\t';

	nuclei_feature_vector << std::endl;

	for(size_t i = 0; i < this->NucleiObjects.size(); i++){

		NucleiObject nuc = this->NucleiObjects[i];
		
		//if(nuc.associativeFeatures.nRoots != 0){
			nuclei_feature_vector << nuc.intrinsicFeatures.ID << '\t' << nuc.intrinsicFeatures.centroid.ndx[0] << '\t' << nuc.intrinsicFeatures.centroid.ndx[1] << '\t' << nuc.intrinsicFeatures.centroid.ndx[2] << '\t';
			nuclei_feature_vector << nuc.intrinsicFeatures.volume << '\t' << nuc.intrinsicFeatures.integratedIntensity << '\t' << nuc.intrinsicFeatures.eccentricity << '\t' << nuc.intrinsicFeatures.elongation << '\t';
			nuclei_feature_vector << nuc.intrinsicFeatures.orientation << '\t' << nuc.intrinsicFeatures.meanIntensity << '\t' << nuc.intrinsicFeatures.minIntensity << '\t' << nuc.intrinsicFeatures.maxIntensity << '\t';
			nuclei_feature_vector << nuc.intrinsicFeatures.varianceIntensity << '\t' << nuc.intrinsicFeatures.meanSurfaceGradient << '\t' << nuc.intrinsicFeatures.interiorGradient << '\t' << nuc.intrinsicFeatures.surfaceIntensity << '\t';
			nuclei_feature_vector << nuc.intrinsicFeatures.interiorIntensity << '\t' <<nuc.intrinsicFeatures.intensityRatio << '\t' << nuc.intrinsicFeatures.convexity << '\t' << nuc.intrinsicFeatures.radiusVariation << '\t';
			nuclei_feature_vector << nuc.intrinsicFeatures.surfaceArea << '\t' << nuc.intrinsicFeatures.shapeMeasure << '\t' << nuc.intrinsicFeatures.sharedBoundary << '\t' << nuc.intrinsicFeatures.energy << '\t';
			nuclei_feature_vector << nuc.intrinsicFeatures.entropy << '\t' << nuc.intrinsicFeatures.inverseDiffMoment << '\t' << nuc.intrinsicFeatures.inertia << '\t' << nuc.intrinsicFeatures.clusterShade << '\t' << nuc.intrinsicFeatures.clusterProminence << '\t'; 
			nuclei_feature_vector << nuc.associativeFeatures.astro_total << '\t' << nuc.associativeFeatures.astro_avg << '\t' << nuc.associativeFeatures.astro_surr << '\t';
			nuclei_feature_vector << nuc.associativeFeatures.micro_total << '\t' << nuc.associativeFeatures.micro_avg << '\t' << nuc.associativeFeatures.micro_surr << '\t';
			nuclei_feature_vector << nuc.associativeFeatures.neuro_total << '\t' << nuc.associativeFeatures.neuro_avg << '\t' << nuc.associativeFeatures.neuro_surr << '\t';
			//nuclei_feature_vector << nuc.associativeFeatures.vessel_total << '\t' << nuc.associativeFeatures.vessel_avg << '\t' << nuc.associativeFeatures.vessel_surr << '\t';
			nuclei_feature_vector << nuc.associativeFeatures.minRootDist << '\t' << nuc.associativeFeatures.maxRootDist << '\t' << nuc.associativeFeatures.meanRootDist << '\t';
			nuclei_feature_vector << nuc.associativeFeatures.varRootDist << '\t' << nuc.associativeFeatures.nRoots << '\t';

			//TEMP: DARPA HACK
			nuclei_feature_vector << nuc.classValue << '\t' << nuc.confidenceMeasure << '\t' << nuc.distanceToElectrode << '\t';

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
			nuclei_object.intrinsicFeatures.centroid = HeapNode_astro(idx, 0);

			nuclei_object.intrinsicFeatures.volume = atof(str_vec[4].c_str());
			nuclei_object.intrinsicFeatures.integratedIntensity = atof(str_vec[5].c_str());
			
			nuclei_object.intrinsicFeatures.meanIntensity = atof(str_vec[6].c_str());
			nuclei_object.intrinsicFeatures.varianceIntensity = atof(str_vec[7].c_str());

			nuclei_object.intrinsicFeatures.eccentricity = atof(str_vec[8].c_str());
			nuclei_object.intrinsicFeatures.elongation = atof(str_vec[9].c_str());
			nuclei_object.intrinsicFeatures.boundingBoxVolume = nuclei_object.intrinsicFeatures.volume; //atof(str_vec[9].c_str());

			nuclei_object.intrinsicFeatures.meanSurfaceGradient = atof(str_vec[10].c_str());
			nuclei_object.intrinsicFeatures.radiusVariation = atof(str_vec[11].c_str());
			nuclei_object.intrinsicFeatures.shapeMeasure = atof(str_vec[12].c_str());
			
			nuclei_object.intrinsicFeatures.energy = atof(str_vec[13].c_str());
			nuclei_object.intrinsicFeatures.entropy = atof(str_vec[14].c_str());
			nuclei_object.intrinsicFeatures.inverseDiffMoment = atof(str_vec[15].c_str());
			nuclei_object.intrinsicFeatures.inertia = atof(str_vec[16].c_str());
			nuclei_object.intrinsicFeatures.clusterShade = atof(str_vec[17].c_str());
			nuclei_object.intrinsicFeatures.clusterProminence = atof(str_vec[18].c_str());

			
			nuclei_object.classValue = atof(str_vec[33].c_str());
			nuclei_object.confidenceMeasure = atof(str_vec[34].c_str());


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


void AstroTracer::Classification_Roots(std::vector< vtkSmartPointer<vtkTable> >& roots_table_vec, std::vector< LabelImageType3D::Pointer >& root_images, std::string TrainingFileName, std::string rootsTableName, std::string rootsImageName, const bool writeResult, bool normalize_from_model)
{
	std::cout<<"Classification Started Using MCLR\n";
	if(!features_table)
		return;

	vtkSmartPointer<vtkTable> roots_table = features_table;
	vtkSmartPointer<vtkTable> active_model_table = ftk::LoadTable(TrainingFileName);

	// to generate the Active Learning Matrix
	vnl_matrix<double> act_learn_matrix;
	act_learn_matrix.set_size((int)active_model_table->GetNumberOfColumns() , (int)active_model_table->GetNumberOfRows() - 2);
	for(int row = 2; row<(int)active_model_table->GetNumberOfRows(); ++row)
	{
		for(int col=0; col<(int)active_model_table->GetNumberOfColumns(); ++col)
		{
			act_learn_matrix.put(col, row-2, active_model_table->GetValue(row,col).ToDouble());
		}
	}

	//to generate the std_deviation and the mean vectors
	vnl_vector<double> std_dev_vec, mean_vec; 
	std_dev_vec.set_size((int)active_model_table->GetNumberOfColumns() - 1);
	mean_vec.set_size((int)active_model_table->GetNumberOfColumns() - 1);
	for(int col=1; col<(int)active_model_table->GetNumberOfColumns(); ++col)
	{
		std_dev_vec.put(col-1, active_model_table->GetValue(0,col).ToDouble());
		mean_vec.put(col-1, active_model_table->GetValue(1,col).ToDouble());
	}
	active_model_table->RemoveRow(0);
	active_model_table->RemoveRow(0);
	active_model_table->RemoveColumn(0);

	std::string classification_name = "astro_root";
	double confidence_thresh = 0.5;

	MCLR* mclr = new MCLR();
	//Number of features and classes needed in "add_bias" fuction of MCLR
	mclr->Set_Number_Of_Classes((int)active_model_table->GetNumberOfRows());
	mclr->Set_Number_Of_Features((int)active_model_table->GetNumberOfColumns());

	vtkSmartPointer<vtkTable> test_table  = vtkSmartPointer<vtkTable>::New();
	test_table->Initialize();
	//test_table->SetNumberOfRows(roots_table->GetNumberOfRows());
	for(int col=0; col<(int)active_model_table->GetNumberOfColumns(); ++col)
	{
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName(active_model_table->GetColumnName(col));
		test_table->AddColumn(column);	
	}
	for(int row = 0; row < (int)roots_table->GetNumberOfRows(); ++row)
	{	
		vtkSmartPointer<vtkVariantArray> model_data1 = vtkSmartPointer<vtkVariantArray>::New();
		for(int col=0; col<(int)active_model_table->GetNumberOfColumns();++col)
		{
			std::string column_name = active_model_table->GetColumnName(col);
			model_data1->InsertNextValue(roots_table->GetValueByName(row,column_name.c_str()));
		}
		test_table->InsertNextRow(model_data1);
	}	

	////// Final Data  to classify from the model
	vnl_matrix<double> data_classify;
	//if(normalize_from_model)
	//	data_classify =  mclr->Normalize_Feature_Matrix_w(mclr->tableToMatrix_w(test_table), std_dev_vec, mean_vec);
	//else
	//	data_classify =  mclr->Normalize_Feature_Matrix(mclr->tableToMatrix_w(test_table));
	data_classify = mclr->Normalize_Feature_Matrix(mclr->tableToMatrix_w(test_table));
	data_classify = data_classify.transpose();

	vnl_matrix<double> currprob;
	currprob = mclr->Test_Current_Model_w(data_classify, act_learn_matrix);

	std::string prediction_col_name = "prediction_active_" + classification_name;
	std::string confidence_col_name = "confidence_" + classification_name;

	//// Add the Prediction Column 
	std::vector< std::string > prediction_column_names = ftk::GetColumsWithString(prediction_col_name, roots_table );
	if(prediction_column_names.size() == 0)
	{
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName(prediction_col_name.c_str());
		column->SetNumberOfValues( roots_table->GetNumberOfRows() );
		roots_table->AddColumn(column);
	}

	// Add the confidence column
	std::vector< std::string > confidence_column_names = ftk::GetColumsWithString(confidence_col_name, roots_table );
	if(confidence_column_names.size() == 0)
	{
		vtkSmartPointer<vtkDoubleArray> column_confidence = vtkSmartPointer<vtkDoubleArray>::New();
		column_confidence->SetName(confidence_col_name.c_str());
		column_confidence->SetNumberOfValues( roots_table->GetNumberOfRows() );
		roots_table->AddColumn(column_confidence);
	}

	for(int row = 0; row<(int)roots_table->GetNumberOfRows(); ++row)  
	{
		vnl_vector<double> curr_col = currprob.get_column(row);
		roots_table->SetValueByName(row, confidence_col_name.c_str(), vtkVariant(curr_col(curr_col.arg_max())));
		if(curr_col(curr_col.arg_max()) > confidence_thresh) 
		{
			roots_table->SetValueByName(row, prediction_col_name.c_str(), vtkVariant(curr_col.arg_max()+1));						
		}
		else
		{
			roots_table->SetValueByName(row, prediction_col_name.c_str(), vtkVariant(0));
		}
	}
	roots_table_vec[0] = roots_table;

	delete mclr;

	//roots_Image = IDImage;
	roots_Image = LabelImageType3D::New();
	roots_Image->SetRegions(IDImage->GetBufferedRegion());
	roots_Image->Allocate();
	//roots_Image->SetSpacing(PaddedCurvImage->GetSpacing());
	roots_Image->FillBuffer(0);

	LabelImageType3D::PixelType * rootsImageArray = roots_Image->GetBufferPointer();
	itk::Size<3> im_size = roots_Image->GetBufferedRegion().GetSize();
	int slice_size = im_size[1] * im_size[0];
	int row_size = im_size[0];

	for(int row=0; row<(int)roots_table->GetNumberOfRows(); ++row)
	{
		if(roots_table->GetValueByName(row, prediction_col_name.c_str()).ToInt() != 1)
		{
			//int x_pos = roots_table->GetValue(row,1).ToUnsignedInt();
			//int y_pos = roots_table->GetValue(row,2).ToUnsignedInt();
			//int z_pos = roots_table->GetValue(row,3).ToUnsignedInt();
			//unsigned long offset = (z_pos*slice_size)+(y_pos*row_size)+x_pos;
			//rootsImageArray[offset] = 0;
			roots_table->RemoveRow(row);
			--row;
		}
	}

	for(int row=0; row<(int)roots_table->GetNumberOfRows(); ++row)
	{
		int IDIndex = roots_table->GetValue(row,0).ToInt();
		ImageType3D::IndexType current_idx;
		current_idx[0] = roots_table->GetValue(row,1).ToInt();
		current_idx[1] = roots_table->GetValue(row,2).ToInt();
		current_idx[2] = roots_table->GetValue(row,3).ToInt();
		roots_Image->SetPixel(current_idx, IDIndex);
	}

	root_images[0] = roots_Image;

	if(writeResult)
	{
		ftk::SaveTable(rootsTableName, roots_table);
		itk::ImageFileWriter< LabelImageType3D >::Pointer rootsImageWriter = itk::ImageFileWriter< LabelImageType3D >::New();
		rootsImageWriter->SetFileName(rootsImageName.c_str());
		rootsImageWriter->SetInput(roots_Image);
		try
		{
			rootsImageWriter->Update();
		}
		catch(itk::ExceptionObject e)
		{
			std::cout << e << std::endl;
			//return EXIT_FAILURE;
		}
	}

	std::cout<< "Classification done" << std::endl;
}

