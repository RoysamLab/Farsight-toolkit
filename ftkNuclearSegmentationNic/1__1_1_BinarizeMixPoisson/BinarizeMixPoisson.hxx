
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType >
void ftk::nucSecNic::BinarizeMixPoisson< inputPixelType, binaryPixelType >::setParameters( unsigned int numberBins_mixPoisson, bool getResultImg_mixPoisson )
{
	_numberBins_mixPoisson = numberBins_mixPoisson;
	_getResultImg_mixPoisson = getResultImg_mixPoisson;
	_sigmaNeighCost = 100; // !! This should be changed
	_wNeigh = 15;
	
	_spacing_XY = 1;
	_spacing_Z = 2;
	
	_connComponentsConnectivity = 26;
	_minObjectSize = 10;
	
	
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType >
void ftk::nucSecNic::BinarizeMixPoisson< inputPixelType, binaryPixelType >::setInput( const ftk::Image::Info* info, typename itk::Image< inputPixelType, 3 >::Pointer inputImage )
{
	_inputImage = inputImage;
	_inputImageArray = _inputImage->GetBufferPointer();
	
	_info = info;

 	_numRows = _info->numRows;
	_numColumns = _info->numColumns;
	_numStacks = _info->numZSlices;	
	_totNumPixels = (long long)_numRows*(long long)_numColumns*(long long)_numStacks;
	
	std::cout << std::endl << "BinPOIs Set up and input image of size, Rows: " << _numRows << ", Col: " << _numColumns << ", Slices: " << _numStacks << ", Voxels: " << _totNumPixels;
	
	_maxValueInputPixelType = std::numeric_limits<inputPixelType>::max();
	
	_maxValueBinaryPixelType = std::numeric_limits<binaryPixelType>::max();
	
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType >
void ftk::nucSecNic::BinarizeMixPoisson< inputPixelType, binaryPixelType >::runBinarization()
{
	std::cout << std::endl << "Running Binarization using mixture of possion";
	std::cout << std::endl << "	Parameters:	numbins: " << _numberBins_mixPoisson;
	
	std::cout << std::endl << "	Allocating:	" << _totNumPixels << " " << _totNumPixels*sizeof(binaryPixelType) / (1024.0 * 1024.0) << " MB of memory.";
	_binaryImage = binaryImageType::New();
	typename binaryImageType::IndexType start;
	start[0] = 0;
	start[1] = 0;
	start[2] = 0;
	typename binaryImageType::SizeType size;
	size[0] = _numColumns;
	size[1] = _numRows;
	size[2] = _numStacks;
	typename binaryImageType::RegionType region;
	region.SetIndex( start );
	region.SetSize( size );
	
	double spacing[3];
	spacing[0] = 1; //spacing along x
	spacing[1] = 1; //spacing along y
	spacing[2] = _spacing_Z; //spacing along z

	_binaryImage->SetRegions( region );
	_binaryImage->SetSpacing( spacing );
	
	try{
		_binaryImage->Allocate();
	}
	catch( itk::ExceptionObject & err ){
		std::cerr << "ExceptionObject caught!" << std::endl;
		std::cerr << err << std::endl;
	}
	const typename  binaryImageType::PixelType ceros = 0;
	_binaryImage->FillBuffer( ceros );
	_binaryImage->Update();
	
	_binaryImageArray = _binaryImage->GetBufferPointer();
	
	
	// 1. Run Binarization
	this->runMinErrorThresholding();
	
	// 2. Refine binarization using graph cuts
	graphCuts_3D();
// 	if( _numPoissonDist == 2 )
// 	{
// 		Seg_GC_Full_2D(imgIn, R, C, alpha_F, alpha_B, P_I, &n_nodes, &n_edges, imgOut);
// 	}
	
	
	
	// 3. write result
	typedef  itk::ImageFileWriter< binaryImageType > WriterType;
	typename WriterType::Pointer writer = WriterType::New();
	writer->SetFileName("/media/sf_11_SharedWithUbuntu/Results/imageBinarizedPoisson.mhd");
	writer->SetInput(_binaryImage);
	try{
		writer->Update();
	}
	catch( itk::ExceptionObject & err ){
		std::cerr << "ExceptionObject caught!" << std::endl;
		std::cerr << err << std::endl;
	}
	
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType >
void ftk::nucSecNic::BinarizeMixPoisson< inputPixelType, binaryPixelType >::runMinErrorThresholding()
{
	// Create a normalized histogram
	double histoGram[_numberBins_mixPoisson];
	for( unsigned int i=0; i<_numberBins_mixPoisson; ++i )
	{
		histoGram[i] = 0.0;
	}
	
	// Create histogram
	double factorBins = ((double)_maxValueInputPixelType+1)/((double)_numberBins_mixPoisson); // +1 because is unsigned char (255), unsigned short (65535), for other data types have to be different (unsigned in general +1)
	unsigned int binIndex;
	for( unsigned int i=0; i<_totNumPixels; ++i )
	{
		binIndex = floor(((double)_inputImageArray[i]/factorBins));
		if( binIndex >= _numberBins_mixPoisson )
		{
			binIndex = _numberBins_mixPoisson-1;
		}
		histoGram[binIndex]++;
	}
	
	// Normalize histogram
	for( unsigned int i=0; i<_numberBins_mixPoisson; ++i )
	{
		histoGram[i] /= _totNumPixels;
	}
	
	
	
	//The three-level min error thresholding algorithm
	float P0, U0, P1, U1, P2, U2, U, J, min_J;
	min_J = 1000000.0;

// 	// Try this: we need to define a penalty term that depends on the number of parameters
// 	//The penalty term is given as 0.5*k*ln(n) where k is the number of parameters of the model and n is the number of samples
// 	//In this case, k=6 and n=_numberBins_mixPoisson
// 	double PenTerm3 =  sqrt(6.0)*log((double)_numberBins_mixPoisson);
// 	
// 	for(int i=0; i<_numberBins_mixPoisson-1; ++i)//to set the first threshold
// 	{
// 		//compute the current parameters of the first component
// 		P0 = U0 = 0.0;		
// 		for(int l=0; l<=i; ++l)
// 		{
// 			P0+=histoGram[l];
// 			U0+=(l+1)*histoGram[l];
// 		}
// 		U0 /= P0;
// 
// 		for(int j=i+1; j<_numberBins_mixPoisson; j++)//to set the second threshold
// 		{
// 			//compute the current parameters of the second component
// 			P1 = U1 = 0.0;		
// 			for(int l=i+1; l<=j; l++)
// 			{
// 				P1+=histoGram[l];
// 				U1+=(l+1)*histoGram[l];
// 			}
// 			U1 /= P1;
// 
// 			//compute the current parameters of the third component
// 			P2 = U2 = 0.0;		
// 			for(int l=j+1; l<=_numberBins_mixPoisson; l++)
// 			{
// 				P2+=histoGram[l];
// 				U2+=(l+1)*histoGram[l];
// 			}
// 			U2 /= P2;
// 
// 			//compute the overall mean
// 			U = P0*U0 + P1*U1 + P2*U2;
// 
// 			//Compute the current value of the error criterion function
// 			J =  U - (P0*(log(P0)+U0*log(U0))+ P1*(log(P1)+U1*log(U1)) + P2*(log(P2)+U2*log(U2)));
// 			//Add the penalty term
// 			J +=PenTerm3;
// 
// 			if(J<min_J)
// 			{
// 				min_J = J;
// 				_alpha_1 = U0;
// 				_P_I = P0;
// 				_alpha_2 = U1;
// 				_P_I2 = P1;
// 				_alpha_3 = U2;
// 				_numPoissonDist = 3;
// 			}
// 		}
// 	}
	
	
	//try this: see if using two components is better
	//The penalty term is given as sqrt(k)*ln(n)	
	//In this case, k=4 and n=numberBins_mixPoisson
// 	double PenTerm2 =  2*log((double)_numberBins_mixPoisson);
	
	for(int i=0; i<_numberBins_mixPoisson-1; ++i)//to set the first threshold
	{
		//compute the current parameters of the first component
		P0 = U0 = 0.0;		
		for(int l=0; l<=i; ++l)
		{
			P0+=histoGram[l];
			U0+=(l+1)*histoGram[l];
		}
		U0 /= P0;

		P1 = U1 = 0.0;		
		for(int j=i+1; j<_numberBins_mixPoisson; ++j)//to set the second threshold
		{
			//compute the current parameters of the second component
			
// 			for(int l=j; l<=_numberBins_mixPoisson; ++l)
// 			{
				P1+=histoGram[j];
				U1+=(j+1)*histoGram[j];
// 			}
			
		}
		U1 /= P1;

		//compute the overall mean
		U = P0*U0 + P1*U1;

		//Compute the current value of the error criterion function
		J =  U - (P0*(log(P0)+U0*log(U0))+ P1*(log(P1)+U1*log(U1)));
		//Add the penalty term
// 		J +=PenTerm2;
		if(J<min_J)
		{
			min_J = J;
			_alpha_1 = U0;
			_P_I = P0;
			_alpha_2 = U1;
			_P_I2 = P1;
			_alpha_3 = U2; //Just a negative number to let the program knows that two levels will be used		
			_numPoissonDist = 2;
		}
	}
	
	
	
	
	
	
	
	
	
	std::ofstream myfile;
	//myfile.open ("/home/nicolasreyv/farsight/src/farsight-src/ftkNuclearSegmentationNic/Results/out",ios::out); ask kedar
	myfile.open ("/media/sf_11_SharedWithUbuntu/Results/out_Histogram",ios::out);

	if (!myfile) {
		std::cerr << "Can't open output file " << std::endl;
		std::exit(1);
	}
	for( unsigned int i=0; i<_numberBins_mixPoisson; ++i )
	{
		myfile << histoGram[i] << "\n";
	}
	myfile.close();
	
	
	
	
	
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType >
double ftk::nucSecNic::BinarizeMixPoisson< inputPixelType, binaryPixelType >::computePoissonProb(int intensity, double alpha)
{
	/*this is the equation
	P = (alpha^intensity)*exp(-alpha)/factorial(intensity);
	however, since alpha and the intensity could be large, computing P in that
	way will result in infinity values from (alpha^intensity) and
	factorial(intensity) as a result of matlab's limitations of the data types*/

	//here is the solution
	double A, P;
	A = exp(-alpha);
	P = 1;
	for (int i=1; i<= intensity; ++i)
	{
		P = P * (alpha/i);
	}

	P = P*A;

	if(P < std::numeric_limits<long double>::epsilon())
		P = std::numeric_limits<long double>::epsilon();

	return P;
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType >
void ftk::nucSecNic::BinarizeMixPoisson< inputPixelType, binaryPixelType >::graphCuts_3D()
{    
	
	
	int num_nodes;
	int num_edges;
	
	long long curr_node_z;
	long long curr_node_yz;
	
	long long curr_node_xy;
	long long curr_node_xyz;
	
	long long rght_node;
	long long down_node;
	long long diag_node;
	
	long long diag_node_1;
	long long rght_node_1;
	long long down_node_1;
	long long diag_node_2;
	
	
	double Df;
	double Db;
	double Dr;
	double Dd; 
	double Dg; 
	double F_H[_maxValueInputPixelType+1];
	double B_H[_maxValueInputPixelType+1];
	
	
	typedef Graph<int,int,int> GraphType;
	// No much diference in the result when using double, the memory increases significatively. 
// 	typedef Graph<double,double,double> GraphType;

// 	Set the number of edges and the number of nodes
	num_nodes = _numRows*_numColumns*_numStacks;
	if( _numStacks == 1 )
	{
		num_edges = (_numRows-1)*(_numColumns-1)*3;
	}
	else
	{
		num_edges = (_numRows-1)*(_numColumns-1)*(_numStacks-1)*3;
		// In case of using 26 neighborhood
// 		num_edges = (_numRows-1)*(_numColumns-1)*(_numStacks-1)*7;
	}


	//Before entering the loop, compute the poisson probs 
	for(int i=0; i<=_maxValueInputPixelType; ++i)
	{
		if(i>=_alpha_1)
			F_H[i] = (1-_P_I)*computePoissonProb((int)_alpha_1,_alpha_1);
		else
			F_H[i] = (1-_P_I)*computePoissonProb(i,_alpha_1);
		if(i<=_alpha_2)
			B_H[i] = _P_I*computePoissonProb(int(_alpha_2),_alpha_2);
		else
			B_H[i] = _P_I*computePoissonProb(i,_alpha_2);
	}
	
	std::cout << std::endl << _P_I;
	std::cout << std::endl << _alpha_1;
	std::cout << std::endl << _alpha_2;
	
	
	std::ofstream myfile;
	myfile.open ("/media/sf_11_SharedWithUbuntu/Results/out_Histogram_F_H",ios::out);

	if (!myfile) {
		std::cerr << "Can't open output file " << std::endl;
		std::exit(1);
	}
	for(int i=0; i<=_maxValueInputPixelType; ++i)
	{
		myfile << F_H[i] << "\n";
	}
	myfile.close();
	
	std::ofstream myfile2;
	myfile2.open ("/media/sf_11_SharedWithUbuntu/Results/out_Histogram_B_H",ios::out);

	if (!myfile2) {
		std::cerr << "Can't open output file " << std::endl;
		std::exit(1);
	}
	for(int i=0; i<=_maxValueInputPixelType; ++i)
	{
		myfile2 << B_H[i] << "\n";
	}
	myfile2.close();
	
	
// 	std::cout << std::endl << "	step 4444235";
// 	std::cout << std::endl << "	step 4444235";
// 	std::cout << std::endl << "	step 4444235";
	
	int nic;
// 	std::cin >> nic;

// 	//Here is the main loop.. 
// 	//For each point, compute the terminal and neighbor edge weights
	
	std::cout << std::endl << "	Allocating:	" << num_nodes << " " << num_nodes*sizeof(int) / (1024.0 * 1024.0) << " MB of memory, for nodes.";
	std::cout << std::endl << "	And, allocating:	" << num_edges << " " << num_edges*sizeof(int) / (1024.0 * 1024.0) << " MB of memory, for edges.";
	
	GraphType *g = new GraphType(/*estimated # of nodes*/ num_nodes, /*estimated # of edges*/ num_edges); 
	
	//std::cin >> nic;
	
// 	try{
// 		g = new GraphType(/*estimated # of nodes*/ num_nodes, /*estimated # of edges*/ num_edges); 
// 	}
// 	catch (std::exception& e)
// 	{
// 		cerr << "exception caught: " << e.what() << endl;
// 	}
	
// 	std::cout << std::endl << "	step 111111111";
// 	std::cout << std::endl << "	step 111111111";
// 	std::cout << std::endl << "	step 111111111";
	
	
	int intenPixel;
	for( int k=0; k<_numStacks; ++k )
	{
		curr_node_z = k*_numColumns*_numRows;
		for(int i=0; i<_numRows; ++i)
		{
			curr_node_yz = (i*_numColumns)+curr_node_z;
			for(int j=0; j<_numColumns; ++j)
			{	
				curr_node_xyz = j + curr_node_yz;
			 
				intenPixel = (int)_inputImageArray[curr_node_xyz];
				
				// This hard contraints are not a good idea, improve this !
				if(intenPixel == _maxValueInputPixelType)
				{
					Df = 0;
					Db = 1000;
				}
				else if(intenPixel == 0)
				{
					Df = 1000;
					Db = 0;
				}
				else
				{                
					Df = -log(F_H[intenPixel]);
					if(Df>1000.0)
						Df = 1000;
					Db = -log(B_H[intenPixel]);
					if(Db>1000.0)
						Db=1000;
				}         
				
				g -> add_node();
				g -> add_tweights( curr_node_xyz,   /* capacities */ Df,Db);  
				
// 				std::cout << std::endl << curr_node_xyz << " " << rght_node << " " << down_node << " " << diag_node;
// 				std::cout << std::endl << curr_node_xyz ;
			}
		}
	}

// 	std::cout << std::endl << "	step 444";
// 	std::cout << std::endl << "	step 444";
// 	std::cout << std::endl << "	step 444";

	//std::cin >> nic;
	
// 	double w = 30.0;
	for( int k=0; k<_numStacks-1; ++k )
	{
		curr_node_z = k*_numColumns*_numRows;
		for(int i=0; i<_numRows-1; ++i)
		{
			curr_node_yz = (i*_numColumns)+curr_node_z;
			for(int j=0; j<_numColumns-1; ++j)
			{	
				curr_node_xyz = j + curr_node_yz+1;
				
				rght_node = curr_node_xyz+1;
				down_node = curr_node_xyz+_numColumns;
				diag_node = curr_node_xyz+_numColumns*_numRows;
				


				Dr = _wNeigh*exp(-pow((double)_inputImageArray[curr_node_xyz]-(double)_inputImageArray[rght_node],2)/(2*pow(_sigmaNeighCost,2)));
				
// 				std::cout << std::endl << curr_node_xyz << " " << rght_node << " " << down_node << " " << diag_node << " " << Dr;
// 				std::cout << std::endl << curr_node_xyz << " " << rght_node << " " << down_node << " " << diag_node << " " << Dr;
				
				g->add_edge( curr_node_xyz, rght_node,    /* capacities */  Dr, Dr );	
				
// 				std::cout << std::endl << curr_node_xyz << " " << rght_node << " " << down_node << " " << diag_node;
// 				std::cout << std::endl << curr_node_xyz << " " << rght_node << " " << down_node << " " << diag_node;
				
				Dd = _wNeigh*exp(-pow((double)_inputImageArray[curr_node_xyz]-(double)_inputImageArray[down_node],2)/(2*pow(_sigmaNeighCost,2)));
				g->add_edge( curr_node_xyz, down_node,    /* capacities */  Dd, Dd );
				
				Dg = _wNeigh*exp(-pow((double)_inputImageArray[curr_node_xyz]-(double)_inputImageArray[diag_node],2)/(2*pow(_sigmaNeighCost,2)));
				g->add_edge( curr_node_xyz, diag_node,    /* capacities */  Dg, Dg );  
				
				
				// In case of using 26 neighborhood system
// 				diag_node_1 = down_node+1;
// 				rght_node_1 = rght_node+_numColumns*_numRows;
// 				down_node_1 = down_node+_numColumns*_numRows;
// 				diag_node_2 = diag_node_1+_numColumns*_numRows;
// 				
// 				Dr = _wNeigh*exp(-pow((double)_inputImageArray[curr_node_xyz]-(double)_inputImageArray[diag_node_1],2)/(2*pow(_sigmaNeighCost,2)));
// 				g->add_edge( curr_node_xyz, diag_node_1,    /* capacities */  Dr, Dr );	
// 				
// 				Dr = _wNeigh*exp(-pow((double)_inputImageArray[curr_node_xyz]-(double)_inputImageArray[rght_node_1],2)/(2*pow(_sigmaNeighCost,2)));
// 				g->add_edge( curr_node_xyz, rght_node_1,    /* capacities */  Dr, Dr );	
// 				
// 				Dr = _wNeigh*exp(-pow((double)_inputImageArray[curr_node_xyz]-(double)_inputImageArray[down_node_1],2)/(2*pow(_sigmaNeighCost,2)));
// 				g->add_edge( curr_node_xyz, down_node_1,    /* capacities */  Dr, Dr );	
// 				
// 				Dr = _wNeigh*exp(-pow((double)_inputImageArray[curr_node_xyz]-(double)_inputImageArray[diag_node_2],2)/(2*pow(_sigmaNeighCost,2)));
// 				g->add_edge( curr_node_xyz, diag_node_2,    /* capacities */  Dr, Dr );	
			}
		}
	} 

// 	std::cout << std::endl << "	step 55";
// 	std::cout << std::endl << "	step 55";
// 	std::cout << std::endl << "	step 55";

	//Compute the maximum flow:
	g->maxflow();

	//std::cin >> nic;
	
// 	std::cout << std::endl << "	step 66";
// 	std::cout << std::endl << "	step 66";
// 	std::cout << std::endl << "	step 66";

	int RR,CC;
	for(int i=0; i<num_nodes; ++i)
	{
		CC = ((long)i)%_numColumns;
		RR = (i-CC)/_numColumns;
		if(g->what_segment(i) == GraphType::SOURCE)
		{
			_binaryImageArray[RR*_numColumns + CC] = 0;
		}
		else
		{
			_binaryImageArray[RR*_numColumns + CC] = _maxValueBinaryPixelType;
		}
	}
// 	std::cin >> nic;
 	delete g;	
// 	std::cin >> nic;
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType >
typename itk::Image< binaryPixelType, 3 >::Pointer ftk::nucSecNic::BinarizeMixPoisson< inputPixelType, binaryPixelType >::getBinarizedImage()
{
	return _binaryImage;
}




// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType >
ftk::nucSecNic::ConnComp* ftk::nucSecNic::BinarizeMixPoisson< inputPixelType, binaryPixelType >::getConnComponents()
{
	typedef itk::ConnectedComponentImageFilter< binaryImageType, binaryImageType > FilterType;
	typename FilterType::Pointer filter = FilterType::New();
	typedef itk::RelabelComponentImageFilter< binaryImageType, binaryImageType > RelabelType;
	typename RelabelType::Pointer relabel = RelabelType::New();

	filter->SetInput ( _binaryImage );
	filter->SetFullyConnected( _connComponentsConnectivity );
	relabel->SetInput( filter->GetOutput() );
	relabel->SetMinimumObjectSize( _minObjectSize );

	try
	{
		relabel->Update();
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Relabel: exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}
	
	_numConnectedComponents = relabel->GetNumberOfObjects();
	
	_binaryImage = relabel->GetOutput();
	
	
	
	
	_binaryImageArray = _binaryImage->GetBufferPointer();
	
	int val; 
	int curNode;
	
// 	ConnComp* _myConnComp
	
	_myConnComp = new ConnComp[_numConnectedComponents];
	for( int i=0; i<_numConnectedComponents; ++i)
	{
		_myConnComp[i].x1 =  _numColumns+1;
		_myConnComp[i].y1 =  _numRows+1;
		_myConnComp[i].x2 =  _myConnComp[i].y2 = -1;
		if (_numStacks == 1)
		{
			_myConnComp[i].z1 = 1;
			_myConnComp[i].z2 = 1;
		}
		else
		{
			_myConnComp[i].z1 = _numStacks+1;
			_myConnComp[i].z2 = -1;
		}
	}
	if (_numStacks == 1)
	{
		for ( int j=0; j<_numRows; ++j )
		{
			for ( int i=0; i<_numColumns; ++i )
			{
				curNode = (j*_numColumns)+i;
				val = _binaryImageArray[curNode];
				if(val>0)
				{
					val = val-1;
					_myConnComp[val].y1 = (int) std::min((double)j,(double)_myConnComp[val].y1);
					_myConnComp[val].y2 = (int) std::max((double)j,(double)_myConnComp[val].y2);
					_myConnComp[val].x1 = (int) std::min((double)i,(double)_myConnComp[val].x1);
					_myConnComp[val].x2 = (int) std::max((double)i,(double)_myConnComp[val].x2);
				}
			}
		}
	}
	else
	{
		for (int k=0; k<_numStacks; ++k)
		{
			for (int j=0; j<_numRows; ++j)
			{
				for (int i=0; i<_numColumns; ++i)
				{
					curNode = (k*_numRows*_numColumns)+(j*_numColumns)+i;
					val = _binaryImageArray[curNode];
					if(val>0)
					{
						val = val-1;
						_myConnComp[val].y1 = (int) std::min((double)j,(double)_myConnComp[val].y1);
						_myConnComp[val].y2 = (int) std::max((double)j,(double)_myConnComp[val].y2);
						_myConnComp[val].x1 = (int) std::min((double)i,(double)_myConnComp[val].x1);
						_myConnComp[val].x2 = (int) std::max((double)i,(double)_myConnComp[val].x2);
						_myConnComp[val].z1 = (int) std::min((double)k,(double)_myConnComp[val].z1);
						_myConnComp[val].z2 = (int) std::max((double)k,(double)_myConnComp[val].z2);
					}
				}
			}
		}
	}
	
	std::cout << std::endl << "The number of connected components is: " << _numConnectedComponents;
	
}








