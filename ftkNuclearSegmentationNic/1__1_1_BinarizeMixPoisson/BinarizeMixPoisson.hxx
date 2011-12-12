
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType >
void ftk::nucSecNic::BinarizeMixPoisson< inputPixelType, binaryPixelType >::setParameters( unsigned int numberBins_mixPoisson, bool getResultImg_mixPoisson )
{
	_numberBins_mixPoisson = numberBins_mixPoisson;
	_getResultImg_mixPoisson = getResultImg_mixPoisson;
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType >
void ftk::nucSecNic::BinarizeMixPoisson< inputPixelType, binaryPixelType >::setInput( const ftk::Image::Info* info, typename itk::Image< inputPixelType, 3 >::Pointer inputImage )
{
	_inputImage = inputImage;
	
	_info = info;
	_numRows = _info->numRows;
	_numColumns = _info->numColumns;
	_numStacks = _info->numZSlices;
	
	_totNumPixels = (long long)_numRows*(long long)_numColumns*(long long)_numStacks;
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType >
void ftk::nucSecNic::BinarizeMixPoisson< inputPixelType, binaryPixelType >::runBinarization()
{
	std::cout << std::endl << "Running Binarization using mixture of possion";
	std::cout << std::endl << "	Parameters:	numbins: " << _numberBins_mixPoisson;
	
	std::cout << std::endl << "	Allocating:	" << _numRows*_numColumns*_numStacks*sizeof(binaryPixelType) / (1024.0 * 1024.0) << " MB of memory.";
	_binaryImage = binaryImageType::New();
	typename binaryImageType::IndexType start;
	start[0] = 0;
	start[1] = 0;
	start[2] = 0;
	typename binaryImageType::SizeType size;
	size[0] = _numRows;
	size[1] = _numColumns;
	size[2] = _numStacks;
	typename binaryImageType::RegionType region;
	region.SetIndex( start );
	region.SetSize( size );
	_binaryImage->SetRegions( region );
	_binaryImage->Allocate();
	const typename  binaryImageType::PixelType ceros = 0;
	_binaryImage->FillBuffer( ceros );
	try{
		_binaryImage->Update();
	}
	catch( itk::ExceptionObject & err ){
		std::cerr << "ExceptionObject caught!" << std::endl;
		std::cerr << err << std::endl;
	}
	
	this->runMinErrorThresholding();
	
	
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
	
	
	
	
	
}