
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType, typename seedDetecPixelType, typename loGResponsePixelType >
void ftk::nucSecNic::AlphaExpansion< inputPixelType, binaryPixelType, seedDetecPixelType, loGResponsePixelType >::setParameters( bool getResultImg_AphaExpansion )
{
	_getResultImg_AphaExpansion = getResultImg_AphaExpansion;
	
	_spacing_XY = 3;
	_spacing_Z = 4;
	
}


// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType, typename seedDetecPixelType, typename loGResponsePixelType >
void ftk::nucSecNic::AlphaExpansion< inputPixelType, binaryPixelType, seedDetecPixelType, loGResponsePixelType >::setInput( const ftk::Image::Info* info, typename inputImageType::Pointer inputImage, typename binaryImageType::Pointer binaryImage, typename seedDetectImageType::Pointer seedDetectImage, typename loGResponseImageType::Pointer loGResponseImage, typename maxClustImageType::Pointer maxClustImage, ConnComp myConnComp )
{
	_inputImage = inputImage;
	_inputImageArray = _inputImage->GetBufferPointer();

	_binaryImage = binaryImage;
	_binaryImageArray = _binaryImage->GetBufferPointer();
	
	_seedDetectImage = seedDetectImage;
	_seedDetectImageArray = _seedDetectImage->GetBufferPointer();
	
	_loGResponseImage = loGResponseImage;
	_loGResponseImageArray = _loGResponseImage->GetBufferPointer();
	
	_maxClustImage = maxClustImage;
	_maxClustImageArray = _maxClustImage->GetBufferPointer();
	
	_myConnComp = myConnComp;
	
	_info = info;

 	_numRows = _info->numRows;
	_numColumns = _info->numColumns;
	_numStacks = _info->numZSlices;	
	_totNumPixels = _numRows*_numColumns*_numStacks;
	
	std::cout << std::endl << "BinPOIs Set up and input image of size, Rows: " << _numRows << ", Col: " << _numColumns << ", Slices: " << _numStacks << ", Voxels: " << _totNumPixels;
	
	_maxValueInputPixelType = std::numeric_limits<inputPixelType>::max();
	
	_maxValueBinaryPixelType = std::numeric_limits<binaryPixelType>::max();
	
	_maxValueSeedDetectPixelType = std::numeric_limits<seedDetecPixelType>::max();
	
	_maxValueMaxClustPixelType = std::numeric_limits<maxClustPixelType>::max();
	
	typedef itk::MinimumMaximumImageCalculator < loGResponseImageType > ImageCalculatorFilterType;
	typename ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New ();
	imageCalculatorFilter->SetImage(_loGResponseImage);
	imageCalculatorFilter->Compute();
	
	_minLoGResponse = imageCalculatorFilter->GetMinimum();
	
	
	
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType, typename seedDetecPixelType, typename loGResponsePixelType >
void ftk::nucSecNic::AlphaExpansion< inputPixelType, binaryPixelType, seedDetecPixelType, loGResponsePixelType >::runAlphaExpansion()
{
	int maxNumColors = 0;
	
	
	
	// is there any problem with this ?
	if(_minLoGResponse<=0)
	{
		_minLoGResponse = -_minLoGResponse;
		for(int i=0; i<_totNumPixels; ++i)
		{
			_loGResponseImageArray[i] += _minLoGResponse+1;
		}
	}
	
	
	
}




// // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// template < typename inputPixelType, typename binaryPixelType, typename seedDetecPixelType, typename loGResponsePixelType >
// typename itk::Image< seedDetecPixelType, 3 >::Pointer ftk::nucSecNic::AlphaExpansion< inputPixelType, binaryPixelType, seedDetecPixelType, loGResponsePixelType >::getMaxClustImage()
// {
// 	return _maxClustImage;
// }



// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType, typename seedDetecPixelType, typename loGResponsePixelType >
double ftk::nucSecNic::AlphaExpansion< inputPixelType, binaryPixelType, seedDetecPixelType, loGResponsePixelType >::MultivaNormal_3D( double X, double Y, double Z, double Ux, double Uy, double Uz, double S00, double S01, double S11, double S02, double S12, double S22 )
{
// 	//A 3-D multivariate gaussian
// 	double Multivar_Norm( double X, double Y, double Z, double Ux, double Uy, double Uz, double S00, double S01, double S11, double S02, double S12, double S22)
// 	{
	double det_segma = S00*(S11*S22-S12*S12) - S01*(S01*S22-S12*S02) + S02*(S01*S12-S11*S02);
	double Sinv00 = (S22*S11-S12*S12)/det_segma;
	double Sinv01 = -(S22*S01-S12*S02)/det_segma;
	double Sinv02 = (S12*S01-S11*S02)/det_segma;
	double Sinv11 = (S22*S00-S02*S02)/det_segma;
	double Sinv12 = -(S12*S00-S01*S02)/det_segma;
	double Sinv22 = (S11*S00-S01*S01)/det_segma;
	X = X-Ux;
	Y = Y-Uy;
	Z = Z-Uz;
	double Mah_Dist = X*(X*Sinv00+Y*Sinv01+Z*Sinv02) + Y*(X*Sinv01+Y*Sinv11+Z*Sinv12) + Z*(X*Sinv02+Y*Sinv12+Z*Sinv22);
	double Ex = std::exp(-.5*Mah_Dist);
	double A = 1/(std::pow((2*3.1416),1.5)*sqrt(det_segma));
	
	return (A*Ex);
}