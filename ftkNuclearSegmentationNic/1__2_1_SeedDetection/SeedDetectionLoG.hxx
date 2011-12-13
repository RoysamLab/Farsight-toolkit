
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType, typename seedDetecPixelType >
void ftk::nucSecNic::SeedDetectionLoG< inputPixelType, binaryPixelType, seedDetecPixelType >::setParameters( long long seedDetectMinScale, long long seedDetectMaxScale, bool getResultImg_seedDetect )
{
	_seedDetectMinScale = seedDetectMinScale;
	_seedDetectMaxScale = seedDetectMaxScale;
	_getResultImg_seedDetect = getResultImg_seedDetect;
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType, typename seedDetecPixelType >
void ftk::nucSecNic::SeedDetectionLoG< inputPixelType, binaryPixelType, seedDetecPixelType >::setInput( const ftk::Image::Info* info, typename itk::Image< inputPixelType, 3 >::Pointer inputImage )
{
	_inputImage = inputImage;
	_inputImageArray = _inputImage->GetBufferPointer();
	
	_info = info;

 	_numRows = _info->numRows;
	_numColumns = _info->numColumns;
	_numStacks = _info->numZSlices;	
	_totNumPixels = _numRows*_numColumns*_numStacks;
	
	std::cout << std::endl << "BinPOIs Set up and input image of size, Rows: " << _numRows << ", Col: " << _numColumns << ", Slices: " << _numStacks << ", Voxels: " << _totNumPixels;
	
	_maxValueInputPixelType = std::numeric_limits<inputPixelType>::max();
	
	_maxValueBinaryPixelType = std::numeric_limits<binaryPixelType>::max();
	
	_maxValueSeedDetectPixelType = std::numeric_limits<seedDetecPixelType>::max();
	
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType, typename seedDetecPixelType >
void ftk::nucSecNic::SeedDetectionLoG< inputPixelType, binaryPixelType, seedDetecPixelType >::runSeedDetection()
{

	
}




// // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// template < typename inputPixelType, typename binaryPixelType, typename seedDetecPixelType >
// typename itk::Image< binaryPixelType, 3 >::Pointer ftk::nucSecNic::SeedDetectionLoG< inputPixelType, binaryPixelType, seedDetecPixelType >::getBinarizedImage()
// {
// 	return _binaryImage;
// }

