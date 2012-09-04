
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType, typename seedDetecPixelType, typename loGResponsePixelType >
void ftk::nucSecNic::SeedDetectionLoG< inputPixelType, binaryPixelType, seedDetecPixelType, loGResponsePixelType >::setParameters( long long seedDetectMinScale, long long seedDetectMaxScale, bool getResultImg_seedDetect )
{
	_seedDetectMinScale = seedDetectMinScale;
	_seedDetectMaxScale = seedDetectMaxScale;
	_getResultImg_seedDetect = getResultImg_seedDetect;
	
	_spacing_XY = 3;
	_spacing_Z = 4;
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType, typename seedDetecPixelType, typename loGResponsePixelType >
void ftk::nucSecNic::SeedDetectionLoG< inputPixelType, binaryPixelType, seedDetecPixelType, loGResponsePixelType >::setInput( const ftk::Image::Info* info, typename itk::Image< inputPixelType, 3 >::Pointer inputImage, typename itk::Image< binaryPixelType, 3 >::Pointer binaryImage )
{
	_inputImage = inputImage;
	_inputImageArray = _inputImage->GetBufferPointer();

	_binaryImage = binaryImage;
	_binaryImageArray = _binaryImage->GetBufferPointer();
	
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
template < typename inputPixelType, typename binaryPixelType, typename seedDetecPixelType, typename loGResponsePixelType >
void ftk::nucSecNic::SeedDetectionLoG< inputPixelType, binaryPixelType, seedDetecPixelType, loGResponsePixelType >::runSeedDetection()
{
	
	_seedDetectImage = seedDetectImageType::New();
	typename seedDetectImageType::IndexType start;
	start[0] = 0;
	start[1] = 0;
	start[2] = 0;
	typename seedDetectImageType::SizeType size;
	size[0] = _numColumns;
	size[1] = _numRows;
	size[2] = _numStacks;
	typename seedDetectImageType::RegionType region;
	region.SetIndex( start );
	region.SetSize( size );
	
	double spacing[3];
	spacing[0] = 1; //spacing along x
	spacing[1] = 1; //spacing along y
	spacing[2] = _spacing_Z; //spacing along z

	_seedDetectImage->SetRegions( region );
	_seedDetectImage->SetSpacing( spacing );
	
	try{
		_seedDetectImage->Allocate();
	}
	catch( itk::ExceptionObject & err ){
		std::cerr << "ExceptionObject caught!" << std::endl;
		std::cerr << err << std::endl;
	}
	const typename  seedDetectImageType::PixelType ceros = 0;
	_seedDetectImage->FillBuffer( ceros );
	_seedDetectImage->Update();
	
	_seedDetectImageArray = _seedDetectImage->GetBufferPointer();
	
	
	
	
	
	if( _numStacks != 1 )
	{
		distMap_2D( );
	}
	
	std::cout << std::endl << "Running LoG Response";
	std::cout << std::endl << "	Parameters:	minScale: " << _seedDetectMinScale;
	std::cout << std::endl << "			maxScale: " << _seedDetectMaxScale;
	
	std::cout << std::endl << "	Allocating:	" << _totNumPixels << " " << 2*_totNumPixels*sizeof(loGResponsePixelType) / (1024.0 * 1024.0) << " MB of memory.";
	
	_inputLoGImage = inputLoGImageType::New();
	typename inputLoGImageType::IndexType start_3;
	start_3[0] = 0;
	start_3[1] = 0;
	start_3[2] = 0;
	typename inputLoGImageType::SizeType size_3;
	size_3[0] = _numColumns;
	size_3[1] = _numRows;
	size_3[2] = _numStacks;
	typename inputLoGImageType::RegionType region_3;
	region_3.SetIndex( start_3 );
	region_3.SetSize( size_3 );
	_inputLoGImage->SetRegions( region_3 );
	
	
	try{
		_inputLoGImage->Allocate();
	}
	catch( itk::ExceptionObject & err ){
		std::cerr << "ExceptionObject caught!" << std::endl;
		std::cerr << err << std::endl;
	}
	const typename  inputLoGImageType::PixelType ceros_3 = 0;
	_inputLoGImage->FillBuffer( ceros_3 );
	_inputLoGImage->Update();
	
// 	
// 	_loGResponseImageMax = loGResponseImageType::New();
// 	_loGResponseImageMax->SetRegions( region );
// 	
// 	try{
// 		_loGResponseImageMax->Allocate();
// 	}
// 	catch( itk::ExceptionObject & err ){
// 		std::cerr << "ExceptionObject caught!" << std::endl;
// 		std::cerr << err << std::endl;
// 	}
// 	_loGResponseImageMax->FillBuffer( ceros );
// 	_loGResponseImageMax->Update();
// 
// 	
// 	if( _seedDetectMinScale > _seedDetectMaxScale )
// 	{
// 		exit(1);
// 	}
	
	
	
// 	_loGResponseImageMaxArray = _loGResponseImageMax->GetBufferPointer();
// 	_loGResponseImageArray = _loGResponseImage->GetBufferPointer();
	
	_distMapImageArray = _distMapImage->GetBufferPointer();
	_inputLoGImageArray = _inputLoGImage->GetBufferPointer();
	_inputImageArray = _inputImage->GetBufferPointer();
	double mult = 1.0;
	// Apply Yousef Transformation
	for( int j = 0; j<_totNumPixels; ++j )
	{
		mult = 1 + _distMapImageArray[j]/(2*_maxMapDistance);
		if( _binaryImageArray[j] > 0 )
		{
			_inputLoGImageArray[j] = ((double)_inputImageArray[j]/((double)(_maxValueInputPixelType)))/mult;
		}
		else
		{
			_inputLoGImageArray[j] = 1;
		}
	}
	
	

	std::cout << std::endl << "	Allocating:	" << _totNumPixels << " " << 2*_totNumPixels*sizeof(loGResponsePixelType) / (1024.0 * 1024.0) << " MB of memory.";
	
	_laplacian = LoGFilterType::New();
	
// 			typedef  itk::ImageFileWriter< binaryImageType > WriterType_bin;
// 			typename WriterType_bin::Pointer writer_bin = WriterType_bin::New();
// 			writer_bin->SetFileName("/media/sf_11_SharedWithUbuntu/Results/imageBinarizedLoG.tif");
// 			writer_bin->SetInput(_binaryImage);
// 			try{
// 				writer_bin->Update();
// 			}
// 			catch( itk::ExceptionObject & err ){
// 				std::cerr << "ExceptionObject caught!" << std::endl;
// 				std::cerr << err << std::endl;
// 			}
	
	
	std::cout << std::endl << "Calculating Scale: " << _seedDetectMinScale;
// 	std::cout << std::endl << "Calculating Scale: " << _seedDetectMinScale;
	loGResponse( _seedDetectMinScale, _loGResponseImageMax ); // _loGResponseImagel
	
				// Write absolute value of the distance map transform
				typedef  itk::ImageFileWriter< loGResponseImageType > WriterType_log;
				typename WriterType_log::Pointer writer_log = WriterType_log::New();
				std::stringstream outas45;
				outas45 << _seedDetectMinScale;
				std::string sas45 = outas45.str();
				std::string filenameCones_45as = "/media/sf_11_SharedWithUbuntu/Results/imageLoG_"+sas45+".mhd";
				writer_log->SetFileName(filenameCones_45as.c_str());
				writer_log->SetInput(_loGResponseImageMax);
				try{
					writer_log->Update();
				}
				catch( itk::ExceptionObject & err ){
					std::cerr << "ExceptionObject caught!" << std::endl;
					std::cerr << err << std::endl;
				}
	
	
	_distMapImageArray = _distMapImage->GetBufferPointer();
	for( int actualScale = _seedDetectMinScale+1; actualScale<=_seedDetectMaxScale; ++actualScale )
	{
		std::cout << std::endl << "Calculating Scale: " << actualScale;
// 		std::cout << std::endl << "Calculating Scale: " << actualScale;
		loGResponse( actualScale, _loGResponseImage ); // _loGResponseImage
		
				std::stringstream outas46;
				outas46 << actualScale;
				std::string sas46 = outas46.str();
				std::string filenameCones_46as = "/media/sf_11_SharedWithUbuntu/Results/imageLoG_"+sas46+".mhd";
				writer_log->SetFileName(filenameCones_46as.c_str());
				writer_log->SetInput(_loGResponseImage);
				try{
					writer_log->Update();
				}
				catch( itk::ExceptionObject & err ){
					std::cerr << "ExceptionObject caught!" << std::endl;
					std::cerr << err << std::endl;
				}
		
		//#pragma omp parallel for
		_loGResponseImageMaxArray = _loGResponseImageMax->GetBufferPointer();
		_loGResponseImageArray = _loGResponseImage->GetBufferPointer();
		#pragma omp parallel for
		for( int j = 0; j<_totNumPixels; ++j )
		{
// 			std::cout << std::endl << j;
// 			std::cout << std::endl << j;
// 			if( actualScale <= _distMapImageArray[j] )
// 			{
// 				if( _loGResponseImageMaxArray[j]>=_loGResponseImageArray[j] )
// 				{
// 					_loGResponseImageMaxArray[j] = _loGResponseImageArray[j];
// 				}
				_loGResponseImageMaxArray[j] = (_loGResponseImageMaxArray[j]>=_loGResponseImageArray[j]) ? _loGResponseImageMaxArray[j] : _loGResponseImageArray[j];
// 			}
			

		}
	}
	
// 				std::stringstream outas47;
// 				outas47 << _seedDetectMaxScale+1;
// 				std::string sas47 = outas47.str();
				std::string filenameCones_47as = "/media/sf_11_SharedWithUbuntu/Results/imageLoG_maxLoG.mhd";
				writer_log->SetFileName(filenameCones_47as.c_str());
				writer_log->SetInput(_loGResponseImageMax);
				try{
					writer_log->Update();
				}
				catch( itk::ExceptionObject & err ){
					std::cerr << "ExceptionObject caught!" << std::endl;
					std::cerr << err << std::endl;
				}
				
	
	
	double scale_xy = _spacing_XY;
	double scale_z = _spacing_Z;
	
	_binaryImageArray = _binaryImage->GetBufferPointer();
	_loGResponseImageMaxArray = _loGResponseImageMax->GetBufferPointer();
	
	seedDetecPixelType seedID = 1;
	
//#pragma omp parallel for 
	for( int k=0; k<_numStacks; ++k )
	{
// 		std::cout << std::endl << "Max Stack: " << k << " " << _spacing_Z << " " << _spacing_XY;
		long long curr_node_z = k*_numColumns*_numRows;
		for(int i=0; i<_numRows; ++i)
		{
			long long curr_node_yz = (i*_numColumns)+curr_node_z;
			for(int j=0; j<_numColumns; ++j)
			{	
				long long curr_node_xyz = j + curr_node_yz;
				
				
				
				
				
					loGResponsePixelType maxLoGResponse = _loGResponseImageMaxArray[curr_node_xyz];
					
					int min_r = (int) std::max((double)0.0,(double)i-scale_xy);
					int min_c = (int) std::max((double)0.0,(double)j-scale_xy);
					int min_s = (int) std::max((double)0.0,(double)k-scale_z);
					int max_r = (int) std::min((double)_numRows-1,(double)i+scale_xy); // In case _numRows was unsigned int
					int max_c = (int) std::min((double)_numColumns-1,(double)j+scale_xy);                         
					int max_s = (int) std::min((double)_numStacks-1,(double)k+scale_z);     
						
					for( int k_in=min_s; k_in<=max_s; ++k_in )
					{
						long long curr_node_z_in = k_in*_numColumns*_numRows;
						for(int i_in=min_r; i_in<=max_r; ++i_in)
						{
							long long  curr_node_yz_in = (i_in*_numColumns)+curr_node_z_in;
							for(int j_in=min_c; j_in<=max_c; ++j_in)
							{	
								long long  curr_node_xyz_in = j_in + curr_node_yz_in;
								
								if( _loGResponseImageMaxArray[curr_node_xyz_in] > maxLoGResponse)
								{
									maxLoGResponse = _loGResponseImageMaxArray[curr_node_xyz_in];
								}
							}
						}
					}
					
					if( _loGResponseImageMaxArray[curr_node_xyz] == maxLoGResponse)
					{
						if( _binaryImageArray[curr_node_xyz] == 0 ){ // THIS WILL NOT ALLOW SEEDS IN THE BINARY IMAGE (in case the seed points are in the background)
							_seedDetectImageArray[curr_node_xyz] = _maxValueSeedDetectPixelType;
						}
						else
						{
							_seedDetectImageArray[curr_node_xyz] = seedID;
							seedID++;
						}
					}
					else
					{
						_seedDetectImageArray[curr_node_xyz] = 0;
					}
			}
		}
	}
	
	
				// Write absolute value of the distance map transform
				typedef  itk::ImageFileWriter< seedDetectImageType > WriterType_log3;
				typename WriterType_log3::Pointer writer_log3 = WriterType_log3::New();
// 				std::stringstream outas49;
// 				outas49 << _seedDetectMaxScale+3;
// 				std::string sas49 = outas49.str();
				std::string filenameCones_49as = "/media/sf_11_SharedWithUbuntu/Results/imageLoG_seed.mhd";
				writer_log3->SetFileName(filenameCones_49as.c_str());
				writer_log3->SetInput(_seedDetectImage);
				try{
					writer_log3->Update();
				}
				catch( itk::ExceptionObject & err ){
					std::cerr << "ExceptionObject caught!" << std::endl;
					std::cerr << err << std::endl;
				}
}
	
	
	
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType, typename seedDetecPixelType, typename loGResponsePixelType >
typename itk::Image< seedDetecPixelType, 3 >::Pointer ftk::nucSecNic::SeedDetectionLoG< inputPixelType, binaryPixelType, seedDetecPixelType, loGResponsePixelType >::getSeedDetectImage()
{
	return _seedDetectImage;
}
	
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType, typename seedDetecPixelType, typename loGResponsePixelType >
typename itk::Image< loGResponsePixelType, 3 >::Pointer ftk::nucSecNic::SeedDetectionLoG< inputPixelType, binaryPixelType, seedDetecPixelType, loGResponsePixelType >::getLoGResponseImage()
{
	return _loGResponseImageMax;
}
	



	


// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType, typename seedDetecPixelType, typename loGResponsePixelType >
void ftk::nucSecNic::SeedDetectionLoG< inputPixelType, binaryPixelType, seedDetecPixelType, loGResponsePixelType >::distMap_2D()
{
	// Distance Map Transform
	typedef itk::SignedMaurerDistanceMapImageFilter<  binaryImageType, distMapImageType >  DTFilter_2D;
// 	#include <itkDanielssonDistanceMapImageFilter.h>
	
	typename DTFilter_2D::Pointer objDTFilter_2D= DTFilter_2D::New();
	objDTFilter_2D->SetInput(_binaryImage);
	objDTFilter_2D->SetSquaredDistance( false );      
	objDTFilter_2D->SetInsideIsPositive( true );
	
	try{
		objDTFilter_2D->Update() ;
	}
	catch( itk::ExceptionObject & err ){
		std::cerr << "Error calculating distance transform: " << err << endl ;
	}
	
	// In case of testing the output without absoulte value
	_distMapImage = objDTFilter_2D->GetOutput();
	
	_distMapImageArray = _distMapImage->GetBufferPointer();
	
	_maxMapDistance = 0;
	for( int i=0; i<_totNumPixels; ++i )
	{
		double ds = _distMapImageArray[i];//*100;
		if( ds <= 0 )
		{
			_distMapImageArray[i] = 0;
		}
		else
		{
			_distMapImageArray[i] = (distMapPixelType) ds;
		}
		if( _maxMapDistance < _distMapImageArray[i] )
		{
			_maxMapDistance = _distMapImageArray[i];
		}
		
		
	}
	
			// Write absolute value of the distance map transform
			typedef  itk::ImageFileWriter< distMapImageType > WriterType;
			typename WriterType::Pointer writer = WriterType::New();
			writer->SetFileName("/media/sf_11_SharedWithUbuntu/Results/imageDistanceMap.mhd");
			writer->SetInput(_distMapImage);
			try{
				writer->Update();
			}
			catch( itk::ExceptionObject & err ){
				std::cerr << "ExceptionObject caught!" << std::endl;
				std::cerr << err << std::endl;
			}
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType, typename seedDetecPixelType, typename loGResponsePixelType >
void ftk::nucSecNic::SeedDetectionLoG< inputPixelType, binaryPixelType, seedDetecPixelType, loGResponsePixelType >::loGResponse( long long actualScale, typename loGResponseImageType::Pointer &loGResponse  )
{
	_laplacian->SetNormalizeAcrossScale( true );
	_laplacian->SetInput(_inputLoGImage);
	_laplacian->SetSigma( actualScale );
	
	try
	{
		_laplacian->Update();
	}
	catch( itk::ExceptionObject & err ) 
	{ 
		std::cout << "ExceptionObject caught !" << std::endl; 
		std::cout << err << std::endl; 
	} 
	
	
	typedef itk::ImageDuplicator< loGResponseImageType > DuplicatorType;
	typename DuplicatorType::Pointer duplicatorFilter = DuplicatorType::New();

	duplicatorFilter->SetInputImage( _laplacian->GetOutput() );
	
	try{
		duplicatorFilter->Update();}
	catch( itk::ExceptionObject & err ){
		std::cerr << "ExceptionObject caught!" << std::endl;
		std::cerr << err << std::endl;
	}
	
	loGResponse = duplicatorFilter->GetOutput();
	
}


// // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// template < typename inputPixelType, typename binaryPixelType, typename seedDetecPixelType, typename loGResponsePixelType >
// typename itk::Image< binaryPixelType, 3 >::Pointer ftk::nucSecNic::SeedDetectionLoG< inputPixelType, binaryPixelType, seedDetecPixelType, loGResponsePixelType >::getBinarizedImage()
// {
// 	return _binaryImage;
// }

