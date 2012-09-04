
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType, typename seedDetecPixelType, typename loGResponsePixelType >
void ftk::nucSecNic::MaxClustering< inputPixelType, binaryPixelType, seedDetecPixelType, loGResponsePixelType >::setParameters( bool getResultImg_MaxClustering )
{
// 	_seedDetectMinScale = seedDetectMinScale;
// 	_seedDetectMaxScale = seedDetectMaxScale;
	_getResultImg_MaxClustering = getResultImg_MaxClustering;
	
	_spacing_XY = 3;
	_spacing_Z = 4;
	
	_maxNumberIterations = 10;
}


// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType, typename seedDetecPixelType, typename loGResponsePixelType >
void ftk::nucSecNic::MaxClustering< inputPixelType, binaryPixelType, seedDetecPixelType, loGResponsePixelType >::setInput( const ftk::Image::Info* info, typename inputImageType::Pointer inputImage, typename binaryImageType::Pointer binaryImage, typename seedDetectImageType::Pointer seedDetectImage, typename loGResponseImageType::Pointer loGResponseImage )
{
	_inputImage = inputImage;
	_inputImageArray = _inputImage->GetBufferPointer();

	_binaryImage = binaryImage;
	_binaryImageArray = _binaryImage->GetBufferPointer();
	
	_seedDetectImage = seedDetectImage;
	_seedDetectImageArray = _seedDetectImage->GetBufferPointer();
	
	_loGResponseImage = loGResponseImage;
	_loGResponseImageArray = _loGResponseImage->GetBufferPointer();
	
	
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
void ftk::nucSecNic::MaxClustering< inputPixelType, binaryPixelType, seedDetecPixelType, loGResponsePixelType >::runMaxClustering()
{
	
	
				// Write absolute value of the distance map transform
				typedef  itk::ImageFileWriter< loGResponseImageType > WriterType_log;
				typename WriterType_log::Pointer writer_log = WriterType_log::New();
// 				std::stringstream outas45;
// 				outas45 << _seedDetectMinScale;
// 				std::string sas45 = outas45.str();
				std::string filenameCones_45as = "/media/sf_11_SharedWithUbuntu/Results/imageLoG_maxClus.mhd";
				writer_log->SetFileName(filenameCones_45as.c_str());
				writer_log->SetInput(_loGResponseImage);
				try{
					writer_log->Update();
				}
				catch( itk::ExceptionObject & err ){
					std::cerr << "ExceptionObject caught!" << std::endl;
					std::cerr << err << std::endl;
				}
	
	
	
	_seedDetectImageArray = _seedDetectImage->GetBufferPointer();
	
	_loGResponseImageArray = _loGResponseImage->GetBufferPointer();
	
	double scale_xy = _spacing_XY;
	double scale_z = _spacing_Z;
	
	int*** maxClustImage;
    
	
	// I CAN JUST CREATE AN ITK IMAGE, AND ACTUALLY
	//create maxClustImage and initialize it with its index (node) value
	maxClustImage = (int ***) malloc(_numRows*sizeof(int**)); 
    
	#pragma omp parallel for
	for(int i=0; i<_numRows; i++)
	{        
		maxClustImage[i] = (int **) malloc(_numColumns*sizeof(int*));
		for(int j=0; j<_numColumns; j++)
		{			
			maxClustImage[i][j] = (int *) malloc(_numStacks*sizeof(int));
			for(int k=0; k<_numStacks; k++)
			{				
				maxClustImage[i][j][k] = (k*_numRows*_numColumns)+(i*_numColumns)+j;//LMX;
			}
		}
	}
    
	cerr << "maxClustImage initialized" << endl;
	
	unsigned short *maxClustImageRows = new unsigned short[_numRows*_numColumns*_numStacks];
	unsigned short *maxClustImageColu = new unsigned short[_numRows*_numColumns*_numStacks];
	unsigned short *maxClustImageStac = new unsigned short[_numRows*_numColumns*_numStacks];

	

	#pragma omp parallel for
	for(int i=0; i<_numRows; i++)
	{
		for(int j=0; j<_numColumns; j++)
		{
			for(int k=0; k<_numStacks; k++)
			{        				
				int min_r = (int) std::max((double)(0.0),(double)(i-scale_xy));
				int min_c = (int) std::max((double)(0.0),(double)(j-scale_xy));
				int min_z = (int) std::max((double)(0.0),(double)(k-scale_z));
				int max_r = (int) std::min((double)(_numRows-1),(double)(i+scale_xy));
				int max_c = (int) std::min((double)(_numColumns-1),(double)(j+scale_xy));                         
				int max_z = (int) std::min((double)(_numStacks-1),(double)(k+scale_z));

				if(_seedDetectImageArray[(k*_numRows*_numColumns)+(i*_numColumns)+j] !=0 )
					continue;					
				else
				{                                              
					loGResponsePixelType maxLoGResponse = _loGResponseImageArray[min_z*(_numRows*_numColumns)+(min_r*_numColumns)+min_c];
					int maxCol = min_c;
					int maxRow = min_r;
					int maxSta = min_z;
					
					for( int k_in=min_z; k_in<=max_z; ++k_in )
					{
						long long curr_node_z_in = k_in*_numColumns*_numRows;
						for(int i_in=min_r; i_in<=max_r; ++i_in)
						{
							long long  curr_node_yz_in = (i_in*_numColumns)+curr_node_z_in;
							for(int j_in=min_c; j_in<=max_c; ++j_in)
							{	
								long long  curr_node_xyz_in = j_in + curr_node_yz_in;
								
								if( _loGResponseImageArray[curr_node_xyz_in] >= maxLoGResponse)
								{
									maxLoGResponse = _loGResponseImageArray[curr_node_xyz_in];
									maxCol = j_in;
									maxRow = i_in;
									maxSta = k_in;
								}
							}
						}
					}
					
					maxClustImageRows[i * (_numColumns * _numStacks) + j * _numStacks + k] = maxRow;
					maxClustImageColu[i * (_numColumns * _numStacks) + j * _numStacks + k] = maxCol;
					maxClustImageStac[i * (_numColumns * _numStacks) + j * _numStacks + k] = maxSta;
				}				
			}
		}
	}
	
	std::cout << "Max_response array done" << endl;

	for(int i=0; i<_numRows; i++)
	{
		for(int j=0; j<_numColumns; j++)
		{
			for(int k=0; k<_numStacks; k++)
			{
				if(_seedDetectImageArray[(k*_numRows*_numColumns)+(i*_numColumns)+j] !=0 )
					continue;
				else
					maxClustImage[i][j][k] =	maxClustImage[	maxClustImageRows[i * (_numColumns * _numStacks) + j * _numStacks + k]	]
														[	maxClustImageColu[i * (_numColumns * _numStacks) + j * _numStacks + k]	]
														[	maxClustImageStac[i * (_numColumns * _numStacks) + j * _numStacks + k]	];
			}
		}
	}
	
	delete [] maxClustImageRows;
	delete [] maxClustImageColu;
	delete [] maxClustImageStac;

	
				typedef itk::Image< int, 3 > maxClusterImageType;
				typename maxClusterImageType::Pointer maxClusterImage = maxClusterImageType::New();
				typename maxClusterImageType::IndexType start_2;
				start_2[0] = 0;
				start_2[1] = 0;
				start_2[2] = 0;
				typename maxClusterImageType::SizeType size_2;
				size_2[0] = _numColumns;
				size_2[1] = _numRows;
				size_2[2] = _numStacks;
				typename maxClusterImageType::RegionType region_2;
				region_2.SetIndex( start_2 );
				region_2.SetSize( size_2 );
				
				double spacing_2[3];
				spacing_2[0] = 1; //spacing along x
				spacing_2[1] = 1; //spacing along y
				spacing_2[2] = _spacing_Z; //spacing along z

				maxClusterImage->SetRegions( region_2 );
				maxClusterImage->SetSpacing( spacing_2 );
				
				try{
					maxClusterImage->Allocate();
				}
				catch( itk::ExceptionObject & err ){
					std::cerr << "ExceptionObject caught!" << std::endl;
					std::cerr << err << std::endl;
				}
				const typename  maxClusterImageType::PixelType ceros_2 = 0;
				maxClusterImage->FillBuffer( ceros_2 );
				maxClusterImage->Update();
	
	int change = 1;
	double LM;
	cerr << "Entering main Clustering Loop" << endl;
	//Now continue to update until no more changes occur, eventually will have clusters pointing to seeds	
	int iterr = 0;
	while(change)	
	{
		
				typename maxClusterImageType::PixelType *maxClusterImageArray = maxClusterImage->GetBufferPointer();
				for(int i=0; i<_numRows; ++i)
				{
					for(int j=0; j<_numColumns; ++j)
					{
						for(int k=0; k<_numStacks; ++k)
						{
							long long curr_node_xyz = k*(_numRows*_numColumns) + i*(_numColumns) + j;
							maxClusterImageArray[curr_node_xyz] = maxClustImage[i][j][k];
						}
					}
				}
				// Write absolute value of the distance map transform
				typedef  itk::ImageFileWriter< maxClusterImageType > WriterType_log3;
				typename WriterType_log3::Pointer writer_log4 = WriterType_log3::New();
				std::stringstream outas499;
				outas499 << iterr+1;
				std::string sas499 = outas499.str();
				std::string filenameCones_499as = "/media/sf_11_SharedWithUbuntu/Results/imageMaxCluster_"+sas499+".mhd";
				writer_log4->SetFileName(filenameCones_499as.c_str());
				writer_log4->SetInput(maxClusterImage);
				try{
					writer_log4->Update();
				}
				catch( itk::ExceptionObject & err ){
					std::cerr << "ExceptionObject caught!" << std::endl;
					std::cerr << err << std::endl;
				}
	    
		//For now, limit it to a maximum of 10 iterations
		iterr++;
		if(iterr == _maxNumberIterations)
			break;		
		change=0;
		
		
		for(int i=0; i<_numRows; i++)
		{
			for(int j=0; j<_numColumns; j++)
			{
				for(int k=0; k<_numStacks; k++)
				{                   
					LM = maxClustImage[i][j][k];
					if(LM==0)
						continue;													    
						
					//Calculate coordinates of local maximum based on its index
					int rem = ((long)LM) % (_numRows*_numColumns);
					int Z = ((int)LM-rem) / (_numRows*_numColumns); 
					int C = ((long)rem) % _numColumns;
					int R = (rem-C)/_numColumns;
						
			
					//Check for seed (already a local max value) or connected to seed (my local maximum is a seed point)					
					if(_seedDetectImageArray[(k*_numRows*_numColumns)+(i*_numColumns)+j] !=0 || _seedDetectImageArray[(Z*_numRows*_numColumns)+(R*_numColumns)+C]!=0 ) 
						continue;
					else
					{
						change++;
						maxClustImage[i][j][k]=maxClustImage[R][C][Z];
					}
				}
			}
		}
		cerr<<"change="<<change<<endl;
	}
	
	
	// STORE THE RESULT IN THE ITK IMAGE
	
	_maxClustImage = maxClustImageType::New();
	typename maxClustImageType::IndexType start;
	start[0] = 0;
	start[1] = 0;
	start[2] = 0;
	typename maxClustImageType::SizeType size;
	size[0] = _numColumns;
	size[1] = _numRows;
	size[2] = _numStacks;
	typename maxClustImageType::RegionType region;
	region.SetIndex( start );
	region.SetSize( size );
	
	double spacing[3];
	spacing[0] = 1; //spacing along x
	spacing[1] = 1; //spacing along y
	spacing[2] = _spacing_Z; //spacing along z

	_maxClustImage->SetRegions( region );
	_maxClustImage->SetSpacing( spacing );
	
	try{
		_maxClustImage->Allocate();
	}
	catch( itk::ExceptionObject & err ){
		std::cerr << "ExceptionObject caught!" << std::endl;
		std::cerr << err << std::endl;
	}
	const typename  maxClustImageType::PixelType ceros = 0;
	_maxClustImage->FillBuffer( ceros_2 );
	_maxClustImage->Update();
	
	_maxClustImageArray = _maxClustImage->GetBufferPointer();
	
	
	
	
	for(int i=0; i<_numRows; ++i)
	{        
		for(int j=0; j<_numColumns; ++j)
		{
			for(int k=0; k<_numStacks; ++k)
			{
				LM = maxClustImage[i][j][k];
				if(_seedDetectImageArray[(int)LM] == _maxValueSeedDetectPixelType/*65535*/ || _binaryImageArray[(k*_numRows*_numColumns)+(i*_numColumns)+j]==0)
					_maxClustImageArray[(k*_numRows*_numColumns)+(i*_numColumns)+j] = 0;
				else
				{
					//modified by Yousef on 8/21/2009
					//if the distance between me and my seed is more than a threshold.. then ignore me
					/*int rem = ((long)LM) % (r*c);
					int Z = (LM-rem) / (r*c); 
					int C = ((long)rem) % c;
					int R = (rem-C)/c;
					double d = (i-R)*(i-R) + (j-C)*(j-C) + 3*(k-Z)*(k-Z);
					d = sqrt(d);
					if(d>10)
						out1[(k*r*c)+(i*c)+j] = 0;
					else*/
					_maxClustImageArray[(k*_numRows*_numColumns)+(i*_numColumns)+j] = _seedDetectImageArray[(int)LM];
				}
			}
		}
	}    

	
	for(int i=0; i<_numRows; i++)
	{        
		for(int j=0; j<_numColumns; j++)
		{			
			free(maxClustImage[i][j]);
		}
		free(maxClustImage[i]);
	}
	free(maxClustImage);
}



// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType, typename seedDetecPixelType, typename loGResponsePixelType >
typename itk::Image< seedDetecPixelType, 3 >::Pointer ftk::nucSecNic::MaxClustering< inputPixelType, binaryPixelType, seedDetecPixelType, loGResponsePixelType >::getMaxClustImage()
{
	return _maxClustImage;
}

