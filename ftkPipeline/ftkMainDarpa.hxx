
// #include "ftkMainDarpa.h"

template<typename TGFP, typename TSOMA >
void ftkMainDarpa::projectImageRGB( std::string inputImageNameGFP /*GFP*/, std::string inputImageNameSOMA /*Soma*/, std::string inputImageName3 /*Other */, std::string outputPath, std::string NAME )
{
	typedef itk::ImageRegionIterator< TGFP > IteratorType_GFP;
	typedef itk::ImageRegionIterator< TSOMA > IteratorType_SOMA;

	typename TGFP::Pointer inputImageGFP = readImage< TGFP >(inputImageNameGFP.c_str());
	typename TSOMA::Pointer inputImageSOMA = readImage< TSOMA >(inputImageNameSOMA.c_str());
// 	typename TGFP::Pointer inputImage3 = readImage< TGFP >(inputImageName3.c_str()); // Nothing

	int found=inputImageNameGFP.find(".");
	std::string inputImageNameLocal = inputImageNameGFP.substr(0,found);
	
	found = inputImageNameLocal.find_last_of("/\\");
	inputImageNameLocal = inputImageNameLocal.substr(found+1);
	
	std::vector< typename TGFP::Pointer > inputImageGFPProjected(3);
	inputImageGFPProjected = getProjectImage< TGFP, TGFP >( inputImageGFP, "ORG" );
	
	std::vector< typename TSOMA::Pointer > inputImageSOMAProjected(3);
	if( NAME.find("GF_PDAP_I") )
		inputImageSOMAProjected = getProjectImage< TSOMA, TSOMA >( inputImageSOMA, "ORG" );
	else
		inputImageSOMAProjected = getProjectImage< TSOMA, TSOMA >( inputImageSOMA, "BIN" );
	
	typedef typename itk::RGBPixel< typename TGFP::PixelType > RGBPixelType;
	typedef itk::Image<RGBPixelType,3> RGBImageType;
	typedef itk::ImageRegionIterator<RGBImageType> RGBIteratorType;
	
// 	Z Image
	for( int channel=0; channel<3; ++channel )
	{
		typename RGBImageType::Pointer imageProjectZ = RGBImageType::New();
		typename RGBImageType::IndexType imdexZ;
		imdexZ.Fill(0);
		typename RGBImageType::SizeType sizeZ;
		sizeZ = inputImageGFPProjected[channel]->GetLargestPossibleRegion().GetSize();
		typename RGBImageType::RegionType regionZ;
		regionZ.SetIndex(imdexZ);
		regionZ.SetSize(sizeZ);
		imageProjectZ->SetRegions(regionZ);
		try
		{
			imageProjectZ->Allocate();
		}
		catch(itk::ExceptionObject &err)
		{
			std::cerr << "ExceptionObject caught!" <<std::endl;
			std::cerr << err << std::endl;
		}
		IteratorType_GFP itinputImageGFP(inputImageGFPProjected[channel], inputImageGFPProjected[channel]->GetRequestedRegion()); itinputImageGFP.GoToBegin();
		IteratorType_SOMA itinputImageSOMA(inputImageSOMAProjected[channel], inputImageSOMAProjected[channel]->GetRequestedRegion()); itinputImageSOMA.GoToBegin();
		RGBIteratorType itimageProjectZ(imageProjectZ, imageProjectZ->GetRequestedRegion()); itimageProjectZ.GoToBegin();
		
		for (itimageProjectZ.GoToBegin(); !itimageProjectZ.IsAtEnd(); ++itimageProjectZ, ++itinputImageGFP, ++itinputImageSOMA) 
		{
			typename RGBImageType::PixelType newPixel;
			
			int foundPro=NAME.find("GFP");
			if (foundPro!=std::string::npos)
			{
				newPixel.SetRed( itinputImageGFP.Get() );
				newPixel.SetGreen( itinputImageGFP.Get() );
				
				if( itinputImageSOMA.Get() != 0 )
					newPixel.SetBlue( std::numeric_limits<typename TGFP::PixelType>::max() );
				else
					newPixel.SetBlue( 0 );
			}
			foundPro=NAME.find("DAPI");
			if (foundPro!=std::string::npos)
			{
// 				newPixel.SetBlue( itinputImageGFP.Get()*4 );
// 				newPixel.SetRed( 0 );
// 				
// 				if( itinputImageSOMA.Get() != 0 )
// 					newPixel.SetGreen( std::numeric_limits<typename TGFP::PixelType>::max()/4 );
// 				else
// 					newPixel.SetGreen( 0 );

				newPixel.SetRed( itinputImageGFP.Get() );
				newPixel.SetGreen( itinputImageGFP.Get() );
				
				if( itinputImageSOMA.Get() != 0 )
					newPixel.SetBlue( std::numeric_limits<typename TGFP::PixelType>::max() );
				else
					newPixel.SetBlue( 0 );
			}
			foundPro=NAME.find("GF_PDAP_I");
			if (foundPro!=std::string::npos)
			{
				newPixel.SetRed( itinputImageGFP.Get() );
				newPixel.SetGreen( itinputImageGFP.Get() );
				newPixel.SetBlue( itinputImageSOMA.Get() );
			}
			itimageProjectZ.Set(newPixel);
		}
		std::string temp5b;
		if( channel== 0 )
			temp5b = outputPath + "/zRGB_"+NAME+"Pro_Z.tif";
		if( channel== 1 )
			temp5b = outputPath + "/zRGB_"+NAME+"Pro_Y.tif";
		if( channel== 2 )
			temp5b = outputPath + "/zRGB_"+NAME+"Pro_X.tif";
		writeImage< RGBImageType >(imageProjectZ,temp5b.c_str());
	}
	
}




template<typename TINPUT, typename TOUTPUT >
void ftkMainDarpa::rescaleImage( std::string inputImageName, std::string outputImageName )
{
	typename TINPUT::Pointer inputImage = readImage< TINPUT >(inputImageName.c_str());
	std::string temp1a = outputImageName;
	writeImageRescaled< TINPUT, TOUTPUT >(inputImage,temp1a.c_str());
	
}

template<typename TINPUT >
void ftkMainDarpa::saveNRRD( std::string inputImageName )
{
	typename TINPUT::Pointer inputImage = readImage< TINPUT >(inputImageName.c_str());
	
	int found=inputImageName.find(".");
	std::string inputImageNameLocal = inputImageName.substr(0,found);
	
	std::string temp3 = inputImageNameLocal + ".nrrd";
	writeImage< TINPUT >(inputImage,temp3.c_str());
}


template<typename TINPUT, typename TOUTPUT >
std::vector< typename TOUTPUT::Pointer > ftkMainDarpa::getProjectImage( typename TINPUT::Pointer inputImage, std::string projectOptions )
{
// 	typename TINPUT::Pointer inputImage = readImage< TINPUT >(inputImageName.c_str());

// 	int found=inputImageName.find(".");
// 	std::string inputImageNameLocal = inputImageName.substr(0,found);
// 	std::cout<<std::endl<<inputImageNameLocal;
	
// 	found = inputImageNameLocal.find_last_of("/\\");
// 	inputImageNameLocal = inputImageNameLocal.substr(found+1);
// 	std::cout<<std::endl<<inputImageNameLocal;
	
	std::vector< typename TOUTPUT::Pointer > vectorOut(3);
	
	typename TINPUT::PixelType * inputImageArray = inputImage->GetBufferPointer();
	itk::Size<3> inputImage_sizez = inputImage->GetLargestPossibleRegion().GetSize();
	unsigned long long inputImage_slice_size = inputImage_sizez[1] * inputImage_sizez[0];
	unsigned long long inputImage_row_size = inputImage_sizez[0];	

	// Create z projecImage
	typename TINPUT::Pointer zProjectImage = TINPUT::New();
	typename TINPUT::PointType originz;
	originz[0] = 0; 
	originz[1] = 0;
	originz[2] = 0;
	zProjectImage->SetOrigin( originz );
	typename TINPUT::IndexType startz;
	startz[0] = 0;
	startz[1] = 0;
	startz[2] = 0;
	typename TINPUT::SizeType sizez;
	sizez[0] = inputImage_sizez[0];
	sizez[1] = inputImage_sizez[1];
	sizez[2] = 1;
	typename TINPUT::RegionType regionz;
	regionz.SetSize ( sizez  );
	regionz.SetIndex( startz );
	zProjectImage->SetRegions( regionz );
	zProjectImage->Allocate();
	zProjectImage->FillBuffer(0);
	zProjectImage->Update();
	typename TINPUT::PixelType * zProjectImageArray = zProjectImage->GetBufferPointer();
	
	// Create Y projecImage
	typename TINPUT::Pointer yProjectImage = TINPUT::New();
	typename TINPUT::PointType originy;
	originy[0] = 0; 
	originy[1] = 0;
	originy[2] = 0;
	yProjectImage->SetOrigin( originy );
	typename TINPUT::IndexType starty;
	starty[0] = 0;
	starty[1] = 0;
	starty[2] = 0;
	typename TINPUT::SizeType sizey;
	sizey[0] = inputImage_sizez[0];
	sizey[1] = inputImage_sizez[2];
	sizey[2] = 1;
	typename TINPUT::RegionType regiony;
	regiony.SetSize ( sizey  );
	regiony.SetIndex( starty );
	yProjectImage->SetRegions( regiony );
	yProjectImage->Allocate();
	yProjectImage->FillBuffer(0);
	yProjectImage->Update();
	typename TINPUT::PixelType * yProjectImageArray = yProjectImage->GetBufferPointer();
	
	// Create X projecImage
	typename TINPUT::Pointer xProjectImage = TINPUT::New();
	typename TINPUT::PointType originx;
	originx[0] = 0; 
	originx[1] = 0;
	originx[2] = 0;
	xProjectImage->SetOrigin( originx );
	typename TINPUT::IndexType startx;
	startx[0] = 0;
	startx[1] = 0;
	startx[2] = 0;
	typename TINPUT::SizeType sizex;
	sizex[0] = inputImage_sizez[2];
	sizex[1] = inputImage_sizez[1];
	sizex[2] = 1;
	typename TINPUT::RegionType regionx;
	regionx.SetSize ( sizex  );
	regionx.SetIndex( startx );
	xProjectImage->SetRegions( regionx );
	xProjectImage->Allocate();
	xProjectImage->FillBuffer(0);
	xProjectImage->Update();
	typename TINPUT::PixelType * xProjectImageArray = xProjectImage->GetBufferPointer();
	
// 	std::cout << std::endl << "NOW IS GOING TO PROJECT"<<std::flush;
	
		
	#if _OPENMP < 200805L 
		#pragma omp parallel for //collapse(2) //TEST
	#else
		#pragma omp parallel for collapse(2)
	#endif

	for( long long x=0; x<inputImage_sizez[0]; ++x)
	{
		for( long long y=0; y<inputImage_sizez[1]; ++y)
		{
			typename TINPUT::PixelType max_val = 0;
			for( long long z=0; z<inputImage_sizez[2]; ++z)
			{
				typename TINPUT::PixelType value = inputImageArray[(inputImage_slice_size*z) + (inputImage_row_size*y) + (x)];
				if( max_val < value )
					max_val = value;
			}
			zProjectImageArray[(inputImage_row_size*y) + (x)] = max_val;
		}
	}
	
	#if _OPENMP < 200805L 
		#pragma omp parallel for //collapse(2) //TEST
	#else
		#pragma omp parallel for collapse(2)
	#endif
	for( long long z=0; z<inputImage_sizez[2]; ++z)
	{
		for( long long x=0; x<inputImage_sizez[0]; ++x)
		{
			typename TINPUT::PixelType max_val = 0;
			for( long long y=0; y<inputImage_sizez[1]; ++y)
			{
				typename TINPUT::PixelType value = inputImageArray[(inputImage_slice_size*z) + (inputImage_row_size*y) + (x)];
				if( max_val < value )
					max_val = value;
			}
			yProjectImageArray[(inputImage_sizez[0]*z) + (x)] = max_val;
		}
	}
	
	#if _OPENMP < 200805L 
		#pragma omp parallel for //collapse(2) //TEST
	#else
		#pragma omp parallel for collapse(2)
	#endif
	for( long long z=0; z<inputImage_sizez[2]; ++z)
	{
		for( long long y=0; y<inputImage_sizez[1]; ++y)
		{
			typename TINPUT::PixelType max_val = 0;
			for( long long x=0; x<inputImage_sizez[0]; ++x)
			{
				typename TINPUT::PixelType value = inputImageArray[(inputImage_slice_size*z) + (inputImage_row_size*y) + (x)];
				if( max_val < value )
					max_val = value;
			}
			xProjectImageArray[(inputImage_sizez[2]*y) + (z)] = max_val;
		}
	}


	int foundPro=projectOptions.find("ORG");
	if (foundPro!=std::string::npos)
	{
		vectorOut[0] = zProjectImage;
		vectorOut[1] = yProjectImage;
		vectorOut[2] = xProjectImage;
	}
	
	int foundBin=projectOptions.find("BIN");
	if (foundBin!=std::string::npos)
	{
		#if _OPENMP < 200805L 
			#pragma omp parallel for //collapse(2) //TEST
		#else
			#pragma omp parallel for collapse(2)
		#endif
		for( long long x=0; x<inputImage_sizez[0]; ++x)
		{
			for( long long y=0; y<inputImage_sizez[1]; ++y)
			{
				if( zProjectImageArray[(inputImage_row_size*y) + (x)] != 0 )
					zProjectImageArray[(inputImage_row_size*y) + (x)] = std::numeric_limits<typename TOUTPUT::PixelType>::max();
			}
		}
		#if _OPENMP < 200805L 
			#pragma omp parallel for //collapse(2) //TEST
		#else
			#pragma omp parallel for collapse(2)
		#endif
		for( long long z=0; z<inputImage_sizez[2]; ++z)
		{
			for( long long x=0; x<inputImage_sizez[0]; ++x)
			{
				if( yProjectImageArray[(inputImage_sizez[0]*z) + (x)] != 0 )
					yProjectImageArray[(inputImage_sizez[0]*z) + (x)] = std::numeric_limits<typename TOUTPUT::PixelType>::max();
			}
		}
		
		#if _OPENMP < 200805L 
			#pragma omp parallel for //collapse(2) //TEST
		#else
			#pragma omp parallel for collapse(2)
		#endif
		for( long long z=0; z<inputImage_sizez[2]; ++z)
		{
			for( long long y=0; y<inputImage_sizez[1]; ++y)
			{
				if( xProjectImageArray[(inputImage_sizez[2]*y) + (z)] != 0 )
					xProjectImageArray[(inputImage_sizez[2]*y) + (z)] = std::numeric_limits<typename TOUTPUT::PixelType>::max();
			}
		}
		vectorOut[0] = zProjectImage;
		vectorOut[1] = yProjectImage;
		vectorOut[2] = xProjectImage;
	}
	return vectorOut;
}



template<typename TINPUT, typename TOUTPUT >
void ftkMainDarpa::projectImage( std::string inputImageName, std::string outputPath, std::string projectOptions, std::string imageType )
{
// 	std::cout << std::endl << "KK" << inputImageName<< " " << outputPath << " " << projectOptions;
	
	std::string tipoImagen;
	int foundType=imageType.find("TIFF");
	if (foundType!=std::string::npos)
	{
		tipoImagen = ".tif";
	}
	foundType=imageType.find("NRRD");
	if (foundType!=std::string::npos)
	{
		tipoImagen = ".nrrd";
	}
	
	typename TINPUT::Pointer inputImage = readImage< TINPUT >(inputImageName.c_str());

	int found=inputImageName.find(".");
	std::string inputImageNameLocal = inputImageName.substr(0,found);
	
	found = inputImageNameLocal.find_last_of("/\\");
	inputImageNameLocal = inputImageNameLocal.substr(found+1);
	
	typename TINPUT::PixelType * inputImageArray = inputImage->GetBufferPointer();
	itk::Size<3> inputImage_sizez = inputImage->GetLargestPossibleRegion().GetSize();
	unsigned long long inputImage_slice_size = inputImage_sizez[1] * inputImage_sizez[0];
	unsigned long long inputImage_row_size = inputImage_sizez[0];	

	// Create z projecImage
	typename TINPUT::Pointer zProjectImage = TINPUT::New();
	typename TINPUT::PointType originz;
	originz[0] = 0; 
	originz[1] = 0;
	originz[2] = 0;
	zProjectImage->SetOrigin( originz );
	typename TINPUT::IndexType startz;
	startz[0] = 0;
	startz[1] = 0;
	startz[2] = 0;
	typename TINPUT::SizeType sizez;
	sizez[0] = inputImage_sizez[0];
	sizez[1] = inputImage_sizez[1];
	sizez[2] = 1;
	typename TINPUT::RegionType regionz;
	regionz.SetSize ( sizez  );
	regionz.SetIndex( startz );
	zProjectImage->SetRegions( regionz );
	zProjectImage->Allocate();
	zProjectImage->FillBuffer(0);
	zProjectImage->Update();
	typename TINPUT::PixelType * zProjectImageArray = zProjectImage->GetBufferPointer();
	
	// Create Y projecImage
	typename TINPUT::Pointer yProjectImage = TINPUT::New();
	typename TINPUT::PointType originy;
	originy[0] = 0; 
	originy[1] = 0;
	originy[2] = 0;
	yProjectImage->SetOrigin( originy );
	typename TINPUT::IndexType starty;
	starty[0] = 0;
	starty[1] = 0;
	starty[2] = 0;
	typename TINPUT::SizeType sizey;
	sizey[0] = inputImage_sizez[0];
	sizey[1] = inputImage_sizez[2];
	sizey[2] = 1;
	typename TINPUT::RegionType regiony;
	regiony.SetSize ( sizey  );
	regiony.SetIndex( starty );
	yProjectImage->SetRegions( regiony );
	yProjectImage->Allocate();
	yProjectImage->FillBuffer(0);
	yProjectImage->Update();
	typename TINPUT::PixelType * yProjectImageArray = yProjectImage->GetBufferPointer();
	
	// Create X projecImage
	typename TINPUT::Pointer xProjectImage = TINPUT::New();
	typename TINPUT::PointType originx;
	originx[0] = 0; 
	originx[1] = 0;
	originx[2] = 0;
	xProjectImage->SetOrigin( originx );
	typename TINPUT::IndexType startx;
	startx[0] = 0;
	startx[1] = 0;
	startx[2] = 0;
	typename TINPUT::SizeType sizex;
	sizex[0] = inputImage_sizez[2];
	sizex[1] = inputImage_sizez[1];
	sizex[2] = 1;
	typename TINPUT::RegionType regionx;
	regionx.SetSize ( sizex  );
	regionx.SetIndex( startx );
	xProjectImage->SetRegions( regionx );
	xProjectImage->Allocate();
	xProjectImage->FillBuffer(0);
	xProjectImage->Update();
	typename TINPUT::PixelType * xProjectImageArray = xProjectImage->GetBufferPointer();
	
// 	std::cout << std::endl << "NOW IS GOING TO PROJECT"<<std::flush;
	
	//#pragma omp parallel for collapse(2)
	for( long long x=0; x<inputImage_sizez[0]; ++x)
	{
		for( long long y=0; y<inputImage_sizez[1]; ++y)
		{
			typename TINPUT::PixelType max_val = 0;
			for( long long z=0; z<inputImage_sizez[2]; ++z)
			{
				typename TINPUT::PixelType value = inputImageArray[(inputImage_slice_size*z) + (inputImage_row_size*y) + (x)];
				if( max_val < value )
					max_val = value;
			}
			zProjectImageArray[(inputImage_row_size*y) + (x)] = max_val;
		}
	}
	
	#if _OPENMP < 200805L 
		#pragma omp parallel for //collapse(2) //TEST
	#else
		#pragma omp parallel for collapse(2)
	#endif
	for(long long z=0; z<inputImage_sizez[2]; ++z)
	{
		for( long long x=0; x<inputImage_sizez[0]; ++x)
		{
			typename TINPUT::PixelType max_val = 0;
			for( long long y=0; y<inputImage_sizez[1]; ++y)
			{
				typename TINPUT::PixelType value = inputImageArray[(inputImage_slice_size*z) + (inputImage_row_size*y) + (x)];
				if( max_val < value )
					max_val = value;
			}
			yProjectImageArray[(inputImage_sizez[0]*z) + (x)] = max_val;
		}
	}
	
	#if _OPENMP < 200805L 
		#pragma omp parallel for //collapse(2) //TEST
	#else
		#pragma omp parallel for collapse(2)
	#endif
	for(long long z=0; z<inputImage_sizez[2]; ++z)
	{
		for( long long y=0; y<inputImage_sizez[1]; ++y)
		{
			typename TINPUT::PixelType max_val = 0;
			for( long long x=0; x<inputImage_sizez[0]; ++x)
			{
				typename TINPUT::PixelType value = inputImageArray[(inputImage_slice_size*z) + (inputImage_row_size*y) + (x)];
				if( max_val < value )
					max_val = value;
			}
			xProjectImageArray[(inputImage_sizez[2]*y) + (z)] = max_val;
		}
	}


	int foundPro=projectOptions.find("ORG");
	if (foundPro!=std::string::npos)
	{
		std::string temp3 = outputPath + "/" + inputImageNameLocal + "zPro_Z"+tipoImagen;
		writeImage< TINPUT >(zProjectImage,temp3.c_str());
	
		std::string temp4 = outputPath + "/" + inputImageNameLocal + "zPro_Y"+tipoImagen;
		writeImage< TINPUT >(yProjectImage,temp4.c_str());
		
		std::string temp5 = outputPath + "/" + inputImageNameLocal + "zPro_X"+tipoImagen;
		writeImage< TINPUT >(xProjectImage,temp5.c_str());
	
	}
	int foundRes=projectOptions.find("RES");
	if (foundRes!=std::string::npos)
	{
		std::string temp3a = outputPath + "/" + inputImageNameLocal + "zPro_Z_Re"+tipoImagen;
		writeImageRescaled< TINPUT, TOUTPUT >(zProjectImage,temp3a.c_str());
	
		std::string temp4a = outputPath + "/" + inputImageNameLocal + "zPro_Y_Re"+tipoImagen;
		writeImageRescaled< TINPUT, TOUTPUT >(yProjectImage,temp4a.c_str());
		
		std::string temp5a = outputPath + "/" + inputImageNameLocal + "zPro_X_Re"+tipoImagen;
		writeImageRescaled< TINPUT, TOUTPUT >(xProjectImage,temp5a.c_str());
	}
		
	int foundBin=projectOptions.find("BIN");
	if (foundBin!=std::string::npos)
	{
		#if _OPENMP < 200805L 
			#pragma omp parallel for //collapse(2) //TEST
		#else
			#pragma omp parallel for collapse(2)
		#endif
		for(long long x=0; x<inputImage_sizez[0]; ++x)
		{
			for( long long y=0; y<inputImage_sizez[1]; ++y)
			{
				if( zProjectImageArray[(inputImage_row_size*y) + (x)] != 0 )
					zProjectImageArray[(inputImage_row_size*y) + (x)] = std::numeric_limits<typename TINPUT::PixelType>::max();
			}
		}
		#if _OPENMP < 200805L 
			#pragma omp parallel for //collapse(2) //TEST
		#else
			#pragma omp parallel for collapse(2)
		#endif
		for(long long z=0; z<inputImage_sizez[2]; ++z)
		{
			for( long long x=0; x<inputImage_sizez[0]; ++x)
			{
				if( yProjectImageArray[(inputImage_sizez[0]*z) + (x)] != 0 )
					yProjectImageArray[(inputImage_sizez[0]*z) + (x)] = std::numeric_limits<typename TINPUT::PixelType>::max();
			}
		}
		
		#if _OPENMP < 200805L 
			#pragma omp parallel for //collapse(2) //TEST
		#else
			#pragma omp parallel for collapse(2)
		#endif
		for(long long z=0; z<inputImage_sizez[2]; ++z)
		{
			for(long long y=0; y<inputImage_sizez[1]; ++y)
			{
				if( xProjectImageArray[(inputImage_sizez[2]*y) + (z)] != 0 )
					xProjectImageArray[(inputImage_sizez[2]*y) + (z)] = std::numeric_limits<typename TINPUT::PixelType>::max();
			}
		}

		std::string temp3b = outputPath + "/" + inputImageNameLocal + "zPro_Z_Bin"+tipoImagen;
		writeImage< TINPUT >(zProjectImage,temp3b.c_str());
		std::string temp4b = outputPath + "/" + inputImageNameLocal + "zPro_Y_Bin"+tipoImagen;
		writeImage< TINPUT >(yProjectImage,temp4b.c_str());
		std::string temp5b = outputPath + "/" + inputImageNameLocal + "zPro_X_Bin"+tipoImagen;
		writeImage< TINPUT >(xProjectImage,temp5b.c_str());
	}
	
	int foundHisto=projectOptions.find("HISTO");
	if (foundHisto!=std::string::npos)
	{
		std::string temp3 = outputPath + "/" + inputImageNameLocal + "zHisto.txt";
		std::vector< std::vector< unsigned long long > > histoGram(inputImage_sizez[2]);
		for(long long z=0; z<inputImage_sizez[2]; ++z)
		{
			std::vector< unsigned long long > temp((unsigned long long)std::numeric_limits<typename TINPUT::PixelType>::max()+1,0);
			histoGram.at(z) = temp;
		}
		#pragma omp parallel for
		for( long long z=0; z<inputImage_sizez[2]; ++z)
		{
			for(long long x=0; x<inputImage_sizez[0]; ++x)
			{
				for(long long y=0; y<inputImage_sizez[1]; ++y)
				{
					typename TINPUT::PixelType value = inputImageArray[(inputImage_slice_size*z) + (inputImage_row_size*y) + (x)];
					histoGram.at(z).at(value)++;
				}
			}
		}
		std::vector< unsigned long long > result((unsigned long long)std::numeric_limits<typename TINPUT::PixelType>::max()+1,0);
		for( long long num=0; num<(unsigned long long)std::numeric_limits<typename TINPUT::PixelType>::max()+1; ++num)
		{
			unsigned long long maxMax = 0;
			for( long long z=0; z<inputImage_sizez[2]; ++z)
			{
				maxMax = maxMax + histoGram.at(z).at(num);
			}
			result.at(num) = maxMax;
		}
		
		ofstream myfile;
		myfile.open (temp3.c_str());
		myfile << (unsigned long long)std::numeric_limits<typename TINPUT::PixelType>::max() << "\n";
		for( long long num=0; num<(unsigned long long)std::numeric_limits<typename TINPUT::PixelType>::max()+1; ++num)
		{
			myfile << result.at(num) << "\n";
		}
		myfile.close();
	}
}

template<typename TINPUT, typename TOUTPUT >
void ftkMainDarpa::projectImage( typename TINPUT::Pointer inputImage, std::string inputImageName, std::string outputPath, std::string projectOptions, std::string imageType )
{
// 	typename TINPUT::Pointer inputImage = readImage< TINPUT >(inputImageName.c_str());

	std::string tipoImagen;
	int foundType=imageType.find("TIFF");
	if (foundType!=std::string::npos)
	{
		tipoImagen = ".tif";
	}
	foundType=imageType.find("NRRD");
	if (foundType!=std::string::npos)
	{
		tipoImagen = ".nrrd";
	}

	int found=inputImageName.find(".");
	std::string inputImageNameLocal = inputImageName.substr(0,found);
	
	found = inputImageNameLocal.find_last_of("/\\");
	inputImageNameLocal = inputImageNameLocal.substr(found+1);
	
	typename TINPUT::PixelType * inputImageArray = inputImage->GetBufferPointer();
	itk::Size<3> inputImage_sizez = inputImage->GetLargestPossibleRegion().GetSize();
	unsigned long long inputImage_slice_size = inputImage_sizez[1] * inputImage_sizez[0];
	unsigned long long inputImage_row_size = inputImage_sizez[0];	

	// Create z projecImage
	typename TINPUT::Pointer zProjectImage = TINPUT::New();
	typename TINPUT::PointType originz;
	originz[0] = 0; 
	originz[1] = 0;
	originz[2] = 0;
	zProjectImage->SetOrigin( originz );
	typename TINPUT::IndexType startz;
	startz[0] = 0;
	startz[1] = 0;
	startz[2] = 0;
	typename TINPUT::SizeType sizez;
	sizez[0] = inputImage_sizez[0];
	sizez[1] = inputImage_sizez[1];
	sizez[2] = 1;
	typename TINPUT::RegionType regionz;
	regionz.SetSize ( sizez  );
	regionz.SetIndex( startz );
	zProjectImage->SetRegions( regionz );
	zProjectImage->Allocate();
	zProjectImage->FillBuffer(0);
	zProjectImage->Update();
	typename TINPUT::PixelType * zProjectImageArray = zProjectImage->GetBufferPointer();
	
	// Create Y projecImage
	typename TINPUT::Pointer yProjectImage = TINPUT::New();
	typename TINPUT::PointType originy;
	originy[0] = 0; 
	originy[1] = 0;
	originy[2] = 0;
	yProjectImage->SetOrigin( originy );
	typename TINPUT::IndexType starty;
	starty[0] = 0;
	starty[1] = 0;
	starty[2] = 0;
	typename TINPUT::SizeType sizey;
	sizey[0] = inputImage_sizez[0];
	sizey[1] = inputImage_sizez[2];
	sizey[2] = 1;
	typename TINPUT::RegionType regiony;
	regiony.SetSize ( sizey  );
	regiony.SetIndex( starty );
	yProjectImage->SetRegions( regiony );
	yProjectImage->Allocate();
	yProjectImage->FillBuffer(0);
	yProjectImage->Update();
	typename TINPUT::PixelType * yProjectImageArray = yProjectImage->GetBufferPointer();
	
	// Create X projecImage
	typename TINPUT::Pointer xProjectImage = TINPUT::New();
	typename TINPUT::PointType originx;
	originx[0] = 0; 
	originx[1] = 0;
	originx[2] = 0;
	xProjectImage->SetOrigin( originx );
	typename TINPUT::IndexType startx;
	startx[0] = 0;
	startx[1] = 0;
	startx[2] = 0;
	typename TINPUT::SizeType sizex;
	sizex[0] = inputImage_sizez[2];
	sizex[1] = inputImage_sizez[1];
	sizex[2] = 1;
	typename TINPUT::RegionType regionx;
	regionx.SetSize ( sizex  );
	regionx.SetIndex( startx );
	xProjectImage->SetRegions( regionx );
	xProjectImage->Allocate();
	xProjectImage->FillBuffer(0);
	xProjectImage->Update();
	typename TINPUT::PixelType * xProjectImageArray = xProjectImage->GetBufferPointer();
	
// 	std::cout << std::endl << "NOW IS GOING TO PROJECT"<<std::flush;
	
	
	#if _OPENMP < 200805L 
		#pragma omp parallel for //collapse(2) //TEST
	#else
		#pragma omp parallel for collapse(2)
	#endif

	for(long long x=0; x<inputImage_sizez[0]; ++x)
	{
		for(long long y=0; y<inputImage_sizez[1]; ++y)
		{
			typename TINPUT::PixelType max_val = 0;
			for(long long z=0; z<inputImage_sizez[2]; ++z)
			{
				typename TINPUT::PixelType value = inputImageArray[(inputImage_slice_size*z) + (inputImage_row_size*y) + (x)];
				if( max_val < value )
					max_val = value;
			}
			zProjectImageArray[(inputImage_row_size*y) + (x)] = max_val;
		}
	}
	
	#if _OPENMP < 200805L 
		#pragma omp parallel for //collapse(2) //TEST
	#else
		#pragma omp parallel for collapse(2)
	#endif

	for( long long z=0; z<inputImage_sizez[2]; ++z)
	{
		for(long long x=0; x<inputImage_sizez[0]; ++x)
		{
			typename TINPUT::PixelType max_val = 0;
			for(long long y=0; y<inputImage_sizez[1]; ++y)
			{
				typename TINPUT::PixelType value = inputImageArray[(inputImage_slice_size*z) + (inputImage_row_size*y) + (x)];
				if( max_val < value )
					max_val = value;
			}
			yProjectImageArray[(inputImage_sizez[0]*z) + (x)] = max_val;
		}
	}
	
	#if _OPENMP < 200805L 
		#pragma omp parallel for //collapse(2) //TEST
	#else
		#pragma omp parallel for collapse(2)
	#endif

	for( long long z=0; z<inputImage_sizez[2]; ++z)
	{
		for( long long y=0; y<inputImage_sizez[1]; ++y)
		{
			typename TINPUT::PixelType max_val = 0;
			for( long long x=0; x<inputImage_sizez[0]; ++x)
			{
				typename TINPUT::PixelType value = inputImageArray[(inputImage_slice_size*z) + (inputImage_row_size*y) + (x)];
				if( max_val < value )
					max_val = value;
			}
			xProjectImageArray[(inputImage_sizez[2]*y) + (z)] = max_val;
		}
	}


	int foundPro=projectOptions.find("ORG");
	if (foundPro!=std::string::npos)
	{
		std::string temp3 = outputPath + "/" + inputImageNameLocal + "zPro_Z"+tipoImagen;
		writeImage< TINPUT >(zProjectImage,temp3.c_str());
	
		std::string temp4 = outputPath + "/" + inputImageNameLocal + "zPro_Y"+tipoImagen;
		writeImage< TINPUT >(yProjectImage,temp4.c_str());
		
		std::string temp5 = outputPath + "/" + inputImageNameLocal + "zPro_X"+tipoImagen;
		writeImage< TINPUT >(xProjectImage,temp5.c_str());
	
	}
	int foundRes=projectOptions.find("RES");
	if (foundRes!=std::string::npos)
	{
		std::string temp3a = outputPath + "/" + inputImageNameLocal + "zPro_Z_Re"+tipoImagen;
		writeImageRescaled< TINPUT, TOUTPUT >(zProjectImage,temp3a.c_str());
	
		std::string temp4a = outputPath + "/" + inputImageNameLocal + "zPro_Y_Re"+tipoImagen;
		writeImageRescaled< TINPUT, TOUTPUT >(yProjectImage,temp4a.c_str());
		
		std::string temp5a = outputPath + "/" + inputImageNameLocal + "zPro_X_Re"+tipoImagen;
		writeImageRescaled< TINPUT, TOUTPUT >(xProjectImage,temp5a.c_str());
	}
		
	int foundBin=projectOptions.find("BIN");
	if (foundBin!=std::string::npos)
	{
		#if _OPENMP < 200805L 
			#pragma omp parallel for //collapse(2) //TEST
		#else
			#pragma omp parallel for collapse(2)
		#endif

		for( long long x=0; x<inputImage_sizez[0]; ++x)
		{
			for(long long y=0; y<inputImage_sizez[1]; ++y)
			{
				if( zProjectImageArray[(inputImage_row_size*y) + (x)] != 0 )
					zProjectImageArray[(inputImage_row_size*y) + (x)] = std::numeric_limits<typename TINPUT::PixelType>::max();
			}
		}
		#if _OPENMP < 200805L 
			#pragma omp parallel for //collapse(2) //TEST
		#else
			#pragma omp parallel for collapse(2)
		#endif
		for(long long z=0; z<inputImage_sizez[2]; ++z)
		{
			for(long long x=0; x<inputImage_sizez[0]; ++x)
			{
				if( yProjectImageArray[(inputImage_sizez[0]*z) + (x)] != 0 )
					yProjectImageArray[(inputImage_sizez[0]*z) + (x)] = std::numeric_limits<typename TINPUT::PixelType>::max();
			}
		}
		
		#if _OPENMP < 200805L 
			#pragma omp parallel for //collapse(2) //TEST
		#else
			#pragma omp parallel for collapse(2)
		#endif
		for( long long z=0; z<inputImage_sizez[2]; ++z)
		{
			for( long long y=0; y<inputImage_sizez[1]; ++y)
			{
				if( xProjectImageArray[(inputImage_sizez[2]*y) + (z)] != 0 )
					xProjectImageArray[(inputImage_sizez[2]*y) + (z)] = std::numeric_limits<typename TINPUT::PixelType>::max();
			}
		}

		std::string temp3b = outputPath + "/" + inputImageNameLocal + "zPro_Z_Bin"+tipoImagen;
		writeImage< TINPUT >(zProjectImage,temp3b.c_str());
		std::string temp4b = outputPath + "/" + inputImageNameLocal + "zPro_Y_Bin"+tipoImagen;
		writeImage< TINPUT >(yProjectImage,temp4b.c_str());
		std::string temp5b = outputPath + "/" + inputImageNameLocal + "zPro_X_Bin"+tipoImagen;
		writeImage< TINPUT >(xProjectImage,temp5b.c_str());
	}
	
	int foundHisto=projectOptions.find("HISTO");
	if (foundHisto!=std::string::npos)
	{
		std::string temp3 = outputPath + "/" + inputImageNameLocal + "zHisto.txt";
		std::vector< std::vector< unsigned long long > > histoGram(inputImage_sizez[2]);
		for( long long z=0; z<inputImage_sizez[2]; ++z)
		{
			std::vector< unsigned long long > temp((unsigned long long)std::numeric_limits<typename TINPUT::PixelType>::max()+1,0);
			histoGram.at(z) = temp;
		}
		#pragma omp parallel for
		for( long long z=0; z<inputImage_sizez[2]; ++z)
		{
			for( long long x=0; x<inputImage_sizez[0]; ++x)
			{
				for( long long y=0; y<inputImage_sizez[1]; ++y)
				{
					typename TINPUT::PixelType value = inputImageArray[(inputImage_slice_size*z) + (inputImage_row_size*y) + (x)];
					histoGram.at(z).at(value)++;
				}
			}
		}
		std::vector< unsigned long long > result((unsigned long long)std::numeric_limits<typename TINPUT::PixelType>::max()+1,0);
		for( long long num=0; num<(unsigned long long)std::numeric_limits<typename TINPUT::PixelType>::max()+1; ++num)
		{
			unsigned long long maxMax = 0;
			for( long long z=0; z<inputImage_sizez[2]; ++z)
			{
				maxMax = maxMax + histoGram.at(z).at(num);
			}
			result.at(num) = maxMax;
		}
		
		ofstream myfile;
		myfile.open (temp3.c_str());
		myfile << std::numeric_limits<typename TINPUT::PixelType>::max() << "\n";
		for( long long num=0; num<(unsigned long long)std::numeric_limits<typename TINPUT::PixelType>::max()+1; ++num)
		{
			myfile << result.at(num) << "\n";
		}
		myfile.close();
	}
}

template<typename TINPUT, typename TOUTPUT >
void ftkMainDarpa::computeDistMap( std::string inputImageName, std::string outputPath, std::string imageType )
{
	std::string tipoImagen;
	int foundType=imageType.find("TIFF");
	if (foundType!=std::string::npos)
	{
		tipoImagen = ".tif";
	}
	foundType=imageType.find("NRRD");
	if (foundType!=std::string::npos)
	{
		tipoImagen = ".nrrd";
	}

	//int found=inputImageName.find(".");
	//std::string inputImageNameLocal = inputImageName.substr(0,found);
	//
	//found = inputImageNameLocal.find_last_of("/\\");
	//inputImageNameLocal = inputImageNameLocal.substr(found+1);

	typename TINPUT::Pointer inputImage = readImage< TINPUT >(inputImageName.c_str());
	//std::string temp1a = outputPath + "/" + inputImageNameLocal + "_dist_map" + tipoImagen;
	
	typedef itk::BinaryThresholdImageFilter<TINPUT, rawImageType_8bit> ThresholdFilterType;
	typename ThresholdFilterType::Pointer threshold_filter = ThresholdFilterType::New();
	threshold_filter->SetLowerThreshold(1);
	threshold_filter->SetInsideValue(255);
	threshold_filter->SetOutsideValue(0);
	threshold_filter->SetInput(inputImage);
	//threshold_filter->Update();

	typedef itk::SignedMaurerDistanceMapImageFilter<rawImageType_8bit, TOUTPUT> SignedMaurerDistanceMapImageFilterType;
	typename SignedMaurerDistanceMapImageFilterType::Pointer MaurerFilter = SignedMaurerDistanceMapImageFilterType::New();
	MaurerFilter->SetInput(threshold_filter->GetOutput());
	MaurerFilter->SetSquaredDistance(false);
	MaurerFilter->SetUseImageSpacing(false);
	MaurerFilter->SetInsideIsPositive(false);
	MaurerFilter->Update();

	writeImage< TOUTPUT >(MaurerFilter->GetOutput(),outputPath.c_str());

}


template<typename TINPUT, typename TOUTPUT >
void ftkMainDarpa::computeMedianFilter( std::string inputImageName, std::string outputImageName, std::string imageType )
{
	std::string tipoImagen;
	int foundType=imageType.find("TIFF");
	if (foundType!=std::string::npos)
	{
		tipoImagen = ".tif";
	}
	foundType=imageType.find("NRRD");
	if (foundType!=std::string::npos)
	{
		tipoImagen = ".nrrd";
	}

	//int found=inputImageName.find(".");
	//std::string inputImageNameLocal = inputImageName.substr(0,found);
	//
	//found = inputImageNameLocal.find_last_of("/\\");
	//inputImageNameLocal = inputImageNameLocal.substr(found+1);

	typename TINPUT::Pointer inputImage = readImage< TINPUT >(inputImageName.c_str());
	//std::string temp1a = outputPath + "/" + inputImageNameLocal + "_dist_map" + tipoImagen;
	
	typedef itk::MedianImageFilter<TINPUT, TOUTPUT> MedianFilterType;
	typename MedianFilterType::Pointer median_filter = MedianFilterType::New();
	typename MedianFilterType::InputSizeType radius;
	radius.Fill(2);
   
	median_filter->SetRadius(radius);
	median_filter->SetInput(inputImage);
	median_filter->Update();
	
	std::cout<<"HERE";
	std::cout<<"HERE";
	std::cout<<"HERE";

	writeImage< TOUTPUT >(median_filter->GetOutput(),outputImageName.c_str());

}
template<typename TINPUT >
void ftkMainDarpa::cropImageDarpa( std::string imageInput, std::string tableInput, std::string coordinate )
{
	
	std::cout<<"reading the file 2"<<std::endl;

	/// read region of interest ///
	std::ifstream coordinate_file;
	coordinate_file.open( coordinate.c_str() );
	itk::SizeValueType Xmin;
	coordinate_file >> Xmin;

	itk::SizeValueType Xsize;
	coordinate_file >> Xsize;
	itk::SizeValueType Ymin;
	coordinate_file >> Ymin;
	itk::SizeValueType Ysize;
	coordinate_file >> Ysize;
	itk::SizeValueType Zmin;
	coordinate_file >> Zmin;
	itk::SizeValueType Zsize;
	coordinate_file >> Zsize;

	std::cout<<Xmin<<std::endl;
	std::cout<<Xsize<<std::endl;
	std::cout<<Ymin<<std::endl;
	std::cout<<Ysize<<std::endl;
	std::cout<<Zmin<<std::endl;
	std::cout<<Zsize<<std::endl;
	//scanf("%d");

	/// read table ///
	vtkSmartPointer< vtkTable > feature_table = ftk::LoadTable(tableInput);

	/// read label image ///
	typename TINPUT::Pointer labelImage = readImage< TINPUT >(imageInput.c_str());

	/// start cleaning the labels ///

	/// crop the label image ///
	itk::Image<unsigned short,3>::Pointer clean_labelImage = itk::Image<unsigned short,3>::New();
	//itk::Image<unsigned short,3>
	itk::Size<3> im_size;
	im_size[0] = Xsize;
	im_size[1] = Ysize;
	im_size[2] = Zsize;

	itk::Index<3> im_index;
	im_index[0] = 0;
	im_index[1] = 0;
	im_index[2] = 0;

	std::cout<<"finished reading image, started allocation..."<<std::endl << std::flush;
	typename TINPUT::RegionType region;
	region.SetSize( im_size );
	region.SetIndex( im_index );
	clean_labelImage->SetRegions( region );
	clean_labelImage->Allocate();
	clean_labelImage->FillBuffer(0);
	try
	{
		clean_labelImage->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "ExceptionObject caught!" <<std::endl;
		std::cerr << err << std::endl;
	}



	itk::SizeValueType Xmax =  Xmin + Xsize;
	itk::SizeValueType Ymax =  Ymin + Ysize;
	itk::SizeValueType Zmax =  Zmin + Zsize;
	unsigned int id = -1;
	std::map<unsigned int, bool> classMap;
	std::map<unsigned int, unsigned int> idMap;
	unsigned short new_id =1;
	std::cout<<"starting.."<<std::endl << std::flush;

	for(int row=0; row<(int)feature_table->GetNumberOfRows(); ++row)
	{
		itk::SizeValueType x_centroid = feature_table->GetValue(row,1).ToUnsignedInt();
		itk::SizeValueType y_centroid = feature_table->GetValue(row,2).ToUnsignedInt();
		itk::SizeValueType z_centroid = feature_table->GetValue(row,3).ToUnsignedInt();
		unsigned int id = feature_table->GetValue(row,0).ToUnsignedInt();


		if((x_centroid > Xmin && x_centroid < Xmax) && (y_centroid > Ymin && y_centroid < Ymax) && (z_centroid >=Zmin && z_centroid < Zmax) )
		{
			classMap[id] = true;
			idMap[id] = new_id;
		///	std::cout<<new_id<<std::endl;
			new_id++;
		}
		else
		{
			classMap[id] = false;
		}

	}
	std::cout<<"finished looking for inside cells.."<<std::endl << std::flush;
	itk::Image<unsigned short,3>::PixelType * clean_labelArray = clean_labelImage->GetBufferPointer();
	typename TINPUT::PixelType * labelArray = labelImage->GetBufferPointer();

	itk::SizeValueType sizeXYsmall = im_size[1] * im_size[0];
	itk::SizeValueType sizeXsmall = im_size[0];


	itk::Size<3> montage_size = labelImage->GetBufferedRegion().GetSize();
	itk::SizeValueType sizeXYlarge= montage_size[1] * montage_size[0];
	itk::SizeValueType sizeXlarge = montage_size[0];

#if _OPENMP < 200805L
	#pragma omp parallel
#else
	#pragma omp parallel for collapse(3)
#endif
	for(int i=0; i<im_size[2]; ++i)
	{
		//itk::SizeValueType sm01 = i*sizeXYsmall;
		//itk::SizeValueType la01 = (i+Xmin)*sizeXYlarge;

		for(int j=0; j<im_size[1]; ++j)
		{
			//itk::SizeValueType  sm02 = sm01 +(j*sizeXsmall);
			//itk::SizeValueType  la02 = la01 +(j+Ymin)*sizeXlarge;

			for(int k=0; k<im_size[0]; ++k)
			{
		
				//itk::SizeValueType small_offset = sm02 + k;
				//itk::SizeValueType large_offset = la02 +(k+Zmin);

				itk::SizeValueType small_offset = (i*sizeXYsmall)+(j*sizeXsmall)+k;
				itk::SizeValueType large_offset = ((i+Zmin)*sizeXYlarge)+((j+Ymin)*sizeXlarge)+(k+Xmin);
				///std::cout<<labelArray[large_offset]<<std::endl;
				if( labelArray[large_offset] != 0 )
					if(classMap[labelArray[large_offset]])
					{
//						clean_labelArray[small_offset] = labelArray[large_offset];
						clean_labelArray[small_offset] = idMap[labelArray[large_offset]];
						//std::cout<<idMap[labelArray[large_offset]]<<std::endl;
					}
			}
		}
	}	
	std::cout<<"finished copying the image.."<<std::endl;
	/// clean up the table ///
	for( unsigned int row = 0; row < feature_table->GetNumberOfRows(); ++row)
	{
		if(!classMap[feature_table->GetValue(row,0).ToUnsignedInt()])
		{
			feature_table->RemoveRow(row);
			--row;
		}
		else
		{
			unsigned int newId = idMap[feature_table->GetValue(row,0).ToUnsignedInt()];
			feature_table->SetValue(row,0,newId);

			itk::SizeValueType x_centroid = feature_table->GetValue(row,1).ToUnsignedInt();
			itk::SizeValueType y_centroid = feature_table->GetValue(row,2).ToUnsignedInt();
			itk::SizeValueType z_centroid = feature_table->GetValue(row,3).ToUnsignedInt();

			feature_table->SetValue(row,1,x_centroid-Xmin);
			feature_table->SetValue(row,2,y_centroid-Ymin);
			feature_table->SetValue(row,3,z_centroid-Zmin);


		}
	}
	std::cout<<"finished relabeling  the table.."<<std::endl;

	std::stringstream out_x;
  std::stringstream out_y;
	std::stringstream out_z;
	out_x<<Xmin;
	out_y<<Ymin;
  out_z<<Zmin;
	

	std::string curr_path = ftk::GetFilePath(imageInput);
	std::string imfname = ftk::GetFilenameFromFullPath(imageInput);
	std::string imoutfname = curr_path + "/cropped_"+out_x.str()+"_"+out_y.str()+"_"+out_z.str()+"_"+imfname;
	writeImage< itk::Image<unsigned short,3> >(clean_labelImage,imoutfname.c_str());


	std::string tabfname = ftk::GetFilenameFromFullPath(tableInput);
	ftk::SaveTable(curr_path + "/cropped_"+out_x.str()+"_"+out_y.str()+"_"+out_z.str()+"_"+tabfname, feature_table);


}
