/*=========================================================================
Copyright 2009 Rensselaer Polytechnic Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
=========================================================================*/
#include "ftkMultipleImageHandler.h"


namespace ftk
{

//constructor
MultipleImageHandler::MultipleImageHandler()
{
	outputDirectory = "";
	outputBaseString = "b";
}

void MultipleImageHandler::SeriesToBlocks(std::string seriesFormat, int startIndex, int endIndex, int dx, int dy, int dz, int color)
{
	typedef itk::NumericSeriesFileNames NameGeneratorType;
	NameGeneratorType::Pointer nameGenerator = NameGeneratorType::New();
	nameGenerator->SetSeriesFormat( seriesFormat );
	nameGenerator->SetStartIndex( startIndex );
	nameGenerator->SetEndIndex( endIndex );
	nameGenerator->SetIncrementIndex( 1 );

	this->SeriesToBlocks( nameGenerator->GetFileNames(), dx, dy, dz, color );
}


//This function will take a list of input files - each 2D, starting with z=0.
//The ouput is a new set of 3D image files of the image made into blocks
//dx,dy,dz specify the number of divisions in the x,y,z directions.
void MultipleImageHandler::SeriesToBlocks(StrVector inFiles, int dx, int dy, int dz, int color)
{
	int zSize = (int)inFiles.size();
	if(zSize == 0)
		return;

	// Find out the pixel type of the image in file
	itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO( inFiles.at(0).c_str(), itk::ImageIOFactory::ReadMode );
	if( !imageIO )
    {
		std::cerr << "NO IMAGEIO WAS FOUND" << std::endl;
		std::cerr << inFiles.at(0) << std::endl;
		return;
    }

	// Now that we found the appropriate ImageIO class, ask it to 
	// read the meta data from the image file.
	imageIO->SetFileName( inFiles.at(0) );
	imageIO->ReadImageInformation();

	int numDimensions = imageIO->GetNumberOfDimensions();
	if(numDimensions > 2)
		return;

	int numComponents =  imageIO->GetNumberOfComponents();
	if( numComponents != 3 && numComponents != 1 )
		return;

	itk::ImageIOBase::IOComponentType dataType = imageIO->GetComponentType();
	if( dataType != itk::ImageIOBase::UCHAR && dataType != itk::ImageIOBase::USHORT )
		return;

	int xSize = imageIO->GetDimensions(0);
	int ySize = imageIO->GetDimensions(1);

	imageIO = 0; //free the IO

	PairVector xPairs = this->CreateOutputRegions(xSize, dx);
	PairVector yPairs = this->CreateOutputRegions(ySize, dy);
	PairVector zPairs = this->CreateOutputRegions(zSize, dz);

	std::cerr << "SeriesToBlocks..." << std::endl;

	ProjectManager * projF = new ProjectManager();

	int totalBlocks = dx*dy*dz;
	int counter = 1;

	for(int zp=0; zp<(int)zPairs.size(); ++zp)
	{
		for(int yp=0; yp<(int)yPairs.size(); ++yp)
		{
			for(int xp=0; xp<(int)xPairs.size(); ++xp)
			{
				std::string filename = outputBaseString;
				if(totalBlocks > 1)
					filename += "_" + NumToString(counter);
				filename += ".tif";

				UCharImageType3D::RegionType region = CreateRegion( xPairs.at(xp), yPairs.at(yp), zPairs.at(zp) );

				projF->addFile(filename, "image", region.GetIndex(0), region.GetIndex(1), region.GetIndex(2));

				std::cerr << "Processing " << counter++ << " of " << totalBlocks << "\r";

				if(numComponents == 3) //must be RGB:
				{
					ExtractRegionColor(inFiles, region, color, filename);
				}
				else if(dataType == itk::ImageIOBase::UCHAR)	//1 component 8-bit grayscale
				{
					ExtractRegion(inFiles, region, false, filename );
				}
				else											//1 component 16-bit grayscale
				{
					ExtractRegion(inFiles, region, true, filename );
				}
			}
		}
	}

	std::string outFname = outputDirectory + "bFiles.xml";
	if( projF->size() > 1 )
		projF->writeProject(outFname.c_str());
	delete projF;

	std::cerr << std::endl << "...Done" << std::endl;
}

MultipleImageHandler::UCharImageType3D::Pointer MultipleImageHandler::ExtractRegion
	(StrVector inFiles, UCharImageType3D::RegionType region, bool rescale, std::string fname)
{
	UCharImageType3D::Pointer img = NULL;

	if(rescale)
	{
		typedef itk::ImageSeriesReader< UShortImageType3D > SeriesReaderType;
		SeriesReaderType::Pointer reader = SeriesReaderType::New();
		reader->SetFileNames(inFiles);

		typedef itk::ExtractImageFilter< UShortImageType3D, UShortImageType3D > ExtractFilterType;
		ExtractFilterType::Pointer extract = ExtractFilterType::New();
		extract->SetInput( reader->GetOutput() );
		extract->SetExtractionRegion( region );

		typedef itk::RescaleIntensityImageFilter< UShortImageType3D, UCharImageType3D > RescaleType;
		RescaleType::Pointer rescale = RescaleType::New();
		rescale->SetInput( extract->GetOutput() );
		rescale->SetOutputMinimum(0);
		rescale->SetOutputMaximum(255);

		try
		{
			rescale->Update();
		}
		catch( itk::ExceptionObject & err )
		{
			std::cerr << "EXTRACT FAILED: " << err << std::endl;
			return NULL;
		}
		
		img = rescale->GetOutput();
	}
	else
	{
		typedef itk::ImageSeriesReader< UCharImageType3D > SeriesReaderType;
		SeriesReaderType::Pointer reader = SeriesReaderType::New();
		reader->SetFileNames(inFiles);

		typedef itk::ExtractImageFilter< UCharImageType3D, UCharImageType3D > ExtractFilterType;
		ExtractFilterType::Pointer extract = ExtractFilterType::New();
		extract->SetInput( reader->GetOutput() );
		extract->SetExtractionRegion( region );

		try
		{
			extract->Update();
		}
		catch( itk::ExceptionObject & err )
		{
			std::cerr << "EXTRACT FAILED: " << err << std::endl;
			return NULL;
		}

		img = extract->GetOutput();
	}

	if(fname != "")
	{
		//Setup the writer
		typedef itk::ImageFileWriter< UCharImageType3D > ImageWriterType;
		ImageWriterType::Pointer writer = ImageWriterType::New();
		writer->SetFileName( fname.c_str() );
		writer->SetInput( img );
		//This line has no effect
		//writer->SetNumberOfStreamDivisions( 10 );
		try
		{
			writer->Update();
		}
		catch( itk::ExceptionObject & err )
		{
			std::cerr << "WRITER FAILED: " << err << std::endl;
		}
	}

	return img;
}

//Extract a region from a series of color images:
MultipleImageHandler::UCharImageType3D::Pointer MultipleImageHandler::ExtractRegionColor
	(StrVector inFiles, UCharImageType3D::RegionType region, int color, std::string fname)
{
		typedef itk::ImageSeriesReader< RGBImageType3D > SeriesReaderType;
		SeriesReaderType::Pointer reader = SeriesReaderType::New();
		reader->SetFileNames(inFiles);

		typedef itk::ExtractImageFilter< RGBImageType3D, RGBImageType3D > ExtractFilterType;
		ExtractFilterType::Pointer extract = ExtractFilterType::New();
		extract->SetInput( reader->GetOutput() );
		extract->SetExtractionRegion( region );
		try
		{
			extract->Update();
		}
		catch( itk::ExceptionObject & err )
		{
			std::cerr << "EXTRACT FAILED: " << err << std::endl;
			return NULL;
		}

		RGBImageType3D::Pointer extracted = extract->GetOutput();
		UCharImageType3D::Pointer img;

		if(color==0 || color==1 || color==2)
		{
			int size1 = (int)extracted->GetLargestPossibleRegion().GetSize()[0];
			int size2 = (int)extracted->GetLargestPossibleRegion().GetSize()[1];
			int size3 = (int)extracted->GetLargestPossibleRegion().GetSize()[2];

			std::cout << "Extracting color " << color << " from Image of size " << size1 << " " << size2 << " " << size3 << std::endl;

			img = UCharImageType3D::New();
			UCharImageType3D::PointType origin;
			origin[0] = 0;
			origin[1] = 0;
			origin[2] = 0;
			img->SetOrigin( origin );
			UCharImageType3D::IndexType start = {{ 0,0,0 }};
			UCharImageType3D::SizeType  size = {{ size1, size2, size3 }};
			UCharImageType3D::RegionType region;
			region.SetSize( size );
			region.SetIndex( start );
			img->SetRegions( region ); 
			img->Allocate();

			typedef itk::ImageRegionIterator< RGBImageType3D > rgbIteratorType;
			rgbIteratorType iteratorRGB ( extracted, extracted->GetLargestPossibleRegion() );

			typedef itk::ImageRegionIterator< UCharImageType3D > myIteratorType;
			myIteratorType iteratorIn( img, img->GetLargestPossibleRegion() );

			for( iteratorRGB.GoToBegin(), iteratorIn.GoToBegin();
				!iteratorRGB.IsAtEnd(), !iteratorIn.IsAtEnd();
				++iteratorRGB, ++iteratorIn )
			{ 
				RGBPixelType p_rgb = iteratorRGB.Get();
				iteratorIn.Set( p_rgb[color] );
			}

			std::cout << "Done color extract\n";
		}
		else
		{
			typedef itk::RGBToLuminanceImageFilter< RGBImageType3D, UCharImageType3D > ConvertFilterType;
			ConvertFilterType::Pointer convert = ConvertFilterType::New();
			convert->SetInput( extracted );
			try
			{
				convert->Update();
			}
			catch( itk::ExceptionObject & err )
			{
				std::cerr << "EXTRACT FAILED: " << err << std::endl;
				return NULL;
			}
			img = convert->GetOutput();
		}

	if(fname != "")
	{
		//Setup the writer
		typedef itk::ImageFileWriter< UCharImageType3D > ImageWriterType;
		ImageWriterType::Pointer writer = ImageWriterType::New();
		writer->SetFileName( fname.c_str() );
		writer->SetInput( img );
		//This line has no effect
		//writer->SetNumberOfStreamDivisions( 10 );
		try
		{
			writer->Update();
		}
		catch( itk::ExceptionObject & err )
		{
			std::cerr << "WRITER FAILED: " << err << std::endl;
		}
	}

	std::cout << "Done ExtractRegionColor\n";

	return img;

}


MultipleImageHandler::UCharImageType3D::RegionType MultipleImageHandler::CreateRegion(PairType xPair, PairType yPair, PairType zPair)
{
	//Create output region:
	UCharImageType3D::RegionType block;

	block.SetIndex( 0, xPair.first );
	block.SetSize( 0, xPair.second );
	block.SetIndex( 1, yPair.first );
	block.SetSize( 1, yPair.second );
	block.SetIndex( 2, zPair.first );
	block.SetSize( 2, zPair.second );

	return block;
}

MultipleImageHandler::PairVector MultipleImageHandler::CreateOutputRegions(int size, int divs)
{
	int maxPerBlock = (int)ceil( (double)size / (double)divs );

	PairVector pairs;	//includes startZ and sizeZ
	for( int z=0; z<divs; z++ )
	{
		int start = ceil( (double)z / (double)divs * (double)size );
		int end = start + maxPerBlock;
		if( end  > size) end = size;
		int size = end-start;

		PairType p(start, size);
		pairs.push_back(p);
	}

	return pairs;
}

//NOTE: the series projection function does not seem 
//to use multiple threads and it does not only load up 1 slice at a time.
//If we need this capability on large images we should re-implement it to 
//just load up on image at a time and compare it to the output image pixel by pixel
MultipleImageHandler::UCharImageType2D::Pointer MultipleImageHandler::SeriesProjection(std::string seriesFormat, int startIndex, int endIndex, std::string outName)
{
	typedef itk::NumericSeriesFileNames NameGeneratorType;
	NameGeneratorType::Pointer nameGenerator = NameGeneratorType::New();
	nameGenerator->SetSeriesFormat( seriesFormat );
	nameGenerator->SetStartIndex( startIndex );
	nameGenerator->SetEndIndex( endIndex );
	nameGenerator->SetIncrementIndex( 1 );

	return this->SeriesProjection( nameGenerator->GetFileNames(), outName );
}

MultipleImageHandler::UCharImageType2D::Pointer MultipleImageHandler::SeriesProjection(StrVector inFiles, std::string outName)
{
	typedef itk::ImageSeriesReader< UCharImageType3D > SeriesReaderType;
	SeriesReaderType::Pointer reader = SeriesReaderType::New();
	reader->SetFileNames(inFiles);

	typedef itk::MaximumProjectionImageFilter< UCharImageType3D, UCharImageType2D > ProjectionType;
	ProjectionType::Pointer projection = ProjectionType::New();
	projection->SetInput( reader->GetOutput() );
	projection->SetProjectionDimension( 2 );
	projection->SetNumberOfThreads( 4 );

	try
	{
		projection->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "PROJECTION FAILED: " << err << std::endl;
		return NULL;
	}

	UCharImageType2D::Pointer img = projection->GetOutput();

	if(outName != "")
	{
		//Setup the writer
		typedef itk::ImageFileWriter< UCharImageType2D > ImageWriterType;
		ImageWriterType::Pointer writer = ImageWriterType::New();
		writer->SetFileName( outName.c_str() );
		writer->SetInput( img );
		//This line has no effect for most image types
		//writer->SetNumberOfStreamDivisions( 10 );
		try
		{
			writer->Update();
		}
		catch( itk::ExceptionObject & err )
		{
			std::cerr << "WRITER FAILED: " << err << std::endl;
		}
	}

	return img;
}

MultipleImageHandler::UCharImageType2D::Pointer MultipleImageHandler::ImageProjection(std::string inFile, std::string outName)
{
	typedef itk::ImageFileReader< UCharImageType3D > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(inFile);

	typedef itk::MaximumProjectionImageFilter< UCharImageType3D, UCharImageType2D > ProjectionType;
	ProjectionType::Pointer projection = ProjectionType::New();
	projection->SetInput( reader->GetOutput() );
	projection->SetProjectionDimension( 2 );
	projection->SetNumberOfThreads( 4 );

	try
	{
		projection->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "PROJECTION FAILED: " << err << std::endl;
		return NULL;
	}

	UCharImageType2D::Pointer img = projection->GetOutput();

	if(outName != "")
	{
		//Setup the writer
		typedef itk::ImageFileWriter< UCharImageType2D > ImageWriterType;
		ImageWriterType::Pointer writer = ImageWriterType::New();
		writer->SetFileName( outName.c_str() );
		writer->SetInput( img );
		//This line has no effect for most image types
		//writer->SetNumberOfStreamDivisions( 10 );
		try
		{
			writer->Update();
		}
		catch( itk::ExceptionObject & err )
		{
			std::cerr << "WRITER FAILED: " << err << std::endl;
		}
	}

	return img;
}


MultipleImageHandler::UCharImageType2D::Pointer MultipleImageHandler::ReadImage2D(std::string filename)
{

	typedef itk::ImageFileReader< UShortImageType2D > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(filename);

	typedef itk::RescaleIntensityImageFilter< UShortImageType2D, UCharImageType2D > RescaleType;
	RescaleType::Pointer rescale = RescaleType::New();
	rescale->SetInput( reader->GetOutput() );
	rescale->SetOutputMinimum(0);
	rescale->SetOutputMaximum(255);

	try
	{
		rescale->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "READER FAILED: " << err << std::endl ;
		return NULL;
	}

	return rescale->GetOutput();

}

} //end namespace











