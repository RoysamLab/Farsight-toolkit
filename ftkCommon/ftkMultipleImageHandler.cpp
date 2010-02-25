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

void MultipleImageHandler::SeriesToBlocks(std::string seriesFormat, int startIndex, int endIndex, int dx, int dy, int dz)
{
	typedef itk::NumericSeriesFileNames NameGeneratorType;
	NameGeneratorType::Pointer nameGenerator = NameGeneratorType::New();
	nameGenerator->SetSeriesFormat( seriesFormat );
	nameGenerator->SetStartIndex( startIndex );
	nameGenerator->SetEndIndex( endIndex );
	nameGenerator->SetIncrementIndex( 1 );

	this->SeriesToBlocks( nameGenerator->GetFileNames(), dx, dy, dz );
}


//This function will take a list of input files - each 2D, starting with z=0.
//The ouput is a new set of 3D image files of the image made into blocks
//dx,dy,dz specify the number of divisions in the x,y,z directions.
void MultipleImageHandler::SeriesToBlocks(StrVector inFiles, int dx, int dy, int dz)
{
	int zSize = (int)inFiles.size();
	if(zSize == 0)
		return;

	// Find out the pixel type of the image in file
	itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO( inFiles.at(0).c_str(), itk::ImageIOFactory::ReadMode );
	if( !imageIO )
    {
		std::cerr << "NO IMAGEIO WAS FOUND" << std::endl;
		return ;
    }

	// Now that we found the appropriate ImageIO class, ask it to 
	// read the meta data from the image file.
	imageIO->SetFileName( inFiles.at(0) );
	imageIO->ReadImageInformation();

	int numComponents =  imageIO->GetNumberOfComponents();
	itk::ImageIOBase::IOComponentType dataType = imageIO->GetComponentType();
	int numDimensions = imageIO->GetNumberOfDimensions();

	int xSize = imageIO->GetDimensions(0);
	int ySize = imageIO->GetDimensions(1);

	imageIO = 0; //free the IO

	PairVector xPairs = this->CreateOutputRegions(xSize, dx);
	PairVector yPairs = this->CreateOutputRegions(ySize, dy);
	PairVector zPairs = this->CreateOutputRegions(zSize, dz);

	for(int zp=0; zp<(int)zPairs.size(); ++zp)
	{
		for(int yp=0; yp<(int)yPairs.size(); ++yp)
		{
			for(int xp=0; xp<(int)xPairs.size(); ++xp)
			{
				std::string filename;
				filename = "b" + NumToString(zp) + NumToString(yp) + NumToString(xp) + ".tif";
				ExtractBlock(inFiles,xPairs.at(xp),yPairs.at(yp),zPairs.at(zp), filename );
			}
		}
	}
}


void MultipleImageHandler::ExtractBlock(StrVector inFiles,  PairType xPair, PairType yPair, PairType zPair, std::string fname)
{
	//Create output region:
	UCharImageType3D::RegionType block;
	block.SetIndex( 0, xPair.first );
	block.SetSize( 0, xPair.second );
	block.SetIndex( 1, yPair.first );
	block.SetSize( 1, yPair.second );
	block.SetIndex( 2, zPair.first );
	block.SetSize( 2, zPair.second );

	//Create the output image
	typedef itk::ImageSeriesReader< UCharImageType3D > SeriesReaderType;
	SeriesReaderType::Pointer reader = SeriesReaderType::New();
	reader->SetFileNames(inFiles);

	typedef itk::ExtractImageFilter< UCharImageType3D, UCharImageType3D > ExtractFilterType;
	ExtractFilterType::Pointer extract = ExtractFilterType::New();
	extract->SetInput( reader->GetOutput() );
	extract->SetExtractionRegion( block );

	typedef itk::ImageFileWriter< UCharImageType3D > ImageWriterType;
	ImageWriterType::Pointer writer = ImageWriterType::New();
	writer->SetFileName( fname.c_str() );
	writer->SetInput( extract->GetOutput() );

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











