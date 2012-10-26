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

#ifndef _ftkImage_txx
#define _ftkImage_txx
//#include "ftkImage.h"	//This .txx file is included from ftkImage.h!!!
//USAGE: see ftkImage.h

#include <itkFixedArray.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>
#include <itkImportImageContainer.h>

namespace ftk
{

//I NEVER RELEASE MEMORY FOR SLICES!!!
template <typename pixelType> pixelType * Image::GetSlicePtr(itk::SizeValueType T, itk::SizeValueType CH, itk::SizeValueType Z, PtrMode mode)
{
	if( T >= m_Info.numTSlices || CH >= m_Info.numChannels || Z >= m_Info.numZSlices )
		return NULL;
		
	if( !IsMatch<pixelType>(m_Info.dataType) )
		return NULL;

	if( mode == RELEASE_CONTROL)
		return NULL;

	itk::SizeValueType numPix = (m_Info.numColumns)*(m_Info.numRows);
	pixelType * stack = static_cast<pixelType *>(imageDataPtrs[T][CH].mem);
	pixelType * slice = stack + Z*numPix;
	pixelType * mem;
	if( mode == DEEP_COPY)
	{
		//mem = new pixelType[m_Info.numColumns*m_Info.numRows];
		mem = (pixelType*)malloc(numPix * m_Info.bytesPerPix);
		if(mem == NULL)
			return (pixelType *)NULL;
		memcpy(mem, slice, numPix * m_Info.bytesPerPix);
	}
	else
	{
		mem = slice;
	}
	return mem;
}


template <typename pixelType> typename itk::Image<pixelType, 3>::Pointer Image::GetItkPtr(itk::SizeValueType T, itk::SizeValueType CH, PtrMode mode)
{
	if( !IsMatch<pixelType>(m_Info.dataType) ) //Forced to duplicate image because the datatypes are different
		mode = DEEP_COPY;

	if( T >= m_Info.numTSlices || CH >= m_Info.numChannels )
		return NULL;

	bool makeCopy;
	bool itkManageMemory;

	if(mode == RELEASE_CONTROL)
	{
		if(imageDataPtrs[T][CH].manager != FTK)	//I can't release control because I don't have it
			return NULL;

		makeCopy = false;
		itkManageMemory = true;
	}
	else if(mode == DEEP_COPY)
	{
		makeCopy = true;
		itkManageMemory = true;
	}
	else		//DEFAULT
	{
		makeCopy = false;
		itkManageMemory = false;
	}

	itk::SizeValueType numPixels = m_Info.numColumns * m_Info.numRows * m_Info.numZSlices;
	itk::SizeValueType numBytes = sizeof(pixelType) * numPixels;
	
	void * mem = NULL;	//This will point to the data I want to create the array from;
	if(makeCopy)		//Make a copy of the data	
	{
		//mem = (void *)new unsigned char[numBytes];
		mem = malloc(numBytes);
		if(mem == NULL)
			return NULL;
		if( !IsMatch<pixelType>(m_Info.dataType) ){		//Datatypes are not the same; use C-style cast and copy all pixels
			pixelType *out_pixel_container = static_cast<pixelType *>(mem);
			for(itk::SizeValueType k=0; k<m_Info.numZSlices; ++k)
				for(itk::SizeValueType j=0; j<m_Info.numRows; ++j)
					for(itk::SizeValueType i=0; i<m_Info.numColumns; ++i){
						//*out_pixel_container = static_cast<pixelType>( this->GetPixel(T,CH,k,j,i) );
						*out_pixel_container = this->GetPixelT<pixelType>(T,CH,k,j,i);
						++out_pixel_container;
					}
		}
		else											//Datatypes are the same; use compiler optimizations to memcpy 
			memcpy(mem,imageDataPtrs[T][CH].mem,numBytes);
	}
	else
	{
		mem = imageDataPtrs[T][CH].mem;
	}

	bool letItkManageMemory = false;			//itk DOES NOT manage the memory (default)
	if( itkManageMemory )	//DEEP COPY or RELEASE CONTROL
	{
		if(mode == RELEASE_CONTROL)
			imageDataPtrs[T][CH].manager = ITK;
		letItkManageMemory = true;	//itk DOES manage the memory
	}

	typedef itk::Image< pixelType, 3 > OutputImageType;
	typedef typename OutputImageType::PixelContainer ImageContainerType;
	typename ImageContainerType::Pointer container = ImageContainerType::New();

	container->Initialize();
	container->Reserve(numPixels);
	container->SetImportPointer( static_cast<pixelType *>(mem), numPixels, letItkManageMemory );

	
	typename OutputImageType::Pointer image = OutputImageType::New();

	typename OutputImageType::PointType origin;
   	origin[0] = 0; 
	origin[1] = 0;    
	origin[2] = 0; 
    image->SetOrigin( origin );
    typename OutputImageType::IndexType start;
    start[0] = 0;  // first index on X
    start[1] = 0;  // first index on Y    
	start[2] = 0;  // first index on Z    
    typename OutputImageType::SizeType  size;
	size[0] = m_Info.numColumns;  // size along X
    size[1] = m_Info.numRows;  // size along Y
	size[2] = m_Info.numZSlices;  // size along Z
    typename OutputImageType::RegionType region;
    region.SetSize( size );
    region.SetIndex( start );
    image->SetRegions( region );
    typename OutputImageType::SpacingType spacing;
    spacing[0] = m_Info.spacing.at(0);
    spacing[1] = m_Info.spacing.at(1);
    spacing[2] = m_Info.spacing.at(2);
    image->SetSpacing(spacing);
	image->Allocate();

	image->SetPixelContainer(container);
	image->Update();

	return image;
}

template <typename newType> void Image::Cast()
{
	//If I match do nothing
	if( IsMatch<newType>(m_Info.dataType) )
		return;
	
	itk::SizeValueType x = m_Info.numColumns;
	itk::SizeValueType y = m_Info.numRows;
	itk::SizeValueType z = m_Info.numZSlices;
	int n = m_Info.bytesPerPix;
	int numBytes = x*y*z*n;
	
	for (itk::SizeValueType t=0; t<m_Info.numTSlices; ++t)
	{
		for (int ch=0; ch<m_Info.numChannels; ++ch)
		{
			//Cast existing data to a char for access
			//char *p = static_cast<char *>(imageDataPtrs[t][ch].mem);	//any 8-bit type works here
			
			//Create a new array of the new type to put data into
			//newType *newArray = new newType[numBytes];
			newType *newArray = (newType *)malloc(numBytes);
			if(newArray == NULL)
				return;
			
			for(itk::SizeValueType k=0; k<z; ++k)
			{
				for(itk::SizeValueType j=0; j<y; ++j)
				{
					for(itk::SizeValueType i=0; i<x; ++i)
					{
						newType pix = this->GetPixelT<newType>(t,ch,k,j,i);
						newArray[k*y*x + j*x + i] = pix;
					}
				}
			}
			//delete[] imageDataPtrs[t][ch].mem;
			free( imageDataPtrs[t][ch].mem );
			imageDataPtrs[t][ch].mem = (void *)newArray;
			
		}	//end channels loop
	}	//end t-slices loop
	
	m_Info.bytesPerPix = sizeof(newType);		//Number of bytes per pixel (UCHAR - 1 or USHORT - 2)		
	m_Info.dataType = GetDataType<newType>();	//From enum ImgDataType;
	
}

//IF MULTIPLE T, SAVE THEM IN THEIR OWN FILE:
template <typename TPixel> bool Image::WriteImageITK(itk::SizeValueType channel, std::string baseName, std::string ext)
{
	if(m_Info.numTSlices == 1)
	{
		std::string fullname = baseName + "." + ext;
		if(!WriteImageITK<TPixel>(fullname, 0, channel))
			return false;
	}
	else
	{
		for(itk::SizeValueType t=0; t<m_Info.numTSlices; ++t)
		{
			std::string fullname = baseName + "_t" + itoa(t) + "." + ext;
			if(!WriteImageITK<TPixel>(fullname, t, channel))
				return false;
		}
	}
	return true;
}

template <typename TPixel> bool Image::WriteImageITK(std::string fileName, itk::SizeValueType T, itk::SizeValueType CH)
{
	//Will not work if wrong pixel type is requested
	if(!IsMatch<TPixel>(m_Info.dataType))
		return false;
		
	if( T >= m_Info.numTSlices || CH >= m_Info.numChannels )
		return false;
		
	typedef itk::ImageBase< 3 >                 ImageBaseType;
	typedef itk::Image< TPixel, 3 >             ImageType;
	typedef itk::ImageFileWriter< ImageType >   WriterType;
	typedef typename WriterType::Pointer        WriterPointer;
	
	WriterPointer writer = WriterType::New(); 
    writer->SetFileName( fileName );
    writer->SetInput( this->GetItkPtr<TPixel>(T, CH) );
    
    try
	{
		writer->Update();
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cerr << excp << std::endl;
		writer = 0;
		return false;
	}
	
	writer = 0;
	return true;
}

template <typename pixelType1> bool Image::IsMatch(DataType pixelType2)
{
	bool retVal = false;
	
	switch(pixelType2)
	{
		case itk::ImageIOBase::CHAR:
			if( typeid(char) == typeid(pixelType1) ) retVal = true;
		break;
		case itk::ImageIOBase::UCHAR:
			if( typeid(unsigned char) == typeid(pixelType1) ) retVal = true;
		break;
		case itk::ImageIOBase::SHORT:
			if( typeid(short) == typeid(pixelType1) ) retVal = true;
		break;
		case itk::ImageIOBase::USHORT:
			if( typeid(unsigned short) == typeid(pixelType1) ) retVal = true;
		break;
		case itk::ImageIOBase::INT:
			if( typeid(int) == typeid(pixelType1) ) retVal = true;
		break;
		case itk::ImageIOBase::UINT:
			if( typeid(unsigned int) == typeid(pixelType1) ) retVal = true;
		break;
		case itk::ImageIOBase::LONG:
			if( typeid(long) == typeid(pixelType1) ) retVal = true;
		break;
		case itk::ImageIOBase::ULONG:
			if( typeid(unsigned long) == typeid(pixelType1) ) retVal = true;
		break;
		case itk::ImageIOBase::FLOAT:
			if( typeid(float) == typeid(pixelType1) ) retVal = true;
		break;
		case itk::ImageIOBase::DOUBLE:
			if( typeid(double) == typeid(pixelType1) ) retVal = true;
		break;
    //just silencing a warning for now...
    case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
    break;
	}
	
	return retVal;	
}

template <typename pixelType> Image::DataType Image::GetDataType()
{
	Image::DataType retVal;

	if( typeid(char) == typeid(pixelType) ) retVal = itk::ImageIOBase::CHAR;
	else if( typeid(unsigned char) == typeid(pixelType) ) retVal = itk::ImageIOBase::UCHAR;
	else if( typeid(short) == typeid(pixelType) ) retVal = itk::ImageIOBase::SHORT;
	else if( typeid(unsigned short) == typeid(pixelType) ) retVal = itk::ImageIOBase::USHORT;
	else if( typeid(int) == typeid(pixelType) ) retVal = itk::ImageIOBase::INT;
	else if( typeid(unsigned int) == typeid(pixelType) ) retVal = itk::ImageIOBase::UINT;
	else if( typeid(long) == typeid(pixelType) ) retVal = itk::ImageIOBase::LONG;
	else if( typeid(unsigned long) == typeid(pixelType) ) retVal = itk::ImageIOBase::ULONG;
	else if( typeid(float) == typeid(pixelType) ) retVal = itk::ImageIOBase::FLOAT;
	else if( typeid(double) == typeid(pixelType) ) retVal = itk::ImageIOBase::DOUBLE;
	else retVal = itk::ImageIOBase::UNKNOWNCOMPONENTTYPE;
	
	return retVal;	
}

//pType is the pixel type to cast to, to retrieve the data from the pointer
//rType is the type to cast to, to return the data to user
template <typename pType, typename rType> rType Image::GetPixelValue(void * p)
{
	pType * pval = static_cast<pType *>(p);
	return static_cast<rType>(*pval);
}

//Casts the value to rType and returns it
template <typename rType> rType Image::GetPixelT(itk::SizeValueType T, itk::SizeValueType CH, itk::SizeValueType Z, itk::SizeValueType R, itk::SizeValueType C)
{
	if( T >= m_Info.numTSlices || CH >= m_Info.numChannels || Z >= m_Info.numZSlices \
		|| R >= m_Info.numRows || C >= m_Info.numColumns )
		return static_cast<rType>(0.0);

	if( T < 0 || CH < 0 || Z < 0 || R < 0 || C < 0 )
		return static_cast<rType>(0.0);

	itk::SizeValueType x = m_Info.numColumns;
	itk::SizeValueType y = m_Info.numRows;
	unsigned int n = m_Info.bytesPerPix;
	char *p = static_cast<char *>(imageDataPtrs[T][CH].mem) + Z*y*x*n + R*x*n + C*n;	//any 8-bit type works here

	double value = 0; 
	switch(m_Info.dataType)
	{
		case itk::ImageIOBase::CHAR:
			value = GetPixelValue<char, rType>(p);
		break;
		case itk::ImageIOBase::UCHAR:
			value = GetPixelValue<unsigned char, rType>(p);
		break;
		case itk::ImageIOBase::SHORT:
			value = GetPixelValue<short, rType>(p);
		break;
		case itk::ImageIOBase::USHORT:
			value = GetPixelValue<unsigned short, rType>(p);
		break;
		case itk::ImageIOBase::INT:
			value = GetPixelValue<int, rType>(p);
		break;
		case itk::ImageIOBase::UINT:
			value = GetPixelValue<unsigned int, rType>(p);
		break;
		case itk::ImageIOBase::LONG:
			value = GetPixelValue<long, rType>(p);
		break;
		case itk::ImageIOBase::ULONG:
			value = GetPixelValue<unsigned long, rType>(p);
		break;
		case itk::ImageIOBase::FLOAT:
			value = GetPixelValue<float, rType>(p);
		break;
		case itk::ImageIOBase::DOUBLE:
			value = GetPixelValue<double, rType>(p);
		break;
		default:
			//These are not supported!
		break;
	}
	return value;
}

template<typename TComp> void Image::LoadImageITK(std::string fileName, itk::SizeValueType numChannels, itkPixelType pixType, bool stacksAreForTime, bool appendChannels)
{
	if(imageDataPtrs.size() > 0)
	{
		if(m_Info.bytesPerPix != sizeof(TComp))
		{
			itk::ExceptionObject excp;
			excp.SetDescription("Component Types do not match");
			throw excp;
			return;
		}
	}

	switch( pixType)
	{
		//All of these types are based on itk::FixedArray (or can be) so I can load them.
		case itk::ImageIOBase::SCALAR:
		case itk::ImageIOBase::RGB:
		case itk::ImageIOBase::RGBA:
		case itk::ImageIOBase::VECTOR:
		case itk::ImageIOBase::POINT:
		case itk::ImageIOBase::COVARIANTVECTOR:
		case itk::ImageIOBase::SYMMETRICSECONDRANKTENSOR:
		case itk::ImageIOBase::DIFFUSIONTENSOR3D:
			LoadImageITK<TComp>( fileName, numChannels, stacksAreForTime, appendChannels );
		break;

		case itk::ImageIOBase::OFFSET:
		case itk::ImageIOBase::COMPLEX:
		case itk::ImageIOBase::MATRIX:
		default:
			itk::ExceptionObject excp;
			excp.SetDescription("Pixel type is not supported in this application");
			throw excp;
			return;
		break;
	}
}

//THIS IS AN INTERMEDIATE STEP WHEN A TYPE WITH VARIABLE LENGTH FIELDS IS USED:
template< typename TComp> void Image::LoadImageITK(std::string fileName, itk::SizeValueType numChannels, bool stacksAreForTime, bool appendChannels)
{
	switch(numChannels)
	{
	case 1:
		LoadImageITK< TComp, 1 >( fileName, stacksAreForTime, appendChannels );
	break;
	case 2:
		LoadImageITK< TComp, 2 >( fileName, stacksAreForTime, appendChannels );
	break;
	case 3:
		LoadImageITK< TComp, 3 >( fileName, stacksAreForTime, appendChannels );
	break;
	case 4:
		LoadImageITK< TComp, 4 >( fileName, stacksAreForTime, appendChannels );
	break;
	case 5:
		LoadImageITK< TComp, 5 >( fileName, stacksAreForTime, appendChannels );
	break;
	case 6:
		LoadImageITK< TComp, 6 >( fileName, stacksAreForTime, appendChannels );
	break;
		//NOTE 6 IS MAXIMUM ITK WILL HANDLE
	}
}

//THIS ASSUMES THAT THE IMAGE CAN BE LOADED AS A FIXED ARRAY!!!! 
//IF DATA ALREADY EXITS WILL ATTEMPT TO APPEND THE NEW IMAGE
template< typename TComp, itk::SizeValueType channels > void Image::LoadImageITK(std::string fileName, bool stacksAreForTime, bool appendChannels)
{
	typedef itk::FixedArray<TComp,channels>		PixelType;
	typedef itk::Image< PixelType, 5 >			ImageType;
	typedef typename ImageType::Pointer			ImagePointer;		
	typedef itk::ImageFileReader< ImageType >   ReaderType;
	typedef typename ReaderType::Pointer        ReaderPointer;

	ReaderPointer reader = ReaderType::New();
	reader->SetFileName( fileName );
	reader->Update();
	ImagePointer img = reader->GetOutput();

	//Set up the size info of the new image I am about to load:
	typename ImageType::RegionType region = img->GetBufferedRegion();
	typename ImageType::RegionType::IndexType start = region.GetIndex();
	typename ImageType::RegionType::SizeType size = region.GetSize();
	Info nInfo;											//Info of this new image!!
	nInfo.numColumns = size[0] - start[0];				//x-dimension
	nInfo.numRows = size[1] - start[1];					//y-dimension
	nInfo.numZSlices = size[2] - start[2];				//z-dimension
	nInfo.numTSlices = size[3] - start[3];		  	//t-dimension
	nInfo.bytesPerPix = sizeof(TComp);					//Already know bytes per pixel by component type
	nInfo.dataType = m_Info.dataType;					//Preserve (set earlier)
	nInfo.numChannels = channels;						//Part of template definition
	nInfo.spacing = m_Info.spacing;						//Preserve (set earlier)

	if(stacksAreForTime)	//Switch the Z and T numbers:
	{
    int tmp = nInfo.numTSlices;
		nInfo.numTSlices = nInfo.numZSlices;
		nInfo.numZSlices  = tmp;
	}

	unsigned short startingCH = 0;
	unsigned short startingT = 0;

	//If data is already loaded we must make sure that this new image will be appended correctly:
	if(imageDataPtrs.size() > 0)		//Already have data
	{
		//Check Size to be sure match:
		if(m_Info.BytesPerChunk() != nInfo.BytesPerChunk())
		{
			itk::ExceptionObject excp;
			excp.SetDescription("Image Sizes do not match");
			throw excp;
			return;
		}

		if(appendChannels)		//T must also match!!!
		{
			if(m_Info.numTSlices != nInfo.numTSlices)
			{
				itk::ExceptionObject excp;
				excp.SetDescription("Image T Slices does not match");
				throw excp;
				return;
			}
			startingCH = m_Info.numChannels;	//Number of channels changes
		}

		if(stacksAreForTime)	//Channels must also match!!!
		{
			if(m_Info.numChannels != nInfo.numChannels)		//Make sure number of components matches:
			{
				itk::ExceptionObject excp;
				excp.SetDescription("Number of Channels does not match");
				throw excp;
				return;
			}
			startingT = m_Info.numTSlices;		//Number of existing T slices changes
		}
	}
	else	//Don't already have data
	{
		m_Info = nInfo;		//Set the Image Info
	}
	
	//Change the number of Time slices if I'm adding some, and resize the vector:
	m_Info.numTSlices = startingT + nInfo.numTSlices;
	imageDataPtrs.resize(m_Info.numTSlices);

	m_Info.numChannels = startingCH + nInfo.numChannels;

	itk::SizeValueType numBytesPerChunk = m_Info.BytesPerChunk();

	//Create a pointer for each channel and time slice (allocate memory)
	ImageMemoryBlock block;
	block.manager = FTK;
	for(itk::SizeValueType t=startingT; t<m_Info.numTSlices; ++t)
	{
		for(itk::SizeValueType c=startingCH; c<m_Info.numChannels; ++c)
		{
			//block.mem = (void *)(new TComp[numBytesPerChunk]);
			void * mem = malloc(numBytesPerChunk);
			if(mem == NULL)
				return;
			block.mem = mem;
			imageDataPtrs[t].push_back(block);
		}
	}

	//Iterate through the input image and extract time and channel images:
	typedef itk::ImageRegionConstIterator< ImageType > IteratorType;
	IteratorType it( img, img->GetRequestedRegion() );
	itk::SizeValueType b = 0;
	itk::SizeValueType t = startingT;
	for ( it.GoToBegin(); !it.IsAtEnd(); ++it) 
	{
		//create a pixel object & Get each channel value
		typename ImageType::PixelType pixelValue = it.Get();
		unsigned short ch = startingCH;
		for(int c=0; c<nInfo.numChannels; ++c)
		{
			TComp *toLoc = ((TComp*)(imageDataPtrs[t][ch++].mem));
			toLoc[b] = (TComp)pixelValue[c];
		}

		//Update pointers
		b++;
		if(b>=numBytesPerChunk)	//I've finished this time slice
		{
			b=0;
			t++;
		}
	}
	img = 0;		//itk smartpointer cleans itself	
}

}  // end namespace ftk
#endif
