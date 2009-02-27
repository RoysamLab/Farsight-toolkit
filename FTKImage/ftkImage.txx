#ifndef _ftkImage_txx
#define _ftkImage_txx
#include "ftkImage.h"

namespace ftk
{

template <typename pixelType> pixelType * Image::GetSlicePtr(int T, int CH, int Z)
{
	if( T >= imageInfo.numTSlices || CH >= imageInfo.numChannels || Z >= imageInfo.numZSlices )
		return NULL;
		
	if( !IsMatch<pixelType>(imageInfo.dataType) )
		return NULL;

	pixelType * stack = static_cast<pixelType *>(imageDataPtrs[T][CH]);
	pixelType * slice = stack + Z*(imageInfo.numColumns)*(imageInfo.numRows);
	return ( slice );
}

template <typename pixelType> typename itk::Image<pixelType, 3>::Pointer Image::GetItkPtr(int T, int CH)
{
	if( !IsMatch<pixelType>(imageInfo.dataType) )
		return NULL;
		
	if( T >= imageInfo.numTSlices || CH >= imageInfo.numChannels )
		return NULL;

	typedef itk::ImportImageContainer<unsigned long, pixelType> ImageContainerType;
	ImageContainerType::Pointer container = ImageContainerType::New();
	
	int numPixels = imageInfo.numColumns * imageInfo.numRows * imageInfo.numZSlices;
	
	container->Initialize();
	container->Reserve(numPixels);
	container->SetImportPointer( static_cast<pixelType *>(imageDataPtrs[T][CH]), numPixels, false );
	
	typedef itk::Image< pixelType, 3 > OutputImageType;
	OutputImageType::Pointer image = OutputImageType::New();
	
	OutputImageType::PointType origin;
   	origin[0] = 0; 
    origin[1] = 0;    
	origin[2] = 0;    
    image->SetOrigin( origin );
    OutputImageType::IndexType start;
    start[0] = 0;  // first index on X
    start[1] = 0;  // first index on Y    
	start[2] = 0;  // first index on Z    
    OutputImageType::SizeType  size;
	size[0] = imageInfo.numColumns;  // size along X
    size[1] = imageInfo.numRows;  // size along Y
	size[2] = imageInfo.numZSlices;  // size along Z
    OutputImageType::RegionType region;
    region.SetSize( size );
    region.SetIndex( start );
    image->SetRegions( region );
    OutputImageType::SpacingType spacing;
    spacing[0] = imageInfo.spacing.at(0);
    spacing[1] = imageInfo.spacing.at(1);
    spacing[2] = imageInfo.spacing.at(2);
    image->SetSpacing(spacing);
	image->Allocate();
	
	image->SetPixelContainer(container);
	image->Update();
	
	return image;
}

template <typename newType> void Image::Cast()
{
	//If I match do nothing
	if( IsMatch<newType>(imageInfo.dataType) )
		return;
	
	int x = imageInfo.numColumns;
	int y = imageInfo.numRows;
	int z = imageInfo.numZSlices;
	int n = imageInfo.bytesPerPix;
	int numBytes = x*y*z*n;
	
	for (int t=0; t<imageInfo.numTSlices; ++t)
	{
		for (int ch=0; ch<imageInfo.numChannels; ++ch)
		{
			//Cast existing data to a char for access
			char *p = static_cast<char *>(imageDataPtrs[t][ch]);		//any 8-bit type works here
			
			//Create a new array of the new type to put data into
			newType *newArray = new newType[numBytes];
			
			for(int k=0; k<z; ++k)
			{
				for(int j=0; j<y; ++j)
				{
					for(int i=0; i<x; ++i)
					{
						newType pix = GetPixel<newType>(t,ch,k,j,i);
						newArray[k*y*x + j*x + i] = pix;
					}
				}
			}
			delete[] imageDataPtrs[t][ch];
			imageDataPtrs[t][ch] = (void *)newArray;
			
		}	//end channels loop
	}	//end t-slices loop
	
	imageInfo.bytesPerPix = sizeof(newType);		//Number of bytes per pixel (UCHAR - 1 or USHORT - 2)		
	imageInfo.dataType = GetDataType<newType>();	//From enum ImgDataType;
	
}

template <typename TPixel> bool Image::Write(std::string fullFilename, int T, int CH)
{
	//Will not work if wrong pixel type is requested
	if(!IsMatch<TPixel>(imageInfo.dataType))
		return false;
		
	if( T >= imageInfo.numTSlices || CH >= imageInfo.numChannels )
		return false;
		
	typedef itk::ImageBase< 3 >                 ImageBaseType;
	typedef itk::Image< TPixel, 3 >             ImageType;
	typedef itk::ImageFileWriter< ImageType >   WriterType;
	typedef typename WriterType::Pointer        WriterPointer;
	
	WriterPointer writer = WriterType::New(); 
    writer->SetFileName( fullFilename );
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

template <typename TPixel> bool Image::WriteAll(std::string path, std::string baseName, std::string ext)
{
	if(!IsMatch<TPixel>(imageInfo.dataType))
		return false;
	
	for (int t=0; t<imageInfo.numTSlices; ++t)
	{
		for (int ch=0; ch<imageInfo.numChannels; ++ch)
		{
			std::string fullName = path + "/" + baseName;
			if(imageInfo.numTSlices > 1)
				fullName += "_T" + itoa(t);
			if(imageInfo.numChannels > 1)
				fullName += "_CH" + itoa(ch);
				
			fullName += ("." + ext);
				
			if(!this->Write<TPixel>(fullName, t, ch))
			{
				std::cerr << "FAILED WRITE ON T=" << t << " & CH=" << ch << std::endl;
				return false;
			}
			
		}//end channels loop
	}//end t-slices loop
	
	return true;
}

template <typename pixelType1> bool Image::IsMatch(ImageDataType pixelType2)
{
	bool retVal = false;
	
	switch(pixelType2)
	{
		case CHAR:
			if( typeid(char) == typeid(pixelType1) ) retVal = true;
		break;
		case UCHAR:
			if( typeid(unsigned char) == typeid(pixelType1) ) retVal = true;
		break;
		case SHORT:
			if( typeid(short) == typeid(pixelType1) ) retVal = true;
		break;
		case USHORT:
			if( typeid(unsigned short) == typeid(pixelType1) ) retVal = true;
		break;
		case INT:
			if( typeid(int) == typeid(pixelType1) ) retVal = true;
		break;
		case UINT:
			if( typeid(unsigned int) == typeid(pixelType1) ) retVal = true;
		break;
		case LONG:
			if( typeid(long) == typeid(pixelType1) ) retVal = true;
		break;
		case ULONG:
			if( typeid(unsigned long) == typeid(pixelType1) ) retVal = true;
		break;
		case FLOAT:
			if( typeid(float) == typeid(pixelType1) ) retVal = true;
		break;
		case DOUBLE:
			if( typeid(double) == typeid(pixelType1) ) retVal = true;
		break;
	}
	
	return retVal;	
}

template <typename pixelType> ImageDataType Image::GetDataType()
{
	ImageDataType retVal;

	if( typeid(char) == typeid(pixelType) ) retVal = CHAR;
	else if( typeid(unsigned char) == typeid(pixelType) ) retVal = UCHAR;
	else if( typeid(short) == typeid(pixelType) ) retVal = SHORT;
	else if( typeid(unsigned short) == typeid(pixelType) ) retVal = USHORT;
	else if( typeid(int) == typeid(pixelType) ) retVal = INT;
	else if( typeid(unsigned int) == typeid(pixelType) ) retVal = UINT;
	else if( typeid(long) == typeid(pixelType) ) retVal = LONG;
	else if( typeid(unsigned long) == typeid(pixelType) ) retVal = ULONG;
	else if( typeid(float) == typeid(pixelType) ) retVal = FLOAT;
	else if( typeid(double) == typeid(pixelType) ) retVal = DOUBLE;
	else retVal = IVOID;
	
	return retVal;	
}

//pType is the pixel type to cast to, to retrieve the data from the pointer
//rType is the type to cast to, to return the data to user
template <typename pType, typename rType> rType Image::GetPixelValue(void * p)
{
	pType * pval = static_cast<pType *>(p);
	return static_cast<rType>(*pval);
}

template <typename inType, typename outType> outType Image::CastValue(inType inVal)
{
	return static_cast<outType>(inVal);
}

//*****************************************************************************************
// GetPixel
//*****************************************************************************************
template <typename rType> rType Image::GetPixel(int T, int CH, int Z, int R, int C)
{
	if( T >= imageInfo.numTSlices || CH >= imageInfo.numChannels || Z >= imageInfo.numZSlices \
		|| R >= imageInfo.numRows || C >= imageInfo.numColumns )
		return NULL;

	unsigned int x = imageInfo.numColumns;
	unsigned int y = imageInfo.numRows;
	unsigned int n = imageInfo.bytesPerPix;
	char *p = static_cast<char *>(imageDataPtrs[T][CH]) + Z*y*x*n + R*x*n + C*n;	//any 8-bit type works here

	rType value = 0; 
	switch(imageInfo.dataType)
	{
		case CHAR:
			value = GetPixelValue<char, rType>(p);
		break;
		case UCHAR:
			value = GetPixelValue<unsigned char, rType>(p);
		break;
		case SHORT:
			value = GetPixelValue<short, rType>(p);
		break;
		case USHORT:
			value = GetPixelValue<unsigned short, rType>(p);
		break;
		case INT:
			value = GetPixelValue<int, rType>(p);
		break;
		case UINT:
			value = GetPixelValue<unsigned int, rType>(p);
		break;
		case LONG:
			value = GetPixelValue<long, rType>(p);
		break;
		case ULONG:
			value = GetPixelValue<unsigned long, rType>(p);
		break;
		case FLOAT:
			value = GetPixelValue<float, rType>(p);
		break;
		case DOUBLE:
			value = GetPixelValue<double, rType>(p);
		break;
		case IVOID:
		case BIT:
		default:
			//These are not supported!
		break;
	}
	return value;
}

template <typename pixelType> void Image::SetPixel(int T, int CH, int Z, int R, int C, pixelType newValue)
{
	if( T >= imageInfo.numTSlices || CH >= imageInfo.numChannels || Z >= imageInfo.numZSlices \
		|| R >= imageInfo.numRows || C >= imageInfo.numColumns )
		return;
		
	unsigned int x = imageInfo.numColumns;
	unsigned int y = imageInfo.numRows;

	if (imageInfo.dataType == CHAR)
	{
			char value = CastValue<pixelType, char>(newValue);
			char * p = static_cast<char *>(imageDataPtrs[T][CH]) + Z*y*x + R*x + C;
			(*p) = value;
	}
	else if (imageInfo.dataType == UCHAR)
	{
			unsigned char value = CastValue<pixelType, unsigned char>(newValue);
			unsigned char * p = static_cast<unsigned char *>(imageDataPtrs[T][CH]) + Z*y*x + R*x + C;
			(*p) = value;
	}
	else if (imageInfo.dataType == SHORT)
	{
			short value = CastValue<pixelType, short>(newValue);
			short * p = static_cast<short *>(imageDataPtrs[T][CH]) + Z*y*x + R*x + C;
			(*p) = value;
	}
	else if (imageInfo.dataType == USHORT)
	{
			unsigned short value = CastValue<pixelType, unsigned short>(newValue);
			unsigned short * p = static_cast<unsigned short *>(imageDataPtrs[T][CH]) + Z*y*x + R*x + C;
			(*p) = value;
	}
	else if (imageInfo.dataType == INT)
	{
			int value = CastValue<pixelType, int>(newValue);
			int * p = static_cast<int *>(imageDataPtrs[T][CH]) + Z*y*x + R*x + C;
			(*p) = value;
	}
	else if (imageInfo.dataType == UINT)
	{
			unsigned int value = CastValue<pixelType, unsigned int>(newValue);
			unsigned int * p = static_cast<unsigned int *>(imageDataPtrs[T][CH]) + Z*y*x + R*x + C;
			(*p) = value;
	}
	else if (imageInfo.dataType == LONG)
	{
			long value = CastValue<pixelType, long>(newValue);
			long * p = static_cast<long *>(imageDataPtrs[T][CH]) + Z*y*x + R*x + C;
			(*p) = value;
	}
	else if (imageInfo.dataType == ULONG)
	{
			unsigned long value = CastValue<pixelType, unsigned long>(newValue);
			unsigned long * p = static_cast<unsigned long *>(imageDataPtrs[T][CH]) + Z*y*x + R*x + C;
			(*p) = value;
	}
	else if (imageInfo.dataType == FLOAT)
	{
			float value = CastValue<pixelType, float>(newValue);
			float * p = static_cast<float *>(imageDataPtrs[T][CH]) + Z*y*x + R*x + C;
			(*p) = value;
	}
	else if (imageInfo.dataType == DOUBLE)
	{
			double value = CastValue<pixelType, double>(newValue);
			double * p = static_cast<double *>(imageDataPtrs[T][CH]) + Z*y*x + R*x + C;
			(*p) = value;
	}
}

}  // end namespace ftk
#endif
