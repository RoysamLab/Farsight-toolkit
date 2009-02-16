#ifndef _ftkImage_txx
#define _ftkImage_txx
#include "ftkImage.h"

namespace ftk
{

//pType is the pixel type to cast to, to retrieve the data from the pointer
//rType is the type to cast to, to return the data to user
template <typename pType, typename rType> rType GetPixelValue(void * p)
{
	pType * pval = static_cast<pType *>(p);
	return static_cast<rType>(*pval);
}

//*****************************************************************************************
// GetPixel
//*****************************************************************************************
template <typename rType> rType Image::GetPixel(int T, int CH, int Z, int R, int C)
{
	if( T > imageInfo.numTSlices || CH > imageInfo.numChannels || Z > imageInfo.numZSlices \
		|| R > imageInfo.numRows || C > imageInfo.numColumns )
		return NULL;

	unsigned short x = imageInfo.numColumns;
	unsigned short y = imageInfo.numRows;
	unsigned short n = imageInfo.bytesPerPix;
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

}  // end namespace ftk
#endif
