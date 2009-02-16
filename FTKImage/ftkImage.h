/*=========================================================================

  Program:   Farsight Biological Image Segmentation and Visualization Toolkit
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision: 0.00 $

=========================================================================*/
#ifndef __ftkImage_h
#define __ftkImage_h

//ITK includes:
#include <itkImage.h>
#include <itkRGBPixel.h>
#include <itkSmartPointer.h>
#include <itkImageRegionConstIterator.h>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

//VTK includes:
#include "vtkImageData.h"

//Local includes:
#include "vtkKWImage.h"
#include "vtkKWImageIO.h"
#include "vtkLSMReader.h"

//Std includes:
#include <string>


namespace ftk
{

//THESE ARE THE 4 TYPES OF ITK IMAGES THAT FTK KNOWS ABOUT GLOBALY (WHEN INCLUDING FTK IMAGE)
//YOU, OF COURSE, MAY USE OTHERS IN YOUR CODE
typedef itk::Image< unsigned char, 3 > UcharImage3DType;
typedef itk::Image< unsigned short, 3 > UshortImage3DType;

typedef itk::Image< unsigned char, 2 > UcharImage2DType;
typedef itk::Image< unsigned short, 2 > UshortImage2DType;

typedef itk::ImageBase< 3 > ImageBaseType;
typedef ImageBaseType::ConstPointer ImageBaseConstPtr;
typedef ImageBaseType::Pointer ImageBasePtr;

typedef  enum {IVOID,BIT,CHAR,UCHAR,SHORT,USHORT,INT,UINT,LONG,ULONG,FLOAT,DOUBLE} ImageDataType;

//**************************************************************************************************************
//This Image class can load a single image file as an image or multiple image files that should be associated 
//as a single image in memory.
//**************************************************************************************************************
class Image
{
public:

	Image();
	~Image();

	bool LoadFile( std::string fName );
	void LoadFiles( std::vector< std::string > fNames );

	std::vector< unsigned short > Size(void);

	unsigned char* GetSlicePtr(int T, int CH, int Z);	//IF SCALAR TYPE = UCHAR, RETURN POINTER
	void * GetDataPtr(int T, int CH);

	//Also have a templated function (must be included here)
	template <typename rType> rType GetPixel(int T, int CH, int Z, int R, int C);	//converts the value to rType and returns it

	typedef struct 
	{
		std::string path;				//Path to this image file
		std::string filename;			//Filename of this image
		unsigned short numColumns;		//Number of Columns in Image (x)
		unsigned short numRows;			//Number of Rows in Image (y)
		unsigned short numZSlices;		//Number of Z Slices in Image (z)
		unsigned short numTSlices;		//Number of Time Slices in Image (t)
		unsigned short numChannels;		//Number of Channels in this Image (ch)
		unsigned char bytesPerPix;		//Number of bytes per pixel (UCHAR - 1 or USHORT - 2)
		unsigned char dataType;			//From enum ImgDataType;

		std::vector< std::vector <unsigned char> > channelColors;	//Holds the color components of each channel
		std::vector< std::string > channelNames;					//Holds the name of each channel

	} Info;

	Info * GetImageInfo(void) { return &(this->imageInfo); };

private:
	//Private Variables:
	Info imageInfo;											//Hold Image info of the 'image file' that make up this 'image'
	std::vector< std::vector< void * > > imageDataPtrs;

	//Private Functions:
	std::string GetFileExtension(std::string);
	std::string GetFilename(std::string);
	std::string GetPath(std::string);
	std::string itoa(const int x);

	bool LoadStandardImage( std::string filename );
	bool LoadLSMImage( std::string fileName );

};

}  // end namespace ftk

#include "ftkImage.txx"

#endif	//end __ftkImage_h

