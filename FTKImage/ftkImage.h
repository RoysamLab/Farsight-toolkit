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
#include <itkImportImageContainer.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

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

	bool LoadFile( std::string fName, bool castToUchar = false );
	void LoadFiles( std::vector< std::string > fNames );
	bool SaveAs( std::string path, std::string fName, std::string ext );

	bool ImageFromData3D(void *dptr, ImageDataType dataType, int bpPix, int cs, int rs, int zs);
	void SetSpacing(int x, int y, int z);

	std::vector< unsigned short > Size(void);

	void * GetDataPtr(int T, int CH);

	//Also have templated functions
	template <typename rType> rType GetPixel(int T, int CH, int Z, int R, int C);		// Casts the value to rType and returns it
	template <typename pixelType> void SetPixel(int T, int Ch, int Z, int R, int C, pixelType newValue);// Casts from pixelType to image pixel type and sets pixel
	template <typename pixelType> pixelType * GetSlicePtr(int T, int CH, int Z);		// IF pixelType agrees with image pixel type
	template <typename pixelType> typename itk::Image<pixelType, 3>::Pointer GetItkPtr(int T, int CH); //IF pixelType agrees with image pixel type
	template <typename newType> void Cast();	//Cast the Image to newType (does not scale)

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
		ImageDataType dataType;			//From enum ImgDataType;

		std::vector< std::vector <unsigned char> > channelColors;	//Holds the color components of each channel
		std::vector< std::string > channelNames;					//Holds the name of each channel

		std::vector<int> spacing;		//Holds the spacing of the image (defaults to 1,1,1 (x,y,z) )

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

	bool LoadStandardImage( std::string filename, bool forDisplay = false );
	bool LoadLSMImage( std::string fileName );

	template<typename pType, typename rType> rType GetPixelValue(void * p);
	template <typename inType, typename outType> outType CastValue(inType inVal);
	template<typename pixelType1> bool IsMatch(ImageDataType pixelType2);
	template<typename pixelType> ImageDataType GetDataType();
	template<typename TPixel> bool Write(std::string fullFilename, int T, int CH);
	template<typename TPixel> bool WriteAll(std::string path, std::string baseName, std::string ext);

};

}  // end namespace ftk

#include "ftkImage.txx"

#endif	//end __ftkImage_h

