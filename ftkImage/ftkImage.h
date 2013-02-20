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

/*=========================================================================

  Program:   Farsight Biological Image Segmentation and Visualization Toolkit
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision: 0.00 $
  USAGE:	 SEE BOTTOM OF FILE

=========================================================================*/
#ifndef __ftkImage_h
#define __ftkImage_h

//ITK includes:
#include <itkImage.h>
#include <itkImageIOBase.h>
#include <itkLightObject.h>
#include <itkObjectFactory.h>
#include <itkSmartPointer.h>

//VTK includes:
#include "vtkImageData.h"
#include "vtkSmartPointer.h"
#include "vtkDataArray.h"

//Std includes:
#include <string>

namespace ftk
{

//**************************************************************************************************************
//This Image class can load a single image file as an image or multiple image files that should be associated 
//as a single image in memory.
//**************************************************************************************************************
class Image : public itk::LightObject
{
public:

	Image();
	~Image();

	/** Smart pointer typedef support. */
	typedef Image Self;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	/** Methods for creation through the object factory. */
	itkNewMacro(Self);								
	itkTypeMacro(Image,itk::LightObject);										

	typedef vtkSmartPointer<vtkImageData> VtkImagePtr;
	typedef itk::ImageIOBase::IOComponentType DataType;
	typedef itk::ImageIOBase::IOPixelType itkPixelType;
	typedef enum { FTK, VTK, ITK, OTHER } WhoManageMemory;
	typedef enum { DEFAULT, RELEASE_CONTROL, DEEP_COPY } PtrMode;

	//Each of these LoadFile commands will clear any previous image data
	bool LoadFile( std::string fName ); //Load 1 file (multi-page assumed to be z-direction)
	bool LoadFileAsTimeSeries( std::string fName ); //Load 1 file (multi-page assumed to be time series)
	bool LoadFileSeries( std::string arg, int start, int end, int step ); //Always assume each file contains a new Z
	bool LoadFilesAsMultipleChannels(std::vector<std::string> fnames, std::vector<std::string> channelnames, std::vector<unsigned char> colors);

	bool SaveChannelAs( int channel, std::string baseName, std::string ext );

	bool AppendChannelFromData3D(void *dptr, DataType dataType, int bpPix, itk::SizeValueType cs, itk::SizeValueType rs, itk::SizeValueType zs, std::string name, std::vector<unsigned char> color, bool copy);
	bool AppendImageFromData3D(void *dptr, DataType dataType, int bpPix, itk::SizeValueType cs, itk::SizeValueType rs, itk::SizeValueType zs, std::string name, bool copy);
	bool AppendImage(ftk::Image::Pointer img, PtrMode mode);	//Will add the image data as a new time slice or slices if all other sizes match.
	bool AppendImage(ftk::Image::Pointer img, PtrMode mode, bool isforOneTime);    // overloaded function	
	void SetSpacing(float x, float y, float z);

	std::vector< itk::SizeValueType > Size(void);
	std::vector< std::string > GetFilenames(void){ return filenames; };
	std::vector< std::vector <std::string> > GetTimeChannelFilenames(void){ return this->FileNames; };
	void SetTimeChannelFilenames(std::vector< std::vector <std::string> > filenames){this->FileNames.clear();this->FileNames = filenames;};

	void * GetDataPtr(itk::SizeValueType T, itk::SizeValueType CH, PtrMode mode = DEFAULT);			//Returns void * to this 3D stack using 1 of 3 modes, PtrMode defaults to DEFAULT
	VtkImagePtr GetVtkPtr(itk::SizeValueType T, itk::SizeValueType CH, PtrMode mode = DEFAULT);		//Returns vtkSmartPointer of vtkImageData at this T and CH, PtrMode defaults to DEFAULT
	void SetPixel(itk::SizeValueType T, itk::SizeValueType Ch, itk::SizeValueType Z, itk::SizeValueType R, itk::SizeValueType C, double newValue); // Casts from double to image pixel type and sets pixel
	double GetPixel(itk::SizeValueType T, itk::SizeValueType CH, itk::SizeValueType Z, itk::SizeValueType R, itk::SizeValueType C);				// Casts the value to double and returns it
	std::vector< std::string > GetChannelNames(void){ return m_Info.channelNames; };

	//Also have templated functions
	template <typename rType> rType GetPixelT(itk::SizeValueType T, itk::SizeValueType CH, itk::SizeValueType Z, itk::SizeValueType R, itk::SizeValueType C);	//Casts the value to rType and returns it
	template <typename newType> void Cast();	//Cast the Image to newType (does not scale)
	template <typename pixelType> typename itk::Image<pixelType, 3>::Pointer GetItkPtr(itk::SizeValueType T, itk::SizeValueType CH, PtrMode mode = DEFAULT);	//IF pixelType agrees with image pixel type, PtrMode defaults to DEFAULT
	template <typename pixelType> pixelType * GetSlicePtr(itk::SizeValueType T, itk::SizeValueType CH, itk::SizeValueType Z,PtrMode mode = DEFAULT);	// IF pixelType agrees with image pixel type (NOTE MEMORY MANAGER DOES NOT CHANGE)

	typedef struct 
	{
		itk::SizeValueType numColumns;		//Number of Columns in Image (x)
		itk::SizeValueType numRows;		//Number of Rows in Image (y)
		itk::SizeValueType numZSlices;		//Number of Z Slices in Image (z)
		itk::SizeValueType numTSlices;		//Number of Time Slices in Image (t)
		itk::SizeValueType numChannels;		//Number of Channels in this Image (ch)
		int bytesPerPix;		//Number of bytes per pixel (UCHAR - 1 or USHORT - 2)		
		DataType dataType;				//From enum ftk::Image::DataType;

		std::vector< std::vector <unsigned char> > channelColors;	//Holds the color components of each channel
		std::vector< std::string > channelNames;					//Holds the name of each channel

		std::vector<float> spacing;		//Holds the spacing of the image (defaults to 1,1,1 (x,y,z) )

		itk::SizeValueType BytesPerChunk(void)
		{
			return numZSlices*numRows*numColumns*bytesPerPix;
		};

	} Info;

	const Info * GetImageInfo(void) { return &(this->m_Info); };	//Returns a pointer to constant data (these values cannot be changed from outside

protected:


private:
	Image(const Self&);				//purposely not implemented
	void operator=(const Self&);	//purposely not implemented

	typedef struct
	{	void * mem;						//A Memory Block Ptr
		WhoManageMemory manager;		//Who Manages this memory block?
	} ImageMemoryBlock;

	//Private Variables:
	Info m_Info;									//Hold Image info of the 'image file' that make up this 'image'
	std::string path;								//Path to this image file
	std::vector< std::string > filenames;			//Filenames of this image
	std::vector< std::vector< ImageMemoryBlock > > imageDataPtrs;		//Pointers to all of the data
	std::vector< std::vector< std::string > > FileNames; // Filenames of this 5D image (including both time and channel names)

	//Private Functions:
	void DeleteData();
	std::string GetFileExtension(std::string);
	std::string GetFilename(std::string);
	std::string GetPath(std::string);
	std::string itoa(const int x);

	bool LoadStandardImage( std::string fileName, bool stacksAreForTime = false, bool appendChannels = false, bool readRGBasSingleChannel = true );
	void SetDefaultColors(void);
#if VTK_MAJOR_VERSION <= 5
	bool LoadLSMImage( std::string fileName );
#endif

	vtkSmartPointer<vtkDataArray> GetVtkDataArray(itk::SizeValueType T, itk::SizeValueType CH, bool makeCopy, bool vtkManageMemory); //Returns vtkSmartPointer at this T and CH
	int GetDataTypeVTK(DataType itk);
	DataType GetDataTypeITK(int vtk_type);

	template<typename pType, typename rType> rType GetPixelValue(void * p);
	template<typename pixelType1> bool IsMatch(DataType pixelType2);
	template<typename pixelType> DataType GetDataType();

	template<typename TPixel> bool WriteImageITK(itk::SizeValueType channel, std::string baseName, std::string ext);
	template<typename TPixel> bool WriteImageITK(std::string fullFilename, itk::SizeValueType T, itk::SizeValueType CH);

	template<typename TComp> void LoadImageITK(std::string fileName, itk::SizeValueType numChannels, itkPixelType pixType, bool stacksAreForTime, bool appendChannels);
	template<typename TComp> void LoadImageITK(std::string filename, itk::SizeValueType numChannels, bool stacksAreForTime, bool appendChannels);
	template<typename TComp, itk::SizeValueType channels> void LoadImageITK(std::string fileName, bool stacksAreForTime, bool appendChannels);

};

}  // end namespace ftk

#include "ftkImage.txx"

#endif	//end __ftkImage_h

/**********************************************************************************************************************
HOW TO USE FTK::IMAGE

-Intro

This class has been created to handle five dimensional image loading/creation/manipulation/writing.
The dimensions included are 3 spatial, 1 time, and 1 channel.  All images are stored in memory using
void* to a 3D image block, this means that each time and channel component of the image may have its
own memory location and are stored separately.  This is done because it is most common for segmentation
algorithms to be executed on one channel at a time - therefore it would be required that the channels be
split at some point anyway.

-Loading From File(s)

Three methods for loading and image from file have been created.
LoadFile() will load one file and assume that multiple pages make up a 3D spatial dimension.
LoadFileAsTimeSeries() will load one file and assume that multiple pages make up the time dimension.
LoadFileSeries() will create an image from multiple files following the filename pattern provided.
This method assumes each file contains all spatial dimensions and separate files make up the time dimension

These three methods for loading from file use one of two imageIO techniques.  The first is a special class for 
loading Ziess images (.lsm), the second is the standard ITK file reader.  ftk::Image should be able to handle all
image types the itk can load.

-Dynamically Creating Images

Images may also be created from existing data buffers using AppendChannelFromData3D().  This method will allow
you to create a 4D ftk image (no time dimension) by combining 3D channels.  To make 5D image you must create
a 4D ftk image for each t and then use AppendImage() to combine them.  This method will add the image data 
as a new time slice or slices if all other sizes match.

-Saving To File

Currently one method exists to save images to file: SaveChannelAs().  You must specify the time and channel you would
like to save.

-Accessing Image Data

One of the main strengths of ftk image is the ease at which the data can be accessed and manipulated.
It is important to remember that this flexibility also places a greater responsibility on you to make 
sure that you are using the correct memory modes (see below).

Pointers to the 3D images can be requested in ITK, VTK, or raw pointer form. When requesting the pointer
you may decide how the memory buffer for that image should be managed - no change, new pointer manages, or
a copy is made (see below).

Note that the request for an ITK image is templated, and will return a NULL pointer if the wrong pixel type
is requested.  One way to be assured that you will get the desired pixel type is to cast the ftk image to
that type by using Cast().  Note that this function could result in a loss of data.

FUTURE IMPROVEMENT: if an itk pointer is requested using DEEP_COPY the new image will be cast to the template
type and be returned appropriately.

Convenience functions GetPixel() and SetPixel() have been provided to allow for pixel access.
For speed, it is recommended that one uses itk iterators for changing the data by:
1. GetItkPtr() in DEFAULT memory mode
2. Use itk iterators to manipulate the data through the itk pointer

GetSlicePtr() is also provided to get a pointer to the beginning of a 2D slice. This method only works in DEFAULT memory
management mode.

-Memory Management

ftk::Image takes advantage of ITK's smart pointers for automatic cleanup.  When no references to the image are left,
the image will automatically delete itself and free the memory of the image data.  This is great, unless you have 
requested a pointer to one of the 3D chunks of data and are still using that data.  To solve this problem, 3 pointer
modes have been created:

1. DEFAULT
2. RELEASE_CONTROL
3. DEEP_COPY

The DEFAULT mode will pass a pointer to the data, but ftk image continues to manage the memory (i.e. cleanup the memory).
RELEASE_CONTROL will pass a pointer to the data, and allow this new pointer to clean up the data when it is deleted.
DEEP_COPY will create a new copy of the data, and both ftk image and the new pointer will have control of their own memory.

Once ftk image releases control of a memory location it cannot get it back and it cannot release it again.  This means
that if you have a created an itk image using RELEASE_CONTROL that itk image now controls the memory.  If you want another
ITK image you should use the pointer you have already created, or you can still get a DEEP_COPY.  This same principle
applies to VTK image data.


*/

