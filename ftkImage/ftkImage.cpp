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

#include "ftkImage.h"

//itk includes:
#include <itkImageRegionConstIterator.h>
#include <itkNumericSeriesFileNames.h>
#include <itkRescaleIntensityImageFilter.h>

//vtk includes:
#include <vtkCharArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkShortArray.h>
#include <vtkUnsignedShortArray.h>
#include <vtkIntArray.h>
#include <vtkUnsignedIntArray.h>
#include <vtkLongArray.h>
#include <vtkUnsignedLongArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include "vtkPointData.h"

//Local includes:
#if VTK_MAJOR_VERSION <= 5
#include "vtkLSMReader.h"
#endif

//****************************************************************************************
// FTKImage is a class the handles the storage of multi-dimensional images.
// 
// The class is not complete, and should be updated and modified as needs arise.
// 
// images are stored in memory using void * for each channel and time-slice.
// 
// USAGE: see ftkImage.h
//****************************************************************************************

namespace ftk
{

//Constructor
Image::Image()
{
	path.clear();
	filenames.clear();
	imageDataPtrs.clear();
}

//Destructor
Image::~Image()
{
	DeleteData();
}

void Image::DeleteData()
{
	//Delete Image Memory:
	for(int i=0; i<(int)imageDataPtrs.size(); ++i)
	{
		for(int j=0; j<(int)imageDataPtrs[i].size(); ++j)
		{
			if(imageDataPtrs[i][j].manager == FTK)
			{
				//void *mem = imageDataPtrs[i][j].mem;
				//delete[] mem;
				free( imageDataPtrs[i][j].mem );
			}
		}
	}
	imageDataPtrs.clear();
}

void Image::SetSpacing(float x, float y, float z)
{
	m_Info.spacing.at(0) = x;
	m_Info.spacing.at(1) = y;
	m_Info.spacing.at(2) = z;
}

//Each file contains 1 or more channels that will be appended.  All other sizes must match in each file.
bool Image::LoadFilesAsMultipleChannels(std::vector<std::string> fnames, std::vector<std::string> channelnames, std::vector<unsigned char> colors)
{
	int count = (int)fnames.size();
	if( (int)channelnames.size() != count || (int)colors.size() != 3*count )
		return false;

	DeleteData();
	m_Info.channelColors.clear();
	m_Info.channelNames.clear();

	int counter = 0;
	for(int i=0; i<count; ++i)
	{

		if( colors.at(i*3+0) == 255 && colors.at(i*3+1) == 255 && colors.at(i*3+2) == 255 ){
			if( !this->LoadStandardImage( fnames.at(i), false, true, false ) ) //All input colors needed, load as RGB
				return false;
		}
		else if( !this->LoadStandardImage( fnames.at(i), false, true, true ) ) //Not All colors are needed, read as a single channel grayscale and scale to desired colorspace
			return false;

		std::string chname = channelnames.at(i);

		if( colors.at(i*3+0) == 255 && colors.at(i*3+1) == 255 && colors.at(i*3+2) == 255 && m_Info.numChannels == (counter+3) ){
			//If loaded as RGB and display with full colorscale
			std::vector< unsigned char > color(3);
			color[0] = (unsigned char)( colors.at(i*3 + 0) ); color[1] = 0;  color[2] = 0; m_Info.channelColors.push_back(color);
			color[0] = 0; color[1] = (unsigned char)( colors.at(i*3 + 1) );  color[2] = 0; m_Info.channelColors.push_back(color);
			color[0] = 0; color[1] = 0;  color[2] = (unsigned char)( colors.at(i*3 + 2) ); m_Info.channelColors.push_back(color);

			m_Info.channelNames.push_back(chname); m_Info.channelNames.push_back(chname); m_Info.channelNames.push_back(chname);

			counter += 3;
		} else {
			//Parse the color info for the scalar channel
			std::vector< unsigned char > color(3);
			color[0] = (unsigned char)( colors.at(i*3 + 0) );
			color[1] = (unsigned char)( colors.at(i*3 + 1) );
			color[2] = (unsigned char)( colors.at(i*3 + 2) );

			m_Info.channelColors.push_back(color);	//Use provided color
			m_Info.channelNames.push_back(chname);	//Use provided name

			++counter;
		}
	}


	this->filenames = fnames;
	//std::vector<std::string> fnamestmp;
	//for(int i=0; i<count; ++i)
	//	fnamestmp.push_back(this->GetFilename(fnames.at(i)));
	//this->FileNames.push_back(fnames);

	return true;
}

//Each file contains a 2D image (may be multiple channels) that will be a new Z:
//Always assume each file contains a new Z
bool Image::LoadFileSeries( std::string arg, int start, int end, int step)
{	
	DeleteData();

	typedef itk::NumericSeriesFileNames NameGeneratorType;
	NameGeneratorType::Pointer nameGenerator = NameGeneratorType::New();
	nameGenerator->SetSeriesFormat( arg );
	nameGenerator->SetStartIndex( start );
	nameGenerator->SetEndIndex( end );
	nameGenerator->SetIncrementIndex( step );
	
	std::vector<std::string> names = nameGenerator->GetFileNames();
	for(int i=0; i<(int)names.size(); ++i)
	{
		if( !this->LoadStandardImage( names.at(i), false, false, false ) )
			return false;
	}
	this->SetDefaultColors();
	return true;
}

////Load 1 file normally (2D/3D - multi-page assumed to be z-direction)
bool Image::LoadFile( std::string fName)
{
	DeleteData();

#if VTK_MAJOR_VERSION <= 5
	if( GetFileExtension(fName) == "lsm" )
		return this->LoadLSMImage( fName );
	else
#endif
	{
		if(!this->LoadStandardImage( fName, false, false, false))
			return false;
		this->SetDefaultColors();
	}
	return true;
}

//Attempts to load a multi-page image and treat each page as a separate time slice.
//Does not work for LSM images
bool Image::LoadFileAsTimeSeries( std::string fName)
{
	DeleteData();

	if( GetFileExtension(fName) == "lsm" )
		return false;
	else
	{
		if(!this->LoadStandardImage( fName, true, false, false))
			return false;
		this->SetDefaultColors();
	}
	return true;
}

//**********************************************************************************************************
// IF IMAGE HAS MULTIPLE PAGES IT IS ASSUMED THAT THEY REPRESENT A 3D IMAGE, BUT THEY CAN BE FORCED TO BE T
//
//***********************************************************************************************************
bool Image::LoadStandardImage( std::string fileName, bool stacksAreForTime, bool appendChannels, bool readRGBasSingleChannel )
{
	// Find out the pixel type of the image in file
	itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO( fileName.c_str(), itk::ImageIOFactory::ReadMode );
	if( !imageIO )
	{
		std::cerr << "NO IMAGEIO WAS FOUND" << std::endl;
		return false;
	}
	// Now that we found the appropriate ImageIO class, ask it to 
	// read the meta data from the image file.
	// Now that we found the appropriate ImageIO class, ask it to 
	// read the meta data from the image file.
	imageIO->SetFileName( fileName.c_str() );
	imageIO->ReadImageInformation();

	int numComponents;
	itkPixelType pixelType = imageIO->GetPixelType();
	if( readRGBasSingleChannel && ( pixelType == itk::ImageIOBase::RGB || pixelType == itk::ImageIOBase::RGBA ) ){
		pixelType = itk::ImageIOBase::SCALAR;
		numComponents = 1;
	}
	else
		numComponents =  imageIO->GetNumberOfComponents();
	DataType dataType = imageIO->GetComponentType();
	int numDimensions = imageIO->GetNumberOfDimensions();
	int dim3size = 1;

	if( numComponents > 12 )	//Can't handle more than 12 channels
	{
		itk::ExceptionObject excp;
		excp.SetDescription("TOO MANY COMPONENTS");
		throw excp;
		return false;
	}

	if(imageDataPtrs.size() > 0)	//Already have data, make sure I can load this new image:
	{
		if(m_Info.dataType != dataType)
			return false;

		for(int d=0; d<numDimensions; ++d)
		{
			if(d==0 && m_Info.numColumns != imageIO->GetDimensions(d))
				return false;
			if(d==1 && m_Info.numRows != imageIO->GetDimensions(d))
				return false;
			if(d==2)
				dim3size = imageIO->GetDimensions(d);
			if(d>=3)
				return false;
		}
		if(stacksAreForTime)
		{
			if(m_Info.numChannels != numDimensions)
				return false;
			if(appendChannels && m_Info.numTSlices != dim3size)
				return false;
		}
		else
		{
			if(appendChannels && m_Info.numZSlices != dim3size)
				return false;
		}
	}
	else
	{
		//My first image so set some stuff (rest will be set later):
		m_Info.dataType = dataType;
		m_Info.spacing.assign(3,1);
	}

	switch( m_Info.dataType )
    {
		case itk::ImageIOBase::UCHAR:
			LoadImageITK<unsigned char>(fileName, numComponents, pixelType, stacksAreForTime, appendChannels);
		break;
		case itk::ImageIOBase::CHAR:
			LoadImageITK<char>(fileName, numComponents, pixelType, stacksAreForTime, appendChannels);
		break;
		case itk::ImageIOBase::USHORT:
			LoadImageITK<unsigned short>(fileName, numComponents, pixelType, stacksAreForTime, appendChannels);
		break;
		case itk::ImageIOBase::SHORT:
			LoadImageITK<short>(fileName, numComponents, pixelType, stacksAreForTime, appendChannels);
		break;
		case itk::ImageIOBase::UINT:
			LoadImageITK<unsigned int>(fileName, numComponents, pixelType, stacksAreForTime, appendChannels);
		break;
		case itk::ImageIOBase::INT:
			LoadImageITK<int>(fileName, numComponents, pixelType, stacksAreForTime, appendChannels);
		break;
		case itk::ImageIOBase::ULONG:
			LoadImageITK<unsigned long>(fileName, numComponents, pixelType, stacksAreForTime, appendChannels);
		break;
		case itk::ImageIOBase::LONG:
			LoadImageITK<long>(fileName, numComponents, pixelType, stacksAreForTime, appendChannels);
		break;
		case itk::ImageIOBase::FLOAT:
			LoadImageITK<float>(fileName, numComponents, pixelType, stacksAreForTime, appendChannels);
		break;
		case itk::ImageIOBase::DOUBLE:
			LoadImageITK<double>(fileName, numComponents, pixelType, stacksAreForTime, appendChannels);
		break;
		default:
			itk::ExceptionObject excp;
			excp.SetDescription("Pixel type is not supported in this application");
			throw excp;
			return false;
		break;
    }

	path = this->GetPath(fileName);							//Path to this image file
	filenames.push_back( this->GetFilename(fileName) );		//Filename of this image

	return true;
}


#if VTK_MAJOR_VERSION <= 5
bool Image::LoadLSMImage( std::string fileName )
{ 
	DeleteData();

	vtkSmartPointer<vtkLSMReader> lsmR = vtkSmartPointer<vtkLSMReader>::New();
	lsmR->SetFileName(fileName.c_str());
	lsmR->OpenFile();
	lsmR->Update();

	m_Info.spacing.clear();
	double *sp = lsmR->GetDataSpacing();
	for(int i=0; i<3; ++i)
	{
		m_Info.spacing.push_back( (float)sp[0] );
	}
  
	m_Info.numChannels = lsmR->GetNumberOfChannels();
	m_Info.numTSlices = lsmR->GetNumberOfTimePoints();
	m_Info.channelColors.clear();
	m_Info.channelNames.clear();

	vtkSmartPointer<vtkImageData> vimdata;
	int extent[6];
	int numBytes = 0;

	imageDataPtrs.clear();
	imageDataPtrs.resize( m_Info.numTSlices );		//First level is Time
	for (itk::SizeValueType t=0; t < m_Info.numTSlices; ++t)
	{
		for (int ch=0; ch < this->m_Info.numChannels; ch++)
		{
			//First get the time point and channel image
  			vimdata = lsmR->GetTimePointOutput(t, ch);
			lsmR->Update();
			vimdata->Update();

			if(t == 0)	//Only do this on the first time through t
  			{
				if(ch == 0)
				{
					//std::cerr << "Scalar Type = " << vimdata->GetScalarTypeAsString() << std::endl;
					vimdata->GetExtent(extent); 
					m_Info.numZSlices = extent[5]-extent[4]+1;
					m_Info.numRows = extent[3]-extent[2]+1;
					m_Info.numColumns = extent[1]-extent[0]+1;
					m_Info.bytesPerPix = vimdata->GetScalarSize();	//number of Bytes per pixel
					m_Info.dataType = GetDataTypeITK( vimdata->GetScalarType() );
					numBytes = (m_Info.numZSlices)*(m_Info.numRows)*(m_Info.numColumns)*(m_Info.bytesPerPix);
				}
				//Assign RGB values for channel
				std::vector< unsigned char > color(3,0);
				color[0] = (unsigned char)(lsmR->GetChannelColorComponent(ch,0));
				color[1] = (unsigned char)(lsmR->GetChannelColorComponent(ch,1));
				color[2] = (unsigned char)(lsmR->GetChannelColorComponent(ch,2));
				m_Info.channelColors.push_back(color);

				//Get the name of each channel
				char *name = lsmR->GetChannelName(ch);
				m_Info.channelNames.push_back(name);
			}

			//Get a scalar pointer to the image data and copy it
			void * p = vimdata->GetScalarPointer();
			//void * mem = new unsigned char [numBytes];
			void * mem = malloc(numBytes);
			if(mem == NULL)
				return false;
			memcpy(mem, p,numBytes);
			ImageMemoryBlock block;
			block.manager = FTK;
			block.mem = mem;
			imageDataPtrs[t].push_back( block );
		}
	}

	//Now delete unneeded data (smart pointers do the work)
	lsmR = 0;
	vimdata = 0;

	path = this->GetPath(fileName);							//Path to this image file
	filenames.push_back( this->GetFilename(fileName) );		//Filename of this image

	return true;
}
#endif

std::vector< itk::SizeValueType > Image::Size(void)
{
	std::vector< itk::SizeValueType > rVal;
	rVal.push_back( m_Info.numTSlices );
	rVal.push_back( m_Info.numZSlices );
	rVal.push_back( m_Info.numRows );
	rVal.push_back( m_Info.numColumns );
	return rVal;
}

void * Image::GetDataPtr(itk::SizeValueType T, itk::SizeValueType CH, PtrMode mode)
{
	if( T > m_Info.numTSlices || CH > m_Info.numChannels )
		return NULL;

	void * mem = NULL;	//This will point to the data I want to pass back;
	if(mode == DEFAULT)
	{
		mem = imageDataPtrs[T][CH].mem;
	}
	else if(mode == RELEASE_CONTROL)
	{
		if(imageDataPtrs[T][CH].manager == FTK)	//I am able to release management
		{
			imageDataPtrs[T][CH].manager = OTHER;
			mem = imageDataPtrs[T][CH].mem;
		}
	}
	else if(mode == DEEP_COPY)
	{
		itk::SizeValueType numBytes = m_Info.BytesPerChunk();
		//mem = (void *)new unsigned char[numBytes];
		mem = malloc( numBytes );
		if(mem == NULL)
			return (void *)NULL;
		memcpy(mem,imageDataPtrs[T][CH].mem,numBytes);
	}
	return mem;
}

Image::VtkImagePtr Image::GetVtkPtr(itk::SizeValueType T, itk::SizeValueType CH, PtrMode mode)
{
	if( T >= m_Info.numTSlices || CH >= m_Info.numChannels )
		return NULL;

	bool makeCopy;
	bool vtkManageMemory;

	if(mode == RELEASE_CONTROL)
	{
		makeCopy = false;
		vtkManageMemory = true;
	}
	else if(mode == DEEP_COPY)
	{
		makeCopy = true;
		vtkManageMemory = true;
	}
	else
	{
		makeCopy = false;
		vtkManageMemory = false;
	}

	vtkSmartPointer<vtkDataArray> m_array = this->GetVtkDataArray(T,CH, makeCopy, vtkManageMemory);
	if (m_array == NULL)
	{
		//FAILED (PROBABLY BECAUSE MEMORY HAS ALREADY BEEN RELEASED
		return NULL;
	}

	VtkImagePtr imageData = VtkImagePtr::New();
#if VTK_MAJOR_VERSION <= 5
	imageData->SetScalarType( GetDataTypeVTK(m_Info.dataType) );
#else
	imageData->SetScalarType( GetDataTypeVTK(m_Info.dataType), NULL );
#endif
	imageData->GetPointData()->SetScalars( m_array );
	imageData->SetDimensions(m_Info.numColumns, m_Info.numRows, m_Info.numZSlices);
	imageData->SetSpacing(m_Info.spacing.at(0), m_Info.spacing.at(1), m_Info.spacing.at(2));
	imageData->SetOrigin(0.0, 0.0, 0.0);

	return imageData;
}

//Add new time slices to myself if all other sizes match.
//Note that any pointer modes can be applied:
// 1. DEFAULT: old image (img) still manages the memory
// 2. RELEASE_CONTROL: this image (this) gets control of the memory
// 3. DEEP_COPY: both images have their own copy of the memory
// It is up to the programmer to use these correctly!!!
bool Image::AppendImage( ftk::Image::Pointer img, PtrMode mode, bool isforOneTime)// overloaded function
{
	if (!isforOneTime)
		return false;

	const Image::Info * in_size = img->GetImageInfo();

	//First check sizes to be sure that they match!!!
	if( m_Info.numChannels != in_size->numChannels ||
		m_Info.numZSlices != in_size->numZSlices ||
		m_Info.numRows != in_size->numRows ||
		m_Info.numColumns != in_size->numColumns ||
		m_Info.dataType != in_size->dataType )
	{
		return false;
	}

	ImageMemoryBlock block;
	if(mode == DEFAULT)
		block.manager = OTHER;
	else
		block.manager = FTK;

	//Resize the data pointers:
	m_Info.numTSlices += in_size->numTSlices;

	imageDataPtrs.resize(m_Info.numTSlices);

	itk::SizeValueType t = m_Info.numTSlices - 1;
	if(m_Info.numTSlices==2)
		FileNames.push_back(this->filenames);
	FileNames.push_back(img->GetFilenames());
	//filenames.push_back( this->GetFilename(fileName.at(0)));		//Filename of this image
	for (int ch=0; ch<m_Info.numChannels; ++ch)
	{
		block.mem = img->GetDataPtr(0,ch,mode);
		imageDataPtrs.at(t).push_back(block);
	}
	return true;
}

bool Image::AppendImage(ftk::Image::Pointer img, PtrMode mode)
{
	const Image::Info * in_size = img->GetImageInfo();

	//First check sizes to be sure that they match!!!
	if( m_Info.numChannels != in_size->numChannels ||
		m_Info.numZSlices != in_size->numZSlices ||
		m_Info.numRows != in_size->numRows ||
		m_Info.numColumns != in_size->numColumns ||
		m_Info.dataType != in_size->dataType )
	{
		return false;
	}

	ImageMemoryBlock block;
	if(mode == DEFAULT)
		block.manager = OTHER;
	else
		block.manager = FTK;

	//Resize the data pointers:
	m_Info.numTSlices += in_size->numTSlices;
	imageDataPtrs.resize(m_Info.numTSlices);

	for( itk::SizeValueType t=0; t<in_size->numTSlices; ++t)
	{
		for (int ch=0; ch<m_Info.numChannels; ++ch)
		{
			block.mem = img->GetDataPtr(0,ch,mode);
			imageDataPtrs.at(t).push_back(block);
		}
	}
	return true;
}

//*********************************************************************
//Will append a new CH to existing data
//Assumes only one 1 time point in the image
// dptr points to the data
// dataType comes from DataType and tells type of pixels
// cs is number of columns (x)
// rs is number of rows (y)
// zs is number of stacks (z)
// name is the name of this channel
// color is a 3 element vector assigning RGB values to this channel
// copy = true will make a copy of the data, copy = false will just assign pointer
// ftk::Image will manage the memory from this point on.
//*********************************************************************
bool Image::AppendChannelFromData3D(void *dptr, DataType dataType, int bpPix, itk::SizeValueType cs, itk::SizeValueType rs, itk::SizeValueType zs, std::string name, std::vector<unsigned char> color, bool copy)
{
	if( imageDataPtrs.size() == 0 )	//This is my first image, so set the size info
	{
		m_Info.numColumns = cs;		//Number of Columns in Image (x)
		m_Info.numRows = rs;			//Number of Rows in Image (y)
		m_Info.numZSlices = zs;		//Number of Z Slices in Image (z)
		m_Info.bytesPerPix = bpPix;	//Number of bytes per pixel (UCHAR - 1 or USHORT - 2)
		m_Info.dataType = dataType;		//See ENUM DataType
		m_Info.numTSlices = 1;		//1 time point assumed
		m_Info.numChannels = 0;		//Start at 0, will be incremented later
		m_Info.spacing.assign(3,1);			//Init the spacing
		imageDataPtrs.resize(1);
	}
	else	//Check to be sure sizes match!
	{
		if(m_Info.numColumns != cs || m_Info.numRows != rs || m_Info.numZSlices != zs)
			return false;
		if(m_Info.bytesPerPix != bpPix || m_Info.dataType != dataType || m_Info.numTSlices != 1)
			return false;
	}

	m_Info.numChannels++;					//Increment the Number of Channels in this Image (ch)
	m_Info.channelColors.push_back(color);	//Use provided color
	m_Info.channelNames.push_back(name);		//Use provided name

	ImageMemoryBlock block;
	block.manager = FTK;
	if(copy)
	{
		//unsigned char *imgdata = new unsigned char[zs*cs*rs*bpPix];
		void *imgdata = malloc(zs*cs*rs*bpPix);
		if(!imgdata) return false;
		memcpy(imgdata,dptr,zs*cs*rs*bpPix);
		block.mem = imgdata;
	}
	else
	{
		block.mem = dptr;
	}
	imageDataPtrs[0].push_back( block );

	return true;
}

bool Image::AppendImageFromData3D(void *dptr, DataType dataType, int bpPix, itk::SizeValueType cs, itk::SizeValueType rs, itk::SizeValueType zs, std::string name, bool copy)
{
	if( imageDataPtrs.size() == 0 )	//This is my first image, so set the size info
	{
		m_Info.numColumns = cs;		//Number of Columns in Image (x)
		m_Info.numRows = rs;			//Number of Rows in Image (y)
		m_Info.numZSlices = zs;		//Number of Z Slices in Image (z)
		m_Info.bytesPerPix = bpPix;	//Number of bytes per pixel (UCHAR - 1 or USHORT - 2)
		m_Info.dataType = dataType;		//See ENUM DataType
		m_Info.numTSlices = 1;		//1 time point assumed
		m_Info.numChannels = 1;		//Start at 1, if to be incremented change code below (was meant for label images only)
		m_Info.spacing.assign(3,1);			//Init the spacing
		//imageDataPtrs.clear();
	}
	else	//Check to be sure sizes match!
	{
		if(m_Info.numColumns != cs || m_Info.numRows != rs || m_Info.numZSlices != zs)
			return false;
		if(m_Info.bytesPerPix != bpPix || m_Info.dataType != dataType)
			return false;
	}

	filenames.push_back(name);		//Use provided name

	ImageMemoryBlock block;
	block.manager = FTK;
	if(copy)
	{
		//unsigned char *imgdata = new unsigned char[zs*cs*rs*bpPix];
		void *imgdata = malloc(zs*cs*rs*bpPix);
		if(!imgdata) return false;
		memcpy(imgdata,dptr,zs*cs*rs*bpPix);
		block.mem = imgdata;
	}
	else
	{
		block.mem = dptr;
	}
	std::vector<ImageMemoryBlock> ChannelData;
	ChannelData.push_back(block);
	imageDataPtrs.push_back(ChannelData);
	m_Info.numTSlices = imageDataPtrs.size();			//Increment the Number of Time Slices in this Image (T)

	return true;
}
bool Image::SaveChannelAs( int channel, std::string baseName, std::string ext )
{
	if( imageDataPtrs.size() == 0 ) return false;

	bool retVal = false;
	switch(m_Info.dataType)
	{
		case itk::ImageIOBase::CHAR:
			if( this->WriteImageITK<char>(channel, baseName, ext) ) retVal = true;
		break;
		case itk::ImageIOBase::UCHAR:
			if( this->WriteImageITK<unsigned char>(channel, baseName, ext) ) retVal = true;
		break;
		case itk::ImageIOBase::SHORT:
			if( this->WriteImageITK<short>(channel, baseName, ext) ) retVal = true;
		break;
		case itk::ImageIOBase::USHORT:
			if( this->WriteImageITK<unsigned short>(channel, baseName, ext) ) retVal = true;
		break;
		case itk::ImageIOBase::INT:
			if( this->WriteImageITK<int>(channel, baseName, ext) ) retVal = true;
		break;
		case itk::ImageIOBase::UINT:
			if( this->WriteImageITK<unsigned int>(channel, baseName, ext) ) retVal = true;
		break;
		case itk::ImageIOBase::LONG:
			if( this->WriteImageITK<long>(channel, baseName, ext) ) retVal = true;
		break;
		case itk::ImageIOBase::ULONG:
			if( this->WriteImageITK<unsigned long>(channel, baseName, ext) ) retVal = true;
		break;
		case itk::ImageIOBase::FLOAT:
			if( this->WriteImageITK<float>(channel, baseName, ext) ) retVal = true;
		break;
		case itk::ImageIOBase::DOUBLE:
			if( this->WriteImageITK<double>(channel, baseName, ext) ) retVal = true;
		break;
    //just silencing a warning for now...
    case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
    break;
	}
	return retVal;
}

//************************************************************************************
//Private Utilities:
std::string Image::GetFileExtension(std::string f)
{
	std::string ext;
	size_t found;
	found = f.find_last_of(".");
	ext = f.substr(found + 1);
	return ext;
}

std::string Image::GetPath(std::string f)
{
	std::string ext;
	size_t found;
	found = f.find_last_of("/\\");
	ext = f.substr(0,found);
	return ext;
}

std::string Image::GetFilename(std::string f)
{
	std::string ext;
	size_t found;
	found = f.find_last_of("/\\");
	ext = f.substr(found + 1);
	return ext;
}

std::string Image::itoa(const int x)
{
  std::ostringstream o;
  if (!(o << x)) return "ERROR";
  return o.str();
}

int Image::GetDataTypeVTK(DataType itk)
{
	int retVal;

	if( itk == itk::ImageIOBase::CHAR ) retVal = VTK_CHAR;
	else if( itk == itk::ImageIOBase::UCHAR ) retVal = VTK_UNSIGNED_CHAR;
	else if( itk == itk::ImageIOBase::SHORT ) retVal = VTK_SHORT;
	else if( itk == itk::ImageIOBase::USHORT ) retVal = VTK_UNSIGNED_SHORT;
	else if( itk == itk::ImageIOBase::INT ) retVal = VTK_INT;
	else if( itk == itk::ImageIOBase::UINT ) retVal = VTK_UNSIGNED_INT;
	else if( itk == itk::ImageIOBase::LONG ) retVal = VTK_LONG;
	else if( itk == itk::ImageIOBase::ULONG ) retVal = VTK_UNSIGNED_LONG;
	else if( itk == itk::ImageIOBase::FLOAT ) retVal = VTK_FLOAT;
	else if( itk == itk::ImageIOBase::DOUBLE ) retVal = VTK_DOUBLE;
	else retVal = VTK_VOID;
	
	return retVal;	
}

Image::DataType Image::GetDataTypeITK(int vtk_type)
{
	DataType retVal;

	if( vtk_type == VTK_CHAR ) retVal = itk::ImageIOBase::CHAR;
	else if( vtk_type == VTK_UNSIGNED_CHAR ) retVal = itk::ImageIOBase::UCHAR;
	else if( vtk_type == VTK_SHORT ) retVal = itk::ImageIOBase::SHORT;
	else if( vtk_type == VTK_UNSIGNED_SHORT ) retVal = itk::ImageIOBase::USHORT;
	else if( vtk_type == VTK_INT ) retVal = itk::ImageIOBase::INT;
	else if( vtk_type == VTK_UNSIGNED_INT ) retVal = itk::ImageIOBase::UINT;
	else if( vtk_type == VTK_LONG ) retVal = itk::ImageIOBase::LONG;
	else if( vtk_type == VTK_UNSIGNED_LONG ) retVal = itk::ImageIOBase::ULONG;
	else if( vtk_type == VTK_FLOAT ) retVal = itk::ImageIOBase::FLOAT;
	else if( vtk_type == VTK_DOUBLE ) retVal = itk::ImageIOBase::DOUBLE;
	else retVal = itk::ImageIOBase::UNKNOWNCOMPONENTTYPE;
	
	return retVal;
}

//This function will return a vtkSmartPointer to a vtkDataArray.  The Array will contain image data.
//There are three valid parameter options:
// 1. makeCopy = false, vtkManageMemory = false (default)
//    ftk::Image continues to manage the memory, and vtkDataArray is wrapped around the same memory block.
// 2. makeCopy = false, vtkManageMemory = true
//	  ftk::Image will no longer manage the memory, management is released to vtkImage, still only one memory block.
//    NOTE: If memory has already been release to someone else (ITK) then this operation will return a NULL pointer.
// 3. makeCopy = true, vtkManageMemory = true
//	  we copy the image data to create a new image and vtkImage handles the new data, ftk::Image still handles original data
// 4. makeCopy = true, vtkManageMemory = false (INVALID OPTION)
//	  Option 3 will be used in this case
vtkSmartPointer<vtkDataArray> Image::GetVtkDataArray(itk::SizeValueType T, itk::SizeValueType CH, bool makeCopy = false, bool vtkManageMemory = false)
{
	if( T >= m_Info.numTSlices || CH >= m_Info.numChannels )
		return NULL;

	int numPixels = m_Info.numColumns * m_Info.numRows * m_Info.numZSlices;
	int numBytes = m_Info.bytesPerPix * numPixels;

	void * mem = NULL;	//This will point to the data I want to create the array from;
	if(makeCopy)		//Make a copy of the data	
	{
		//mem = (void *)new unsigned char[numBytes];
		mem = malloc(numBytes);
		if(mem == NULL)
			return NULL;
		memcpy(mem,imageDataPtrs[T][CH].mem,numBytes);
		vtkManageMemory = true;
	}
	else
	{
		mem = imageDataPtrs[T][CH].mem;
	}

	int save = 1;				//vtk DOES NOT manage the memory (default)
	if( vtkManageMemory )
	{
		if( imageDataPtrs[T][CH].manager != FTK )	//I can't manage it because someone else already does
		{
			return NULL;
		}
		else
		{
			imageDataPtrs[T][CH].manager = VTK;
			save = 0;	//vtk DOES manage the memory
		}
	}

	vtkSmartPointer<vtkDataArray> d_array = NULL;

	if (m_Info.dataType == itk::ImageIOBase::CHAR)
	{
		vtkSmartPointer<vtkCharArray> arr = vtkSmartPointer<vtkCharArray>::New();
		arr->SetArray( static_cast<char*>(mem), numPixels, save ); //Last parameter (1) says vtk Image doesn't manange memory
		d_array = static_cast< vtkSmartPointer<vtkDataArray> >(arr);
	}
	else if (m_Info.dataType == itk::ImageIOBase::UCHAR)
	{
		vtkSmartPointer<vtkUnsignedCharArray> arr = vtkSmartPointer<vtkUnsignedCharArray>::New();
		arr->SetArray( static_cast<unsigned char*>(mem), numPixels, save ); //Last parameter (1) says vtk Image doesn't manange memory
		d_array = static_cast< vtkSmartPointer<vtkDataArray> >(arr);
	}
	else if (m_Info.dataType == itk::ImageIOBase::SHORT)
	{
		vtkSmartPointer<vtkShortArray> arr = vtkSmartPointer<vtkShortArray>::New();
		arr->SetArray( static_cast<short*>(mem), numPixels, save ); //Last parameter (1) says vtk Image doesn't manange memory
		d_array = static_cast< vtkSmartPointer<vtkDataArray> >(arr);
	}
	else if (m_Info.dataType == itk::ImageIOBase::USHORT)
	{
		vtkSmartPointer<vtkUnsignedShortArray> arr = vtkSmartPointer<vtkUnsignedShortArray>::New();
		arr->SetArray( static_cast<unsigned short*>(mem), numPixels, save ); //Last parameter (1) says vtk Image doesn't manange memory
		d_array = static_cast< vtkSmartPointer<vtkDataArray> >(arr);
	}
	else if (m_Info.dataType == itk::ImageIOBase::INT)
	{
		vtkSmartPointer<vtkIntArray> arr = vtkSmartPointer<vtkIntArray>::New();
		arr->SetArray( static_cast<int*>(mem), numPixels, save ); //Last parameter (1) says vtk Image doesn't manange memory
		d_array = static_cast< vtkSmartPointer<vtkDataArray> >(arr);
	}
	else if (m_Info.dataType == itk::ImageIOBase::UINT)
	{
		vtkSmartPointer<vtkUnsignedIntArray> arr = vtkSmartPointer<vtkUnsignedIntArray>::New();
		arr->SetArray( static_cast<unsigned int*>(mem), numPixels, save ); //Last parameter (1) says vtk Image doesn't manange memory
		d_array = static_cast< vtkSmartPointer<vtkDataArray> >(arr);
	}
	else if (m_Info.dataType == itk::ImageIOBase::LONG)
	{
		vtkSmartPointer<vtkLongArray> arr = vtkSmartPointer<vtkLongArray>::New();
		arr->SetArray( static_cast<long*>(mem), numPixels, save ); //Last parameter (1) says vtk Image doesn't manange memory
		d_array = static_cast< vtkSmartPointer<vtkDataArray> >(arr);
	}
	else if (m_Info.dataType == itk::ImageIOBase::ULONG)
	{
		vtkSmartPointer<vtkUnsignedLongArray> arr = vtkSmartPointer<vtkUnsignedLongArray>::New();
		arr->SetArray( static_cast<unsigned long*>(mem), numPixels, save ); //Last parameter (1) says vtk Image doesn't manange memory
		d_array = static_cast< vtkSmartPointer<vtkDataArray> >(arr);
	}
	else if (m_Info.dataType == itk::ImageIOBase::FLOAT)
	{
		vtkSmartPointer<vtkFloatArray> arr = vtkSmartPointer<vtkFloatArray>::New();
		arr->SetArray( static_cast<float*>(mem), numPixels, save ); //Last parameter (1) says vtk Image doesn't manange memory
		d_array = static_cast< vtkSmartPointer<vtkDataArray> >(arr);
	}
	else if (m_Info.dataType == itk::ImageIOBase::DOUBLE)
	{
		vtkSmartPointer<vtkDoubleArray> arr = vtkSmartPointer<vtkDoubleArray>::New();
		arr->SetArray( static_cast<double*>(mem), numPixels, save ); //Last parameter (1) says vtk Image doesn't manange memory
		d_array = static_cast< vtkSmartPointer<vtkDataArray> >(arr);
	}
	
	return d_array;
}

void Image::SetDefaultColors(void)
{
	m_Info.channelColors.clear();
	m_Info.channelNames.clear();

	std::vector<unsigned char> gray(3,255);

	std::vector<unsigned char> red(3,0);
	std::vector<unsigned char> green(3,0);
	std::vector<unsigned char> blue(3,0);
	std::vector<unsigned char> cyan(3,255);
	std::vector<unsigned char> magenta(3,255);
	std::vector<unsigned char> yellow(3,255);
	red[0] = 255;
	green[1] = 255;
	blue[2] = 255;
	cyan[0] = 0;
	magenta[1] = 0;
	yellow[2] = 0;

	std::vector< std::vector<unsigned char> > colorWheel;
	std::vector< std::string > colorNames;
	colorWheel.push_back(red);
	colorNames.push_back("red");
	colorWheel.push_back(green);
	colorNames.push_back("green");
	colorWheel.push_back(blue);
	colorNames.push_back("blue");
	colorWheel.push_back(cyan);
	colorNames.push_back("cyan");
	colorWheel.push_back(magenta);
	colorNames.push_back("magenta");
	colorWheel.push_back(yellow);
	colorNames.push_back("yellow");

	if(m_Info.numChannels == 1)
	{
		m_Info.channelColors.push_back(gray);
		m_Info.channelNames.push_back("gray");
	}
	else
	{
		for(int i=0; i<m_Info.numChannels; ++i)
		{
			int c = i%5;
			m_Info.channelColors.push_back( colorWheel.at(c) );
			m_Info.channelNames.push_back( colorNames.at(c) );
		}
	}
}

void Image::SetPixel(itk::SizeValueType T, itk::SizeValueType CH, itk::SizeValueType Z, itk::SizeValueType R, itk::SizeValueType C, double newValue)
{
	if( T >= m_Info.numTSlices || CH >= m_Info.numChannels || Z >= m_Info.numZSlices \
		|| R >= m_Info.numRows || C >= m_Info.numColumns )
		return;
		
	unsigned int x = m_Info.numColumns;
	unsigned int y = m_Info.numRows;

	if (m_Info.dataType == itk::ImageIOBase::CHAR)
	{
			char value = static_cast<char>(newValue);
			char * p = static_cast<char *>(imageDataPtrs[T][CH].mem) + Z*y*x + R*x + C;
			(*p) = value;
	}
	else if (m_Info.dataType == itk::ImageIOBase::UCHAR)
	{
			unsigned char value = static_cast<unsigned char>(newValue);
			unsigned char * p = static_cast<unsigned char *>(imageDataPtrs[T][CH].mem) + Z*y*x + R*x + C;
			(*p) = value;
	}
	else if (m_Info.dataType == itk::ImageIOBase::SHORT)
	{
			short value = static_cast<short>(newValue);
			short * p = static_cast<short *>(imageDataPtrs[T][CH].mem) + Z*y*x + R*x + C;
			(*p) = value;
	}
	else if (m_Info.dataType == itk::ImageIOBase::USHORT)
	{
			unsigned short value = static_cast<unsigned short>(newValue);
			unsigned short * p = static_cast<unsigned short *>(imageDataPtrs[T][CH].mem) + Z*y*x + R*x + C;
			(*p) = value;
	}
	else if (m_Info.dataType == itk::ImageIOBase::INT)
	{
			int value = static_cast<int>(newValue);
			int * p = static_cast<int *>(imageDataPtrs[T][CH].mem) + Z*y*x + R*x + C;
			(*p) = value;
	}
	else if (m_Info.dataType == itk::ImageIOBase::UINT)
	{
			unsigned int value = static_cast<unsigned int>(newValue);
			unsigned int * p = static_cast<unsigned int *>(imageDataPtrs[T][CH].mem) + Z*y*x + R*x + C;
			(*p) = value;
	}
	else if (m_Info.dataType == itk::ImageIOBase::LONG)
	{
			long value = static_cast<long>(newValue);
			long * p = static_cast<long *>(imageDataPtrs[T][CH].mem) + Z*y*x + R*x + C;
			(*p) = value;
	}
	else if (m_Info.dataType == itk::ImageIOBase::ULONG)
	{
			unsigned long value = static_cast<unsigned long>(newValue);
			unsigned long * p = static_cast<unsigned long *>(imageDataPtrs[T][CH].mem) + Z*y*x + R*x + C;
			(*p) = value;
	}
	else if (m_Info.dataType == itk::ImageIOBase::FLOAT)
	{
			float value = static_cast<float>(newValue);
			float * p = static_cast<float *>(imageDataPtrs[T][CH].mem) + Z*y*x + R*x + C;
			(*p) = value;
	}
	else if (m_Info.dataType == itk::ImageIOBase::DOUBLE)
	{
			double value = static_cast<double>(newValue);
			double * p = static_cast<double *>(imageDataPtrs[T][CH].mem) + Z*y*x + R*x + C;
			(*p) = value;
	}
}

//*****************************************************************************************
// GetPixel
//*****************************************************************************************
double Image::GetPixel(itk::SizeValueType T, itk::SizeValueType CH, itk::SizeValueType Z, itk::SizeValueType R, itk::SizeValueType C)
{
	return this->GetPixelT<double>(T,CH,Z,R,C);
}

}  // end namespace ftk
