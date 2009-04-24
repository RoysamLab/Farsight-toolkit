#include "ftkImage.h"

//itk includes:
#include <itkImageRegionConstIterator.h>
#include <itkNumericSeriesFileNames.h>

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
#include "vtkLSMReader.h"

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
				void *mem = imageDataPtrs[i][j].mem;
				delete[] mem;
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

bool Image::LoadFileSeries( std::string arg, int start, int end, int step)
{	
	typedef itk::NumericSeriesFileNames NameGeneratorType;
	NameGeneratorType::Pointer nameGenerator = NameGeneratorType::New();
	nameGenerator->SetSeriesFormat( arg );
	nameGenerator->SetStartIndex( start );
	nameGenerator->SetEndIndex( end );
	nameGenerator->SetIncrementIndex( step );
	
	std::vector<std::string> names = nameGenerator->GetFileNames();
	for(int i=0; i<(int)names.size(); ++i)
	{
		if( !this->LoadFile( names.at(i) ) )
			return false;
	}
	return true;
}

//Loads a File normally (2D/3D)
bool Image::LoadFile( std::string fName)
{
	if( GetFileExtension(fName) == "lsm" )
		return this->LoadLSMImage( fName );
	else
		return this->LoadStandardImage( fName, false);
}

//Attempts to load a multi-page image and treat each page as a separate time slice.
//Does not work for LSM images
bool Image::LoadFileAsTimeSeries( std::string fName)
{
	if( GetFileExtension(fName) == "lsm" )
		return false;
	else
		return this->LoadStandardImage( fName, true);
}

//**********************************************************************************************************
// IF IMAGE HAS MULTIPLE PAGES IT IS ASSUMED THAT THEY REPRESENT A 3D IMAGE, BUT THEY CAN BE FORCED TO BE T
//
//***********************************************************************************************************
bool Image::LoadStandardImage( std::string fileName, bool stacksAreForTime = false )
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
	imageIO->SetFileName( fileName.c_str() );
	imageIO->ReadImageInformation();

	itkPixelType pixelType = imageIO->GetPixelType();
	DataType dataType = imageIO->GetComponentType();
	int numComponents = imageIO->GetNumberOfComponents();

	if( numComponents > 6 )	//Can't handle more than 6 channels
	{
		itk::ExceptionObject excp;
		excp.SetDescription("TOO MANY COMPONENTS");
		throw excp;
		return false;
	}

	if(imageDataPtrs.size() > 0)
	{
		if(m_Info.dataType != dataType || m_Info.numChannels != numComponents)
			return false;
	}
	else
	{
		//My first image so set some stuff:
		m_Info.numChannels = numComponents;
		m_Info.dataType = dataType;
		m_Info.spacing.assign(3,1);
	}

	switch( m_Info.dataType )
    {
		case itk::ImageIOBase::UCHAR:
			LoadImageITK<unsigned char>(fileName, pixelType, stacksAreForTime);
		break;
		case itk::ImageIOBase::CHAR:
			LoadImageITK<char>(fileName, pixelType, stacksAreForTime);
		break;
		case itk::ImageIOBase::USHORT:
			LoadImageITK<unsigned short>(fileName, pixelType, stacksAreForTime);
		break;
		case itk::ImageIOBase::SHORT:
			LoadImageITK<short>(fileName, pixelType, stacksAreForTime);
		break;
		case itk::ImageIOBase::UINT:
			LoadImageITK<unsigned int>(fileName, pixelType, stacksAreForTime);
		break;
		case itk::ImageIOBase::INT:
			LoadImageITK<int>(fileName, pixelType, stacksAreForTime);
		break;
		case itk::ImageIOBase::ULONG:
			LoadImageITK<unsigned long>(fileName, pixelType, stacksAreForTime);
		break;
		case itk::ImageIOBase::LONG:
			LoadImageITK<long>(fileName, pixelType, stacksAreForTime);
		break;
		case itk::ImageIOBase::FLOAT:
			LoadImageITK<float>(fileName, pixelType, stacksAreForTime);
		break;
		case itk::ImageIOBase::DOUBLE:
			LoadImageITK<double>(fileName, pixelType, stacksAreForTime);
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
	for (int t=0; t < m_Info.numTSlices; ++t)
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
			void * mem = new unsigned char [numBytes];
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

std::vector< unsigned short > Image::Size(void)
{
	std::vector< unsigned short > rVal;
	rVal.push_back( m_Info.numTSlices );
	rVal.push_back( m_Info.numZSlices );
	rVal.push_back( m_Info.numRows );
	rVal.push_back( m_Info.numColumns );
	return rVal;
}

void * Image::GetDataPtr(int T, int CH, PtrMode mode)
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
		int numBytes = m_Info.numColumns * m_Info.numRows * m_Info.numZSlices * m_Info.bytesPerPix;
		mem = (void *)new unsigned char[numBytes];
		memcpy(mem,imageDataPtrs[T][CH].mem,numBytes);
	}
	return mem;
}

Image::VtkImagePtr Image::GetVtkPtr(int T, int CH, PtrMode mode)
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
	imageData->SetScalarType( GetDataTypeVTK(m_Info.dataType) );
	imageData->GetPointData()->SetScalars( m_array );
	imageData->SetDimensions(m_Info.numColumns, m_Info.numRows, m_Info.numZSlices);
	imageData->SetSpacing(m_Info.spacing.at(0), m_Info.spacing.at(1), m_Info.spacing.at(2));
	imageData->SetOrigin(0.0, 0.0, 0.0);

	int type = imageData->GetScalarType();

	return imageData;
}

//Add new time slices to myself if all other sizes match.
//Note that any pointer modes can be applied:
// 1. DEFAULT: old image (img) still manages the memory
// 2. RELEASE_CONTROL: this image (this) gets control of the memory
// 3. DEEP_COPY: both images have their own copy of the memory
// It is up to the programmer to use these correctly!!!
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

	for( int t=0; t<in_size->numTSlices; ++t)
	{
		for (int ch=0; ch<m_Info.numChannels; ++ch)
		{
			block.mem = img->GetDataPtr(t,ch,mode);
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
bool Image::AppendChannelFromData3D(void *dptr, DataType dataType, int bpPix, int cs, int rs, int zs, std::string name, std::vector<unsigned char> color, bool copy)
{
	if( imageDataPtrs.size() == 0 )	//This is my first image, so set the size info
	{
		m_Info.numColumns = cs;		//Number of Columns in Image (x)
		m_Info.numRows = rs;			//Number of Rows in Image (y)
		m_Info.numZSlices = zs;		//Number of Z Slices in Image (z)
		m_Info.bytesPerPix = bpPix;	//Number of bytes per pixel (UCHAR - 1 or USHORT - 2)
		m_Info.dataType = dataType;		//See ENUM DataType
		m_Info.numTSlices = 1;		//1 time point assumed
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
		unsigned char *imgdata = new unsigned char[zs*cs*rs*bpPix];
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
vtkSmartPointer<vtkDataArray> Image::GetVtkDataArray(int T, int CH, bool makeCopy = false, bool vtkManageMemory = false)
{
	if( T >= m_Info.numTSlices || CH >= m_Info.numChannels )
		return NULL;

	int numPixels = m_Info.numColumns * m_Info.numRows * m_Info.numZSlices;
	int numBytes = m_Info.bytesPerPix * numPixels;

	void * mem = NULL;	//This will point to the data I want to create the array from;
	if(makeCopy)		//Make a copy of the data	
	{
		mem = (void *)new unsigned char[numBytes];
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

void Image::SetPixel(int T, int CH, int Z, int R, int C, double newValue)
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
double Image::GetPixel(int T, int CH, int Z, int R, int C)
{
	if( T >= m_Info.numTSlices || CH >= m_Info.numChannels || Z >= m_Info.numZSlices \
		|| R >= m_Info.numRows || C >= m_Info.numColumns )
		return NULL;

	unsigned int x = m_Info.numColumns;
	unsigned int y = m_Info.numRows;
	unsigned int n = m_Info.bytesPerPix;
	char *p = static_cast<char *>(imageDataPtrs[T][CH].mem) + Z*y*x*n + R*x*n + C*n;	//any 8-bit type works here

	double value = 0; 
	switch(m_Info.dataType)
	{
		case itk::ImageIOBase::CHAR:
			value = GetPixelValue<char, double>(p);
		break;
		case itk::ImageIOBase::UCHAR:
			value = GetPixelValue<unsigned char, double>(p);
		break;
		case itk::ImageIOBase::SHORT:
			value = GetPixelValue<short, double>(p);
		break;
		case itk::ImageIOBase::USHORT:
			value = GetPixelValue<unsigned short, double>(p);
		break;
		case itk::ImageIOBase::INT:
			value = GetPixelValue<int, double>(p);
		break;
		case itk::ImageIOBase::UINT:
			value = GetPixelValue<unsigned int, double>(p);
		break;
		case itk::ImageIOBase::LONG:
			value = GetPixelValue<long, double>(p);
		break;
		case itk::ImageIOBase::ULONG:
			value = GetPixelValue<unsigned long, double>(p);
		break;
		case itk::ImageIOBase::FLOAT:
			value = GetPixelValue<float, double>(p);
		break;
		case itk::ImageIOBase::DOUBLE:
			value = GetPixelValue<double, double>(p);
		break;
		default:
			//These are not supported!
		break;
	}
	return value;
}

}  // end namespace ftk
