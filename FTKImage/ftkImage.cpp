#include "ftkImage.h"

//****************************************************************************************
// FTKImage is a class the handles the storage of multi-dimensional images.
// 
// The class is not complete, and should be updated and modified as needs arise.
//
// It currently loads an image from a file using the vtkKW Image Loader or the LSM
// image loader.  A new image can also be created from data with a single channel
// 
// The images are stored in memory using 3D ITK Images for each channel and time-slice".
// If only a single 2D slice is desired, the image will be cast into a 2D image to be returned.
// Currently Images can either be 16 bit or 8 bit.  They will be stored in one of these two
// formats. 
//****************************************************************************************

namespace ftk
{

//Constructor
Image::Image()
{
}

//Destructor
Image::~Image()
{
}

void Image::LoadFiles( std::vector< std::string > fNames )
{
	for (int i = 0; i < (int)fNames.size(); ++i)
	{
		this->LoadFile( fNames.at(i) );
	}
}

bool Image::LoadFile( std::string fName )
{
	if( GetFileExtension(fName) == "lsm" )
	{
		return this->LoadLSMImage( fName );
	}
	else
	{
		return this->LoadStandardImage( fName );
	}
}

bool Image::LoadStandardImage( std::string fileName )
{
	//Create a new reader to automatically read the image file.
	//It can be extended to read a list of image files as time points.
	vtkKWImageIO *reader = vtkKWImageIO::New();
	reader->SetFileName(fileName);
	try
	{
		reader->ReadImage();
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cerr << excp << std::endl;
		return false;
	}

	itk::ImageIOBase::IOPixelType pixtype = reader->GetImagePixelType();
    itk::ImageIOBase::IOComponentType comptype = reader->GetImageComponentType();

	vtkKWImage *kwImage = reader->HarvestReadImage();
	vtkImageData *vtkImage = kwImage->GetVTKImage();

	//Fill in Common Image Info:
	int extent[6];
	vtkImage->GetExtent(extent);
	imageInfo.path = this->GetPath( fileName );
	imageInfo.filename = this->GetFilename( fileName );
	imageInfo.numColumns = extent[1]-extent[0]+1;		//Number of Columns in Image (x)
	imageInfo.numRows  = extent[3]-extent[2]+1;			//Number of Rows in Image (y)
	imageInfo.numZSlices  = extent[5]-extent[4]+1;		//Number of Z Slices in Image (z)
	imageInfo.numTSlices = 1;							//Number of Time Slices in Image (t)
	imageInfo.bytesPerPix = vtkImage->GetScalarSize();	//Number of bytes per pixel (UCHAR - 1 or USHORT - 2)
	imageInfo.dataType = vtkImage->GetScalarType();		//See ENUM ImageDataType
	imageInfo.channelColors.clear();
	imageInfo.channelNames.clear();
	imageDataPtrs.clear();
	imageDataPtrs.resize(1);

	int numBytes = (imageInfo.numZSlices)*(imageInfo.numRows)*(imageInfo.numColumns)*(imageInfo.bytesPerPix);

	if( pixtype == itk::ImageIOBase::SCALAR )
	{
		imageDataPtrs[0].push_back( vtkImage->GetScalarPointer() );
		imageInfo.numChannels = 1;							//Number of Channels in this Image (ch)
		std::vector<unsigned char> gray(3,255);
		imageInfo.channelColors.push_back(gray);
		imageInfo.channelNames.push_back("gray");
	}
	else if( pixtype == itk::ImageIOBase::RGB && comptype == itk::ImageIOBase::UCHAR)
	{
		//This is more complicated because I need to split the color channels:
		//First Get ITK Pointer:
		typedef itk::RGBPixel< unsigned char > RGBPixelType;
		typedef itk::Image< RGBPixelType, 3 > CImageType3D;
		typedef itk::ImageRegionConstIterator< CImageType3D > CIteratorType3D;

		unsigned char * R = new unsigned char [numBytes];
		unsigned char * G = new unsigned char [numBytes];
		unsigned char * B = new unsigned char [numBytes];
		int b;

		CImageType3D::ConstPointer ptr = static_cast<const CImageType3D *>( kwImage->GetITKImageBase() );
		CIteratorType3D it( ptr, ptr->GetBufferedRegion() );
		for ( it.GoToBegin(), b=0; !it.IsAtEnd(); ++it, ++b ) 
		{
			//create a pixel object
			CImageType3D::PixelType pixelValue = it.Get();

			R[b] = (unsigned char)pixelValue[0]; //Red
			G[b] = (unsigned char)pixelValue[1]; //Green
			B[b] = (unsigned char)pixelValue[2]; //Blue
		}

		imageDataPtrs[0].push_back( (void *)R );
		imageDataPtrs[0].push_back( (void *)G );
		imageDataPtrs[0].push_back( (void *)B );
		imageInfo.numChannels = 3;
				
		//Assign RGB values for channel
		std::vector<unsigned char> red(3,0);
		std::vector<unsigned char> green(3,0);
		std::vector<unsigned char> blue(3,0);
		red[0] = 255;
		green[1] = 255;
		blue[2] = 255;

		imageInfo.channelColors.push_back(red);
		imageInfo.channelColors.push_back(green);
		imageInfo.channelColors.push_back(blue);

		//Get the name of each channel
		imageInfo.channelNames.push_back("red");
		imageInfo.channelNames.push_back("green");
		imageInfo.channelNames.push_back("blue");
	}
	else
	{
		std::cerr << "Unrecognized pixel type" << std::endl;
		return false;
	}

	//vtkImage->Delete();		//Caused program to crash, so removed
	return true;
}

bool Image::LoadLSMImage( std::string fileName )
{ 
	vtkLSMReader *lsmR = vtkLSMReader::New();
	lsmR->SetFileName(fileName.c_str());
	lsmR->OpenFile();
	lsmR->Update();
  
	this->imageInfo.numChannels = lsmR->GetNumberOfChannels();
	this->imageInfo.numTSlices = lsmR->GetNumberOfTimePoints();
	this->imageInfo.channelColors.clear();
	this->imageInfo.channelNames.clear();

	vtkImageData* vimdata;
	int extent[6];
	int numBytes = 0;

	this->imageDataPtrs.resize( this->imageInfo.numTSlices );		//First level is Time
	for (int t=0; t < this->imageInfo.numTSlices; ++t)
	{
		for (int ch=0; ch < this->imageInfo.numChannels; ch++)
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
					imageInfo.numZSlices = extent[5]-extent[4]+1;
					imageInfo.numRows = extent[3]-extent[2]+1;
					imageInfo.numColumns = extent[1]-extent[0]+1;
					imageInfo.bytesPerPix = vimdata->GetScalarSize();	//number of Bytes per pixel
					imageInfo.dataType = vimdata->GetScalarType();
					numBytes = (imageInfo.numZSlices)*(imageInfo.numRows)*(imageInfo.numColumns)*(imageInfo.bytesPerPix);
				}
				//Assign RGB values for channel
				std::vector< unsigned char > color(3,0);
				color[0] = (unsigned char)(lsmR->GetChannelColorComponent(ch,0));
				color[1] = (unsigned char)(lsmR->GetChannelColorComponent(ch,1));
				color[2] = (unsigned char)(lsmR->GetChannelColorComponent(ch,2));
				imageInfo.channelColors.push_back(color);

				//Get the name of each channel
				char *name = lsmR->GetChannelName(ch);
				imageInfo.channelNames.push_back(name);
			}

			//Get a scalar pointer to the image data and copy it
			void * p = vimdata->GetScalarPointer();
			void * mem = new unsigned char [numBytes];
			memcpy(mem, p,numBytes);
			imageDataPtrs[t].push_back( mem );
		}
	}

	//Now delete unneeded data
	lsmR->Delete();

	return true;
}

std::vector< unsigned short > Image::Size(void)
{
	std::vector< unsigned short > rVal;
	rVal.push_back( imageInfo.numTSlices );
	rVal.push_back( imageInfo.numZSlices );
	rVal.push_back( imageInfo.numRows );
	rVal.push_back( imageInfo.numColumns );
	return rVal;
}

unsigned char* Image::GetSlicePtr(int T, int CH, int Z)
{
	if( T > imageInfo.numTSlices || CH > imageInfo.numChannels || Z > imageInfo.numZSlices )
		return NULL;
	if(imageInfo.dataType != UCHAR)
		return NULL;

	unsigned char* stack = static_cast<unsigned char *>(imageDataPtrs[T][CH]);
	return ( stack + Z*(imageInfo.numColumns)*(imageInfo.numRows)*(imageInfo.bytesPerPix) );
}

void * Image::GetDataPtr(int T, int CH)
{
	if( T > imageInfo.numTSlices || CH > imageInfo.numChannels )
		return NULL;

	return imageDataPtrs[T][CH];
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

}  // end namespace ftk
