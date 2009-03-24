/*=========================================================================

  Program:   KWImage - Kitware Image IO Library
  Module:    $RCSfile: vtkKWImageIO.cxx,v $

  Copyright (c) Kitware, Inc., Insight Consortium.  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "vtkKWImageIO.h"

#include "vtkImageData.h"
#include "vtkObjectFactory.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"

#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"

#include "vtkKWImage.h"

#include <vector>
#include <string>

/** \class ReaderCreator
 *  This helper class will take care of instantiating the appropriate ITK
 *  ImageFileReader class corresponding to the actual pixel type of the input
 *  image. */
template <class TPixel >
class ReaderCreator
{
public:
  
  typedef itk::Image< TPixel, 3 >             ImageType;
  typedef itk::ImageFileReader< ImageType >   ReaderType;
  typedef typename ReaderType::Pointer        ReaderPointer;
  typedef vtkKWImage *                        ImagePointer;

  static ImagePointer
  CreateAndRead( const std::string & filename ) 
    {
    ImagePointer kwImage = vtkKWImage::New();
    ReaderPointer reader = ReaderType::New();
   
    reader->SetFileName( filename );
    reader->Update();

    kwImage->SetITKImageBase( reader->GetOutput() );

    return kwImage;
    }
};

/** \class SeriesReaderCreator
 *  This helper class will take care of instantiating the appropriate ITK
 *  ImageFileReader class corresponding to the actual pixel type of the input
 *  image. */
template <class TPixel >
class SeriesReaderCreator
{
public:
  
  typedef itk::Image< TPixel, 3 >             ImageType;
  typedef itk::ImageSeriesReader< ImageType > ReaderType;
  typedef typename ReaderType::Pointer        ReaderPointer;
  typedef vtkKWImage *                        ImagePointer;

  static ImagePointer
  CreateAndRead( const std::vector<std::string> & filenames ) 
    {
    ImagePointer kwImage = vtkKWImage::New();
    ReaderPointer reader = ReaderType::New();
   
    reader->SetFileNames( filenames );
    reader->Update();

    kwImage->SetITKImageBase( reader->GetOutput() );

    return kwImage;
    }
};

/** \class ReadAndCastImageCreator
  read and cast the image to unsigned char type */
template <class TPixel >
class ReadAndCastImageCreator
{
public:
  
  typedef itk::Image< TPixel, 3 >                InputImageType;
  typedef itk::Image< unsigned char, 3>          OutputImageType;
  typedef itk::RescaleIntensityImageFilter< 
              InputImageType, OutputImageType >  RescaleFilterType;
  typedef typename RescaleFilterType::Pointer    RescaleFilterPointer;
  typedef itk::ImageFileReader< InputImageType > ReaderType;
  typedef typename ReaderType::Pointer           ReaderPointer;
  typedef vtkKWImage *                           ImagePointer;

  static ImagePointer
  CreateAndRead( const std::string & filename ) 
    {
    ImagePointer kwImage = vtkKWImage::New();
    ReaderPointer reader = ReaderType::New();

    RescaleFilterPointer rescaleFilter = RescaleFilterType::New();

    rescaleFilter->SetOutputMinimum( 0 );
    rescaleFilter->SetOutputMaximum( 255 );
   
    reader->SetFileName( filename );
    rescaleFilter->SetInput ( reader->GetOutput());
    rescaleFilter->Update();

    kwImage->SetITKImageBase( rescaleFilter->GetOutput() );

    return kwImage;
    }
};

/** \class WriteCreator
 *  This helper class will take care of instantiating the appropriate ITK
 *  ImageFileWriter class corresponding to the actual pixel type of the input
 *  image. */
template <class TPixel >
class WriteCreator
{
public:
  
  typedef itk::ImageBase< 3 >                 ImageBaseType;
  typedef itk::Image< TPixel, 3 >             ImageType;
  typedef itk::ImageFileWriter< ImageType >   WriterType;
  typedef typename WriterType::Pointer        WriterPointer;

  static void
  CreateAndWrite( const std::string & filename, const vtkKWImage * image ) 
    {
    WriterPointer writer = WriterType::New();
   
    const ImageBaseType * imageBase = image->GetITKImageBase();

    typename ImageType::ConstPointer typedImage = 
      dynamic_cast< const ImageType * >( imageBase );

    if( typedImage.IsNull() )
      {
      itk::ExceptionObject excp;
      excp.SetDescription("Error extracting ITK typed image"
          " from the vtkKWImage");
      throw excp;
      }
    
    writer->SetFileName( filename );
    writer->SetInput( typedImage ); 
    writer->Update();
    }
};


/** This helper macro will instantiate the pipeline for a reader creator for a
 * particular pixel type */
#define ReadMacro( PixelType ) \
  ReaderCreator< PixelType >::CreateAndRead( this->FileName );

/** This helper macro will instantiate the pipeline for a reader creator for a
 * particular pixel type */
#define SeriesReaderCreator( PixelType ) \
  SeriesReaderCreator< PixelType >::CreateAndRead( this->SeriesFileNames );

/** This helper macro will instantiate the pipeline for a reader creator for a
 * particular pixel type and cast it to unsigned char type*/
#define ReadAndCastMacro( PixelType ) \
  ReadAndCastImageCreator< PixelType >::CreateAndRead( this->FileName );

/** This helper macro will instantiate the pipeline fora writer creator for a
 * particular pixel type */
#define WriteMacro( PixelType ) \
  WriteCreator< PixelType >::CreateAndWrite( this->FileName, \
                                             this->ImageToBeWritten );


//----------------------------------------------------------------------------
vtkStandardNewMacro( vtkKWImageIO );
vtkCxxRevisionMacro(vtkKWImageIO, "$Revision: 1.1 $");

//----------------------------------------------------------------------------
vtkKWImageIO::vtkKWImageIO()
{
  this->ImageThatHasBeenRead = NULL;
  this->ImageToBeWritten = NULL;
}

//----------------------------------------------------------------------------
vtkKWImageIO::~vtkKWImageIO()
{

}

//-----------------------------------------------------------------------------
void vtkKWImageIO::SetFileName( const std::string & fileName )
{
  this->FileName = fileName;
}

//-----------------------------------------------------------------------------
void 
vtkKWImageIO::SetSeriesFileNames( const std::vector<std::string> & fileNames )
{
  this->SeriesFileNames = fileNames;
}

//-----------------------------------------------------------------------------
void vtkKWImageIO::SetDirectory( const std::string & directory )
{
  this->Directory = directory;
}

// Read Image 
void vtkKWImageIO::ReadImage()
{
  // Find out the pixel type of the image in file
  // CHANGED 1-30-2009 BY ISAAC ABBOTT:
  //typedef itk::ImageIOBase::IOComponentType  PixelType;
	typedef itk::ImageIOBase::IOComponentType  ComponentType;
	typedef itk::ImageIOBase::IOPixelType  PixelType;

  itk::ImageIOBase::Pointer imageIO = 
    itk::ImageIOFactory::CreateImageIO( this->FileName.c_str(), 
                                   itk::ImageIOFactory::ReadMode );

  if( !imageIO )
    {
    std::cerr << "NO IMAGEIO WAS FOUND" << std::endl;
    return;
    }

  // Now that we found the appropriate ImageIO class, ask it to 
  // read the meta data from the image file.
  imageIO->SetFileName( this->FileName.c_str() );
  imageIO->ReadImageInformation();

  //Modified by ISAAC ABBOTT 1-30-2009
  //REMOVED:
  //	PixelType pixelType = imageIO->GetComponentType();
  //ADDED:
  PixelType pixelType = imageIO->GetPixelType();
  ComponentType componentType = imageIO->GetComponentType();

  // Use the componentType and pixelType to instantiate the appropriate reader
  if ( pixelType == itk::ImageIOBase::SCALAR )
  {
   this->ImagePixelType = pixelType;
   switch( componentType )
    {
    case itk::ImageIOBase::UCHAR:
      {
      this->ImageThatHasBeenRead = ReadMacro( unsigned char ); 
	  this->ImageComponentType = componentType;
      break;
      }
    case itk::ImageIOBase::CHAR:
      {
      this->ImageThatHasBeenRead = ReadMacro( char );
      this->ImageComponentType = componentType;
      break;
      }
    case itk::ImageIOBase::USHORT:
      {
      this->ImageThatHasBeenRead = ReadMacro( unsigned short );
      this->ImageComponentType = componentType;
      break;
      }
    case itk::ImageIOBase::SHORT:
      {
      this->ImageThatHasBeenRead = ReadMacro( short );
      this->ImageComponentType = componentType;
      break;
      }
    case itk::ImageIOBase::UINT:
      {
      this->ImageThatHasBeenRead = ReadMacro( unsigned int );
      this->ImageComponentType = componentType;
      break;
      }
    case itk::ImageIOBase::INT:
      {
      this->ImageThatHasBeenRead = ReadMacro( int );
      this->ImageComponentType = componentType;
      break;
      }
    case itk::ImageIOBase::ULONG:
      {
      this->ImageThatHasBeenRead = ReadMacro( unsigned long );
      this->ImageComponentType = componentType;
      break;
      }
    case itk::ImageIOBase::LONG:
      {
      this->ImageThatHasBeenRead = ReadMacro( long );
      this->ImageComponentType = componentType;
      break;
      }
    case itk::ImageIOBase::FLOAT:
      {
      this->ImageThatHasBeenRead = ReadMacro( float );
      this->ImageComponentType = componentType;
      break;
      }
    case itk::ImageIOBase::DOUBLE:
      {
      this->ImageThatHasBeenRead = ReadMacro( double );
      this->ImageComponentType = componentType;
      break;
      }
    default:
      {
      itk::ExceptionObject excp;
      excp.SetDescription("The file uses a pixel type"
            "that is not supported in this application");
      throw excp;
      }
    }
  }
  else if( pixelType == itk::ImageIOBase::RGB )
  {
	this->ImagePixelType = pixelType;
    switch( componentType )
    {
    case itk::ImageIOBase::UCHAR:
      {
      this->ImageThatHasBeenRead = ReadMacro( itk::RGBPixel<unsigned char> ); 
	  this->ImageComponentType = componentType;
      break;
      }
    case itk::ImageIOBase::CHAR:
      {
      this->ImageThatHasBeenRead = ReadMacro( itk::RGBPixel<char> );
      this->ImageComponentType = componentType;
      break;
      }
    case itk::ImageIOBase::USHORT:
      {
      this->ImageThatHasBeenRead = ReadMacro( itk::RGBPixel<unsigned short> );
      this->ImageComponentType = componentType;
      break;
      }
    case itk::ImageIOBase::SHORT:
      {
      this->ImageThatHasBeenRead = ReadMacro( itk::RGBPixel<short> );
      this->ImageComponentType = componentType;
      break;
      }
    case itk::ImageIOBase::UINT:
      {
      this->ImageThatHasBeenRead = ReadMacro( itk::RGBPixel<unsigned int> );
      this->ImageComponentType = componentType;
      break;
      }
    case itk::ImageIOBase::INT:
      {
      this->ImageThatHasBeenRead = ReadMacro( itk::RGBPixel<int> );
      this->ImageComponentType = componentType;
      break;
      }
    case itk::ImageIOBase::ULONG:
      {
      this->ImageThatHasBeenRead = ReadMacro( itk::RGBPixel<unsigned long> );
      this->ImageComponentType = componentType;
      break;
      }
    case itk::ImageIOBase::LONG:
      {
      this->ImageThatHasBeenRead = ReadMacro( itk::RGBPixel<long> );
      this->ImageComponentType = componentType;
      break;
      }
    case itk::ImageIOBase::FLOAT:
      {
      this->ImageThatHasBeenRead = ReadMacro( itk::RGBPixel<float> );
      this->ImageComponentType = componentType;
      break;
      }
    case itk::ImageIOBase::DOUBLE:
      {
      this->ImageThatHasBeenRead = ReadMacro( itk::RGBPixel<double> );
      this->ImageComponentType = componentType;
      break;
      }
    default:
      {
      itk::ExceptionObject excp;
      excp.SetDescription("The file uses a pixel type"
            "that is not supported in this application");
      throw excp;
      }
    }
  }
  else
  {
      itk::ExceptionObject excp;
      excp.SetDescription("The file uses a pixel type"
            "that is not supported in this application");
      throw excp;
  }
}


// Read Image 
void vtkKWImageIO::ReadImageSeries()
{
  //Extract the series UID of the selected dicom image
  // CHANGED 1-30-2009 BY ISAAC ABBOTT:
  //typedef itk::ImageIOBase::IOComponentType  PixelType;
  typedef itk::ImageIOBase::IOComponentType  ComponentType;
  typedef itk::ImageIOBase::IOPixelType  PixelType;

  /***************************************************************************
  // REMOVE DICOM/GDCM BECAUSE I WANT TO OPEN ANY SERIES OF IMAGES:
  // BY ISAAC ABBOTT3/24/2009
  itk::GDCMImageIO::Pointer imageIO = itk::GDCMImageIO::New();

  if( !imageIO )
    {
    std::cerr << "Failure generating GDCM image IO " << std::endl;
    return;
    }

  // Now that we found the appropriate ImageIO class, ask it to 
  // read the meta data from the image file.
  imageIO->SetFileName( this->FileName.c_str() );
  imageIO->ReadImageInformation();

  itk::MetaDataDictionary & dict = imageIO->GetMetaDataDictionary();
 
  std::string  tagkey = "0020|000e";
  std::string seriesIdentifier;

  if( ! itk::ExposeMetaData<std::string>(dict,tagkey, seriesIdentifier ) )
    {
    std::cerr << "Series UID not found: " << std::endl;
    return; 
    }

  //Generate filenames of the DICOM series that belong
  //to the DICOM volume of the selected image  
  typedef itk::GDCMSeriesFileNames NamesGeneratorType;
  NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();

  nameGenerator->SetInputDirectory( this->Directory );

  typedef std::vector< std::string >   FileNamesContainerType;
  FileNamesContainerType fileNames;
  fileNames = nameGenerator->GetFileNames( seriesIdentifier );
  this->SetSeriesFileNames( fileNames );
  ****************************************************************************/

  itk::ImageIOBase::Pointer imageIO = 
    itk::ImageIOFactory::CreateImageIO( this->SeriesFileNames.at(0).c_str(), 
                                   itk::ImageIOFactory::ReadMode );

  if( !imageIO )
  {
	  std::cerr << "NO IMAGEIO WAS FOUND" << std::endl;
	  return;
  }


  //Modified by ISAAC ABBOTT 1-30-2009
  //REMOVED:
  //	PixelType pixelType = imageIO->GetComponentType();
  //ADDED:
  PixelType pixelType = imageIO->GetPixelType();
  ComponentType componentType = imageIO->GetComponentType();

  if( pixelType != itk::ImageIOBase::SCALAR )
  {
      itk::ExceptionObject excp;
      excp.SetDescription("The file uses a pixel type"
            "that is not supported in this application");
      throw excp;
  }
  
  this->ImagePixelType = itk::ImageIOBase::SCALAR;

  // Use the pixel type to instantiate the appropriate reader
  switch( componentType )
    {
    case itk::ImageIOBase::UCHAR:
      {
      this->ImageThatHasBeenRead = SeriesReaderCreator( unsigned char ); 
      this->ImageComponentType = componentType;
      break;
      }
    case itk::ImageIOBase::CHAR:
      {
      this->ImageThatHasBeenRead = SeriesReaderCreator( char );
      this->ImageComponentType = componentType;
      break;
      }
    case itk::ImageIOBase::USHORT:
      {
      this->ImageThatHasBeenRead = SeriesReaderCreator( unsigned short );
      this->ImageComponentType = componentType;
      break;
      }
    case itk::ImageIOBase::SHORT:
      {
      this->ImageThatHasBeenRead = SeriesReaderCreator( short );
      this->ImageComponentType = componentType;
      break;
      }
    case itk::ImageIOBase::UINT:
      {
      this->ImageThatHasBeenRead = SeriesReaderCreator( unsigned int );
      this->ImageComponentType = componentType;
      break;
      }
    case itk::ImageIOBase::INT:
      {
      this->ImageThatHasBeenRead = SeriesReaderCreator( int );
      this->ImageComponentType = componentType;
      break;
      }
    case itk::ImageIOBase::ULONG:
      {
      this->ImageThatHasBeenRead = SeriesReaderCreator( unsigned long );
      this->ImageComponentType = componentType;
      break;
      }
    case itk::ImageIOBase::LONG:
      {
      this->ImageThatHasBeenRead = SeriesReaderCreator( long );
      this->ImageComponentType = componentType;
      break;
      }
    case itk::ImageIOBase::FLOAT:
      {
      this->ImageThatHasBeenRead = SeriesReaderCreator( float );
      this->ImageComponentType = componentType;
      break;
      }
    case itk::ImageIOBase::DOUBLE:
      {
      this->ImageThatHasBeenRead = SeriesReaderCreator( double );
      this->ImageComponentType = componentType;
      break;
      }
    default:
      {
      itk::ExceptionObject excp;
      excp.SetDescription("The file uses a pixel type"
            "that is not supported in this application");
      throw excp;
      }
    }
}

// Read and cast Image to unsigned char 
void vtkKWImageIO::ReadAndCastImage()
{
  // Find out the pixel type of the image in file
  // CHANGED 1-30-2009 BY ISAAC ABBOTT:
  //typedef itk::ImageIOBase::IOComponentType  PixelType;
  typedef itk::ImageIOBase::IOComponentType  ComponentType;
  typedef itk::ImageIOBase::IOPixelType  PixelType;

  itk::ImageIOBase::Pointer imageIO = 
    itk::ImageIOFactory::CreateImageIO( this->FileName.c_str(), 
                                   itk::ImageIOFactory::ReadMode );

  if( !imageIO )
    {
    std::cerr << "NO IMAGEIO WAS FOUND" << std::endl;
    return;
    }

  // Now that we found the appropriate ImageIO class, ask it to 
  // read the meta data from the image file.
  imageIO->SetFileName( this->FileName.c_str() );
  imageIO->ReadImageInformation();

  //Modified by ISAAC ABBOTT 1-30-2009
  //REMOVED:
  //	PixelType pixelType = imageIO->GetComponentType();
  //ADDED:
  PixelType pixelType = imageIO->GetPixelType();
  ComponentType componentType = imageIO->GetComponentType();

  // Use the pixel type to instantiate the appropriate reader
  // Use the componentType and pixelType to instantiate the appropriate reader
  if ( pixelType == itk::ImageIOBase::SCALAR )
  {
   this->ImagePixelType = pixelType;
   switch( componentType )
    {
    case itk::ImageIOBase::UCHAR:
      {
      //this->ImageThatHasBeenRead = ReadAndCastMacro( unsigned char );
      this->ImageThatHasBeenRead = ReadMacro( unsigned char );
	  this->ImageComponentType = itk::ImageIOBase::UCHAR;
      break;
      }
    case itk::ImageIOBase::CHAR:
      {
      this->ImageThatHasBeenRead = ReadAndCastMacro( char );
      this->ImageComponentType = itk::ImageIOBase::UCHAR;
      break;
      }
    case itk::ImageIOBase::USHORT:
      {
      this->ImageThatHasBeenRead = ReadAndCastMacro( unsigned short );
      this->ImageComponentType = itk::ImageIOBase::UCHAR;
      break;
      }
    case itk::ImageIOBase::SHORT:
      {
      this->ImageThatHasBeenRead = ReadAndCastMacro( short );
      this->ImageComponentType = itk::ImageIOBase::UCHAR;
      break;
      }
    case itk::ImageIOBase::UINT:
      {
      this->ImageThatHasBeenRead = ReadAndCastMacro( unsigned int );
      this->ImageComponentType = itk::ImageIOBase::UCHAR;
      break;
      }
    case itk::ImageIOBase::INT:
      {
      this->ImageThatHasBeenRead = ReadAndCastMacro( int );
      this->ImageComponentType = itk::ImageIOBase::UCHAR;
      break;
      }
    case itk::ImageIOBase::ULONG:
      {
      this->ImageThatHasBeenRead = ReadAndCastMacro( unsigned long );
      this->ImageComponentType = itk::ImageIOBase::UCHAR;
      break;
      }
    case itk::ImageIOBase::LONG:
      {
      this->ImageThatHasBeenRead = ReadAndCastMacro( long );
      this->ImageComponentType = itk::ImageIOBase::UCHAR;
      break;
      }
    case itk::ImageIOBase::FLOAT:
      {
      this->ImageThatHasBeenRead = ReadAndCastMacro( float );
      this->ImageComponentType = itk::ImageIOBase::UCHAR;
      break;
      }
    case itk::ImageIOBase::DOUBLE:
      {
      this->ImageThatHasBeenRead = ReadAndCastMacro( double );
      this->ImageComponentType = itk::ImageIOBase::UCHAR;
      break;
      }
    default:
      {
      itk::ExceptionObject excp;
      excp.SetDescription("The file uses a pixel type"
            "that is not supported in this application");
      throw excp;
      }
    }
  }
  else if ( pixelType == itk::ImageIOBase::RGB && componentType == itk::ImageIOBase::UCHAR )
  {
      this->ImageThatHasBeenRead = ReadMacro( itk::RGBPixel<unsigned char> ); 
	  this->ImageComponentType = componentType;
	  this->ImagePixelType = pixelType;
  }
  else
  {
      itk::ExceptionObject excp;
      excp.SetDescription("The file uses a pixel type"
            "that is not supported in this application");
      throw excp;
  }
}

// return pixel type
itk::ImageIOBase::IOPixelType vtkKWImageIO::GetImagePixelType()
{
  return this->ImagePixelType;
}

itk::ImageIOBase::IOComponentType vtkKWImageIO::GetImageComponentType()
{
  return this->ImageComponentType;
}


// Write Image 
void vtkKWImageIO::WriteImage()
{
  // CHANGED 1-30-2009 BY ISAAC ABBOTT:
  //typedef itk::ImageIOBase::IOComponentType  PixelType;
  typedef itk::ImageIOBase::IOComponentType  ComponentType;
  typedef itk::ImageIOBase::IOPixelType  PixelType;

  PixelType pixelType = this->ImageToBeWritten->GetITKPixelType();
  ComponentType componentType = this->ImageToBeWritten->GetITKComponentType();

  if ( pixelType == itk::ImageIOBase::SCALAR )
  {
   this->ImagePixelType = pixelType;
   switch( componentType )
    {
    case itk::ImageIOBase::UCHAR:
      {
      WriteMacro( unsigned char ); 
      break;
      }
    case itk::ImageIOBase::CHAR:
      {
      WriteMacro( char );
      break;
      }
    case itk::ImageIOBase::USHORT:
      {
      WriteMacro( unsigned short );
      break;
      }
    case itk::ImageIOBase::SHORT:
      {
      WriteMacro( short );
      break;
      }
    case itk::ImageIOBase::UINT:
      {
      WriteMacro( unsigned int );
      break;
      }
    case itk::ImageIOBase::INT:
      {
      WriteMacro( int );
      break;
      }
    case itk::ImageIOBase::ULONG:
      {
      WriteMacro( unsigned long );
      break;
      }
    case itk::ImageIOBase::LONG:
      {
      WriteMacro( long );
      break;
      }
    case itk::ImageIOBase::FLOAT:
      {
      WriteMacro( float );
      break;
      }
    case itk::ImageIOBase::DOUBLE:
      {
      WriteMacro( double );
      break;
      }
    default:
      {
      itk::ExceptionObject excp;
      excp.SetDescription("The file uses a pixel type"
            "that is not supported in this application");
      throw excp;
      }
    }
  }
  else if( pixelType == itk::ImageIOBase::RGB )
  {
	this->ImagePixelType = pixelType;
    switch( componentType )
    {
    case itk::ImageIOBase::UCHAR:
      {
      WriteMacro( itk::RGBPixel<unsigned char> ); 
      break;
      }
    case itk::ImageIOBase::CHAR:
      {
      WriteMacro( itk::RGBPixel<char> );
      break;
      }
    case itk::ImageIOBase::USHORT:
      {
      WriteMacro( itk::RGBPixel<unsigned short> );
      break;
      }
    case itk::ImageIOBase::SHORT:
      {
      WriteMacro( itk::RGBPixel<short> );
      break;
      }
    case itk::ImageIOBase::UINT:
      {
      WriteMacro( itk::RGBPixel<unsigned int> );
      break;
      }
    case itk::ImageIOBase::INT:
      {
      WriteMacro( itk::RGBPixel<int> );
      break;
      }
    case itk::ImageIOBase::ULONG:
      {
      WriteMacro( itk::RGBPixel<unsigned long> );
      break;
      }
    case itk::ImageIOBase::LONG:
      {
      WriteMacro( itk::RGBPixel<long> );
      break;
      }
    case itk::ImageIOBase::FLOAT:
      {
      WriteMacro( itk::RGBPixel<float> );
      break;
      }
    case itk::ImageIOBase::DOUBLE:
      {
      WriteMacro( itk::RGBPixel<double> );
      break;
      }
    default:
      {
      itk::ExceptionObject excp;
      excp.SetDescription("The file uses a pixel type"
            "that is not supported in this application");
      throw excp;
      }
    }
  }
  else
  {
      itk::ExceptionObject excp;
      excp.SetDescription("The file uses a pixel type"
            "that is not supported in this application");
      throw excp;
  }
}

//-----------------------------------------------------------------------------
vtkKWImage * vtkKWImageIO::HarvestReadImage()
{
  vtkKWImage * image = this->ImageThatHasBeenRead;
  this->ImageThatHasBeenRead = NULL;
  return image;
}


//-----------------------------------------------------------------------------
void vtkKWImageIO::SetImageToBeWritten( const vtkKWImage * imageToBeWritten )
{
  this->ImageToBeWritten = imageToBeWritten;
}
