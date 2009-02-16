/*=========================================================================

  Program:   KWImage - Kitware Image IO Library
  Module:    $RCSfile: vtkKWImage.cxx,v $

  Copyright (c) Kitware, Inc., Insight Consortium.  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "vtkKWImage.h"

#include "itkImage.h"
#include "itkRGBPixel.h"
#include "itkVTKImageExport.h"

#include "vtkObjectFactory.h"
#include "vtkImageData.h"
#include "vtkImageImport.h"


/** \class PipelineCreator
 *  This helper class will take care of instantiating the appropriate
 *  ITK Export class corresponding to the actual pixel type of the 
 *  input image. */
template < class TPixel >
class PipelineCreator
{
public:
  
  typedef itk::ImageBase<3>           ImageBaseType;
  typedef ImageBaseType::Pointer      ImageBasePointer;
  typedef itk::ProcessObject          ExporterBaseType;
  typedef itk::ProcessObject::Pointer ExporterBasePointer;
  typedef itk::Image< TPixel, 3 >     ImageType;

  static void 
  CreateExporter( ImageBasePointer    & imageBase, 
                  ExporterBasePointer & exporter,
                  vtkImageImport      * importer  )
    {
    ImageType * image = 
      dynamic_cast< ImageType * >( imageBase.GetPointer() );

    if( image )
      {
      typedef itk::VTKImageExport< ImageType >   ExportFilterType;
      typedef typename ExportFilterType::Pointer ExportFilterPointer;
      ExportFilterPointer itkExporter = ExportFilterType::New();
      itkExporter->SetInput( image );

      exporter = itkExporter;

      importer->SetUpdateInformationCallback(
        itkExporter->GetUpdateInformationCallback());
      importer->SetPipelineModifiedCallback(
        itkExporter->GetPipelineModifiedCallback());
      importer->SetWholeExtentCallback(
        itkExporter->GetWholeExtentCallback());
      importer->SetSpacingCallback(
        itkExporter->GetSpacingCallback());
      importer->SetOriginCallback(
        itkExporter->GetOriginCallback());
      importer->SetScalarTypeCallback(
        itkExporter->GetScalarTypeCallback());
      importer->SetNumberOfComponentsCallback(
        itkExporter->GetNumberOfComponentsCallback());
      importer->SetPropagateUpdateExtentCallback(
        itkExporter->GetPropagateUpdateExtentCallback());
      importer->SetUpdateDataCallback(
        itkExporter->GetUpdateDataCallback());
      importer->SetDataExtentCallback(
        itkExporter->GetDataExtentCallback());
      importer->SetBufferPointerCallback(
        itkExporter->GetBufferPointerCallback());
      importer->SetCallbackUserData(
        itkExporter->GetCallbackUserData());
      }
    }
};


/** This helper macro will instantiate the pipeline creator for a particular
 * pixel type */
#define CreatePipelineMacro( PixelType ) \
  PipelineCreator< PixelType >::CreateExporter( \
      this->ItkImage, this->Exporter, this->Importer );

//----------------------------------------------------------------------------
vtkStandardNewMacro( vtkKWImage );
vtkCxxRevisionMacro( vtkKWImage, "$Revision: 1.1 $" );

//----------------------------------------------------------------------------
vtkKWImage::vtkKWImage()
{
  this->Importer = vtkImageImport::New();
}

//----------------------------------------------------------------------------
vtkKWImage::~vtkKWImage()
{
  if( this->Importer )
    { 
    this->Importer->Delete();
    }
}

//----------------------------------------------------------------------------
void vtkKWImage::SetITKImageBase( ImageBaseType * image )
{
  if( this->ItkImage.GetPointer() == image )
    {
    return;
    }

  this->ItkImage = image;
  this->Modified();
  
  CreatePipelineMacro( unsigned char );
  CreatePipelineMacro( char );
  CreatePipelineMacro( unsigned short );
  CreatePipelineMacro( short );
  CreatePipelineMacro( unsigned int );
  CreatePipelineMacro( int );
  CreatePipelineMacro( unsigned long );
  CreatePipelineMacro( long );
  CreatePipelineMacro( float );
  CreatePipelineMacro( double );

  //**********************************************************
  // ADDED 1-30-2009 BY ISAAC ABBOTT
  CreatePipelineMacro( itk::RGBPixel< unsigned char > );
  CreatePipelineMacro( itk::RGBPixel< char > );
  CreatePipelineMacro( itk::RGBPixel< unsigned short > );
  CreatePipelineMacro( itk::RGBPixel< short > );
  CreatePipelineMacro( itk::RGBPixel< unsigned int > );
  CreatePipelineMacro( itk::RGBPixel< int > );
  CreatePipelineMacro( itk::RGBPixel< unsigned long > );
  CreatePipelineMacro( itk::RGBPixel< long > );
  CreatePipelineMacro( itk::RGBPixel< float > );
  CreatePipelineMacro( itk::RGBPixel< double > );
  //***********************************************************

  this->Importer->Update();
}


//----------------------------------------------------------------------------
vtkImageData * vtkKWImage::GetVTKImage()
{
  return this->Importer->GetOutput();
}

//----------------------------------------------------------------------------
const vtkKWImage::ImageBaseType * vtkKWImage::GetITKImageBase() const
{
  return this->ItkImage;
}
/******************************************************************************
// REMOVED 1-30-2009 BY ISAAC ABBOTT. 
//----------------------------------------------------------------------------
vtkKWImage::ITKScalarPixelType vtkKWImage::GetITKScalarPixelType() const
{
  ITKScalarPixelType pixelType = itk::ImageIOBase::UCHAR;

  ImageBaseType * itkImageBase = this->ItkImage.GetPointer();

  if( dynamic_cast< itk::Image< unsigned char, 3> * >( itkImageBase ) )
    {
    pixelType = itk::ImageIOBase::UCHAR;
    }
  else if( dynamic_cast< itk::Image< char, 3> * >( itkImageBase ) )
    {
    pixelType = itk::ImageIOBase::CHAR;
    }
  else if( dynamic_cast< itk::Image< short, 3> * >( itkImageBase ) )
    {
    pixelType = itk::ImageIOBase::SHORT;
    }
  else if( dynamic_cast< itk::Image< unsigned short, 3> * >( itkImageBase ) )
    {
    pixelType = itk::ImageIOBase::USHORT;
    }
  else if( dynamic_cast< itk::Image< int, 3> * >( itkImageBase ) )
    {
    pixelType = itk::ImageIOBase::INT;
    }
  else if( dynamic_cast< itk::Image< unsigned int, 3> * >( itkImageBase ) )
    {
    pixelType = itk::ImageIOBase::UINT;
    }
  else if( dynamic_cast< itk::Image< long, 3> * >( itkImageBase ) )
    {
    pixelType = itk::ImageIOBase::LONG;
    }
  else if( dynamic_cast< itk::Image< unsigned long, 3> * >( itkImageBase ) )
    {
    pixelType = itk::ImageIOBase::ULONG;
    }
  else if( dynamic_cast< itk::Image< float, 3> * >( itkImageBase ) )
    {
    pixelType = itk::ImageIOBase::FLOAT;
    }
  else if( dynamic_cast< itk::Image< double, 3> * >( itkImageBase ) )
    {
    pixelType = itk::ImageIOBase::DOUBLE;
    }

  return pixelType;
}
************************************************************************************/
//**********************************************************************************
//ADDED 1-30-2009 BY ISAAC ABBOTT
vtkKWImage::ITKComponentType vtkKWImage::GetITKComponentType() const
{
  ITKComponentType componentType = itk::ImageIOBase::UCHAR;

  ImageBaseType * itkImageBase = this->ItkImage.GetPointer();

  if( dynamic_cast< itk::Image< unsigned char, 3> * >( itkImageBase ) ||
	  dynamic_cast< itk::Image< itk::RGBPixel<unsigned char>,  3 > * >( itkImageBase ) )
    {
    componentType = itk::ImageIOBase::UCHAR;
    }
  else if( dynamic_cast< itk::Image< char, 3> * >( itkImageBase ) ||
	  dynamic_cast< itk::Image< itk::RGBPixel<char>,  3 > * >( itkImageBase ) )
    {
    componentType = itk::ImageIOBase::CHAR;
    }
  else if( dynamic_cast< itk::Image< short, 3> * >( itkImageBase ) ||
	  dynamic_cast< itk::Image< itk::RGBPixel<short>,  3 > * >( itkImageBase ) )
    {
    componentType = itk::ImageIOBase::SHORT;
    }
  else if( dynamic_cast< itk::Image< unsigned short, 3> * >( itkImageBase ) ||
	  dynamic_cast< itk::Image< itk::RGBPixel<unsigned short>,  3 > * >( itkImageBase ) )
    {
    componentType = itk::ImageIOBase::USHORT;
    }
  else if( dynamic_cast< itk::Image< int, 3> * >( itkImageBase ) ||
	  dynamic_cast< itk::Image< itk::RGBPixel<int>,  3 > * >( itkImageBase ) )
    {
    componentType = itk::ImageIOBase::INT;
    }
  else if( dynamic_cast< itk::Image< unsigned int, 3> * >( itkImageBase ) ||
	  dynamic_cast< itk::Image< itk::RGBPixel<unsigned int>,  3 > * >( itkImageBase ) )
    {
    componentType = itk::ImageIOBase::UINT;
    }
  else if( dynamic_cast< itk::Image< long, 3> * >( itkImageBase ) ||
	  dynamic_cast< itk::Image< itk::RGBPixel<long>,  3 > * >( itkImageBase ) )
    {
    componentType = itk::ImageIOBase::LONG;
    }
  else if( dynamic_cast< itk::Image< unsigned long, 3> * >( itkImageBase ) ||
	  dynamic_cast< itk::Image< itk::RGBPixel<unsigned long>,  3 > * >( itkImageBase ) )
    {
    componentType = itk::ImageIOBase::ULONG;
    }
  else if( dynamic_cast< itk::Image< float, 3> * >( itkImageBase ) ||
	  dynamic_cast< itk::Image< itk::RGBPixel<float>,  3 > * >( itkImageBase ) )
    {
    componentType = itk::ImageIOBase::FLOAT;
    }
  else if( dynamic_cast< itk::Image< double, 3> * >( itkImageBase ) ||
	  dynamic_cast< itk::Image< itk::RGBPixel<double>,  3 > * >( itkImageBase ) )
    {
    componentType = itk::ImageIOBase::DOUBLE;
    }

  return componentType;
}

vtkKWImage::ITKPixelType vtkKWImage::GetITKPixelType() const
{
  ITKPixelType pixelType = itk::ImageIOBase::SCALAR;

  ImageBaseType * itkImageBase = this->ItkImage.GetPointer();

  if( dynamic_cast< itk::Image< unsigned char, 3> * >( itkImageBase ) ||
	  dynamic_cast< itk::Image< char, 3> * >( itkImageBase ) ||
	  dynamic_cast< itk::Image< unsigned short, 3> * >( itkImageBase ) ||
	  dynamic_cast< itk::Image< short, 3> * >( itkImageBase ) ||
	  dynamic_cast< itk::Image< unsigned int, 3> * >( itkImageBase ) ||
	  dynamic_cast< itk::Image< int, 3> * >( itkImageBase ) ||
	  dynamic_cast< itk::Image< unsigned long, 3> * >( itkImageBase ) ||
	  dynamic_cast< itk::Image< long, 3> * >( itkImageBase ) ||
	  dynamic_cast< itk::Image< float, 3> * >( itkImageBase ) ||
	  dynamic_cast< itk::Image< double, 3> * >( itkImageBase ) )
    {
    pixelType = itk::ImageIOBase::SCALAR;
    }
  else if( dynamic_cast< itk::Image< itk::RGBPixel<unsigned char>,  3 > * >( itkImageBase ) ||
	  dynamic_cast< itk::Image< itk::RGBPixel<char>,  3 > * >( itkImageBase ) ||
	  dynamic_cast< itk::Image< itk::RGBPixel<unsigned short>,  3 > * >( itkImageBase ) ||
	  dynamic_cast< itk::Image< itk::RGBPixel<short>,  3 > * >( itkImageBase ) ||
	  dynamic_cast< itk::Image< itk::RGBPixel<unsigned int>,  3 > * >( itkImageBase ) ||
	  dynamic_cast< itk::Image< itk::RGBPixel<int>,  3 > * >( itkImageBase ) ||
	  dynamic_cast< itk::Image< itk::RGBPixel<unsigned long>,  3 > * >( itkImageBase ) ||
	  dynamic_cast< itk::Image< itk::RGBPixel<long>,  3 > * >( itkImageBase ) ||
	  dynamic_cast< itk::Image< itk::RGBPixel<float>,  3 > * >( itkImageBase ) ||
	  dynamic_cast< itk::Image< itk::RGBPixel<double>,  3 > * >( itkImageBase ) )
    {
    pixelType = itk::ImageIOBase::RGB;
    }
  
  return pixelType;
}
//**********************************************************************************
//**********************************************************************************


/*********************************************************************************
// REMOVED 1-30-2009 BY ISAAC ABBOTT
//----------------------------------------------------------------------------
int vtkKWImage::GetVTKScalarPixelType()
{
  return this->GetVTKImage()->GetScalarType();
}
**********************************************************************************/
/**********************************************************************************
// ADDED 1-30-2009 BY ISAAC ABBOTT												 */
int vtkKWImage::GetVTKPixelType()
{
	return this->GetVTKImage()->GetNumberOfScalarComponents();
}
int vtkKWImage::GetVTKComponentType()
{
	return this->GetVTKImage()->GetScalarType();
}
//**********************************************************************************