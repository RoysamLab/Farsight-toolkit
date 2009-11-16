/* 
 * Copyright 2009 Rensselaer Polytechnic Institute
 * This program is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
 * for more details.
 * 
 * You should have received a copy of the GNU General Public License along 
 * with this program; if not, write to the Free Software Foundation, Inc., 
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkImageToVTKImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2005-06-12 01:23:44 $
  Version:   $Revision: 1.7 $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkImageToVTKImageFilter_txx
#define _itkImageToVTKImageFilter_txx

#include "itkImageToVTKImageFilter.h"

namespace itk
{



/**
 * Constructor
 */
template <class TInputImage>
ImageToVTKImageFilter<TInputImage>
::ImageToVTKImageFilter()
{

  m_Importer = vtkImageImport::New();

  m_Exporter = ExporterFilterType::New();

  m_Importer->SetUpdateInformationCallback(m_Exporter->GetUpdateInformationCallback());
  m_Importer->SetPipelineModifiedCallback(m_Exporter->GetPipelineModifiedCallback());
  m_Importer->SetWholeExtentCallback(m_Exporter->GetWholeExtentCallback());
  m_Importer->SetSpacingCallback(m_Exporter->GetSpacingCallback());
  m_Importer->SetOriginCallback(m_Exporter->GetOriginCallback());
  m_Importer->SetScalarTypeCallback(m_Exporter->GetScalarTypeCallback());
  m_Importer->SetNumberOfComponentsCallback(m_Exporter->GetNumberOfComponentsCallback());
  m_Importer->SetPropagateUpdateExtentCallback(m_Exporter->GetPropagateUpdateExtentCallback());
  m_Importer->SetUpdateDataCallback(m_Exporter->GetUpdateDataCallback());
  m_Importer->SetDataExtentCallback(m_Exporter->GetDataExtentCallback());
  m_Importer->SetBufferPointerCallback(m_Exporter->GetBufferPointerCallback());
  m_Importer->SetCallbackUserData(m_Exporter->GetCallbackUserData());

}




/**
 * Destructor
 */
template <class TInputImage>
ImageToVTKImageFilter<TInputImage>
::~ImageToVTKImageFilter()
{
  if( m_Importer )
    {
    m_Importer->Delete();
    m_Importer = 0;
    }
}



/**
 * Set an itk::Image as input 
 */
template <class TInputImage>
void
ImageToVTKImageFilter<TInputImage>
::SetInput( const InputImageType * inputImage )
{
  m_Exporter->SetInput( inputImage );
}



/**
 * Get a vtkImage as output
 */
template <class TInputImage>
vtkImageData *
ImageToVTKImageFilter<TInputImage>
::GetOutput() const
{
  return m_Importer->GetOutput();
}




/**
 * Get the importer filter
 */
template <class TInputImage>
vtkImageImport *
ImageToVTKImageFilter<TInputImage>
::GetImporter() const
{
  return m_Importer;
}



/**
 * Get the exporter filter
 */
template <class TInputImage>
typename ImageToVTKImageFilter<TInputImage>::ExporterFilterType *
ImageToVTKImageFilter<TInputImage>
::GetExporter() const
{
  return m_Exporter.GetPointer();
}



/**
 * Delegate the Update to the importer
 */
template <class TInputImage>
void
ImageToVTKImageFilter<TInputImage>
::Update()
{
  m_Importer->Update();
}





} // end namespace itk

#endif

