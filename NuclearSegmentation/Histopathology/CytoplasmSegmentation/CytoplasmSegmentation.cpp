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

  Program:   Farsight Biological Image Segmentation and Visualization Toolkit
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision: 0.00 $

=========================================================================*/
#include "CytoplasmSegmentation.h"

bool CytoplasmSegmentation::Run(std::string cytoFname, std::string nucLabelFname, std::string outLabelFname)
{
	typedef itk::Image< unsigned short, 2 > UShortImageType;
	typedef itk::ImageFileReader< UShortImageType >  ReaderType;

	ReaderType::Pointer cytoReader = ReaderType::New();
	cytoReader->SetFileName( cytoFname.c_str() );
	try
	{
		cytoReader->Update();
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
		return false;
	}

	ReaderType::Pointer nucReader = ReaderType::New();
	nucReader->SetFileName( nucLabelFname.c_str() );
	try
	{
		nucReader->Update();
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
		return false;
	}

	WholeCellSeg *whole_cell = new WholeCellSeg;
	whole_cell->set_nuc_img( nucReader->GetOutput() );
	whole_cell->set_cyt_img( cytoReader->GetOutput() );
	whole_cell->RunBinarization();
	whole_cell->RunSegmentation();

	typedef itk::ImageFileWriter < UShortImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName( outLabelFname.c_str() );
	writer->SetInput( whole_cell->getSegPointer() );
	try
	{
		writer->Update();
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
		return false;
	}

	delete whole_cell;
	return true;
}