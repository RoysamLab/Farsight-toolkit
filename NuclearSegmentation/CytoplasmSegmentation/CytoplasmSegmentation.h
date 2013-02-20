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
#ifndef __CytoplasmSegmentation_h
#define __CytoplasmSegmentation_h

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <iostream>
#include "whole_cell.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"
#include "itkExtractImageFilter.h"

#include "ftkImage/ftkImage.h"

namespace ftk
{

class CytoplasmSegmentation
{
public:
	CytoplasmSegmentation();

	bool LoadInputs(std::string datafname, std::string nucfname, int nucDataChNumber = 0, int nucLabChNumber = 0);	//Load the input image from a file
	bool SetDataInput(ftk::Image::Pointer inImg, std::string fname, int cytNumber = -1, int memNumber = -1 );	//Pass a pointer to the already loaded image
	std::vector<std::string> GetParameterNames(){ return paramNames; };
	void SetParameter(std::string name, int value);
	bool SetNucleiInput(ftk::Image::Pointer lbImg, std::string fname, int chNumber = 0);
	bool Run();
	int GetParameter(std::string name);

	bool SaveOutputImage(std::string fname = "");			//Save the output image

private:
	std::string dataFilename;				//Name of the file that the data came from
	ftk::Image::Pointer dataImage;			//The data image
	int cytchannelNumber,memchannelNumber;	//Use this channel from the dataImage for segmentation
	std::string nucLabelFilename;			//Name of the file that is the nuclear label image
	ftk::Image::Pointer labelImage;			//my label image (I will append the cyto label image)
	int nucChannelNumber;					//Input Nucleus Labeled Image
	int cytoChannelNumber;					//Output Labeled Image
	std::string cytoLabelFilename;			//The label image of cytoplasm channel
	std::vector<std::string> paramNames;
	int params[6];
};

}  // end namespace ftk

#endif