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

namespace ftk 
{

CytoplasmSegmentation::CytoplasmSegmentation()
{
	dataFilename.clear();				//Name of the file that the data came from
	dataImage = NULL;					//The data image
	channelNumber = 0;					//Use this channel from the dataImage for segmentation
	nucLabelFilename.clear();			//Name of the file that is the nuclear label image
	labelImage = NULL;					//my label image (I will append the cytoplasm segmentation)
	nucChannelNumber = 0;
	cytoChannelNumber = 0;
	cytoLabelFilename.clear();			//The label image of cytoplasm channel
}

//Load the input images from files
bool CytoplasmSegmentation::LoadInputs(std::string datafname, std::string nucfname, int nucDataChNumber, int nucLabChNumber )
{
	ftk::Image::Pointer tmpDImg = ftk::Image::New();
	if(!tmpDImg->LoadFile(datafname))	//Load for display
		return false;
	if(!this->SetDataInput(tmpDImg,datafname,nucDataChNumber))
		return false;

	ftk::Image::Pointer tmpLImg = ftk::Image::New();
	if(!tmpLImg->LoadFile(nucfname))
		return false;
	if(!this->SetNucleiInput(tmpLImg, nucfname, nucLabChNumber))
		return false;

	return true;
}

//Pass a pointer to the already loaded data image
bool CytoplasmSegmentation::SetDataInput(ftk::Image::Pointer inImg, std::string fname, int chNumber)
{
	if(chNumber > inImg->GetImageInfo()->numChannels)
		return false;
	if(inImg->GetImageInfo()->numZSlices != 1)
		return false;
	if(inImg->GetImageInfo()->dataType != itk::ImageIOBase::UCHAR)
		return false;

	this->dataFilename = fname;
	this->dataImage = inImg;
	this->channelNumber = chNumber;
	return true;

}

//Pass a pointer to the already loaded label image
bool CytoplasmSegmentation::SetNucleiInput(ftk::Image::Pointer lbImg, std::string fname, int chNumber)
{
	if(chNumber > lbImg->GetImageInfo()->numChannels)
		return false;

	if(lbImg->GetImageInfo()->numZSlices != 1)
		return false;
	if(lbImg->GetImageInfo()->dataType != itk::ImageIOBase::USHORT)
		return false;

	this->nucLabelFilename = fname;
	this->labelImage = lbImg;
	this->nucChannelNumber = chNumber;
	return true;
}

bool CytoplasmSegmentation::Run()
{
	typedef itk::Image< unsigned short, 3 > UShortImageType3D;
	typedef itk::Image< unsigned char, 3 > UCharImageType3D;
	typedef itk::Image< unsigned short, 2 > UShortImageType2D;
	typedef itk::Image< unsigned char, 2 > UCharImageType2D;
	typedef itk::ExtractImageFilter< UCharImageType3D, UCharImageType2D > DataExtractType;
	typedef itk::ExtractImageFilter< UShortImageType3D, UShortImageType2D > LabelExtractType; 
	typedef itk::CastImageFilter< UCharImageType2D, UShortImageType2D > DataCastType;
	typedef itk::CastImageFilter< UShortImageType2D, UShortImageType2D > LabelCastType;

	if(!dataImage)
		return false;
	if(!labelImage)
		return false;

	UCharImageType3D::Pointer data = dataImage->GetItkPtr<unsigned char>(0,channelNumber);
	UShortImageType3D::Pointer nuc = labelImage->GetItkPtr<unsigned short>(0,nucChannelNumber);

	DataExtractType::Pointer deFilter = DataExtractType::New();
	UCharImageType3D::RegionType dRegion = data->GetLargestPossibleRegion();
	dRegion.SetSize(2,0);
	deFilter->SetExtractionRegion(dRegion);
	deFilter->SetInput( data );

	LabelExtractType::Pointer leFilter = LabelExtractType::New();
	UShortImageType3D::RegionType lRegion = nuc->GetLargestPossibleRegion();
	lRegion.SetSize(2,0);
	leFilter->SetExtractionRegion(lRegion);
	leFilter->SetInput( nuc );

	DataCastType::Pointer dFilter = DataCastType::New();
	dFilter->SetInput( deFilter->GetOutput() );
	LabelCastType::Pointer lFilter = LabelCastType::New();
	lFilter->SetInput( leFilter->GetOutput() );
	try
	{
		dFilter->Update();
		lFilter->Update();
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
		return false;
	}

	WholeCellSeg *whole_cell = new WholeCellSeg;
	whole_cell->set_nuc_img( lFilter->GetOutput() );
	whole_cell->set_cyt_img( dFilter->GetOutput() );
	whole_cell->RunBinarization();
	whole_cell->RunSegmentation();

	UShortImageType2D::Pointer cyto = whole_cell->getSegPointer();

	ftk::Image::DataType dType = labelImage->GetImageInfo()->dataType;
	unsigned char bpPix = labelImage->GetImageInfo()->bytesPerPix;
	unsigned short cs = labelImage->GetImageInfo()->numColumns;
	unsigned short rs = labelImage->GetImageInfo()->numRows;
	unsigned short zs = 1;
	std::string name = "cyto";
	std::vector<unsigned char> color(3,255);
	labelImage->AppendChannelFromData3D(cyto->GetBufferPointer(), dType, bpPix, cs, rs, zs, name, color, true);
	cytoChannelNumber = labelImage->GetImageInfo()->numChannels - 1;

	delete whole_cell;
	return true;
}

bool CytoplasmSegmentation::SaveOutputImage(std::string fname)
{
	if(!labelImage)
		return false;

	if(fname.size() == 0)
	{
		size_t pos = dataFilename.find_last_of(".");
		std::string base = dataFilename.substr(0,pos);
		std::string tag = "_cyto";
		std::string ext = "tif";
		cytoLabelFilename = base + tag + "." + ext;
	}

	size_t pos = cytoLabelFilename.find_last_of(".");
	std::string base = cytoLabelFilename.substr(0,pos);
	std::string ext = cytoLabelFilename.substr(pos+1);

	if(!labelImage->SaveChannelAs(cytoChannelNumber, base, ext ))
		return false;

	return true;
}

//*******************************************************************************************************
//*******************************************************************************************************
//*******************************************************************************************************
//*******************************************************************************************************
//*******************************************************************************************************
} //END NAMESPACE FTK
