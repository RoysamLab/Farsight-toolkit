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
	cytchannelNumber  = -1;				//Use this channel from the dataImage for segmentation
	memchannelNumber = -1;				//Use this channel from the dataImage for segmentation
	nucLabelFilename.clear();			//Name of the file that is the nuclear label image
	labelImage = NULL;					//my label image (I will append the cytoplasm segmentation)
	nucChannelNumber = -1;
	cytoChannelNumber = -1;
	cytoLabelFilename.clear();			//The label image of cytoplasm channel

	paramNames.push_back("draw_real_boundaries");				//0
	paramNames.push_back("remove_stromal_cell_boundaries");		//1
	paramNames.push_back("draw_synthetic_boundaries");			//2
	paramNames.push_back("radius_of_synthetic_boundaries");		//3
	paramNames.push_back("number_of_levels");					//4
	paramNames.push_back("number_of_levels_in_foreground");		//5

	//Defaults
	params[3] = 12;
	params[4] = 1;
	params[5] = 1;

}

void CytoplasmSegmentation::SetParameter(std::string name, int value)
{
	for( int i=0; i<(int)paramNames.size(); ++i )
		if( paramNames.at(i) == name )
			params[i] = value;
}

int CytoplasmSegmentation::GetParameter(std::string name)
{
	for(int i=0; i<(int)paramNames.size(); ++i)
		if(paramNames.at(i) == name)
			return params[i];
	return -1;
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
bool CytoplasmSegmentation::SetDataInput(ftk::Image::Pointer inImg, std::string fname, int cytNumber , int memNumber )
{
	if(cytNumber > inImg->GetImageInfo()->numChannels)
		return false;
	if(memNumber > inImg->GetImageInfo()->numChannels)
		return false;
	if(inImg->GetImageInfo()->numZSlices != 1)
		return false;
	if(inImg->GetImageInfo()->dataType != itk::ImageIOBase::UCHAR)
		return false;

	this->dataFilename = fname;
	this->dataImage = inImg;
	this->cytchannelNumber = cytNumber;
	this->memchannelNumber = memNumber;
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

	WholeCellSeg *whole_cell = new WholeCellSeg;

	/*List from constructor
	paramNames.push_back("draw_real_boundaries");				//0
	paramNames.push_back("remove_stromal_cell_boundaries");		//1
	paramNames.push_back("draw_synthetic_boundaries");			//2
	paramNames.push_back("radius_of_synthetic_boundaries");		//3
	paramNames.push_back("number_of_levels");					//4
	paramNames.push_back("number_of_levels_in_foreground");		//5*/

	if( cytchannelNumber == -1 ){
		params[0] = 0;
		params[2] = 1;
	}

	if( cytchannelNumber > -1 ){
		UCharImageType3D::Pointer data = dataImage->GetItkPtr<unsigned char>(0,cytchannelNumber);
		DataExtractType::Pointer deFilter = DataExtractType::New();
		UCharImageType3D::RegionType dRegion = data->GetLargestPossibleRegion();
		dRegion.SetSize(2,0);
		deFilter->SetExtractionRegion(dRegion);
		deFilter->SetDirectionCollapseToIdentity();
		deFilter->SetInput( data );

		DataCastType::Pointer dFilter = DataCastType::New();
		dFilter->SetInput( deFilter->GetOutput() );

		try{
			dFilter->Update();
		}
		catch( itk::ExceptionObject & excep )
		{
				std::cerr << "Exception caught !" << std::endl;
			std::cerr << excep << std::endl;
			return false;
		}
		whole_cell->set_cyt_img( dFilter->GetOutput() );
	}

	if( memchannelNumber > -1 ){
		UCharImageType3D::Pointer data = dataImage->GetItkPtr<unsigned char>(0,memchannelNumber);
		DataExtractType::Pointer deFilter = DataExtractType::New();
		UCharImageType3D::RegionType dRegion = data->GetLargestPossibleRegion();
		dRegion.SetSize(2,0);
		deFilter->SetExtractionRegion(dRegion);
		deFilter->SetDirectionCollapseToIdentity();
		deFilter->SetInput( data );

		DataCastType::Pointer dFilter = DataCastType::New();
		dFilter->SetInput( deFilter->GetOutput() );

		try{
			dFilter->Update();
		}
		catch( itk::ExceptionObject & excep )
		{
				std::cerr << "Exception caught !" << std::endl;
			std::cerr << excep << std::endl;
			return false;
		}
		whole_cell->set_mem_img( dFilter->GetOutput() );
	}

	UShortImageType3D::Pointer nuc = labelImage->GetItkPtr<unsigned short>(0,nucChannelNumber);

	LabelExtractType::Pointer leFilter = LabelExtractType::New();
	UShortImageType3D::RegionType lRegion = nuc->GetLargestPossibleRegion();
	lRegion.SetSize(2,0);
	leFilter->SetExtractionRegion(lRegion);
	leFilter->SetDirectionCollapseToIdentity();
	leFilter->SetInput( nuc );

	LabelCastType::Pointer lFilter = LabelCastType::New();
	lFilter->SetInput( leFilter->GetOutput() );
	try
	{
		lFilter->Update();
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
		return false;
	}

	whole_cell->set_nuc_img( lFilter->GetOutput() );
	whole_cell->set_parameters(params);
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
