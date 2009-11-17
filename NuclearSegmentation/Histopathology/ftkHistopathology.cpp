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
#include "ftkHistopathology.h"
#include <ctime>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>

namespace ftk 
{

//Constructor
Histopathology::Histopathology()
{
}

Histopathology::~Histopathology()
{
}

//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************
// LOAD UP FORMER RESULTS FROM FILES:
//*******************************************************************************************
bool Histopathology::LoadAll(std::string filename)
{
	TiXmlDocument doc;
	if ( !doc.LoadFile( filename.c_str() ) )
	{
		errorMessage = "Unable to load XML File";
		return 0;
	}

	TiXmlElement* rootElement = doc.FirstChildElement();
	const char* docname = rootElement->Value();
	if ( strcmp( docname, "Histopathology" ) != 0 )
	{
		errorMessage = "Incorrect XML root Element: ";
		errorMessage.append(rootElement->Value());
		return 0;
	}

	std::vector<std::string> dataChName;
	std::vector<unsigned char> dataColor;
	std::vector<std::string> lablChName;
	std::vector<unsigned char> lablColor;

	//Parents we know of: datafilename,resultfilename,object,parameter
	TiXmlElement* parentElement = rootElement->FirstChildElement();
	while (parentElement)
	{
		const char * parent = parentElement->Value();

		if ( strcmp( parent, "datafile" ) == 0 )
		{
			dataChName.push_back( parentElement->Attribute("chname") );
			dataColor.push_back( atoi(parentElement->Attribute("r")) );
			dataColor.push_back( atoi(parentElement->Attribute("g")) );
			dataColor.push_back( atoi(parentElement->Attribute("b")) );
			dataFilename.push_back( parentElement->GetText() );
		}
		else if ( strcmp( parent, "resultfile" ) == 0 )
		{
			lablChName.push_back( parentElement->Attribute("chname") );
			lablColor.push_back( atoi(parentElement->Attribute("r")) );
			lablColor.push_back( atoi(parentElement->Attribute("g")) );
			lablColor.push_back( atoi(parentElement->Attribute("b")) );
			labelFilename.push_back( parentElement->GetText() );
		}
		else if ( strcmp( parent, "featurefile" ) == 0 )
		{
			featureFilename = parentElement->GetText();
		}
		else if ( strcmp( parent, "headerfile" ) == 0 )
		{
			headerFilename = parentElement->GetText();
		}
		else if ( strcmp( parent, "paramfile" ) == 0 )
		{
			paramFilename = parentElement->GetText();
		}
		else if ( strcmp( parent, "editfile" ) == 0 )
		{
			editFilename = parentElement->GetText();
		}
		else
		{
			errorMessage = "Unrecognized parent element: ";
			errorMessage.append(parent);
			return 0;
		}
		parentElement = parentElement->NextSiblingElement();
	} // end while(parentElement)

	//doc.close();

	dataImage = ftk::Image::New();
	if(!dataImage->LoadFilesAsMultipleChannels(dataFilename,dataChName,dataColor))	//Load for display
	{
		errorMessage = "Data Image failed to load";
		dataImage = 0;
		return false;
	}

	labelImage = ftk::Image::New();
	if(!labelImage->LoadFilesAsMultipleChannels(labelFilename,lablChName,lablColor))	//Load for display
	{
		errorMessage = "Label Image failed to load";
		labelImage = 0;
		return false;
	}

	return true;
}

} //END NAMESPACE FTK