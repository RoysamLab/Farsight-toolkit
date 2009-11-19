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

namespace ftk 
{

//Constructor
Histopathology::Histopathology()
{
	dataImage = NULL;
	labelImage = NULL;
}

Histopathology::~Histopathology()
{
	dataImage = NULL;
	labelImage = NULL;
}

//*************************************************************************
//
// The project definition is found in an XML file.  This function loads up
// all of the information from this XML file.  It does not actually load
// any data.
//
//*************************************************************************
bool Histopathology::LoadProject(std::string xmlfname)
{
	xmlFilename = xmlfname;
	input.clear();

	TiXmlDocument doc;
	if ( !doc.LoadFile( xmlFilename.c_str() ) )
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

	//Parents we know of: input
	TiXmlElement* parentElement = rootElement->FirstChildElement();
	while (parentElement)
	{
		const char * parent = parentElement->Value();

		if ( strcmp( parent, "input" ) == 0 )	//Get the new input:
		{
			input.push_back( parseInputElement(parentElement) );
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

	return true;
}

Input Histopathology::parseInputElement(TiXmlElement *inputElement)
{
	Input input;

	if ( strcmp( inputElement->Value(), "input" ) != 0 )
		return input;

	TiXmlElement *member = inputElement->FirstChildElement();
	while(member)
	{
		const char* memberName = member->Value();
		if ( strcmp( memberName, "filename" ) == 0 )
		{
			input.filename = member->GetText();
		}
		else if ( strcmp( memberName, "channel" ) == 0 )
		{
			input.channel = atoi(member->GetText());
		}
		else if ( strcmp( memberName, "name" ) == 0 )
		{
			input.chname = member->GetText();
		}
		else if ( strcmp( memberName, "color" ) == 0 )
		{
			input.color[0] = atoi(member->Attribute("r"));
			input.color[1] = atoi(member->Attribute("g"));
			input.color[2] = atoi(member->Attribute("b"));
		}
		else if ( strcmp( memberName, "type") == 0 )
		{
			input.type = member->GetText();
		}
		else if ( strcmp( memberName, "operation") == 0 )
		{
			input.operation = parseOperationElement(member);
		}
		member = member->NextSiblingElement();
	}
	return input;
}

Operation Histopathology::parseOperationElement(TiXmlElement *operationElement)
{
	Operation op;
	op.done = false;
	op.type = "none";
	op.outfile = "";

	if ( strcmp( operationElement->Value(), "operation" ) != 0 )
		return op;

	op.type = operationElement->Attribute("type");

	TiXmlElement *member = operationElement->FirstChildElement();
	while(member)
	{
		const char* memberName = member->Value();
		if ( strcmp( memberName, "resultfile" ) == 0 )
		{
			op.outfile = member->GetText();
		}
		if( strcmp( memberName, "nucfile" ) == 0 )
		{
			op.infile = member->GetText();
		}
		member = member->NextSiblingElement();
	}

	if(FileExists(op.outfile))
		op.done = true;

	return op;
}

bool Histopathology::ProcessInputs(void)
{
	//First do the nuclear segmentations:
	for(int i=0; i<input.size(); ++i)
	{
		if( input.at(i).operation.type == "NucSeg" && input.at(i).operation.done == false )
		{
			if(!RunNuclearSegmentation(input.at(i)))
				return false;
		}
	}
	//Then do the Cytoplasm Segmentations:
	for(int i=0; i<input.size(); ++i)
	{
		if( input.at(i).operation.type == "CytoSeg" && input.at(i).operation.done == false )
		{
			if(!RunCytoplasmSegmentation(input.at(i)))
				return false;
		}
	}
	return true;
}

bool Histopathology::RunNuclearSegmentation(Input in)
{
	ftk::NuclearSegmentation * seg = new ftk::NuclearSegmentation();
	seg->SetInputs(in.filename,"");
	seg->SetChannel(in.channel);
	seg->LoadData();
	seg->Binarize();
	seg->DetectSeeds(false);
	seg->RunClustering();
	seg->Finalize();
	//seg->ComputeFeatures();
	//seg->ReleaseSegMemory();
	seg->SaveResultImage();
	in.operation.outfile = seg->GetLabelFilename();
	in.operation.done = true;

	delete seg;
	return true;
}

bool Histopathology::RunCytoplasmSegmentation(Input in)
{
	CytoplasmSegmentation * seg = new CytoplasmSegmentation();
	seg->Run(in.filename,in.operation.infile,in.operation.outfile);
	in.operation.done = true;

	delete seg;
	return true;
}

//This function loads the data images and the result (label) images into an ftk::Image 
bool Histopathology::LoadImages()
{
	std::vector<std::string> dataChName;
	std::vector<unsigned char> dataColor;
	std::vector<std::string> lablChName;
	std::vector<unsigned char> lablColor;

	//First do the nuclear segmentations:
	for(int i=0; i<input.size(); ++i)
	{
		if( input.at(i).operation.type == "NucSeg" || input.at(i).operation.type == "CytoSeg")
		{
			dataChName.push_back( input.at(i).chname );
			dataColor.push_back(input.at(i).color[0]);
			dataColor.push_back(input.at(i).color[1]);
			dataColor.push_back(input.at(i).color[2]);
			dataFilename.push_back( input.at(i).filename );
			if( input.at(i).operation.done == true )
			{
				lablChName.push_back( input.at(i).chname );
				lablColor.push_back(input.at(i).color[0]);
				lablColor.push_back(input.at(i).color[1]);
				lablColor.push_back(input.at(i).color[2]);
				labelFilename.push_back( input.at(i).operation.outfile );
			}
		}
		else
		{
			dataChName.push_back( input.at(i).chname );
			dataColor.push_back(input.at(i).color[0]);
			dataColor.push_back(input.at(i).color[1]);
			dataColor.push_back(input.at(i).color[2]);
			dataFilename.push_back( input.at(i).filename );
		}
	}

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