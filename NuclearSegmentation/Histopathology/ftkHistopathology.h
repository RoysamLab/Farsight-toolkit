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
#ifndef __ftkHistopathology_h
#define __ftkHistopathology_h

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "ftkImage/ftkImage.h"
#include "tinyxml/tinyxml.h"
#include "ftkCommon/ftkUtils.h"
#include "NuclearSegmentation/ftkNuclearSegmentation.h"
#include "CytoplasmSegmentation/CytoplasmSegmentation.h"

namespace ftk
{ 

class Operation;
class Input;

/** \class Histopathology
 *  \brief For storage of a complete Histopathology segmentation 
 *   
 *  Handles the execution, result, and editing of a Histopathology segmentation
 *  
 */
class Histopathology
{
public:
	Histopathology();
	~Histopathology();

	bool LoadProject(std::string xmlfname);		//Load from the xml file containing the project definition
	bool ProcessInputs();						//Process the inputs from the beginning
	bool LoadImages();
	bool LoadAll(std::string filename);			//Complete Restore from files

	//Misc string Gets
	std::string GetErrorMessage() { return errorMessage; };

	//Get Data:
	ftk::Image::Pointer GetDataImage(void){ return dataImage; };	
	ftk::Image::Pointer GetLabelImage(void){ return labelImage; };
	//*********************************************************************************************

protected:
	Input parseInputElement(TiXmlElement *objectElement);
	Operation parseOperationElement(TiXmlElement *operationElement);

	bool RunNuclearSegmentation(Input in);
	bool RunCytoplasmSegmentation(Input in);

	std::string xmlFilename;				//The project file (full path)
	std::vector<Input> input;				//The input information defining the project

	std::vector<std::string> dataFilename;	//the filename of the data image		(full path)
	std::vector<std::string> labelFilename;	//the filename of the label image		(full path)
	std::string paramFilename;				//the filename of the parameter file	(full path)
	std::string featureFilename;			//the filename of the feature txt file	(full path)
	std::string headerFilename;				//the filename of the feature names txt (full path)
	std::string editFilename;				//the filename of the edit record txt	(full path)
	std::string wrkDir;						//Working Directory

	std::string errorMessage;

	ftk::Image::Pointer dataImage;		//The data image
	ftk::Image::Pointer labelImage;		//My label image

}; // end Histopathology

class Operation
{
public:
	Operation(){done = false;};
	std::string type;
	std::string infile;
	std::string outfile;
	bool done;
};

class Input
{
public:
	Input(){channel = 0;};
	std::string filename;
	int channel;
	std::string chname;
	unsigned char color[3];
	std::string type;
	Operation operation;
};



}  // end namespace ftk

#endif	// end __ftkHistopathology_h

