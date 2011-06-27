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
#ifndef __ftkProjectFiles_h
#define __ftkProjectFiles_h

#include <tinyxml/tinyxml.h>
#include <ftkCommon/ftkUtils.h>

#include <vector>
#include <string>
#include <cstdio>
#include <sstream>
#include <iomanip>

namespace ftk
{

class ProjectFiles
{
public:
	ProjectFiles();
	//FUNCTIONS:
	bool Read(std::string filename);
	bool Write();
	std::string GetFullName();
	std::string GetFullInput();
	std::string GetFullOutput();
	std::string GetFullLog();
	std::string GetFullDef();
	std::string GetFullTable();
	std::string GetFullAdjTables();
	void ClearAll();

	//VARIABLES
	std::string path;		//Includes last separator character
	std::string name;		//Name of the project (also used as filename)
	std::string type;		// If it is a single image project -> "single"
							// or multiImage/Time Series project -> "multi"

	bool nucSegValidated;	//Nuclear Segmentation is validated

	//Filenames:
	std::string input;
	std::string output;
	std::string log;
	std::string definition;
	std::string table;
	std::string adjTables;

	bool inputSaved;
	bool outputSaved;
	//bool logSaved; //Log is continuously saved
	bool definitionSaved;
	bool tableSaved;
	bool adjTablesSaved;

protected:
private:
};

}  // end namespace ftk

#endif	// end __ftkProjectFiles_h
