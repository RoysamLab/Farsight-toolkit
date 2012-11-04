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

#include "ftkNuclearSegmentation.h"

#include <iostream>

int main(int argc, char* argv[])
{
	if(argc != 3)
	{
		std::cout<<"Usage: classify_nuclei <xmlResultFile> <ClassFile>\n";
		return 0;
	}

	std::string xmlFullName = argv[1];

	size_t found = xmlFullName.find_last_of(".");
	std::string xmlBaseName = xmlFullName.substr(0,found);

	std::string classFile = argv[2];

	ftk::NuclearSegmentation *segmentation = new ftk::NuclearSegmentation();
	//segmentation->RestoreFromXML(xmlFullName);
	//segmentation->LoadClassInfoFromFile(classFile);
	//segmentation->SaveChanges(xmlFullName);

	//Now load up the label image and split it into a separate image for each class
	//segmentation->SaveLabelByClass();
	
	delete segmentation;

	return 1;
}
