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

#include "nuclear_segmentation/yousef_core/yousef_seg.h"
#include "nuclear_segmentation/gui/nuclei.h"

#include <iostream>

using namespace std;

int main(unsigned int argc, char* argv[])
{
	if(argc < 3)
	{
		std::cout<<"Usage: nucseg2 <InputImageFileName> <ParametersFileName>\n";
		return 0;
	}
	
	FTKAbstractSegmentation *segmentation = new nuclei();
	segmentation->setup(argv[1], argv[2]);
	segmentation->executeModule(0);	
	segmentation->executeModule(1);
	segmentation->executeModule(2);
	segmentation->save();
	segmentation->generateXML();

	segmentation->initMetaNeural();

	delete segmentation;

	return 1;
}
