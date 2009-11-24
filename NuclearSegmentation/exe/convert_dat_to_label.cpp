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
	if(argc < 3)
	{
		std::cout<<"Usage: convert_dat_to_label <InputImageFileName> <InputdatlFileName>\n";
		return 0;
	}
    
	std::string imageName = argv[1];
	std::string datName = argv[2];
	
	ftk::NuclearSegmentation *segmentation = new ftk::NuclearSegmentation();	

	std::cout<<" Reading from dat segmentation output...";
	segmentation->LoadInput(imageName);
	segmentation->LoadFromDAT(imageName);
	std::cout<<"done"<<std::endl;
	
	segmentation->SaveLabelImage();
	
	delete segmentation;

	return 1;
}
