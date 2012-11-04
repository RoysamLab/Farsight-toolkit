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

#include "ftkImage.h"
#include "ftkLabelImageToFeatures.h"
#include "ftkUtils.h"

#include <iostream>

int main(int argc, char* argv[])
{
	if(argc < 4)
	{
		std::cout<<"Usage: compute_nuclei_features <InputImageFileName> <InputLabelFileName> <SegmentationResultsFileName>\n";
		return 0;
	}
    
	std::string imageName = argv[1];
	std::string labelName = argv[2];
	std::string resultsName = argv[3];
	
	ftk::Image::Pointer inpImg = ftk::Image::New();
	if(!inpImg->LoadFile(imageName))
		return 0;

	ftk::Image::Pointer labImg = ftk::Image::New();
	if(!labImg->LoadFile(labelName))
		return 0;

	ftk::IntrinsicFeatureCalculator *calc = new ftk::IntrinsicFeatureCalculator();
	calc->SetInputImages(inpImg, labImg);

	vtkSmartPointer<vtkTable> table = calc->Compute();
	ftk::SaveTable(resultsName, table);

	delete calc;

	return 1;
}
