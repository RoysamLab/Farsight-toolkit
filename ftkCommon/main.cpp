/*=========================================================================
Copyright 2009 Rensselaer Polytechnic Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 
=========================================================================*/
#include "ftkMultipleImageHandler.h"

int main(int argc, char *argv[])
{
	if(argc < 4)
    {
		std::cerr << "Usage: " << argv[0] << "pattern(string) start(int) end(int) <outputBase(string)> <outputDir(string)> <x(int)> <y(int)> <z(int)> <color>" << std::endl;
		return EXIT_FAILURE;
    }

	ftk::MultipleImageHandler * iHandle = new ftk::MultipleImageHandler();
	
	std::string inPattern = argv[1];
	int start = atoi(argv[2]);
	int end = atoi(argv[3]);
	std::string outBase = "b";
	std::string outDir = "";
	int x = 1;
	int y = 1;
	int z = 1;
	int color = 1;

	if(argc >= 5)
		outBase = argv[4];

	if(argc >= 6)
		outDir = argv[5];
	if(argc >= 7)
		x = atoi(argv[6]);
	if(argc >= 8)
		y = atoi(argv[7]);
	if(argc >= 9)
		z = atoi(argv[8]);
	if(argc >= 10)
		color = atoi(argv[9]);

	iHandle->SetOutputDirectory(outDir);
	iHandle->SetOutputBase(outBase);
	iHandle->SeriesToBlocks(inPattern, start, end, x, y, z, color);
	//iHandle->SeriesProjection("%03d.tif", 001, 110, "p.tif");

	delete iHandle;


	return EXIT_SUCCESS;
}
