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

#include <iostream>
#include <string>
#include <vector>
#include "View3D.h"

struct SeedPT{
	double x,y,z;
	int id;
};
std::vector<SeedPT> seeds;
bool readPts(char*filename)
{
	FILE * fp = fopen(filename, "r");
	if(fp==NULL)
	{
		printf("Couldn't open file %s for parsing\n",filename);
		return false;
	}
	char buff[1024];
	int id, parent;
	double x,y,z;
	while(!feof(fp))
	{
		if(fgets(buff,1024,fp)==NULL)
			break;
		int pc = 0;
		while(buff[pc]==' '&&(pc<1023))
			pc++;
		if(buff[pc]=='#') // ignoring comment lines for now. We have to read scale from it! TODO
			continue;
		sscanf(buff,"%d %lf %lf %lf  %d",&id,&x,&y,&z,&parent);
		SeedPT newSeed;
		newSeed.x=x;	newSeed.y=y;	newSeed.z=z;	newSeed.id=id;
		seeds.push_back(newSeed);
	}
	fclose(fp);
	printf("size of list= %i", seeds.size());
	return true;
};
int main (int argc, char* argv[])	{
	View3d view; //initalizes the 3d viewer
	view.RenderWin();
	//view.readImg(argv[1]);
	view.rayCast(argv[1]);
	view.AddVolumeSliders();
	readPts(argv[2]);
	//view.AddPointsAsPoints(seeds);
	view.renWin->Render();
	view.iren->Start();

}