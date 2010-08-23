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
#if defined(_MSC_VER)
#pragma warning(disable : 4996)
#endif

#include "mdlUtils.h"

namespace mdl
{

vtkFileHandler::vtkFileHandler()
{
	nodesPtr = NULL;
	linesPtr = NULL;
}

void vtkFileHandler::SetNodes(std::vector<fPoint3D> * nodes)
{
	nodesPtr = nodes;
}
void vtkFileHandler::SetLines(std::vector<pairE> * lines)
{
	linesPtr = lines;
}

bool vtkFileHandler::Write(std::string filename)
{
	if(!nodesPtr)
		return false;

	FILE * fout = fopen(filename.c_str(), "w");
	if (fout == NULL)
		return false;

	fprintf(fout, "# vtk DataFile Version 3.0\n");
	fprintf(fout,"Tracing Result\n");
	fprintf(fout,"ASCII\n");
	fprintf(fout,"DATASET POLYDATA\n");

	int num_nodes = (int)nodesPtr->size();
	fprintf(fout,"POINTS %d float\n",num_nodes);
	for(int i=0; i<num_nodes; ++i)
	{
		mdl::fPoint3D nd = nodesPtr->at(i);
		fprintf(fout,"%f %f %f\n", (float)nd.x, (float)nd.y, (float)nd.z);
	}

	if(linesPtr)
	{
		int num_lines = (int)linesPtr->size();
		fprintf(fout,"LINES %d %d\n", num_lines, num_lines*3);
		for(int i=0; i<num_lines; ++i)
		{
			pairE e = linesPtr->at(i);
			fprintf(fout, "2 %d %d\n", e.first, e.second);
		}
	}
	else
	{
		fprintf(fout,"VERTICES %d %d\n", num_nodes, num_nodes*2);
		for(int i=0; i<num_nodes; ++i)
		{
			fprintf(fout,"%d %d\n", 1, i);
		}
	}

	fclose(fout);
	
	return true;
}

bool vtkFileHandler::Read(std::string filename)
{
	if(!nodesPtr)
		return false;

	FILE * infile = fopen(filename.c_str(),"rb");
	if(infile == NULL)
		return false;

	//FIND NUMBER OF POINTS IN FILE:
	char str[200];
	int num_nodes = 0;
	while(num_nodes == 0)
    {
		if( fscanf(infile,"%s",str) == EOF )
		{
			std::cerr << "end of file!" << std::endl;
			return false;
		}
		if( strcmp(str, "POINTS") == 0) //found POINTS
		{
			fscanf(infile,"%s",str);	//get next string (number)
			num_nodes = atoi(str);
			fscanf(infile,"%s",str);	//get last string ( type )
		}
    }

	//Get the points
	nodesPtr->clear();

	float temp;
	for(int i=0; i<num_nodes; i++)
	{
		fPoint3D p;
		if( fscanf (infile,"%f",&temp) == EOF )
		{
			std::cerr << "end of file!" << std::endl;
			return false;
		}
		p.x = temp;
		if( fscanf (infile,"%f",&temp) == EOF )
		{
			std::cerr << "end of file!" << std::endl;
			return false;
		}
		p.y = temp;
		if( fscanf (infile,"%f",&temp) == EOF )
		{
			std::cerr << "end of file!" << std::endl;
			return false;
		}
		p.z = temp;
		nodesPtr->push_back(p);
	}


	if(!linesPtr)
		return true;

	//FIND NUMBER OF LINES IN FILE:
	int num_lines = 0;
	while(num_lines == 0)
    {
		if( fscanf(infile,"%s",str) == EOF )
		{
			std::cerr << "end of file!" << std::endl;
			return true;
		}
		if( strcmp(str, "LINES") == 0) //found LINES
		{
			fscanf(infile,"%s",str);	//get next string (number)
			num_lines = atoi(str);
			fscanf(infile,"%s",str);	//get last string ( type )
		}
    }

	//Get the lines
	linesPtr->clear();

	int tmp;
	for(int i=0; i<num_lines; i++)
	{
		pairE p;
		if( fscanf (infile,"%d",&tmp) == EOF )
		{
			std::cerr << "end of file!" << std::endl;
			return true;
		}
		if( fscanf (infile,"%d",&tmp) == EOF )
		{
			std::cerr << "end of file!" << std::endl;
			return true;
		}
		p.first = tmp;
		if( fscanf (infile,"%d",&tmp) == EOF )
		{
			std::cerr << "end of file!" << std::endl;
			return true;
		}
		p.second = tmp;
		linesPtr->push_back(p);
	}
	
	fclose(infile);

	return true;
}

/*
	void SetNodes(std::vector<fPoint3D> * nodes);
	void SetLines(std::vector<pairE> * lines);
	bool Write(std::string filename);
	*/

void SWCFileWriter::SetNodes(std::vector<fPoint3D> *nodes)
{
	nodesPtr = nodes;
}

void SWCFileWriter::SetLines(std::vector<pairE> *lines)
{
	linesPtr = lines;
}

}  // end namespace mdl
