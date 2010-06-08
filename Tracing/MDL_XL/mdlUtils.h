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
#ifndef __mdlUtils_h
#define __mdlUtils_h

#include "mdlTypes.h"
#include <string>
#include <vector>
#include <iostream>
#include <stdio.h>

namespace mdl
{

class vtkFileHandler
{
public:
	vtkFileHandler();
	void SetNodes(std::vector<fPoint3D> * nodes);
	void SetLines(std::vector<pairE> * lines);

	bool Write(std::string filename);
	bool Read(std::string filename);
	
	bool GetNodesandLinesFromVtkfile(std::string filename);

    std::vector<fPoint3D>getNodes();
    std::vector<pairE>getLines();

private:

	std::vector<fPoint3D> Nodes;
	std::vector<pairE> Lines;	

	std::vector<fPoint3D> * nodesPtr;
	std::vector<pairE> * linesPtr;

};

}  // end namespace mdl

#endif

