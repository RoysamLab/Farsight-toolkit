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
#ifndef __BRANCHPT_H
#define __BRANCHPT_H
#include <vector>
#include "TraceLine.h"
#include "TraceBit.h"

class branchPT
{
public:
	branchPT();
  ~branchPT();
	branchPT(TraceBit branchBit, std::vector<TraceLine*> connected);
	bool SeekParent(int id);
	TraceLine* getParent();
	TraceBit GetBit();
	void SetBit(TraceBit b);
	bool state();
	std::vector<TraceLine*> GetChildren();
	std::vector<int> childIDS();
  void AddConnection(TraceLine *line);

private:
	TraceBit branchBit;
	TraceLine* branch;
	std::vector<TraceLine*> connected;
	std::vector<TraceLine*> children;
	bool solved;
};
#endif
