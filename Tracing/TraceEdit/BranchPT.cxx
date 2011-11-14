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

#include <vector>
#include "branchPT.h"

branchPT::branchPT()
{
	this->solved = false;
}

branchPT::branchPT(TraceBit branchBit, std::vector<TraceLine*> connected)
{
	this->branchBit = branchBit;
	this->connected = connected;
	this->solved = false;
}

//it's the TraceObject's responsibility to delete these lines.
branchPT::~branchPT()
{
  this->connected.clear();
  this->children.clear();
}

bool branchPT::SeekParent(int id)
{	
	if (!this->solved)
	{
		this->children.clear();
		for (unsigned int i = 0; i < this->connected.size(); i++)
		{
			if (id== this->connected.at(i)->GetId())
			{
				this->branch = this->connected.at(i);			
				this->solved = true;			
			}//end if
			else
			{
				this->children.push_back(this->connected.at(i));
			}//end else
		}//end while
	}
	return this->solved;
}
TraceBit branchPT::GetBit()
{
	return this->branchBit;
}
void branchPT::SetBit(TraceBit b)
{
  this->branchBit = b;
}
TraceLine* branchPT::getParent()
{
	return this->branch;
}
std::vector<TraceLine*> branchPT::GetChildren()
{
	return this->children;
}
bool branchPT::state()
{
	return this->solved;
}
std::vector<int> branchPT::childIDS()
{
	std::vector<int> ids;
	unsigned int i;
	for ( i = 0; i < this->children.size(); i++)
	{
		ids.push_back( this->children.at(i)->GetId());
	}
	return ids;
}

void branchPT::AddConnection(TraceLine *line)
{
  this->connected.push_back(line);
}

