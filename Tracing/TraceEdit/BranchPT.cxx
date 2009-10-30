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
#include "TraceLine.h"
#include "TraceBit.h"
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
void branchPT::parent(int id)
{
	unsigned int found = 0, i = 0;
	while ( !found && (i < this->connected.size() ))
	{
		if (id== this->connected.at(i)->GetId())
		{
			this->branch = this->connected.at(i);
			this->connected.erase(this->connected.begin() + i - 1);
			this->solved = true;
			found = true;
		}//end if
		else
		{
			i++;
		}//end else
	}//end while
}
bool branchPT::state()
{
	return this->solved;
}